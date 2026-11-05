// Magma side of the Magma/Oscar benchmark: a server that bench.py drives over a pipe.
//
//   magma -b cases:=CASES spec:=/path/to/lattices.spec [memlimit:=GB] bench_magma.m
//
// Reads the cases file (bench.py sample), warms up every task on small fixed
// lattices in this process (the TimeoutCall children inherit what Magma initialises
// lazily, see WarmLatticeCaches), prints "@@READY", then answers one job per stdin
// line
//     JOB <case id> <task> <variant> <timeout> <min_time> <max_reps>
// with one line
//     @@RESULT <case id> <task> <variant> <status> <K> <wall,...> <cpu,...> <extra> <value>
// (tab separated; status ok | timeout | crash | error | unsupported; times in
// microseconds per call).  Anything else on stdout is log output, which bench.py
// ignores.  Each job runs in a HardKill TimeoutCall child, so a timeout or a memory
// overrun (memlimit) kills only that job.
//
// Timing: Magma's clocks tick in 10 ms steps, so a call is timed in blocks of K
// calls, with the task's untimed setup done beforehand; K grows tenfold until a
// block takes BLOCK seconds (or the untimed preparation of a block gets long).
// Blocks of the final K are repeated until min_time seconds or max_reps blocks; the
// times reported are per call.  Every call in a job gets its own Gram matrix U G U^t
// (see scramble and next_grams): while a lattice is referenced, LatticeWithGram
// returns that same object for an equal Gram matrix, so the K lattices of a block
// would share their cached minimum, short vectors etc.
//
// Feed the jobs through a pipe, not a file: a TimeoutCall child's exit rewinds a
// shared file offset, and the server then reads its input again (the fork gotcha in
// lattices/CLAUDE.md).  Run it from an empty directory: the orbit-method enumeration
// reads and writes a disk cache under orth_cache/ in the working directory, and a
// warm cache would make class numbers look cheaper than they are from scratch.

if not assigned cases then error "usage: magma -b cases:=CASES spec:=SPEC bench_magma.m"; end if;
if not assigned spec then spec := "../lattices.spec"; end if;
T0 := Realtime();
SetColumns(0);         // no line wrapping in the replies
AttachSpec(spec);
if not assigned memlimit then memlimit := "0"; end if;
if StringToInteger(memlimit) gt 0 then SetMemoryLimit(StringToInteger(memlimit) * 2^30); end if;

Z := Integers();
BLOCK := 0.2;          // seconds a timed block should take
PREP_LIMIT := 5.0;     // stop growing K once preparing a block takes this long
// Largest K: rank 2 has few distinct bases with small entries (see next_grams).
KMAX := func< n | n le 2 select 1000 else 10000 >;
MAX_ATTEMPTS := 200;

CaseFmt := recformat< id, n, nplus, gram, gram_iso, gram_other, prec, mass >;

function parse_gram(n, s)
    if s eq "-" then return 0; end if;
    return Matrix(Z, n, n, [StringToInteger(x) : x in Split(s, ",")]);
end function;

function read_cases(path)
    lines := [l : l in Split(Read(path), "\n") | #l gt 0];
    // Split drops empty fields; bench.py writes "-" for an empty one, never "".
    header := Split(lines[1], "\t");
    col := func< name | Index(header, name) >;
    cs := AssociativeArray();
    for line in lines[2..#lines] do
        f := Split(line, "\t");
        error if #f ne #header, "malformed case line", line;
        n := StringToInteger(f[col("rank")]);
        m := f[col("mass")];
        cs[f[col("id")]] := rec< CaseFmt | id := f[col("id")], n := n,
            nplus := StringToInteger(f[col("nplus")]), gram := parse_gram(n, f[col("gram")]),
            gram_iso := parse_gram(n, f[col("gram_iso")]), gram_other := parse_gram(n, f[col("gram_other")]),
            prec := StringToInteger(f[col("theta_prec")]), mass := m eq "-" select 0 else eval m >;
    end for;
    return cs;
end function;

function lattice(G)
    return LatticeWithGram(G : CheckPositive := false);
end function;

// U G U^t for the unimodular U made by max(2n, 12) row operations U[i] +:= e U[j] drawn from
// the high bits of a linear congruential generator seeded with k (its low bits have
// short periods); k = 0 gives G itself.
// bench_oscar.jl has the same function, so both systems see the same bases.
function scramble(G, k)
    n := Nrows(G);
    if k eq 0 or n lt 2 then return G; end if;
    U := IdentityMatrix(Z, n);
    s := k;
    for t in [1..Max(2*n, 12)] do
        s := (1103515245 * s + 12345) mod 2^31; i := (s div 2^16) mod n + 1;
        s := (1103515245 * s + 12345) mod 2^31; j := (s div 2^16) mod (n - 1) + 1;
        if j ge i then j +:= 1; end if;
        s := (1103515245 * s + 12345) mod 2^31; e := IsEven(s div 2^16) select 1 else -1;
        AddRow(~U, e, j, i);
    end for;
    return U * G * Transpose(U);
end function;

// The Gram matrices for the next call of a job: scramble(c`gram, k) for the next k
// whose matrix (and its partner's, for the isometry tasks) this job has not used
// yet, or a repeated one after MAX_ATTEMPTS tries (a small lattice with a large
// automorphism group has few); used holds their Eltseqs.  Partner M is 0 for the
// other tasks.
procedure next_grams(task, c, ~k, ~used, ~G, ~M)
    attempt := 0;
    repeat
        k +:= 1;
        attempt +:= 1;
        G := scramble(c`gram, k);
        M := task eq "isometric_same" select scramble(c`gram_iso, k) else
             (task eq "isometric_other" select scramble(c`gram_other, k) else 0);
    until attempt ge MAX_ATTEMPTS or
          (Eltseq(G) notin used and (M cmpeq 0 or Eltseq(M) notin used));
    Include(~used, Eltseq(G));
    if M cmpne 0 then Include(~used, Eltseq(M)); end if;
end procedure;

function definite(c) return c`nplus eq c`n; end function;

// ---------------------------------------------------------------------------
// Tasks.  For each task and variant: applies(c), an untimed setup on the fresh
// lattice returning a state, and the timed run returning <value, extra> as strings
// in the format the Oscar server uses (value "@@UNSUPPORTED" when that is decided
// only at run time, "@@TIMEOUT" when an inner time limit ran out).  Variant
// "default" is what the pipeline does; "stock" is Magma's own intrinsic where the
// pipeline replaced it.

VARIANTS := AssociativeArray();
VARIANTS["class_number"] := ["default", "stock"];
VARIANTS["mass"] := ["default"];
VARIANTS["aut_order"] := ["default", "stock"];
VARIANTS["aut_label"] := ["default"];
VARIANTS["minimum"] := ["default"];
VARIANTS["kissing"] := ["default"];
VARIANTS["successive_minima"] := ["default"];
VARIANTS["root_system"] := ["default"];
VARIANTS["theta"] := ["default"];
VARIANTS["orth_decomp"] := ["default", "stock"];
VARIANTS["isometric_same"] := ["default"];
VARIANTS["isometric_other"] := ["default"];
VARIANTS["disc_image"] := ["default"];
VARIANTS["oq_order"] := ["default"];      // Oscar only: Magma has no O(q) for torsion forms

function applies(task, variant, c)
    if task eq "oq_order" then return false; end if;
    if task eq "class_number" then return true; end if;
    if task in ["isometric_same", "isometric_other"] then
        G := task eq "isometric_same" select c`gram_iso else c`gram_other;
        // Magma's IsIsometric needs a definite form; the pipeline's binary test
        // (genus_reps.m) covers indefinite rank 2.  Indefinite rank >= 3 has none.
        return G cmpne 0 and (definite(c) or c`n eq 2);
    end if;
    return definite(c);
end function;

// The pipeline runs AutomorphismGroupFaster on an LLL-reduced Gram matrix (fill's
// canonical or LLL gram), and it is very slow on unreduced bases: E8+E8 in bases
// with entries up to ~100 took over 60 s in 9 of 10 cases, against 0.07 s once
// reduced (stock AutomorphismGroup: 0.02 s either way).  So the default variant
// reduces first, as fill does.
function reduced(L)
    return LatticeWithGram(LLLGram(GramMatrix(L)));
end function;

// M: the partner's Gram matrix for the isometry tasks (see next_grams)
function setup(task, variant, L, c, M)
    if task in ["aut_label", "disc_image"] then
        Lr := reduced(L);
        return <Lr, AutomorphismGroupFaster(Lr)>;
    elif task in ["isometric_same", "isometric_other"] then
        return lattice(M);
    end if;
    return 0;
end function;

function seq_string(s) return Join([Sprint(x) : x in s], ","); end function;

function run(task, variant, L, st, c, timeout)
    case task:
    when "class_number":
        if variant eq "stock" then
            return <Sprint(#GenusRepresentatives(L)), "">;
        end if;
        // The pipeline's dispatch (fill_genus.m), with its own time limit.
        ok, reps := GenusReps(L : Timeout := timeout,
                              UseOrth := definite(c) and UseOrthHeuristic(c`n, c`mass));
        return <ok select Sprint(#reps) else "@@TIMEOUT", "">;
    when "mass":
        return <Sprint(Mass(L)), "">;
    when "aut_order":
        G := variant eq "stock" select AutomorphismGroup(L) else AutomorphismGroupFaster(reduced(L));
        return <Sprint(#G), "">;
    when "aut_label":
        if not CanIdentifyGroup(#st[2]) then return <"@@UNSUPPORTED", "">; end if;
        id := IdentifyGroup(st[2]);
        return <Sprintf("%o.%o", id[1], id[2]), "">;
    when "minimum":
        return <Sprint(Minimum(L)), "">;
    when "kissing":
        return <Sprint(KissingNumber(L)), "">;
    when "successive_minima":
        m, _ := SuccessiveMinima(L);
        return <seq_string(m), "">;
    when "root_system":
        return <RootString(L : IsRootSystem := false), "">;
    when "theta":
        f := ThetaSeries(L, c`prec);
        return <seq_string([Coefficient(f, i) : i in [0..c`prec]]), "">;
    when "orth_decomp":
        D := variant eq "stock" select OrthogonalDecomposition(L) else OrthogonalDecompositionFaster(L);
        return <Join([Sprintf("%o:%o", t[1], t[2]) : t in Sort([<Rank(M), Abs(Determinant(M))> : M in D])], ","), "">;
    when "isometric_same", "isometric_other":
        iso := definite(c) select IsIsometric(L, st) else IsIsometricBinary(L, st);
        return <iso select "true" else "false", "">;
    when "disc_image":
        // As aut_data in fill_genus.m: the image of O(L) in the automorphism group of
        // the discriminant group.  extra records |Aut(L^#/L)|, which fill divides by
        // the image order for festi_veniani_index (Oscar's extra is |O(q_L)|).
        Lr, aut := Explode(st);           // the group acts in Lr's basis
        LD := Dual(Lr : Rescale := false);
        A, quo := LD / Lr;
        AA := AutomorphismGroup(A);
        gens := [A.i : i in [1..Ngens(A)]];
        im := [hom< A -> A | [quo(x @@ quo * ChangeRing(alpha, Rationals())) : x in gens] > : alpha in Generators(aut)];
        return <Sprint(#sub< AA | im >), Sprintf("AutA=%o", #AA)>;
    end case;
    error "unknown task", task;
end function;

usec := func< x | Sprint(Round(10^6 * x)) >;

// CPU seconds (user + system) of this process, from /proc: Cputime() stays at 0 in a
// TimeoutCall child, where every job runs.  Ticks are 1/100 s, as on Linux.
function proc_cpu()
    s := Read("/proc/self/stat");
    i := Index(s, ")");                   // the command name may contain spaces
    f := Split(s[i+2..#s], " ");          // f[1] is field 3 (state); utime, stime are 14, 15
    return (StringToInteger(f[12]) + StringToInteger(f[13])) / 100;
end function;

// One block: K lattices on fresh Gram matrices (next_grams; k and used carry over
// between the blocks of a job) with their setups, untimed, then the K timed calls.
// Sets the values, the wall and CPU seconds of the timed part, and the preparation
// time.
procedure block(task, variant, c, K, timeout, ~k, ~used, ~vals, ~wall, ~cpu, ~prep)
    t := Realtime();
    Ls := [];
    sts := [* *];
    for i in [1..K] do
        next_grams(task, c, ~k, ~used, ~G, ~M);
        L := lattice(G);
        Append(~Ls, L);
        Append(~sts, setup(task, variant, L, c, M));
    end for;
    prep := Realtime(t);
    t := Realtime(); u := proc_cpu();
    vals := [run(task, variant, Ls[i], sts[i], c, timeout) : i in [1..K]];
    wall := Realtime(t); cpu := proc_cpu() - u;
end procedure;


function measure(task, variant, c, timeout, min_time, max_reps)
    if not applies(task, variant, c) then return "unsupported\t0\t-\t-\t-\t-"; end if;
    K := 1;
    k := -1;               // the first call gets k = 0, the Gram matrix as given
    used := {PowerSequence(Z) | };
    block(task, variant, c, K, timeout, ~k, ~used, ~vals, ~wall, ~cpu, ~prep);
    value := vals[1][1];
    extra := vals[1][2];
    if value eq "@@UNSUPPORTED" then return "unsupported\t0\t-\t-\t-\t-"; end if;
    if value eq "@@TIMEOUT" then return "timeout\t0\t-\t-\t-\t-"; end if;
    // Grow the block until it is long enough to time with a 10 ms clock.
    while wall lt BLOCK and K lt KMAX(c`n) and prep lt PREP_LIMIT do
        K := Min(KMAX(c`n), 10 * K);
        block(task, variant, c, K, timeout, ~k, ~used, ~vals, ~wall, ~cpu, ~prep);
        error if exists{v : v in vals | v[1] ne value}, "value changed between calls";
    end while;
    walls := [wall / K];
    cpus := [cpu / K];
    while #walls lt max_reps and K * &+walls lt min_time do
        block(task, variant, c, K, timeout, ~k, ~used, ~vals, ~wall, ~cpu, ~prep);
        error if exists{v : v in vals | v[1] ne value}, "value changed between calls";
        Append(~walls, wall / K);
        Append(~cpus, cpu / K);
    end while;
    return Sprintf("ok\t%o\t%o\t%o\t%o\t%o", K, Join([usec(w) : w in walls], ","),
                   Join([usec(u) : u in cpus], ","), extra eq "" select "-" else extra, value);
end function;

function clean(s)
    if #s eq 0 then return "-"; end if;
    s := &cat[ch in ["\t", "\n", "\r"] select " " else ch : ch in Eltseq(s)];
    return #s gt 300 select s[1..300] else s;
end function;

// Runs in the TimeoutCall child; returns the reply fields after <task> <variant>.
// Errors are caught here: Call (Parallel.magma) would print them to stdout.
function job(task, variant, c, timeout, min_time, max_reps)
    try
        return measure(task, variant, c, timeout, min_time, max_reps);
    catch e
        return "error\t0\t-\t-\t-\t" * clean(e`Object);
    end try;
end function;

// One untimed call, to initialise in this process what the task initialises lazily.
procedure warm_once(task, variant, c)
    if applies(task, variant, c) then
        L := lattice(c`gram);
        M := task eq "isometric_same" select c`gram_iso else c`gram_other;
        _ := run(task, variant, L, setup(task, variant, L, c, M), c, 60);
    end if;
end procedure;

procedure reply(id, task, variant, fields)
    printf "@@RESULT\t%o\t%o\t%o\t%o\n", id, task, variant, fields;
end procedure;

// Small fixed lattices covering the code paths (definite rank 2, 3, 5 and 8, odd and
// even, class numbers 1 and 3; indefinite rank 2 and 3), the same ones
// bench_oscar.jl warms up on.
function warm_cases()
    mk := func< id, n, s, g, gi, go, p | rec< CaseFmt | id := id, n := n, nplus := s,
        gram := Matrix(Z, n, n, g), gram_iso := gi cmpeq 0 select 0 else Matrix(Z, n, n, gi),
        gram_other := go cmpeq 0 select 0 else Matrix(Z, n, n, go), prec := p, mass := 0 > >;
    E8 := Eltseq(CartanMatrix("E8"));
    return [
        mk("warm1", 3, 3, [2,-1,0,-1,2,0,0,0,6], [2,1,0,1,2,0,0,0,6], [2,0,0,0,2,0,0,0,6], 20),
        mk("warm2", 8, 8, E8, E8, 0, 6),
        mk("warm3", 2, 2, [1,0,0,5], [1,1,1,6], [2,1,1,3], 20),
        mk("warm4", 3, 3, [1,0,0,0,1,0,0,0,3], 0, 0, 20),
        mk("warm5", 2, 1, [2,1,1,-2], [2,-1,-1,-2], 0, 0),
        mk("warm6", 3, 2, [2,1,0,1,2,0,0,0,-18], 0, 0, 0),
        mk("warm7", 3, 3, [2,0,0,0,2,-1,0,-1,16], [2,0,0,0,2,1,0,1,16], 0, 40),
        mk("warm8", 5, 5, [2,-1,1,-1,1,-1,2,-1,0,-1,1,-1,4,1,2,-1,0,1,4,-1,1,-1,2,-1,6], 0, 0, 12)];
end function;

procedure() // read the whole script before the first fork (see lattices/CLAUDE.md)
    cs := read_cases(cases);
    t_load := Realtime(T0);
    t := Realtime();
    WarmLatticeCaches();
    for c in warm_cases() do
        for task in Keys(VARIANTS) do
            for variant in VARIANTS[task] do
                try
                    warm_once(task, variant, c);
                catch e
                    printf "warm-up %o %o/%o: %o\n", c`id, task, variant, clean(e`Object);
                end try;
            end for;
        end for;
    end for;
    t_warm := Realtime(t);
    a, b, cc := GetVersion();
    printf "@@READY\t%o\t%o\tMagma V%o.%o-%o\n", t_load, t_warm, a, b, cc;
    F := Open("/dev/stdin", "r");
    while true do
        line := Gets(F);
        if IsEof(line) then break; end if;
        f := Split(line, "\t");
        if #f lt 7 or f[1] ne "JOB" then continue; end if;
        id := f[2]; task := f[3]; variant := f[4];
        timeout := StringToInteger(f[5]);
        min_time := eval f[6];
        max_reps := StringToInteger(f[7]);
        if not IsDefined(cs, id) or not IsDefined(VARIANTS, task) or variant notin VARIANTS[task] then
            reply(id, task, variant, "error\t0\t-\t-\t-\tunknown case, task or variant");
            continue;
        end if;
        // The child gets a little more than the job's timeout, which also bounds the
        // inner TimeoutCall of GenusReps.
        t := Realtime();
        ok, out := TimeoutCall(timeout + 5, job, <task, variant, cs[id], timeout, min_time, max_reps>, 1 : HardKill := true);
        if ok then
            reply(id, task, variant, out[1]);
        elif Realtime(t) ge timeout then
            reply(id, task, variant, "timeout\t0\t-\t-\t-\t-");
        else
            // the child died early: a segfault, or the memory limit
            reply(id, task, variant, "crash\t0\t-\t-\t-\t-");
        end if;
    end while;
    exit 0;
end procedure();
