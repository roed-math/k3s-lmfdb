declare verbose FillGenus, 1;

function hecke_primes(rank)
    if rank lt 8 then
        return [2,3,5];
    else
        return [2];
    end if;
end function;

intrinsic StringToReal(s::MonStgElt) -> RngIntElt
{ Converts a decimal string (like 123.456 or 1.23456e40 or 1.23456e-10) to a real number at default precision. }
    if #s eq 0 then return 0.0; end if;
    if "e" in s then
        t := Split(s,"e");
        require #t eq 2: "Input should have the form 123.456e20 or 1.23456e-10";
        return StringToReal(t[1])*10.0^StringToInteger(t[2]);
    end if;
    t := Split(s,".");
    require #t le 2: "Input should have the form 123 or 123.456 or 1.23456e-10";
    n := StringToInteger(t[1]);  s := t[1][1] eq "-" select -1 else 1;
    return #t eq 1 select RealField()!n else RealField()!n + s*RealField()!StringToInteger(t[2])/10^#t[2];
end intrinsic;

intrinsic ThetaSeriesIncremental(L::Lat, target_prec::RngIntElt, timeout::RngIntElt : Min := 0) -> SeqEnum, RngIntElt, Assoc
{Compute the theta series coefficients of L, doubling the precision until either target_prec is reached or the time budget timeout (in seconds) is exhausted. Returns the coefficient sequence obtained, the precision actually reached, and an associative array mapping each precision computed to the time (in seconds) it took to compute the theta series to that precision.  Min is the minimum of L if known (0 if not); it only sets the starting precision, and is a parameter so that this never computes the minimum itself outside a TimeoutCall.}
    best_theta := [];
    best_prec := 0;
    theta_elapsed := AssociativeArray();
    remaining := timeout;
    prec := Maximum(16, Min + 4);
    while prec le target_prec and remaining gt 0 do
        current_prec := Minimum(prec, target_prec);
        success, theta, elapsed := TimeoutCall(remaining, ThetaSeries, <L, current_prec - 1>, 1 : HardKill := true);
        if not success then
            vprintf FillGenus, 1 : "Theta series timed out at precision %o\n", current_prec;
            break;
        end if;
        best_theta := Eltseq(theta[1]);
        best_prec := current_prec;
        // ThetaSeries recomputes the whole series each time, so this elapsed time
        // is the cost of computing theta to current_prec (which is what SetHashes
        // uses to cost theta-based hashing).
        theta_elapsed[current_prec] := StringToReal(elapsed);
        vprintf FillGenus, 1 : "Theta series to precision %o in %o s\n", current_prec, elapsed;
        if current_prec ge target_prec then break; end if;
        remaining -:= Ceiling(StringToReal(elapsed));
        prec *:= 2;
    end while;
    return best_theta, best_prec, theta_elapsed;
end intrinsic;

function dict_to_jsonb(dict)
    return "{" * Join([Sprintf("\"%o\":%o", key, dict[key]) : key in Keys(dict)], ",") * "}";
end function;

intrinsic to_postgres(val::Any : jsonb_val := false) -> Any
{Serialise a value (matrix, sequence, tuple, associative array or scalar) into the Postgres array / JSON text format used for the database export.}
    delims := jsonb_val select "[]" else "{}";
    if ISA(Type(val),Mtrx) then
        return to_postgres(Eltseq(val) : jsonb_val:=jsonb_val);
    elif val cmpeq "\\N" then
        return val;
    // I think this is unnecessary, and used to cause problems, so removing for now.
    //elif Type(val) eq MonStgElt then
    //    return "\"" * val * "\""; // This will fail if the string has quotes, but I don't think that's ever true for us.
    elif Type(val) in [SeqEnum, Tup] then
        return delims[1] * Join([Sprintf("%o",to_postgres(x : jsonb_val:=jsonb_val)) : x in val],",") * delims[2];
    elif Type(val) eq Assoc then
        val_prime := AssociativeArray();
        for key in Keys(val) do
            val_prime[to_postgres(key)] := to_postgres(val[key] : jsonb_val:=true);
        end for;
        return dict_to_jsonb(val_prime);
    else
        return val;
    end if;
end intrinsic;

function RescaledDualNF(L)
    Q := Rationals();
    K := BaseRing(L);
    B := ChangeRing(BasisMatrix(L),Q);
    M := ChangeRing(InnerProductMatrix(L),Q);
    F := ChangeRing(GramMatrix(L),Q);
    B := F^-1 * B;
    B := IntegralMatrix(B);
    B div:= GCD(Eltseq(B));
    F := B * M * Transpose(B);
    F, d := IntegralMatrix(F);
    g := GCD(Eltseq(F));
    F div:= g;
    M := (d/g) * M;
    return NumberFieldLattice(Rows(ChangeRing(B, K)) : InnerProduct := ChangeRing(M,K));
end function;

function hecke_data_of(G, p)
    // Wrapped so it can run under TimeoutCall: AdjacencyMatrix walks the p-neighbour graph
    // in a C-level builtin that does not yield to Magma's Alarm, so it needs the HardKill
    // watchdog to be bounded (it can take >75s for a single rank-8 genus).  The
    // characteristic polynomial and its factorisation are dense in the class number, so
    // they run in the same child.
    Ap := AdjacencyMatrix(G, p);
    fpf := Factorization(CharacteristicPolynomial(Ap));
    return Ap, [(<Coefficients(pair[1]), pair[2]>) : pair in fpf];
end function;

// Per-lattice computations that run in HardKill children (see the loop in FillGenus).
// The batch process itself must make no call that can run unbounded: nothing outside
// can interrupt it, so one stuck builtin there hangs the whole batch.  The children
// are grouped because at low rank the fork, not the computation, is the cost: a
// HardKill TimeoutCall costs ~7-15 ms of wall-clock in a fill batch, against well
// under 1 ms for each of these at rank 4.

// The minimum and what is read off it (Magma caches it on L, so the kissing number,
// density and Hermite number reuse it), the root system, and the minimum of the
// rescaled dual D, which sets where its theta series starts.
function minimum_data(L, D)
    m := Minimum(L);
    roots := m ge 3 select "0" else RootString(L : IsRootSystem := false);
    return m, KissingNumber(L), Density(L), HermiteNumber(L), roots, Minimum(D);
end function;

function root_string_of(L)
    return RootString(L : IsRootSystem := false);
end function;

function successive_minima_of(L)
    minima, _ := SuccessiveMinima(L);
    return minima;
end function;

// The order of the image of O(L) in the automorphism group of the discriminant group
// L^#/L, and the SmallGroups label of the automorphism group ("\N" when IdentifyGroup
// cannot handle the order).  The image lies in O(q_L), and the Festi-Veniani index is
// [O(q_L) : image]; Magma has no orthogonal groups of torsion quadratic forms, so
// oscar/fix_festi_veniani.jl computes |O(q_L)| and the index from the image order,
// which fill records in festi_veniani/ (the index used to divide #Aut(L^#/L), the
// automorphism group of the abelian group, instead of |O(q_L)|).  The discriminant
// group's automorphism group order and invariants are returned so that the caller can
// check them against the genus: an assert in the child would only make the call fail
// quietly.
function aut_data(Lcanon, aut_group, disc_aut_size)
    LD := Dual(Lcanon : Rescale:=false);
    discL, quo := LD/Lcanon;
    disc_aut := AutomorphismGroup(discL);
    gens_disc := [discL.i : i in [1..Ngens(discL)]];
    im_aut := [hom< discL -> discL | [quo(x@@quo*ChangeRing(alpha, Rationals())): x in gens_disc]> : alpha in Generators(aut_group)];
    aut_label := "\\N";
    if CanIdentifyGroup(#aut_group) then
        Aid := IdentifyGroup(aut_group);
        aut_label := Sprintf("%o.%o", Aid[1], Aid[2]);
    end if;
    return #sub<disc_aut | im_aut>, #disc_aut, Invariants(discL), aut_label;
end function;

function SphereVolume(n)
    RR := RealField();
    pi := Pi(RR);
    m := n div 2;
    if IsEven(n) then
        return pi^m / Factorial(m);
    else
        return 2^n * pi^m * Factorial(m) / Factorial(n);
    end if;
end function;

// Coefficient i (1-based) of a stored theta series.  ThetaSeriesIncremental stores
// Eltseq of the power series, which drops trailing zero coefficients, so a list
// can be shorter than its precision (common for rescaled duals, whose scale makes
// most coefficients zero); the missing coefficients are 0, as FirstDiff and
// ThetaHash in hash.m already assume.  Indexing past the end crashed the sort:
// "Sequence index 182 should be in the range 1 to 181" on 8.8.504.001.01.3827.
function coeff(f, i)
    return i le #f select f[i] else 0;
end function;

// Once per process, before the first genus (FillGenus and ConnectGenus both call it):
// run each kind of per-lattice computation on two fixed small lattices, here in the
// batch process.  Magma initialises a lot lazily on first use (the SmallGroups
// database, real constants, factorisation and enumeration machinery), and the
// children FillGenus forks inherit whatever the batch process has initialised, but
// discard what they initialise themselves.  When the
// minimum, density, IdentifyGroup and the Hecke polynomials were computed in the
// batch process this happened as a side effect; once they moved into HardKill
// children, every child paid it again, and fill took ~1.5x as long at rank <= 4.
// The inputs are fixed and tiny, so this is bounded (~0.1 s).
fill_warm := NewStore();

intrinsic WarmLatticeCaches()
{Initialise, once per process, what the per-lattice computations of fill and connect
initialise lazily, so that the TimeoutCall children they fork inherit it.}
    if StoreIsDefined(fill_warm, "done") then return; end if;
    E8 := LatticeWithGram(CartanMatrix("E8"));
    W2 := LatticeWithGram(DiagonalJoin(CartanMatrix("A2"), Matrix(Integers(), 1, 1, [6])));
    for W in [E8, W2] do
        _ := Minimum(W); _ := KissingNumber(W); _ := Density(W); _ := HermiteNumber(W);
        _ := SuccessiveMinima(W); _ := RootString(W : IsRootSystem := false);
        _ := ThetaSeries(W, 20); _ := ThetaSeries(Dual(W), 20);
        _ := CanonicalForm(GramMatrix(W));
        _ := GroupToString(AutomorphismGroupFaster(W) : use_id := false);
        // connect
        _ := IsIsomorphic(W, LatticeWithGram(GramMatrix(W)));
        _ := OrthogonalDecompositionFaster(W);
        _ := BVHash(W, 1, 3);
        _ := EvenSublattice(W);
    end for;
    _ := IdentifyGroup(SmallGroup(12, 3));
    _ := Factorization(CharacteristicPolynomial(Matrix(Integers(), 2, 2, [1, 2, 3, 4])));
    _ := #AutomorphismGroup(AbelianGroup([2, 6]));
    StoreSet(fill_warm, "done", true);
end intrinsic;

// Two optional integer keys, larger first.  A key that was not computed ("\N": its
// TimeoutCall overran or failed) sorts after every computed one.  Comparing only when both are
// present, as cmp_lat used to, is not transitive once some keys are missing
// (a ~ b ~ c with a < c), and Sort then returns an order that depends on the order
// the classes were enumerated in.
function cmp_int_desc(a, b)
    ka := Type(a) eq RngIntElt;
    kb := Type(b) eq RngIntElt;
    if ka and kb then return b - a; end if;
    if ka ne kb then return ka select -1 else 1; end if;
    return 0;
end function;

// Two stored theta series, coefficient by coefficient, smaller first.  A coefficient
// at or past a series' precision is unknown and sorts after every known one, so
// that two series are compared as fixed sequences and the order stays transitive
// when the precisions differ.
function cmp_theta(f1, p1, f2, p2)
    if Type(f1) ne SeqEnum then p1 := 1; end if;
    if Type(f2) ne SeqEnum then p2 := 1; end if;
    for i in [1..Maximum(p1, p2) - 1] do
        k1 := i lt p1;
        k2 := i lt p2;
        if k1 ne k2 then return k1 select -1 else 1; end if;
        d := coeff(f1, i) - coeff(f2, i);
        if d ne 0 then return d; end if;
    end for;
    return 0;
end function;

// Sort order of the lattices of a genus (records as built in FillGenus); file level
// so tests/test_fill_sort.m can exercise it.  The order fixes the labels.  Every key
// but the last can be missing when its computation timed out, so labels are fixed by
// the run that enumerates the genus, not by the genus alone; a missing key sorts
// after every present one, which keeps this a total preorder.
function cmp_lat(L1, L2)
    n := Isqrt(#L1["gram"]);
    // 1. automorphism group size (larger first)
    d := cmp_int_desc(L1["aut_size"], L2["aut_size"]);
    if d ne 0 then return d; end if;
    // 2. density (denser first); within a genus the determinant is fixed, so we
    //    sort by the minimum as a proxy (larger minimum = denser)
    d := cmp_int_desc(L1["minimum"], L2["minimum"]);
    if d ne 0 then return d; end if;
    // 3. root system (larger first, see CompareRootStrings)
    r1 := L1["root_sublattice"];
    r2 := L2["root_sublattice"];
    if r1 ne "\\N" and r2 ne "\\N" then
        d := CompareRootStrings(r1, r2);
        if d ne 0 then return d; end if;
    elif r1 ne r2 then
        return r1 ne "\\N" select -1 else 1;
    end if;
    // 4. theta series
    d := cmp_theta(L1["theta_series"], L1["theta_prec"], L2["theta_series"], L2["theta_prec"]);
    if d ne 0 then return d; end if;
    // 5. dual theta series
    d := cmp_theta(L1["dual_theta_series"], L1["dual_theta_prec"], L2["dual_theta_series"], L2["dual_theta_prec"]);
    if d ne 0 then return d; end if;
    // 6. tiebreaker: the Gram matrix entries (canonical when CanonicalForm finished,
    //    else the LLL-reduced Gram matrix, which depends on the enumeration)
    for i in [1..n^2] do
        d := L1["gram"][i] - L2["gram"][i];
        if (d ne 0) then return d; end if;
    end for;
    return 0;
end function;

intrinsic FillGenus(label::MonStgElt : timeout := 1800, masslimit := 0, sizelimit := 0, timelimit := 0, adjlimit := 0, perlattice := 0, use_orth := true, force_orth := false, stats := false, oscar := "oscar_genus")
{Fill the data for a genus and its lattice representatives, given files in the genera_basic format.

Enumeration guards (0 = unlimited): masslimit skips enumerating a definite genus whose
mass exceeds it (its class number is too large to enumerate); adjlimit is a work budget
for the adjacency (Hecke) matrix -- skip it when the estimated work, ~class_number times
sum_p p^(rank-2) over the Hecke primes p, exceeds adjlimit;
sizelimit records the genus-level data but skips storing individual lattices past that
class number; timelimit is a per-genus wall-clock cap on the per-lattice loop.  In every
case the genus-level record is still written, so a guarded genus is bounded rather than
either running for hours or going missing.

perlattice is the wall-clock each individual lattice gets for its own computations
(canonical form, automorphism group, theta series).  It used to be derived as
timeout/class_number, which starved exactly the genera that needed it most -- at
timeout 300 a class-number-55 genus gave each lattice 6s, while CanonicalForm alone
measures 8-66s at rank 12, so those fields were silently written as "\N".  Note this
trades against timelimit: the per-lattice loop costs roughly class_number*perlattice,
and a genus that overruns timelimit has ALL of its lattices discarded, so raising
perlattice without raising timelimit turns partial data into no data.  Left at 0 it
falls back to the old timeout/class_number budget.

use_orth selects the orbit-method
enumeration for definite genera of rank >= 3 (see GenusReps in genus_reps.m).
stats writes one fill_stats/ record per genus (fill_stats.format): which guard, if any,
decided the outcome, the class number, per-phase wall-clock and the slowest single
lattice -- the data for calibrating the guards (calibrate_fill.py).}
    WarmLatticeCaches();
    genus_t0 := Realtime();
    data := Split(Split(Read(LabelPath("genera_basic", label)), "\n")[1], "|");
    basic_format := Split(Read("genera_basic.format"), "|");
    advanced_format := Split(Read("genera_advanced.format"), "|");
    hash_format := Split(Read("lat_hash.format"), "|");
    // This function only fills in basic lattice entries (essentially those that don't require interactions between different genera)
    lat_format := Split(Split(Read("lat_basic.format"), "\n")[1], "|");
    assert #data eq #basic_format;
    basics := AssociativeArray();
    for i in [1..#data] do
        basics[basic_format[i]] := data[i];
        if data[i] eq "None" then basics[basic_format[i]] := "\\N"; end if;
    end for;
    advanced := AssociativeArray();
    lats := [];

    n := StringToInteger(basics["rank"]);
    s := StringToInteger(basics["nplus"]);
    K := Rationals();
    LWG := LatticeWithGram;
    DualLat := Dual;
    rep := basics["rep"];
    // Switch to square brackets
    rep := "[" * rep[2..#rep - 1] * "]"; // Switch to square brackets
    gram0 := Matrix(K, n, eval rep);
    L0 := LWG(gram0 : CheckPositive := false);
    vprintf FillGenus, 1 : "Computing genus representatives...";
    reps := [];
    // masslimit: a definite genus of mass M has class number on the order of M, so a
    // very large mass means enumeration is hopeless -- skip it and keep only the
    // (cheap) genus-level record.  Mass is only defined for definite genera; indefinite
    // ones store "\N" here, so only parse it when definite.
    genus_mass := 0;
    if (n eq s) and (basics["mass"] ne "\\N") then
        mass_pieces := Split(basics["mass"][2..#basics["mass"]-1], ",");
        genus_mass := StringToInteger(mass_pieces[1]) / StringToInteger(mass_pieces[2]);
    end if;
    spinor := 0;   // the spinor generators, for indefinite rank >= 3 (see GenusReps)
    // Classes enumerated with Oscar (oscar/pool.jl --task enumerate, run by
    // run_all.py before fill): indefinite genera of rank >= 3, whose Magma spinor
    // kernel can overcount, and definite genera of small rank, where Oscar is faster.
    // An ok record replaces GenusReps; its automorphism group orders and root systems
    // are kept per class (oscar_aut, oscar_root) to cross-check Magma's and to stand in
    // when Magma's computation overruns.  A failed record falls back to Magma for a
    // definite genus but not for an indefinite one of rank >= 3, whose Magma class
    // number may be wrong.
    oscar_status := "none";
    oscar_aut := [];
    oscar_root := [];
    opath := LabelPath(oscar, label);
    if OpenTest(opath, "r") then
        olines := Split(Read(opath), "\n");
        ohead := AssociativeArray();
        for kv in Split(olines[1], " ") do
            pos := Index(kv, "=");
            if pos gt 0 then ohead[kv[1..pos-1]] := kv[pos+1..#kv]; end if;
        end for;
        oscar_status := IsDefined(ohead, "status") select ohead["status"] else "bad";
        if oscar_status eq "ok" then
            oreps := [];
            for line in olines[2..#olines] do
                f := Split(line, "|");
                Append(~oreps, GramStringToLat(f[1], n));
                Append(~oscar_aut, f[2] eq "\\N" select 0 else StringToInteger(f[2]));
                Append(~oscar_root, f[3]);
            end for;
            error if #oreps ne StringToInteger(ohead["classes"]), Sprintf("%o: %o classes in %o, header says %o", label, #oreps, opath, ohead["classes"]);
            if ohead["spinor_generators"] ne "\\N" then
                spinor := [StringToInteger(x) : x in Split(ohead["spinor_generators"][2..#ohead["spinor_generators"]-1], ",")];
            end if;
        end if;
    end if;
    from_oscar := oscar_status eq "ok";
    no_magma := (not from_oscar) and oscar_status ne "none" and (n ne s) and (n ge 3);
    skip_massive := (not from_oscar) and (masslimit gt 0) and (n eq s) and (genus_mass gt masslimit);
    if from_oscar then
        vprintf FillGenus, 1 : "Classes from %o\n", opath;
        genus_success, reps, orth_prov := Explode(<true, oreps, 0>);
    elif no_magma then
        vprintf FillGenus, 1 : "Oscar enumeration %o; no Magma fallback for an indefinite genus of rank >= 3\n", oscar_status;
        genus_success := false;
    elif skip_massive then
        vprintf FillGenus, 1 : "Skipping enumeration: mass exceeds masslimit %o\n", masslimit;
        genus_success := false;
    else
        // Method selection, per-method timing and fallbacks live in GenusReps
        // (genus_reps.m); use_orth toggles the orbit method for definite rank >= 3,
        // refined by the rank/mass heuristic (the neighbour walk stays the default
        // for small class numbers at low rank, where it wins).
        genus_success, reps, orth_prov, spinor := GenusReps(L0 : Timeout := timeout,
            UseOrth := use_orth and (force_orth or UseOrthHeuristic(n, genus_mass)));
    end if;
    t_enum := Realtime(genus_t0);
    enum_method := from_oscar select "oscar" else no_magma select "oscar_" * oscar_status else
                   skip_massive select "none" else ((n eq s and use_orth and n ge 3 and (force_orth or UseOrthHeuristic(n, genus_mass))) select "orth" else "default");
    adj_status := "none";
    advanced["class_number"] := "\\N";
    // Indefinite rank >= 3: the classes are labelled by their spinor coordinates
    // relative to the published rep (genus_reps_Spinor), recorded per lattice in
    // spinor_coordinates; this records the primes they refer to.
    advanced["spinor_generators"] := (spinor cmpne 0) select to_postgres(spinor) else "\\N";
    advanced["adjacency_matrix"] := "\\N";
    advanced["adjacency_polynomials"] := "\\N";
    // ambient_genus: the genus of the parent lattices from which this genus's
    // classes were obtained as orthogonal complements (orbit-method provenance);
    // the per-lattice counter/vector references live in lat_lattices.
    if not assigned orth_prov then orth_prov := 0; end if;
    advanced["ambient_genus"] := (orth_prov cmpeq 0 or orth_prov[1] eq "") select "\\N" else orth_prov[1];
    have_adjacency := false;   // set below iff the adjacency (Hecke) matrix is computed
    if genus_success then
        vprintf FillGenus, 1 : "Number of genus representatives: %o\n", #reps;
        advanced["class_number"] := #reps;
        vprintf FillGenus, 1 : "Computing adjacency matrix for p = ";
        hecke_mats := AssociativeArray();
        hecke_polys := AssociativeArray();
        G := Genus(L0);
        // This works for 2.28 - should be replaced by SetGenus in 2.29
        G`Representatives := reps;
        G`IsNatural := true;
        // The adjacency (Hecke) matrix walks the p-neighbour graph: each of the #reps
        // representatives has ~p^(n-2) p-neighbours (the isotropic lines in L/pL), each
        // needing an isometry test.  So the work is ~ #reps * sum_p p^(n-2), which grows
        // as p^(n-2) with rank -- a single rank-12 genus costs 16x a rank-8 one even at
        // class number 1.  We estimate that work and skip the matrix past adjlimit (a
        // work budget), an adaptive cutoff that captures both the rank and class-number
        // dependence, rather than a crude cap on either alone.  The class number is still
        // recorded, and pneighbors below is emitted only when the matrix exists.
        adj_work := #reps * (&+[Integers() | p^Maximum(n-2, 0) : p in hecke_primes(n)]);
        have_adjacency := (n eq s) and ((adjlimit le 0) or (adj_work le adjlimit));
        if have_adjacency then
          adj_ok := true;
          for p in hecke_primes(n) do
            vprintf FillGenus, 1 : "%o:", p;
            // Bound the adjacency computation: it is an uninterruptible builtin that can
            // run for well over a minute on a single high-rank genus, so it needs the
            // HardKill watchdog.  If it overruns we drop the (optional) Hecke data for
            // this genus -- the class number is already recorded -- rather than let a few
            // genera dominate wall-clock.
            adj_success, adj_out, adj_elapsed := TimeoutCall(timeout, hecke_data_of, <G, p>, 2 : HardKill := true);
            if not adj_success then
                vprintf FillGenus, 1 : "(timed out after ~%os, skipping Hecke data)\n", timeout;
                adj_ok := false;
                break;
            end if;
            vprintf FillGenus, 1 : "computed in %o seconds\n", adj_elapsed;
            hecke_mats[p] := adj_out[1];
            hecke_polys[p] := adj_out[2];
          end for;
          adj_status := adj_ok select "ok" else "timeout";
          if adj_ok then
            vprintf FillGenus, 1 : "Done!\n";
            advanced["adjacency_matrix"] := to_postgres(hecke_mats);
            advanced["adjacency_polynomials"] := to_postgres(hecke_polys);
          else
            have_adjacency := false;   // no Hecke data -> the pneighbors block below is skipped
          end if;
        end if;
    else
        reps := [];
    end if;
    // sizelimit: past this class number, keep the genus-level record but skip the
    // per-lattice loop below (individual automorphism groups, theta series, hashes).
    t_adj := Realtime(genus_t0) - t_enum;
    size_skipped := (sizelimit gt 0) and (#reps gt sizelimit);
    if size_skipped then
        vprintf FillGenus, 1 : "Class number %o exceeds sizelimit %o: storing genus-level data only\n", #reps, sizelimit;
        reps := [];
    end if;
    disc_invs := basics["discriminant_group_invs"];
    disc_invs := "[" * disc_invs[2..#disc_invs-1] * "]"; // Switch to square brackets
    disc_invs := eval disc_invs;
    disc_aut_size := #AutomorphismGroup(AbelianGroup(disc_invs));

    if (n eq s) then
        vprintf FillGenus, 1 : "Computing canonical forms and automorphism groups for representative ";
    end if;

    if (#reps gt 0) then
        // PER-OPERATION timeout, NOT a per-lattice aggregate: to_per_rep is passed
        // independently to CanonicalForm, AutomorphismGroupFaster, the primal theta series,
        // the dual theta series and each of the smaller guarded calls below (minimum,
        // root system, successive minima, ...), so a single lattice can spend several
        // multiples of to_per_rep.  The aggregate wall-clock guard is the per-GENUS
        // timelimit checked at the top of this loop, not this value.
        // perlattice gives every operation the same budget regardless of class number; the
        // old timeout/#reps split a fixed budget across the representatives, which starved
        // precisely the large genera (class number 55 at timeout 300 => 6s per operation,
        // against a measured 8-66s just for CanonicalForm at rank 12) and silently wrote the
        // results as "\N".  Kept as the fallback when perlattice is 0.
        to_per_rep := (perlattice gt 0) select perlattice else timeout div #reps + 1;
    end if;

    // Theta-series timing per precision, aggregated over the genus representatives.
    // DECISION: aggregate by max time per precision (the slowest representative
    // bounds the cost of hashing the whole genus to that precision).  This may
    // change in future -- e.g. to the mean or to per-representative timings.
    theta_elapsed := AssociativeArray();

    // timelimit: a per-genus wall-clock backstop.  If the per-lattice work overruns it
    // (a genus under sizelimit whose individual lattices are still pathologically slow),
    // stop and fall back to the genus-level record rather than grinding for hours.
    genus_timed_out := false;
    reps_done := 0;
    max_lat_time := 0.0;
    t_loop0 := Realtime();
    for Li->L in reps do
        t_lat0 := Realtime();
        if (timelimit gt 0) and (Realtime(genus_t0) gt timelimit) then
            vprintf FillGenus, 1 : "Per-genus timelimit %o s exceeded at rep %o/%o: storing genus-level data only\n", timelimit, Li, #reps;
            genus_timed_out := true;
            break;
        end if;
        lat := AssociativeArray();
        lat["lattice"] := L; // useful for subroutines; removed before saving to disk
        for col in ["rank", "nplus", "nminus", "disc_abs", "disc_sign", "disc_radical", "disc_witt", "disc_geometric", "disc_quadratic", "disc_half", "disc_2adic_unit", "bad_primes", "discriminant_group_invs", "discriminant_group_exponent", "is_even", "level", "scale", "conway_symbol", "dual_conway_symbol"] do
            lat[col] := basics[col];
        end for;
        lat["genus_label"] := basics["label"];
        lat["class_number"] := advanced["class_number"];
        // TODO (Eran) := The code for ConwaySymbol is currently in sage.
        // The magma implemntation is in version 2.29 that has some bugs
        // This is no longer part of the lattice, only of the genus
        // lat["dual_conway"] := "\\N";
        lat["aut_size"] := "\\N";
        lat["festi_veniani_index"] := "\\N";   // filled in by oscar/fix_festi_veniani.jl
        lat["fv_image"] := "\\N";              // |image of O(L) in Aut(L^#/L)|, for festi_veniani/
        lat["aut_label"] := "\\N";
        lat["aut_group"] := "\\N";
        lat["is_chiral"] := "\\N";
        // Orbit-method provenance: this lattice is the orthogonal complement of
        // orthogonal_complement (a vector) inside the lattice with counter
        // ambient_lattice in the genus ambient_genus (see genera_advanced).  The
        // vector lives in the first lattice up the parent chain for which a Gram
        // matrix is available; \N when the class came from a p-neighbour walk or
        // the parent genus has no recorded label/counters.
        lat["ambient_lattice"] := "\\N";
        lat["orthogonal_complement"] := "\\N";
        if orth_prov cmpne 0 and orth_prov[1] ne "" and Li le #orth_prov[2] and orth_prov[2][Li][1] gt 0 then
            lat["ambient_lattice"] := orth_prov[2][Li][1];
            lat["orthogonal_complement"] := orth_prov[2][Li][2];
        end if;
        lat["density"] := "\\N";
        lat["hermite"] := "\\N";
        lat["kissing"] := "\\N";
        lat["minimum"] := "\\N";
        lat["theta_series"] := "\\N";
        lat["theta_prec"] := "\\N";
        lat["dual_theta_series"] := "\\N";
        lat["dual_theta_prec"] := "\\N";
        lat["successive_minima"] := "\\N";
        lat["root_sublattice"] := "\\N";
        // reps are in label order for these genera, so position Li is label .Li
        lat["spinor_coordinates"] := (spinor cmpne 0) select SpinorCoordinates(Li, #spinor) else "\\N";
        // Only set below in the positive-definite (n eq s) branch; default them here
        // so the indefinite case still has every column.
        lat["center_density"] := "\\N";
        lat["canonical_gram"] := "\\N";
        lat["is_universal"] := "\\N";
        lat["is_even_universal"] := "\\N";
        // Trying to reduce the size of the entries in the gram matrix
        gram0 := GramMatrix(L);
        gram := LLLGram(gram0);
        max_abs := Max([Abs(x) : x in Eltseq(gram)]);
        max_abs_0 := Max([Abs(x) : x in Eltseq(gram0)]);
        if max_abs_0 le max_abs then
            gram := gram0;
        end if;
        lat["gram"] := Eltseq(gram);
        lat["gram_is_canonical"] := false;
        lat["gram_others"] := []; // This will be manually set in cases like E8 where we want to store other options
        // At the moment we do not have a notion of a canonical gram in the indefinite case
        if (n eq s) then
            vprintf FillGenus, 1 : "%o", gram;
            // The automorphism group and dual must be expressed in the same basis as
            // the stored gram (ConnectGenus rebuilds the lattice from that gram, and
            // reuses lat["aut_group"] unchanged against it -- see VoronoiData/
            // AutOrbitCoords/tDesign in connect_genus.m). Lcanon must therefore start
            // out matching lat["gram"] (= gram, LLL-reduced above), NOT L (gram0's
            // basis) -- gram0 and gram can differ by a basis change LLLGram picked,
            // in which case an automorphism group computed on L would not preserve
            // gram at all (confirmed: on a scrambled-basis E8 test case, ALL of its
            // generators failed to preserve the stored gram). If CanonicalForm below
            // succeeds, Lcanon and lat["gram"] both get overwritten together to the
            // canonical lattice, staying in sync the same way.
            Lcanon := LatticeWithGram(gram);
            success, canonical_gram, elapsed := TimeoutCall(to_per_rep, CanonicalForm, <gram>, 1 : HardKill := true);
            vprintf FillGenus, 1 : "Canonical form computed in %o seconds\n", elapsed;
            if success then
                canonical_gram := canonical_gram[1];
                lat["gram"] := Eltseq(canonical_gram);
                lat["canonical_gram"] := Eltseq(canonical_gram);
                lat["gram_is_canonical"] := true;
                Lcanon := LatticeWithGram(canonical_gram);
                lat["lattice"] := Lcanon;
            end if;
            success, aut_group, elapsed := TimeoutCall(to_per_rep, AutomorphismGroupFaster, <Lcanon>, 1 : HardKill := true);
            vprintf FillGenus, 1 : "Automorphism group computed in %o seconds\n", elapsed;
            if success then
                aut_group := aut_group[1];
                lat["aut_group"] := GroupToString(aut_group : use_id:=false);
                lat["aut_size"] := #aut_group;
                lat["is_chiral"] := &and[Determinant(g) eq 1 : g in Generators(aut_group)];
                // the discriminant group is checked against the genus here, and is
                // also what the Festi-Veniani index needs
                ok, ad := TimeoutCall(to_per_rep, aut_data, <Lcanon, aut_group, disc_aut_size>, 4 : HardKill := true);
                if ok then
                    assert disc_aut_size eq ad[2];
                    assert disc_invs eq ad[3];
                    lat["fv_image"] := ad[1];
                    lat["aut_label"] := ad[4];
                end if;
            end if;
            // Oscar's automorphism group order (isometry_group): a cross-check, and the
            // value when Magma's computation overran (aut_size is a sort key, cmp_lat)
            if from_oscar and oscar_aut[Li] gt 0 then
                if Type(lat["aut_size"]) eq RngIntElt then
                    error if lat["aut_size"] ne oscar_aut[Li],
                        Sprintf("%o class %o: #Aut %o from Magma, %o from Oscar", label, Li, lat["aut_size"], oscar_aut[Li]);
                else
                    lat["aut_size"] := oscar_aut[Li];
                end if;
            end if;
            D := Dual(L);
            // LLL always finishes (polynomial time), so it is safe here in the batch
            // process, and it caches a reduced basis on L and D that every child below
            // inherits; without it each theta-series child redoes that work (~4 ms per
            // child at rank 4, measured -- it was the unguarded Minimum here that used
            // to warm this cache).
            _ := LLL(L);
            _ := LLL(D);
            m := 0;    // the minimum, 0 while unknown
            dm := 0;   // the dual's minimum, likewise
            ok, md := TimeoutCall(to_per_rep, minimum_data, <L, D>, 6 : HardKill := true);
            if ok then
                m := md[1];
                lat["minimum"] := m;
                lat["kissing"] := md[2];
                lat["density"] := md[3];
                lat["center_density"] := md[3] / SphereVolume(n);
                lat["hermite"] := md[4];
                lat["root_sublattice"] := md[5];
                dm := md[6];
            elif not (from_oscar and oscar_root[Li] ne "\\N") then
                // The root system is a sort key (cmp_lat), so try it on its own:
                // enumerating vectors of norm <= 2 can be much cheaper than finding
                // the minimum.
                ok, rs := TimeoutCall(to_per_rep, root_string_of, <L>, 1 : HardKill := true);
                if ok then lat["root_sublattice"] := rs[1]; end if;
            end if;
            // Oscar's root system (root_lattice_recognition, same string format): a
            // cross-check, and the value when Magma's overran
            if from_oscar and oscar_root[Li] ne "\\N" then
                if lat["root_sublattice"] ne "\\N" then
                    error if lat["root_sublattice"] ne oscar_root[Li],
                        Sprintf("%o class %o: root system %o from Magma, %o from Oscar", label, Li, lat["root_sublattice"], oscar_root[Li]);
                else
                    lat["root_sublattice"] := oscar_root[Li];
                end if;
            end if;
            target_prec := Max(150, m+4);
            theta, theta_prec, rep_theta_elapsed := ThetaSeriesIncremental(L, target_prec, to_per_rep : Min := m);
            for p in Keys(rep_theta_elapsed) do
                theta_elapsed[p] := IsDefined(theta_elapsed, p)
                    // Maximum for worst case Hash computation.
                    // Change to sum if we care about the average case
                    select Maximum(theta_elapsed[p], rep_theta_elapsed[p]) else rep_theta_elapsed[p];
            end for;
            for p in Keys(theta_elapsed) do
                if not IsDefined(rep_theta_elapsed, p) then theta_elapsed[p] := Infinity(); end if;
            end for;
            lat["is_universal"] := "\\N";
            lat["is_even_universal"] := "\\N";
            if theta_prec gt 0 then
                lat["theta_series"] := theta;
                lat["theta_prec"] := theta_prec;
                // 15-theorem for is_universal, doubled 290-theorem for
                // is_even_universal (sources, and why 2*S290 is needed even for odd
                // lattices, in universality.m).  The previous scan certified
                // is_even_universal from theta support on 2..30, which is only valid
                // for doubly even lattices.  Norms beyond theta_prec are checked by
                // a TimeoutCall-guarded first-vector search; undecided stays "\N".
                udec, uval, eudec, euval := UniversalityFlags(L : Theta := theta);
                if udec then
                    lat["is_universal"] := uval;
                end if;
                if eudec then
                    lat["is_even_universal"] := euval;
                end if;
            else
                lat["theta_series"] := [1];
                lat["theta_prec"] := 1;
            end if;
            // Dual theta series (used as a tie-breaker in cmp_lat below, and stored
            // alongside the primal theta series).
            dual_theta, dual_theta_prec := ThetaSeriesIncremental(D, Max(150, dm + 4), to_per_rep : Min := dm);
            if dual_theta_prec gt 0 then
                lat["dual_theta_series"] := dual_theta;
                lat["dual_theta_prec"] := dual_theta_prec;
            else
                lat["dual_theta_series"] := [1];
                lat["dual_theta_prec"] := 1;
            end if;
            ok, minima := TimeoutCall(to_per_rep, successive_minima_of, <L>, 1 : HardKill := true);
            if ok then
                lat["successive_minima"] := minima[1];
            end if;

        end if;
        lat["hash"] := "\\N";

        //lat["level"] := Level(LatticeWithGram(ChangeRing(GramMatrix(L), Integers()) : CheckPositive:=false));

        Append(~lats, lat);
        reps_done +:= 1;
        max_lat_time := Maximum(max_lat_time, Realtime(t_lat0));
    end for;
    t_lattices := Realtime(t_loop0);
    // per-operation timeouts (or failures) inside the loop, counted before a timelimit
    // discard empties lats; definite only (indefinite lattices never get either)
    canon_nulls := (n eq s) select #[lat : lat in lats | not lat["gram_is_canonical"]] else 0;
    aut_nulls := (n eq s) select #[lat : lat in lats | Type(lat["aut_size"]) ne RngIntElt] else 0;

    // If the per-genus timelimit tripped mid-loop, discard the partial lattices so the
    // record is consistent (class_number set, genus-level data only) rather than
    // claiming a class number it did not store.
    if genus_timed_out then lats := []; end if;

    vprintf FillGenus, 1 : "Done!\n";


    // Label order (cmp_lat): automorphism group size, minimum, root system, theta
    // series, dual theta series, then the Gram matrix.  Any of the first five can be
    // missing when its computation timed out, so the labels are fixed by the run that
    // enumerates the genus rather than by the genus alone.
    if (n eq s) then
        Sort(~lats, cmp_lat, ~perm);
    end if;
    // Store the Hecke matrices in label order: row and column i belong to the lattice
    // labelled .i, which is the enumerated class i^perm (Sort gives new lats[i] =
    // old lats[i^perm]).  The enumeration order is random from run to run, so a
    // matrix left in it cannot be matched to the labels.  When no lattices are stored
    // (sizelimit, timelimit) there are no labels, and it stays in enumeration order.
    if have_adjacency and #lats gt 0 then
        for p in Keys(hecke_mats) do
            A := hecke_mats[p];
            N := Nrows(A);
            hecke_mats[p] := Matrix(Integers(), N, N, &cat[[A[i^perm, j^perm] : j in [1..N]] : i in [1..N]]);
        end for;
        advanced["adjacency_matrix"] := to_postgres(hecke_mats);
    end if;

    SetColumns(0);
    for idx->L in lats do
        // Need label for lattice.  counter records this sort position explicitly, so
        // cross-genus references (ambient_lattice in child genera) can point at it.
        lats[idx]["label"] := Sprintf("%o.%o", basics["label"], idx);
        lats[idx]["counter"] := idx;
    end for;
    // Persist the finished genus (reps in counter order) to the on-disk orth cache:
    // in rank-descending orchestration the genera of the next-lower rank descend
    // from these lattices and resolve their ambient_lattice counters against them.
    if (n eq s) and (#lats gt 0) then
        WriteGenusOrthCache(basics["label"], [lat["lattice"] : lat in lats]);
    end if;

    if #lats gt 0 then
        // the Hecke matrices (in label order, above) let the theta-Hecke pass skip
        // recomputing the adjacency of a genus SetHashes leaves undecided
        SetHashes(~lats, ~advanced, theta_elapsed, timeout :
                  Adjacency := have_adjacency select hecke_mats else AssociativeArray(), Sidecar := "theta_hecke");
    else
        // GenusRepresentatives failed (e.g. timed out on a hard genus): record the
        // genus with no stored lattices.  HashGenus factors through the canonical
        // label, so the genus hash is that of the label -- no representative needed.
        advanced["genus_hash"] := CollapseIntList(StringToBytes(basics["label"]));
        advanced["theta_distinguishing_prec"] := "\\N";
        advanced["is_theta_distinguished"] := "\\N";
        advanced["hash_function"] := "\\N";
        advanced["is_hash_distinguished"] := "\\N";
    end if;
    // We need to be able to look up hash functions for lattices that are not in the main
    // genus being processed.  So we write the hash function used to a separate file
    // so that it can be looked up when needed (see lookup_hash_function in connect_genus.m)
    Write(LabelPath("genera_hash", n, s, IntegerToString(advanced["genus_hash"]) : Create), advanced["hash_function"] : Overwrite);

    // festi_veniani/<rank>/<nplus>/<genus>: "<lattice label>|<image order>" for every
    // lattice whose image of O(L) on the discriminant group was computed; the index
    // itself is \N here and oscar/fix_festi_veniani.jl fills it in.
    fv_lines := [Sprintf("%o|%o", lat["label"], lat["fv_image"]) : lat in lats | Type(lat["fv_image"]) eq RngIntElt];
    if #fv_lines gt 0 then
        Write(LabelPath("festi_veniani", basics["label"] : Create), Join(fv_lines, "\n") : Overwrite);
    end if;
    for idx->L in lats do
        lat := L;
        Remove(~lat, "fv_image");
        if have_adjacency then   // pneighbors are read off the Hecke matrix, so only when it was computed
            pNeighbors := AssociativeArray();
            for p in hecke_primes(n) do
                pNeighbors[p] := ["\"" * lats[j]["label"] * "\"" : j in [1..#lats] | hecke_mats[p][idx, j] ne 0];
            end for;
            lat["pneighbors"] := to_postgres(pNeighbors);
        else
            lat["pneighbors"] := "\\N";
        end if;
        Remove(~lat, "lattice");
        error if Keys(lat) ne Set(lat_format), [k : k in lat_format | k notin Keys(lat)], [k : k in Keys(lat) | k notin lat_format];
        output := Join([Sprintf("%o", to_postgres(lat[k])) : k in lat_format], "|");
        Write(LabelPath("lattice_basic_data", lat["label"] : Create), output : Overwrite);

    end for;
    // Now write hash data
    output := Join([Join([Sprintf("%o", to_postgres(lat[k])) : k in hash_format], "|") : lat in lats], "\n");
    Write(LabelPath("lattice_hashes", n, s, IntegerToString(advanced["genus_hash"]) : Create), output : Overwrite);

    error if Keys(basics) ne Set(basic_format), [k : k in basic_format | k notin Keys(basics)], [k : k in Keys(basics) | k notin basic_format];
    error if Keys(advanced) ne Set(advanced_format), [k : k in advanced_format | k notin Keys(advanced)], [k : k in Keys(advanced) | k notin advanced_format];
    output := Join([basics[k] : k in basic_format] cat [Sprintf("%o", advanced[k]) : k in advanced_format], "|");
    Write(LabelPath("genera_advanced", label : Create), output : Overwrite);
    if stats then
        outcome := skip_massive select "mass_skip" else
                   (not genus_success) select "enum_fail" else
                   size_skipped select "size_skip" else
                   genus_timed_out select "time_discard" else "ok";
        st := AssociativeArray();
        st["label"] := label;
        st["rank"] := n;
        st["nplus"] := s;
        st["disc_abs"] := basics["disc_abs"];
        st["mass"] := basics["mass"];
        st["method"] := enum_method;
        st["outcome"] := outcome;
        st["class_number"] := advanced["class_number"];
        st["lattices_stored"] := #lats;
        st["reps_processed"] := reps_done;
        st["canon_nulls"] := canon_nulls;
        st["aut_nulls"] := aut_nulls;
        st["adjacency"] := adj_status;
        st["t_enum"] := t_enum;
        st["t_adj"] := t_adj;
        st["t_lattices"] := t_lattices;
        st["max_lat_time"] := max_lat_time;
        st["t_total"] := Realtime(genus_t0);
        stats_format := Split(Split(Read("fill_stats.format"), "\n")[1], "|");
        Write(LabelPath("fill_stats", label : Create), Join([Sprintf("%o", st[k]) : k in stats_format], "|") : Overwrite);
    end if;
    return;
end intrinsic;
