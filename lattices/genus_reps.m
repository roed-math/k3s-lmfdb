// Genus enumeration methods for the lattice pipeline.
//
// This file collects the per-method wrappers and the method dispatch that used to
// live in fill_genus.m, together with a new enumeration method GenusRepresentativesOrth
// ("the orbit method") following Chenevier-Taibi, "Unimodular lattices of rank 29 and
// related even genera of small determinant" (https://otaibi.perso.math.cnrs.fr/hdiscp.pdf),
// Sect. 2.6 (Prop 2.7, Examples 2.9 and 2.13) and Sect. 7.2 (steps A1-A4).
//
// The orbit method enumerates a positive definite genus by descending from a parent
// genus of rank one higher and strictly smaller determinant: each element x of the
// discriminant group L^#/L of order o with o*(x,x) = -m mod o gives an integral
// overlattice P = < L + <mo>e, xlift + e/o > of L + <mo> with det P = det(L)*m/o.
// Conversely (Prop 2.7 / Ex 2.9), the classes of the target genus are recovered from
// the classes of Genus(P) as v^perp with v running over O(P)-orbits of primitive
// vectors of norm mo and modulus m (equivalently v = m*u with u in P^# of norm o/m).
// The chain of parents terminates in small-determinant (usually unimodular) genera.
// Every level is verified against the Minkowski-Siegel mass, and any shortfall is
// topped up with the p-neighbour walk, so the output is provably complete; the orbit
// method is an accelerator, not a source of truth.

declare verbose GenusReps, 2;
declare attributes RngInt : GenusOrthCache, GenusOrthDeadline;

import "neighbours.mag" : neighbours_old;
import "neighbours_canonical.mag" : list_genus;
import "aut_char.mag" : aut_graph;
import "isom.mag" : isom_graph;

Z := Integers();
Q := Rationals();

/****************************************************************************
 * Wrappers for the individual enumeration methods (moved from fill_genus.m)
 ****************************************************************************/

function genus_reps_Magma(L)
    // The bound is set to infinity to avoid Magma printing an error message
    // without throwing a runtime error.
    if IsPositiveDefinite(GramMatrix(L)) or (Rank(L) eq 2) then
      return GenusRepresentatives(L : Bound := Infinity());
    end if;
    // due to some bugs in Magma, we convert to number field
    LF := NumberFieldLattice(L);
    reps := GenusRepresentatives(LF);
    return [LatticeWithGram(ChangeRing(GramMatrix(r), Integers()) :
			    CheckPositive := false) : r in reps];
end function;

function genus_reps_Logan(L)
    // neighbours was renamed neighbours_old when list_genus took over; keep the
    // wrapper pointing at the surviving implementation.
    return Setseq(neighbours_old(L : thorough));
end function;

function genus_reps_Faster(L)
    // Our own mass-verified p-neighbour enumeration (definite lattices only).
    reps := GenusRepresentativesFaster(L);
    return reps;
end function;

function genus_reps_Spinor(L)
    // Indefinite lattices of rank >= 3: by Eichler each (proper) spinor genus consists
    // of a single isometry class -- but a genus may contain SEVERAL spinor genera, so
    // all of them must be enumerated.  SpinorRepresentatives(L) alone returns only the
    // class of L's own spinor genus and under-enumerates, e.g. diag(1,20,-25), whose
    // genus has two spinor genera hence two classes.  This route also avoids Magma's
    // GenusRepresentatives, which fails on some indefinite genera (e.g. signature (2,2)
    // rank 4 -- "Illegal null sequence" / "cs must be non-empty" in its LatNF code).
    return [Representative(S) : S in SpinorGenera(Genus(L))];
end function;

forward genus_reps_orth;

function genus_reps_negdef(L, use_orth)
    // Negative definite lattices: negate to positive definite, enumerate the full
    // genus there, and negate back.  (The indefinite spinor shortcut does NOT apply
    // to definite lattices, and the definite machinery requires positive gram.)
    Lp := LatticeWithGram(-GramMatrix(L));
    if Rank(Lp) ge 3 then
        reps := use_orth select genus_reps_orth(Lp) else genus_reps_Faster(Lp);
    else
        reps := genus_reps_Magma(Lp);
    end if;
    return [LatticeWithGram(-GramMatrix(R) : CheckPositive := false) : R in reps];
end function;

function genus_reps_negdef_orth(L)
    return genus_reps_negdef(L, true);
end function;

function genus_reps_negdef_walk(L)
    return genus_reps_negdef(L, false);
end function;

// Exact GL2(Z)-isometry test for the canonical rank-2 forms [[0,m],[m,k1]] and
// [[0,m],[m,k2]] of (square) determinant -m^2.  These isotropic forms have finite
// isometry groups, so we can decide isometry directly: an isometry sends a
// primitive isotropic vector of the first to one of the second, and the rest of
// the basis is then determined up to the (finitely many) sign/line choices.
function square_disc_isometric(k1, k2, m)
    G1 := Matrix(Integers(), 2, 2, [0, m, m, k1]);
    G2 := Matrix(Integers(), 2, 2, [0, m, m, k2]);
    g1 := GCD(k1, 2*m);                       // GCD(0, 2m) = 2m
    isotropic := [ Vector(Integers(), [1, 0]),
                   Vector(Integers(), [-k1 div g1, (2*m) div g1]) ];
    for v1 in isotropic, sgn in [1, -1] do
        w1 := sgn * v1;
        row := Vector(Integers(), Eltseq(Matrix(w1) * G1));   // w1^t G1
        d := GCD(row[1], row[2]);
        if d eq 0 or m mod d ne 0 then continue; end if;
        _, p, q := XGCD(row[1], row[2]);                      // p*row[1] + q*row[2] = d
        v20 := (m div d) * Vector(Integers(), [p, q]);        // <w1, v20> = m
        Qv20 := (Matrix(v20) * G1 * Transpose(Matrix(v20)))[1][1];
        if (k2 - Qv20) mod (2*m) ne 0 then continue; end if;  // Q(v20 + t w1) = Qv20 + 2 t m
        v2 := v20 + ((k2 - Qv20) div (2*m)) * w1;
        U := Matrix(Integers(), 2, 2, [w1[1], v2[1], w1[2], v2[2]]);   // columns w1, v2
        if Abs(Determinant(U)) eq 1 and Transpose(U)*G1*U eq G2 then
            return true;
        end if;
    end for;
    return false;
end function;

// Genus representatives of a rank-2 lattice L0 of square determinant -m^2.  Magma's
// GenusRepresentatives fails here (it reduces to binary quadratic forms, which
// reject square discriminants).  Every such lattice is isometric to a canonical
// [[0,m],[m,k]] with 0 <= k < 2m; we keep those in L0's genus (Genus() is
// reliable for these even though the isometry routines are not) and remove
// duplicate isometry classes with the exact test above.
function genus_reps_square_disc(L0)
    m := Isqrt(Integers() ! (-Determinant(L0)));
    G0 := Genus(L0);
    reps := [];
    rep_ks := [];
    for k in [0 .. 2*m - 1] do
        Lk := LatticeWithGram(Matrix(Rationals(), 2, 2, [0, m, m, k]) : CheckPositive := false);
        if Genus(Lk) ne G0 then continue; end if;
        if forall{ kr : kr in rep_ks | not square_disc_isometric(k, kr, m) } then
            Append(~rep_ks, k);
            Append(~reps, Lk);
        end if;
    end for;
    return reps;
end function;

/****************************************************************************
 * The orbit method
 ****************************************************************************/

// Wall-clock budget for orth work (see SetGenusOrthTimeBudget): raise a
// CATCHABLE error once the deadline passes.  Used by the batch pre-warm, which
// cannot run under TimeoutCall (the fork would take the warmed cache with it)
// and cannot rely on Alarm (uncatchable, and not delivered inside C builtins).
procedure check_orth_deadline()
    ZZ := Integers();
    if assigned ZZ`GenusOrthDeadline and Realtime() gt ZZ`GenusOrthDeadline then
        error "orth: time budget exceeded";
    end if;
end procedure;

// Canonical representative of {v, -v}: make the first nonzero entry positive.
function normsign(v)
    for c in Eltseq(v) do
        if c lt 0 then return -v; end if;
        if c gt 0 then return v; end if;
    end for;
    return v;
end function;

// Glue data for ascending from L to a rank+1 parent of smaller determinant.
// Each x in L^#/L of order o > 1 with o*(x,x) =: -m mod o, 0 < m < o, yields the
// integral overlattice P = < L + <mo>e, xlift + e/o > with det P = det(L)*m/o.
// Returns tuples <detP, o/m, o, m, xlift> sorted by detP (strongest determinant
// reduction first, avoiding long chains of weak reductions), then by the descent
// cost o/m (the norm of the dual vectors enumerated when descending back).
function glue_candidates(L)
    D := Determinant(L);
    LD := Dual(L : Rescale := false);
    A, quo := LD / L;
    cands := [];
    seen := { A | };
    for x in A do
        o := Order(x);
        if o eq 1 or x in seen then continue; end if;
        Include(~seen, -x);
        xl := x @@ quo;
        r := Z ! (o * Norm(xl));         // o*(x,x) is an integer since x in L^#
        m := (-r) mod o;
        if m eq 0 then continue; end if; // isotropic: no rank+1 determinant drop
        detP := (D * m) div o;           // o | D since o | #(L^#/L)
        // The parent's parity is known without constructing it: the glue vector has
        // norm Norm(xl) + m/o (an integer), and P is even iff L is even and that
        // norm is even.  Even parents sort first (their neighbour base cases are
        // uniformly fast), so glue_scored's construction window sees them even when
        // many odd candidates precede them in determinant order.
        oddP := (IsEven(L) and IsEven(Z ! (Norm(xl) + m/o))) select 0 else 1;
        Append(~cands, <detP, o/m, o, m, xl, oddP>);
    end for;
    cmp := function(a, b)
        if a[6] ne b[6] then return a[6] lt b[6] select -1 else 1; end if;
        if a[1] ne b[1] then return a[1] lt b[1] select -1 else 1; end if;
        if a[2] ne b[2] then return a[2] lt b[2] select -1 else 1; end if;
        return 0;
    end function;
    Sort(~cands, cmp);
    return cands;
end function;

// Build the parent overlattice for a glue candidate; returns false if the glued
// lattice is not integral (cannot happen for candidates from glue_candidates,
// but kept as a safety check).
function glue_parent(L, o, m, xl)
    n := Rank(L);
    nu := m * o;
    IP := DiagonalJoin(ChangeRing(GramMatrix(L), Q), Matrix(Q, 1, 1, [nu]));
    g := [Q | c : c in Eltseq(xl)] cat [1/o];
    M := VerticalJoin(ScalarMatrix(Q, n + 1, 1), Matrix(Q, 1, n + 1, g));
    H := HermiteForm(ChangeRing(o * M, Z));
    B := ChangeRing(Matrix([H[i] : i in [1 .. n + 1]]), Q) / o;
    P := Lattice(B, IP);
    GP := GramMatrix(P);
    if exists{ e : e in Eltseq(GP) | not IsIntegral(e) } then
        return false, _;
    end if;
    P := LatticeWithGram(ChangeRing(GP, Z));
    assert Determinant(P) eq (Determinant(L) * m) div o;
    return true, P;
end function;

// Representatives for the O(P)-orbits of type-(o, m) vectors in P: primitive v
// with v.v = m*o and v.P = m*Z, i.e. v = m*u with u in P^# of norm o/m, u
// primitive in P^#, and v primitive in P.  Enumeration happens in the dual
// (Remark 7.4 in Chenevier-Taibi), where the norm o/m is small.  Returns
// false if more than max_vecs vectors would have to be considered.
function typed_vector_orbits(P, o, m, max_vecs)
    n1 := Rank(P);
    GP := GramMatrix(P);
    Pd := Dual(P : Rescale := false);
    Bd := ChangeRing(BasisMatrix(Pd), Q);
    GD := GramMatrix(Pd);
    sc := Determinant(P) * m;               // clears denominators of GD and of o/m
    Pdi := LatticeWithGram(ChangeRing(sc * GD, Z));
    bound := Z ! (sc * o / m);
    proc := ShortVectorsProcess(Pdi, bound, bound);
    S := [];
    while true do
        u := NextVector(proc);
        if IsZero(u) then break; end if;
        if #S mod 4096 eq 0 then check_orth_deadline(); end if;
        if #S ge max_vecs then
            vprintf GenusReps, 1 : "orth: more than %o type vectors, giving up on this descent\n", max_vecs;
            return false, _;
        end if;
        cu := Eltseq(u);
        if GCD([Z | c : c in cu]) ne 1 then continue; end if;      // u primitive in P^#
        uamb := ChangeRing(Vector(cu), Q) * Bd;
        vamb := m * uamb;
        if exists{ c : c in Eltseq(vamb) | not IsIntegral(c) } then continue; end if;
        v := Vector([Z ! c : c in Eltseq(vamb)]);
        if GCD(Eltseq(v)) ne 1 then continue; end if;              // v primitive in P
        Append(~S, normsign(v));
    end while;
    if #S eq 0 then return true, []; end if;
    // BFS orbit decomposition under the (sign-normalized) action of O(P).
    gens := [Matrix(Z, g) : g in Generators(AutomorphismGroupFaster(P))];
    idx := AssociativeArray();
    for i -> v in S do idx[Eltseq(v)] := i; end for;
    seen := [false : i in [1 .. #S]];
    reps := [];
    for i in [1 .. #S] do
        if seen[i] then continue; end if;
        seen[i] := true;
        orb := [S[i]];
        head := 1;
        while head le #orb do
            w := orb[head];
            head +:= 1;
            for g in gens do
                wn := normsign(w * g);
                key := Eltseq(wn);
                error if not IsDefined(idx, key),
                    "orth: type-vector set not closed under O(P) -- this should not happen";
                j := idx[key];
                if not seen[j] then
                    seen[j] := true;
                    Append(~orb, wn);
                end if;
            end for;
        end while;
        Append(~reps, S[i]);
    end for;
    vprintf GenusReps, 2 : "orth: %o type vectors in %o orbits\n", #S, #reps;
    return true, reps;
end function;

// Gram matrix of v^perp in P (a saturated sublattice of corank 1).
function orth_complement_gram(P, v)
    GP := GramMatrix(P);
    K := KernelMatrix(GP * Matrix(Z, Rank(P), 1, Eltseq(v)));
    NG := K * GP * Transpose(K);
    NG := LLLGram(MatrixAlgebra(Z, Nrows(NG)) ! NG);
    return NG;
end function;

// Cheap isometry-invariant bucket key used before exact isometry comparison.
function cheap_lat_key(L)
    th := ThetaSeries(L, 8);
    return Sprint([Coefficient(th, i) : i in [0 .. 8]]);
end function;

// Bucket key for the process-wide genus cache.
function cache_bucket_key(L)
    th := ThetaSeries(L, 4);
    return Sprint(<Rank(L), Determinant(L), IsEven(L), [Coefficient(th, i) : i in [0 .. 4]]>);
end function;

intrinsic SetGenusOrthTimeBudget(seconds::RngIntElt)
{Abort orth enumeration work (with a catchable error) once this many seconds of
wall clock have elapsed.  Used to bound the batch cache pre-warm, which must run
in the parent process.  Clear with ClearGenusOrthTimeBudget.}
    ZZ := Integers();
    ZZ`GenusOrthDeadline := Realtime() + seconds;
end intrinsic;

intrinsic ClearGenusOrthTimeBudget()
{Remove the wall-clock budget set by SetGenusOrthTimeBudget.}
    ZZ := Integers();
    if assigned ZZ`GenusOrthDeadline then
        delete ZZ`GenusOrthDeadline;
    end if;
end intrinsic;

// Cache entries are <query lattice, aut_graph, reps, label, counters, prov>:
// label is the genus's pipeline label ("" when unknown), counters[i] the
// database counter of reps[i] within the genus (label sort order; [] when
// unknown), and prov the provenance <ambient genus label, [<ambient counter,
// vector>]> recorded when the genus was produced by a descent (0 when unknown).
// Label and counters are what make cross-genus references (ambient_lattice)
// resolvable when a child genus descends from this one.

forward cache_lookup;

// The disk tier: FillGenus persists each definite genus it finishes (with its
// label and counters) under orth_cache/, so parallel fill processes -- and, in
// rank-descending orchestration, the children of this rank -- reuse hub genera
// across process boundaries.  Files are named <bucket-hash>.<content-hash>; two
// processes writing the same genus from different representatives produce two
// files, which lookup deduplicates by isometry.  Writes go through a temp file
// and mv, so concurrent writers are safe.
function orth_cache_path(gkey, content)
    // gkey must be GENUS-invariant (rank/det/parity): a class-dependent key (e.g.
    // theta coefficients) would make lookups from a different class of the same
    // genus miss the file entirely.
    bh := IntegerToString(CollapseIntList(StringToBytes(gkey)));
    ch := IntegerToString(CollapseIntList(StringToBytes(content)));
    return "orth_cache/" * bh * "." * ch, bh;
end function;

function orth_cache_genus_key(L)
    return Sprint(<Rank(L), Determinant(L), IsEven(L)>);
end function;

procedure orth_cache_write(L, reps, label, counters)
    content := label * "\n" *
        Join([ (#counters eq #reps select IntegerToString(counters[i]) else "0")
               * "|" * Join([IntegerToString(Z ! c) : c in Eltseq(GramMatrix(reps[i]))], " ")
             : i in [1 .. #reps]], "\n");
    path, _ := orth_cache_path(orth_cache_genus_key(L), content);
    System("mkdir -p orth_cache");
    tmp := path * ".tmp" * IntegerToString(Getpid());
    Write(tmp, content : Overwrite);
    System("mv " * tmp * " " * path);
end procedure;

intrinsic WriteGenusOrthCache(label::MonStgElt, reps::SeqEnum)
{Persist a finished genus (reps in counter order, i.e. label sort order) to the
on-disk orth cache, so other processes can use it as a descent parent and
resolve ambient_lattice counters against it.}
    if #reps eq 0 then return; end if;
    orth_cache_write(reps[1], reps, label, [1 .. #reps]);
end intrinsic;

// Try to load a genus isometric to L from the disk tier; returns success and
// the parsed entry-shaped tuple.
function orth_cache_disk_lookup(L)
    _, bh := orth_cache_path(orth_cache_genus_key(L), "");
    // "|| true": Pipe raises on nonzero exit, and ls fails when no genus has been
    // persisted yet (the directory is only created on first write).
    files := Split(Pipe("ls orth_cache 2>/dev/null || true", ""), "\n");
    n := Rank(L);
    GL := 0;   // Genus(L), computed lazily (only when a candidate file exists)
    for f in files do
        if #f lt #bh + 1 or f[1 .. #bh + 1] ne bh * "." or "tmp" in f then continue; end if;
        lines := Split(Read("orth_cache/" * f), "\n");
        label := lines[1];
        reps := [];
        counters := [];
        okfile := true;
        for i in [2 .. #lines] do
            parts := Split(lines[i], "|");
            ents := [StringToInteger(c) : c in Split(parts[2], " ")];
            if #ents ne n^2 then okfile := false; break; end if;
            Append(~counters, StringToInteger(parts[1]));
            Append(~reps, LatticeWithGram(Matrix(Z, n, ents)));
        end for;
        if not okfile or #reps eq 0 then continue; end if;
        // A cache file represents a GENUS, while L is one class of it (typically not
        // the class stored first), so membership is decided by genus equality, not
        // by isometry with the first representative.
        if GL cmpeq 0 then GL := Genus(L); end if;
        if Genus(reps[1]) eq GL then
            return true, <reps[1], aut_graph(reps[1] : orth_bd := 6), reps, label, counters, 0>;
        end if;
    end for;
    return false, <L, 0, [], "", [], 0>;
end function;

// Insert a complete list of genus representatives into the process-wide cache,
// keyed by the query lattice and (for small genera) by every class rep, so a
// later chain gluing into any class of this genus hits.  With guard, skips
// insertion when the genus is already cached.
procedure cache_insert(L, reps : guard := false, label := "", counters := [], prov := 0)
    ZZ := Integers();
    if not assigned ZZ`GenusOrthCache then
        ZZ`GenusOrthCache := AssociativeArray();
    end if;
    if guard and cache_lookup(L) then return; end if;
    keyed := [<cache_bucket_key(L), L>];
    if #reps le 50 then
        keyed cat:= [<cache_bucket_key(R), R> : R in reps | R ne L];
    end if;
    for kr in keyed do
        entry := <kr[2], aut_graph(kr[2] : orth_bd := 6), reps, label, counters, prov>;
        if IsDefined(ZZ`GenusOrthCache, kr[1]) then
            lst := (ZZ`GenusOrthCache)[kr[1]];
            Append(~lst, entry);
            (ZZ`GenusOrthCache)[kr[1]] := lst;
        else
            (ZZ`GenusOrthCache)[kr[1]] := [* entry *];
        end if;
    end for;
end procedure;

// Exhaust the remaining mass of the genus of L with the p-neighbour walk,
// seeded with the representatives already found.
function kneser_topup(L, reps0, mass_left)
    vprintf GenusReps, 1 : "orth: topping up mass %o with p-neighbours\n", mass_left;
    p := 2;
    if IsEven(Determinant(L)) and not IsEven(L) then p := 3; end if;
    silent := GetVerbose("GenusReps") lt 2;
    dset := {@ x : x in reps0 @};
    repeat
        dset, _, _, mass_left := list_genus(L : p := p, done := dset,
                                            mass_left := mass_left, silent := silent);
        p := NextPrime(p);
        error if p gt 100, "orth: failed to exhaust the genus by p-neighbours";
    until mass_left eq 0;
    return [x : x in dset];
end function;

// Descend from representatives of the parent genus via O(P)-orbits of
// type-(o, m) vectors, binning the complements by genus.  targets is a
// sequence of tuples <L, Genus(L), Mass(L)>; complements in other genera of
// the same determinant are collected too (they are free), in bins keyed by
// their order of first appearance.  Returns success, a parallel list of
// <found, mass_left> per target, and a list of <genus, lattices, complete>
// for the non-target genera encountered (complete iff the bin's mass tally
// equals its genus mass, in which case it is a full set of representatives).
// Deduplication buckets by a cheap theta key and decides isometry with
// aut_graph/isom_graph (canonical forms can hang on lattices with large
// Weyl-type automorphism groups, e.g. glued D_n classes).
function orth_descend_multi(parents, o, m, targets, max_vecs)
    nt := #targets;
    Dt := Determinant(targets[1][1]);
    found := [* *];
    for i in [1 .. nt] do Append(~found, []); end for;
    tallies := [Rationals() | 0 : i in [1 .. nt]];
    done := [false : i in [1 .. nt]];
    // per-target dedup state
    graphs := [* [* *] : i in [1 .. nt] *];
    buckets := [* AssociativeArray() : i in [1 .. nt] *];
    // bins for genera nobody asked for: parallel lists
    other_genera := [* *];
    other_lats := [* *];
    npairs := 0;
    provs := [* *];
    for i in [1 .. nt] do Append(~provs, []); end for;
    pidx := 0;
    for P in parents do
        pidx +:= 1;
        check_orth_deadline();
        if forall{ i : i in [1 .. nt] | done[i] } then break; end if;
        ok, vreps := typed_vector_orbits(P, o, m, max_vecs);
        if not ok then
            return false, [<found[i], targets[i][3] - tallies[i], provs[i]> : i in [1 .. nt]], [* *];
        end if;
        for v in vreps do
            npairs +:= 1;
            NG := orth_complement_gram(P, v);
            if Determinant(NG) ne Dt then continue; end if;
            NL := LatticeWithGram(NG);
            GN := Genus(NL);
            ti := 0;
            for i in [1 .. nt] do
                if GN eq targets[i][2] then ti := i; break; end if;
            end for;
            if ti eq 0 then
                // not one of the requested genera: bin it (dedup lazily at the end)
                oi := 0;
                for i in [1 .. #other_genera] do
                    if GN eq other_genera[i] then oi := i; break; end if;
                end for;
                if oi eq 0 then
                    Append(~other_genera, GN);
                    Append(~other_lats, [* NL *]);
                else
                    lst := other_lats[oi];
                    Append(~lst, NL);
                    other_lats[oi] := lst;
                end if;
                continue;
            end if;
            if done[ti] then continue; end if;
            key := cheap_lat_key(NL);
            newlat := true;
            tb := buckets[ti];
            tg := graphs[ti];
            if IsDefined(tb, key) then
                for j in tb[key] do
                    if tg[j] cmpeq 0 then
                        tg[j] := aut_graph(found[ti][j] : orth_bd := 6);
                        graphs[ti] := tg;
                    end if;
                    if isom_graph(tg[j], NL) then newlat := false; break; end if;
                end for;
                if newlat then
                    fl := found[ti]; Append(~fl, NL); found[ti] := fl;
                    tgn := graphs[ti]; Append(~tgn, 0); graphs[ti] := tgn;
                    tb[key] cat:= [#found[ti]];
                    buckets[ti] := tb;
                end if;
            else
                fl := found[ti]; Append(~fl, NL); found[ti] := fl;
                tgn := graphs[ti]; Append(~tgn, 0); graphs[ti] := tgn;
                tb[key] := [#found[ti]];
                buckets[ti] := tb;
            end if;
            if newlat then
                pl := provs[ti]; Append(~pl, <pidx, Eltseq(v)>); provs[ti] := pl;
                tallies[ti] +:= 1 / #AutomorphismGroupFaster(NL);
                error if tallies[ti] gt targets[ti][3],
                    "orth: mass exceeded during descent (isometry-class deduplication failed?)";
                if tallies[ti] eq targets[ti][3] then
                    done[ti] := true;
                    vprintf GenusReps, 1 : "orth: descent target %o complete with %o classes (%o pairs so far)\n",
                        ti, #found[ti], npairs;
                end if;
            end if;
        end for;
    end for;
    for i in [1 .. nt] do
        if not done[i] then
            vprintf GenusReps, 1 : "orth: descent target %o: %o classes, mass shortfall %o\n",
                i, #found[i], targets[i][3] - tallies[i];
        end if;
    end for;
    // Deduplicate and mass-check the free bins; only complete genera are returned
    // as complete (safe to cache).
    others := [* *];
    for i in [1 .. #other_genera] do
        lst := other_lats[i];
        if #lst gt 50 then continue; end if;   // cap the free work
        reps := [];
        rgraphs := [* *];
        for NLx in lst do
            new := true;
            for j in [1 .. #reps] do
                if rgraphs[j] cmpeq 0 then rgraphs[j] := aut_graph(reps[j] : orth_bd := 6); end if;
                if isom_graph(rgraphs[j], NLx) then new := false; break; end if;
            end for;
            if new then Append(~reps, NLx); Append(~rgraphs, 0); end if;
        end for;
        gmass := Mass(reps[1]);
        tally := &+[ Rationals() | 1 / #AutomorphismGroupFaster(R) : R in reps ];
        Append(~others, <other_genera[i], reps, tally eq gmass>);
    end for;
    return true, [<found[i], targets[i][3] - tallies[i], provs[i]> : i in [1 .. nt]], others;
end function;

// Single-target wrapper (the recursion uses this).  Returns success, the
// representatives, the missing mass, and the per-class provenance <parent
// index, vector> list (parallel to the representatives).
function orth_descend(parents, o, m, Ltarget, max_vecs, use_cache)
    ok, res, others := orth_descend_multi(parents, o, m,
        [<Ltarget, Genus(Ltarget), Mass(Ltarget)>], max_vecs);
    // Cache any complete free genera picked up along the way.
    if ok and use_cache then
        for oth in others do
            if oth[3] and #oth[2] gt 0 then
                cache_insert(oth[2][1], oth[2] : guard := true);
            end if;
        end for;
    end if;
    return ok, res[1][1], res[1][2], res[1][3];
end function;

// Process-wide cache lookup with disk fallback.  The in-memory cache maps a
// cheap invariant bucket to entries (see above); membership within a bucket is
// decided by exact isometry (isom_graph), never by canonical forms (see
// orth_descend).  On a memory miss the disk tier is consulted and hits are
// promoted, so genera filled by other processes are reused transparently.
function cache_lookup_full(L)
    ZZ := Integers();
    if not assigned ZZ`GenusOrthCache then
        ZZ`GenusOrthCache := AssociativeArray();
    end if;
    bkey := cache_bucket_key(L);
    if IsDefined(ZZ`GenusOrthCache, bkey) then
        for ent in (ZZ`GenusOrthCache)[bkey] do
            if isom_graph(ent[2], L) then return true, ent; end if;
        end for;
    end if;
    ok, ent := orth_cache_disk_lookup(L);
    if ok then
        cache_insert(ent[1], ent[3] : label := ent[4], counters := ent[5], prov := ent[6]);
        return true, ent;
    end if;
    return false, ent;
end function;

function cache_lookup(L)
    ok, ent := cache_lookup_full(L);
    if ok then return true, ent[3]; end if;
    return false, _;
end function;

// Scored glue candidates for L: up to four constructible parents, preferring
// (1) cached parents (free), then (2) even parents -- the p-neighbour base case
// is fast on even lattices at any level, whereas odd unimodular genera at
// rank >= 17 can force the walk to escalate to 3-neighbours -- then (3) the
// strongest determinant reduction, then (4) the cheapest descent (o/m).  The
// guards are purely structural (rank_cap, dual_norm_bound): filtering on the
// parent's mass instead was tried and causes "determinant creep", since Kneser
// cost grows like 2^rank, which the mass does not see.  Tuples are
// <cached, odd, detP, o/m, o, m, P>.
function glue_scored(L, dual_norm_bound, rank_cap, use_cache)
    cands := [c : c in glue_candidates(L) | c[2] le dual_norm_bound];
    scored := [];
    for c in cands do
        okP, P := glue_parent(L, c[3], c[4], c[5]);
        if not okP then continue; end if;
        incache := use_cache and cache_lookup(P);
        vprintf GenusReps, 2 : "orth: candidate parent rank %o det %o (o=%o, m=%o)%o\n",
            Rank(P), Determinant(P), c[3], c[4], incache select " [cached]" else "";
        // Cold parents above rank_cap are blocked (their eventual neighbour base
        // case would be hopeless), cold ODD parents already above rank 16;
        // cached parents of any rank are free.
        if not incache then
            if Rank(P) gt rank_cap then continue; end if;
            if (not IsEven(P)) and Rank(P) gt 16 then continue; end if;
        end if;
        Append(~scored, <incache select 0 else 1, IsEven(P) select 0 else 1,
                         c[1], c[2], c[3], c[4], P>);
        // Construct a wider window than we keep: the static candidate order cannot
        // see cached parents, so stopping at the first four constructible candidates
        // could discard a cached (or better-scored) parent appearing later.
        if #scored ge 12 then break; end if;
    end for;
    cmp := function(a, b)
        for i in [1 .. 4] do
            if a[i] ne b[i] then return a[i] lt b[i] select -1 else 1; end if;
        end for;
        return 0;
    end function;
    Sort(~scored, cmp);
    if #scored gt 4 then scored := scored[1 .. 4]; end if;
    return scored;
end function;

// Enumerate the genus of the positive definite lattice L, preferring the orbit
// method and falling back to (and topping up with) the p-neighbour walk.
// Returns the representatives and the provenance <ambient genus label,
// [<ambient counter, vector>]> (0 when unknown, e.g. p-neighbour fallback or a
// parent with no recorded label/counters).  The recursion itself only consumes
// the first return value; the provenance is what FillGenus stores as
// ambient_genus / ambient_lattice / orthogonal_complement.
function genus_reps_orth_rec(L, depth, dual_norm_bound, max_vecs, rank_cap, use_cache)
    check_orth_deadline();
    if use_cache then
        hit, ent := cache_lookup_full(L);
        if hit then
            vprintf GenusReps, 1 : "orth: cache hit for rank %o det %o\n", Rank(L), Determinant(L);
            return ent[3], ent[6];
        end if;
    end if;

    mass := Mass(L);
    aut := AutomorphismGroupFaster(L);
    reps := [];
    prov := 0;
    if mass eq 1/#aut then
        // The genus consists of a single class.
        reps := [L];
    else
        done := false;
        if depth gt 0 then
            // A descent can fail only by exceeding max_vecs, in which case the
            // next candidate is tried.
            for cand in glue_scored(L, dual_norm_bound, rank_cap, use_cache) do
                besto := cand[5]; bestm := cand[6]; bestP := cand[7];
                vprintf GenusReps, 1 : "orth: rank %o det %o -> parent rank %o det %o (o=%o, m=%o)\n",
                    Rank(L), Determinant(L), Rank(bestP), Determinant(bestP), besto, bestm;
                tpar := Cputime();
                parent_reps := genus_reps_orth_rec(bestP, depth - 1, dual_norm_bound, max_vecs, rank_cap, use_cache);
                vprintf GenusReps, 1 : "orth: parent genus rank %o det %o: %o classes in %o s\n",
                    Rank(bestP), Determinant(bestP), #parent_reps, Cputime(tpar);
                tdesc := Cputime();
                okd, reps, mass_left, dprov := orth_descend(parent_reps, besto, bestm, L, max_vecs, use_cache);
                vprintf GenusReps, 1 : "orth: descent to rank %o det %o took %o s\n",
                    Rank(L), Determinant(L), Cputime(tdesc);
                if not okd then continue; end if;    // too many type vectors: next candidate
                if mass_left ne 0 then
                    reps := kneser_topup(L, reps, mass_left);
                end if;
                // Resolve the ambient references: the parent's label and per-class
                // counters live on its cache entry (present when the parent came
                // from the disk tier, i.e. was already filled by the pipeline).
                plabel := "";
                pcounters := [];
                if use_cache then
                    okp, pent := cache_lookup_full(bestP);
                    if okp then plabel := pent[4]; pcounters := pent[5]; end if;
                end if;
                if plabel ne "" and #pcounters eq #parent_reps then
                    entries := [<pcounters[pr[1]], pr[2]> : pr in dprov];
                    // classes appended by the top-up have no descent provenance
                    entries cat:= [<0, [Z | ]> : i in [#entries + 1 .. #reps]];
                    prov := <plabel, entries>;
                end if;
                done := true;
                break;
            end for;
        end if;
        if not done then
            vprintf GenusReps, 1 : "orth: falling back to p-neighbour enumeration for rank %o det %o\n",
                Rank(L), Determinant(L);
            dset := GenusRepresentativesFaster(L);
            reps := [x : x in dset];
        end if;
    end if;
    if use_cache then
        // Insert the query rep and (for small genera) every class rep, so a later
        // chain that glues into a different class of this genus still hits.
        cache_insert(L, reps : prov := prov);
    end if;
    return reps, prov;
end function;

intrinsic GenusRepresentativesOrth(L::Lat : DualNormBound := 32, MaxVectors := 100000,
                                            Depth := 16, RankCap := 24, UseCache := true) -> SeqEnum, Any
{Representatives for the isometry classes in the genus of the positive definite
lattice L, computed by the orbit method of Chenevier-Taibi: recursively glue L up to
a parent genus of rank+1 and smaller determinant, and descend back by taking orthogonal
complements of O(P)-orbit representatives of vectors of a fixed norm and modulus.
Every level is verified against the Minkowski-Siegel mass; missing classes (which can
happen only in exotic 2-adic situations) are topped up by the p-neighbour walk, and
lattices for which no cheap descent exists fall back to it entirely.

DualNormBound caps the dual norm o/m of the vectors enumerated in a descent step;
MaxVectors caps the number of type vectors considered per parent representative;
Depth caps the recursion; RankCap is the largest rank for which an uncached parent
genus may be enumerated (cached parents are always used); UseCache memoizes
enumerated genera in the current Magma process (with a disk tier under orth_cache/).

The second return value is the provenance <ambient genus label, [<ambient counter,
vector>]> when the final descent's parent genus carries a label and counters (i.e.
was already filled by the pipeline), and 0 otherwise.}
    G := GramMatrix(L);
    require IsPositiveDefinite(G) : "The lattice must be positive definite";
    require forall{ e : e in Eltseq(G) | IsIntegral(e) } : "The lattice must be integral";
    L0 := LatticeWithGram(MatrixAlgebra(Z, Rank(L)) ! G);
    // Reduce by the content: the genus of c*L' is the c-scaling of the genus of L'.
    c := GCD(Eltseq(GramMatrix(L0)));
    if c gt 1 then
        prim := LatticeWithGram(GramMatrix(L0) div c);
        reps := genus_reps_orth_rec(prim, Depth, DualNormBound, MaxVectors, RankCap, UseCache);
        // provenance refers to the primitive lattices, so do not propagate it
        return [LatticeWithGram(c * GramMatrix(r)) : r in reps], 0;
    end if;
    reps, prov := genus_reps_orth_rec(L0, Depth, DualNormBound, MaxVectors, RankCap, UseCache);
    return reps, prov;
end intrinsic;

function genus_reps_orth(L)
    reps, prov := GenusRepresentativesOrth(L);
    return reps, prov;
end function;

intrinsic GenusRepresentativesOrthBatch(targets::SeqEnum[Lat] : DualNormBound := 32,
        MaxVectors := 100000, Depth := 16, RankCap := 24, UseCache := true) -> List
{Enumerate the genera of a batch of positive definite integral lattices together by
the orbit method.  Targets whose best glue lands in the same parent genus with the
same vector type share a single descent sweep, during which every orthogonal
complement is binned by genus -- so one pass over the parent's vector orbits can
emit representatives for several requested genera at once, and complete
representative lists for genera nobody asked for are cached for later use.
Each target is mass-verified independently (with p-neighbour top-up on any
shortfall), so the output is provably complete per genus.  Returns a list of
representative sequences parallel to the input.}
    nt := #targets;
    results := [* *];
    for i in [1 .. nt] do Append(~results, 0); end for;
    pending := [];
    for i in [1 .. nt] do
        L := targets[i];
        G := GramMatrix(L);
        require IsPositiveDefinite(G) : "All lattices must be positive definite";
        require forall{ e : e in Eltseq(G) | IsIntegral(e) } : "All lattices must be integral";
        L0 := LatticeWithGram(MatrixAlgebra(Z, Rank(L)) ! G);
        if GCD(Eltseq(GramMatrix(L0))) gt 1 then
            // scaled lattices take the single-target path (rare)
            results[i] := GenusRepresentativesOrth(L : DualNormBound := DualNormBound,
                MaxVectors := MaxVectors, Depth := Depth, RankCap := RankCap, UseCache := UseCache);
            continue;
        end if;
        if UseCache then
            hit, cached := cache_lookup(L0);
            if hit then results[i] := cached; continue; end if;
        end if;
        Append(~pending, <i, L0>);
    end for;
    // Group the remaining targets by (o, m, parent genus) of their best candidate.
    groups := [* *];   // entries <o, m, P, GP, [<index, L0, GenusL0, MassL0>]>
    singles := [];
    if Depth le 0 then
        // No glue level available: the single-target path's depth accounting
        // handles the fallback (the audit found the batch path allowed one more
        // level than Depth permits).
        singles := pending;
        pending := [];
    end if;
    for t in pending do
        scored := glue_scored(t[2], DualNormBound, RankCap, UseCache);
        if #scored eq 0 then
            Append(~singles, t);
            continue;
        end if;
        c := scored[1];
        o := c[5]; m := c[6]; P := c[7];
        entry := <t[1], t[2], Genus(t[2]), Mass(t[2])>;
        gi := 0;
        for j in [1 .. #groups] do
            if groups[j][1] eq o and groups[j][2] eq m and groups[j][4] eq Genus(P) then
                gi := j; break;
            end if;
        end for;
        if gi eq 0 then
            Append(~groups, <o, m, P, Genus(P), [entry]>);
        else
            g := groups[gi];
            lst := g[5];
            Append(~lst, entry);
            groups[gi] := <g[1], g[2], g[3], g[4], lst>;
        end if;
    end for;
    for g in groups do
        o := g[1]; m := g[2]; P := g[3]; members := g[5];
        vprintf GenusReps, 1 : "orth batch: %o target(s) share parent rank %o det %o (o=%o, m=%o)\n",
            #members, Rank(P), Determinant(P), o, m;
        // The glue step above already consumed one level, matching the
        // single-target recursion's Depth accounting.
        parent_reps := genus_reps_orth_rec(P, Depth - 1, DualNormBound, MaxVectors, RankCap, UseCache);
        tlist := [<mem[2], mem[3], mem[4]> : mem in members];
        okd, res, others := orth_descend_multi(parent_reps, o, m, tlist, MaxVectors);
        if not okd then
            // vector explosion: fall back to the single-target path per member,
            // which tries the remaining candidates internally
            for mem in members do Append(~singles, <mem[1], mem[2]>); end for;
            continue;
        end if;
        plabel := "";
        pcounters := [];
        if UseCache then
            okp, pent := cache_lookup_full(P);
            if okp then plabel := pent[4]; pcounters := pent[5]; end if;
        end if;
        for k in [1 .. #members] do
            reps := res[k][1];
            if res[k][2] ne 0 then
                reps := kneser_topup(members[k][2], reps, res[k][2]);
            end if;
            prov := 0;
            if plabel ne "" and #pcounters eq #parent_reps then
                entries := [<pcounters[pr[1]], pr[2]> : pr in res[k][3]];
                entries cat:= [<0, [Z | ]> : i in [#entries + 1 .. #reps]];
                prov := <plabel, entries>;
            end if;
            if UseCache then cache_insert(members[k][2], reps : prov := prov); end if;
            results[members[k][1]] := reps;
        end for;
        if UseCache then
            for oth in others do
                if oth[3] and #oth[2] gt 0 then
                    cache_insert(oth[2][1], oth[2] : guard := true);
                end if;
            end for;
        end if;
    end for;
    for t in singles do
        results[t[1]] := genus_reps_orth_rec(t[2], Depth, DualNormBound, MaxVectors, RankCap, UseCache);
    end for;
    return results;
end intrinsic;

/****************************************************************************
 * Method dispatch (moved from fill_genus.m)
 ****************************************************************************/

intrinsic UseOrthHeuristic(n::RngIntElt, mass::Any) -> BoolElt
{Whether the orbit method should be preferred over the p-neighbour walk for a
definite genus of rank n and Minkowski-Siegel mass mass (pass 0 when unknown).

Tuned on benchmarks (2026-07): the orbit method wins from rank 15 (Kneser cost
grows like 2^rank), and at any rank once the genus is genuinely large (27x
faster at rank 6, det 2^11, h = 815).  NB the mass threshold keys on
DETERMINANT scale, not class-number scale: mass ~ h/|Aut| and automorphism
groups at low determinant are enormous (measured on the C=256 census: rank-12
genera with h = 55 have mass ~ 1.4e-4), so small-determinant genera never
trip it -- which matches the benchmarks, where the neighbour walk wins that
regime anyway.  Both measured crossover regimes (high rank at small
determinant, large smooth determinant at low rank with h in the hundreds)
satisfy mass >= 5; be wary of reading the threshold as "class number >= 5".}
    if n ge 15 then return true; end if;
    if Type(mass) in [RngIntElt, FldRatElt] and mass ge 5 then return true; end if;
    return false;
end intrinsic;

intrinsic GenusReps(L0::Lat : Timeout := 1800, UseOrth := true) -> BoolElt, SeqEnum, Any
{Enumerate the genus of L0 with the method dispatch used by the LMFDB lattice
pipeline, under a wall-clock timeout (in seconds).  Returns success, the list of
representatives (empty on failure), and -- for the orbit method -- the provenance
<ambient genus label, [<ambient counter, vector>]> (0 when unavailable).

Rank-2 lattices of square determinant are handled directly (Magma's
GenusRepresentatives rejects square discriminants); positive definite lattices of
rank >= 3 use the orbit method (UseOrth, with internal p-neighbour fallback/top-up)
or the mass-verified p-neighbour enumeration; negative definite lattices are negated,
enumerated as positive definite, and negated back; genuinely indefinite lattices of
rank >= 3 enumerate every spinor genus in the genus (each is a single class by
Eichler); everything else uses Magma's GenusRepresentatives.

Every TimeoutCall here carries HardKill: Magma's Alarm is an uncatchable SIGALRM
that is never delivered inside uninterruptible C-level builtins (exactly what the
enumerators are), so without the watchdog a stuck child never dies and the parent
blocks forever in WaitForAllChildren -- in production this looked like ~112-minute
fill batches that still reported failed=0.  Do not drop it in a refactor.}
    n := Rank(L0);
    Gm := GramMatrix(L0);
    posdef := IsPositiveDefinite(Gm);
    negdef := (not posdef) and IsPositiveDefinite(-Gm);
    prov := 0;
    // Rank-2 lattices of square determinant -m^2 have isotropic (split) forms, on
    // which Magma's GenusRepresentatives fails (it reduces to binary quadratic
    // forms, which reject square discriminants); these are handled directly.  Note
    // the class number is NOT always 1 -- e.g. for m = 5 some genera have two
    // classes.
    if n eq 2 and IsSquare(-Determinant(L0)) then
        genus_success, reps, elapsed := TimeoutCall(Timeout, genus_reps_square_disc, <L0>, 1 : HardKill := true);
        vprintf GenusReps, 1 : "Genus representatives (square discriminant) computed in %o seconds\n", elapsed;
    elif posdef and (n ge 3) then
        // Positive definite rank >= 3: the orbit method (GenusRepresentativesOrth)
        // descends from easier parent genera and mass-verifies, falling back
        // internally to the p-neighbour enumeration (GenusRepresentativesFaster)
        // whenever no cheap descent exists; both are provably complete when they
        // terminate.  Magma's GenusRepresentatives is avoided entirely (slow, and
        // reliably fails at rank >= 7 with the "cs must be non-empty" bug).
        if UseOrth then
            genus_success, reps, elapsed := TimeoutCall(Timeout, genus_reps_orth, <L0>, 2 : HardKill := true);
        else
            genus_success, reps, elapsed := TimeoutCall(Timeout, genus_reps_Faster, <L0>, 1 : HardKill := true);
        end if;
        vprintf GenusReps, 1 : "Genus representatives (%o) computed in %o seconds\n",
            UseOrth select "orbit method" else "p-neighbours", elapsed;
    elif negdef then
        // Negative definite: negate, enumerate the full positive definite genus,
        // negate back.  (Previously these fell into the indefinite branch, whose
        // Eichler justification does not apply to definite lattices.)
        fn := UseOrth select genus_reps_negdef_orth else genus_reps_negdef_walk;
        genus_success, reps, elapsed := TimeoutCall(Timeout, fn, <L0>, 1 : HardKill := true);
        vprintf GenusReps, 1 : "Genus representatives (negated definite) computed in %o seconds\n", elapsed;
    elif (n ge 3) then
        // Genuinely indefinite rank >= 3: one representative per spinor genus of the
        // genus; each spinor genus is a single isometry class by Eichler.
        genus_success, reps, elapsed := TimeoutCall(Timeout, genus_reps_Spinor, <L0>, 1 : HardKill := true);
        vprintf GenusReps, 1 : "Genus representatives (spinor genera) computed in %o seconds\n", elapsed;
    else
        // Rank 1-2: neither the definite rank >= 3 methods nor the spinor route
        // applies, so use Magma's general GenusRepresentatives.
        genus_success, reps, elapsed := TimeoutCall(Timeout, genus_reps_Magma, <L0>, 1 : HardKill := true);
        vprintf GenusReps, 1 : "Genus representatives computed in %o seconds\n", elapsed;
    end if;
    if genus_success then
        if posdef and (n ge 3) and UseOrth then
            prov := reps[2];
        end if;
        return true, reps[1], prov;
    end if;
    return false, [], 0;
end intrinsic;
