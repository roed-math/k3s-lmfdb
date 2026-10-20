// Universality, even-universality and (rank-1) regularity helpers.
//
// Everything these intrinsics report as decided is a theorem: positive claims come
// from published criterion sets (plus the Bhargava--Hanke 290 theorem, whose proof
// is a 2011 preprint that remains unpublished -- flagged in the knowl), negative
// claims from explicit counterexamples (a truant value, or a witness lattice that
// provably does not embed).  When neither is available within budget, the caller
// stores "\N".  Layer architecture and the post-Oscar db-relative refinement:
// oscar/INTEGRATION_PLAN.md section 10.
//
// Criterion sets (positive definite integer-matrix forms = our lattices):
//  * S1 = {1,2,3,5,6,7,10,14,15}: Conway--Schneeberger 15-theorem (proof Bhargava,
//    Contemp. Math. 272 (2000)).  Represents S1 <=> represents all positive integers.
//  * S290 = {1,...,290} criticals (Bhargava--Hanke, for integer-VALUED forms).
//    Halving an even lattice gives an integer-valued form, so an even lattice
//    represents all even integers iff it represents 2*S290.  Every even-norm vector
//    of any lattice L lies in the even sublattice of L, so testing 2*S290 on L
//    itself is correct for odd L too.
//  * S2 = {<1,1>, <2,3>, <3,3>, [2,1,2], [2,1,3], [2,1,4]} where [a,b,c] has Gram
//    [[a,b],[b,c]]: Kim--Kim--Oh (Contemp. Math. 249 (1999)); unique minimal set
//    (Kominers, arXiv:0708.4340).  Represents S2 <=> 2-universal.
//  * S8 = {I8, E8}: B.-K. Oh (Proc. AMS 128 (1999)); unique minimal set (Kominers,
//    arXiv:0807.2099).  Represents S8 <=> 8-universal.
//  * No criterion sets are known for n in 3..7, nor for even n-universality with
//    n >= 2 (Ji--Kim--Oh, J. Korean Math. Soc. 58 (2021), still have 20 unproven
//    rank-5 even-2-universal candidates).  There, an exact value needs the proven
//    lower bound to meet a witness upper bound; otherwise we stay undecided.
//  * Rank facts used as shortcuts: n-universal requires rank >= n+3 for n <= 5
//    (Oh, Proc. AMS 128); Ko (1937): I_{n+3} is n-universal for n <= 5; norm-2
//    vectors of Z^m span root lattices of type A/D only, so E6 never embeds in any
//    I_m and universality(I_m) caps at 5.

function S1Set()
    return [1, 2, 3, 5, 6, 7, 10, 14, 15];
end function;

function S290Set()
    return [1, 2, 3, 5, 6, 7, 10, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31,
            34, 35, 37, 42, 58, 93, 110, 145, 203, 290];
end function;

// Kim--Kim--Oh 2-criterion set; SymmetricMatrix([a,b,c]) is [[a,b],[b,c]]
function S2Grams()
    return [SymmetricMatrix(t) : t in
            [[1, 0, 1], [2, 0, 3], [3, 0, 3], [2, 1, 2], [2, 1, 3], [2, 1, 4]]];
end function;

// Even binary lattices used to certify even_universality <= 1, cheap frequent
// failures first ([2,1,4] already fails for D4: (v,w)=1 with v a D4 root forces an
// odd coordinate sum).
function EvenBinaryWitnessGrams()
    return [SymmetricMatrix(t) : t in
            [[2, 1, 4], [2, 0, 2], [2, 1, 2], [2, 0, 4], [2, 0, 6], [4, 0, 4],
             [2, 0, 14], [4, 1, 4], [2, 0, 30]]];
end function;

// Rank-k lattices used to certify universality <= k-1.  <1,...,1,7> works against
// unimodular-complement lattices: any I_{k-1} embedded in I_m has orthogonal
// complement I_{m-k+1}, and 7 is not a sum of three squares.  (The provable
// 8-universal lower bound is checked directly against S8 = {I8, E8} in
// LatticeUniversalityData, independent of the rest of this pool.)
function WitnessGrams(k)
    pool := [MatrixRing(Integers(), k) ! 1];                                   // I_k
    if k ge 2 then
        Append(~pool, DiagonalMatrix(Integers(), [1 : i in [1..k-1]] cat [7]));
        Append(~pool, DiagonalMatrix(Integers(), [1 : i in [1..k-1]] cat [15]));
    end if;
    Append(~pool, GramMatrix(Lattice("A", k)));
    if k ge 4 then
        Append(~pool, GramMatrix(Lattice("D", k)));
    end if;
    if k in [6, 7, 8] then
        Append(~pool, GramMatrix(Lattice("E", k)));
    end if;
    Append(~pool, DiagonalMatrix(Integers(), [2 : i in [1..k]]));              // A1^k
    return pool;
end function;

intrinsic FirstNonRepresentedValue(L::Lat, values::SeqEnum) -> RngIntElt, RngIntElt
{Scan values for representability by L as vector norms.  Returns <first provably
non-represented value or 0, first undecided value or 0>.  A thin shell [m,m] costs
a Fincke--Pohst traversal of the whole norm-<=m tree -- a single C call that
Magma's Alarm (hence TimeoutCall's child-side deadline) CANNOT interrupt, and
hopeless above rank ~5 anyway -- so every enumeration here is bounded by
construction: positive certificates come from the rank-4 and rank-5 sublattices
spanned by the shortest LLL basis vectors (trees of ~m^2..m^2.5 nodes), and a
conclusive full-L negative is attempted only when rank(L) <= 5.  Values that a
large-rank lattice's small sublattices miss are reported undecided, never guessed.}
    subs := [];
    if Rank(L) gt 5 then
        // gram of the sublattice spanned by the k shortest LLL basis vectors.  LLL
        // does not sort by norm: on E8 perp I8 the leading block is all E8 roots,
        // an even sublattice that misses every odd value.
        GL := GramMatrix(LLL(L));
        idx := [1..Rank(L)];
        Sort(~idx, func<a, b | GL[a][a] - GL[b][b]>);
        for k in [4, 5] do
            Append(~subs, LatticeWithGram(Matrix(Integers(), k, k,
                [GL[idx[i]][idx[j]] : j in [1..k], i in [1..k]]) : CheckPositive := false));
        end for;
    end if;
    for m in values do
        if exists{K : K in subs | #ShortVectors(K, m, m : Max := 1) gt 0} then
            continue;
        end if;
        if Rank(L) le 5 then
            if #ShortVectors(L, m, m : Max := 1) eq 0 then
                return m, 0;
            end if;
        else
            return 0, m;
        end if;
    end for;
    return 0, 0;
end intrinsic;

// Scan `values` for representability by L.  theta covers norms up to prec exactly;
// the values beyond that go through a single budgeted fork (TimeoutCall) rather
// than one per value.  Returns <decided, all_represented>.
function represents_all_values(L, values, theta, prec, tail_seconds)
    values := Sort(values);
    for m in values do
        if m le prec and theta[m+1] eq 0 then
            return true, false;
        end if;
    end for;
    tail := [Integers() | m : m in values | m gt prec];
    if #tail eq 0 then
        return true, true;
    end if;
    success, res, _ := TimeoutCall(tail_seconds * #tail, FirstNonRepresentedValue, <L, tail>, 2 : HardKill := true);
    if not success then
        return false, false;                               // budget exhausted: undecided
    end if;
    missing := res[1];
    undecided := res[2];
    if missing gt 0 then
        return true, false;
    elif undecided gt 0 then
        return false, false;
    end if;
    return true, true;
end function;

// Shared core for the boolean flags.  theta may be [] (prec then -1).
// Returns u_dec, u_val, eu_dec, eu_val.
function flags_internal(L, theta, tail_seconds)
    prec := #theta - 1;
    if IsEven(L) then
        u_dec := true;  u_val := false;                    // an even lattice misses 1
    else
        u_dec, u_val := represents_all_values(L, S1Set(), theta, prec, tail_seconds);
    end if;
    if u_dec and u_val then
        eu_dec := true;  eu_val := true;                   // universal => even-universal
    else
        eu_dec, eu_val := represents_all_values(L, [2*m : m in S290Set()], theta, prec, tail_seconds);
    end if;
    return u_dec, u_val, eu_dec, eu_val;
end function;

// L isometric to I_n?  n mutually orthogonal norm-1 vectors span a unimodular
// sublattice, which therefore splits L; determinant 1 forces equality.
function is_standard_In(L)
    if IsEven(L) or Determinant(L) ne 1 or Minimum(L) ne 1 then
        return false;
    end if;
    return #ShortVectors(L, 1, 1) eq Rank(L);              // up to sign: one per axis
end function;

// Orbit representatives (up to sign) for a matrix group acting on a list of
// vectors that is closed under the action up to sign (a ShortVectors shell).
// Sound to use at the top search level: composing an embedding with an
// automorphism of L (or a global sign flip) moves v_1 onto a representative.
function orbit_reps(A, vecs)
    reps := [];
    seen := {};
    for v in vecs do
        if v in seen then
            continue;
        end if;
        Append(~reps, v);
        orb := Orbit(A, v);
        seen join:= {x : x in orb} join {-x : x in orb};
    end for;
    return reps;
end function;

forward embed_search;

// Depth-first exact search for vectors v_1..v_n in L with prescribed Gram matrix
// MG.  Level j+1 needs (v_i, y) = MG[i][j+1] for i <= j -- an integral affine
// solve -- and (y, y) = MG[j+1][j+1] -- close-vector enumeration on the solve's
// kernel lattice.  Candidates stream lazily (ShortVectorsProcess /
// CloseVectorsProcess), so successful searches return without materializing large
// shells, while running a stream dry is an exhaustion proof.  The zero vector is
// the end-of-process sentinel; the one candidate it could genuinely encode
// (y = sol, when |sol|^2 = t) is checked out-of-band, and a mid-stream zero is
// skipped once in that case.  `nodes` counts candidates tried against maxnodes.
// Returns <status, nodes>: 1 = embedding found, 0 = exhausted (proof of
// non-embedding), -1 = node cap hit (undecided).
function embed_search(L, G, MG, chosen, nodes, maxnodes, deadline)
    j := #chosen;
    n := Nrows(MG);
    N := Nrows(G);
    t := MG[j+1][j+1];
    cands := [];
    proc := 0;
    sol := 0;
    zero_is_solution := false;
    if j eq 0 then
        // the Aut-collapsed top level lives in RepresentsLattice; this path is the
        // aut-free one, streamed lazily
        proc := ShortVectorsProcess(L, t, t);
    else
        V := Matrix(Integers(), j, N, &cat[Eltseq(v) : v in chosen]);
        T := Transpose(V * G);                             // (v_i, y) = (y * T)_i
        c := Vector(Integers(), [MG[i][j+1] : i in [1..j]]);
        ok, sol := IsConsistent(T, c);
        if not ok then
            return 0, nodes;                               // inner products unattainable
        end if;
        if (Matrix(sol) * G * Transpose(Matrix(sol)))[1][1] eq t then
            Append(~cands, sol);                           // the kernel-zero candidate
            zero_is_solution := true;
        end if;
        K := KernelMatrix(T);
        if Nrows(K) gt 0 then
            KL := LatticeWithBasis(K, G);
            proc := CloseVectorsProcess(KL, ChangeRing(-sol, Rationals()), t, t);
        end if;
    end if;
    zero_skipped := false;
    repeat
        if #cands gt 0 then
            y := cands[#cands];
            Prune(~cands);
        elif proc cmpne 0 then
            x := NextVector(proc);
            if IsZero(x) then
                if zero_is_solution and not zero_skipped then
                    zero_skipped := true;                  // the genuine k = 0 solution,
                    continue;                              // already tried via cands
                end if;
                break;                                     // end of process
            end if;
            y := j eq 0 select Vector(Integers(), Eltseq(x))
                        else Vector(Integers(), Eltseq(x)) + sol;
        else
            break;
        end if;
        nodes +:= 1;
        if nodes gt maxnodes or Cputime() gt deadline then
            return -1, nodes;                              // budget hit: undecided
        end if;
        if j + 1 eq n then
            return 1, nodes;
        end if;
        status, nodes := embed_search(L, G, MG, Append(chosen, y), nodes, maxnodes, deadline);
        if status ne 0 then
            return status, nodes;
        end if;
    until false;
    return 0, nodes;
end function;

intrinsic RepresentsLattice(L::Lat, M::Lat : MaxNodes := 10000000, MaxSeconds := 120, Aut := 0) -> BoolElt, BoolElt
{Decide whether the positive definite integral lattice M embeds isometrically in L
(any embedding, not necessarily primitive).  Returns <decided, represents>; both
outcomes are certificates (an explicit embedding / an exhausted exact search).  If
the node or time cap is hit, decided is false and the second value is meaningless.
Aut may be the automorphism group of L (matrices in L's basis coordinates, as
returned by AutomorphismGroup); it collapses the top search level to orbit
representatives, which is what makes exhaustive non-embedding proofs feasible on
vector-rich lattices.}
    require IsPositiveDefinite(GramMatrix(L)) and IsPositiveDefinite(GramMatrix(M)) :
        "universality machinery is for positive definite lattices";
    if Rank(M) gt Rank(L) then
        return true, false;
    end if;
    if IsEven(L) and not IsEven(M) then
        return true, false;
    end if;
    if Minimum(M) lt Minimum(L) then
        return true, false;                                // vectors of M are vectors of L
    end if;
    G := ChangeRing(GramMatrix(L), Integers());
    MG := ChangeRing(GramMatrix(M), Integers());
    // Search in Gram coordinates (identity basis matrix), so that element
    // coordinates and G live in the same frame regardless of L's ambient
    LG := LatticeWithGram(G : CheckPositive := false);
    // LLL-reduce, then order the target basis by increasing norm: small shells at
    // the top of the search tree, and the top shell is the one Aut can collapse
    MG := LLLGram(MG);
    n := Nrows(MG);
    idx := [1..n];
    Sort(~idx, func<a, b | MG[a][a] - MG[b][b]>);
    MG := Matrix(Integers(), n, n, [MG[idx[i]][idx[j]] : j in [1..n], i in [1..n]]);
    deadline := Cputime() + MaxSeconds;
    if Aut cmpne 0 then
        // With the top level materialized as orbit representatives, hunt for an
        // easy success under a small per-representative budget before paying for
        // exhaustion: one representative's dead subtree (e.g. the I8-part roots
        // when searching for E8 inside E8 perp I8) must not starve a sibling
        // whose subtree succeeds immediately.  The probe phase has a global
        // budget of its own so that a large uncollapsed shell cannot multiply it.
        t1 := MG[1][1];
        // un-forked C call: acceptable because the ascending-diagonal ordering puts
        // the SMALLEST norm on top, and minimal shells are the small ones; a target
        // with a large minimal norm would be the case to fork
        shell := [Vector(Integers(), Eltseq(v[1])) : v in ShortVectors(LG, t1, t1)];
        reps := #shell le 20000 select orbit_reps(Aut, shell) else shell;
        if n eq 1 then
            return true, #reps gt 0;
        end if;
        probe_total := MaxNodes div 4;
        probe := Max(200, Min(2000, probe_total div (2 * Max(1, #reps))));
        spent := 0;
        for y in reps do
            if spent ge probe_total or Cputime() gt deadline then
                break;
            end if;
            status, used := embed_search(LG, G, MG, [y], 0, Min(probe, probe_total - spent), deadline);
            spent +:= used;
            if status eq 1 then
                return true, true;
            end if;
        end for;
        nodes := 0;
        for y in reps do
            status, nodes := embed_search(LG, G, MG, [y], nodes, MaxNodes, deadline);
            if status eq 1 then
                return true, true;
            end if;
            if status eq -1 then
                vprintf FillGenus, 1 : "RepresentsLattice: budget hit (%o nodes)\n", nodes;
                return false, false;
            end if;
        end for;
        return true, false;
    end if;
    status, nodes := embed_search(LG, G, MG, [], 0, MaxNodes, deadline);
    if status eq -1 then
        vprintf FillGenus, 1 : "RepresentsLattice: budget hit (%o nodes)\n", nodes;
        return false, false;
    end if;
    return true, status eq 1;
end intrinsic;

intrinsic UniversalityFlags(L::Lat : Theta := [], TailSeconds := 10) -> BoolElt, BoolElt, BoolElt, BoolElt
{The booleans is_universal (15-theorem) and is_even_universal (doubled 290-theorem)
for a positive definite integral lattice, each as a pair <decided, value>.  Theta is
the stored theta coefficient list ([theta_0, ..., theta_prec]); norms beyond its
precision go through one TimeoutCall-guarded first-vector scan, so an exhausted
budget yields undecided rather than a wrong answer.}
    require IsPositiveDefinite(GramMatrix(L)) : "universality machinery is for positive definite lattices";
    u_dec, u_val, eu_dec, eu_val := flags_internal(L, Theta, TailSeconds);
    return u_dec, u_val, eu_dec, eu_val;
end intrinsic;

intrinsic LatticeUniversalityData(L::Lat : Theta := [], TailSeconds := 10, MaxNodes := 10000000,
                                            MaxSeconds := 120, TotalSeconds := 600, Aut := 0)
    -> BoolElt, RngIntElt, BoolElt, BoolElt, RngIntElt, BoolElt
{The universality and even_universality columns (largest n such that every, resp.
every even, positive definite rank-n lattice embeds in L), each as a triple
<defined, value, proved>.  proved means the value is a theorem: the criterion
ladder (S1/S2/S8, plus Ko's I_m values) meets a witness upper bound.  When only
one side is available the value is a heuristic -- the largest k whose entire
witness pool provably embeds -- stored with proved = false, to be refined by the
db-relative closure pass over lat_inclusions.  defined = false means not even a
heuristic is available (budgets exhausted); store null.  Aut is the automorphism
group of L in basis coordinates (connect passes the stored one); if absent it is
computed, budgeted, the first time an embedding search needs it.}
    require IsPositiveDefinite(GramMatrix(L)) : "universality machinery is for positive definite lattices";
    u_dec, u_val, eu_dec, eu_val := flags_internal(L, Theta, TailSeconds);
    ladder_deadline := Cputime() + TotalSeconds;

    // The embedding searches below only run for the rare lattices that survive the
    // flag stage; an automorphism group makes their exhaustive branches feasible.
    A := Aut;
    In_override := u_dec and u_val and is_standard_In(L);
    needs_search := (u_dec and u_val and not In_override and Rank(L) ge 5)
                    or (eu_dec and eu_val);
    if needs_search and A cmpeq 0 then
        okA, resA, _ := TimeoutCall(60, AutomorphismGroupFaster,
                                    <LatticeWithGram(ChangeRing(GramMatrix(L), Integers()) : CheckPositive := false)>, 1
                                    : HardKill := true);
        if okA then
            A := resA[1];
        end if;
    end if;

    Udef := false;  Uval := 0;  Uproved := false;
    if u_dec and not u_val then
        Udef := true;  Uproved := true;                    // not even 1-universal
    elif u_dec and u_val then
        if In_override then
            // Ko: I_{n+3} is n-universal for n <= 5; E6 never embeds in I_m, so 5 caps
            Udef := true;  Uval := Min(Rank(L) - 3, 5);  Uproved := true;
        elif Rank(L) lt 5 then
            Udef := true;  Uval := 1;  Uproved := true;    // 2-universal needs rank >= 5
        else
            lower := 1;
            upper := -1;                                   // -1 = no upper bound yet
            all2 := true;
            unk2 := false;
            for g in S2Grams() do
                if Cputime() gt ladder_deadline then
                    unk2 := true;  break;
                end if;
                dec, rep := RepresentsLattice(L, LatticeWithGram(g) : MaxNodes := MaxNodes,
                                              MaxSeconds := MaxSeconds, Aut := A);
                if dec and not rep then
                    all2 := false;  break;
                elif not dec then
                    unk2 := true;
                end if;
            end for;
            if not all2 then
                Udef := true;  Uval := 1;  Uproved := true;
            elif unk2 then
                // S2 undecidable within budget: universal but 2-universality open
                Udef := true;  Uval := 1;  Uproved := false;
            else
                lower := 2;
                heur := 2;                                 // largest k with a fully-embedded pool
                heur_stuck := false;
                for k in [3..9] do
                    if Cputime() gt ladder_deadline then
                        heur_stuck := true;  break;
                    end if;
                    k_all_decided_true := true;
                    for g in WitnessGrams(k) do
                        if Cputime() gt ladder_deadline then
                            k_all_decided_true := false;  break;
                        end if;
                        dec, rep := RepresentsLattice(L, LatticeWithGram(g) : MaxNodes := MaxNodes,
                                                      MaxSeconds := MaxSeconds, Aut := A);
                        if dec and not rep then
                            upper := k - 1;  break;
                        elif not dec then
                            k_all_decided_true := false;
                        end if;
                    end for;
                    if upper ne -1 then
                        break;
                    end if;
                    if k_all_decided_true and not heur_stuck then
                        heur := k;
                    else
                        heur_stuck := true;                // a gap: stop advancing the heuristic
                    end if;
                    if k eq 8 and upper eq -1 and lower lt 8 then
                        // Oh's criterion needs only S8 = {I8, E8}: check them
                        // directly so an undecided unrelated pool member (A8, D8,
                        // ...) cannot downgrade a provable 8 to a heuristic
                        okI, repI := RepresentsLattice(L, StandardLattice(8) :
                                        MaxNodes := MaxNodes, MaxSeconds := MaxSeconds, Aut := A);
                        okE, repE := RepresentsLattice(L, Lattice("E", 8) :
                                        MaxNodes := MaxNodes, MaxSeconds := MaxSeconds, Aut := A);
                        if okI and repI and okE and repE then
                            lower := 8;
                        end if;
                    end if;
                end for;
                if upper ne -1 then
                    assert upper ge lower;                 // proven-lower witnesses embed
                end if;
                if upper ne -1 and upper eq lower then
                    Udef := true;  Uval := lower;  Uproved := true;
                else
                    // heuristic: every rank-<=heur witness embeds (and heur <= upper
                    // when an upper bound exists, since heur stops before a failure)
                    Udef := true;  Uval := Max(lower, heur);  Uproved := false;
                end if;
            end if;
        end if;
    end if;

    EUdef := false;  EUval := 0;  EUproved := false;
    if eu_dec and not eu_val then
        EUdef := true;  EUproved := true;                  // not even even-1-universal
    elif eu_dec and eu_val then
        // no even 2-criterion set is known (Ji--Kim--Oh 2021), so 1 is the best
        // provable value; if every even binary witness embeds, report the
        // heuristic value 2 for the closure pass to refine
        even_unknown := false;
        even_failed := false;
        for g in EvenBinaryWitnessGrams() do
            if Cputime() gt ladder_deadline then
                even_unknown := true;  break;
            end if;
            dec, rep := RepresentsLattice(L, LatticeWithGram(g) : MaxNodes := MaxNodes,
                                          MaxSeconds := MaxSeconds, Aut := A);
            if dec and not rep then
                even_failed := true;  break;
            elif not dec then
                even_unknown := true;
            end if;
        end for;
        if even_failed then
            EUdef := true;  EUval := 1;  EUproved := true;
        elif even_unknown then
            EUdef := true;  EUval := 1;  EUproved := false;
        else
            EUdef := true;  EUval := 2;  EUproved := false;
        end if;
    end if;

    return Udef, Uval, Uproved, EUdef, EUval, EUproved;
end intrinsic;

intrinsic RankOneIrregularityWitness(thetas::SeqEnum, i::RngIntElt) -> RngIntElt
{Given the theta coefficient lists of every class in one genus (thetas[j][m+1] =
r_j(m)) and a class index i, return the smallest norm within the common precision
that some class represents but class i does not -- an unconditional certificate
that class i is irregular.  Returns 0 if there is none up to the common precision,
which proves nothing (deciding regularity positively is hard: even the classical
ternary classification retains GRH-conditional cases).}
    B := Min([#t : t in thetas]) - 1;
    for m in [1..B] do
        if thetas[i][m+1] eq 0 and exists{t : t in thetas | t[m+1] ne 0} then
            return m;
        end if;
    end for;
    return 0;
end intrinsic;
