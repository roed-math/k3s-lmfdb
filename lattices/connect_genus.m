declare verbose ConnectGenus, 1;

intrinsic GramStringToLat(s::MonStgElt, n::RngIntElt) -> Lat
{Given a string encoding the entries of an n x n Gram matrix, return the corresponding lattice}
    gram := "[" * s[2..#s-1] * "]"; // Switch to square brackets
    gram := Matrix(Rationals(), n, eval gram);
    return LatticeWithGram(gram : CheckPositive := false);
end intrinsic;

function load_genus_data(genus_label)
    genus := AssociativeArray();
    basic_format := Split(Split(Read("genera_basic.format"), "\n")[1], "|");
    adv_format   := Split(Split(Read("genera_advanced.format"), "\n")[1], "|");
    // The genera_advanced file is written by FillGenus with both the basic and
    // advanced genus fields, so read it with the combined format.
    stage_format := AssociativeArray();
    stage_format["basic"] := basic_format;
    stage_format["advanced"] := basic_format cat adv_format;
    for stage in ["basic", "advanced"] do
        genus_data := Split(Split(Read(LabelPath("genera_"*stage, genus_label)), "\n")[1], "|");
        genus_format := stage_format[stage];
        assert #genus_data eq #genus_format;
        for i in [1..#genus_data] do
            genus[genus_format[i]] := genus_data[i];
            if genus_data[i] eq "None" then genus[genus_format[i]] := "\\N"; end if;
        end for;
    end for;
    return genus;
end function;

function lookup_hash_function(genus_hash, rank, nplus)
    fname := LabelPath("genera_hash", rank, nplus, IntegerToString(genus_hash));
    if not OpenTest(fname, "r") then
        return "";
    end if;
    return Read(fname);
end function;

hash_format := Split(Split(Read("lat_hash.format"), "\n")[1], "|");
// Read a genus's lattice-hash file into raw split-field rows.
// NB: these two loaders are kept monomorphic (one returns Assoc, the other
// SeqEnum) on purpose -- a single function with a polymorphic return type
// (as_assoc controlling Assoc vs SeqEnum) makes the Magma package compiler
// reject any intrinsic that calls it more than once with "bad syntax".
function load_hash_rows(genus_hash, rank, nplus)
    fname := LabelPath("lattice_hashes", rank, nplus, IntegerToString(genus_hash));
    if not OpenTest(fname, "r") then return []; end if;
    hash_format := Split(Split(Read("lat_hash.format"), "\n")[1], "|");
    rows := [Split(line, "|") : line in Split(Read(fname), "\n")];
    assert &and[#row eq #hash_format : row in rows];
    return rows;
end function;

// Group the lattices in a genus by their per-lattice hash value.
// Returns Assoc: hash value -> sequence of raw split-field rows.
function load_hash_by_hash(genus_hash, rank, nplus)
    hash_format := Split(Split(Read("lat_hash.format"), "\n")[1], "|");
    hash_pos := Index(hash_format, "hash");
    lats := AssociativeArray();
    for row in load_hash_rows(genus_hash, rank, nplus) do
        // Indefinite lattices have no per-lattice hash (theta/BV hashing does not
        // apply), so their hash column is "\N"; skip them rather than crash.
        if row[hash_pos] eq "\\N" then continue; end if;
        h := StringToInteger(row[hash_pos]);
        if IsDefined(lats, h) then Append(~lats[h], row); else lats[h] := [row]; end if;
    end for;
    return lats;
end function;

// Return the lattices in a genus as a sequence of field-name -> value records.
function load_hash_records(genus_hash, rank, nplus)
    hash_format := Split(Split(Read("lat_hash.format"), "\n")[1], "|");
    recs := [];
    for row in load_hash_rows(genus_hash, rank, nplus) do
        lat := AssociativeArray();
        for i in [1..#row] do
            lat[hash_format[i]] := row[i];
            if row[i] eq "None" then lat[hash_format[i]] := "\\N"; end if;
        end for;
        Append(~recs, lat);
    end for;
    return recs;
end function;

hash_cache := NewStore();

intrinsic HashCache() -> Assoc
{Get an associative array for caching lattices by genus and hash value}
    if not StoreIsDefined(hash_cache, "cache") then
        StoreSet(hash_cache, "cache", AssociativeArray());
    end if;
    return StoreGet(hash_cache, "cache");
end intrinsic;

intrinsic GetHashes(genus_hash::RngIntElt, rank::RngIntElt, nplus::RngIntElt) -> Assoc
{Get an associative array with keys the possible hash values for lattices in the genus and values a sequence of strings matching the hash_format}
    cache := HashCache();
    if IsDefined(cache, genus_hash) then
        return cache[genus_hash];
    end if;
    lats := load_hash_by_hash(genus_hash, rank, nplus);
    cache[genus_hash] := lats;
    StoreSet(hash_cache, "cache", cache);
    return lats;
end intrinsic;

intrinsic FindGenusData(L::Lat) -> Tup
{Given a lattice, find the hash of its genus and the hash function used for that genus}
    genus_hash := HashGenus(L);
    nplus := Signature(GramMatrix(L));
    hash_func := lookup_hash_function(genus_hash, Rank(L), nplus);
    return <genus_hash, hash_func>;
end intrinsic;

intrinsic FindLabel(L::Lat : genus_data:=<>) -> MonStgElt
{Given a lattice, find its label (or null if not in database)}
    if #genus_data eq 0 then
        genus_data := FindGenusData(L);
    end if;
    genus_hash, hash_func := Explode(genus_data);
    by_hash := GetHashes(genus_hash, Rank(L), Signature(GramMatrix(L)));
    if #by_hash eq 0 then
        // No lattices stored for this genus
        return "\\N";
    elif #by_hash eq 1 then
        poss := Values(by_hash)[1];
    else
        h := HashLat(L, genus_hash, hash_func);
        poss := by_hash[h];
    end if;
    if #poss eq 1 then
        found := poss[1];
    else
        n := Rank(L);
        gram_pos := Index(hash_format, "gram");
        found := [];
        for lat in poss do
            M := GramStringToLat(lat[gram_pos], n);
            if IsIsomorphic(L, M) then
                found := lat;
                break;
            end if;
        end for;
        assert #found ne 0;
    end if;
    label_pos := Index(hash_format, "label");
    return found[label_pos];
end intrinsic;

intrinsic LatSortKey(label::MonStgElt) -> Tup
{A tuple that sorts how we want lattices to sort}
    // A summand not found in the database (e.g. out of range) has label "\N";
    // sort those first and keep them as their own group.
    if label eq "\\N" then return <0, 0, 0, label>; end if;
    pieces := Split(label, ".");
    // Sort by: rank, then signature (pos def first), then absolute det, then the label string as tiebreaker
    return <StringToInteger(pieces[1]), -StringToInteger(pieces[2]), StringToInteger(pieces[3]), label>;
end intrinsic;

intrinsic LabelPath(folder::MonStgElt, rank::RngIntElt, nplus::RngIntElt, identifier::MonStgElt : Create := false) -> MonStgElt
{The on-disk path "folder/rank/nplus/identifier" for a lattice identifier where
the lattice has rank and nplus.  Centralises the data directory layout so that a
 future change to the folder scheme only needs editing here.  If Create is true,
 the containing directory is created (mkdir -p) so the path is ready to write to.}
    dir := Sprintf("%o/%o/%o", folder, rank, nplus);
    if Create then
        System("mkdir -p " * dir);
    end if;
    return Sprintf("%o/%o", dir, identifier);
end intrinsic;

intrinsic LabelPath(folder::MonStgElt, label::MonStgElt : Create := false) -> MonStgElt
{The on-disk path "folder/rank/nplus/label" for a lattice or genus label of the
 form rank.nplus.det.... .  Centralises the data directory layout so that a
 future change to the folder scheme only needs editing here.  If Create is true,
 the containing directory is created (mkdir -p) so the path is ready to write to.}
    pieces := Split(label, ".");
    require #pieces ge 2 : "label must have the form rank.nplus....";
    rank := StringToInteger(pieces[1]);
    nplus := StringToInteger(pieces[2]);                   // label is rank.nplus.det.... (see create_genus_label)
    return LabelPath(folder, rank, nplus, label : Create := Create);
end intrinsic;

intrinsic AutOrbits(L::Lat, A::GrpMat, vecs::SeqEnum) -> SeqEnum
{Given a lattice L, an automorphism group A of L and a sequence of vectors which is invariant under the action of A (which could be LatElts or ModTupFldElts), return orbit representatives for the action of A on vecs}
    // This assumes A is in GL_n(Z), with respect to the basis of the lattice
    // created with NaturalAction = false.  vecs may be LatElts (e.g. shortest
    // vectors, whose coordinates come from Coordinates) or ModTupFldElts (e.g.
    // deep holes, which are already coordinate vectors); reduce both to their
    // coordinate vectors and match the group action there, so we never mix types.
    coord := func<v | Type(v) eq LatElt select Vector(Rationals(), Coordinates(v))
                                          else  Vector(Rationals(), Eltseq(v))>;
    cvecs := [coord(v) : v in vecs];
    gens := [ChangeRing(g, Rationals()) : g in Generators(A)];
    Sn := Sym(#vecs);
    perms := [Sn![Index(cvecs, cvecs[i]*g) : i in [1..#vecs]] : g in gens];
    A_perm := sub<Sn|perms>;
    orb_reps := OrbitRepresentatives(A_perm);
    // We throw away the orbit sizes, and replace the index of the
    // representative by the representative itself
    return [vecs[o[2]] : o in orb_reps];
end intrinsic;

intrinsic VoronoiData(L::Lat, A::GrpMat) -> RngIntElt, RngIntElt, SeqEnum[ModTupFldElt], RngIntElt, RngIntElt, RngIntElt
{Given a lattice L and its automorphism group A, find the covering norm, orbit representatives for the deep holes, the number of deep holes, the number of deep hole orbits, and the number of holes}
    cn := CoveringRadius(L);
    cnn := Numerator(cn);
    cnd := Denominator(cn);
    dh := DeepHoles(L);
    reps := AutOrbits(L, A, dh);
    return cnn, cnd, reps, #dh, #reps, #Holes(L);
end intrinsic;

intrinsic IsWellRounded(L::Lat, S::SeqEnum) -> BoolElt
{Given the sequence of shortest vectors in a lattice, test if the lattice is well rounded (ie the rank of their span is equal to the rank of the lattice)}
    return Rank(sub<L | S>) eq Rank(L);
end intrinsic;

intrinsic IsStronglyWellRounded(L::Lat, S::SeqEnum) -> BoolElt
{Given the sequence of shortest vectors in a lattice, test if it contains a basis for the lattice (ie some n of them have coordinate matrix of determinant +-1).}
    n := Rank(L);
    if #S eq 0 then
        return n eq 0;
    end if;
    // Coordinates of each minimal vector wrt a basis of L: integer row vectors.
    C := [ Vector([Integers()| c : c in Coordinates(L, s)]) : s in S ];
    m := #C;
    // Must at least be well rounded.
    if Rank(Matrix(C)) ne n then
        return false;
    end if;
    // Backtracking search for n rows forming a unimodular matrix. We only keep
    // partial sets that are independent and saturated (all elementary divisors
    // 1), since any subset of a basis has this property; a saturated
    // independent set of size n is automatically a basis of L.
    search := function(rows, start)
        k := #rows;
        if k eq n then
            return true;
        end if;
        for i in [start..m] do
            cand := Append(rows, C[i]);
            M := Matrix(cand);
            if Rank(M) eq k + 1 then
                ed := ElementaryDivisors(M);
                if #ed eq k + 1 and &and[ d eq 1 : d in ed ] then
                    if $$(cand, i + 1) then
                        return true;
                    end if;
                end if;
            end if;
        end for;
        return false;
    end function;
    return search([], 1);
end intrinsic;

intrinsic IsMinimalVectorGenerated(L::Lat, S::SeqEnum) -> BoolElt
{Given the sequence of shortest vectors in a lattice, test if the lattice is generated by its minimal vectors}
    return sub<L | S> eq L;
end intrinsic;

intrinsic IsEutactic(L::Lat, S::SeqEnum) -> BoolElt, SeqEnum
{Given the sequence S of shortest vectors in a lattice, test if the quadratic form x.x is a linear combination of the (x.s)^2, for s in S, with positive coefficients.
If it is, return the (exact, rational) eutaxy coefficients in the same order as S; otherwise return an empty sequence.}
    // Writing u_s for the coordinates of the minimal vector s with respect to a
    // basis of L and G for the Gram matrix, the condition
    //     x.x = sum_s c_s (x.s)^2  for all x,  with c_s > 0
    // is equivalent to  sum_s c_s u_s^t u_s = G^{-1}.  This is a linear system in
    // the c_s; the lattice is eutactic iff it has a strictly positive solution.
    n := Rank(L);
    if #S eq 0 then
        return n eq 0, [];
    end if;
    G := ChangeRing(GramMatrix(L), Rationals());
    Ginv := G^-1;
    U := [ Vector(Rationals(), [ c : c in Coordinates(L, s) ]) : s in S ];
    m := #U;
    // Vectorise symmetric matrices on the upper triangle (i <= j); two symmetric
    // matrices are equal iff they agree there.
    idx := [ <i,j> : j in [i..n], i in [1..n] ];
    svec := func< M | Vector(Rationals(), [ M[p[1]][p[2]] : p in idx ]) >;
    A := Matrix([ svec(Transpose(Matrix(u)) * Matrix(u)) : u in U ]);   // m x N
    b := svec(Ginv);                                                    // length N

    // Is the (weak) eutaxy equation c*A = b solvable at all?
    cons, sol := IsConsistent(A, b);
    if not cons then
        return false, [];                  // not even weakly eutactic
    end if;
    K := Basis(Kernel(A));                  // {y : y*A = 0}, the solution directions
    d := #K;
    if d eq 0 then
        // Unique solution: eutactic iff every coefficient is positive (exact).
        if &and[ sol[i] gt 0 : i in [1..m] ] then
            return true, Eltseq(sol);
        end if;
        return false, [];
    end if;

    // General case: the solutions form the affine space sol + <K>.  Decide whether
    // it meets the open positive orthant by maximising the smallest coefficient t
    // via a linear program (c = sol + sum mu_j K_j, with all c_i >= t); the lattice
    // is eutactic iff the optimum t* is positive.  The free variables mu_j and t are
    // modelled as differences of nonnegative ones.  The optimum is finite because
    // the kernel meets the nonnegative orthant only at 0: a positive relation among
    // the rank-one positive-semidefinite matrices u_s^t u_s would force every u_s=0.

    // The decision is exact: t* is rational with denominator at most Dbound (an
    // input-derived Hadamard bound -- an optimal vertex of the (m+1)-variable LP
    // solves a square integer subsystem whose determinant bounds the denominator),
    // so a positive t* is at least 1/Dbound.  We therefore take the working
    // precision, the decision threshold tau and the certificate denominator all
    // from Dbound rather than from fixed constants.  A is integral; clearing the
    // denominators of b leaves every constraint coefficient an integer of size at
    // most Mm, whence any (m+1)-square subdeterminant is at most (sqrt(m+1)*Mm)^(m+1).
    Db := Lcm([ Denominator(x) : x in Eltseq(b) ] cat [1]);
    Mm := Integers() ! Max([ Abs(x) : x in Eltseq(A) ]
                           cat [ Abs(x*Db) : x in Eltseq(b) ] cat [Db, 1]);
    Dbound := Mm^(m+1) * (m+1)^((m+1) div 2 + 1);
    tau := 1/(2*Dbound);
    // Recovering an exact rational of denominator <= Dbound from the LP's floating
    // solution needs ~2*log10(Dbound) digits, not log10(Dbound): two rationals with
    // denominators <= Dbound differ by as little as 1/Dbound^2, so BestApproximation
    // below cannot resolve the vertex coordinate with only log10(Dbound) digits and
    // throws "Overflow or precision loss" once Dbound exceeds ~10^30 (dense rank-6
    // lattices with many minimal vectors).  Use 2x the digits, +30 to guard LP roundoff.
    R := RealField(2*Ceiling(Log(10, 2*Dbound)) + 30);

    nv := 2*d + 2;                          // p_1..p_d, q_1..q_d, tp, tn
    rows := [];  rhs := [];
    for i in [1..m] do
        row := [ R | 0 : k in [1..nv] ];
        for j in [1..d] do
            row[j]   :=  R ! K[j][i];        // coeff of p_j   (mu_j = p_j - q_j)
            row[d+j] := -R ! K[j][i];        // coeff of q_j
        end for;
        row[2*d+1] := -1;  row[2*d+2] := 1;  // -tp + tn   (t = tp - tn)
        Append(~rows, row);
        Append(~rhs, [ -R ! sol[i] ]);       // c_i - t >= 0
    end for;
    LP := LPProcess(R, nv);
    AddConstraints(LP, Matrix(R, m, nv, &cat rows), Matrix(R, m, 1, &cat rhs) : Rel := "ge");
    obj := [ R | 0 : k in [1..nv] ];
    obj[2*d+1] := 1;  obj[2*d+2] := -1;      // maximise t = tp - tn
    SetObjectiveFunction(LP, Matrix(R, 1, nv, obj));
    SetMaximiseFunction(LP, true);
    v := Eltseq(Solution(LP));
    tstar := v[2*d+1] - v[2*d+2];
    if tstar le R ! tau then
        return false, [];                    // t* <= 0 (any positive t* would exceed tau)
    end if;

    // t* > 0: a strictly positive solution exists.  Recover the exact rational one
    // (its denominator is at most Dbound) and verify its positivity, so the YES
    // answer and the returned coefficients are certified exactly.  Any c in
    // sol + <K> automatically satisfies c*A = b.
    lam := sol;
    for j in [1..d] do
        lam +:= BestApproximation(v[j] - v[d+j], Dbound) * Vector(Rationals(), Eltseq(K[j]));
    end for;
    return true, Eltseq(lam);
end intrinsic;

intrinsic PerfectionDefect(L::Lat, S::SeqEnum) -> RngIntElt
{Given a lattice L of rank n, return n(n+1)/2 - rank of span (in M_n(R)) of the rank 1 matrices (s_i s_i^t), with s in S.}
    n := Rank(L);
    N := n*(n+1) div 2;
    if #S eq 0 then
        return N;
    end if;
    U := [ Vector(Rationals(), [ c : c in Coordinates(L, s) ]) : s in S ];
    idx := [ <i,j> : j in [i..n], i in [1..n] ];
    svec := func< M | Vector(Rationals(), [ M[p[1]][p[2]] : p in idx ]) >;
    A := Matrix([ svec(Transpose(Matrix(u)) * Matrix(u)) : u in U ]);
    return N - Rank(A);
end intrinsic;

intrinsic IsAdditivelyDecomposableByRankOne(L::Lat) -> BoolElt
{Returns true if L is provably additively decomposable by splitting off a rank-1
form, i.e. if its dual lattice has a nonzero vector of norm <= 1 (rank 1: Gram
entry >= 2).  This is only a SUFFICIENT condition for additive decomposability,
NOT a characterisation: e.g. the rank-15 complement of (x,y) in E8 + E8 with
Q(x)=2, Q(y)=4 is additively decomposable yet has dual minimum 4/3 > 1.  Hence a
"false" here does NOT imply additively indecomposable; the full decision is much
harder (see Wang, Acta Math. Sinica 41 (2025) 908-924, and Plesken).}
    if Rank(L) eq 1 then
        return GramMatrix(L)[1][1] ge 2;
    end if;
    return Minimum(Dual(L : Rescale := false)) le 1;
end intrinsic;

// ---- helpers for additive (in)decomposability of low-discriminant lattices ----

function lat_adjugate(G)
    GQ := ChangeRing(G, Rationals());
    Gi := GQ^-1;
    return ChangeRing(Determinant(GQ) * Gi, Integers());   // det(G) * G^-1, integral
end function;

// The primitive rank-(n-1) sublattices of L of discriminant t, as a sequence of
// lattices.  Such a sublattice is w^perp cap L for a dual vector w, and its
// discriminant equals w * adj(G) * w^t; so we range over the norm-t vectors of
// the lattice with Gram matrix adj(G).  (For squarefree t every rank-(n-1)
// sublattice of discriminant t is primitive, so this is exhaustive there.)
function rank_nm1_sublattices(L, t)
    G := ChangeRing(GramMatrix(L), Integers());
    AL := LatticeWithGram(lat_adjugate(G));
    res := [];
    for p in ShortVectors(AL, t, t) do
        w := Matrix(Integers(), Ncols(G), 1, Eltseq(p[1]));
        B := Matrix(Integers(), [ Eltseq(b) : b in Basis(Kernel(w)) ]);
        Append(~res, LatticeWithGram(B * G * Transpose(B)));
    end for;
    return res;
end function;

function lat_orth_indecomposable(M)
    return Rank(M) le 1 or #OrthogonalDecomposition(M) eq 1;
end function;

function has_indec_rank_nm1_sublattice(L, t)
    return exists{ M : M in rank_nm1_sublattices(L, t) | lat_orth_indecomposable(M) };
end function;

// Whether L has orthogonal sublattices M1 perp M2 with disc M1 = a, disc M2 = b
// spanning a rank-(n-1) sublattice (which then has discriminant a*b).
function has_orth_pair_sublattice(L, a, b)
    for M in rank_nm1_sublattices(L, a*b) do
        ds := Rank(M) le 1 select [ Determinant(M) ]
              else [ Determinant(C) : C in OrthogonalDecomposition(M) ];
        if a eq b then
            if #[ x : x in ds | x eq a ] ge 2 then return true; end if;
        elif a in ds and b in ds then
            return true;
        end if;
    end for;
    return false;
end function;

// Whether L has a decomposable (necessarily unimodular) integral overlattice of
// index 2 -- the "represented by a decomposable unimodular lattice" condition.
function has_decomposable_unimodular_overlattice(L)
    Ld := Dual(L : Rescale := false);
    A, q := quo< Ld | L >;
    for a in A do
        if a eq A!0 or 2*a ne A!0 then continue; end if;
        M := sub< Ld | Basis(L) cat [a @@ q] >;
        if IsIntegral(M) and Index(M, L) eq 2 and #OrthogonalDecomposition(M) gt 1 then
            return true;
        end if;
    end for;
    return false;
end function;

intrinsic SatisfiesPleskenIII1(L::Lat) -> BoolElt
{Plesken's sufficient condition (Prop. III.1) for additive indecomposability: L is
additively indecomposable if (i) min(L*) > 1 (L is a "block form"), and (ii) the
sublattice L(<=3) generated by the vectors of norm <= 3 is orthogonally
indecomposable with dim L - dim L(<=3) <= 5.  Returns true => additively
indecomposable; false is inconclusive.}//'
    n := Rank(L);
    if n eq 0 then return true; end if;
    if Minimum(Dual(L : Rescale := false)) le 1 then
        return false;
    end if;
    sv := ShortVectors(L, 3);
    if #sv eq 0 then
        return n le 5;                      // L(<=3) = 0, corank = n
    end if;
    R := sub<L | [ v[1] : v in sv ]>;
    return #OrthogonalDecomposition(R) eq 1 and (n - Rank(R)) le 5;
end intrinsic;

intrinsic IsAdditivelyIndecomposable(L::Lat) -> BoolElt, BoolElt
{Decide, when possible, whether L is additively indecomposable (its Gram matrix is
not a sum of two nonzero positive semidefinite integral matrices).  Returns
<is_indecomposable, is_determined>; the second value is false when the answer is
not known.  The decision is complete for rank <= 8 (Mordell: none in ranks 2-5;
Plesken Thm. III.4: exactly E_n in ranks 6,7,8; only (Z,1) in rank 1) and, for
larger rank, whenever the lattice is orthogonally decomposable, has a rank-1
additive split, has discriminant 2-5 (Wang Thms. 2.11/2.12/2.16), or meets
Plesken's sufficient condition III.1.  Discriminants >= 6 are left undetermined,
as is the rare discriminant 4/5 case with a unimodular rank-(n-1) sublattice.} //'
    n := Rank(L);
    if n eq 0 then return true, true; end if;
    d := Determinant(L);
    // Complete classification in low rank.
    if n eq 1 then return d eq 1, true; end if;
    if n in {2,3,4,5} then return false, true; end if;          // Mordell
    if n in {6,7,8} then return IsIsometric(L, Lattice("E", n)), true; end if;  // Plesken III.4
    // Rank >= 9: rigorous partial decision.
    if #OrthogonalDecomposition(L) gt 1 then                    // an orthogonal sum is an additive sum
        return false, true;
    end if;
    if Minimum(Dual(L : Rescale := false)) le 1 then            // splits off a rank-1 form
        return false, true;
    end if;
    if d eq 2 then return true, true; end if;                   // Wang Thm. 2.11
    if d eq 3 then                                              // Wang Thm. 2.12
        return #rank_nm1_sublattices(L, 2) eq 0, true;          //   decomposable iff a rank-(n-1) disc-2 sublattice
    end if;
    if d eq 4 or d eq 5 then                                    // Wang Thm. 2.16
        cond1 := (d eq 4) select has_decomposable_unimodular_overlattice(L)
                            else (#rank_nm1_sublattices(L, 2) gt 0);
        if cond1 or has_indec_rank_nm1_sublattice(L, d-2)
                 or has_indec_rank_nm1_sublattice(L, d-1)
                 or has_orth_pair_sublattice(L, 2, d-2) then
            return false, true;                                 // additively decomposable
        end if;
        // The conditions above use the primitive rank-(n-1) sublattices; the only
        // gap is the non-primitive disc-(d-1) ones, which exist iff L has a
        // unimodular rank-(n-1) sublattice.  When it has none, the analysis is
        // complete and L is additively indecomposable.
        if #rank_nm1_sublattices(L, 1) eq 0 then
            return true, true;
        end if;
    end if;
    if SatisfiesPleskenIII1(L) then return true, true; end if;  // Plesken III.1 (min(L*) > 1 already checked)
    return false, false;   // discriminant >= 6 (or the rare disc 4/5 gap): undetermined
end intrinsic;

intrinsic SaveVdat(lat::Assoc)
{Save the voronoi data to disk for later use in the decomposable case}
    path := LabelPath("voronoi", lat["label"] : Create);
    contents := Sprintf("%o|%o|%o|%o|%o", lat["covering_norm_num"], lat["covering_norm_den"], lat["deep_hole_count"], lat["deep_hole_orbit_count"], lat["hole_count"]);
    Write(path, contents : Overwrite);   // Overwrite, not append: connect must be re-runnable (a stale file would otherwise double up and break LoadVdat's parse)
end intrinsic;

intrinsic LoadVdat(labels::SeqEnum[MonStgElt]) -> SeqEnum[Tup]
{Given a sequence of lattice labels, load Voronoi data as a sequence of tuples (covering norm, num deep holes, num deep hole orbits, num holes).  If any not available, return empty sequence instead.  Connect is run rank-ascending (see run_all.py), so a decomposable lattice's strictly-lower-rank orthogonal factors are always already on disk here.}
    ans := [];
    for label in labels do
        if label eq "\\N" then return []; end if;   // a factor not in the database
        fname := LabelPath("voronoi", label);
        if not OpenTest(fname, "r") then
            return [];
        end if;
        cnn, cnd, ndh, ndho, nh := Explode(Split(Read(fname), "|"));
        Append(~ans, <StringToInteger(cnn), StringToInteger(cnd), StringToInteger(ndh), StringToInteger(ndho), StringToInteger(nh)>);
    end for;
    return ans;
end intrinsic;

// Field order of the shortest/<label> files, matching what ConnectGenus writes.
sv_fields := ["minimum", "shortest", "is_well_rounded", "is_minimal_vector_generated",
              "is_strongly_well_rounded", "is_eutactic", "is_strongly_eutactic",
              "t_design", "perfection_defect", "is_perfect", "is_strongly_perfect"];

function parse_sv_field(s)
    // Recover a stored value: booleans, integers, the null marker, or a raw string.
    if s eq "\\N" then
        return "\\N";
    elif s eq "true" then
        return true;
    elif s eq "false" then
        return false;
    elif #s gt 0 and forall{ i : i in [1..#s] | s[i] in "0123456789-" }
                 and exists{ i : i in [1..#s] | s[i] in "0123456789" } then
        return StringToInteger(s);
    else
        return s;
    end if;
end function;

intrinsic SaveSVdat(lat::Assoc)
{Save the shortest vector data to disk for later use in the decomposable case}
    path := LabelPath("shortest", lat["label"] : Create);
    // to_postgres keeps "shortest" (a sequence of vectors) on a single line;
    // "%o" would print it across several lines and LoadSVdat reads only the first.
    contents := Sprintf("%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o", lat["minimum"], to_postgres(lat["shortest"]), lat["is_well_rounded"], lat["is_minimal_vector_generated"], lat["is_strongly_well_rounded"], lat["is_eutactic"], lat["is_strongly_eutactic"], lat["t_design"], lat["perfection_defect"], lat["is_perfect"], lat["is_strongly_perfect"]);
    Write(path, contents : Overwrite);   // Overwrite, not append: connect must be re-runnable (a stale file would otherwise double up and break LoadSVdat's parse)
end intrinsic;

intrinsic LoadSVdat(labels::SeqEnum[MonStgElt]) -> SeqEnum
{Given a sequence of lattice labels, load short-vector data for each as an associative array keyed by property name (minimum, shortest, is_well_rounded, ...), with booleans and integers parsed and "\N" denoting a missing value.  If any file is not available, return an empty sequence instead.  Connect is run rank-ascending (see run_all.py), so a decomposable lattice's strictly-lower-rank orthogonal factors are always already on disk here.}
    ans := [];
    for label in labels do
        if label eq "\\N" then return []; end if;   // a factor not in the database
        fname := LabelPath("shortest", label);
        if not OpenTest(fname, "r") then
            return [];
        end if;
        vals := Split(Split(Read(fname), "\n")[1], "|");
        assert #vals eq #sv_fields;
        dat := AssociativeArray();
        for i in [1..#sv_fields] do
            dat[sv_fields[i]] := parse_sv_field(vals[i]);
        end for;
        Append(~ans, dat);
    end for;
    return ans;
end intrinsic;

intrinsic TimeoutAssign(~D::Assoc, key::MonStgElt, func::Program, inp::Tup, timeout::RngIntElt : Parameters := [])
{Compute value using given function, then store in D}
    success, out, elapsed := TimeoutCall(timeout, func, inp, 1 : Parameters := Parameters);
    if success then
        D[key] := out[1];
    else
        D[key] := "\\N";
    end if;
end intrinsic;

intrinsic ShortestVectorCoords(L::Lat) -> SeqEnum
{The coordinate sequences of the shortest vectors of L (up to sign). Returns plain
 integer sequences rather than lattice vectors so the result can be serialized
 across a TimeoutCall fork.}
    return [Eltseq(v) : v in ShortestVectors(L)];
end intrinsic;

intrinsic AutOrbitCoords(L::Lat, A::GrpMat, vecs::SeqEnum) -> SeqEnum
{As AutOrbits, but returns the coordinate sequences of the orbit representatives
 (serializable across a TimeoutCall fork).}
    return [Eltseq(v) : v in AutOrbits(L, A, vecs)];
end intrinsic;

intrinsic DeriveDecomposableFields(~lat::Assoc, n::RngIntElt)
{For a decomposable definite lattice of rank n, derive the geometric fields that combine
over an orthogonal direct sum (deep-hole counts and the well-rounded / eutactic / perfect
booleans) from its orthogonal factors' Voronoi and short-vector data, which must already
be on disk.  Fields with no orthogonal-sum rule are left as they are ("\N").}
    vdat := LoadVdat(lat["orthogonal_factors"]);
    if #vdat gt 0 then
        cnorm := 0;
        num_deep_holes := 1;
        num_deep_hole_orbits := 1;
        num_holes := 1;
        for i in [1..#vdat] do
            cnn, cnd, ndh, ndho, nh := Explode(vdat[i]);
            m := lat["orthogonal_multiplicities"][i];
            cnorm +:= m * (cnn/cnd);
            num_deep_holes *:= ndh^m;
            num_deep_hole_orbits *:= Binomial(ndho + m - 1, m);
            num_holes *:= nh^m;
        end for;
        lat["covering_norm"] := cnorm;
        lat["deep_hole_count"] := num_deep_holes;
        lat["deep_hole_orbit_count"] := num_deep_hole_orbits;
        lat["hole_count"] := num_holes;
    end if;

    // Derive the short-vector properties that combine cleanly over an orthogonal
    // direct sum.  Minimal vectors live entirely in the factors of smallest minimum,
    // so a property such as well-roundedness holds for L iff every factor attains that
    // minimum ("active") and has the property itself.
    svdat := LoadSVdat(lat["orthogonal_factors"]);
    if #svdat gt 0 then
        mults := lat["orthogonal_multiplicities"];
        ranks := [ StringToInteger(Split(f, ".")[1]) : f in lat["orthogonal_factors"] ];
        minimum := Minimum([ d["minimum"] : d in svdat ]);
        active := [ d["minimum"] eq minimum : d in svdat ];
        all_active := &and active;
        // "All factors active and each has the property": false if some factor is
        // inactive or lacks it, "\N" if any value is unknown.
        derive := function(key)
            if not all_active then return false; end if;
            vals := [ d[key] : d in svdat ];
            if exists{ v : v in vals | v cmpeq "\\N" } then return "\\N"; end if;
            return &and vals;
        end function;
        lat["is_well_rounded"]             := derive("is_well_rounded");
        lat["is_minimal_vector_generated"] := derive("is_minimal_vector_generated");
        lat["is_strongly_well_rounded"]    := derive("is_strongly_well_rounded");
        lat["is_eutactic"]                 := derive("is_eutactic");
        // Perfection defect: n(n+1)/2 minus the dimension spanned by the minimal
        // vectors, which sit block-diagonally inside the active factors, so the
        // spanned dimension adds up over those factors.
        if forall{ i : i in [1..#svdat] | not active[i]
                       or svdat[i]["perfection_defect"] cmpne "\\N" } then
            perfrank := &+[ Integers() |
                mults[i] * (ranks[i]*(ranks[i]+1) div 2 - svdat[i]["perfection_defect"])
                : i in [1..#svdat] | active[i] ];
            lat["perfection_defect"] := n*(n+1) div 2 - perfrank;
            lat["is_perfect"] := (lat["perfection_defect"] eq 0);
        end if;
        // The remaining shell properties (shortest, is_strongly_eutactic, t_design,
        // is_strongly_perfect) have no simple orthogonal-sum rule and are left "\N".
    end if;
end intrinsic;

function SplitPipe(s)
    // Split on "|" preserving empty fields.  Magma's Split drops empty tokens, which
    // misaligns an advanced record whose columns can be empty strings (e.g. an empty
    // root system serialises to "").
    parts := [];
    cur := "";
    for i in [1..#s] do
        if s[i] eq "|" then Append(~parts, cur); cur := ""; else cur := cur cat s[i]; end if;
    end for;
    Append(~parts, cur);
    return parts;
end function;

function ParsePGList(s)
    // Inverse of to_postgres for a flat sequence: "{a,b,c}" -> ["a","b","c"].  Labels
    // contain "." but never "," so a plain split is safe; "{}" and "\N" -> [].
    if s eq "\\N" or #s lt 2 then return []; end if;
    inner := s[2..#s-1];   // strip the { }
    if inner eq "" then return []; end if;
    return Split(inner, ",");
end function;

intrinsic ConnectGenusDerive(label::MonStgElt)
{Consume pass of the two-pass connect: for each decomposable lattice in this genus,
derive the orthogonal-sum geometric fields the produce pass left "\N" from the factors'
data (now guaranteed on disk) and patch them into the lattice's advanced record.}
    SetColumns(0);
    advanced_format := Split(Split(Read("lat_advanced.format"), "\n")[1], "|");
    sl := Split(label, ".");
    n := StringToInteger(sl[1]);
    s := StringToInteger(sl[2]);
    if s ne n then return; end if;   // indefinite: no definite geometry to derive
    ii := Index(advanced_format, "is_indecomposable");
    of_i := Index(advanced_format, "orthogonal_factors");
    om_i := Index(advanced_format, "orthogonal_multiplicities");
    // The fields DeriveDecomposableFields fills in for a decomposable lattice.
    derived := ["deep_hole_count", "deep_hole_orbit_count", "hole_count",
                "is_well_rounded", "is_minimal_vector_generated", "is_strongly_well_rounded",
                "is_eutactic", "perfection_defect", "is_perfect"];
    lats := load_hash_records(CollapseIntList(StringToBytes(label)), n, s);
    for i in [1..#lats] do
        path := LabelPath("lattice_advanced_data", lats[i]["label"]);
        if not OpenTest(path, "r") then continue; end if;   // produce pass never wrote it (e.g. an upstream fill failure) -- nothing to derive
        fields := SplitPipe(Split(Read(path), "\n")[1]);   // first line, empties preserved
        error if #fields ne #advanced_format, "advanced record field count mismatch for", lats[i]["label"];
        if fields[ii] ne "false" then continue; end if;   // only decomposable lattices
        lat := AssociativeArray();
        lat["orthogonal_factors"] := ParsePGList(fields[of_i]);
        lat["orthogonal_multiplicities"] := [ StringToInteger(x) : x in ParsePGList(fields[om_i]) ];
        for k in derived do lat[k] := "\\N"; end for;
        DeriveDecomposableFields(~lat, n);
        for k in derived do
            fields[Index(advanced_format, k)] := Sprintf("%o", to_postgres(lat[k]));
        end for;
        Write(path, Join(fields, "|") : Overwrite);
    end for;
end intrinsic;

intrinsic ConnectGenus(label::MonStgElt : timeout := 1800, derive_decomposable := true)
{Fill in lattice data that requires working with lattices in different genera.  When
derive_decomposable is false, decomposable lattices' orthogonal-sum-derived geometric
fields are left "\N" for a later ConnectGenusDerive (consume) pass, so that every
factor's data is guaranteed on disk before it is read (the "produce" pass of the
two-pass connect in run_all.py).}
    SetColumns(0);
    advanced_format := Split(Split(Read("lat_advanced.format"), "\n")[1], "|");
    basic_format := Split(Split(Read("lat_basic.format"), "\n")[1], "|");
    atomic_names := LoadAtomicNames();              // stage-4 atomic names (run_basic_names)
    name_i := Index(advanced_format, "name");
    // Everything ConnectGenus needs is determined by the label: rank and nplus are
    // its first two components, and the genus (hence its hash) is recovered by
    // GenusSymbolFromLabel.  HashGenus is canonical, so this reproduces the hash
    // FillGenus used to name the lattice_hashes file; scale is the gcd of each
    // lattice's Gram entries (computed in the loop below).
    sl := Split(label, ".");
    n := StringToInteger(sl[1]);
    s := StringToInteger(sl[2]);
    // Geometric invariants (minimum, shortest vectors, density, Voronoi data, the
    // orthogonal decomposition, root/norm-1 sublattices) require a positive
    // definite form; indefinite lattices (nplus < rank) leave them all "\N".
    definite := (s eq n);
    // The genus hash is HashGenus(genus) = hash(CreateGenusLabel(genus)), and the
    // label IS that canonical label, so hash it directly.  This avoids
    // reconstructing the symbol via GenusSymbolFromLabel (and its round-trip bugs).
    lats := load_hash_records(CollapseIntList(StringToBytes(label)), n, s);
    if #lats gt 0 then
        to_per_rep := timeout div #lats + 1;
    end if;
    // Theta coefficient lists for every class in the genus: the rank-1 regularity
    // scan needs all classmates at once, and the universality ladder reuses the
    // stored series instead of recomputing them.
    genus_thetas := [];
    if definite then
        theta_col := Index(basic_format, "theta_series");
        for j in [1..#lats] do
            pieces := Split(Split(Read(LabelPath("lattice_basic_data", lats[j]["label"])), "\n")[1], "|");
            Append(~genus_thetas, [Integers() | StringToInteger(x) : x in ParsePGList(pieces[theta_col])]);
        end for;
    end if;

    for i in [1..#lats] do
        lat := lats[i];
        // The hash records only carry label|hash|minimum|gram; ConnectGenus needs
        // basic fields (aut_group, is_even, ...) computed in FillGenus, so merge in
        // the full basic record (without clobbering the hash fields).
        basic_pieces := Split(Split(Read(LabelPath("lattice_basic_data", lat["label"])), "\n")[1], "|");
        for j in [1..#basic_format] do
            if not IsDefined(lat, basic_format[j]) then
                lat[basic_format[j]] := basic_pieces[j];
            end if;
        end for;
        L := GramStringToLat(lat["gram"], n);
        scale := GCD([Integers() | x : x in Eltseq(GramMatrix(L))]);   // = gcd(rep), a genus invariant
        D := Dual(L);

        // Default every advanced column to "\N"; the computations below overwrite
        // the ones that apply (most are definite-only, so indefinite lattices keep
        // the defaults).
        for col in advanced_format do
            if not IsDefined(lat, col) then lat[col] := "\\N"; end if;
        end for;

        lat["dual_det"] := Determinant(D);
        lat["dual_label"] := FindLabel(D);
        if definite then
            lat["dual_density"] := Density(D);
            lat["dual_hermite"] := HermiteNumber(D);
            lat["dual_kissing"] := KissingNumber(D);
        end if;
        // dual_theta_series / dual_theta_prec are computed and stored in FillGenus
        // (lat_basic.format), so they are not recomputed here.

        // Criterion-set layer of universality/regularity (definite only); the
        // db-relative layer runs once the Oscar embedding stage has produced
        // lat_inclusions (oscar/INTEGRATION_PLAN.md section 10).  Values the
        // intrinsics cannot decide keep the "\N" defaults set above.
        if definite then
            // The stored automorphism group (computed in FillGenus on the same gram)
            // lets the embedding searches collapse their top level to orbits.
            autg := 0;
            if lat["aut_group"] ne "\\N" then
                try
                    autg := StringToGroup(lat["aut_group"]);
                catch e
                    autg := 0;
                end try;
            end if;
            udec, uval, eudec, euval := LatticeUniversalityData(L : Theta := genus_thetas[i], Aut := autg,
                                                                TotalSeconds := Max(60, to_per_rep));
            if udec then
                lat["universality"] := uval;
            end if;
            if eudec then
                lat["even_universality"] := euval;
            end if;
            if #lats eq 1 then
                // genus-representation means representation by some class, so a
                // class number 1 lattice is regular at every rank up to its own
                lat["is_regular"] := true;
                lat["regularity"] := n;
            else
                wit := RankOneIrregularityWitness(genus_thetas, i);
                if wit gt 0 then
                    lat["is_regular"] := false;
                    lat["regularity"] := 0;
                end if;
            end if;
        end if;

        if scale eq 1 then
            lat["primitive_scaling"] := "\\N";
        else
            lat["primitive_scaling"] := FindLabel(LatticeWithGram(GramMatrix(L) div scale : CheckPositive := false));
        end if;

        if definite then
        summands := OrthogonalDecompositionFaster(L);
        lat["is_indecomposable"] := (#summands eq 1);
        summand_labels := [FindLabel(M) : M in summands];
        collected := CountFibers(summand_labels, func<x | LatSortKey(x)>);   // CountFibers wants a UserProgram, not an intrinsic
        lat["orthogonal_factors"] := [fib[1][4] : fib in collected];
        lat["orthogonal_multiplicities"] := [fib[2] : fib in collected];
        lat["name"] := LatticeName(lat["label"], lat["orthogonal_factors"],
                                   lat["orthogonal_multiplicities"], atomic_names, name_i);

        for col in ["covering_norm_num", "covering_norm_den", "deep_holes", "deep_hole_count", "deep_hole_orbit_count", "hole_count"] do
            lat[col] := "\\N"; // Overwritten below if possible
        end for;
        for col in ["shortest", "is_well_rounded", "is_minimal_vector_generated", "is_strongly_well_rounded", "is_eutactic", "is_strongly_eutactic", "t_design", "perfection_defect", "is_perfect", "is_strongly_perfect"] do
            lat[col] := "\\N"; // Overwritten below if possible
        end for;
        if lat["is_indecomposable"] then
            aut_group := (lat["aut_group"] cmpne "\\N") select StringToGroup(lat["aut_group"]) else 0;
            // aut_group is the integer 0 when the automorphism group is unknown ("\N",
            // e.g. an AutomorphismGroup timeout -- more frequent at high rank).
            // VoronoiData and AutOrbitCoords take a GrpMat, so passing 0 raises "Bad
            // argument types"; guard on the group being present and otherwise leave
            // those fields "\N" (they cannot be computed without it anyway).  (tDesign
            // takes A := 0 as its default and handles the no-group case itself.)
            have_aut := Type(aut_group) eq GrpMat;
            // In addition to the timeout, we may want to impose a rank limit
            if have_aut then
                success, vdat, elapsed := TimeoutCall(timeout, VoronoiData, <L, aut_group>, 6);
                if success then
                    lat["covering_norm_num"], lat["covering_norm_den"], lat["deep_holes"], lat["deep_hole_count"], lat["deep_hole_orbit_count"], lat["hole_count"] := Explode(vdat);
                    // We write the data to a file for loading in the decomposable case
                    SaveVdat(lat);
                end if;
            end if;

            has_sv, S, elapsed := TimeoutCall(timeout, ShortestVectorCoords, <L>, 1);
            if has_sv then
                // magma returns representatives up to +-; ShortestVectorCoords
                // returns coordinate sequences (lattice vectors are not
                // serializable across the TimeoutCall fork), so rebuild them in L
                half := [L ! c : c in S[1]];
                S := half cat [-v : v in half];
                if have_aut then
                    TimeoutAssign(~lat, "shortest", AutOrbitCoords, <L, aut_group, S>, timeout);
                end if;
                TimeoutAssign(~lat, "is_well_rounded", IsWellRounded, <L, S>, timeout);
                TimeoutAssign(~lat, "is_minimal_vector_generated", IsMinimalVectorGenerated, <L, S>, timeout);
                TimeoutAssign(~lat, "is_strongly_well_rounded", IsStronglyWellRounded, <L, S>, timeout);
                has_eutaxy, eutaxy_data, elapsed := TimeoutCall(timeout, IsEutactic, <L, S>, 2);
                if has_eutaxy then
                    lat["is_eutactic"], eutaxy := Explode(eutaxy_data);
                    lat["is_strongly_eutactic"] := lat["is_eutactic"] and (#Set(eutaxy) eq 1);
                end if;
                // tDesign only needs reps up to +/-
                TimeoutAssign(~lat, "t_design", tDesign, <L, half>, timeout : Parameters := [<"A", aut_group>]);
                if lat["t_design"] cmpne "\\N" then
                    if lat["t_design"] ge 2 then
                        // Venkov: minimal vectors forming a 2-design => strongly
                        // eutactic.  This is definitive, so trust it rather than the
                        // direct IsEutactic derivation above, whose eutaxy-coefficient
                        // equality test (#Set(eutaxy) eq 1) can spuriously disagree
                        // (e.g. A4, which is strongly eutactic).
                        lat["is_eutactic"] := true;
                        lat["is_strongly_eutactic"] := true;
                    end if;
                    if lat["t_design"] ge 4 then
                        lat["is_strongly_perfect"] := true;
                        lat["is_perfect"] := true;
                    end if;
                end if;
                TimeoutAssign(~lat, "perfection_defect", PerfectionDefect, <L, S>, timeout);
                if lat["perfection_defect"] cmpne "\\N" then
                    if lat["is_perfect"] cmpeq "\\N" then
                        lat["is_perfect"] := (lat["perfection_defect"] eq 0);
                    else
                        assert lat["is_perfect"] eq (lat["perfection_defect"] eq 0);
                    end if;
                end if;
                SaveSVdat(lat);
            end if;
        elif derive_decomposable then
            // Consumer step: derive this decomposable lattice's fields from its
            // orthogonal factors' data.  Skipped when derive_decomposable is false
            // (the "produce" pass), which leaves those fields "\N" for a later
            // consume pass -- see the two-pass connect in run_all.py.
            DeriveDecomposableFields(~lat, n);
        end if;
        end if;   // definite: orthogonal decomposition / Voronoi / shortest vectors

        // Decided when the theory is complete (rank <= 8, disc 2, Plesken III.1, ...);
        // discriminant 3-5 and the general high-rank case are left as "\N".
        // is_additively_indecomposable is a positive-definite notion (whether the
        // Gram matrix is a sum of two nonzero PSD integral matrices).  For indefinite
        // L the low-rank branches return meaningless values and the rank 6-8 branch
        // feeds L to IsIsometric, which requires positive-definiteness and fails with
        // "Numerical GSO does not seem to be positive definite".  Only decide it in
        // the definite case; leave it "\N" (undetermined) otherwise.
        if definite then
            ai_ok, ai_data := TimeoutCall(timeout, IsAdditivelyIndecomposable, <L>, 2);
            if ai_ok then
                ai_val, ai_known := Explode(ai_data);
                lat["is_additively_indecomposable"] := ai_known select ai_val else "\\N";
            else
                lat["is_additively_indecomposable"] := "\\N";
            end if;
        else
            lat["is_additively_indecomposable"] := "\\N";
        end if;

        if lat["is_even"] eq "T" then
            lat["even_sublattice"] := "\\N";
        else
            lat["even_sublattice"] := FindLabel(EvenSublattice(L));
        end if;

        if definite then
        sv1 := ShortVectors(L, 1, 1);
        R1 := sub<L|[x[1] : x in sv1]>;
        lat["norm1_rank"] := Rank(R1);
        if Rank(R1) eq 0 or Rank(R1) eq n then
            lat["norm1_complement"] := "\\N";
        else
            lat["norm1_complement"] := FindLabel(OrthogonalComplementFaster(L, R1));
        end if;

        sv2 := ShortVectors(L, 1, 2);
        R := sub<L|[x[1] : x in sv2]>;
        lat["root_sublattice"] := RootString(R);

        if Rank(R) eq 0 or Rank(R) eq n then
            lat["root_complement"] := "\\N";
        else
            lat["root_complement"] := FindLabel(OrthogonalComplementFaster(L, R));
        end if;
        end if;   // definite: norm-1 and root sublattices

        lat["is_algebraic"] := "\\N"; // TODO Eran: We are not sure about what to do here - see below.
        /*
        
        1. The definition of an "algebraic" lattice is unclear at the moment.
        Here are several options: A lattice is algebraic iff it is (isometric to)
        Option (A): a free Z-module in a number field K with the bilinear trace form (x,y) -> Tr(xy)
        Option (B): a free Z-module in a number field K with a bilinear trace form (x,y) -> Tr(a xy) with a in K
        Option (C): a free Z-module in a number field K with a bilinear trace form (x,y) -> Tr(x sigma(y)) with sigma in Aut(K)
        
        Option (B) would identify almost all lattices as algebraic, and 
        Option (C) allows for Tr(x xbar) on CM fields

        One can also decide to restrict to the following classes of Z-modules:
        (a) full rank Z-modules (rank(L) = [K:Q])
        (b) fractional R-ideals for some order R in K
        (c) fractional Z_K-ideals
        (d) Orders in K
        (e) Only Z_K

        2. We then have to decide how to check that. One option is to run over 
        number fields in the LMFDB, and run over sublattices up to appropriate index. 

        */


        // lat carries input fields not in lat_advanced.format (the hash record's
        // gram/hash/minimum and the basic fields merged in above); drop them so
        // only the advanced columns remain.
        for k in Keys(lat) do
            if k notin advanced_format then Remove(~lat, k); end if;
        end for;
        error if Keys(lat) ne Set(advanced_format), [k : k in advanced_format | k notin Keys(lat)], [k : k in Keys(lat) | k notin advanced_format];
        output := Join([Sprintf("%o", to_postgres(lat[k])) : k in advanced_format], "|");
        Write(LabelPath("lattice_advanced_data", lat["label"] : Create), output : Overwrite);
    end for;
end intrinsic;
