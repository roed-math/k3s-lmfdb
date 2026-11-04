// Deciding whether classes of a genus have equal theta series from far fewer
// coefficients than the Sturm bound, using the genus's p-neighbour graph and Hecke
// operators.
//
// Notation.  L is positive definite of rank n >= 3 with Gram matrix A.  The pipeline's
// theta series is f_L = sum_v q^(v,v).  For an even lattice we work with
// g_L = sum_v q^((v,v)/2) (so f_L(q) = g_L(q^2)), whose level is that of A; for an odd
// one with f_L, of level N the level of 2A (ThetaFormData).  Either way the form lies in
// M_(n/2)(Gamma_0(N), chi) with
//     chi(d) = ((-1)^(n/2) det(S) / d)   (n even),     chi(d) = (2 det(S) / d)   (n odd),
// S the even Gram matrix (A or 2A); half-integral weight forms transform like theta^n
// times chi, as in Magma's HalfIntegralWeightForms.  For two classes of one genus the
// difference of theta series is a cusp form (Siegel).
//
// Eichler's commutation relation.  For a prime p not dividing N, with
// Adj = AdjacencyMatrix(G, p) (Adj[j, i] = the number of p-neighbours of L_i isometric
// to L_j: the source class indexes the column),
//     sum_j Adj[j, i] theta_j = (T(p^2) + c_n(p)) theta_i,    c_n(p) = p + p^2 + ... + p^(n-3)
// (c_2 = -1, c_3 = 0), where in integral weight k, T(p^2) = T(p)^2 - chi(p) p^(k-1), and in
// weight lambda + 1/2, T(p^2) is Shimura's operator
//     b(m) = a(p^2 m) + chi(p) ((-1)^lambda m / p) p^(lambda-1) a(m) + p^(2 lambda - 1) a(m / p^2).
// The normalisation (the constant, the source-in-column convention and the
// half-integral character) was fitted to Magma's AdjacencyMatrix on genera of ranks
// 2-9 and is checked by tests/test_theta_hecke.m (EichlerRelationHolds).
//
// So theta: Q^h -> M_k, e_i -> theta_i, intertwines alpha_p: v -> v Adj^t with
// T(p^2) + c_n(p), and its kernel K (the linear relations among the theta series) is
// alpha_p-stable.  Two consequences, used by ThetaHeckeDecide:
//   1. K lies in the largest alpha-stable subspace K* of the relations that hold to
//      the known precision (ThetaKernelCore).  A pair with e_i - e_j outside K* has
//      different theta series -- with no modular forms computation.
//   2. theta(K*) lies in W = {f cuspidal : F_p(T(p^2) + c_n(p)) f = 0}, F_p the minimal
//      polynomial of alpha_p on K*.  W is a small piece of the cusp space (the Hecke
//      eigenspaces matching the adjacency eigenvalues), usually determined by a
//      handful of coefficients; once the known coefficients determine W, theta(K*) = 0
//      and every pair with e_i - e_j in K* has equal theta series
//      (ThetaHeckeCertificate).

intrinsic EvenGramLevel(S::Mtrx) -> RngIntElt
{The level of the even symmetric matrix S: the smallest N with N S^-1 integral with even diagonal.}
    B := ChangeRing(S, Rationals())^-1;
    n := Nrows(B);
    return LCM([Integers() | Denominator(B[i, j]) : i, j in [1 .. n]] cat
               [Integers() | Denominator(B[i, i] / 2) : i in [1 .. n]]);
end intrinsic;

intrinsic ThetaFormData(L::Lat) -> BoolElt, RngIntElt, RngIntElt
{For the positive definite lattice L of rank n: halve, N, D.  The theta series used here
is g = sum_v q^((v,v)/2) when halve (L even) and the pipeline's f = sum_v q^(v,v)
otherwise; it lies in M_(n/2)(Gamma_0(N), chi) with chi(d) = (D / d).}
    A := GramMatrix(L);
    n := Nrows(A);
    halve := IsEven(L);
    S := halve select A else 2 * A;
    N := EvenGramLevel(S);
    det := Integers() ! Determinant(S);
    D := IsEven(n) select (-1)^(n div 2) * det else 2 * det;
    return halve, N, D;
end intrinsic;

intrinsic NeighbourHeckeShift(n::RngIntElt, p::RngIntElt) -> RngIntElt
{c_n(p) = p + p^2 + ... + p^(n-3) (-1 for n = 2): the p-neighbour operator on theta series
of rank n is T(p^2) + c_n(p).}
    return (p^(n - 2) - p) div (p - 1);
end intrinsic;

// The theta coefficients used here (see ThetaFormData) from the pipeline's: every
// other one for an even lattice.  Returns the first prec of them.
function form_coefficients(theta, halve, prec)
    coeff := func<m | m lt #theta select theta[m + 1] else 0>;
    return [Integers() | coeff(halve select 2 * m else m) : m in [0 .. prec - 1]];
end function;

// How many of the form's coefficients pipeline precision P (a_0 .. a_(P-1)) gives.
function form_precision(P, halve)
    return halve select (P + 1) div 2 else P;
end function;

// The Hecke image (T(p^2) + c_n(p)) a of the coefficient sequence a = [a_0, a_1, ...]
// of a form of weight n/2 and character (D / .), to precision prec; needs
// #a > p^2 (prec - 1).
function hecke_p2_image(a, n, D, p, prec)
    chip := KroneckerSymbol(D, p);
    c := NeighbourHeckeShift(n, p);
    if IsEven(n) then
        k := n div 2;
        // T(p^2) = T(p)^2 - chi(p) p^(k-1), with T(p): b(m) = a(p m) + chi(p) p^(k-1) a(m / p)
        tp := func<s, m | s[p * m + 1] + (m mod p eq 0 select chip * p^(k - 1) * s[m div p + 1] else 0)>;
        t1 := [tp(a, m) : m in [0 .. p * (prec - 1)]];
        return [tp(t1, m) - chip * p^(k - 1) * a[m + 1] + c * a[m + 1] : m in [0 .. prec - 1]];
    end if;
    lam := (n - 1) div 2;
    return [a[p^2 * m + 1] + chip * KroneckerSymbol((-1)^lam * m, p) * p^(lam - 1) * a[m + 1]
            + (m mod p^2 eq 0 select p^(2 * lam - 1) * a[m div p^2 + 1] else 0) + c * a[m + 1]
            : m in [0 .. prec - 1]];
end function;

intrinsic EichlerRelationHolds(reps::SeqEnum[Lat], Adj::Mtrx, p::RngIntElt, prec::RngIntElt) -> BoolElt
{Whether sum_j Adj[j, i] theta_j = (T(p^2) + c_n(p)) theta_i holds to precision prec for
every class i (reps a complete set of classes of a genus, Adj = AdjacencyMatrix in that
order, p a prime not dividing the level).  For testing the normalisation.}
    halve, N, D := ThetaFormData(reps[1]);
    require N mod p ne 0 : "p divides the level";
    n := Rank(reps[1]);
    depth := p^2 * (prec - 1) + 1;            // coefficients needed
    th := [form_coefficients(Eltseq(ThetaSeries(R, (halve select 2 * depth else depth) - 1)), halve, depth) : R in reps];
    for i in [1 .. #reps] do
        lhs := [&+[Adj[j, i] * th[j][m + 1] : j in [1 .. #reps]] : m in [0 .. prec - 1]];
        if lhs ne hecke_p2_image(th[i], n, D, p, prec) then return false; end if;
    end for;
    return true;
end intrinsic;

// A matrix whose columns are a basis of the column space of C (at least one column).
function column_basis(C)
    E := EchelonForm(Transpose(C));
    r := Rank(E);
    return r eq 0 select ZeroMatrix(BaseRing(C), Nrows(C), 1)
                    else Transpose(Submatrix(E, 1, 1, r, Ncols(E)));
end function;

intrinsic ThetaKernelCore(Theta::Mtrx, adjs::SeqEnum) -> ModTupFld, Mtrx
{Theta: h x M, row i the first M coefficients of theta_i; adjs: adjacency matrices
(AdjacencyMatrix's convention) at primes not dividing the level.  Returns K*, the largest
subspace of the relations (the v with v Theta = 0) that is stable under v -> v Adj^t for
every Adj in adjs, and a matrix C such that K* is the kernel of C.  The true relations
among the theta series lie in K*, so e_i - e_j outside it (rows i and j of C differ)
proves theta_i ne theta_j.}
    C := column_basis(ChangeRing(Theta, Rationals()));
    ops := [Transpose(ChangeRing(A, Rationals())) : A in adjs];
    repeat
        r := Ncols(C);
        for At in ops do
            C := column_basis(HorizontalJoin(C, At * C));
        end for;
    until Ncols(C) eq r;
    return Kernel(C), C;
end intrinsic;

// The smallest subspace of Q^h containing the vectors vs and stable under the maps
// v -> v Op for Op in ops.
function stable_span(vs, ops, h)
    V := VectorSpace(Rationals(), h);
    Z := sub<V | vs>;
    repeat
        d := Dimension(Z);
        Z := sub<V | Basis(Z) cat [b * Op : b in Basis(Z), Op in ops]>;
    until Dimension(Z) eq d;
    return Z;
end function;

// 1 + the largest pivot exponent of the echelon form of Q (rows: q-expansions
// a_0 .. a_(prec-1) of a basis of a space of dimension d), or -1 if Q has rank < d.
function determining_precision(Q, d)
    E := EchelonForm(Q);
    if Rank(E) lt d then return -1; end if;
    return Maximum([Minimum([m : m in [1 .. Ncols(E)] | E[i, m] ne 0]) : i in [1 .. d]]);
end function;

intrinsic ThetaHeckeSpace(n::RngIntElt, N::RngIntElt, D::RngIntElt, p::RngIntElt, prec::RngIntElt) -> Assoc
{The cusp space S_(n/2)(Gamma_0(N), (D / .)) (n >= 3; p a prime not dividing N), with what
ThetaHeckeEigenspace needs to cut out eigenspaces of the p-neighbour operator
T(p^2) + c_n(p) and expand them to precision prec.  Integral weight uses modular symbols;
half-integral weight (Magma has no Hecke operators there) a q-expansion basis known to
p^2 times the precision that determines the space, from which the operator is read off.}
    require n ge 3 : "rank at least 3 (no modular symbols in weight 1)";
    require N mod p ne 0 : "p divides the level";
    space := AssociativeArray();
    space["n"] := n; space["N"] := N; space["D"] := D; space["p"] := p; space["prec"] := prec;
    chi := DirichletGroup(N) ! KroneckerCharacter(D, Rationals());
    if IsEven(n) then
        space["modsym"] := CuspidalSubspace(ModularSymbols(chi, n div 2, +1));
        space["dim"] := Dimension(space["modsym"]);
        return space;
    end if;
    S := CuspidalSubspace(HalfIntegralWeightForms(chi, n / 2));
    e := Dimension(S);
    space["dim"] := e;
    if e eq 0 then return space; end if;
    X0 := PrecisionBound(S);
    X := Maximum(prec, p^2 * (X0 - 1) + 1);
    B := [];
    for b in Basis(S) do
        f := qExpansion(b, X);
        Append(~B, [Rationals() | Coefficient(f, m) : m in [0 .. X - 1]]);
    end for;
    B0 := Matrix(Rationals(), [b[1 .. X0] : b in B]);
    ok, Op := IsConsistent(B0, Matrix(Rationals(), [hecke_p2_image(b, n, D, p, X0) : b in B]));
    error if not ok, "T(p^2) does not preserve the half-integral weight cusp space (normalisation bug)";
    space["op"] := MatrixAlgebra(Rationals(), e) ! Op;
    space["basis"] := Matrix(Rationals(), [b[1 .. prec] : b in B]);
    return space;
end intrinsic;

intrinsic ThetaHeckeEigenbasis(space::Assoc, u::RngUPolElt) -> Mtrx, RngIntElt
{For an irreducible polynomial u over Q: the q-expansions a_0 .. a_(prec-1) (rows; prec
the space's) of a basis of W = ker u(T(p^2) + c_n(p)) in the cusp space of
ThetaHeckeSpace, and the dimension of W (the rows are dependent when prec does not
determine W).}
    n := space["n"]; p := space["p"]; prec := space["prec"];
    if space["dim"] eq 0 then return ZeroMatrix(Rationals(), 0, prec), 0; end if;
    R<x> := PolynomialRing(Rationals());
    if IsEven(n) then
        k := n div 2;
        G := Evaluate(R ! u, x^2 - KroneckerSymbol(space["D"], p) * p^(k - 1) + NeighbourHeckeShift(n, p));
        W := Kernel([<p, G>], space["modsym"]);
        d := Dimension(W);
        if d eq 0 then return ZeroMatrix(Rationals(), 0, prec), 0; end if;
        return Matrix(Rationals(), [[Coefficient(f, m) : m in [0 .. prec - 1]] : f in qExpansionBasis(W, prec)]), d;
    end if;
    Wc := Kernel(Evaluate(R ! u, space["op"]));
    d := Dimension(Wc);
    if d eq 0 then return ZeroMatrix(Rationals(), 0, prec), 0; end if;
    return BasisMatrix(Wc) * space["basis"], d;
end intrinsic;

intrinsic ThetaHeckeEigenspace(space::Assoc, u::RngUPolElt) -> RngIntElt, RngIntElt
{For an irreducible polynomial u over Q: the dimension of W = ker u(T(p^2) + c_n(p)) in the
cusp space of ThetaHeckeSpace, and the precision that determines W -- 1 + the largest
pivot exponent of an echelon basis of q-expansions, 0 if W = 0, -1 if more than the
space's prec.}
    Q, d := ThetaHeckeEigenbasis(space, u);
    if d eq 0 then return 0, 0; end if;
    return d, determining_precision(Q, d);
end intrinsic;

intrinsic ThetaHeckePrime(L::Lat) -> RngIntElt
{The smallest prime not dividing the level of the theta form of L (ThetaFormData): the
prime whose neighbours ThetaHeckeDecide uses.  It is 2 for an even lattice of odd level.}
    _, N, _ := ThetaFormData(L);
    p := 2;
    while N mod p eq 0 do p := NextPrime(p); end while;
    return p;
end intrinsic;

intrinsic ThetaHeckeOpen(thetas::SeqEnum, prec::RngIntElt, halve::BoolElt, adjs::SeqEnum) -> SeqEnum, SeqEnum
{The step that needs no modular forms.  thetas[i]: the pipeline theta coefficients
a_0, a_1, ... of class i of a genus (all known to precision at least prec); halve as in
ThetaFormData; adjs: adjacency matrices at primes not dividing the level, in the same
order, adjs[1] at the prime of ThetaHeckePrime.  Returns the pairs <i, j>, i < j, the
stable core (ThetaKernelCore) leaves open -- every other pair has different theta
series -- and the irreducible factors of the minimal polynomial of v -> v adjs[1]^t on the
relations the open pairs generate (ThetaHeckeEigenspace's u).}
    h := #thetas;
    M := form_precision(prec, halve);
    Theta := Matrix(Rationals(), [form_coefficients(t, halve, M) : t in thetas]);
    _, C := ThetaKernelCore(Theta, adjs);
    open := [<i, j> : i in [1 .. h], j in [1 .. h] | i lt j and C[i] eq C[j]];
    if #open eq 0 or #adjs eq 0 then return open, []; end if;
    V := VectorSpace(Rationals(), h);
    ops := [Transpose(ChangeRing(A, Rationals())) : A in adjs];
    BZ := BasisMatrix(stable_span([V.(t[1]) - V.(t[2]) : t in open], ops, h));
    F := MinimalPolynomial(Solution(BZ, BZ * ops[1]));
    return open, [t[1] : t in Factorization(F)];
end intrinsic;

intrinsic ThetaHeckeVerdict(open::SeqEnum, precisions::SeqEnum, prec::RngIntElt, halve::BoolElt) -> MonStgElt, RngIntElt
{From ThetaHeckeOpen's pairs (for theta series known to pipeline precision prec) and the
determining precisions (ThetaHeckeEigenspace, in the form's coefficients; -1 = more than
was tried) of the eigenspaces of its factors: the status ("distinguished", "equal" or
"unknown") and the pipeline precision the decision needs (0 when no eigenspace was
needed, -1 when unknown).  theta(relations) lies in the sum of the eigenspaces and
vanishes to the known precision, so once each eigenspace is determined by it the
relations hold exactly: the open pairs have equal theta series.}
    if #open eq 0 then return "distinguished", 0; end if;
    M := form_precision(prec, halve);
    if #precisions eq 0 or exists{P : P in precisions | P lt 0 or P gt M} then return "unknown", -1; end if;
    P := Maximum(precisions);
    return "equal", halve select Maximum(2 * P - 1, 0) else P;
end intrinsic;

intrinsic ThetaHeckeDecide(reps::SeqEnum[Lat], thetas::SeqEnum, prec::RngIntElt, adjacency::Assoc) -> MonStgElt, SeqEnum, RngIntElt
{Decide which classes of a genus have equal theta series, building the cusp space on the
spot (the pipeline's pass groups genera by space instead).  reps: the classes (rank >= 3,
positive definite); thetas[i]: the pipeline theta coefficients a_0, a_1, ... of reps[i],
known to precision at least prec; adjacency: p -> AdjacencyMatrix(G, p) in the order of
reps, which must include the prime of ThetaHeckePrime (others not dividing the level
refine the first step).  Returns the status ("distinguished", "equal" or "unknown"), the
pairs <i, j> proved to have equal theta series, and the pipeline precision the decision
needs (ThetaHeckeVerdict).}
    n := Rank(reps[1]);
    halve, N, D := ThetaFormData(reps[1]);
    p := ThetaHeckePrime(reps[1]);
    if not IsDefined(adjacency, p) then return "unknown", [], -1; end if;
    adjs := [adjacency[p]] cat [adjacency[q] : q in Keys(adjacency) | q ne p and N mod q ne 0];
    open, factors := ThetaHeckeOpen(thetas, prec, halve, adjs);
    if #open eq 0 then return "distinguished", [], 0; end if;
    space := ThetaHeckeSpace(n, N, D, p, form_precision(prec, halve));
    precisions := [];
    for u in factors do
        _, P := ThetaHeckeEigenspace(space, u);
        Append(~precisions, P);
    end for;
    status, need := ThetaHeckeVerdict(open, precisions, prec, halve);
    return status, status eq "equal" select open else [], need;
end intrinsic;

intrinsic ThetaHeckeAdjacency(reps::SeqEnum[Lat], p::RngIntElt) -> Mtrx
{AdjacencyMatrix(G, p) for the genus whose classes are reps, in the order of reps.}
    G := Genus(reps[1]);
    G`Representatives := reps;
    G`IsNatural := true;
    return AdjacencyMatrix(G, p);
end intrinsic;

// ---- The pass over the genera fill leaves undecided --------------------------------
//
// SetHashes writes a record for each such genus to theta_hecke/<rank>/<nplus>/<label>:
//   n|N|D|p|halve|h; the adjacency matrix at p (row-major, label order) or \N; then one
//   line per class: gram|theta series.
// ThetaHeckePass reads them, does the step without modular forms per genus, builds each
// cusp space once for all the eigenvalue factors its genera need (in a HardKill child),
// and writes theta_hecke_done/<rank>/<nplus>/<label>:
//   status|pairs proved equal|pipeline precision needed|p|seconds|note
// The eigenspace data is kept in theta_hecke_spaces/<n>/<N>.<D>.<p>, one line per factor u
// of the minimal polynomial of the p-neighbour operator: u's coefficients|dim|precision|
// precision tried (precision -1: not determined by that many coefficients), so a later
// run, or another genus of the space, need not build the space again.

function parse_ints(s)
    s := s[2 .. #s - 1];
    return s eq "" select [Integers() |] else [Integers() | StringToInteger(x) : x in Split(s, ",")];
end function;

function space_path(folder, n, N, D, p)
    return Sprintf("%o/%o", folder, n), Sprintf("%o/%o/%o.%o.%o", folder, n, N, D, p);
end function;

// cache: factor coefficient string -> <dim, precision, precision tried>
function read_space_cache(path)
    cache := AssociativeArray();
    ok, F := OpenTest(path, "r");
    if not ok then return cache; end if;
    for line in Split(Read(path), "\n") do
        f := Split(line, "|");
        cache[f[1]] := <StringToInteger(f[2]), StringToInteger(f[3]), StringToInteger(f[4])>;
    end for;
    return cache;
end function;

procedure write_space_cache(dir, path, cache)
    System("mkdir -p " * dir);
    tmp := path * ".tmp";
    Write(tmp, Join([Sprintf("%o|%o|%o|%o", k, cache[k][1], cache[k][2], cache[k][3]) : k in Keys(cache)], "\n") : Overwrite);
    System(Sprintf("mv %o %o", tmp, path));
end procedure;

function factor_key(u)
    return "{" * Join([Sprint(c) : c in Coefficients(u)], ",") * "}";
end function;

// in a child: the eigenspace data of the factors us of one space
function space_factors(n, N, D, p, prec, us)
    space := ThetaHeckeSpace(n, N, D, p, prec);
    R := PolynomialRing(Rationals());
    out := [];
    for u in us do
        d, P := ThetaHeckeEigenspace(space, R ! u);
        Append(~out, [d, P]);
    end for;
    return out;
end function;

intrinsic ThetaHeckePass(labels::SeqEnum[MonStgElt] : timeout := 600, Sidecar := "theta_hecke", Done := "theta_hecke_done", Spaces := "theta_hecke_spaces")
{Decide the genera with records in Sidecar (see SetHashes), writing a record per genus to
Done and keeping eigenspace data per cusp space in Spaces.  timeout bounds each
adjacency computation and each space (seconds, HardKill).  Genera sharing a space should
be in one call, which builds the space once.}
    R := PolynomialRing(Rationals());
    cols := GetColumns(); SetColumns(0);          // Write wraps long lines otherwise
    genera := [];
    for label in labels do
        t0 := Realtime();
        lines := Split(Read(LabelPath(Sidecar, label)), "\n");
        n, N, D, p, halve, h := Explode([StringToInteger(x) : x in Split(lines[1], "|")]);
        halve := halve eq 1;
        reps := []; thetas := [];
        for i in [1 .. h] do
            f := Split(lines[2 + i], "|");
            Append(~reps, LatticeWithGram(Matrix(Integers(), n, n, parse_ints(f[1]))));
            Append(~thetas, parse_ints(f[2]));
        end for;
        prec := Minimum([#t : t in thetas]);
        g := AssociativeArray();
        g["label"] := label; g["key"] := <n, N, D, p>; g["halve"] := halve; g["prec"] := prec;
        g["M"] := form_precision(prec, halve); g["seconds"] := 0.0; g["note"] := "-";
        A := 0;
        if lines[2] ne "\\N" then
            A := Matrix(Integers(), h, h, parse_ints(lines[2]));
        else
            ok, res := TimeoutCall(timeout, ThetaHeckeAdjacency, <reps, p>, 1 : HardKill := true);
            if ok then A := res[1]; else g["note"] := "adjacency_timeout"; end if;
        end if;
        if Type(A) eq RngIntElt then
            g["open"] := [<i, j> : i, j in [1 .. h] | i lt j]; g["factors"] := []; g["status"] := "unknown";
        else
            g["open"], g["factors"] := ThetaHeckeOpen(thetas, prec, halve, [A]);
            g["status"] := #g["open"] eq 0 select "distinguished" else "pending";
        end if;
        g["seconds"] := Realtime(t0);
        Append(~genera, g);
    end for;
    // each space once, for every factor its pending genera need
    keys := {g["key"] : g in genera | g["status"] eq "pending"};
    caches := AssociativeArray();
    for key in keys do
        n, N, D, p := Explode(key);
        dir, path := space_path(Spaces, n, N, D, p);
        cache := read_space_cache(path);
        mine := [g : g in genera | g["status"] eq "pending" and g["key"] eq key];
        M := Maximum([g["M"] : g in mine]);
        need := {};
        for g in mine do
            for u in g["factors"] do
                k := factor_key(u);
                if not IsDefined(cache, k) or (cache[k][2] lt 0 and cache[k][3] lt g["M"]) then
                    Include(~need, k);
                end if;
            end for;
        end for;
        if #need gt 0 then
            us := [[Rationals() | StringToInteger(c) : c in Split(k[2 .. #k - 1], ",")] : k in need];
            t0 := Realtime();
            ok, res := TimeoutCall(timeout, space_factors, <n, N, D, p, M, us>, 1 : HardKill := true);
            if ok then
                for i in [1 .. #us] do
                    cache[factor_key(R ! us[i])] := <res[1][i][1], res[1][i][2], M>;
                end for;
                write_space_cache(dir, path, cache);
            end if;
            // the space's time is charged to its genera evenly
            for i in [1 .. #genera] do
                if genera[i]["status"] eq "pending" and genera[i]["key"] eq key then
                    genera[i]["seconds"] +:= Realtime(t0) / #mine;
                    if not ok then genera[i]["note"] := "space_timeout"; end if;
                end if;
            end for;
        end if;
        caches[key] := cache;
    end for;
    for g in genera do
        status := g["status"]; pairs := []; needp := -1;
        if status eq "distinguished" then
            needp := 0;
        elif status eq "pending" then
            cache := caches[g["key"]];
            precisions := [];
            for u in g["factors"] do
                k := factor_key(u);
                Append(~precisions, IsDefined(cache, k) select cache[k][2] else -1);
            end for;
            status, needp := ThetaHeckeVerdict(g["open"], precisions, g["prec"], g["halve"]);
            if status eq "equal" then pairs := g["open"]; end if;
        end if;
        Write(LabelPath(Done, g["label"] : Create),
              Sprintf("%o|%o|%o|%o|%.3o|%o", status, to_postgres([[t[1], t[2]] : t in pairs]), needp, g["key"][4], g["seconds"], g["note"]) : Overwrite);
    end for;
    SetColumns(cols);
end intrinsic;
