// ThetaSturmPrecision (hash.m): the precision at which the pipeline's theta series
// sum_v q^((v,v)) determine a modular form in their space, and SetHashes' use of it.
AttachSpec("lattices.spec");
load "tests/assertions.m";
procedure()   // SetHashes forks (TimeoutCall)
results := NewResults();
Lat := func<g | LatticeWithGram(Matrix(Rationals(), Isqrt(#g), Isqrt(#g), g))>;

// E8: theta(q) = E_4(2 tau) in M_4(Gamma_0(2)): N = 2, m = 3, P = Floor(4*3/12) + 1 = 2.
AssertEqual(~results, ThetaSturmPrecision(Lattice("E", 8)), 2, "E8");
// Z: sum q^(n^2) in M_{1/2}(Gamma_0(4)): N = 4, m = 6, P = Floor(6/24) + 1 = 1.
AssertEqual(~results, ThetaSturmPrecision(Lat([1])), 1, "Z");
// Z^2: theta^2 in M_1(Gamma_0(4)): P = Floor(1*6/12) + 1 = 1.
AssertEqual(~results, ThetaSturmPrecision(Lat([1,0,0,1])), 1, "Z^2");
// D4 (even, det 4): Theta_std in M_2(Gamma_0(2)), so theta(q) = Theta_std(2 tau) has
// level 4: m = 6, P = Floor(2*6/12) + 1 = 2.
AssertEqual(~results, ThetaSturmPrecision(Lattice("D", 4)), 2, "D4");

// The bound holds: for the even rank cases, Magma's own precision bound for the space
// M_k(Gamma_0(N), chi), chi(d) = ((-1)^k det(2A) / d), is no larger, and the theta
// series lies in that space (so N and the weight are right).
for t in [<Lattice("D", 4), 4>, <Lat([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,3]), 12>, <Lattice("A", 2), 6>] do
    L, N := Explode(t);
    k := Rank(L) div 2;
    chi := DirichletGroup(N) ! KroneckerCharacter((-1)^k * Integers() ! Determinant(2 * GramMatrix(L)), Rationals());
    M := ModularForms(chi, k);
    P := ThetaSturmPrecision(L);
    AssertTrue(~results, PrecisionBound(M) le P, Sprintf("Magma's precision bound %o <= ours %o for rank %o level %o", PrecisionBound(M), P, Rank(L), N));
    th := ThetaSeries(L, P + 20);
    B := [qExpansion(f, P + 21) : f in Basis(M)];
    rows := [[Coefficient(b, i) : i in [0 .. P + 20]] : b in B];
    S := RowSpace(Matrix(Rationals(), rows));
    AssertTrue(~results, Vector(Rationals(), [Coefficient(th, i) : i in [0 .. P + 20]]) in S,
        Sprintf("theta series of the rank-%o lattice lies in M_%o(Gamma_0(%o), chi)", Rank(L), k, N));
end for;

// SetHashes: two classes with equal theta series are certified (false) once they
// agree to the bound -- here E8+E8 and D16+, theta = E_4(2 tau)^2, P = 3 -- and a
// pair that differs is distinguished (true).
function rec(L, theta, prec)
    r := AssociativeArray();
    r["lattice"] := L; r["theta_series"] := theta; r["theta_prec"] := prec;
    r["rank"] := IntegerToString(Rank(L)); r["nplus"] := IntegerToString(Rank(L));
    r["level"] := "1"; r["disc_abs"] := "1"; r["is_even"] := "T"; r["minimum"] := Minimum(L);
    return r;
end function;
E := DirectSum(Lattice("E", 8), Lattice("E", 8));
// D16+ is the other class of the genus of E8+E8 (even unimodular, rank 16)
D := [R : R in GenusRepresentatives(E) | not IsIsometric(R, E)][1];
th := Eltseq(ThetaSeries(E, 9));
lats := [rec(E, th, 10), rec(D, Eltseq(ThetaSeries(D, 9)), 10)];
G := AssociativeArray(); te := AssociativeArray(); te[10] := 0.01;
SetHashes(~lats, ~G, te, 60);
AssertEqual(~results, G["is_theta_distinguished"], false, "E8+E8 vs D16+: equal theta series certified");
lats := [rec(E, th, 10), rec(Lattice("A", 16), Eltseq(ThetaSeries(Lattice("A", 16), 9)), 10)];
SetHashes(~lats, ~G, te, 60);
AssertEqual(~results, G["is_theta_distinguished"], true, "E8+E8 vs A16: distinguished");
Report(~results, "test_theta_bound");
end procedure();
