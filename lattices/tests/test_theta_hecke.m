// theta_hecke.m: Eichler's commutation relation between p-neighbours and Hecke
// operators, and the decision of which classes of a genus have equal theta series.
AttachSpec("lattices.spec");
load "tests/assertions.m";
procedure()   // GenusRepresentatives / TimeoutCall may fork
results := NewResults();
Lat := func<g | LatticeWithGram(Matrix(Integers(), Isqrt(#g), Isqrt(#g), g))>;
Diag := func<v | LatticeWithGram(DiagonalMatrix(Integers(), v))>;

// c_n(p) = p + ... + p^(n-3)
AssertEqual(~results, [NeighbourHeckeShift(n, 3) : n in [2 .. 6]], [-1, 0, 3, 12, 39], "c_n(3), n = 2..6");
AssertEqual(~results, NeighbourHeckeShift(8, 5), 3905, "c_8(5)");

// The form used: g = sum q^((v,v)/2) of level N(A) for even lattices, the pipeline's
// f = sum q^(v,v) of level N(2A) otherwise; chi = (D / .).
for t in [<Lattice("E", 8), <true, 1, 1>, "E8: even, level 1, trivial character">,
          <Lattice("D", 4), <true, 2, 4>, "D4: even, level 2, D = det">,
          <Diag([1, 1, 1]), <false, 4, 16>, "Z^3: odd, level 4, D = 2 det(2A)">] do
    halve, N, D := ThetaFormData(t[1]);
    AssertEqual(~results, <halve, N, D>, t[2], t[3]);
end for;
AssertEqual(~results, ThetaHeckePrime(DirectSum(Lattice("A", 4), Lattice("A", 4))), 2, "A4+A4: odd level, so p = 2");
AssertEqual(~results, ThetaHeckePrime(Diag([1, 1, 1, 1, 1, 1, 5])), 3, "I6+<5>: level 20, so p = 3");

// Eichler's relation, sum_j Adj[j, i] theta_j = (T(p^2) + c_n(p)) theta_i, to precision
// 8, on genera with distinct theta series: odd ranks 3 and 5 (half-integral weight),
// rank 4 with chi(p) = +1 and -1, and even lattices.
cases := [* <"3.3.10.af", ["{1,0,0,0,1,0,0,0,10}", "{2,0,0,0,2,1,0,1,3}"], 3>,
            <"4.4.10.63", ["{1,0,0,0,0,1,0,0,0,0,2,0,0,0,0,5}", "{1,0,0,0,0,4,3,3,0,3,4,3,0,3,3,4}"], 3>,
            <"4.4.11.1d", ["{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,11}", "{1,0,0,0,0,1,0,0,0,0,3,1,0,0,1,4}"], 3>,
            <"5.5.14.df", ["{1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,3,1,0,0,0,1,5}", "{1,0,0,0,0,0,2,0,0,0,0,0,3,2,2,0,0,2,3,2,0,0,2,2,3}"], 3> *];
for c in cases do
    reps := [Lat([StringToInteger(x) : x in Split(g[2 .. #g - 1], ",")]) : g in c[2]];
    A := ThetaHeckeAdjacency(reps, c[3]);
    AssertTrue(~results, EichlerRelationHolds(reps, A, c[3], 8), Sprintf("%o: Eichler relation at p = %o", c[1], c[3]));
    if c[1] eq "4.4.10.63" then
        // the matrix is not symmetric here, so this checks the source-in-column convention
        AssertTrue(~results, A ne Transpose(A) and not EichlerRelationHolds(reps, Transpose(A), c[3], 8),
                   "4.4.10.63: the transposed adjacency matrix fails the relation");
    end if;
end for;
// even lattices, through g = sum q^((v,v)/2): odd level at p = 2, even level, odd rank
for t in [<"A4+A4", [<"A", 4>, <"A", 4>], 2>, <"A3+A1+A1+A1", [<"A", 3>, <"A", 1>, <"A", 1>, <"A", 1>], 3>,
          <"A1+A1+A1+A2", [<"A", 1>, <"A", 1>, <"A", 1>, <"A", 2>], 5>] do
    reps := GenusRepresentatives(DirectSum([Lattice(x[1], x[2]) : x in t[2]]));
    AssertTrue(~results, EichlerRelationHolds(reps, ThetaHeckeAdjacency(reps, t[3]), t[3], 8), Sprintf("%o: Eichler relation at p = %o", t[1], t[3]));
end for;

// The stable core: with theta_1 = theta_2 to the known precision but an adjacency
// operator that moves e_1 - e_2 out of the relations, the pair is distinguished.
Th := Matrix(Rationals(), [[1, 0], [1, 0], [1, 1]]);
K, C := ThetaKernelCore(Th, [Matrix(Integers(), [[1, 0, 0], [0, 1, 0], [1, 0, 1]])]);   // (e1 - e2) A^t = (1, -1, 1)
AssertEqual(~results, Dimension(K), 0, "stable core of a non-stable relation is 0");
K, C := ThetaKernelCore(Th, [Matrix(Integers(), [[2, 1, 0], [1, 2, 0], [0, 0, 3]])]);
AssertEqual(~results, Dimension(K), 1, "a stable relation stays");

// Equal theta series certified: E8+E8 and D16+ (S_8(SL_2(Z)) = 0, nothing to check),
// and the two classes of the genus of A2^4.
for t in [<"E8+E8", DirectSum(Lattice("E", 8), Lattice("E", 8)), 10>, <"A2^4", DirectSum([Lattice("A", 2) : i in [1 .. 4]]), 100>] do
    reps := GenusRepresentatives(t[2]);
    p := ThetaHeckePrime(reps[1]);
    adjacency := AssociativeArray(); adjacency[p] := ThetaHeckeAdjacency(reps, p);
    th := [Eltseq(ThetaSeries(R, t[3] - 1)) : R in reps];
    status, pairs := ThetaHeckeDecide(reps, th, t[3], adjacency);
    AssertEqual(~results, <status, pairs>, <"equal", [<1, 2>]>, t[1] * ": equal theta series certified");
end for;

// Never "equal" for theta series that differ, at any truncation up to their first
// difference; "distinguished" once it is included.
for c in cases do
    reps := [Lat([StringToInteger(x) : x in Split(g[2 .. #g - 1], ",")]) : g in c[2]];
    p := ThetaHeckePrime(reps[1]);
    adjacency := AssociativeArray(); adjacency[p] := ThetaHeckeAdjacency(reps, p);
    th := [Eltseq(ThetaSeries(R, 49)) : R in reps];
    m := FirstDiff(th[1], th[2], 50);
    st := [ThetaHeckeDecide(reps, th, M, adjacency) : M in [1 .. m]];
    AssertTrue(~results, "equal" notin st, Sprintf("%o: not equal at any truncation <= %o", c[1], m));
    AssertEqual(~results, ThetaHeckeDecide(reps, th, m + 1, adjacency), "distinguished", c[1] * ": distinguished");
end for;

// The pass end to end: SetHashes told not to extend the theta series leaves the
// genus undecided and writes its record; ThetaHeckePass decides it.  A2^4's classes
// agree to precision 3 < Sturm 5 and are equal; 4.4.11.1d's agree to their first
// difference and must stay unknown.
tmp := Sprintf("/tmp/test_theta_hecke_%o", Getpid());
function rec(L, theta, prec, glabel)
    r := AssociativeArray();
    r["lattice"] := L; r["theta_series"] := theta; r["theta_prec"] := prec;
    r["rank"] := IntegerToString(Rank(L)); r["nplus"] := IntegerToString(Rank(L));
    r["minimum"] := Minimum(L); r["genus_label"] := glabel;
    return r;
end function;
reps42 := [Lat([StringToInteger(x) : x in Split(g[2 .. #g - 1], ",")]) : g in cases[3][2]];
th42 := [Eltseq(ThetaSeries(R, 49)) : R in reps42];
m42 := FirstDiff(th42[1], th42[2], 50);
for t in [<"8.8.81.test", GenusRepresentatives(DirectSum([Lattice("A", 2) : i in [1 .. 4]])), 3, "equal">,
          <"4.4.11.test", reps42, m42, "unknown">] do
    lats := [rec(R, Eltseq(ThetaSeries(R, t[3] - 1)), t[3], t[1]) : R in t[2]];
    G := AssociativeArray(); te := AssociativeArray(); te[t[3]] := 0.01;
    SetHashes(~lats, ~G, te, 60 : Sidecar := tmp * "/theta_hecke", Extend := false);
    AssertEqual(~results, G["is_theta_distinguished"], "\\N", t[1] * ": undecided by fill");
    ThetaHeckePass([t[1]] : timeout := 300, Sidecar := tmp * "/theta_hecke", Done := tmp * "/done", Spaces := tmp * "/spaces");
    out := Split(Read(LabelPath(tmp * "/done", t[1])), "|");
    AssertEqual(~results, out[1], t[4], t[1] * ": pass verdict");
end for;
// the space cache: a rerun reads it instead of building the space
cached := Read(Sprintf("%o/spaces/8/3.81.2", tmp));
AssertTrue(~results, #Split(cached, "\n") ge 1, "space cache written");
System("rm -rf " * tmp);

Report(~results, "test_theta_hecke");
end procedure();
