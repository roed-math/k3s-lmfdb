// Minimal-vector lattice invariants from connect_genus.m:
// IsWellRounded, IsStronglyWellRounded, IsMinimalVectorGenerated, IsEutactic,
// PerfectionDefect.
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

// --- IsStronglyWellRounded -------------------------------------------------
for tup in [* <Lattice("A",2), true>, <Lattice("A",3), true>, <Lattice("D",4), true>,
              <Lattice("D",5), true>, <Lattice("E",6), true>, <Lattice("E",8), true>,
              <StandardLattice(3), true>,
              <LatticeWithGram(DiagonalMatrix(Rationals(),[1,2,3])), false> *] do
    L := tup[1];
    AssertEqual(~results, IsStronglyWellRounded(L, ShortestVectors(L)), tup[2],
        Sprintf("IsStronglyWellRounded(%o)", L));
end for;
// rank-6 lattice that is well rounded but NOT strongly well rounded
G6 := SymmetricMatrix(Rationals(), [55, 1,19, 36,0,55, -1,-3,-37,46, 7,3,7,0,6, -20,6,-35,29,1,31]);
L6 := LatticeWithGram(G6);  S6 := ShortestVectors(L6);
AssertTrue(~results, IsWellRounded(L6, S6), "rank-6 example: well rounded");
AssertTrue(~results, not IsStronglyWellRounded(L6, S6), "rank-6 example: not strongly well rounded");

// --- IsEutactic (with exact certificate) -----------------------------------
for L in [ StandardLattice(3), Lattice("A",2), Lattice("A",3), Lattice("D",4),
           Lattice("E",6), Lattice("E",8) ] do
    S := ShortestVectors(L);
    eu, coeffs := IsEutactic(L, S);
    AssertTrue(~results, eu, Sprintf("IsEutactic(%o)", L));
    AssertTrue(~results, #coeffs eq #S and forall{ c : c in coeffs | c gt 0 },
        Sprintf("IsEutactic(%o): positive coefficients, one per shortest vector", L));
    G := ChangeRing(GramMatrix(L), Rationals());
    U := [ Vector(Rationals(), [ x : x in Coordinates(L, s) ]) : s in S ];
    M := &+[ coeffs[i] * (Transpose(Matrix(U[i])) * Matrix(U[i])) : i in [1..#U] ];
    AssertEqual(~results, M, G^-1,
        Sprintf("IsEutactic(%o): sum_s c_s u_s^t u_s = G^-1", L));
end for;
// well rounded but not (weakly) eutactic, and a non-well-rounded lattice
Lne := LatticeWithGram(Matrix(Rationals(),2,2,[2,1/2,1/2,2]));
AssertTrue(~results, not IsEutactic(Lne, ShortestVectors(Lne)),
    "well-rounded but not eutactic example");
Lnw := LatticeWithGram(DiagonalMatrix(Rationals(),[1,2,3]));
AssertTrue(~results, not IsEutactic(Lnw, ShortestVectors(Lnw)),
    "non-well-rounded example: not eutactic");

// --- PerfectionDefect / IsMinimalVectorGenerated ---------------------------
AssertEqual(~results, PerfectionDefect(Lattice("E",8), ShortestVectors(Lattice("E",8))), 0,
    "PerfectionDefect(E8): perfect");
AssertEqual(~results, PerfectionDefect(Lattice("A",2), ShortestVectors(Lattice("A",2))), 0,
    "PerfectionDefect(A2): perfect");
AssertEqual(~results, PerfectionDefect(Lattice("A",3), ShortestVectors(Lattice("A",3))), 0,
    "PerfectionDefect(A3): perfect");
AssertTrue(~results, IsMinimalVectorGenerated(Lattice("E",8), ShortestVectors(Lattice("E",8))),
    "IsMinimalVectorGenerated(E8)");
AssertTrue(~results, IsMinimalVectorGenerated(Lattice("A",3), ShortestVectors(Lattice("A",3))),
    "IsMinimalVectorGenerated(A3)");

Report(~results, "test_shell_invariants");
