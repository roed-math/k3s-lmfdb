// IsAdditivelyIndecomposable and the Plesken III.1 helper from connect_genus.m.
// The intrinsic returns <is_indecomposable, is_determined>.
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

// Plesken's families a;(s)b = diag(1,...,1,-s) + a^t a.
function pleskenLat(av, s)
    m := #av;
    G := DiagonalMatrix(Integers(), [1 : i in [1..m-1]] cat [-s])
         + Transpose(Matrix(Integers(),1,m,av)) * Matrix(Integers(),1,m,av);
    return LatticeWithGram(G);
end function;

// Complete classification in low rank (Mordell; Plesken Thm. III.4 -> E_n).
for tup in [* <StandardLattice(1), true>, <Lattice("A",2), false>, <Lattice("D",5), false>,
              <Lattice("E",6), true>, <Lattice("A",6), false>, <Lattice("E",7), true>,
              <Lattice("E",8), true>, <Lattice("D",8), false>, <StandardLattice(9), false> *] do
    v, known := IsAdditivelyIndecomposable(tup[1]);
    AssertTrue(~results, known, Sprintf("IsAdditivelyIndecomposable(%o): determined", tup[1]));
    AssertEqual(~results, v, tup[2], Sprintf("IsAdditivelyIndecomposable(%o): value", tup[1]));
end for;

// Discriminant 4 and 5 (Wang Thm. 2.16): Plesken's indecomposable examples.
for L in [ pleskenLat([1: i in [1..11]] cat [4], 1),     // 1^11;4, det 4
           pleskenLat([1: i in [1..9]] cat [5], 2) ] do  // 1^9;(2)5, det 5
    v, known := IsAdditivelyIndecomposable(L);
    AssertTrue(~results, known, "Plesken family example: determined");
    AssertTrue(~results, v, "Plesken family example: indecomposable");
end for;

// A det-4 lattice that IS additively decomposable: complement of (x,y) with
// Q(x)=Q(y)=2 in E8 perp E8 (decomposes as E8 + E8).
E8 := Lattice("E",8);  G8 := ChangeRing(GramMatrix(E8), Integers());
wc := Vector(Integers(),
        Coordinates(E8, ShortestVectors(E8)[1]) cat Coordinates(E8, ShortestVectors(E8)[2]));
G16 := DiagonalJoin(G8, G8);
Bc := Matrix(Integers(), [ Eltseq(b) : b in Basis(Kernel(G16 * Transpose(Matrix(wc)))) ]);
Ldec := LatticeWithGram(Bc * G16 * Transpose(Bc));
v, known := IsAdditivelyIndecomposable(Ldec);
AssertTrue(~results, known, "E8+E8 complement: determined");
AssertTrue(~results, not v, "E8+E8 complement: decomposable");

// Plesken III.1 sufficient condition, standalone.
AssertTrue(~results, SatisfiesPleskenIII1(Lattice("E",8)), "SatisfiesPleskenIII1(E8)");
AssertTrue(~results, not SatisfiesPleskenIII1(StandardLattice(3)),
    "not SatisfiesPleskenIII1(Z^3): dual minimum <= 1");

Report(~results, "test_additively_indecomposable");
