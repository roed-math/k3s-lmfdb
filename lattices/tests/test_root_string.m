// Root systems (orth_dec.mag): RootString's canonical form, ParseRootString,
// RootCount, and CompareRootStrings, the root_sublattice key of cmp_lat in
// fill_genus.m.
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

function DS(grams)
    return LatticeWithGram(DiagonalJoin(<Matrix(Integers(), G) : G in grams>));
end function;

// A lattice in a random basis, so the component order orth_dec finds is arbitrary.
function scrambled(L)
    n := Rank(L);
    repeat
        U := RandomSLnZ(n, 4*n, 3*n);
    until Abs(Determinant(U)) eq 1;
    return LatticeWithGram(U * GramMatrix(L) * Transpose(U));
end function;

E8 := CartanMatrix("E8");
A2 := CartanMatrix("A2");

// --- canonical form ---------------------------------------------------------
L := DS([* E8, IdentityMatrix(Integers(), 4), A2, E8 *]);
for i in [1..3] do
    AssertEqual(~results, RootString(scrambled(L)), "I4 A2 2E8",
        "I4 + A2 + E8^2 in a random basis: components in I, A, D, E order");
end for;
AssertEqual(~results, RootString(DS([* CartanMatrix("D5"), CartanMatrix("A3"), CartanMatrix("A1"), CartanMatrix("A3") *])),
    "A1 2A3 D5", "multiplicities, and increasing rank within a type");
AssertEqual(~results, RootString(DS([* CartanMatrix("E6"), CartanMatrix("E7") *])), "E6 E7", "E6 E7");
AssertEqual(~results, RootString(StandardLattice(1)), "I1", "Z is I1");
AssertEqual(~results, RootString(LatticeWithGram(Matrix(Integers(), 2, 2, [4, 2, 2, 4])) : IsRootSystem := false),
    "0", "no vectors of norm 1 or 2: \"0\", not the empty string");
// [[2,1],[1,4]]: only +-e1 has norm <= 2
AssertEqual(~results, RootString(LatticeWithGram(Matrix(Integers(), 2, 2, [2, 1, 1, 4])) : IsRootSystem := false),
    "A1", "root sublattice of a non-root lattice");
// I2 + <3>: the norm-2 vectors e1 +- e2 lie in I2
AssertEqual(~results, RootString(LatticeWithGram(DiagonalMatrix(Integers(), [1, 3, 1])) : IsRootSystem := false),
    "I2", "norm-2 vectors inside the I_k part");
// E8 contains E7 + A1, but its root sublattice is E8
AssertEqual(~results, RootString(scrambled(LatticeWithGram(E8)) : IsRootSystem := false), "E8", "E8");

// --- ParseRootString / RootCount ------------------------------------------------
AssertEqual(~results, ParseRootString("I4 A2 2E8"), [<"I", 4>, <"A", 2>, <"E", 8>, <"E", 8>], "parse");
AssertEqual(~results, ParseRootString("0"), [], "parse the empty root system");
for s in ["X3", "2", "A", "A0", "0A2", "A2  E8", "E8 A2", "E8 E8", "1A2", "", "I2 I1"] do
    ok := true;
    try
        _ := ParseRootString(s);
    catch e
        ok := false;
    end try;
    AssertTrue(~results, not ok, Sprintf("ParseRootString rejects %o", s));
end for;
// RootCount counts vectors of norm 1 and 2 (ShortVectors lists them up to sign)
for M in [DS([* E8, IdentityMatrix(Integers(), 4), A2, E8 *]), DS([* CartanMatrix("D5"), CartanMatrix("A3") *]),
          DS([* CartanMatrix("E6"), CartanMatrix("E7"), IdentityMatrix(Integers(), 3) *]),
          LatticeWithGram(Matrix(Integers(), 3, 3, [2, 1, 0, 1, 2, 1, 0, 1, 5]))] do
    L := M;
    s := RootString(L : IsRootSystem := false);
    AssertEqual(~results, RootCount(s), 2 * #ShortVectors(L, 1, 2), Sprintf("RootCount(%o)", s));
end for;
AssertEqual(~results, RootCount("0"), 0, "RootCount of the empty root system");

// --- CompareRootStrings -----------------------------------------------------------
AssertTrue(~results, CompareRootStrings("A1", "0") lt 0, "any roots sort before none");
AssertTrue(~results, CompareRootStrings("0", "A1") gt 0, "antisymmetric");
AssertTrue(~results, CompareRootStrings("E8", "D8") lt 0, "more roots first (240 > 112)");
AssertTrue(~results, CompareRootStrings("2A2", "A1 A3") gt 0, "more roots first (12 < 14)");
AssertTrue(~results, CompareRootStrings("6A1", "A3") lt 0, "equal roots (12): larger rank first (6 > 3)");
AssertTrue(~results, CompareRootStrings("A1 2A3", "A1 2A3") eq 0, "a root system equals itself");
// same number of roots (480) and rank (16): the larger component first
AssertTrue(~results, CompareRootStrings("D16", "2E8") lt 0, "D16 before 2E8");
// E6 and I6 both have rank 6 and 72 vectors of norm <= 2: type breaks the tie
AssertTrue(~results, CompareRootStrings("E6", "I6") lt 0, "E6 before I6");
// a total order: sorting is independent of the input order and agrees pairwise
pool := ["0", "A1", "2A1", "A2", "I1", "I2", "A3", "3A1", "D4", "4A1", "I4", "A1 A3", "2A2",
         "E6", "I6", "A6", "D6", "E8", "D8", "2E8", "D16", "I4 A2 2E8", "A1 2A3 D5", "6A1"];
cmp := func<a, b | CompareRootStrings(a, b)>;   // Sort takes a function, not an intrinsic
srt := Sort(pool, cmp);
AssertTrue(~results, forall{i : i in [1..#srt-1] | CompareRootStrings(srt[i], srt[i+1]) lt 0},
    "the pool sorts strictly (distinct root systems never compare equal)");
AssertTrue(~results, forall{<a, b> : a, b in pool |
    Sign(CompareRootStrings(a, b)) eq -Sign(CompareRootStrings(b, a))}, "antisymmetric on the pool");
AssertEqual(~results, Sort(Reverse(pool), cmp), srt, "sorting does not depend on the input order");

Report(~results, "test_root_string");
