AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

lats := [LatticeWithGram(Matrix(11,11,i)): i in [
    [ 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 2, 1, -1, 0, 1, 1, 0, -1, 1, 2, 0, 1,
    2, -1, 0, 1, 1, 0, 0, 1, 2, 0, -1, -1, 2, -1, -1, -1, 0, 0, 0, -2, 0, 0, 0,
    -1, 4, 2, 1, 3, 4, 0, 1, 0, 1, 1, -1, 2, 4, 1, 2, 2, 2, 3, 0, 1, 1, -1, 1,
    1, 2, 2, 1, 1, 2, 0, 0, 0, 0, 3, 2, 2, 6, 4, 1, 1, 0, -1, 0, 0, 4, 2, 1, 4,
    6, 0, 0, -1, 1, 1, 0, 0, 2, 1, 1, 0, 4, 2, 0, 2, 2, -2, 1, 3, 2, 1, 0, 2, 5
    ],
    [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -1, -1, 0, -1, 1, -1, 2, 2, 0,
    1, 2, -1, -1, 0, -1, 1, -1, 2, 1, 0, -1, -1, 2, 0, 0, 1, 0, 1, -1, -1, 0,
    -1, -1, 0, 2, 0, 0, -1, 0, -1, -1, 0, 0, 0, 0, 0, 2, -1, 0, -1, 0, 0, 0, -1,
    -1, 1, 0, -1, 3, -1, 2, -2, -2, 0, 1, 1, 0, -1, 0, -1, 3, 0, 3, 3, 0, -1,
    -1, 1, 0, -1, 2, 0, 3, -1, 0, 0, 2, 2, -1, -1, 0, -2, 3, -1, 5, 4, 0, 2, 1,
    -1, -1, 0, -2, 3, 0, 4, 6 ]
]];

AssertTrue(~results, IsIsomorphic(lats[1], lats[2]), "IsIsomorphic(lats[1], lats[2])");
AssertTrue(~results, CanonicalForm(lats[1]) eq CanonicalForm(lats[2]),
    "CanonicalForm(lats[1]) eq CanonicalForm(lats[2])");
AssertTrue(~results, CanonicalForm(GramMatrix(lats[1])) eq CanonicalForm(GramMatrix(lats[2])),
    "CanonicalForm(GramMatrix(lats[1])) eq CanonicalForm(GramMatrix(lats[2]))");

Report(~results, "test_bug2");
