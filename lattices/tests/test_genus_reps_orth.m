// Orbit-method genus enumeration (genus_reps.m).  Class numbers of the genera
// G_{n,p} of even lattices of determinant p (n even) or 2p (n odd) are known
// (Chenevier-Taibi, Tables 1.2/1.3); seeds from their Table 7.2, padded by E8.
AttachSpec("lattices.spec");
load "tests/assertions.m";
results := NewResults();

E8 := LatticeWithGram(CartanMatrix("E8"));
E6 := LatticeWithGram(CartanMatrix("E6"));
A1 := LatticeWithGram(CartanMatrix("A1"));
A2 := LatticeWithGram(CartanMatrix("A2"));
A4 := LatticeWithGram(CartanMatrix("A4"));

// Single-class genera (mass = 1/|Aut| short-circuit).
AssertEqual(~results, #GenusRepresentativesOrth(E8), 1, "G(E8): class number");
AssertEqual(~results, #GenusRepresentativesOrth(A2), 1, "G(A2): class number");
AssertEqual(~results, #GenusRepresentativesOrth(DirectSum(A1, A2)), 1, "G_{3,3}: class number");

// G_{11,3} (det 6): class number 2.
L := DirectSum(DirectSum(A1, A2), E8);
reps := GenusRepresentativesOrth(L);
AssertEqual(~results, #reps, 2, "G_{11,3}: class number");
AssertTrue(~results, forall{ R : R in reps | Genus(R) eq Genus(L) }, "G_{11,3}: reps lie in Genus(L)");
AssertEqual(~results, &+[ 1/#AutomorphismGroup(R) : R in reps ], Mass(L), "G_{11,3}: mass formula");

// G_{12,5} (det 5): class number 2.
L := DirectSum(A4, E8);
reps := GenusRepresentativesOrth(L);
AssertEqual(~results, #reps, 2, "G_{12,5}: class number");
AssertTrue(~results, forall{ R : R in reps | Genus(R) eq Genus(L) }, "G_{12,5}: reps lie in Genus(L)");

// G_{14,3} (det 3): class number 2.  The parent chain passes through unimodular
// lattices of rank 15/16, so this also exercises the recursion and the
// p-neighbour base case.
L := DirectSum(E6, E8);
reps := GenusRepresentativesOrth(L);
AssertEqual(~results, #reps, 2, "G_{14,3}: class number");
AssertTrue(~results, forall{ R : R in reps | Genus(R) eq Genus(L) }, "G_{14,3}: reps lie in Genus(L)");
AssertEqual(~results, &+[ 1/#AutomorphismGroup(R) : R in reps ], Mass(L), "G_{14,3}: mass formula");

// Odd lattices: cross-check against the p-neighbour enumeration.
for G in [ DiagonalJoin(IdentityMatrix(Integers(), 3), CartanMatrix("A2")),
           DiagonalJoin(Matrix(Integers(), 1, 1, [3]), CartanMatrix("D4")) ] do
    L := LatticeWithGram(G);
    reps := GenusRepresentativesOrth(L);
    fast := GenusRepresentativesFaster(L);
    AssertEqual(~results, #reps, #fast, Sprintf("odd lattice %o: orbit vs p-neighbour count", G));
    cfs := { Eltseq(CanonicalForm(GramMatrix(R))) : R in reps };
    AssertEqual(~results, cfs, { Eltseq(CanonicalForm(GramMatrix(R))) : R in fast },
        Sprintf("odd lattice %o: orbit vs p-neighbour canonical forms", G));
end for;

// Content reduction: a scaled lattice has the scaled representatives.
reps := GenusRepresentativesOrth(LatticeWithGram(3 * CartanMatrix("D4")));
AssertTrue(~results,
    #reps eq 1 and GramMatrix(reps[1]) eq 3 * GramMatrix(LatticeWithGram(CartanMatrix("D4"))),
    "content reduction: 3*D4");

// Batch enumeration: same answers as single-target calls; targets of equal
// determinant may share one descent sweep with complements binned by genus.
B := [ DirectSum(DirectSum(A1, A1), A1),
       LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 8])),
       LatticeWithGram(DiagonalMatrix(Integers(), [1, 2, 4])) ];
batch := GenusRepresentativesOrthBatch(B);
for i in [1 .. #B] do
    single := GenusRepresentativesOrth(B[i]);
    AssertEqual(~results, #batch[i], #single, Sprintf("batch[%o]: matches single-target count", i));
    AssertTrue(~results, forall{ R : R in batch[i] | Genus(R) eq Genus(B[i]) },
        Sprintf("batch[%o]: reps lie in Genus(B[i])", i));
    AssertEqual(~results, &+[ 1/#AutomorphismGroup(R) : R in batch[i] ], Mass(B[i]),
        Sprintf("batch[%o]: mass formula", i));
end for;

// Indefinite dispatch enumerates EVERY spinor genus in the genus (audit P1):
// diag(1,20,-25) has two spinor genera, hence two classes by Eichler;
// SpinorRepresentatives alone used to return only one of them.
Lind := LatticeWithGram(DiagonalMatrix(Rationals(), [1, 20, -25]) : CheckPositive := false);
ok, reps := GenusReps(Lind : Timeout := 120);
AssertTrue(~results, ok and #reps eq 2, "indefinite dispatch: diag(1,20,-25), both spinor genera");

// Negative definite lattices are enumerated via negation, not the indefinite
// shortcut (audit P1).
Lneg := LatticeWithGram(DiagonalMatrix(Rationals(), [-1, -1, -8]) : CheckPositive := false);
ok, nreps := GenusReps(Lneg : Timeout := 120);
posreps := GenusRepresentativesOrth(LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 8])));
AssertTrue(~results, ok and #nreps eq #posreps, "negative definite: same count as positive analogue");
AssertTrue(~results, forall{ R : R in nreps | IsPositiveDefinite(-GramMatrix(R)) },
    "negative definite: every rep is actually negative definite");

// Batch with Depth := 0 must not perform any glue step (audit: it used one more
// level than allowed); the single-target fallback still gives correct answers.
batch0 := GenusRepresentativesOrthBatch(B : Depth := 0);
for i in [1 .. #B] do
    AssertEqual(~results, #batch0[i], #batch[i], Sprintf("Depth:=0 batch[%o]: matches Depth-unset", i));
end for;

// Mass-aware routing: the neighbour walk keeps small low-rank genera, the
// orbit method takes high rank or large mass.
AssertTrue(~results, not UseOrthHeuristic(14, 1/2), "UseOrthHeuristic(14, 1/2): false");
AssertTrue(~results, UseOrthHeuristic(15, 1/2), "UseOrthHeuristic(15, 1/2): true");
AssertTrue(~results, UseOrthHeuristic(6, 128), "UseOrthHeuristic(6, 128): true");
AssertTrue(~results, not UseOrthHeuristic(6, 0), "UseOrthHeuristic(6, 0): false");

// The dispatch intrinsic agrees on a definite lattice...
ok, reps := GenusReps(DirectSum(E6, E8) : Timeout := 600);
AssertTrue(~results, ok and #reps eq 2, "GenusReps dispatch: E6+E8");
// ... handles the rank-2 square-discriminant case ...
ok, reps := GenusReps(LatticeWithGram(Matrix(Rationals(),2,2,[0,5,5,1]) : CheckPositive := false) : Timeout := 60);
AssertTrue(~results, ok and #reps ge 1, "GenusReps dispatch: square-discriminant rank 2");
// ... and the indefinite spinor case.
ok, reps := GenusReps(LatticeWithGram(DiagonalMatrix(Rationals(), [1,1,1,-1]) : CheckPositive := false) : Timeout := 60);
AssertTrue(~results, ok and #reps eq 1, "GenusReps dispatch: indefinite spinor case");

Report(~results, "test_genus_reps_orth");
