// Orbit-method genus enumeration (genus_reps.m).  Class numbers of the genera
// G_{n,p} of even lattices of determinant p (n even) or 2p (n odd) are known
// (Chenevier-Taibi, Tables 1.2/1.3); seeds from their Table 7.2, padded by E8.

E8 := LatticeWithGram(CartanMatrix("E8"));
E6 := LatticeWithGram(CartanMatrix("E6"));
A1 := LatticeWithGram(CartanMatrix("A1"));
A2 := LatticeWithGram(CartanMatrix("A2"));
A4 := LatticeWithGram(CartanMatrix("A4"));

// Single-class genera (mass = 1/|Aut| short-circuit).
assert #GenusRepresentativesOrth(E8) eq 1;
assert #GenusRepresentativesOrth(A2) eq 1;
assert #GenusRepresentativesOrth(DirectSum(A1, A2)) eq 1;   // G_{3,3}

// G_{11,3} (det 6): class number 2.
L := DirectSum(DirectSum(A1, A2), E8);
reps := GenusRepresentativesOrth(L);
assert #reps eq 2;
assert forall{ R : R in reps | Genus(R) eq Genus(L) };
assert &+[ 1/#AutomorphismGroup(R) : R in reps ] eq Mass(L);

// G_{12,5} (det 5): class number 2.
L := DirectSum(A4, E8);
reps := GenusRepresentativesOrth(L);
assert #reps eq 2;
assert forall{ R : R in reps | Genus(R) eq Genus(L) };

// G_{14,3} (det 3): class number 2.  The parent chain passes through unimodular
// lattices of rank 15/16, so this also exercises the recursion and the
// p-neighbour base case.
L := DirectSum(E6, E8);
reps := GenusRepresentativesOrth(L);
assert #reps eq 2;
assert forall{ R : R in reps | Genus(R) eq Genus(L) };
assert &+[ 1/#AutomorphismGroup(R) : R in reps ] eq Mass(L);

// Odd lattices: cross-check against the p-neighbour enumeration.
for G in [ DiagonalJoin(IdentityMatrix(Integers(), 3), CartanMatrix("A2")),
           DiagonalJoin(Matrix(Integers(), 1, 1, [3]), CartanMatrix("D4")) ] do
    L := LatticeWithGram(G);
    reps := GenusRepresentativesOrth(L);
    fast := GenusRepresentativesFaster(L);
    assert #reps eq #fast;
    cfs := { Eltseq(CanonicalForm(GramMatrix(R))) : R in reps };
    assert cfs eq { Eltseq(CanonicalForm(GramMatrix(R))) : R in fast };
end for;

// Content reduction: a scaled lattice has the scaled representatives.
reps := GenusRepresentativesOrth(LatticeWithGram(3 * CartanMatrix("D4")));
assert #reps eq 1 and GramMatrix(reps[1]) eq 3 * GramMatrix(LatticeWithGram(CartanMatrix("D4")));

// Batch enumeration: same answers as single-target calls; targets of equal
// determinant may share one descent sweep with complements binned by genus.
B := [ DirectSum(DirectSum(A1, A1), A1),
       LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 8])),
       LatticeWithGram(DiagonalMatrix(Integers(), [1, 2, 4])) ];
batch := GenusRepresentativesOrthBatch(B);
for i in [1 .. #B] do
    single := GenusRepresentativesOrth(B[i]);
    assert #batch[i] eq #single;
    assert forall{ R : R in batch[i] | Genus(R) eq Genus(B[i]) };
    assert &+[ 1/#AutomorphismGroup(R) : R in batch[i] ] eq Mass(B[i]);
end for;

// Indefinite dispatch enumerates EVERY spinor genus in the genus (audit P1):
// diag(1,20,-25) has two spinor genera, hence two classes by Eichler;
// SpinorRepresentatives alone used to return only one of them.
Lind := LatticeWithGram(DiagonalMatrix(Rationals(), [1, 20, -25]) : CheckPositive := false);
ok, reps := GenusReps(Lind : Timeout := 120);
assert ok and #reps eq 2;

// Negative definite lattices are enumerated via negation, not the indefinite
// shortcut (audit P1).
Lneg := LatticeWithGram(DiagonalMatrix(Rationals(), [-1, -1, -8]) : CheckPositive := false);
ok, nreps := GenusReps(Lneg : Timeout := 120);
posreps := GenusRepresentativesOrth(LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 8])));
assert ok and #nreps eq #posreps;
assert forall{ R : R in nreps | IsPositiveDefinite(-GramMatrix(R)) };

// Batch with Depth := 0 must not perform any glue step (audit: it used one more
// level than allowed); the single-target fallback still gives correct answers.
batch0 := GenusRepresentativesOrthBatch(B : Depth := 0);
for i in [1 .. #B] do
    assert #batch0[i] eq #batch[i];
end for;

// Mass-aware routing: the neighbour walk keeps small low-rank genera, the
// orbit method takes high rank or large mass.
assert not UseOrthHeuristic(14, 1/2);
assert UseOrthHeuristic(15, 1/2);
assert UseOrthHeuristic(6, 128);
assert not UseOrthHeuristic(6, 0);

// The dispatch intrinsic agrees on a definite lattice...
ok, reps := GenusReps(DirectSum(E6, E8) : Timeout := 600);
assert ok and #reps eq 2;
// ... handles the rank-2 square-discriminant case ...
ok, reps := GenusReps(LatticeWithGram(Matrix(Rationals(),2,2,[0,5,5,1]) : CheckPositive := false) : Timeout := 60);
assert ok and #reps ge 1;
// ... and the indefinite spinor case.
ok, reps := GenusReps(LatticeWithGram(DiagonalMatrix(Rationals(), [1,1,1,-1]) : CheckPositive := false) : Timeout := 60);
assert ok and #reps eq 1;
