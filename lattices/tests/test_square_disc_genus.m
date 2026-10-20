// Genus representatives for rank-2 lattices of square determinant -m^2, where
// Magma's GenusRepresentatives fails (genus_reps.m).
AttachSpec("lattices.spec");
import "genus_reps.m" : genus_reps_square_disc, square_disc_isometric;
load "tests/assertions.m";
results := NewResults();

// m = 3: every genus has class number 1.
for k in [0..5] do
    L := LatticeWithGram(Matrix(Rationals(),2,2,[0,3,3,k]) : CheckPositive := false);
    AssertEqual(~results, #genus_reps_square_disc(L), 1, Sprintf("m=3, k=%o", k));
end for;

// m = 5: the class number is NOT always 1 -- some genera have two classes,
// others one.
class_numbers := { #genus_reps_square_disc(
        LatticeWithGram(Matrix(Rationals(),2,2,[0,5,5,k]) : CheckPositive := false))
    : k in [0..9] };
AssertTrue(~results, 2 in class_numbers, "m=5: some genus has class number 2");
AssertTrue(~results, 1 in class_numbers, "m=5: some genus has class number 1");

// Non-canonical input (det -25) is handled, and the reps lie in its genus.
Lnc := LatticeWithGram(Matrix(Rationals(),2,2,[2,3,3,-8]) : CheckPositive := false);
reps := genus_reps_square_disc(Lnc);
AssertEqual(~results, #reps, 2, "non-canonical det -25 input: #reps");
AssertTrue(~results, forall{ R : R in reps | Genus(R) eq Genus(Lnc) },
    "non-canonical det -25 input: reps lie in Genus(Lnc)");

// Exact isometry test: distinct canonical forms for m = 3 are non-isometric.
AssertTrue(~results, square_disc_isometric(1, 1, 3), "square_disc_isometric(1,1,3)");
AssertTrue(~results, not square_disc_isometric(0, 1, 3), "not square_disc_isometric(0,1,3)");

Report(~results, "test_square_disc_genus");
