// Genus representatives for rank-2 lattices of square determinant -m^2, where
// Magma's GenusRepresentatives fails (genus_reps.m).
import "genus_reps.m" : genus_reps_square_disc, square_disc_isometric;

// m = 3: every genus has class number 1.
for k in [0..5] do
    L := LatticeWithGram(Matrix(Rationals(),2,2,[0,3,3,k]) : CheckPositive := false);
    assert #genus_reps_square_disc(L) eq 1;
end for;

// m = 5: the class number is NOT always 1 -- some genera have two classes,
// others one.
class_numbers := { #genus_reps_square_disc(
        LatticeWithGram(Matrix(Rationals(),2,2,[0,5,5,k]) : CheckPositive := false))
    : k in [0..9] };
assert 2 in class_numbers;
assert 1 in class_numbers;

// Non-canonical input (det -25) is handled, and the reps lie in its genus.
Lnc := LatticeWithGram(Matrix(Rationals(),2,2,[2,3,3,-8]) : CheckPositive := false);
reps := genus_reps_square_disc(Lnc);
assert #reps eq 2;
assert forall{ R : R in reps | Genus(R) eq Genus(Lnc) };

// Exact isometry test: distinct canonical forms for m = 3 are non-isometric.
assert square_disc_isometric(1, 1, 3);
assert not square_disc_isometric(0, 1, 3);
