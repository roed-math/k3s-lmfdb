// Genus representatives of indefinite binary lattices with non-square
// discriminant must be the ISOMETRY classes of the genus (genus_reps.m,
// genus_reps_Magma).  Unpatched Magma's GenusRepresentatives gets these wrong --
// lattices outside the genus, missing classes, a class next to its mirror image
// -- so this suite fails on it; the fix is magma_patches/.  Class numbers below
// were checked with Sage (BinaryQF_reduced_representatives plus Genus).
AttachSpec("lattices.spec");
import "genus_reps.m" : genus_reps_Magma, binary_isometric, genus_reps_square_disc, square_disc_isometric, square_disc_key;
load "tests/assertions.m";
results := NewResults();

ok, msg := MagmaBinaryGenusPatched();
AssertTrue(~results, ok, "Magma carries magma_patches/BinaryGenusSpinor-2.30-1.patch: " * msg);

Lat := func<g | LatticeWithGram(Matrix(Integers(), 2, 2, g) : CheckPositive := false)>;

// x^2 - 37y^2: the genus has 3 proper classes, {1, g, g^-1}; g and g^-1 are
// mirror images, so it has 2 isometry classes.
reps := genus_reps_Magma(Lat([1, 0, 0, -37]));
AssertEqual(~results, #reps, 2, "x^2 - 37y^2: class number");
AssertTrue(~results, forall{R : R in reps | Genus(R) eq Genus(Lat([1, 0, 0, -37]))},
    "x^2 - 37y^2: reps lie in the genus");
AssertTrue(~results, not binary_isometric(GramMatrix(reps[1]), GramMatrix(reps[2])),
    "x^2 - 37y^2: the two reps are not isometric");

// The mirror test itself: (4, 6, -7) and (-4, 10, 3) are improperly equivalent.
AssertTrue(~results, binary_isometric(Matrix(Integers(), 2, 2, [4, 3, 3, -7]),
                                      Matrix(Integers(), 2, 2, [-4, 5, 5, 3])),
    "(4,6,-7) ~ (-4,10,3) improperly");
AssertTrue(~results, binary_isometric(Matrix(Integers(), 2, 2, [2, 1, 1, -13]),
                                      Matrix(Integers(), 2, 2, [2, -1, -1, -13])),
    "a form is isometric to its mirror");

// x^2 - 79y^2: 3 proper classes, 2 isometry classes; unpatched Magma usually
// returns 3 lattices, one of them from the other genus.
reps := genus_reps_Magma(Lat([1, 0, 0, -79]));
AssertEqual(~results, #reps, 2, "x^2 - 79y^2: class number");
AssertTrue(~results, forall{R : R in reps | Genus(R) eq Genus(Lat([1, 0, 0, -79]))},
    "x^2 - 79y^2: reps lie in the genus");
// x^2 - 34y^2: -x^2 + 34y^2 is a second class in the same genus, which unpatched
// Magma never returned.
AssertEqual(~results, #genus_reps_Magma(Lat([1, 0, 0, -34])), 2, "x^2 - 34y^2: class number");
AssertEqual(~results, #genus_reps_Magma(Lat([2, 1, 1, -13])), 1, "[[2,1],[1,-13]]: class number");
// definite binary genera are untouched (Magma already returns isometry classes)
AssertEqual(~results, #genus_reps_Magma(Lat([2, 1, 1, 12])), 2, "[[2,1],[1,12]]: class number");
AssertEqual(~results, #genus_reps_Magma(Lat([1, 0, 0, 23])), 2, "x^2 + 23y^2: class number");

// square_disc_key is a complete invariant: it agrees with the exact
// square_disc_isometric test on every pair of canonical forms, and is unchanged by
// a change of basis.
agree := true; inv := true;
for m in [1..8] do
    for k1, k2 in [0..2*m-1] do
        G1 := Matrix(Integers(), 2, 2, [0, m, m, k1]);
        G2 := Matrix(Integers(), 2, 2, [0, m, m, k2]);
        if (square_disc_key(G1) eq square_disc_key(G2)) ne square_disc_isometric(k1, k2, m) then
            agree := false;
        end if;
    end for;
    for k in [0..2*m-1], U in [Matrix(Integers(), 2, 2, [2, 1, 1, 1]), Matrix(Integers(), 2, 2, [1, 3, 0, -1]),
                               Matrix(Integers(), 2, 2, [-1, 0, 4, 1])] do
        G := Matrix(Integers(), 2, 2, [0, m, m, k]);
        if square_disc_key(Transpose(U)*G*U) ne square_disc_key(G) then inv := false; end if;
    end for;
end for;
AssertTrue(~results, agree, "square_disc_key agrees with square_disc_isometric (m <= 8)");
AssertTrue(~results, inv, "square_disc_key is basis-independent");

// IsIsometricBinary on both kinds, including a scrambled basis
U := Matrix(Integers(), 2, 2, [2, 1, 1, 1]);
A := Matrix(Integers(), 2, 2, [4, 3, 3, -7]);
AssertTrue(~results, IsIsometricBinary(Lat(Eltseq(A)), Lat(Eltseq(Transpose(U)*Matrix(Integers(),2,2,[-4,5,5,3])*U))),
    "IsIsometricBinary: mirror forms, scrambled basis");
sq := genus_reps_square_disc(LatticeWithGram(Matrix(Rationals(), 2, 2, [2, 3, 3, -8]) : CheckPositive := false));
AssertTrue(~results, not IsIsometricBinary(sq[1], sq[2]), "IsIsometricBinary: the two square-disc classes differ");
AssertTrue(~results, IsIsometricBinary(sq[1], LatticeWithGram(Transpose(U)*GramMatrix(sq[1])*U : CheckPositive := false)),
    "IsIsometricBinary: square-disc class vs itself in another basis");

Report(~results, "test_binary_genus");
