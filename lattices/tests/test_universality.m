// RepresentsLattice, UniversalityFlags, LatticeUniversalityData and
// RankOneIrregularityWitness from universality.m.  Expected values are
// literature-anchored: Ko 1937 (I_{n+3} is n-universal for n <= 5, E6 never
// embeds in Z^m), Mordell 1930 (I5 is 2-universal), Kim--Kim--Oh's S2,
// Oh's S8 = {I8, E8}.
AttachSpec("lattices.spec");
load "tests/assertions.m";

procedure() // forcing magma to read the full input before forking
results := NewResults();

// ---- RepresentsLattice: known embeddings and non-embeddings
E8 := Lattice("E", 8);
d117 := LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 7]));
for tup in [* <E8, Lattice("A", 2), true, "E8 > A2">,
              <E8, StandardLattice(1), false, "E8 > I1 (even vs odd)">,
              <StandardLattice(5), d117, false, "I5 > [1,1,7] (I3 complement misses 7)">,
              <StandardLattice(9), d117, true, "I9 > [1,1,7] (I7 complement has 7)">,
              <Lattice("D", 4), LatticeWithGram(SymmetricMatrix([2, 1, 4])), false,
               "D4 > [2,1,4] (odd coordinate sum)">,
              <StandardLattice(6), Lattice("A", 3), true, "I6 > A3"> *] do
    dec, rep := RepresentsLattice(tup[1], tup[2]);
    AssertEqual(~results, <dec, rep>, <true, tup[3]>, "RepresentsLattice " * tup[4]);
end for;
dec, rep := RepresentsLattice(StandardLattice(9), Lattice("E", 6)
                              : Aut := AutomorphismGroup(StandardLattice(9)));
AssertEqual(~results, <dec, rep>, <true, false>, "RepresentsLattice I9 > E6 (E6 never in Z^m)");

// ---- UniversalityFlags (empty Theta: everything through the guarded tail search)
// I4: Lagrange, and universal => even-universal.  I3 misses 7 and 28 = 2*14.
// E8 is even so misses 1, and represents every even integer.  A2 misses 4.
for tup in [* <StandardLattice(4), <true, true, true, true>, "I4">,
              <StandardLattice(3), <true, false, true, false>, "I3">,
              <E8, <true, false, true, true>, "E8">,
              <Lattice("A", 2), <true, false, true, false>, "A2"> *] do
    udec, uval, eudec, euval := UniversalityFlags(tup[1]);
    AssertEqual(~results, <udec, uval, eudec, euval>, tup[2], "UniversalityFlags " * tup[3]);
end for;

// with a supplied theta series, small truants are read off without any forking
t3 := ThetaSeries(StandardLattice(3), 20);
udec, uval := UniversalityFlags(StandardLattice(3) : Theta := [Coefficient(t3, i) : i in [0..19]]);
AssertEqual(~results, <udec, uval>, <true, false>, "UniversalityFlags I3 with Theta");

// ---- LatticeUniversalityData: universality (Ko ladder for I_m, criterion ladder else)
for tup in [* <StandardLattice(3), 0, "I3">, <StandardLattice(4), 1, "I4">,
              <StandardLattice(5), 2, "I5">, <StandardLattice(6), 3, "I6">,
              <StandardLattice(7), 4, "I7">, <StandardLattice(8), 5, "I8">,
              <StandardLattice(9), 5, "I9">, <StandardLattice(12), 5, "I12">, <E8, 0, "E8"> *] do
    udef, uval, uproved := LatticeUniversalityData(tup[1]);
    AssertEqual(~results, <udef, uval, uproved>, <true, tup[2], true>, "universality " * tup[3]);
end for;

// a universal rank-4 lattice that is not I_4: 2-universality is impossible below
// rank 5, so the value is decided to be exactly 1
udef, uval, uproved := LatticeUniversalityData(LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 1, 3])));
AssertEqual(~results, <udef, uval, uproved>, <true, 1, true>, "universality [1,1,1,3]");

// E8 perp I8 represents S8 = {I8, E8}, hence is 8-universal (Oh); it cannot be
// 9-universal because its only norm-1 axes are the eight in the I8 factor, so the
// I9 witness fails and the value is decided to be exactly 8
udef, uval, uproved := LatticeUniversalityData(DirectSum(E8, StandardLattice(8)));
AssertEqual(~results, <udef, uval, uproved>, <true, 8, true>, "universality E8 + I8");

// ---- even_universality
// A2 misses 4; D4 fails the [2,1,4] witness; norm-2 vectors of I4 pair evenly with
// norm-4.  E8 embeds every witness in the pool, and no even 2-criterion set exists,
// so its value is the heuristic 2 with the proved flag off
for tup in [* <Lattice("A", 2), 0, true, "A2">, <Lattice("D", 4), 1, true, "D4">,
              <StandardLattice(4), 1, true, "I4">, <E8, 2, false, "E8"> *] do
    _, _, _, eudef, euval, euproved := LatticeUniversalityData(tup[1]);
    AssertEqual(~results, <eudef, euval, euproved>, <true, tup[2], tup[3]>,
        "even_universality " * tup[4]);
end for;

// ---- RankOneIrregularityWitness on synthetic genus theta data
thetas := [[1, 0, 1, 0, 1], [1, 0, 1, 1, 1]];
// class 2 has norm 3, class 1 misses it
AssertEqual(~results, RankOneIrregularityWitness(thetas, 1), 3, "RankOneIrregularityWitness class 1");
AssertEqual(~results, RankOneIrregularityWitness(thetas, 2), 0, "RankOneIrregularityWitness class 2");

Report(~results, "test_universality");
end procedure();
