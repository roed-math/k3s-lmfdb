// RepresentsLattice, UniversalityFlags, LatticeUniversalityData and
// RankOneIrregularityWitness from universality.m.  Expected values are
// literature-anchored: Ko 1937 (I_{n+3} is n-universal for n <= 5, E6 never
// embeds in Z^m), Mordell 1930 (I5 is 2-universal), Kim--Kim--Oh's S2,
// Oh's S8 = {I8, E8}.  TimeoutCall forks inside the eval may print a harmless
// post-success "bad syntax" (same as the fill tests).

// ---- RepresentsLattice: known embeddings and non-embeddings
E8 := Lattice("E", 8);
dec, rep := RepresentsLattice(E8, Lattice("A", 2));
assert dec and rep;
dec, rep := RepresentsLattice(E8, StandardLattice(1));            // even vs odd
assert dec and not rep;
dec, rep := RepresentsLattice(StandardLattice(9), Lattice("E", 6) : Aut := AutomorphismGroup(StandardLattice(9)));
assert dec and not rep;                                           // E6 never in Z^m
d117 := LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 7]));
dec, rep := RepresentsLattice(StandardLattice(5), d117);          // I3 complement misses 7
assert dec and not rep;
dec, rep := RepresentsLattice(StandardLattice(9), d117);          // I7 complement has 7
assert dec and rep;
dec, rep := RepresentsLattice(Lattice("D", 4), LatticeWithGram(SymmetricMatrix([2, 1, 4])));
assert dec and not rep;                                           // odd coordinate sum
dec, rep := RepresentsLattice(StandardLattice(6), Lattice("A", 3));
assert dec and rep;

// ---- UniversalityFlags (empty Theta: everything through the guarded tail search)
udec, uval, eudec, euval := UniversalityFlags(StandardLattice(4));
assert udec and uval and eudec and euval;                         // Lagrange; universal => even-universal
udec, uval, eudec, euval := UniversalityFlags(StandardLattice(3));
assert udec and not uval;                                         // misses 7
assert eudec and not euval;                                       // misses 28 = 2*14
udec, uval, eudec, euval := UniversalityFlags(E8);
assert udec and not uval;                                         // even lattice misses 1
assert eudec and euval;                                           // E8 represents every even integer
udec, uval, eudec, euval := UniversalityFlags(Lattice("A", 2));
assert udec and not uval and eudec and not euval;                 // A2 misses 4

// with a supplied theta series, small truants are read off without any forking
t3 := ThetaSeries(StandardLattice(3), 20);
udec, uval, eudec, euval := UniversalityFlags(StandardLattice(3) : Theta := [Coefficient(t3, i) : i in [0..19]]);
assert udec and not uval;

// ---- LatticeUniversalityData: universality (Ko ladder for I_m, criterion ladder else)
for tup in [* <StandardLattice(3), 0>, <StandardLattice(4), 1>, <StandardLattice(5), 2>,
             <StandardLattice(6), 3>, <StandardLattice(7), 4>, <StandardLattice(8), 5>,
             <StandardLattice(9), 5>, <StandardLattice(12), 5>, <E8, 0> *] do
    udec, uval, eudec, euval := LatticeUniversalityData(tup[1]);
    assert udec and uval eq tup[2];
end for;

// a universal rank-4 lattice that is not I_4: 2-universality is impossible below
// rank 5, so the value is decided to be exactly 1
udec, uval := LatticeUniversalityData(LatticeWithGram(DiagonalMatrix(Integers(), [1, 1, 1, 3])));
assert udec and uval eq 1;

// E8 perp I8 represents S8 = {I8, E8}, hence is 8-universal (Oh); it cannot be
// 9-universal because its only norm-1 axes are the eight in the I8 factor, so the
// I9 witness fails and the value is decided to be exactly 8
udec, uval := LatticeUniversalityData(DirectSum(E8, StandardLattice(8)));
assert udec and uval eq 8;

// ---- even_universality
_, _, eudec, euval := LatticeUniversalityData(Lattice("A", 2));
assert eudec and euval eq 0;                                      // misses 4
_, _, eudec, euval := LatticeUniversalityData(Lattice("D", 4));
assert eudec and euval eq 1;                                      // [2,1,4] witness fails
_, _, eudec, euval := LatticeUniversalityData(StandardLattice(4));
assert eudec and euval eq 1;                                      // norm-2 vectors of I4 pair evenly with norm-4
// E8 embeds every witness in the pool, and no even 2-criterion set exists, so the
// value must stay undecided
_, _, eudec, euval := LatticeUniversalityData(E8);
assert not eudec;

// ---- RankOneIrregularityWitness on synthetic genus theta data
thetas := [[1, 0, 1, 0, 1], [1, 0, 1, 1, 1]];
assert RankOneIrregularityWitness(thetas, 1) eq 3;                // class 2 has norm 3, class 1 misses it
assert RankOneIrregularityWitness(thetas, 2) eq 0;
