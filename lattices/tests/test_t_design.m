// The three shell-design-strength methods in tdesign.m must agree (the cutoffs
// in tDesign only choose the fastest, so correctness is method-independent).
AttachSpec("lattices.spec");
import "tdesign.m" : shell_design_strength,
                     shell_design_strength_harmonic,
                     shell_design_strength_Molien;
load "tests/assertions.m";
results := NewResults();

for tup in [* <Lattice("A",2), 5>, <Lattice("A",3), 3>, <Lattice("D",4), 5> *] do
    L := tup[1];  half := ShortestVectors(L);  A := AutomorphismGroup(L);
    label := Sprintf("shell_design_strength(%o)", tup);
    AssertEqual(~results, shell_design_strength(L, half), tup[2], label);
    AssertEqual(~results, shell_design_strength_harmonic(L, half), tup[2],
        "harmonic " * label);
    AssertEqual(~results, shell_design_strength_Molien(L, half, A), tup[2],
        "Molien " * label);
end for;

// The intrinsic itself (direct path): E8 minimal vectors are a 7-design.
AssertEqual(~results, tDesign(Lattice("E",8), ShortestVectors(Lattice("E",8))), 7,
    "tDesign(E8 minimal vectors)");

Report(~results, "test_t_design");
