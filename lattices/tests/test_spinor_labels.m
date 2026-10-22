// Labels of indefinite genera of rank >= 3 (genus_reps_Spinor in genus_reps.m):
// the classes are the spinor genera, labelled from the genus's rep by their spinor
// coordinates e, at position 1 + sum e_i 2^(i-1).
AttachSpec("lattices.spec");
load "tests/assertions.m";

procedure() // GenusReps forks (TimeoutCall): parse everything first
results := NewResults();

// --- SpinorCoordinates ---------------------------------------------------------
AssertEqual(~results, [SpinorCoordinates(i, 2) : i in [1..4]], [[0, 0], [1, 0], [0, 1], [1, 1]],
    "SpinorCoordinates: binary counting, first coordinate least significant");
AssertEqual(~results, SpinorCoordinates(1, 0), [], "SpinorCoordinates: class number 1");

// Cases: diag(1, 20, -25) (one generator; it was the example that showed
// SpinorRepresentatives alone under-enumerates), and two census genera's reps
// with one and two generators.
cases := [* <"diag(1,20,-25)", DiagonalMatrix(Rationals(), [1, 20, -25]), [3]>,
            <"3.2.5408.00001.2.8fb", Matrix(Rationals(), 3, 3, [-9, -11, 1, -11, 14, 7, 1, 7, 23]), 0>,
            <"3.2.10368.100001.0001.1847", Matrix(Rationals(), 3, 3, [0, 0, 8, 0, 162, 4, 8, 4, 7]), 0> *];
for c in cases do
    L0 := LatticeWithGram(c[2] : CheckPositive := false);
    ok, reps, prov, gens := GenusReps(L0 : Timeout := 300);
    AssertTrue(~results, ok, c[1] * ": enumerated");
    if not ok then continue; end if;
    if c[3] cmpne 0 then
        AssertEqual(~results, gens, c[3], c[1] * ": spinor generators");
    end if;
    AssertEqual(~results, gens, SpinorGenerators(Genus(L0) : Proper := false), c[1] * ": generators are Magma's");
    AssertEqual(~results, #reps, 2^#gens, c[1] * ": one class per spinor coordinate vector");
    AssertTrue(~results, GramMatrix(reps[1]) eq GramMatrix(L0), c[1] * ": class .1 is the rep itself");
    AssertTrue(~results, forall{R : R in reps | Genus(R) eq Genus(L0)}, c[1] * ": every class lies in the genus");
    // Magma's SpinorGenera builds the same classes in the same order, but which
    // neighbour Neighbors returns first can depend on the process's state, so the
    // lattices themselves can differ; oscar/test_spinor_labels.jl compares the
    // classes with an isometry test.
    AssertEqual(~results, #reps, #SpinorGenera(Genus(L0)), c[1] * ": as many classes as SpinorGenera");
end for;

// Content 2 (3.2.2000.21.11.26, the even sublattice of diag(1,20,-25)): Magma's
// Neighbours, and so its SpinorGenera, fail on it ("not positive definite").
L0 := LatticeWithGram(Matrix(Rationals(), 3, 3, [4, 2, 0, 2, -4, 0, 0, 0, 100]) : CheckPositive := false);
ok, reps, prov, gens := GenusReps(L0 : Timeout := 300);
AssertTrue(~results, ok and gens eq [3] and #reps eq 2 and GramMatrix(reps[1]) eq GramMatrix(L0)
    and forall{R : R in reps | Genus(R) eq Genus(L0)}, "content 2: two classes in the genus, the first the rep");

// An indefinite genus with one class has no generators: label .1 = the rep.
L0 := LatticeWithGram(DiagonalMatrix(Rationals(), [1, 1, -1]) : CheckPositive := false);
ok, reps, prov, gens := GenusReps(L0 : Timeout := 60);
AssertTrue(~results, ok and #reps eq 1 and gens eq [], "I_{2,1}: one class, no spinor generators");

Report(~results, "test_spinor_labels");
end procedure();
