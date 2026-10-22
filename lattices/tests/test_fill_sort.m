// cmp_lat (fill_genus.m) orders the lattices of a genus, and so fixes their labels.
// Stored theta series are Eltseq of a power series, which drops trailing zero
// coefficients, so a list can be shorter than its precision; the comparator must
// read the missing coefficients as 0.  It used to index past the end:
// "Sequence index 182 should be in the range 1 to 181" on 8.8.504.001.01.3827.
AttachSpec("lattices.spec");
import "fill_genus.m" : cmp_lat;
load "tests/assertions.m";
results := NewResults();

function rec(theta, tprec, dual_theta, dprec, gram)
    r := AssociativeArray();
    r["aut_size"] := 48;
    r["minimum"] := 1;
    r["theta_series"] := theta;
    r["theta_prec"] := tprec;
    r["dual_theta_series"] := dual_theta;
    r["dual_theta_prec"] := dprec;
    r["gram"] := gram;
    return r;
end function;

g1 := [1, 0, 0, 2];
g2 := [1, 0, 0, 3];
theta := [1, 4, 4];

// dual theta of A is [1,0,0,2] at precision 8 (coefficients 5..8 are zero and were
// dropped); B has a nonzero coefficient 7.  Equal up to A's stored length, so the
// comparison has to read past it.
A := rec(theta, 3, [1, 0, 0, 2], 8, g1);
B := rec(theta, 3, [1, 0, 0, 2, 0, 0, 3], 8, g2);
AssertTrue(~results, cmp_lat(A, B) lt 0, "short dual theta: A (coefficient 7 = 0) sorts before B (= 3)");
AssertTrue(~results, cmp_lat(B, A) gt 0, "short dual theta: antisymmetric");

// the same for the primal theta series
C := rec([1, 2], 6, [1], 1, g1);
D := rec([1, 2, 0, 0, 5], 6, [1], 1, g2);
AssertTrue(~results, cmp_lat(C, D) lt 0, "short theta: C sorts before D");

// equal series (one stored short, one padded) fall through to the Gram tiebreaker
E := rec(theta, 3, [1, 0, 0, 2], 8, g1);
F := rec(theta, 3, [1, 0, 0, 2, 0, 0, 0], 8, g2);
AssertTrue(~results, cmp_lat(E, F) lt 0, "equal series up to precision: Gram tiebreak");
AssertEqual(~results, cmp_lat(E, E), 0, "a record equals itself");

Report(~results, "test_fill_sort");
