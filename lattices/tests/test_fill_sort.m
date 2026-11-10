// cmp_lat (fill_genus.m) orders the lattices of a genus, and so fixes their labels.
// Stored theta series are Eltseq of a power series, which drops trailing zero
// coefficients, so a list can be shorter than its precision; the comparator must
// read the missing coefficients as 0.  It used to index past the end:
// "Sequence index 182 should be in the range 1 to 181" on 8.8.504.001.01.3827.
// Keys whose computation timed out are "\N" and must sort after every computed
// value; comparing only when both are present is not transitive.
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
    r["root_sublattice"] := "0";
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

// --- root system: the key between the minimum and the theta series ------------
G := rec(theta, 3, [1], 1, g2);  G["root_sublattice"] := "A1";
H := rec([1, 0, 0, 9], 3, [1], 1, g1);  H["root_sublattice"] := "2A1";
AssertTrue(~results, cmp_lat(H, G) lt 0, "more roots first, before the theta series is consulted");
AssertTrue(~results, cmp_lat(G, H) gt 0, "root system: antisymmetric");
H["minimum"] := 2;
AssertTrue(~results, cmp_lat(H, G) lt 0, "the minimum decides before the root system");
H["minimum"] := 1;  H["root_sublattice"] := "A1";
AssertTrue(~results, cmp_lat(H, G) lt 0, "equal root systems fall through to the theta series (H's is smaller)");
H["root_sublattice"] := "\\N";
AssertTrue(~results, cmp_lat(G, H) lt 0 and cmp_lat(H, G) gt 0, "a missing root system sorts after a computed one");

// --- missing keys sort last, so the order is transitive --------------------------
// aut sizes 2, \N, 4: comparing only when both are present made a ~ b ~ c with c < a,
// so the result of Sort depended on the input order.
a := rec(theta, 3, [1], 1, [1, 0, 0, 1]);  a["aut_size"] := 2;
b := rec(theta, 3, [1], 1, [1, 0, 0, 2]);  b["aut_size"] := "\\N";
c := rec(theta, 3, [1], 1, [1, 0, 0, 3]);  c["aut_size"] := 4;
for input in [[a, b, c], [b, c, a], [c, a, b], [a, c, b], [b, a, c], [c, b, a]] do
    srt := Sort(input, cmp_lat);
    AssertEqual(~results, [r["gram"][4] : r in srt], [3, 1, 2], "aut sizes 4, 2, then the missing one, from any input order");
end for;
b["aut_size"] := 2;  b["minimum"] := "\\N";
AssertTrue(~results, cmp_lat(a, b) lt 0, "a missing minimum sorts after a computed one");
// theta series to different precisions: a coefficient past one series' precision is
// unknown and sorts after a known one (it used to stop at the shorter precision)
P := rec([1, 2, 0, 4], 5, [1], 1, g2);
Q := rec([1, 2], 3, [1], 1, g1);
AssertTrue(~results, cmp_lat(P, Q) lt 0 and cmp_lat(Q, P) gt 0, "a series known further sorts first when they agree as far as both go");
R := rec([1, 2, 0, 9], 5, [1], 1, [1, 0, 0, 5]);
AssertTrue(~results, cmp_lat(P, R) lt 0 and cmp_lat(P, Q) lt 0 and cmp_lat(Q, R) gt 0,
    "theta comparisons at mixed precisions are transitive (P < R < Q)");
AssertTrue(~results, cmp_lat(R, Q) lt 0, "theta comparisons at mixed precisions: R < Q");

Report(~results, "test_fill_sort");
