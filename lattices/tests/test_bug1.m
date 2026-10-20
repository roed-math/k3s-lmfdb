AttachSpec("lattices.spec");
import "canonical_form.m" : test_V, test_canonical, V_cvp;
load "tests/assertions.m";
results := NewResults();

M:=-Matrix(Rationals(), 19, 19 ,[-2, 0, 0, 0, 0, -1, -1, 0, 1, -1, 1, 1, 1, \
-1, 1, -1, 1, 1, -1, 0, -2, -1, -1, -1, 0, 0, 0, -1, 1, -1, -1, -1, 1, 1, 1, 1\
, -1, 1, 0, -1, -2, -1, -1, 0, 0, 0, -1, 1, -1, -1, -1, 0, 1, 1, 1, 0, 1, 0, -\
1, -1, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -2, 0\
, 0, 0, -1, 1, -1, 0, -1, 1, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, -2, -1, 0, 1, -1, \
1, 0, 1, -1, 1, 0, 0, 1, -1, -1, 0, 0, 0, 0, -1, -2, 0, 1, -1, 1, 1, 1, 0, 1, \
-1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, -2, 1, 1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 1, \
-1, -1, 0, -1, 1, 1, 1, -4, 1, -1, -1, -1, 1, 0, 1, 1, -1, 2, -1, 1, 1, 0, 1, \
-1, -1, 1, 1, -4, 3, 2, 3, 0, 0, -2, 0, 1, -2, 1, -1, -1, 0, -1, 1, 1, -1, -1,\
 3, -4, -1, -3, 1, -1, 1, 1, -1, 1, 1, -1, -1, 0, 0, 0, 1, 0, -1, 2, -1, -4, -\
1, -1, 1, 3, -1, -1, 2, 1, -1, -1, 0, -1, 1, 1, -1, -1, 3, -3, -1, -4, 1, -1, \
1, 1, -1, 2, -1, 1, 0, 1, 1, -1, 0, 0, 1, 0, 1, -1, 1, -4, 1, 1, -1, 1, 0, 1, \
1, 1, 1, 1, 1, 1, -1, 0, 0, -1, 1, -1, 1, -4, -1, 0, 0, 0, -1, 1, 1, 0, 0, 0, \
-1, 0, 1, -2, 1, 3, 1, 1, -1, -4, 1, 1, -2, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, -\
1, 1, -1, 0, 1, -4, -1, 0, 1, -1, 0, -1, 0, 1, 1, 0, -1, 1, -1, -1, -1, 1, 0, \
1, -1, -4, 1, -1, 1, 1, 0, 0, -1, -1, 0, 2, -2, 1, 2, 2, 0, 0, -2, 0, 1, -4]);

// test_V and test_canonical (canonical_form.m) assert internally and raise
// (fatal, via SetQuitOnError) rather than returning a pass/fail -- so a
// completed call is itself the check; register it explicitly to keep it in
// the assertion count.
test_V(M);
AssertTrue(~results, true, "test_V(M) completed without raising");
test_canonical(M);
AssertTrue(~results, true, "test_canonical(M) completed without raising");

L:=LatticeWithGram(M);
G:=AutomorphismGroup(L);
s:=V_cvp(M);
print #s;
AssertTrue(~results, &and[i*ChangeRing(Matrix(G.1),Rationals()) in s : i in s],
    "V_cvp(M) is closed under the first automorphism generator");

Report(~results, "test_bug1");
