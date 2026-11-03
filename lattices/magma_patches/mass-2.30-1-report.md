# Bug report: `Mass` / `DedekindZetaExact` wrong for lattices of rank 2 mod 4

Fields for Magma's bug report form.  Fix: `mass-2.30-1.patch` in this directory.

## Title

`Mass(L)` is wrong (or errors) for many lattices of rank 2 mod 4: `DedekindZetaExact`
evaluates L(chi_d, 1-k) numerically instead of using generalized Bernoulli numbers

## Magma version

V2.30-1 (avx2-linux).  Also reproduced with V2.29-10 (same wrong value for the
example below).

## Minimal reproducible example

```
K := QuadraticField(-2127);
DedekindZetaExact(K, -6 : Relative);
-BernoulliNumber(7, KroneckerCharacter(-2127, Rationals())) / 7;

L := LatticeWithGram(DiagonalMatrix(Integers(), [1 : i in [1..13]] cat [2127]));
Mass(L);
```

## Expected and actual behavior

Expected: `DedekindZetaExact(K, -6 : Relative)` is L(chi_{-2127}, -6) = -B_{7,chi}/7,
an integer, so the first two lines print the same number:

```
-16032206109030617744
-16032206109030617744
```

and `Mass(L)` is the mass of the genus of I_13 + <2127>, which Sage
(`Genus(diagonal_matrix(ZZ, [1]*13 + [2127])).mass()`) and Oscar (`mass(L)`) both give
as

```
12561949209912499297859/54655451136000
```

Actual (full output):

```
Warning: in L-series calculation, loss of 2 digits due to cancellation
-1101724814752318821679547421329/68719476737
-16032206109030617744
Warning: in L-series calculation, loss of 2 digits due to cancellation
96684064568219242834132043053569053/420660128323746329001984000
```

The L-value is a wrong rational (denominator 2^36 + 1), and the mass is off by the
factor N/(N - 1), N = 1101724814752318821679547421329: a relative error of about
1e-30.  No error is raised, only the warning.

## Suggestions and additional details

**Cause.**  In `package/Lattice/Lat/mass.m`, `standard_mass` multiplies by
`DedekindZetaExact(E, 1-n : Relative)` for a lattice of even rank m = 2n with
non-square discriminant, where E is the quadratic field of the discriminant.
`DedekindZetaExact` returns an exact value via Bernoulli numbers only when z is odd;
when n is odd (rank 2 mod 4) z = 1 - n is even, so it evaluates `LSeries(...)`
numerically and accepts the rational that `BestApproximation(x, 10^i)` returns for
i = 12, 14, 16, 18 when all four agree.  For large L-values the precision is not enough
and all four agree on a wrong rational.  Higher ranks need larger values, and from
rank 26 `BestApproximation` itself fails ("Overflow or precision loss in truncation").

**Scale.**  Comparing `DedekindZetaExact(QuadraticField(D), z : Relative)` with
-B_{k,chi_D}/k for 183 negative fundamental discriminants D with |D| <= 3000 (every
fifth one):

| z | rank | wrong | runtime error |
|---|---|---|---|
| 0, -2, -4 | 2, 6, 10 | 0 | 0 |
| -6 | 14 | 25 | 0 |
| -8 | 18 | 26 | 0 |
| -10 | 22 | 123 | 0 |
| -12 | 26 | 76 | 87 |
| -14 | 30 | 22 | 152 |

The error is "BestApproximation: Overflow or precision loss in truncation" (mass.m line
632).  We found it through a lattice database whose genus enumeration checks
completeness against `Mass`: an inexact mass breaks that check.

**Fix.**  For a quadratic field K and every k = 1 - z >= 1,
L(chi_d, 1 - k) = -B_{k,chi_d}/k (Washington, *Introduction to Cyclotomic Fields*,
2nd ed., GTM 83, Theorem 4.2), and it is 0 when chi_d(-1) <> (-1)^k.  So the relative
quadratic case needs no numerical L-series for either parity of z.  Two details:
`AbsoluteDiscriminant(Integers(K))` returns |d| (2127 for Q(sqrt(-2127))), so the sign
must be restored before building `KroneckerCharacter(d, Rationals())`, and
`BernoulliNumber(k, chi)` asserts chi(-1) = (-1)^k, so the vanishing case must be
handled before calling it.  Patch against V2.30-1:

```
--- package/Lattice/Lat/mass.m
+++ package/Lattice/Lat/mass.m
@@ in DedekindZetaExact
   require (Relative and z eq 0) or z lt 0 : "The argument must be a negative integer";

+  k:= 1-z;
+  if AbsoluteDegree(K) eq 2 and Relative then
+    // AbsoluteDiscriminant is |d|; chi_d needs the sign of d.
+    d:= AbsoluteDiscriminant(Integers(K)) * (IsTotallyReal(K) select 1 else -1);
+    chi:= KroneckerCharacter(d, Rationals());
+    if Evaluate(chi, -1) ne (-1)^k then return 0; end if;
+    return BernoulliNumber(k, chi) / -k;
+  end if;
+
   if IsOdd(z) then
-    k:= 1-z;
     if AbsoluteDegree(K) eq 1 then
       return BernoulliNumber(k)/-k;
     elif AbsoluteDegree(K) eq 2 and IsTotallyReal(K) then
       d:= AbsoluteDiscriminant(Integers(K));
-      if Relative then
-        return BernoulliNumber(k, KroneckerCharacter(d, Rationals())) / -k;
-      else
-        return BernoulliNumber(k) * BernoulliNumber(k, KroneckerCharacter(d, Rationals())) / k^2;
-      end if;
+      return BernoulliNumber(k) * BernoulliNumber(k, KroneckerCharacter(d, Rationals())) / k^2;
     end if;
   end if;
```

(the full diff, with a comment, is `mass-2.30-1.patch`).  With it, the example gives the
expected values without warnings, `Mass` agrees with Sage's mass on all 545 positive
definite genera of ranks 2-16 we tested (one of which, rank 14, was wrong before), and
the scan above has no errors at any rank; it is also 3-4 times faster.  The numerical
fallback that remains (fields of higher degree) has the same weakness: requiring the
value to be stable under an increase of the L-series precision, not just under the
`BestApproximation` bound, would catch it.

**Random seed.**  Deterministic: the same wrong values with every seed tried
(e.g. [Seed = 1861081374], [Seed = 4134879535]).

**Resources.**  The example takes about 1 s and 32 MB.

**Operating system and binary.**  Ubuntu 26.04.1 LTS, x86_64; `GetBinaryVariant();`
returns `avx2-linux`.
