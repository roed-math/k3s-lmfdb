# Bug report: `SpinorCharacters` / `SpinorGenerators` wrong for indefinite lattices

Fields for Magma's bug report form.  Fix: `diagonal-2.30-1.patch` in this directory.

## Title

`SpinorCharacters` and `SpinorGenerators` return wrong spinor characters for indefinite
lattices of rank >= 3 (three bugs in `package/Lattice/Lat/diagonal.m`)

## Magma version

V2.30-1 (avx2-linux), with an unmodified package directory.  We have not run the
example on other versions.

## Minimal reproducible example

Each case is an indefinite ternary lattice L with a rational isometry S of L tensor Q
such that N = L*S is a q-neighbour of L ([L : L meet N] = [N : L meet N] = q, N
integral).  N is then isometric to L, so q is improperly automorphous and every spinor
character of the genus of L (Proper := false) must take the value 1 at q.

```
procedure Check(g, q, s)
  G := Matrix(Rationals(), 3, 3, g);
  S := Matrix(Rationals(), 3, 3, s);
  L := LatticeWithBasis(IdentityMatrix(Rationals(), 3), G : CheckPositive := false);
  N := LatticeWithBasis(S, G : CheckPositive := false);
  M := L meet N;
  X := SpinorCharacters(Genus(L) : Proper := false);
  printf "isometry %o, N integral %o, indices %o %o; SpinorGenerators %o, characters at %o: %o; #GenusRepresentatives %o\n",
    S * G * Transpose(S) eq G, IsIntegral(N), Index(L, M), Index(N, M),
    SpinorGenerators(Genus(L) : Proper := false), q, [Evaluate(x, q) : x in X], #GenusRepresentatives(L);
end procedure;
Check([4,1,0, 1,10,0, 0,0,-54], 5, [1,-1/5,0, 0,-1,0, 0,0,-1]);
Check([-9,0,0, 0,18,0, 0,0,64], 7, [-9/7,-4/7,0, 8/7,9/7,0, 0,0,-1]);
Check([1,0,0, 0,-48,24, 0,24,213], 11, [-13/11,1/11,0, 48/11,-13/11,0, -24/11,12/11,1]);
```

## Expected and actual behavior

Expected: every character is 1 at q, and, since the classes of an indefinite genus of
rank >= 3 are its spinor genera (Eichler), 2^#SpinorGenerators equals
`#GenusRepresentatives`:

```
isometry true, N integral true, indices 5 5; SpinorGenerators [], characters at 5: []; #GenusRepresentatives 1
isometry true, N integral true, indices 7 7; SpinorGenerators [ 5 ], characters at 7: [ 1 ]; #GenusRepresentatives 2
isometry true, N integral true, indices 11 11; SpinorGenerators [ 7 ], characters at 11: [ 1 ]; #GenusRepresentatives 2
```

(this is the output with the patch below).  Actual (full output):

```
isometry true, N integral true, indices 5 5; SpinorGenerators [ 5 ], characters at 5: [ -1 ]; #GenusRepresentatives 1
isometry true, N integral true, indices 7 7; SpinorGenerators [ 5, 7 ], characters at 7: [ -1, 1 ]; #GenusRepresentatives 2
isometry true, N integral true, indices 11 11; SpinorGenerators [ 11 ], characters at 11: [ -1 ]; #GenusRepresentatives 2
Total time: 0.200 seconds, Total memory usage: 32.09MB
```

In each case a spinor character is -1 at an automorphous prime.  In the first two
cases `SpinorGenerators` also contradicts `GenusRepresentatives` (2 spinor genera but 1
class; 4 spinor genera but 2 classes).  In the third the count is right but the
character is wrong: 11 is returned as a spinor generator, but L is isometric to its
11-neighbour, while the patched version finds the actual generator 7.  A program that
builds the classes of such a genus as neighbours for the primes `SpinorGenerators`
returns (as ours did) stores isometric lattices as different classes and misses
classes.

## Suggestions and additional details

**Cause.**  Three independent bugs, one per example (each example is fixed by the
corresponding change alone, given the earlier ones):

1. `SpinorKernel` does not use the signature.  For an indefinite form the real place
   has every spinor norm, so the idele that is -1 at infinity (and 1 elsewhere) is
   automorphous; `SpinorKernel` omits it, which treats the form as definite at
   infinity.  Fix: when the form is indefinite, add the vector of coordinate-character
   values at -1 to the kernel generators (or remove the coordinate if only one
   character is -1 there), as the existing code does for other automorphous classes.
2. In `SpinorCharacters`, when 2 is in M but -1 is not, a coordinate at -1 is added
   (`M := [-1] cat M`) and every kernel vector is extended by 1 there.  The kernel puts
   no condition on that coordinate, so it must stay free: the vector that is -1 at the
   new coordinate and 1 elsewhere has to be added to S.  Without it a spurious
   character survives.
3. `SpinorCharacters` takes `W := Complement(V, sub< V | B >)`.  The spinor characters
   are the products of coordinate characters that are trivial on the kernel, i.e. the
   exponent vectors orthogonal to every vector of B: the annihilator of B, not an
   arbitrary complement of its span.  Both have the right dimension, which is why the
   count is often right while the characters are not.

**Patch** against V2.30-1 (`diagonal-2.30-1.patch`):

```
@@ in SpinorCharacters
    if 2 in M and -1 notin M then 
       M := [-1] cat M; 
-      S := [ [1] cat s : s in S ];
+      S := [ [1] cat s : s in S ] cat [ [-1] cat [ 1 : i in [1..#M-1] ] ];
    end if;
@@
    B := [ V | [ (1-s[i]) div 2 : i in [1..#s] ] : s in S ];
-   W := Complement(V,sub< V | B >);
+   W := #B eq 0 select V else Nullspace(Transpose(Matrix(GF(2), #B, Ngens(X), [ Eltseq(b) : b in B ])));
@@ in SpinorKernel, before the reflection for improper spinor genera
+   G := GramMatrix(L);
+   if not IsPositiveDefinite(G) and not IsPositiveDefinite(-G) then
+      g := [ q eq -1 select -1 else KroneckerSymbol(-1, q) : q in chars ];
+      I := [ k : k in [1..#chars] | g[k] eq -1 ];
+      if #I eq 1 then
+         k := I[1];
+         Remove(~chars,k);
+         gens := { Remove(h,k) : h in gens };
+      elif #I ge 2 then
+         Include(~gens,g);
+      end if;
+   end if;
```

**Testing.**  On 19162 indefinite genera of ranks 3 to 5 (determinant with square
factors, 18 <= |det| <= 10890) the patched `SpinorGenerators` gives the same number of
spinor genera as Hecke/Oscar (`improper_spinor_generators`, with the fix in Hecke.jl
pull request #2611, a separate bug there) for every genus; unpatched it differs in
99.  For all 39 of these genera with more than one spinor genus, the patched
characters define the same map on primes < 200 (not dividing 2 det) as Hecke's
spinor-norm map.  Every explicit isometry we found between classes of these genera is
consistent with the patched counts.

**References.**  The classes of an indefinite genus of rank >= 3 are its spinor genera
(Eichler; see e.g. J. W. S. Cassels, *Rational Quadratic Forms*, Academic Press 1978).
For a prime p not dividing 2 det(L), the p-neighbours of L lie in the spinor genus of L
translated by the idele that is p at p (Kneser's theory of neighbours), so p is
automorphous if and only if a p-neighbour of L lies in the spinor genus of L.  Spinor
kernels and automorphous numbers as used in `diagonal.m`: Conway and Sloane, *Sphere
Packings, Lattices and Groups*, Chapter 15, Section 9.

**Random seed.**  Deterministic; the run above had [Seed = 2871785296].

**Resources.**  The example takes 0.6 s and 34 MB.

**Operating system and binary.**  Ubuntu 26.04.1 LTS, x86_64; `GetBinaryVariant();`
returns `avx2-linux`.
