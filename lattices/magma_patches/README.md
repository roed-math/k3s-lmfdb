# Patches to Magma's package code

Apply from Magma's root directory (the one containing `package/`), e.g.
`patch -p0 < .../magma_patches/neighbors-2.30-1.patch`.  Written against
V2.30-1; for other versions check that each hunk applies, and port by hand if not.

- **`neighbors-2.30-1.patch` (required).**  `package/Lattice/Lat/neighbors.m`:
  - `GenusRepresentatives` on indefinite binary lattices worked in the wide class
    group, so it returned lattices from other genera and missed classes, and did
    not identify mirror images.  `SpinorRepresentatives` had the same bug.
  - `Neighbours` rejected every indefinite binary lattice and every indefinite
    lattice of content > 1.
  `run_all.py` runs `check_magma.m` before the fill stage and stops if the first
  fix is missing.  Reported to Magma on 2026-10-22.
- **`mass-2.30-1.patch` (required).**  `package/Lattice/Lat/mass.m`:
  `DedekindZetaExact` computed L(chi_d, 1-n), which `Mass` needs for every lattice of
  rank 2 mod 4 with non-square discriminant, by evaluating the L-series numerically
  and recognising the rational with `BestApproximation`.  That gives wrong masses from
  rank 14 on (I_13 + <2127>: off by a relative 1e-30) and runtime errors from rank 26
  on; in a scan of 183 discriminants it was wrong or failed for 14% at rank 14, 67% at
  rank 22, 95% at rank 30.  The patch uses -B_{k,chi_d}/k, exact for every k.  With it
  Magma's masses agree with Sage's on all 545 definite genera of the benchmark sample.
  The orbit method and the p-neighbour walk stop on Magma's `Mass`, so a wrong mass
  breaks enumeration.  `check_magma.m` requires it (`MagmaMassPatched`).  Bug report:
  `mass-2.30-1-report.md`.
- **`diagonal-2.30-1-indefinite.patch` (partial, not yet required).**
  `package/Lattice/Lat/diagonal.m`: `SpinorKernel` ignored the signature, so it
  treated an indefinite form as definite at infinity and doubled the number of
  spinor genera (hence of classes) of many indefinite genera of rank >= 3.  This
  adds the idele that is -1 at infinity.  It brings 82 of 99 overcounted census
  genera in a sample into agreement with Hecke, but 17 with complicated 2- and
  3-adic structure still disagree (explicit isometries show Magma still
  overcounts in 6); see "Indefinite class numbers" in `../CLAUDE.md`.
