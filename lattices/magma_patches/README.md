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
- **The three `diagonal-2.30-1-*.patch` (spinor genera; all three together, not yet
  required by `check_magma.m`).**  `package/Lattice/Lat/diagonal.m`; they apply in any
  order.  Magma counts the classes of an indefinite genus of rank >= 3 as its spinor
  genera via `SpinorCharacters`, which had three bugs:
  - `diagonal-2.30-1-indefinite.patch`: `SpinorKernel` ignored the signature, treating
    an indefinite form as definite at infinity; adds the idele that is -1 at infinity.
  - `diagonal-2.30-1-characters.patch` (found by David, 2026-11-03): when
    `SpinorCharacters` adds a coordinate at -1 (2 in M, -1 not), that coordinate must
    stay free; adds the kernel vector [-1, 1, ..., 1].
  - `diagonal-2.30-1-complement.patch`: the characters were an arbitrary complement
    of the kernel's span (`Complement`), which has the right dimension but can contain
    characters nontrivial on the kernel; uses the annihilator instead.  With the other
    two patches but not this one, Magma's count was right but its character (and so
    the spinor generator, and the pipeline's second class) was wrong in 13 genera.
  On 19162 census genera of ranks 3-5 with square factors in det, against stock Hecke
  0.39.21, stock Magma disagrees on 99 class numbers, with the indefinite patch on 17,
  with indefinite + characters on 13: exactly the 13 genera that Hecke.jl PR #2611 (a
  Hecke bug in `_automorphous_numbers`) moves from 1 class to 2.  Against Hecke with
  that PR, Magma with indefinite + characters agrees on all 19162 class numbers, but its
  spinor characters differ in 13 of the 39 multi-class genera; with all three patches
  they agree in all 39, at every prime < 200, and the pipeline's representatives of the
  17 genera Magma and Hecke used to dispute are pairwise non-isometric (Hecke).  Without
  the complement patch, in 3.2.10800.0001.11.01.293b and .493b the neighbour Magma picks
  as the second class is explicitly isometric to the first.  The characters patch alone
  changes only 4 genera (the det-10368 ones), and indefinite + complement without it
  leaves those 4 wrong.  Evidence and scripts: `/scratch/claudeF/work/spinor_fix/` on
  mimir.
