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
- **`diagonal-2.30-1.patch` (required).**  `package/Lattice/Lat/diagonal.m`: Magma counts
  the classes of an indefinite genus of rank >= 3 as its (improper) spinor genera, through
  `SpinorCharacters`, which had three bugs; `genus_reps_Spinor` builds the classes from
  them.  The patch fixes all three:
  1. `SpinorKernel` ignored the signature, treating an indefinite form as definite at
     infinity; it now adds the idele that is -1 at infinity.
  2. When `SpinorCharacters` adds a coordinate at -1 (2 in M, -1 not), that coordinate
     must stay free: it adds the kernel vector [-1, 1, ..., 1] (found by David,
     2026-11-03).
  3. The characters were an arbitrary complement of the kernel's span (`Complement`),
     which has the right dimension but can contain characters that are nontrivial on
     the kernel; it now takes the annihilator.  With fixes 1 and 2 but not this one,
     the count was right but the character, and so the spinor generator and the
     pipeline's second class, was wrong in 13 genera.
  On 19162 census genera of ranks 3-5 with square factors in det, against stock Hecke
  0.39.21, stock Magma disagrees on 99 class numbers, with fix 1 on 17, with fixes 1-2
  on 13: exactly the 13 genera that Hecke.jl PR #2611 (a Hecke bug in
  `_automorphous_numbers`) moves from 1 class to 2.  Against Hecke with that PR, fixes
  1-2 agree on all 19162 class numbers but the spinor characters differ in 13 of the 39
  multi-class genera; with all three they agree in all 39, at every prime < 200, and
  the pipeline's representatives of the 17 genera Magma and Hecke used to dispute are
  pairwise non-isometric (Hecke).  Without fix 3, in 3.2.10800.0001.11.01.293b and .493b
  the neighbour Magma picks as the second class is explicitly isometric to the first.
  Fix 2 alone changes only 4 genera (the det-10368 ones), and fixes 1 and 3 without it
  leave those 4 wrong.  `check_magma.m` requires the patch (`MagmaSpinorPatched`, one
  example per bug).  Evidence and scripts: `/scratch/claudeF/work/spinor_fix/` on mimir.
  Bug report: `diagonal-2.30-1-report.md`.
