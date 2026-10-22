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
- **`diagonal-2.30-1-indefinite.patch` (partial, not yet required).**
  `package/Lattice/Lat/diagonal.m`: `SpinorKernel` ignored the signature, so it
  treated an indefinite form as definite at infinity and doubled the number of
  spinor genera (hence of classes) of many indefinite genera of rank >= 3.  This
  adds the idele that is -1 at infinity.  It brings 82 of 99 overcounted census
  genera in a sample into agreement with Hecke, but 17 with complicated 2-adic
  structure still disagree; see "Indefinite class numbers" in `../CLAUDE.md`.
