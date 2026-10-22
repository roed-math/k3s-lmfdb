#!/usr/bin/env python3
"""
Check that a regenerated genera_basic keeps the rep of every indefinite genus of
rank >= 3.

Usage: python3 check_spinor_base.py OLD_GENERA_BASIC NEW_GENERA_BASIC

The classes of such a genus are labelled by their spinor coordinates relative to
the genus's rep (genus_reps_Spinor in genus_reps.m), so a different rep relabels
them.  The rep comes from Sage's genus_symbol.representative() in the listing
stage, which can change between versions (the Hanke wiring changed it), so after
listing again, compare with the genera_basic the published labels came from.
Exits 1 if any rep changed or any genus disappeared.
"""
import os
import sys


def reps(root):
    """label -> rep for the indefinite genera of rank >= 3 under root/<rank>/<nplus>/."""
    fmt = open(os.path.join(root, "..", "genera_basic.format")).read().strip().split("|")
    li, repi = fmt.index("label"), fmt.index("rep")
    out = {}
    for rank in os.listdir(root):
        if not rank.isdigit() or int(rank) < 3:
            continue
        for nplus in os.listdir(os.path.join(root, rank)):
            if not 0 < int(nplus) < int(rank):
                continue
            d = os.path.join(root, rank, nplus)
            for f in os.listdir(d):
                with open(os.path.join(d, f)) as fh:
                    pieces = fh.readline().rstrip("\n").split("|")
                out[pieces[li]] = pieces[repi]
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    old, new = reps(sys.argv[1]), reps(sys.argv[2])
    changed = sorted(l for l in old if l in new and old[l] != new[l])
    missing = sorted(l for l in old if l not in new)
    print(f"{len(old)} indefinite genera of rank >= 3 in {sys.argv[1]}: "
          f"{len(changed)} with a different rep, {len(missing)} missing from {sys.argv[2]}")
    for l in changed[:20]:
        print(f"  rep changed: {l}: {old[l]} -> {new[l]}")
    for l in missing[:20]:
        print(f"  missing: {l}")
    sys.exit(1 if changed or missing else 0)


if __name__ == "__main__":
    main()
