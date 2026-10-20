"""
Large-prime regression sweep for hanke_full.maximal_overlattice_2.
Run: sage tests/test_hanke_large_primes.py   (from lattices/)

Not a systematic sweep. The 184-case correctness sweep documented in
install_fast_maximal_overlattice (hanke_full.py) covers rank 3-16, det
1-196, but only SMALL primes -- and it never exercised the bug that was
actually found: find_isotrop_fp used 200 random trials over F_p and
inferred anisotropy from failure, which is fine for small p but silently
wrong over a large one. maximal_overlattice_2(1009*A_2) returned a
non-maximal lattice of det 3*1009^2 (correct answer: 3) with no
exception. find_isotrop_fp is deterministic now, but nothing in the repo
checks that stays true, so this exists to catch it if it regresses.

Oracle: Sage's own IntegralLattice.maximal_overlattice, not the
Genus(rep) == genus check the small-prime sweep used. That check only
exercises representative() end-to-end, calling maximal_overlattice_2
many times per case -- a single bad prime among many can go unnoticed.
Comparing directly against Sage's trusted implementation, one lattice at
a time, is the tighter test for this specific failure mode.

Sage's maximal_overlattice returns the maximal EVEN overlattice (see the
"which notion" note in maximal_overlattice_2's docstring), which is what
we get here too since every case below is run with p=None, the case that
sets did_finish=True and calls finish() to restore evenness. This does
NOT exercise the single-prime p=<prime> path used inside
local_modification -- that path deliberately skips finish() and isn't
comparable to Sage's maximal_overlattice.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from sage.all import IntegralLattice, matrix, ZZ
from sage.quadratic_forms.genera.genus import Genus
from hanke_full import maximal_overlattice_2

# Small named lattices, not a systematic rank sweep -- enough spread (rank
# 1, 2, 4, 8) to see the bug isn't rank-specific, cheap enough that Sage's
# own maximal_overlattice finishes quickly even at det ~ p^2 * base_det.
# E8 (rank 8, det 1) is here specifically because the small-prime sweep's
# rank 3-16 coverage was never combined with a large prime -- wiring this
# into the live pipeline needs both, not either alone.
BASE_LATTICES = {
    "A1": matrix(ZZ, [[2]]),
    "A2": matrix(ZZ, [[2, 1], [1, 2]]),
    "D4": matrix(ZZ, [[2, 0, 1, 0], [0, 2, 1, 0], [1, 1, 2, 1], [0, 0, 1, 2]]),
    "E8": IntegralLattice("E8").gram_matrix(),
}

# 1009 is the documented failure. 97 and 10007 bracket it (one below,
# one an order of magnitude above); 100003 checks the deterministic
# replacement doesn't degrade two orders of magnitude out.
PRIMES = [97, 1009, 10007, 100003]


def run():
    failures = []
    total = 0
    for name, gram0 in BASE_LATTICES.items():
        for p in PRIMES:
            total += 1
            label = f"{p}*{name}"
            L = IntegralLattice(p * gram0)
            try:
                ours = maximal_overlattice_2(L, do_asserts=True)
            except Exception as e:
                failures.append((label, f"raised {type(e).__name__}: {e}"))
                print(f"  ERROR {label}: {type(e).__name__}: {e}")
                continue
            theirs = L.maximal_overlattice()
            ours_det = ours.gram_matrix().det()
            theirs_det = theirs.gram_matrix().det()
            if ours_det != theirs_det:
                failures.append((label, f"det mismatch: ours={ours_det} sage={theirs_det}"))
                print(f"  FAIL {label}: det mismatch: ours={ours_det} sage={theirs_det}")
                continue
            if Genus(ours.gram_matrix()) != Genus(theirs.gram_matrix()):
                failures.append((label, "same det, different genus"))
                print(f"  FAIL {label}: same determinant ({ours_det}) but different genus")
                continue
            print(f"  ok   {label}: det {theirs_det}")

    print()
    if failures:
        print(f"FAIL: {len(failures)}/{total} cases disagree with Sage's maximal_overlattice")
        for label, msg in failures:
            print(f"  - {label}: {msg}")
        sys.exit(1)
    else:
        print(f"PASS: {total}/{total} cases match Sage's maximal_overlattice")


if __name__ == "__main__":
    run()
