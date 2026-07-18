#!/usr/bin/env python3
"""Synthetic-graph test for refine_universality.py (run: python3 test_refine_universality.py).

Fake database: rank-1 lattices <1>, <2>, <3>, <4> and rank-2 lattices I2 plus a
fake class-number-2 genus {P, Q}.  Edges mirror true facts about those lattices:
I2 represents 1 and 2 but not 3; <4> sits inside <1> with index 2 (an overlattice
edge, det 4 -> 1); P represents 3 but its genus-mate Q does not.
"""

import io
import tempfile
from pathlib import Path

from refine_universality import (LatMeta, analyze, apply_to_advanced, closure_flags,
                                 derive_certificates, read_inclusions, read_meta_tsv,
                                 resolve_candidates, write_sidecars)

META = """\
1.1.1.a.1\t1\t1\t1\tF
1.1.2.a.1\t1\t1\t2\tT
1.1.3.a.1\t1\t1\t3\tF
1.1.4.a.1\t1\t1\t4\tT
2.2.1.a.1\t2\t2\t1\tF
2.2.9.a.1\t2\t2\t9\tF
2.2.64.a.1\t2\t2\t64\tF
2.2.64.a.2\t2\t2\t64\tF
"""

TBL = """\
domain|codomain|primitive|quotient_invs|domain_rank|codomain_rank|quotient_rank|multiplicity|orthogonal_complement|quotient_gens|images
text|text|boolean|integer[]|smallint|smallint|smallint|integer|numeric[]|numeric[]|numeric[]

1.1.1.a.1|2.2.1.a.1|T|{}|1|2|1|1|{1}|{}|{}
1.1.2.a.1|2.2.1.a.1|T|{}|1|2|1|1|{2}|{}|{}
1.1.4.a.1|1.1.1.a.1|F|{2}|1|1|0|1|{}|{1}|{}
1.1.1.a.1|2.2.64.a.1|T|{}|1|2|1|1|{64}|{}|{}
1.1.3.a.1|2.2.64.a.1|T|{}|1|2|1|1|{}|{}|{}
1.1.1.a.1|2.2.64.a.2|T|{}|1|2|1|1|{64}|{}|{}
1.1.1.a.1|2.2.9.a.1|T|{}|1|2|1|1|{9}|{}|{}
1.1.2.a.1|2.2.9.a.1|T|{}|1|2|1|1|{}|{}|{}
1.1.3.a.1|2.2.9.a.1|T|{}|1|2|1|1|{}|{}|{}
1.1.4.a.1|2.2.9.a.1|T|{}|1|2|1|1|{}|{}|{}
"""

ADV_FMT = "label|is_regular|universality|universality_proved|even_universality|even_universality_proved|regularity|regularity_proved"

def adv_record(label, **kw):
    vals = {c: kw.get(c, "\\N") for c in ADV_FMT.split("|")[1:]}
    return "|".join([label] + [vals[c] for c in ADV_FMT.split("|")[1:]]) + "\n"


def main():
    with tempfile.TemporaryDirectory() as d:
        meta_path = Path(d) / "meta.tsv"
        tbl_path = Path(d) / "lat_inclusions.tbl"
        meta_path.write_text(META)
        tbl_path.write_text(TBL)

        lats = read_meta_tsv(meta_path)
        assert len(lats) == 8
        edges = read_inclusions(tbl_path, lats)

        lines = []
        results, candidates = analyze(lats, edges, report=lines.append)

        # I2 reaches <1>, <2> directly and <4> through the overlattice chain
        # <4> -> <1> -> I2; <3> (det 3) is missing, so the coverage bound is 2.
        univ = results["2.2.1.a.1"][("univ", 1)]
        assert univ == (2, 3, 4), univ
        # even sources are <2> and <4>, both reached: full even coverage up to det 4
        even = results["2.2.1.a.1"][("even_univ", 1)]
        assert even == (4, 2, 2), even

        # P reaches everything except... P gets <1>, <3>, and <4> via the chain;
        # <2> is missing, so its bound stops at det 1.
        p = results["2.2.64.a.1"][("univ", 1)]
        assert p == (1, 3, 4), p

        # Q misses <3>, which reaches its genus mate P: exactly one regularity
        # candidate, to be confirmed by a RepresentsLattice exhaustion in Magma.
        assert candidates == [("2.2.64.a.2", 1, "1.1.3.a.1")], candidates

        # sanity on the report text
        assert any("univ[1]: det<=2 (3/4)" in l for l in lines), lines

        # ---- bitset flags + apply-to-advanced
        flags, cands2, mmd = closure_flags(lats, edges)
        assert flags["2.2.9.a.1"][1] == {"all_full": True, "even_full": True,
                                         "genus_full": True}, flags["2.2.9.a.1"]
        assert not flags["2.2.1.a.1"][1]["all_full"]           # misses <3>
        assert flags["2.2.1.a.1"][1]["even_full"]              # has <2>, <4>
        assert flags["2.2.64.a.1"][1]["genus_full"]            # P holds the genus OR
        assert not flags["2.2.64.a.2"][1]["genus_full"]        # Q misses <3>
        assert ("2.2.64.a.2", 1, "1.1.3.a.1") in cands2
        assert mmd[("2.2.1.a.1", 1)] == 3              # smallest missing source of I2 is <3>

        fmt_path = Path(d) / "lat_advanced.format"
        fmt_path.write_text(ADV_FMT + "\n")
        adv = Path(d) / "advanced"
        for lab, kw in [
            # Y: universality and regularity should be raised; even_universality is
            # proved and must NOT be touched despite the even_full flag
            ("2.2.9.a.1", dict(even_universality="0", even_universality_proved="true")),
            # P: heuristic regularity should be raised to 1
            ("2.2.64.a.1", dict()),
            # Q: nothing to raise (genus_full false, all_full false)
            ("2.2.64.a.2", dict()),
        ]:
            m = lats[lab]
            pdir = adv / str(m.rank) / str(m.nplus)
            pdir.mkdir(parents=True, exist_ok=True)
            (pdir / lab).write_text(adv_record(lab, **kw))
        seen, changed = apply_to_advanced(lats, flags, adv, fmt_path, report=lines.append)
        assert (seen, changed) == (3, 2), (seen, changed)      # X has no record; Y and P change
        fields = dict(zip(ADV_FMT.split("|"),
                          (adv / "2" / "2" / "2.2.9.a.1").read_text().strip().split("|")))
        assert fields["universality"] == "1" and fields["universality_proved"] == "false"
        assert fields["even_universality"] == "0" and fields["even_universality_proved"] == "true"
        assert fields["regularity"] == "1" and fields["is_regular"] == "true"
        fields = dict(zip(ADV_FMT.split("|"),
                          (adv / "2" / "2" / "2.2.64.a.1").read_text().strip().split("|")))
        assert fields["regularity"] == "1" and fields["universality"] == "\\N"
        q = (adv / "2" / "2" / "2.2.64.a.2").read_text().strip().split("|")
        assert all(x == "\\N" for x in q[1:]), q

        # ---- negative-certificate path: certificates, candidate resolution,
        # ---- and the certified universality upper bound in the sidecars
        emb = Path(d) / "embed_labels"
        man = Path(d) / "data"
        emb.mkdir()
        for job in ("overlat_rank1_all", "corank1_rank1_all"):
            (man / job).mkdir(parents=True)
            (man / job / "manifest.txt").write_text("task=x loaded=4 skipped=0\n")
        C = 64                                             # window: rank-1 dets <= 64

        # complete universe + completed jobs: everything <= 64 is certified
        certs = derive_certificates(emb, man, C, [1, 2])
        assert certs[("overlat", 1)] == 64 and certs[("c1dom", 1)] == 64
        assert certs[("c1cod", 2)] == 32                   # rank-2 window is C//2
        resolved = resolve_candidates(cands2, lats, certs)
        assert resolved["2.2.64.a.2"] == [(1, "1.1.3.a.1", "u")], resolved
        # det(Q) = 64 > 32 = certified codomain bound, so Q's candidate stays "u";
        # shrink Q's det into the certified window and it flips to proven
        lats2 = dict(lats); lats2["2.2.64.a.2"] = LatMeta("2.2.64.a.2", 2, 2, 20, False)
        assert resolve_candidates(cands2, lats2, certs)["2.2.64.a.2"][0][2] == "p"

        # an exclusion at det 3 caps the certified bound at 2: the det-3 source
        # is no longer certified, so the candidate degrades back to "u"
        (emb / "rank1_excluded.txt").write_text("1.1.3.a.1|3|classlimit\n")
        certs3 = derive_certificates(emb, man, C, [1, 2])
        assert certs3[("overlat", 1)] == 2 and certs3[("c1dom", 1)] == 2
        assert resolve_candidates(cands2, lats2, certs3)["2.2.64.a.2"][0][2] == "u"
        (emb / "rank1_excluded.txt").unlink()

        # a loader skip recorded in the manifest caps the bound the same way
        (man / "corank1_rank1_all" / "manifest.txt").write_text(
            "task=corank1 loaded=3 skipped=1\n1.1.2.a.1\n")
        certs4 = derive_certificates(emb, man, C, [1, 2])
        assert certs4[("c1dom", 1)] == 1 and certs4[("overlat", 1)] == 64
        # and a missing manifest certifies nothing
        (man / "overlat_rank1_all" / "manifest.txt").unlink()
        certs5 = derive_certificates(emb, man, C, [1, 2])
        assert certs5[("overlat", 1)] == -1

        # sidecars: with full certificates, I2 (rank 2, misses <3> at det 3 <= 64,
        # det(I2) = 1 <= 32) gets the proven upper bound "not 1-universal"
        sc = Path(d) / "sidecar"
        nw = write_sidecars(sc, lats, flags, mmd, resolved, certs)
        rec = (sc / "2" / "2" / "2.2.1.a.1").read_text().strip().split("|")
        assert rec[4] == "1", rec                          # universality_upper = 1
        assert nw >= 4
        # without certificates the same lattice carries no upper bound
        sc2 = Path(d) / "sidecar2"
        write_sidecars(sc2, lats, flags, mmd, resolved, certs5)
        rec2 = (sc2 / "2" / "2" / "2.2.1.a.1").read_text().strip().split("|")
        assert rec2[4] == "\\N", rec2

    print("refine_universality synthetic test: OK")


if __name__ == "__main__":
    main()
