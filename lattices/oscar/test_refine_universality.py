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

from refine_universality import LatMeta, analyze, read_inclusions, read_meta_tsv

META = """\
1.1.1.a.1\t1\t1\t1\tF
1.1.2.a.1\t1\t1\t2\tT
1.1.3.a.1\t1\t1\t3\tF
1.1.4.a.1\t1\t1\t4\tT
2.2.1.a.1\t2\t2\t1\tF
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
"""


def main():
    with tempfile.TemporaryDirectory() as d:
        meta_path = Path(d) / "meta.tsv"
        tbl_path = Path(d) / "lat_inclusions.tbl"
        meta_path.write_text(META)
        tbl_path.write_text(TBL)

        lats = read_meta_tsv(meta_path)
        assert len(lats) == 7
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

    print("refine_universality synthetic test: OK")


if __name__ == "__main__":
    main()
