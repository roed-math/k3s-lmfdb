#!/usr/bin/env python3
"""Layer-3 (db-relative) universality/regularity pass over lat_inclusions.

Once the Oscar embedding stage has produced inclusion data, every embedding
between database lattices factors as prime-index overlattice steps followed by a
chain of corank-1 primitive embeddings, so the transitive closure of those two
edge types is exactly db-representability -- except that a chain's intermediate
lattices can have determinants beyond the enumeration window.  Closure therefore
yields SOUND "represents" facts and coverage bounds, but NO sound negatives.
Everything this script emits is accordingly a report for curation / follow-up:

  * db_universality:  for each definite lattice L and each n < rank(L), the
    largest D such that every definite rank-n lattice in the database with
    |det| <= D reaches L in the closure ("L represents every rank-n lattice of
    det <= D"), plus the largest n whose full db det range is covered.
  * db_even_universality: the same restricted to even sources.
  * regularity candidates: for h >= 2 genera, sources M that reach some class of
    the genus but not L.  The genus-representation half is rigorous (the path is
    a witness embedding); "L misses M" is only closure-relative and must be
    finished off by a RepresentsLattice(L, M) exhaustion in Magma before it is a
    theorem.  See oscar/INTEGRATION_PLAN.md section 10.

Inputs: a lat_inclusions table (psycodict copy format) and lattice metadata,
either walked from lattice_basic_data/<rank>/<nplus>/<label> files or supplied as
a simple TSV (label, rank, nplus, disc_abs, is_even) for tests.  Memory scales as
(#lattices of source rank) bits per node, one source rank at a time.
"""

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path


@dataclass
class LatMeta:
    label: str
    rank: int
    nplus: int
    disc_abs: int
    is_even: bool

    @property
    def definite(self):
        return self.rank == self.nplus

    @property
    def genus_label(self):
        return self.label.rsplit(".", 1)[0]


def read_meta_tsv(path):
    """label<tab>rank<tab>nplus<tab>disc_abs<tab>is_even(T/F) rows."""
    lats = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            label, rank, nplus, disc, is_even = row[:5]
            lats[label] = LatMeta(label, int(rank), int(nplus), int(disc),
                                  is_even in ("T", "t", "True", "true", "1"))
    return lats


def read_meta_lattice_dir(root):
    """Walk lattice_basic_data/<rank>/<nplus>/<label> pipe-separated records."""
    root = Path(root)
    fmt = (root.parent / "lat_basic.format").read_text().strip().split("|")
    li, ri, si, di, ei = (fmt.index(k) for k in
                          ("label", "rank", "nplus", "disc_abs", "is_even"))
    lats = {}
    for p in root.glob("*/*/*"):
        pieces = p.read_text().split("\n")[0].split("|")
        lats[pieces[li]] = LatMeta(pieces[li], int(pieces[ri]), int(pieces[si]),
                                   int(pieces[di]), pieces[ei] == "T")
    return lats


def read_inclusions(path, lats):
    """Edges domain -> codomain from a lat_inclusions .tbl (header: 3 lines).

    Rows whose endpoints are unknown or indefinite are skipped (the K3 lattice is
    indefinite, so k3_embed rows drop out here).
    """
    edges = defaultdict(set)          # codomain -> set of domains (reverse edges)
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    cols = lines[0].split("|")
    dom_i, cod_i = cols.index("domain"), cols.index("codomain")
    for line in lines[2:]:
        pieces = line.split("|")
        dom, cod = pieces[dom_i], pieces[cod_i]
        if dom == cod:
            continue
        if dom in lats and cod in lats and lats[dom].definite and lats[cod].definite:
            edges[cod].add(dom)
    return edges


def toposort(lats):
    """Edges preserve (rank asc, det desc within rank): overlattice steps keep the
    rank and divide the determinant by p^2, corank-1 steps raise the rank.  Sorting
    by (rank, -det) is therefore a topological order for the closure DP."""
    return sorted((m for m in lats.values() if m.definite),
                  key=lambda m: (m.rank, -m.disc_abs, m.label))


def reach_from_rank(n, order, edges, lats, even_only=False):
    """For every definite lattice, the set of rank-n sources reaching it.

    Returns {label: frozen bitset (python int)} plus the source index map.
    Processes in topological order; memory is one |sources|-bit int per node.
    """
    sources = [m.label for m in order if m.rank == n and (not even_only or m.is_even)]
    src_bit = {lab: 1 << i for i, lab in enumerate(sources)}
    reach = {}
    for m in order:
        bits = src_bit.get(m.label, 0)    # a source reaches itself
        for dom in edges.get(m.label, ()):
            bits |= reach.get(dom, 0)
            if dom in src_bit:
                bits |= src_bit[dom]
        reach[m.label] = bits
    return reach, sources


def analyze(lats, edges, max_source_rank=None, report=print):
    order = toposort(lats)
    by_genus = defaultdict(list)
    for m in order:
        by_genus[m.genus_label].append(m.label)
    ranks = sorted({m.rank for m in order})
    max_source_rank = max_source_rank or (max(ranks) - 1 if len(ranks) > 1 else 0)

    results = defaultdict(dict)       # label -> {(kind, n): (D, covered, total)}
    candidates = []                   # (label, n, missing source reaching a genus mate)

    for n in [r for r in ranks if r <= max_source_rank]:
        for even_only, kind in ((False, "univ"), (True, "even_univ")):
            reach, sources = reach_from_rank(n, order, edges, lats, even_only)
            if not sources:
                continue
            full = (1 << len(sources)) - 1
            src_index = {lab: i for i, lab in enumerate(sources)}
            det_sorted = sorted(sources, key=lambda l: (lats[l].disc_abs, l))
            for m in order:
                if m.rank <= n:
                    continue
                bits = reach[m.label]
                # coverage bound: largest D with EVERY det <= D source reached
                # (a det level counts only when complete); covered counts all
                # reached sources, including scattered ones past the first gap
                best = 0
                for det, group in groupby(det_sorted, key=lambda l: lats[l].disc_abs):
                    if all(bits >> src_index[lab] & 1 for lab in group):
                        best = det
                    else:
                        break
                covered = bin(bits).count("1")
                results[m.label][(kind, n)] = (best, covered, len(sources))
                # regularity candidates: sources reaching a genus mate but not m
                if not even_only and covered < len(sources):
                    mates = [x for x in by_genus[m.genus_label] if x != m.label]
                    if mates:
                        mate_bits = 0
                        for x in mates:
                            mate_bits |= reach[x]
                        missing = mate_bits & full & ~bits
                        while missing:
                            i = (missing & -missing).bit_length() - 1
                            candidates.append((m.label, n, sources[i]))
                            missing &= missing - 1

    for lab in sorted(results, key=lambda l: (lats[l].rank, lats[l].disc_abs, l)):
        cells = []
        for (kind, n), (D, covered, total) in sorted(results[lab].items()):
            cells.append(f"{kind}[{n}]: det<={D} ({covered}/{total})")
        report(f"{lab}\t" + "\t".join(cells))
    if candidates:
        report("# regularity candidate counterexamples "
               "(verify each with RepresentsLattice before storing anything):")
        for lab, n, src in candidates:
            report(f"# candidate\t{lab}\trank {n}\tmisses {src}")
    return results, candidates


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inclusions", help="lat_inclusions .tbl (psycodict copy format)")
    meta = ap.add_mutually_exclusive_group(required=True)
    meta.add_argument("--lattice-dir", help="lattice_basic_data directory")
    meta.add_argument("--meta-tsv", help="TSV: label rank nplus disc_abs is_even")
    ap.add_argument("--max-source-rank", type=int, default=None)
    args = ap.parse_args()

    lats = (read_meta_lattice_dir(args.lattice_dir) if args.lattice_dir
            else read_meta_tsv(args.meta_tsv))
    edges = read_inclusions(args.inclusions, lats)
    analyze(lats, edges, args.max_source_rank)


if __name__ == "__main__":
    main()
