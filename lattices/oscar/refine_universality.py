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
    if len(lines) < 2:
        return edges                  # empty table (headers only, or nothing)
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
    by (rank, -det) is therefore a topological order for the closure DP.  A
    corrupted edge violating this order would merely fail to propagate its mask --
    a conservative under-approximation, never an over-claim."""
    return sorted((m for m in lats.values() if m.definite),
                  key=lambda m: (m.rank, -m.disc_abs, m.label))


def reach_from_rank(n, order, edges, lats, even_only=False, det_cap=None):
    """For every definite lattice, the set of rank-n sources reaching it.

    Sources are the definite rank-n lattices, each assigned a distinct power of
    two (optionally only even ones, optionally capped at det <= det_cap so the
    tests speak about exactly the range the embedding stage covered); masks are
    OR-ed up the embedding poset in topological order.  Returns
    {label: bitset (python int)} plus the source list.  Memory is one
    |sources|-bit int per node, one source rank at a time.
    """
    sources = [m.label for m in order if m.rank == n and (not even_only or m.is_even)
               and (det_cap is None or m.disc_abs <= det_cap)]
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


def closure_flags(lats, edges, max_source_rank=None, disc_ratio=None):
    """The three bitset tests, per definite lattice and source rank n:

      all_full   -- mask over the rank-n sources is all ones: every rank-n lattice
                    in the universe embeds (the db-relative n-universality test)
      even_full  -- the same restricted to even sources
      genus_full -- mask equals the bitwise OR of the masks over the lattice's
                    genus: everything the genus receives, the lattice receives
                    (the db-relative n-regularity test; the genus side is
                    witness-backed, so a strict inequality names candidate
                    counterexamples for a RepresentsLattice confirmation)

    disc_ratio caps the rank-n source universe at det <= disc_ratio // n, matching
    the embedding stage's coverage convention.  Returns ({label: {n: flag dict}},
    candidates) with candidates = [(label, n, missing source label)].
    """
    order = toposort(lats)
    by_genus = defaultdict(list)
    for m in order:
        by_genus[m.genus_label].append(m.label)
    ranks = sorted({m.rank for m in order})
    max_source_rank = max_source_rank or (max(ranks) - 1 if len(ranks) > 1 else 0)
    flags = defaultdict(dict)
    candidates = []
    missing_min_det = {}
    for n in [r for r in ranks if r <= max_source_rank]:
        cap = None if disc_ratio is None else disc_ratio // n
        reach, sources = reach_from_rank(n, order, edges, lats, det_cap=cap)
        if not sources:
            continue
        full = (1 << len(sources)) - 1
        even_mask = 0
        for i, lab in enumerate(sources):
            if lats[lab].is_even:
                even_mask |= 1 << i
        det_order = sorted(range(len(sources)), key=lambda i: lats[sources[i]].disc_abs)
        genus_or = {}
        for g, labs in by_genus.items():
            acc = 0
            for x in labs:
                acc |= reach[x]
            genus_or[g] = acc
        for m in order:
            if m.rank <= n:
                continue
            bits = reach[m.label]
            flags[m.label][n] = {
                "all_full": bits & full == full,
                "even_full": even_mask != 0 and bits & even_mask == even_mask,
                "genus_full": bits == genus_or[m.genus_label],
            }
            if bits & full != full:
                for i in det_order:
                    if not bits >> i & 1:
                        missing_min_det[(m.label, n)] = lats[sources[i]].disc_abs
                        break
            missing = genus_or[m.genus_label] & ~bits
            while missing:
                i = (missing & -missing).bit_length() - 1
                candidates.append((m.label, n, sources[i]))
                missing &= missing - 1
    return flags, candidates, missing_min_det


def label_det(label):
    """Genus and lattice labels are r.nplus.det.... , so the det is the third piece."""
    return int(label.split(".")[2])


def derive_certificates(embed_labels_dir, manifests_dir, disc_ratio, ranks):
    """Completeness certificates from job manifests + exclusion lists.

    A driver job's committed directory (inclusions/data/<jobtag>/, containing
    manifest.txt) appears atomically on clean completion, and its
    computation enumerates ALL embeddings for the pairs it loaded (Nikulin
    machinery), so the certified bound for a task at source rank n is the
    largest B such that every definite rank-n genus with det <= B was included
    (not excluded by build_embedding_inputs, not skipped by the loader) and the
    covering job completed.  Returns {("overlat", n) | ("c1dom", n) |
    ("c1cod", n): B} with B = -1 when nothing is certified.
    """
    certs = {}
    excl_min = {}
    for n in ranks:
        cap = disc_ratio // n
        m = cap + 1
        p = Path(embed_labels_dir) / f"rank{n}_excluded.txt"
        if p.is_file():
            for line in p.read_text().splitlines():
                lab, det, _ = line.split("|")
                m = min(m, int(det))
        excl_min[n] = m                              # smallest hole in the universe
    def manifest_bound(name, want_ranks):
        p = Path(manifests_dir) / name / "manifest.txt"
        if not p.is_file():
            return {r: -1 for r in want_ranks}       # job did not complete: nothing certified
        bound = {r: excl_min.get(r, 0) - 1 for r in want_ranks}
        for line in p.read_text().splitlines()[1:]:
            lab = line.strip()
            if not lab:
                continue
            r = int(lab.split(".")[0])
            if r in bound:
                bound[r] = min(bound[r], label_det(lab) - 1)   # loader skipped it
        return bound
    for n in ranks:
        certs[("overlat", n)] = manifest_bound(f"overlat_rank{n}_all", [n])[n]
        b = manifest_bound(f"corank1_rank{n}_all", [n, n + 1])
        certs[("c1dom", n)] = b[n]
        certs[("c1cod", n + 1)] = b[n + 1]
    return certs


def resolve_candidates(candidates, lats, certs):
    """Split candidates into proven-missing (certificate-backed) and unresolved.

    Sound only for corank 1: the closure detects M -> L through overlattice steps
    (det only decreases, so chains stay inside a certified window) followed by one
    corank-1 edge; if the rank-n overlattice data and the (n, n+1) corank-1 data
    are certified past det(M) and det(L), a missing closure path proves M does not
    embed in L.  Higher-corank candidates keep flag "u" for the Magma
    RepresentsLattice confirmation in connect.
    """
    out = defaultdict(list)                          # label -> [(n, src, flag)]
    for lab, n, src in candidates:
        flag = "u"
        L, M = lats[lab], lats[src]
        if L.rank == n + 1:
            b = min(certs.get(("overlat", n), -1), certs.get(("c1dom", n), -1))
            bc = certs.get(("c1cod", n + 1), -1)
            if M.disc_abs <= b and L.disc_abs <= bc:
                flag = "p"
        out[lab].append((n, src, flag))
    return out


def write_sidecars(out_root, lats, flags, missing_min_det, resolved, certs, max_cands=8):
    """One record per definite lattice with anything to say, LabelPath-sharded:

        label|db_universality|db_even_universality|db_regularity|universality_upper|cands

    where cands = "n,source,flag;..." (flag p = proven missing by a certificate,
    u = needs a RepresentsLattice confirmation), sorted by (n, source det) and
    capped at max_cands, and universality_upper = n means "provably not
    n-universal" (a certified-complete rank-n slice with a missing source of
    certified det).  Consumed by ConnectGenus, which is the single writer of the
    final columns.
    """
    def contiguous(fl, key):
        n = 0
        while fl.get(n + 1, {}).get(key):
            n += 1
        return n
    root = Path(out_root)
    written = 0
    for m in lats.values():
        if not m.definite:
            continue
        fl = flags.get(m.label, {})
        cands = sorted(resolved.get(m.label, ()),
                       key=lambda c: (c[0], label_det(c[1])))[:max_cands]
        db_u, db_e, db_r = (contiguous(fl, k) for k in ("all_full", "even_full", "genus_full"))
        upper = "\\N"
        n = m.rank - 1
        if n >= 1 and n in fl and not fl[n]["all_full"]:
            b = min(certs.get(("overlat", n), -1), certs.get(("c1dom", n), -1))
            md = missing_min_det.get((m.label, n))
            if md is not None and md <= b and m.disc_abs <= certs.get(("c1cod", m.rank), -1):
                upper = str(n)
        if not (db_u or db_e or db_r or cands or upper != "\\N"):
            continue
        d = root / str(m.rank) / str(m.nplus)
        d.mkdir(parents=True, exist_ok=True)
        cstr = ";".join(f"{n},{src},{flag}" for n, src, flag in cands)
        (d / m.label).write_text(
            f"{m.label}|{db_u}|{db_e}|{db_r}|{upper}|{cstr}\n")
        written += 1
    return written


APPLY_COLS = ("is_regular", "universality", "universality_proved", "even_universality",
              "even_universality_proved", "regularity", "regularity_proved")


def apply_to_advanced(lats, flags, advanced_root, format_path, report=print):
    """Patch the db-relative heuristics into lattice_advanced_data records.

    Raise-only, and proved values are never touched: a `*_proved = true` value is a
    theorem from the criterion layer, and a missing bit in the closure is not a
    disproof (chain intermediates can leave the det window), so this pass only
    strengthens `proved = false` (or null) values.  The value written is the largest
    n such that the test holds for every source rank 1..n.  Returns
    (records seen, records changed).
    """
    fmt = Path(format_path).read_text().strip().split("|")
    idx = {c: fmt.index(c) for c in APPLY_COLS}
    root = Path(advanced_root)

    def contiguous(fl, key):
        n = 0
        while fl.get(n + 1, {}).get(key):
            n += 1
        return n

    seen = changed = 0
    for m in lats.values():
        if not m.definite:
            continue
        p = root / str(m.rank) / str(m.nplus) / m.label
        if not p.is_file():
            continue
        seen += 1
        fields = p.read_text().split("\n")[0].split("|")
        if len(fields) != len(fmt):
            report(f"# skipping {m.label}: advanced record has {len(fields)} fields, format has {len(fmt)}")
            continue
        fl = flags.get(m.label, {})
        this = False

        def raise_value(val_col, proved_col, db_val):
            nonlocal this
            if db_val <= 0 or fields[idx[proved_col]] == "true":
                return
            cur = fields[idx[val_col]]
            cur_num = int(cur) if cur not in ("\\N", "") else -1
            if db_val > cur_num:
                fields[idx[val_col]] = str(db_val)
                fields[idx[proved_col]] = "false"
                this = True

        raise_value("universality", "universality_proved", contiguous(fl, "all_full"))
        raise_value("even_universality", "even_universality_proved", contiguous(fl, "even_full"))
        db_r = contiguous(fl, "genus_full")
        raise_value("regularity", "regularity_proved", db_r)
        if this and db_r >= 1 and fields[idx["regularity_proved"]] == "false" \
                and int(fields[idx["regularity"]]) >= 1:
            fields[idx["is_regular"]] = "true"
        if this:
            p.write_text("|".join(fields) + "\n")
            changed += 1
    return seen, changed


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
    ap.add_argument("--apply", action="store_true",
                    help="patch the db-relative heuristics into --advanced-dir records "
                         "(raise-only; proved values untouched) instead of only reporting")
    ap.add_argument("--advanced-dir", default=None, help="lattice_advanced_data directory")
    ap.add_argument("--format", default=None, help="lat_advanced.format path")
    ap.add_argument("--disc-ratio", type=int, default=None,
                    help="cap the rank-n source universe at det <= C/n (the embedding "
                         "stage's coverage convention)")
    ap.add_argument("--candidates-out", default=None,
                    help="write regularity candidate counterexamples (label, rank, missing "
                         "source) to this file for the Magma confirmation pass")
    args = ap.parse_args()

    lats = (read_meta_lattice_dir(args.lattice_dir) if args.lattice_dir
            else read_meta_tsv(args.meta_tsv))
    edges = read_inclusions(args.inclusions, lats)
    if args.apply:
        flags, cands, _ = closure_flags(lats, edges, args.max_source_rank, args.disc_ratio)
        adv = args.advanced_dir or str(Path(args.lattice_dir).parent / "lattice_advanced_data")
        fmt = args.format or str(Path(args.lattice_dir).parent / "lat_advanced.format")
        seen, changed = apply_to_advanced(lats, flags, adv, fmt)
        print(f"updated {changed} of {seen} advanced records")
        if args.candidates_out:
            Path(args.candidates_out).write_text(
                "".join(f"{l}\t{n}\t{s}\n" for l, n, s in cands))
            print(f"{len(cands)} regularity candidates -> {args.candidates_out}")
    else:
        analyze(lats, edges, args.max_source_rank)


if __name__ == "__main__":
    main()
