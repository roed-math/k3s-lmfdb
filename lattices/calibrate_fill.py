#!/usr/bin/env python3
"""Calibrate the fill stage's guards and batch weights against a real census.

The fill defaults (--enum-masslimit/-timeout/-timelimit/-sizelimit/-per-lattice-timeout,
and genus_work's per-rank costs in run_all.py) were tuned at C=256.  This draws a
stratified sample of genera from a listed census, fills it through the real pipeline
path, and reports per-stratum cost and guard outcomes, extrapolated to the census.

  1. sample: pick genera and copy their genera_basic records into a calibration tree
       python3 calibrate_fill.py sample --census /path/to/census --tree . --per-stratum 8

     Strata are (rank, nminus, floor(log2 det)).  Sampling uniformly over labels would
     almost only see the largest determinants (the genus count grows with det), so each
     stratum gets the same number of draws, and the report weights strata by their
     census population instead.  The census population per stratum is saved to
     calibration_population.tsv.

  2. fill: run the pipeline's fill stage on the sample with per-genus stats on, e.g.
       sage -python run_all.py -r 9 -R 12 --skip-list-genera --fill-stats -j 8 \\
           --skip-embeddings --skip-tensor --skip-names --skip-connect \\
           --enum-timelimit 3600 --enum-masslimit 100000 ...
     Use deliberately generous guards for a cost measurement (a guard that trips is a
     censored observation: the report counts it, but its time is a lower bound), or the
     production values to measure what coverage they would give.

  3. report: summarise fill_stats/ against the sample
       python3 calibrate_fill.py report --tree .

Caveat: production fill runs rank-descending with an on-disk orth cache, so a definite
genus's parent (rank+1) is usually already filled.  A sampled genus's parents mostly are
not, so its orbit-method enumeration starts cold.  Enumeration times here are therefore
upper bounds; the per-lattice loop (canonical form, Aut, theta) is unaffected.
"""
import argparse
import csv
import math
import os
import random
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

OUTCOMES = ["ok", "mass_skip", "enum_fail", "size_skip", "time_discard", "missing"]


def parse_label(label):
    p = label.split(".")
    rank, nplus, det = int(p[0]), int(p[1]), int(p[2])
    return rank, nplus, det


def stratum(label):
    rank, nplus, det = parse_label(label)
    return rank, rank - nplus, int(math.log2(det)) if det > 0 else 0


def iter_census_labels(census):
    """All labels of a listed census, from genus_jobs.txt if present (one read of a
    single file), otherwise by walking genera_basic/ (slow at census scale)."""
    jobs = Path(census, "genus_jobs.txt")
    if jobs.exists():
        with open(jobs) as F:
            for line in F:
                for label in line.strip().split(":"):
                    if label:
                        yield label
    else:
        for dp, _, fns in os.walk(Path(census, "genera_basic")):
            yield from fns


def cmd_sample(a):
    rng = random.Random(a.seed)
    lo, _, hi = a.ranks.partition("-")
    ranks = set(range(int(lo), int(hi or lo) + 1))
    nminus_ok = None if a.nminus is None else {int(x) for x in a.nminus.split(",")}
    # Reservoir-sample per stratum in one pass, so the census label list is never held
    # in memory.
    pop = Counter()
    res = defaultdict(list)
    for label in iter_census_labels(a.census):
        st = stratum(label)
        if st[0] not in ranks or (nminus_ok is not None and st[1] not in nminus_ok):
            continue
        pop[st] += 1
        k = pop[st]
        if k <= a.per_stratum:
            res[st].append(label)
        else:
            j = rng.randrange(k)
            if j < a.per_stratum:
                res[st][j] = label
    tree = Path(a.tree)
    dest = tree / "genera_basic"
    if dest.exists() and any(dest.iterdir()):
        if not a.force:
            sys.exit(f"{dest} is not empty; use --force to replace it")
        shutil.rmtree(dest)
    n = 0
    with open(tree / "calibration_sample.txt", "w") as S:
        for st in sorted(res):
            for label in sorted(res[st]):
                rank, nplus, _ = parse_label(label)
                src = Path(a.census, "genera_basic", str(rank), str(nplus), label)
                d = dest / str(rank) / str(nplus)
                d.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, d / label)
                S.write(label + "\n")
                n += 1
    with open(tree / "calibration_population.tsv", "w") as P:
        P.write("rank\tnminus\tlog2det\tpopulation\tsampled\n")
        for st in sorted(pop):
            P.write("%d\t%d\t%d\t%d\t%d\n" % (*st, pop[st], len(res[st])))
    print(f"sampled {n} genera from {sum(pop.values())} in {len(pop)} strata into {dest}")


def read_stats(tree):
    fmt = Path(tree, "fill_stats.format").read_text().strip().split("|")
    stats = {}
    for dp, _, fns in os.walk(Path(tree, "fill_stats")):
        for fn in fns:
            vals = Path(dp, fn).read_text().strip().split("|")
            if len(vals) == len(fmt):
                r = dict(zip(fmt, vals))
                stats[r["label"]] = r
    return stats


# Where a genus's records go: the LMFDB table, and the format file naming its columns
# (oscar_genus/ is pipeline-internal and not uploaded).
RECORD_DIRS = [("genera_basic", "lat_genera", "genera_basic.format"),
               ("genera_advanced", "lat_genera", "genera_advanced.format"),
               ("oscar_genus", "not uploaded", None),
               ("lattice_basic_data", "lat_lattices", "lat_basic.format"),
               ("lattice_advanced_data", "lat_lattices", "lat_advanced.format")]


def genus_footprint(tree, label, stored, formats):
    """Files and bytes the pipeline keeps for one genus in the calibration tree: its
    genera_basic, genera_advanced and oscar_genus records and, per stored lattice,
    its lattice_basic_data and lattice_advanced_data records (the per-genus hash
    files are not counted: they are named by hash, not label).  Also the bytes per
    (table, column), as text (columns from the first line of each record)."""
    rank, nplus, _ = parse_label(label)
    sub = Path(str(rank), str(nplus))
    files, total, cols = 0, 0, Counter()
    for d, table, fmt in RECORD_DIRS:
        names = [label] if d in ("genera_basic", "genera_advanced", "oscar_genus") else \
                [f"{label}.{i}" for i in range(1, stored + 1)]
        for name in names:
            path = Path(tree, d, sub, name)
            if not path.is_file():
                continue
            text = path.read_text()
            files += 1
            total += len(text.encode())
            first, _, rest = text.partition("\n")
            fields = first.split("|")
            if fmt in formats and len(fields) == len(formats[fmt]):
                for c, v in zip(formats[fmt], fields):
                    cols[(table, c)] += len(v.encode())
                cols[(table, d + " (separators, further lines)")] += len(text.encode()) - sum(len(v.encode()) for v in fields)
            else:
                cols[(table, d)] += len(text.encode())
    return files, total, cols


def read_connect_stats(tree):
    """connect_stats/ (run_all.py --connect-stats): label -> (seconds over both
    passes, any error)."""
    out = {}
    for dp, _, fns in os.walk(Path(tree, "connect_stats")):
        for fn in fns:
            f = Path(dp, fn).read_text().strip().split("|")
            if len(f) != 4:
                continue
            label, _, status, secs = f
            t, err = out.get(label, (0.0, False))
            out[label] = (t + float(secs), err or status != "ok")
    return out


def q(xs, p):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


def cmd_report(a):
    tree = Path(a.tree)
    sample = Path(tree, "calibration_sample.txt").read_text().split()
    pop = {}
    with open(tree / "calibration_population.tsv") as P:
        for row in csv.DictReader(P, delimiter="\t"):
            pop[(int(row["rank"]), int(row["nminus"]), int(row["log2det"]))] = int(row["population"])
    stats = read_stats(tree)
    by = defaultdict(list)
    for label in sample:
        by[stratum(label)].append(stats.get(label))
    rows = []
    hdr = ("rank nminus log2det  pop       n  " + " ".join("%-6s" % o[:6] for o in OUTCOMES)
           + "  t_med   t_p90   t_max   h_med  h_max  maxlat_p90  est_cpu_h")
    print(hdr)
    rank_tot = defaultdict(lambda: [0.0, 0, 0.0])   # (rank, definite) -> [cpu_h, pop, lost pop]
    costs = []                                       # (rank, nminus, log2det, mean seconds)
    for st in sorted(by):
        recs = by[st]
        oc = Counter((r["outcome"] if r else "missing") for r in recs)
        times = [float(r["t_total"]) for r in recs if r]
        hs = [int(r["class_number"]) for r in recs if r and r["class_number"] not in ("\\N", "")]
        ml = [float(r["max_lat_time"]) for r in recs if r]
        mean_t = statistics.mean(times) if times else float("nan")
        if times:
            costs.append((*st, mean_t))
        P = pop.get(st, 0)
        cpu_h = mean_t * P / 3600 if times else float("nan")
        lost = sum(oc[o] for o in OUTCOMES if o != "ok") / len(recs)
        rt = rank_tot[(st[0], st[1] == 0)]
        rt[0] += 0 if math.isnan(cpu_h) else cpu_h
        rt[1] += P
        rt[2] += lost * P
        print("%4d %6d %7d %8d %4d  " % (*st, P, len(recs))
              + " ".join("%-6d" % oc[o] for o in OUTCOMES)
              + "  %6.1f  %6.1f  %6.1f  %6s %6s  %10.1f  %9.1f" % (
                  q(times, 0.5), q(times, 0.9), max(times) if times else float("nan"),
                  q(hs, 0.5) if hs else "-", max(hs) if hs else "-", q(ml, 0.9), cpu_h))
    print("\nper rank (population-weighted extrapolation; censored times are lower bounds):")
    print("rank  type   population  est_cpu_hours  est_genera_not_ok  mean_s_per_genus")
    tot = [0.0, 0, 0.0]
    for (r, definite) in sorted(rank_tot):
        cpu_h, P, lost = rank_tot[(r, definite)]
        tot = [tot[0] + cpu_h, tot[1] + P, tot[2] + lost]
        print("%4d  %-5s  %10d  %13.1f  %17d  %16.2f" % (r, "def" if definite else "indef", P, cpu_h,
                                                        round(lost), cpu_h * 3600 / P if P else 0))
    print("total        %10d  %13.1f  %17d" % (tot[1], tot[0], round(tot[2])))
    # genus_work in run_all.py batches fill jobs by these per-genus costs
    def_cost = {r: v[0] * 3600 / v[1] for (r, d), v in rank_tot.items() if d and v[1]}
    indef = [v for (r, d), v in rank_tot.items() if not d and v[1]]
    print("\nmeasured per-genus costs for genus_work (run_all.py):")
    print("_DEF_COST = {" + ", ".join("%d: %.1f" % (r, c) for r, c in sorted(def_cost.items())) + "}")
    if indef:
        print("indefinite: %.2f s per genus" % (sum(v[0] for v in indef) * 3600 / sum(v[1] for v in indef)))
    if a.costs_out:
        with open(a.costs_out, "w") as F:
            F.write("# rank nminus log2det mean_seconds_per_genus (run_all.py --genus-costs)\n")
            for r, nm, b, t in costs:
                F.write("%d %d %d %.3f\n" % (r, nm, b, t))
        print(f"\nwrote {len(costs)} stratum costs to {a.costs_out} (use with run_all.py --genus-costs)")
    cstats = read_connect_stats(tree)
    if cstats:
        ctot = defaultdict(lambda: [0.0, 0, 0, 0])   # (rank, definite) -> [cpu_h, pop, measured, errors]
        for st in sorted(by):
            ts = [cstats[l][0] for l in sample if stratum(l) == st and l in cstats]
            if not ts:
                continue
            c = ctot[(st[0], st[1] == 0)]
            c[0] += statistics.mean(ts) * pop.get(st, 0) / 3600
            c[1] += pop.get(st, 0)
            c[2] += len(ts)
            c[3] += sum(1 for l in sample if stratum(l) == st and l in cstats and cstats[l][1])
        print("\nconnect (both passes; run_all.py --connect-stats), extrapolated like fill:")
        print("rank  type   population  measured  errors  est_cpu_hours  mean_s_per_genus")
        tot = 0.0
        for (r, definite) in sorted(ctot):
            cpu_h, P, m, e = ctot[(r, definite)]
            tot += cpu_h
            print("%4d  %-5s  %10d  %8d  %6d  %13.1f  %16.2f" % (r, "def" if definite else "indef", P, m, e,
                                                                cpu_h, cpu_h * 3600 / P if P else 0))
        print("total connect cpu hours: %.1f" % tot)
    # storage: lattices, files and bytes per genus, extrapolated like the costs
    store = defaultdict(lambda: [0, 0.0, 0.0, 0.0])    # rank -> [pop, lattices, files, bytes]
    colbytes = Counter()                                # (table, column) -> extrapolated bytes
    formats = {fmt: Path(tree, fmt).read_text().strip().split("|") for _, _, fmt in RECORD_DIRS
               if fmt and Path(tree, fmt).exists()}
    for st in sorted(by):
        recs = [(l, stats[l]) for l in sample if stratum(l) == st and l in stats]
        if not recs:
            continue
        P = pop.get(st, 0)
        lat = [int(r["lattices_stored"]) for _, r in recs]
        fp = [genus_footprint(tree, l, int(r["lattices_stored"]), formats) for l, r in recs]
        for _, _, c in fp:
            for k, v in c.items():
                colbytes[k] += v / len(fp) * P
        s_ = store[st[0]]
        s_[0] += P
        s_[1] += statistics.mean(lat) * P
        s_[2] += statistics.mean(f for f, _, _ in fp) * P
        s_[3] += statistics.mean(b for _, b, _ in fp) * P
    if store:
        print("\nstorage (records in this tree after fill and connect, extrapolated):")
        print("rank  genera(pop)  est_lattices  est_files  est_GB")
        tot = [0, 0.0, 0.0, 0.0]
        for r in sorted(store):
            P, nl, nf, nb = store[r]
            tot = [tot[0] + P, tot[1] + nl, tot[2] + nf, tot[3] + nb]
            print("%4d  %11d  %12.0f  %9.0f  %6.2f" % (r, P, nl, nf, nb / 1e9))
        print("total %11d  %12.0f  %9.0f  %6.2f" % (tot[0], tot[1], tot[2], tot[3] / 1e9))
        # Text bytes, as in the flat files; Postgres differs per type (numeric arrays
        # are larger than their text, TOAST compresses values over ~2 kB) and adds indexes.
        tables = Counter()
        for (t, c), b in colbytes.items():
            tables[t] += b
        print("\nstorage by LMFDB table (text bytes, extrapolated):")
        for t, b in tables.most_common():
            print("  %-13s %7.2f GB" % (t, b / 1e9))
        print("largest columns:")
        for (t, c), b in colbytes.most_common(15):
            print("  %-13s %-45s %7.2f GB  %4.1f%%" % (t, c, b / 1e9, 100 * b / sum(tables.values())))
    missing = [l for l in sample if l not in stats]
    if missing:
        print(f"\n{len(missing)} sampled genera have no fill_stats record (crash, error, or not run): "
              + " ".join(missing[:20]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("--census", required=True, help="directory holding genera_basic/ (and ideally genus_jobs.txt)")
    s.add_argument("--tree", default=".", help="lattices/ directory of the calibration checkout")
    s.add_argument("--per-stratum", type=int, default=8)
    s.add_argument("--ranks", default="1-32")
    s.add_argument("--nminus", default=None, help="comma list of nminus values to keep (default all)")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--force", action="store_true")
    r = sub.add_parser("report")
    r.add_argument("--tree", default=".")
    r.add_argument("--costs-out", default="", help="write per-stratum mean seconds per genus for run_all.py --genus-costs")
    a = ap.parse_args()
    {"sample": cmd_sample, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
