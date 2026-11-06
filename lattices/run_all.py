#!/usr/bin/env -S sage -python

import os
import sys
import time
import shutil
import argparse
import hashlib
import subprocess
import contextlib
from collections import defaultdict
from pathlib import Path
# Fully initialise Sage: under "sage -python" the lazy imports are not pre-resolved,
# which breaks genus enumeration in two ways on Sage 10.7 -- a missing
# laurent_polynomial_ring attribute (via discriminant_group(...).normal_form()),
# and a sage.symbolic.function circular import (via the genus mass computation).
import sage.all  # noqa: F401
from genus import write_all_of_sig_between_genera_basic
from sage.parallel.decorate import parallel as sage_parallel

@contextlib.contextmanager
def timed(label):
    print(f"=== {label} ...", flush=True)
    t0 = time.monotonic()
    yield
    STAGE_TIMES.append((label, round(time.monotonic() - t0, 1)))
    print(f"=== {label}: done in {time.monotonic() - t0:.1f}s", flush=True)

STAGE_TIMES = []
RUN_ID = time.strftime("%Y%m%dT%H%M%S") + f"-{os.getpid()}"

def provenance(event, **fields):
    """Append a record to run_provenance.jsonl: which code, tools and options produced
    the data in this directory.  Labels and several invariants depend on the Magma
    version (CanonicalForm, patched intrinsics) and on the enumeration route, and a
    run that pulls new code partway mixes versions, so record them per run."""
    import json, socket
    rec = {"run": RUN_ID, "event": event, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    if event == "start":
        def git(*a):
            r = subprocess.run(["git", *a], capture_output=True, text=True)
            return r.stdout.strip() if r.returncode == 0 else None
        rec.update(host=socket.gethostname(), argv=sys.argv, commit=git("rev-parse", "HEAD"),
                   dirty=bool(git("status", "--porcelain", "--untracked-files=no")),
                   submodules=git("submodule", "status"), python=sys.version.split()[0])
        try:
            import sage.version
            rec["sage"] = sage.version.version
        except Exception:
            pass
    rec.update(fields)
    with open("run_provenance.jsonl", "a") as F:
        F.write(json.dumps(rec) + "\n")

parser = argparse.ArgumentParser("run_all", description="Run the whole pipeline for generating lattice data")

# Bounds on rank and discriminant
parser.add_argument("-r", "--min-rank", type=int, default=1, help="Lower bound on rank")
parser.add_argument("-R", "--max-rank", type=int, default=32, help="Upper bound on rank")
parser.add_argument("-d", "--min-disc-ratio", type=int, default=0, metavar="C", help="Enumerate genera of discriminant more than C/rank")
parser.add_argument("-D", "--max-disc-ratio", type=int, default=32768, help="Enumerate genera of discriminant at most C/rank")
parser.add_argument("-k", "--nokthree", action="store_true", help="By default, we also enumerate genera of discriminant up to d/(22-n) for nminus=1 or 2.  This option turns that off")

# Parallelization
parser.add_argument("-j", "--jobs", type=int, default=128, help="Number of parallel processes to use")
parser.add_argument("-b", "--batch-mass", type=int, default=128, help="DEPRECATED, ignored.  Batches are now sized by estimated wall-clock work (--batch-work-target), not mass.  Mass is a terrible proxy for cost: high-rank genera have huge automorphism groups, so their mass is ~1e-6, and hundreds of expensive genera packed into one mass-bounded batch that then ran for tens of minutes.  Kept only so existing invocations don't error.")
parser.add_argument("--batch-work-target", type=float, default=120.0, help="Batch genus-enumeration instances together until their estimated fill cost reaches this many seconds, so every batch is ~equal wall-clock regardless of rank/signature mix.  Costs come from a measured per-genus calibration (see genus_work): rank<=8 and all indefinite genera are seconds, while a single definite rank 9-12 genus costs minutes and so lands in a batch of its own.  This is the primary batching control.")
parser.add_argument("--genus-costs", default="", help="TSV of measured per-genus fill costs (rank, nminus, log2det, seconds), as written by `calibrate_fill.py report --costs-out`.  Batching then uses the cost of each genus's stratum (or the nearest determinant stratum of its rank and signature) instead of the built-in C=256 per-rank table.")
parser.add_argument("--retry-missing", action="store_true", help="Refill only the genera listed in fill_missing_<rank>.txt (left by an earlier run: errors and crashes), one per job, with this run's guards; also retries failed Oscar enumeration records.  Implies --skip-list-genera.")
parser.add_argument("--batch-size", type=int, default=32, help="Backstop cap on genera per batch, regardless of work target.  Mainly bounds the per-job command-line length (one label each); for cheap genera the work target rarely fills a batch, so this keeps those batches from growing unwieldy.")

# Don't store too many lattices
parser.add_argument("--enum-masslimit", type=int, default=1000, help="If the mass of a genus is larger than this threshold, don't even try to enumerate lattices within")
parser.add_argument("--enum-per-lattice-timeout", type=int, default=0, help="PER-OPERATION wall-clock cap (seconds) applied INDEPENDENTLY to each of a lattice's own computations -- canonical form, automorphism group, primal theta, dual theta, and the smaller guarded calls (minimum, root system, successive minima, ...) -- so a single lattice can spend several times this value; it is NOT a per-lattice aggregate.  The aggregate wall-clock guard is the per-genus --enum-timelimit, not this option.  0 (the default) instead splits enum-timeout across the genus, giving each operation enum-timeout/class_number.  That split looks like starvation -- at enum-timeout 300 a class-number-55 genus allows each operation only 6s, and CanonicalForm alone measures 8-66s at rank 12, so ~1%% of lattices get \\N in those columns -- but it keeps the genus finishing under --enum-timelimit and storing its lattices.  Setting a fixed budget here loosens that, and since a genus overrunning --enum-timelimit has ALL of its lattices discarded, it can convert partial data into no data unless --enum-timelimit is raised to match.  Measured on a rank-12 class-number-55 genus: the default stored 55 lattices (4 with a null canonical_gram) in 813s, while perlattice=120 at timelimit=900 stored 0 lattices in 1007s.  Erasing that ~1%% of nulls costs roughly 8x the runtime, so raise this only deliberately, together with --enum-timelimit.")
parser.add_argument("--enum-timeout", type=int, default=300, help="Seconds allowed for the genus-representatives computation of a single genus.  This was previously not passed through at all, so it silently fell back to run_fill_genus.m's 60s default and 1.8%% of definite genera (all rank 11-12) failed to enumerate entirely -- no class number and no lattices, a hole that null-scanning cannot see since the genus is simply absent from the per-lattice tables.  The computation is hard-killed past this, so it still bounds the worst case.")
parser.add_argument("--enum-timelimit", type=int, default=900, help="Maximum wall-clock seconds for a genus's per-lattice loop.  Overrunning it does NOT keep partial data: FillGenus discards every lattice of that genus (a partial set would contradict the recorded class_number) and keeps only the genus-level record.  So this is a coverage switch, not a speed knob -- and because it is wall-clock, it makes database contents depend on machine load, which is why --enum-masslimit is the better place to express a real coverage policy.  Calibrated at C=256: definite rank 9-12 genera need 97-1005s to COMPLETE (medians 221/460/384/693s by rank), so 900 completes ~97%% of them; 300 (the old default) silently discarded every lattice of most rank 10-12 genera, and 60 discarded essentially all of rank 9-12.")
parser.add_argument("--enum-sizelimit", type=int, default=1000, help="For genera with class number larger than this, do not store individual lattices within the genus")
parser.add_argument("--enum-adjlimit", type=int, default=0, help="Work budget for the adjacency (Hecke) matrix: skip it when the estimated work (~class_number * sum_p p^(rank-2) over the Hecke primes p) exceeds this; 0 = no budget, always attempt it.  Defaults to 0 because the estimate mis-prices reality badly -- at rank 7 the sum_p p^5 term crosses a 20000 budget at class number 6, yet a whole rank-7 genus fills in ~4s, so a budget of 20000 silently dropped Hecke data for 65 rank-7 genera (4.4%% of definite genera overall) that were never expensive.  AdjacencyMatrix now runs under a hard-killed TimeoutCall, which bounds it by measured time instead of a bad model; prefer that.  Set this only if you specifically want to skip adjacency without even trying.")
parser.add_argument("--no-orth", action="store_true", help="Disable the orbit-method genus enumeration (use the p-neighbour walk only)")
parser.add_argument("--warm-limit", type=int, default=0, help="Wall-clock budget (seconds) for each fill job's batch cache pre-warm (GenusRepresentativesOrthBatch over the batch's definite genera); 0 (the default) skips the pre-warm.  It runs in the batch's own Magma process, so forked per-genus enumerations inherit the warmed cache, and that is also why it cannot be bounded by TimeoutCall: the budget is only checked at safe points in the orth code, and a long call inside a C builtin never reaches one.  On the C=32768 calibration a rank-8 batch spent 25 h at 100%% CPU in it without finishing its first genus (each of whose fills takes ~2 min cold).  Its benefit is sharing descents within a batch; the rank-descending disk cache already supplies the parents, so it is off by default.  Set it (e.g. 300) only where batches hold many related definite genera.")

# Embedding stage (Julia/Oscar; see oscar/INTEGRATION_PLAN.md and the README for setup)
parser.add_argument("--embed-tasks", default="overlat,corank1,k3", help="Comma list of embedding computations: overlat (prime-index overlattices), corank1 (rank n -> n+1 primitive embeddings), k3 (primitive-embedding counts into the K3 lattice)")
parser.add_argument("--embed-jobs", type=int, default=8, help="Oscar workers for the embedding stage (a pool, oscar/pool.jl --task embed).  Separate from -j: each worker loads Oscar once (~30 s, more under load) and needs a few GB of memory")
parser.add_argument("--embed-classlimit", type=int, default=200, help="Skip genera with class number above this in the embedding stage (its cost grows with the class lists on both sides)")
parser.add_argument("--embed-max-disc-ratio", type=int, default=0, help="Determinant cap C for the embedding stage, in the enumeration convention (rank n up to C/n); 0 = same as --max-disc-ratio")
parser.add_argument("--embed-k3-min-rank", type=int, default=12, help="Lowest rank for the k3 task: the orthogonal complement inside the rank-22 K3 lattice has rank 22-n, and Oscar's embedding machinery has to enumerate its genus, which explodes as n decreases")
parser.add_argument("--embed-chunk", type=int, default=200, help="Domain genera per embedding job.  Each task's domain genera are det-sorted and split into chunks of this size, one GNU parallel job each, so a timeout or OOM loses one chunk (and caps the completeness certificate just below that chunk's smallest det) instead of a whole rank or rank pair.  Finished chunks are committed atomically and skipped on a rerun without --clean.  Each job pays a fixed cost -- Julia/Oscar startup and compilation, plus loading its codomain or overlattice pool -- measured at ~270s per job on a loaded 32-core machine (94 chunks of 25 at R=4, D=120 cost 3.8x the 7 unchunked jobs in total job time), so keep chunks at well over that much work.")
parser.add_argument("--embed-memlimit", type=float, default=16, help="Memory cap in GB per Oscar worker in the embedding stage; a job over it has its worker killed and replaced and its chunk stays uncommitted.")
parser.add_argument("--embed-max-rows-per-domain", type=int, default=10000, help="Volume guard: a domain genus whose embedding task produces more rows than this has its rows dropped and is listed as overflowed (in its job's units.jsonl and manifest, which also stops the universality certificates below its det) rather than writing an unbounded payload; 0 = no cap.  Changing it reruns every embedding job.")
parser.add_argument("--embed-timelimit", type=int, default=3600, help="Wall-clock seconds per embedding job; the Oscar pool kills and replaces a worker whose job overruns it, and only that chunk stays uncommitted")

# Skip stages
parser.add_argument("--memlimit", type=int, default=32, help="Memory cap in GB for each Magma process in the fill and connect stages (SetMemoryLimit; 0 = none).  A job can hold two processes at once (the batch and a TimeoutCall child), so budget ~2x this per -j slot.  A child over the cap dies and that computation counts as failed; a batch process over it ends, and run_fill_genus.sh reports its unfinished genera.  Without a cap, rank-17 genera on the C=32768 census ran Magma processes up to 220 GB.")
parser.add_argument("--connect-stats", action="store_true", help="Have each connect pass write connect_stats/<rank>/<nplus>/<label>.<phase> (label|phase|status|seconds), so a calibration sample prices connect too (calibrate_fill.py report).")
parser.add_argument("--fill-stats", action="store_true", help="Have FillGenus write one fill_stats/<rank>/<nplus>/<label> record per genus (fields in fill_stats.format): the outcome -- ok, or which guard decided it (mass_skip, enum_fail, size_skip, time_discard) -- the class number, per-phase wall-clock and the slowest single lattice.  Off by default (it is another file per genus); calibrate_fill.py turns it on.")
parser.add_argument("--skip-list-genera", action="store_true", help="Assume that genera have already been listed")
parser.add_argument("--no-oscar-enum", action="store_true", help="Enumerate every genus with Magma.  By default the genera Oscar handles better -- indefinite of rank >= --oscar-indefinite-min-rank (Magma's spinor kernel overcounts some of them) and definite of rank <= --oscar-definite-max-rank (Oscar is faster there) -- are enumerated first by a pool of Oscar workers (oscar/pool.jl --task enumerate, records in oscar_genus/), and fill takes their classes from those records.")
parser.add_argument("--oscar-jobs", type=int, default=0, help="Oscar workers for the enumeration pool (0 = -j).  Each loads Oscar once (~30 s, more under load) and then serves many genera.")
parser.add_argument("--oscar-timeout", type=int, default=600, help="Seconds per genus in the Oscar enumeration pool; a genus that overruns has its worker killed and replaced and gets a 'timeout' record.  A definite genus whose Oscar enumeration failed is enumerated by Magma in fill; an indefinite one of rank >= 3 is not (its Magma class number may be wrong) and gets no class number.")
parser.add_argument("--oscar-memlimit", type=float, default=16, help="Memory cap in GB per Oscar worker process; a worker over it is killed and replaced, and its genus gets a 'memory' record.")
parser.add_argument("--oscar-batch", type=int, default=20, help="Genera per remote call in the Oscar pool.  Batching keeps cheap genera from making the master a bottleneck; a genus that fails is charged alone and the rest of its batch is requeued.")
parser.add_argument("--oscar-definite-max-rank", type=int, default=8, help="Highest rank of definite genera enumerated with Oscar (the Magma/Oscar benchmark: Oscar faster up to rank 8, Magma from rank 9).")
parser.add_argument("--oscar-indefinite-min-rank", type=int, default=3, help="Lowest rank of indefinite genera enumerated with Oscar.  Rank 2 stays with (patched) Magma.")
parser.add_argument("--allow-unpatched-magma", action="store_true", help="Fill even if Magma fails check_magma.m.  Unpatched Magma gets indefinite binary genera wrong (classes outside the genus, missing classes), masses in ranks 2 mod 4 from rank 14 on, and the spinor characters of indefinite genera of rank >= 3; the census must not be run with it.  See magma_patches/.")
parser.add_argument("--skip-enumerate-genera", action="store_true", help="Skip the fill stage (genus enumeration).  This does NOT imply --skip-list-genera: listing still runs unless that is also given, which is how the census listing was run (listing only, with --skip-embeddings).")
parser.add_argument("--skip-embeddings", action="store_true", help="Do not compute primitive and K3 embedding data")
parser.add_argument("--skip-festi-veniani", action="store_true", help="Do not fill in festi_veniani_index after the fill stage (it stays null; run oscar/fix_festi_veniani.jl later)")
parser.add_argument("--fv-jobs", type=int, default=8, help="Parallel Julia workers for festi_veniani_index (each loads Oscar: ~30s and a few GB)")
parser.add_argument("--fv-chunk", type=int, default=5000, help="Genera per festi_veniani_index job; |O(q_L)| takes milliseconds for most genera and under 0.5s for the hardest census discriminant forms, so jobs are dominated by Oscar's startup")
parser.add_argument("--fv-timelimit", type=int, default=3600, help="GNU parallel timeout (seconds) per festi_veniani_index job; genera of a job that times out keep a null index")
parser.add_argument("--skip-theta-hecke", action="store_true", help="Do not run the theta-Hecke pass after fill: is_theta_distinguished stays null for the genera fill could not decide (some theta series agree to the known precision, and extending them to the Sturm bound was too expensive).  The pass decides them with Eichler's commutation relation between p-neighbours and Hecke operators (theta_hecke.m).")
parser.add_argument("--theta-hecke-jobs", type=int, default=8, help="Parallel Magma jobs for the theta-Hecke pass; genera sharing a cusp space go to one job, which builds the space once")
parser.add_argument("--theta-hecke-chunk", type=int, default=200, help="About this many genera per theta-Hecke job (whole cusp-space groups are never split)")
parser.add_argument("--theta-hecke-timeout", type=int, default=600, help="Seconds (HardKill) per adjacency matrix and per cusp space in the theta-Hecke pass; over it the genus stays undecided")
parser.add_argument("--skip-tensor", action="store_true", help="Do not compute tensor product decompositions")
parser.add_argument("--skip-names", action="store_true", help="Do not find lattice names")
parser.add_argument("--skip-connect", action="store_true", help="Do not run_connect_genus.m")

parser.add_argument("--clean", action="store_true", help="Delete the pipeline's output directories before running, so the result cannot silently blend this run's output with stale files from earlier runs.  Every stage writes per-file with Overwrite and never deletes, so anything this run does NOT produce keeps whatever an earlier run left there -- and that mixture is invisible to null-scanning, since the stale row is present rather than \\N.  Use this whenever settings changed and you want the output to mean exactly one run.")

## TODO: Support building on partial progress by checking what files already exist, or by writing a json file that is updated with state

args = parser.parse_args()

def ranks():
    yield from range(args.min_rank, args.max_rank + 1)

def signatures():
    for r in ranks():
        yield from [(r - n_minus, n_minus) for n_minus in range(r//2+1)]

def split_interval(a, b, m, prefix=()):
    n = b - a + 1
    q, r = divmod(n, m)
    bounds = [a]
    for i in range(m):
        bounds.append(bounds[-1] + q + (1 if i < r else 0))
    return [prefix + (bounds[i], bounds[i+1] - 1) for i in range(m)]

def Dbound(sig):
    c = args.min_disc_ratio
    C = args.max_disc_ratio
    r = sum(sig)
    if not args.nokthree and (12 <= r <= 20) and (sig[1] in [1,2]):
        min_det = c // (22 - r) + 1
        max_det = C // (22 - r)
    else:
        min_det = c // r + 1
        max_det = C // r
    return min_det, max_det

def build_genus_inputs():
    inputs = []
    for sig in signatures():
        min_det, max_det = Dbound(sig)
        inputs_sig = split_interval(min_det, max_det, args.jobs, sig)
        inputs.extend(inputs_sig)
    return inputs

# Calibrated per-genus fill cost in seconds (from a C=256 rank-vs-wall-clock calibration).
# Definite genera dominate, and their cost cliffs at rank 9 where the per-lattice loop
# saturates enum_timelimit; indefinite genera (SpinorRepresentatives, small class number)
# are uniformly cheap (~0.5s) at every rank.  We batch by summing this cost to a target so
# every batch is roughly equal wall-clock.
_DEF_COST = {1: 1.1, 2: 1.7, 3: 2.2, 4: 2.4, 5: 2.9, 6: 4.1, 7: 3.8, 8: 11.0,
             9: 221.0, 10: 460.0, 11: 384.0, 12: 693.0}
def genus_work(r, definite):
    if not definite:
        return 0.5
    if r in _DEF_COST:
        return _DEF_COST[r]
    # Rank >= 13 was never calibrated -- the C=256 census stops at 12 for lack of
    # determinant room -- so fall back to the timelimit, the only bound we have.  These
    # weights are also C=256-specific: at larger determinants automorphism groups shrink
    # and class numbers grow, so the rank 9-12 costs will not transfer to the C=32768
    # run.  Recalibrate against the actual census before trusting them there.
    return float(args.enum_timelimit)

_GENUS_COSTS = None
def measured_cost(r, nminus, det):
    """Cost of a genus from --genus-costs: its (rank, nminus, log2 det) stratum, else
    the nearest determinant stratum of the same rank and signature; None if neither."""
    global _GENUS_COSTS
    if not args.genus_costs:
        return None
    if _GENUS_COSTS is None:
        _GENUS_COSTS = defaultdict(dict)      # (rank, nminus) -> {log2det: seconds}
        with open(args.genus_costs) as F:
            for line in F:
                f = line.split()
                if len(f) == 4 and f[0].isdigit():
                    _GENUS_COSTS[(int(f[0]), int(f[1]))][int(f[2])] = float(f[3])
    by_det = _GENUS_COSTS.get((r, nminus))
    if not by_det:
        return None
    b = det.bit_length() - 1
    return by_det[min(by_det, key=lambda k: (abs(k - b), k))]

def build_enumeration_inputs(fname, rank=None):
    """Write fill batches (colon-joined labels, one batch per line) to fname.
    Genera are packed until their estimated cost reaches --batch-work-target (or
    --batch-size labels); a genus that alone reaches the target gets a batch of its
    own, so it never holds up cheaper genera queued behind it.  Batches are written
    most expensive first: GNU parallel starts them in file order, and a rank cannot
    finish before its longest batch (on the C=32768 calibration, rank-8 batches of 11
    serial genera ran for hours on 3 of 24 slots)."""
    T = args.batch_work_target
    batches = []                    # (estimated seconds, labels)
    for sig in signatures():
        r = sum(sig)
        if rank is not None and r != rank:
            continue
        nplus = sig[0]
        base = Path("genera_basic", str(r), str(nplus))
        if not base.is_dir():
            continue            # no genera enumerated for this signature
        definite = (r == nplus)
        W, labels = 0.0, []
        for genus in base.iterdir():
            w = measured_cost(r, r - nplus, int(genus.name.split(".")[2]))
            if w is None:
                w = genus_work(r, definite)
            if w >= T:
                batches.append((w, [genus.name]))
                continue
            W += w
            labels.append(genus.name)
            if W >= T or len(labels) >= args.batch_size:
                batches.append((W, labels))
                W, labels = 0.0, []
        if labels:
            batches.append((W, labels))
    batches.sort(key=lambda b: -b[0])
    with open(fname, "w") as Fout:
        for _, labels in batches:
            _ = Fout.write(":".join(labels) + "\n")
    return len(batches)

def build_tensor_inputs(fname):
    total = 0
    with open(fname, "w") as Fout:
        for sig in signatures():
            r = sum(sig)
            nplus = sig[0]
            d, D = Dbound(sig)
            # TODO: run_tensor.m only takes D right now
            _ = Fout.write(f"{r}.{nplus}:{D}\n")
            total += 1
    return total

def build_embedding_inputs(fname):
    """Write embed_jobs.txt (lines "task:labelsfile:labelsfile2", with "-" when a
    task has no second file) plus the per-task genus-label lists under embed_labels/.

    A genus qualifies when fill produced its genera_advanced record with
    1 <= class_number <= --embed-classlimit and its determinant is within the
    embedding cap (C/rank, with the K3 extension C/(22-rank) for the k3 lists);
    the Julia loader re-verifies completeness (every class file on disk) before
    using a genus, so a qualifying-but-partially-stored genus is skipped there.
    Each task's domain genera are split into det-sorted chunks of --embed-chunk,
    one job each (corank-1 jobs load the whole rank+1 pool as codomains); the
    expected jobtags go to embed_labels/expected_jobs.txt for the certificate
    derivation.  Returns (jobs to run, jobs expected): jobs already committed with
    the same inputs are not rerun."""
    C = args.embed_max_disc_ratio or args.max_disc_ratio
    tasks = [t.strip() for t in args.embed_tasks.split(",") if t.strip()]
    basic_fmt = Path("genera_basic.format").read_text().strip().split("|")
    combined_fmt = basic_fmt + Path("genera_advanced.format").read_text().strip().split("|")
    by_rank = defaultdict(list)            # rank -> [(label, nminus, is_even, disc)]
    for r in ranks():
        base = Path("genera_advanced", str(r))
        if not base.is_dir():
            continue
        for f in base.glob("*/*"):
            pieces = f.read_text().split("\n")[0].split("|")
            if len(pieces) != len(combined_fmt):
                continue
            rec = dict(zip(combined_fmt, pieces))
            cn = rec.get("class_number", "\\N")
            if cn in ("\\N", "None", "") or not (1 <= int(cn) <= args.embed_classlimit):
                continue
            by_rank[r].append((rec["label"], int(rec["nminus"]), rec["is_even"] == "T",
                               int(rec["disc_abs"])))
    labels_dir = Path("embed_labels")
    if labels_dir.is_dir():
        shutil.rmtree(labels_dir)
    labels_dir.mkdir()
    # Exclusion lists for the completeness certificates (oscar/refine_universality.py):
    # a definite genus inside the det window that is NOT part of the embedding
    # universe caps the certified range below its det.
    included = {r: {lab for lab, _, _, _ in recs} for r, recs in by_rank.items()}
    for r in ranks():
        base = Path("genera_basic", str(r), str(r))       # definite: nplus == rank
        if not base.is_dir():
            continue
        base_cap = C // r
        lines = []
        for f in base.glob("*"):
            pieces = f.read_text().split("\n")[0].split("|")
            rec = dict(zip(basic_fmt, pieces))
            det = abs(int(rec["disc_abs"]))
            if det > base_cap:
                continue
            lab = rec["label"]
            if lab in included.get(r, ()):
                continue
            reason = ("no_advanced" if not Path("genera_advanced", str(r), str(r), lab).exists()
                      else "classlimit")
            lines.append(f"{lab}|{det}|{reason}")
        if lines:
            (labels_dir / f"rank{r}_excluded.txt").write_text("\n".join(lines) + "\n")
    files = {}                             # (task list, rank) -> [(chunk file, min det)]
    pools = {}                             # rank -> full rank{r}_all.txt (corank-1 codomain pool)

    def write_chunks(stem, labs):
        # det-sorted chunks, so a chunk that does not finish caps the completeness
        # certificate (refine_universality.derive_certificates) just below its own
        # smallest det rather than voiding the whole rank
        labs = sorted(labs, key=lambda lab: (int(lab.split(".")[2]), lab))
        out = []
        for k in range(0, len(labs), args.embed_chunk):
            chunk = labs[k:k + args.embed_chunk]
            p = labels_dir / f"{stem}.c{k // args.embed_chunk + 1:04d}.txt"
            p.write_text("\n".join(chunk) + "\n")
            out.append((p, int(chunk[0].split(".")[2])))
        return out

    for r, recs in by_rank.items():
        base_cap = C // r
        labels = [lab for lab, _, _, d in recs if d <= base_cap]
        if labels:
            p = labels_dir / f"rank{r}_all.txt"
            p.write_text("\n".join(labels) + "\n")
            pools[r] = p
            files[("all", r)] = write_chunks(f"rank{r}_all", labels)
        if "k3" in tasks and args.embed_k3_min_rank <= r <= 21:
            k3_cap = max(base_cap, C // (22 - r)) if (12 <= r <= 20 and not args.nokthree) else base_cap
            k3_labels = [lab for lab, nminus, is_even, d in recs
                         if is_even and nminus in (1, 2) and d <= k3_cap]
            if k3_labels:
                files[("k3", r)] = write_chunks(f"rank{r}_k3", k3_labels)
    jobs = []                              # (task, chunk file, second file or "-", min det)
    for r in sorted(by_rank):
        if "overlat" in tasks:
            # overlattices are looked up in the whole rank, not just the chunk
            jobs += [("overlat", p, pools[r], m) for p, m in files.get(("all", r), [])]
        if "k3" in tasks:
            jobs += [("k3", p, "-", m) for p, m in files.get(("k3", r), [])]
        if "corank1" in tasks and (r + 1) in pools:
            jobs += [("corank1", p, pools[r + 1], m) for p, m in files.get(("all", r), [])]
    jobtag = lambda task, p: f"{task}_{p.stem}"
    (labels_dir / "expected_jobs.txt").write_text(
        "".join(f"{jobtag(t, p)}\t{m}\n" for t, p, _, m in jobs))
    # every job's inputs, for the pool's embed task (LMFDBLattices/src/embed.jl)
    (labels_dir / "jobs.tsv").write_text(
        "".join(f"{jobtag(t, p)}\t{t}\t{p}\t{p2}\t{args.embed_max_rows_per_domain}\n" for t, p, p2, _ in jobs))
    # Resume: the driver commits each finished job atomically to
    # inclusions/data/<jobtag>/ with the sha1 of its inputs in the manifest.  Keep
    # those, drop committed jobs this run does not expect (a different chunking or
    # det cap would otherwise duplicate rows in lat_inclusions.tbl), and do not
    # launch jobs that are already done -- each would pay ~30s of Oscar startup
    # only for the driver to skip it.
    expected = {jobtag(t, p) for t, p, _, _ in jobs}
    data = Path("inclusions", "data")
    if data.is_dir():
        for d in data.iterdir():
            if d.name not in expected:
                shutil.rmtree(d)
    total = done = 0
    with open(fname, "w") as Fout:
        for t, p, p2, _ in jobs:
            # the same digest as embed.jl's inputs_digest: the label lists and the row cap
            digest = hashlib.sha1(("\0".join(Path(f).read_text() for f in [p, p2] if f != "-")
                                   + f"\0rowcap={args.embed_max_rows_per_domain}").encode()).hexdigest()
            man = data / jobtag(t, p) / "manifest.txt"
            if man.is_file() and f"inputs={digest}" in man.read_text().split("\n")[0].split():
                done += 1
                continue
            Fout.write(jobtag(t, p) + "\n")
            total += 1
    if done:
        print(f"  {done} embedding jobs already committed with the same inputs; skipping them", flush=True)
    return total, len(jobs)

def julia_env():
    """A minimal environment for Julia subprocesses.  Running under `sage -python`,
    `import sage.all` exports Sage's compiler/library variables (CPATH, LIBRARY_PATH,
    LDFLAGS, PERL5LIB, SINGULAR_BIN_DIR, ...), which stall Oscar's startup in the
    child from ~25s to beyond 600s; with only the basics kept it loads normally."""
    keep = ("PATH", "HOME", "USER", "LOGNAME", "SHELL", "TMPDIR", "LANG", "LC_ALL", "TERM")
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env.update({k: v for k, v in os.environ.items() if k.startswith("JULIA_")})
    return env

def oscar_enumerate():
    """Enumerate the genera Oscar handles better (see --no-oscar-enum) on a pool of
    Oscar workers, before fill.  Records land in oscar_genus/; FillGenus reads them.
    Already finished records are kept (rerun without --clean resumes)."""
    labels = []
    for sig in signatures():
        r, nplus = sum(sig), sig[0]
        definite = (r == nplus)
        if (definite and r > args.oscar_definite_max_rank) or (not definite and r < args.oscar_indefinite_min_rank):
            continue
        base = Path("genera_basic", str(r), str(nplus))
        if base.is_dir():
            labels.extend(g.name for g in base.iterdir())
    if not labels:
        return
    Path("oscar_labels.txt").write_text("\n".join(labels) + "\n")
    if not check_oscar():
        indefinite = sum(1 for l in labels if l.split(".")[0] != l.split(".")[1])
        print(f"WARNING: Julia/Oscar environment not ready; Magma will enumerate all {len(labels)} "
              f"genera meant for Oscar, including {indefinite} indefinite ones of rank >= "
              f"{args.oscar_indefinite_min_rank}, where its spinor kernel can overcount classes", flush=True)
        return
    with timed(f"Enumerating {len(labels)} genera with Oscar"):
        subprocess.run(["julia", "--project=oscar", "--startup-file=no", "oscar/pool.jl",
                        "--task", "enumerate", "--labels-file", "oscar_labels.txt",
                        "--jobs", str(args.oscar_jobs or args.jobs), "--timeout", str(args.oscar_timeout),
                        "--memlimit", str(args.oscar_memlimit), "--batch", str(args.oscar_batch),
                        "--data-root", ".", "--out-root", "oscar_genus", "--ledger", "oscar_enum.ledger.tsv"]
                       + (["--retry-failed"] if args.retry_missing else []),
                       env=julia_env())

def check_oscar():
    """The embedding stage needs the pinned dev-Oscar project under oscar/ (see the
    README); a broken or missing environment skips the stage rather than failing it,
    and prints why so the skip is diagnosable."""
    try:
        # LMFDBLattices (oscar/LMFDBLattices) is precompiled with a workload the first
        # time it is loaded on a machine, ~5 min at low load; allow for that
        r = subprocess.run(["julia", "--project=oscar", "-e",
                            "using LMFDBLattices, Oscar; @assert isdefined(Oscar, :_primitive_embeddings_in_unimodular_safe); "
                            "println(\"OSCAR_VERSION \", Oscar.VERSION_NUMBER, \" julia \", VERSION)"],
                           capture_output=True, timeout=3600, env=julia_env())
        v = [l for l in r.stdout.decode(errors="replace").splitlines() if l.startswith("OSCAR_VERSION")]
        if v:
            provenance("oscar", version=v[0].split(" ", 1)[1])
        if r.returncode != 0:
            tail = r.stderr.decode(errors="replace").strip().splitlines()[-3:]
            print("  oscar preflight failed: " + " | ".join(tail), flush=True)
        return r.returncode == 0
    except FileNotFoundError:
        print("  oscar preflight failed: no `julia` on PATH", flush=True)
        return False
    except subprocess.TimeoutExpired:
        print("  oscar preflight failed: loading LMFDBLattices exceeded 3600s", flush=True)
        return False

def parallel(jobfile, joblog, prog_args, colsep=None, prog=None, jobs=None, timeout=None, env=None):
    cmd = ["parallel", "-j", str(jobs or args.jobs), "--joblog", joblog, "--eta", "-a", jobfile]
    if timeout is not None:
        cmd += ["--timeout", str(timeout)]
    if colsep is not None:
        cmd += ["--colsep", colsep]
    cmd += (prog or ["magma", "-b"]) + prog_args
    # GNU parallel returns the number of failed jobs; don't abort the pipeline on
    # per-job failures -- log a warning and continue (the joblog records which
    # jobs failed so they can be retried).
    result = subprocess.run(cmd, env=env)
    if result.returncode != 0:
        print(f"  WARNING: {result.returncode} job(s) failed (see {joblog}); continuing", flush=True)

def audit_fill(jobfile, missing_file):
    """Every genus fill attempted must leave a genera_advanced record -- even a genus a
    guard stopped (masslimit, enum timeout, sizelimit, timelimit) gets its genus-level
    row.  So a missing record means the genus was never written: a Magma error or a
    crash.  That is invisible to null-scanning (there is no row to contain \\N), so list
    those genera after each rank rather than leave them to be noticed downstream."""
    missing = []
    with open(jobfile) as F:
        for line in F:
            for label in line.strip().split(":"):
                if not label:
                    continue
                rank, nplus = label.split(".")[:2]
                if not Path("genera_advanced", rank, nplus, label).exists():
                    missing.append(label)
    if missing:
        Path(missing_file).write_text("\n".join(missing) + "\n")
        print(f"  WARNING: {len(missing)} genera have no genera_advanced record "
              f"(listed in {missing_file}; retry with run_fill_genus.sh)", flush=True)
    elif Path(missing_file).exists():
        Path(missing_file).unlink()

def fill_festi_veniani():
    """festi_veniani_index = [O(q_L) : image of O(L)].  Fill records the image orders in
    festi_veniani/ and leaves the column null; oscar/fix_festi_veniani.jl computes
    |O(q_L)| per genus (Magma has no orthogonal groups of torsion quadratic forms) and
    writes the index.  Chunks of genera run under GNU parallel; every job rewrites only
    finished genera, so a rerun is safe and a timed-out chunk leaves its rest null."""
    base = Path("festi_veniani")
    labels = sorted(f.name for f in base.rglob("*") if f.is_file() and not f.name.endswith(".tmp")) if base.is_dir() else []
    if not labels:
        return
    if not check_oscar():
        print("  WARNING: no working Oscar; festi_veniani_index stays null "
              "(run oscar/fix_festi_veniani.jl from lattices/ once Oscar works)", flush=True)
        return
    jobdir = Path("festi_veniani_jobs")
    if jobdir.exists():
        shutil.rmtree(jobdir)
    jobdir.mkdir()
    chunks = []
    for k in range(0, len(labels), args.fv_chunk):
        f = jobdir / f"chunk_{k // args.fv_chunk:05d}.txt"
        f.write_text("\n".join(labels[k:k + args.fv_chunk]) + "\n")
        chunks.append(str(f))
    Path("festi_veniani_jobs.txt").write_text("\n".join(chunks) + "\n")
    with timed(f"Festi-Veniani indices (Oscar): {len(labels)} genera in {len(chunks)} jobs"):
        parallel("festi_veniani_jobs.txt", "festi_veniani.joblog", ["oscar/fix_festi_veniani.jl", ".", "{1}"],
                 prog=["julia", "--project=oscar"], jobs=args.fv_jobs, timeout=args.fv_timelimit, env=julia_env())

def embed_units_summary():
    """Merge the committed embedding jobs' units.jsonl ledgers (one line per domain
    genus: task, domain, classes, seconds, rows, status ok/overflow/skipped) into
    embed_units.jsonl, the input for calibrating the stage, and print per-task totals
    and the overflowed domains (their rows are not in lat_inclusions.tbl)."""
    import json
    data = Path("inclusions", "data")
    if not data.is_dir():
        return
    units = []
    for f in sorted(data.glob("*/units.jsonl")):
        units += [json.loads(line) for line in f.read_text().splitlines() if line.strip()]
    with open("embed_units.jsonl", "w") as F:
        for u in units:
            F.write(json.dumps(u) + "\n")
    tot = defaultdict(lambda: defaultdict(float))
    for u in units:
        t = tot[u["task"]]
        t[u["status"]] += 1
        t["seconds"] += u["seconds"]
        t["rows"] += u["rows"] if u["status"] == "ok" else 0
    for task, t in sorted(tot.items()):
        print(f"  {task}: {int(t['ok'])} domain genera, {int(t['rows'])} rows, {t['seconds'] / 3600:.2f} h"
              + (f", {int(t['overflow'])} overflowed" if t["overflow"] else "")
              + (f", {int(t['skipped'])} not loadable" if t["skipped"] else ""), flush=True)
    over = [u["domain"] for u in units if u["status"] == "overflow"]
    if over:
        print(f"  WARNING: {len(over)} domain genera exceeded --embed-max-rows-per-domain "
              f"{args.embed_max_rows_per_domain} and are left out of lat_inclusions (e.g. {', '.join(over[:5])}; "
              f"all in embed_units.jsonl)", flush=True)

def theta_hecke_pass():
    """is_theta_distinguished for the genera fill left undecided.  SetHashes writes
    theta_hecke/<r>/<nplus>/<label> for each (the classes' Gram matrices and theta
    series, the adjacency matrix when fill computed it at the right prime); then
    run_theta_hecke.m decides them with the p-neighbour graph and the Hecke eigenspaces
    of the cusp space (theta_hecke.m) and writes theta_hecke_done/.  A job holds whole
    groups of genera sharing a cusp space, so each space is built once per run, and
    theta_hecke_spaces/ keeps the eigenspace data across runs (it depends only on the
    space).  Genera with a verdict are skipped on a rerun; the verdicts are patched into
    genera_advanced, where the field is still null."""
    base = Path("theta_hecke")
    if not base.is_dir():
        return
    groups = defaultdict(list)
    for f in base.rglob("*"):
        if not f.is_file() or ".tmp" in f.name:
            continue
        if not Path("theta_hecke_done", *f.parts[1:]).exists():
            with open(f) as F:
                groups[tuple(F.readline().split("|")[:4])].append(f.name)
    if groups:
        jobdir = Path("theta_hecke_jobs")
        if jobdir.exists():
            shutil.rmtree(jobdir)
        jobdir.mkdir()
        jobs, cur = [], []
        for key in sorted(groups, key=lambda k: -len(groups[k])):
            cur += groups[key]
            if len(cur) >= args.theta_hecke_chunk:
                jobs.append(cur)
                cur = []
        if cur:
            jobs.append(cur)
        names = []
        for k, labs in enumerate(jobs):
            f = jobdir / f"job_{k:05d}.txt"
            f.write_text("\n".join(labs) + "\n")
            names.append(str(f))
        Path("theta_hecke_jobs.txt").write_text("\n".join(names) + "\n")
        ngen = sum(len(v) for v in groups.values())
        with timed(f"Theta-Hecke pass: {ngen} genera, {len(groups)} cusp spaces, {len(jobs)} jobs"):
            parallel("theta_hecke_jobs.txt", "theta_hecke.joblog",
                     [f"timeout:={args.theta_hecke_timeout}", f"memlimit:={args.memlimit}", "labels:={1}", "run_theta_hecke.m"],
                     jobs=args.theta_hecke_jobs)
    patch_theta_hecke()

def patch_theta_hecke():
    """Write the theta-Hecke verdicts into genera_advanced (true = distinguished,
    false = some classes have equal theta series) where fill left is_theta_distinguished
    null.  A verdict contradicting fill's is a bug in one of them: reported, not written."""
    fmt = Path("genera_basic.format").read_text().strip().split("|") + Path("genera_advanced.format").read_text().strip().split("|")
    col = fmt.index("is_theta_distinguished")
    counts = defaultdict(int)
    for f in Path("theta_hecke_done").rglob("*") if Path("theta_hecke_done").is_dir() else []:
        if not f.is_file() or ".tmp" in f.name:
            continue
        status = f.read_text().split("|")[0]
        counts[status] += 1
        if status not in ("distinguished", "equal"):
            continue
        value = "true" if status == "distinguished" else "false"
        adv = Path("genera_advanced", *f.parts[1:])
        if not adv.exists():
            continue
        lines = adv.read_text().split("\n")
        fields = lines[0].split("|")
        if len(fields) != len(fmt):
            print(f"  ERROR: {adv}: {len(fields)} fields, expected {len(fmt)}", flush=True)
            continue
        if fields[col] == value:
            continue
        if fields[col] != "\\N":
            print(f"  ERROR: {f.name}: fill says is_theta_distinguished={fields[col]}, the theta-Hecke pass {value}", flush=True)
            continue
        fields[col] = value
        lines[0] = "|".join(fields)
        tmp = adv.with_name(adv.name + ".tmp")
        tmp.write_text("\n".join(lines))
        os.replace(tmp, adv)
        counts["patched"] += 1
    if counts:
        print("  theta-Hecke verdicts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), flush=True)

def clean_outputs():
    """Delete pipeline outputs so this run's data cannot be a mixture of several runs.

    Every stage writes per-file with Overwrite and never removes anything, so a genus or
    lattice this run does not produce silently keeps the file an earlier run wrote --
    e.g. a genus whose per-lattice loop overruns --enum-timelimit stores no lattices at
    all (FillGenus drops them rather than record a set that contradicts class_number),
    yet the lattices a previous, more generous run stored are still on disk.  Null-scans
    cannot see this (the stale row is present, not \\N), so coverage looks better than it
    is -- it previously made an audit blame connect for a gap that fill had created.
    """
    targets = ["genera_advanced", "genera_hash", "lattice_basic_data",
               "lattice_advanced_data", "lattice_hashes", "voronoi", "shortest",
               "orth_cache", "fill_stats", "connect_stats", "oscar_genus",   # stale cached genera would leak old counters/labels into provenance
               "inclusions", "embed_labels", "embed_pool", "refine_sidecar", "label_lookups", "label_lookups_done",
               "festi_veniani", "festi_veniani_jobs", "theta_hecke", "theta_hecke_done", "theta_hecke_jobs"]
    # theta_hecke_spaces/ is kept: eigenspace data of cusp spaces, a function of the space alone
    if not args.skip_list_genera:
        targets.append("genera_basic")   # regenerated below; otherwise it is this run's input
    removed = []
    for t in targets:
        p = Path(t)
        if p.is_dir():
            removed.append(f"{t} ({sum(1 for q in p.rglob('*') if q.is_file())} files)")
            shutil.rmtree(p)
    for f in ([Path("atomic_names"), Path("lat_inclusions.tbl")] + sorted(Path(".").glob("atomic_names_partial_*"))
              + sorted(Path(".").glob("fill_missing_*.txt"))
              + [Path("oscar_labels.txt"), Path("oscar_enum.ledger.tsv"), Path("lookup_labels.txt"),
                 Path("label_lookups.ledger.tsv"), Path("embed_units.jsonl")]):
        if f.exists():
            f.unlink()
            removed.append(f.name)
    print("Cleaned: " + (", ".join(removed) if removed else "(nothing to remove)"), flush=True)


def main():
    if args.retry_missing:
        args.skip_list_genera = True
    provenance("start")
    if args.clean:
        clean_outputs()
    if args.skip_embeddings and Path("refine_sidecar").is_dir():
        print("WARNING: --skip-embeddings, but a refine_sidecar/ directory from an "
              "earlier run exists; connect will consume it.  Use --clean (or delete "
              "it) if this run's data differs from the run that produced it.", flush=True)
    if not args.skip_list_genera:
        inputs = build_genus_inputs()
        total = len(inputs)
        start_time = time.monotonic()
        print(f"Listing genera ({total} tasks)...")

        # genus.py decorates the lister with @parallel(ncpus=128); rewrap the undecorated
        # function so the Sage worker pool honours -j instead (on a 32-core machine the
        # fixed 128 oversubscribed it 4x, and nothing else in the run bounded it).
        lister = sage_parallel(ncpus=args.jobs)(write_all_of_sig_between_genera_basic.func)
        for i, res in enumerate(lister(inputs), 1):
            elapsed = time.monotonic() - start_time
            avg_time = elapsed / i
            eta = avg_time * (total - i)
            percent = (i / total) * 100

            # Build the bar string
            bar = "#" * int(percent // 5) + "-" * (20 - int(percent // 5))

            # Update the same line in stdout
            sys.stdout.write(f"\r|{bar}| {percent:.1f}% ({i}/{total}) | ETA: {eta:.1f}s")
            sys.stdout.flush()

        print(f"\nDone listing genera! ({time.monotonic() - start_time:.1f}s)")

    gcount = build_enumeration_inputs("genus_jobs.txt")   # full list, used by the connect stages
    if not args.skip_enumerate_genera:
        # Unpatched Magma returns wrong classes for indefinite binary genera.
        # Require the success line: an error inside Magma can still exit 0.
        check = subprocess.run(["magma", "-b", "check_magma.m"], capture_output=True, text=True)
        print(check.stdout.strip(), flush=True)
        provenance("magma", check_magma=check.stdout.strip().splitlines()[-1:] or None)
        if "MAGMA_OK" not in check.stdout and not args.allow_unpatched_magma:
            sys.exit("Stopping: patch Magma (magma_patches/) or pass --allow-unpatched-magma.")
        if not args.no_oscar_enum:
            oscar_enumerate()
        # Work DOWNWARD by rank, with a barrier between ranks.  The orbit method
        # enumerates a definite genus by descending from a parent genus of rank+1
        # and smaller determinant, and FillGenus persists every finished definite
        # genus (with its label and per-lattice counters) to the on-disk orth
        # cache (orth_cache/).  Since a parent's determinant is at most half the
        # child's, every parent of a rank-r genus lies inside the rank-(r+1) box,
        # so completing rank r+1 first turns parent lookups into disk hits:
        # one-step descents against already-filled genera (which also makes the
        # ambient_genus / ambient_lattice provenance resolvable) instead of cold
        # recursive chains.  Within a rank, batches run in parallel as before;
        # only the ranks are sequenced.
        for r in sorted(set(ranks()), reverse=True):
            if args.retry_missing:
                missing = Path(f"fill_missing_{r}.txt")
                if not missing.exists():
                    continue
                todo = [l for l in missing.read_text().split() if l]
                Path(f"genus_jobs_{r}.txt").write_text("\n".join(todo) + "\n")   # one genus per job
                rcount = len(todo)
            else:
                rcount = build_enumeration_inputs(f"genus_jobs_{r}.txt", rank=r)
            if rcount == 0:
                continue
            with timed(f"Enumerating genera of rank {r} ({rcount} batches)"):
                # Pass the enumeration guards through to FillGenus so a few pathological
                # genera can't dominate wall-clock: skip enumerating past enum-masslimit,
                # store genus-level data only past enum-sizelimit, and cap per-genus
                # wall-clock at enum-timelimit.
                parallel(f"genus_jobs_{r}.txt", f"fill_{r}.joblog",
                         [f"timeout:={args.enum_timeout}", f"masslimit:={args.enum_masslimit}",
                          f"sizelimit:={args.enum_sizelimit}", f"timelimit:={args.enum_timelimit}",
                          f"adjlimit:={args.enum_adjlimit}",
                          f"perlattice:={args.enum_per_lattice_timeout}",
                          f"useorth:={0 if args.no_orth else 1}",
                          f"warmlimit:={args.warm_limit}",
                          f"stats:={1 if args.fill_stats else 0}",
                          f"memlimit:={args.memlimit}",
                          f"oscar:={'oscar_genus_disabled' if args.no_oscar_enum else 'oscar_genus'}",
                          "labels:={1}"],
                         # The wrapper adds done:=1 and checks that every label in the
                         # batch printed DONE or ERROR, so a Magma segfault (which can
                         # exit 0) shows up as ERROR lines and a failed job instead of
                         # silently missing genera.
                         prog=["./run_fill_genus.sh"])
            audit_fill(f"genus_jobs_{r}.txt", f"fill_missing_{r}.txt")
        if not args.skip_festi_veniani:
            fill_festi_veniani()
        if not args.skip_theta_hecke:
            theta_hecke_pass()

    if not args.skip_embeddings:
        # Julia/Oscar embedding stage (oscar/INTEGRATION_PLAN.md): overlattices,
        # corank-1 primitive embeddings and K3-embedding counts, assembled into
        # lat_inclusions.tbl.  Needs only the fill outputs, and runs last so the
        # slowest optional stage cannot delay the core tables.
        ecount, etotal = build_embedding_inputs("embed_jobs.txt")
        if etotal == 0:
            print("No embedding jobs in range; skipping embeddings", flush=True)
        elif ecount > 0 and not check_oscar():
            print("WARNING: Julia/Oscar environment not ready (setup: lattices/README.md); skipping embeddings", flush=True)
        else:
            if ecount > 0:
                # On the Oscar worker pool: each worker loads Oscar once and the
                # package's precompiled embedding code, instead of ~270 s of startup
                # and compilation per job (GNU parallel launching oscar/driver.jl).
                # A job over --embed-timelimit or --embed-memlimit has its worker
                # killed and replaced; its chunk stays uncommitted (see below).
                if Path("embed_pool").is_dir():
                    shutil.rmtree("embed_pool")
                with timed(f"Computing embeddings ({ecount} of {etotal} jobs)"):
                    subprocess.run(["julia", "--project=oscar", "--startup-file=no", "oscar/pool.jl",
                                    "--task", "embed", "--labels-file", "embed_jobs.txt",
                                    "--jobs", str(args.embed_jobs), "--timeout", str(args.embed_timelimit),
                                    "--memlimit", str(args.embed_memlimit), "--batch", "1",
                                    "--data-root", ".", "--out-root", "embed_pool",
                                    "--ledger", "embed.ledger.tsv"], env=julia_env())
            unfinished = [line.split("\t")[0] for line in
                          Path("embed_labels", "expected_jobs.txt").read_text().splitlines()
                          if line and not Path("inclusions", "data", line.split("\t")[0], "manifest.txt").is_file()]
            if unfinished:
                print(f"  WARNING: {len(unfinished)} of {etotal} embedding jobs did not finish "
                      f"(e.g. {', '.join(unfinished[:3])}); their domains are left out of "
                      f"lat_inclusions and of the certificates.  Rerunning without --clean "
                      f"(e.g. with --skip-list-genera --skip-enumerate-genera) resumes them.", flush=True)
            embed_units_summary()
            with timed("Writing lat_inclusions.tbl"):
                wt = subprocess.run(["julia", "--startup-file=no", "oscar/write_table.jl",
                                     "--inclusions-root", "inclusions",
                                     "--out", "lat_inclusions.tbl",
                                     "--lattice-root", "lattice_basic_data"], env=julia_env())
                if wt.returncode != 0:
                    # quarantine the partial table so the refine step (and anything
                    # else keyed on its existence) does not consume a bad file
                    if Path("lat_inclusions.tbl").exists():
                        Path("lat_inclusions.tbl").rename("lat_inclusions.tbl.failed")
                    print("  WARNING: lat_inclusions writer failed "
                          "(partial output in lat_inclusions.tbl.failed); "
                          "sidecar refinement will be skipped", flush=True)

    # Bottom-up bitset pass over the embedding poset (oscar/refine_universality.py):
    # give the rank-n lattices distinct powers of two, OR the masks up the corank-1
    # and overlattice edges, then record per-lattice sidecars for connect: db-relative
    # heuristic values (all-ones = n-universal over the coverage window, mask equals
    # the genus OR = n-regular), certified universality upper bounds, and regularity
    # candidate counterexamples.  Certificates come from the drivers' completion
    # manifests plus the exclusion lists: a cleanly finished job certifies its whole
    # det window, making closure-absence SOUND for corank-1 pairs inside it (their
    # chains cannot leave the window, since overlattice steps only lower the det), so
    # those candidates skip the Magma confirmation that connect runs for the rest.
    if not args.skip_embeddings and Path("lat_inclusions.tbl").exists():
        with timed("Deriving embedding-poset sidecars for connect"):
            sys.path.insert(0, "oscar")
            import refine_universality as _refine
            _lats = _refine.read_meta_lattice_dir("lattice_basic_data")
            _edges = _refine.read_inclusions("lat_inclusions.tbl", _lats)
            _C = args.embed_max_disc_ratio or args.max_disc_ratio
            _flags, _cands, _mmd = _refine.closure_flags(_lats, _edges, disc_ratio=_C)
            _ranks = sorted({m.rank for m in _lats.values() if m.definite})
            _certs = _refine.derive_certificates("embed_labels", "inclusions/data", _C, _ranks)
            _resolved = _refine.resolve_candidates(_cands, _lats, _certs)
            if Path("refine_sidecar").is_dir():
                shutil.rmtree("refine_sidecar")
            nsc = _refine.write_sidecars("refine_sidecar", _lats, _flags, _mmd, _resolved, _certs)
            Path("regularity_candidates.txt").write_text(
                "".join(f"{l}\t{n}\t{s0}\n" for l, n, s0 in _cands))
            _np = sum(1 for v in _resolved.values() for c in v if c[2] == "p")
            print(f"  wrote {nsc} sidecar records; {len(_cands)} candidates "
                  f"({_np} certificate-proven)", flush=True)

    if not args.skip_tensor:
        tcount = build_tensor_inputs("tensor_jobs.txt")
        with timed(f"Computing tensor products ({tcount} tasks)"):
            parallel("tensor_jobs.txt", "tensor.joblog",
                     ["sig:={1}", "Dbound:={2}", "run_tensor.m"], colsep=":")

    if not args.skip_names:
        with timed("Computing basic names"):
            # Cap the named-lattice scalings at the determinant we actually enumerated
            # (otherwise NameAtomicLattices scans scalings up to its default 32768,
            # which for low ranks is tens of thousands of wasted FindLabel calls), and
            # split the catalog across workers -- naming is otherwise the only serial
            # stage.  Each worker writes atomic_names_partial_<k>; merge them keeping
            # the lowest-priority name per label (the same collision rule the intrinsic
            # applies in-process).
            with open("name_jobs.txt", "w") as Fout:
                for w in range(args.jobs):
                    Fout.write(f"{w}:{args.jobs}\n")
            parallel("name_jobs.txt", "names.joblog",
                     [f"DETCAP:={args.max_disc_ratio}", "WORKER:={1}", "NWORKERS:={2}", "run_basic_names.m"],
                     colsep=":")
            best = {}   # label -> (priority, name)
            for w in range(args.jobs):
                partial = Path(f"atomic_names_partial_{w}")
                if not partial.exists():
                    continue
                for line in partial.read_text().splitlines():
                    if not line:
                        continue
                    label, name, prio = line.rsplit("|", 2)
                    # match NameAtomicLattices' collision rule: (priority, name),
                    # lower priority wins, ties to the lexicographically smaller name
                    key = (int(prio), name)
                    if label not in best or key < best[label]:
                        best[label] = key
                partial.unlink()
            with open("atomic_names", "w") as Fout:
                Fout.write("\n".join(f"{label}|{name}" for label, (_, name) in best.items()))
            print(f"Named {len(best)} atomic lattices (merged from {args.jobs} workers).", flush=True)

    if not args.skip_connect:
        # Connect in two passes.  A decomposable lattice derives its geometric fields
        # (eutaxy, covering radius, deep holes, ...) from its orthogonal factors --
        # strictly lower-rank lattices in other genera whose data is written during
        # connect.  A single flat pass that processed a composite before its factors
        # silently dropped those fields to \N.  So: pass 1 ("produce") runs the full
        # connect for every genus but skips the decomposable derivations, which writes
        # every indecomposable factor's data to disk; pass 2 ("consume") then derives
        # the decomposable fields with all factor data guaranteed present.  Both passes
        # are flat and fully parallel; the only barrier is between them.
        with timed("Connecting genera (pass 1: produce)"):
            parallel("genus_jobs.txt", "connect1.joblog",
                     ["labels:={1}", "phase:=1", f"memlimit:={args.memlimit}",
                      f"stats:={1 if args.connect_stats else 0}", "run_connect_genus.m"])
        with timed("Connecting genera (pass 2: consume)"):
            parallel("genus_jobs.txt", "connect2.joblog",
                     ["labels:={1}", "phase:=2", f"memlimit:={args.memlimit}",
                      f"stats:={1 if args.connect_stats else 0}", "run_connect_genus.m"])
        # Lookups FindLabel could not decide -- an indefinite genus of rank >= 3 with
        # several classes, or a Magma overrun hashing a lattice or testing isometry --
        # are recorded in label_lookups/ by connect; the Oscar pool's resolve task
        # tests them with is_isometric and patches the advanced records.
        lookups = Path("label_lookups")
        sources = sorted(f.name for f in lookups.rglob("*") if f.is_file()) if lookups.is_dir() else []
        if sources:
            Path("lookup_labels.txt").write_text("\n".join(sources) + "\n")
            if not check_oscar():
                print(f"WARNING: Julia/Oscar environment not ready; {len(sources)} lattices keep "
                      "null labels connect could not decide (label_lookups/)", flush=True)
            else:
                # resolution records go to a fresh directory: the patched columns are
                # rewritten by every connect run, so earlier records do not carry over
                if Path("label_lookups_done").is_dir():
                    shutil.rmtree("label_lookups_done")
                with timed(f"Resolving {len(sources)} lattices' deferred label lookups (Oscar)"):
                    subprocess.run(["julia", "--project=oscar", "--startup-file=no", "oscar/pool.jl",
                                    "--task", "resolve", "--labels-file", "lookup_labels.txt",
                                    "--jobs", str(args.oscar_jobs or args.jobs), "--timeout", str(args.oscar_timeout),
                                    "--memlimit", str(args.oscar_memlimit), "--batch", str(args.oscar_batch),
                                    "--data-root", ".", "--out-root", "label_lookups_done",
                                    "--ledger", "label_lookups.ledger.tsv"], env=julia_env())

if __name__ == "__main__":
    try:
        main()
        provenance("end", status="ok", stages=STAGE_TIMES)
    except SystemExit as e:
        provenance("end", status=f"exit {e.code}", stages=STAGE_TIMES)
        raise
    except BaseException as e:
        provenance("end", status=f"{type(e).__name__}: {e}", stages=STAGE_TIMES)
        raise
