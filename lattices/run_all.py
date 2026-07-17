#!/usr/bin/env -S sage -python

import sys
import time
import shutil
import argparse
import subprocess
import contextlib
from pathlib import Path
# Fully initialise Sage: under "sage -python" the lazy imports are not pre-resolved,
# which breaks genus enumeration in two ways on Sage 10.7 -- a missing
# laurent_polynomial_ring attribute (via discriminant_group(...).normal_form()),
# and a sage.symbolic.function circular import (via the genus mass computation).
import sage.all  # noqa: F401
from genus import write_all_of_sig_between_genera_basic

@contextlib.contextmanager
def timed(label):
    print(f"=== {label} ...", flush=True)
    t0 = time.monotonic()
    yield
    print(f"=== {label}: done in {time.monotonic() - t0:.1f}s", flush=True)

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
parser.add_argument("--batch-size", type=int, default=32, help="Backstop cap on genera per batch, regardless of work target.  Mainly bounds the per-job command-line length (one label each); for cheap genera the work target rarely fills a batch, so this keeps those batches from growing unwieldy.")

# Don't store too many lattices
parser.add_argument("--enum-masslimit", type=int, default=1000, help="If the mass of a genus is larger than this threshold, don't even try to enumerate lattices within")
parser.add_argument("--enum-per-lattice-timeout", type=int, default=0, help="PER-OPERATION wall-clock cap (seconds) applied INDEPENDENTLY to each of a lattice's own computations -- canonical form, automorphism group, primal theta, dual theta -- so a single lattice can spend up to ~4x this value plus unbounded minimum/kissing work; it is NOT a per-lattice aggregate.  The aggregate wall-clock guard is the per-genus --enum-timelimit, not this option.  0 (the default) instead splits enum-timeout across the genus, giving each operation enum-timeout/class_number.  That split looks like starvation -- at enum-timeout 300 a class-number-55 genus allows each operation only 6s, and CanonicalForm alone measures 8-66s at rank 12, so ~1%% of lattices get \\N in those columns -- but it keeps the genus finishing under --enum-timelimit and storing its lattices.  Setting a fixed budget here loosens that, and since a genus overrunning --enum-timelimit has ALL of its lattices discarded, it can convert partial data into no data unless --enum-timelimit is raised to match.  Measured on a rank-12 class-number-55 genus: the default stored 55 lattices (4 with a null canonical_gram) in 813s, while perlattice=120 at timelimit=900 stored 0 lattices in 1007s.  Erasing that ~1%% of nulls costs roughly 8x the runtime, so raise this only deliberately, together with --enum-timelimit.")
parser.add_argument("--enum-timeout", type=int, default=300, help="Seconds allowed for the genus-representatives computation of a single genus.  This was previously not passed through at all, so it silently fell back to run_fill_genus.m's 60s default and 1.8%% of definite genera (all rank 11-12) failed to enumerate entirely -- no class number and no lattices, a hole that null-scanning cannot see since the genus is simply absent from the per-lattice tables.  The computation is hard-killed past this, so it still bounds the worst case.")
parser.add_argument("--enum-timelimit", type=int, default=900, help="Maximum wall-clock seconds for a genus's per-lattice loop.  Overrunning it does NOT keep partial data: FillGenus discards every lattice of that genus (a partial set would contradict the recorded class_number) and keeps only the genus-level record.  So this is a coverage switch, not a speed knob -- and because it is wall-clock, it makes database contents depend on machine load, which is why --enum-masslimit is the better place to express a real coverage policy.  Calibrated at C=256: definite rank 9-12 genera need 97-1005s to COMPLETE (medians 221/460/384/693s by rank), so 900 completes ~97%% of them; 300 (the old default) silently discarded every lattice of most rank 10-12 genera, and 60 discarded essentially all of rank 9-12.")
parser.add_argument("--enum-sizelimit", type=int, default=1000, help="For genera with class number larger than this, do not store individual lattices within the genus")
parser.add_argument("--enum-adjlimit", type=int, default=0, help="Work budget for the adjacency (Hecke) matrix: skip it when the estimated work (~class_number * sum_p p^(rank-2) over the Hecke primes p) exceeds this; 0 = no budget, always attempt it.  Defaults to 0 because the estimate mis-prices reality badly -- at rank 7 the sum_p p^5 term crosses a 20000 budget at class number 6, yet a whole rank-7 genus fills in ~4s, so a budget of 20000 silently dropped Hecke data for 65 rank-7 genera (4.4%% of definite genera overall) that were never expensive.  AdjacencyMatrix now runs under a hard-killed TimeoutCall, which bounds it by measured time instead of a bad model; prefer that.  Set this only if you specifically want to skip adjacency without even trying.")
parser.add_argument("--no-orth", action="store_true", help="Disable the orbit-method genus enumeration (use the p-neighbour walk only)")
parser.add_argument("--warm-limit", type=int, default=300, help="Wall-clock budget (seconds) for each fill job's batch cache pre-warm.  The pre-warm runs in the parent Magma process (so forked per-genus enumerations inherit the warmed cache) and therefore cannot be bounded by TimeoutCall; past this budget it aborts at a safe point and the job continues with a partial cache plus the on-disk tier.")

# Skip stages
parser.add_argument("--skip-list-genera", action="store_true", help="Assume that genera have already been listed")
parser.add_argument("--skip-enumerate-genera", action="store_true", help="Assume that genera have already been enumerated (implies skip-list-genera)")
parser.add_argument("--skip-embeddings", action="store_true", help="Do not compute primitive and K3 embedding data")
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

def build_enumeration_inputs(fname, rank=None):
    W = 0.0
    T = args.batch_work_target
    labels = []
    total = 0
    with open(fname, "w") as Fout:
        for sig in signatures():
            r = sum(sig)
            if rank is not None and r != rank:
                continue
            nplus = sig[0]
            base = Path("genera_basic", str(r), str(nplus))
            if not base.is_dir():
                continue            # no genera enumerated for this signature
            definite = (r == nplus)
            for genus in base.iterdir():
                W += genus_work(r, definite)
                labels.append(genus.name)
                if W >= T or len(labels) >= args.batch_size:
                    _ = Fout.write(":".join(labels) + "\n")
                    W = 0.0
                    labels = []
                    total += 1
        if labels:
            _ = Fout.write(":".join(labels) + "\n")
            total += 1
    return total

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

def parallel(jobfile, joblog, magma_args, colsep=None):
    cmd = ["parallel", "-j", str(args.jobs), "--joblog", joblog, "--eta", "-a", jobfile]
    if colsep is not None:
        cmd += ["--colsep", colsep]
    cmd += ["magma", "-b"] + magma_args
    # GNU parallel returns the number of failed jobs; don't abort the pipeline on
    # per-job failures -- log a warning and continue (the joblog records which
    # jobs failed so they can be retried).
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"  WARNING: {result.returncode} job(s) failed (see {joblog}); continuing", flush=True)

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
               "orth_cache"]   # stale cached genera would leak old counters/labels into provenance
    if not args.skip_list_genera:
        targets.append("genera_basic")   # regenerated below; otherwise it is this run's input
    removed = []
    for t in targets:
        p = Path(t)
        if p.is_dir():
            removed.append(f"{t} ({sum(1 for q in p.rglob('*') if q.is_file())} files)")
            shutil.rmtree(p)
    for f in [Path("atomic_names")] + sorted(Path(".").glob("atomic_names_partial_*")):
        if f.exists():
            f.unlink()
            removed.append(f.name)
    print("Cleaned: " + (", ".join(removed) if removed else "(nothing to remove)"), flush=True)


def main():
    if args.clean:
        clean_outputs()
    if not args.skip_list_genera:
        inputs = build_genus_inputs()
        total = len(inputs)
        start_time = time.monotonic()
        print(f"Listing genera ({total} tasks)...")

        for i, res in enumerate(write_all_of_sig_between_genera_basic(inputs), 1):
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
                          "labels:={1}", "run_fill_genus.m"])

    if not args.skip_embeddings:
        print("Finding lattice embeddings (TODO: Oscar embedding code)")

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
                     ["labels:={1}", "phase:=1", "run_connect_genus.m"])
        with timed("Connecting genera (pass 2: consume)"):
            parallel("genus_jobs.txt", "connect2.joblog",
                     ["labels:={1}", "phase:=2", "run_connect_genus.m"])

if __name__ == "__main__":
    main()
