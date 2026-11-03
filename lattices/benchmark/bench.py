#!/usr/bin/env python3
"""Benchmark Magma against Oscar on the lattice quantities both can compute.

Run from lattices/ (plain python3; Sage is not needed):

  1. sample: draw genera from a census tree and write a cases file
       python3 benchmark/bench.py sample --data-root . --ranks 2-16 --per-stratum 3 \\
           -o bench_cases.tsv
     Strata are (rank, definite or indefinite, floor(log2 det)), with the same number
     of draws from each, as in calibrate_fill.py.  A case is one genus: its
     representative from genera_basic (or, if the genus is filled, its class .1 from
     lattice_basic_data, with class .2 as a non-isometric partner), a copy of it in a
     random basis (for the isometry tests) and a theta-series precision chosen so that
     about --theta-vectors vectors are enumerated.

  2. run: run every task on every case in both systems
       python3 benchmark/bench.py run bench_cases.tsv -o bench_results -j 8 --timeout 600
     -j servers per system, Magma and Oscar at the same time, so that both see the
     same machine load.  Results are appended to bench_results/results.tsv as they
     arrive; rerunning the same command resumes (finished jobs are skipped).  Within
     each (system, task, variant, rank, definite/indefinite) group the cases run in
     increasing determinant, and after --max-timeouts timeouts the rest of the group
     is recorded as skipped.

  3. report: compare
       python3 benchmark/bench.py report bench_results
     writes bench_results/report.md (per task and rank: success counts, median times,
     the Oscar/Magma time ratio, value mismatches) and mismatches.tsv.

The tasks (same names in bench_magma.m and bench_oscar.jl, which hold the actual
calls) are what both systems can compute that the pipeline stores: class_number,
mass, aut_order, aut_label, minimum, kissing, successive_minima, root_system, theta,
orth_decomp, isometric_same, isometric_other, disc_image, and (Oscar only) oq_order,
|O(q_L)| for the Festi-Veniani index.  Variant "default" is the
pipeline's Magma code path and Oscar's documented function; other variants are
alternatives (Magma's stock intrinsics, Oscar's enumeration algorithms) run with
--variants all.  Not covered, because only Magma has them: canonical forms, Hecke
(adjacency) matrices, hashes, universality, Voronoi cells and t-designs.
"""
import argparse
import json
import math
import os
import random
import select
import shutil
import signal
import socket
import statistics
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
LATTICES = HERE.parent
sys.path.insert(0, str(LATTICES))
from calibrate_fill import iter_census_labels, parse_label  # noqa: E402

TASKS = ["class_number", "mass", "aut_order", "aut_label", "minimum", "kissing",
         "successive_minima", "root_system", "theta", "orth_decomp", "isometric_same",
         "isometric_other", "disc_image", "oq_order"]
# Variants other than "default", per system (must match the servers).
VARIANTS = {
    "magma": {"class_number": ["stock"], "aut_order": ["stock"], "orth_decomp": ["stock"]},
    "oscar": {"class_number": ["orbit", "random"], "aut_order": ["order_only"]},
}
SYSTEMS = ["magma", "oscar"]
CASE_COLUMNS = ["id", "rank", "nplus", "det", "class_number", "mass", "theta_prec",
                "gram", "gram_iso", "gram_other"]
RESULT_COLUMNS = ["system", "case", "task", "variant", "status", "K", "wall_us", "cpu_us",
                  "extra", "value", "finished", "job_s"]   # job_s: the whole job, as bench.py saw it
# statuses that count against a group's --max-timeouts (a crash is usually the memory limit)
HARD_FAIL = {"timeout", "oom", "crash"}


# ---------------------------------------------------------------------------
# sample

def read_format(path):
    return Path(path).read_text().splitlines()[0].strip().split("|")


def label_path(root, folder, label):
    p = label.split(".")
    return Path(root, folder, p[0], p[1], label)


def read_record(path, fmt):
    fields = path.read_text().splitlines()[0].split("|")
    if len(fields) != len(fmt):
        return None
    return dict(zip(fmt, fields))


def parse_braces(s):
    return [int(x) for x in s.strip().strip("{}").split(",") if x.strip()]


def random_unimodular(n, rng):
    """A random matrix in GL_n(Z): elementary row operations, then a signed permutation."""
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(3 * n):
        i, j = rng.sample(range(n), 2) if n > 1 else (0, 0)
        if i == j:
            continue
        s = rng.choice((-1, 1))
        U[i] = [a + s * b for a, b in zip(U[i], U[j])]
    rng.shuffle(U)
    signs = [rng.choice((-1, 1)) for _ in range(n)]
    return [[x * e for x in row] for row, e in zip(U, signs)]


def transform(G, U):
    n = len(U)
    UG = [[sum(U[i][k] * G[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return [[sum(UG[i][k] * U[j][k] for k in range(n)) for j in range(n)] for i in range(n)]


def det_int(G):
    """The determinant of an integer matrix (exact Gaussian elimination over Q)."""
    from fractions import Fraction
    A = [[Fraction(x) for x in row] for row in G]
    n, d = len(A), Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return int(d)


def theta_prec(n, det, target):
    """Largest T in [2, 150] with about `target` vectors of norm <= T, by the Gaussian
    heuristic vol(ball of radius sqrt(T)) / sqrt(det)."""
    vn = math.pi ** (n / 2) / math.gamma(n / 2 + 1)
    T = (target * math.sqrt(det) / vn) ** (2 / n)
    return max(2, min(150, int(T)))


def cmd_sample(a):
    rng = random.Random(a.seed)
    lo, _, hi = a.ranks.partition("-")
    ranks = set(range(int(lo), int(hi or lo) + 1))
    kinds = set(a.kinds.split(","))
    pop = Counter()
    res = defaultdict(list)
    for label in iter_census_labels(a.data_root):
        rank, nplus, det = parse_label(label)
        if rank not in ranks or det < a.min_det or (a.max_det and det > a.max_det):
            continue
        kind = "def" if nplus == rank else ("indef" if 0 < nplus < rank else None)
        if kind not in kinds:
            continue
        st = (rank, kind, int(math.log2(det)) if det > 0 else 0)
        pop[st] += 1
        k = pop[st]
        if k <= a.per_stratum:
            res[st].append(label)
        else:
            j = rng.randrange(k)
            if j < a.per_stratum:
                res[st][j] = label
    root = a.data_root
    gfmt = read_format(Path(root, "genera_basic.format"))
    afmt = gfmt + read_format(Path(root, "genera_advanced.format")) \
        if Path(root, "genera_advanced.format").exists() else None
    lfmt = read_format(Path(root, "lat_basic.format")) if Path(root, "lat_basic.format").exists() else None
    rows = []
    for st in sorted(res):
        for label in sorted(res[st]):
            rank, nplus, det = parse_label(label)
            rec = read_record(label_path(root, "genera_basic", label), gfmt)
            if rec is None:
                print(f"  {label}: malformed genera_basic record, skipped", file=sys.stderr)
                continue
            gram = parse_braces(rec["rep"])
            other = None
            cn = "-"
            if afmt and label_path(root, "genera_advanced", label).exists():
                adv = read_record(label_path(root, "genera_advanced", label), afmt)
                if adv and adv.get("class_number") not in (None, "\\N", "None", ""):
                    cn = adv["class_number"]
            if lfmt:
                l1, l2 = (label_path(root, "lattice_basic_data", f"{label}.{i}") for i in (1, 2))
                if l1.exists():
                    r1 = read_record(l1, lfmt)
                    if r1:
                        gram = parse_braces(r1["gram"])
                    if l2.exists():
                        r2 = read_record(l2, lfmt)
                        if r2:
                            other = parse_braces(r2["gram"])
            if len(gram) != rank * rank:
                print(f"  {label}: gram has {len(gram)} entries, skipped", file=sys.stderr)
                continue
            G = [gram[i * rank:(i + 1) * rank] for i in range(rank)]
            Giso = transform(G, random_unimodular(rank, rng))
            assert det_int(Giso) == det_int(G), label
            mass = "-"
            if rec.get("mass") not in (None, "\\N", "None", ""):
                num, den = parse_braces(rec["mass"])
                mass = f"{num}/{den}" if den != 1 else str(num)    # as both systems print it
            rows.append({
                "id": label, "rank": rank, "nplus": nplus, "det": det, "class_number": cn,
                "mass": mass,
                "theta_prec": theta_prec(rank, det, a.theta_vectors) if nplus == rank else 0,
                "gram": ",".join(map(str, gram)),
                "gram_iso": ",".join(str(x) for row in Giso for x in row),
                "gram_other": ",".join(map(str, other)) if other else "-",
            })
    with open(a.output, "w") as F:
        F.write("\t".join(CASE_COLUMNS) + "\n")
        for r in rows:
            F.write("\t".join(str(r[c]) for c in CASE_COLUMNS) + "\n")
    with open(Path(a.output).with_suffix(".population.tsv"), "w") as P:
        P.write("rank\tkind\tlog2det\tpopulation\tsampled\n")
        for st in sorted(pop):
            P.write("%d\t%s\t%d\t%d\t%d\n" % (*st, pop[st], len(res[st])))
    print(f"wrote {len(rows)} cases from {sum(pop.values())} genera in {len(pop)} strata to {a.output}")


def read_cases(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    cases = {}
    for line in lines[1:]:
        if line:
            r = dict(zip(header, line.split("\t")))
            r["rank"], r["nplus"], r["det"] = int(r["rank"]), int(r["nplus"]), int(r["det"])
            r["kind"] = "def" if r["nplus"] == r["rank"] else "indef"
            cases[r["id"]] = r
    return cases


# ---------------------------------------------------------------------------
# run

def julia_env():
    """A minimal environment for Julia (as run_all.py's julia_env: Sage's compiler
    variables stall Oscar's startup), single-threaded like Magma."""
    keep = ("PATH", "HOME", "USER", "LOGNAME", "SHELL", "TMPDIR", "LANG", "LC_ALL", "TERM")
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env.update({k: v for k, v in os.environ.items() if k.startswith("JULIA_")})
    env.update(JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    return env


def rss_bytes(pid):
    try:
        with open(f"/proc/{pid}/status") as F:
            for line in F:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    return 0


class ServerDied(Exception):
    pass


class Server:
    """One Magma or Julia server process speaking the JOB / @@RESULT line protocol
    (see bench_magma.m).  Lines that are not protocol lines go to the log."""

    def __init__(self, system, slot, a, out):
        self.system, self.slot, self.a = system, slot, a
        self.proc = None
        self.buf = b""
        self.workdir = out / "work" / f"{system}-{slot}"
        self.logpath = out / "logs" / f"{system}-{slot}.log"
        self.ready = None

    def launch(self):
        a = self.a
        if self.workdir.exists():
            shutil.rmtree(self.workdir)      # no orth_cache/ left over from a previous server
        self.workdir.mkdir(parents=True)
        if getattr(self, "log", None):
            self.log.close()
        self.log = open(self.logpath, "a")
        self.log.write(f"=== start {time.strftime('%F %T')}\n")
        self.log.flush()
        cases = str(Path(a.cases).resolve())
        if self.system == "magma":
            cmd = [a.magma, "-b", f"cases:={cases}", f"spec:={LATTICES / 'lattices.spec'}",
                   f"memlimit:={a.memlimit}", str(HERE / "bench_magma.m")]
            env = None
        else:
            cmd = [a.julia, f"--project={LATTICES / 'oscar'}", "-t", "1", "--gcthreads=1",
                   "--startup-file=no", str(HERE / "bench_oscar.jl"), cases]
            env = julia_env()
        self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=self.log, cwd=self.workdir, env=env,
                                     start_new_session=True)
        self.buf = b""

    def wait_ready(self):
        a = self.a
        line = self._read_line("@@READY", a.startup_timeout, watch_memory=False)
        if line is None:
            self.kill()
            raise ServerDied(f"{self.system} server did not start within {a.startup_timeout}s; see {self.logpath}")
        f = line.split("\t")
        self.ready = {"load_s": float(f[1]), "warmup_s": float(f[2]), "version": f[3]}

    def kill(self):
        proc, self.proc = self.proc, None
        if proc is not None:
            try:
                os.killpg(proc.pid, signal.SIGKILL)   # the server and its children
            except (ProcessLookupError, PermissionError):
                pass
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                pass

    def _read_line(self, prefix, timeout, watch_memory=True):
        """The next line containing prefix, from the prefix on; None on timeout, "@@OOM"
        past the memory limit, ServerDied if the process exits or is killed."""
        deadline = time.time() + timeout
        proc = self.proc                   # kill() from another thread clears self.proc
        if proc is None:
            raise ServerDied(f"{self.system} server was killed")
        fd = proc.stdout.fileno()
        limit = self.a.memlimit * 2**30 if (watch_memory and self.system == "oscar" and self.a.memlimit) else 0
        while True:
            while b"\n" in self.buf:
                raw, self.buf = self.buf.split(b"\n", 1)
                line = raw.decode(errors="replace")
                i = line.find(prefix)      # progress output can leave a "\r..." before it
                if i >= 0:
                    if i > 0:
                        self.log.write(line[:i] + "\n")
                    return line[i:]
                self.log.write(line + "\n")
            self.log.flush()
            remaining = deadline - time.time()
            if remaining <= 0:
                return None
            r, _, _ = select.select([fd], [], [], min(remaining, 1.0))
            if r:
                chunk = os.read(fd, 1 << 16)
                if not chunk:
                    raise ServerDied(f"{self.system} server exited (code {proc.wait()})")
                self.buf += chunk
            elif limit and rss_bytes(proc.pid) > limit:
                return "@@OOM"

    def request(self, case, task, variant):
        a = self.a
        msg = f"JOB\t{case}\t{task}\t{variant}\t{a.timeout}\t{a.min_time}\t{a.max_reps}\n"
        try:
            self.proc.stdin.write(msg.encode())
            self.proc.stdin.flush()
        except (AttributeError, OSError) as e:          # killed, or exited
            raise ServerDied(f"{self.system} server is gone ({e})")
        # Magma bounds each job itself (HardKill TimeoutCall); this is only a backstop.
        # Julia cannot, so for Oscar this is the timeout.
        grace = 120 if self.system == "magma" else 0
        line = self._read_line("@@RESULT", a.timeout + grace)
        if line is None or line == "@@OOM":
            self.kill()
            return {"status": "timeout" if line is None else "oom"}
        f = line.split("\t")
        if f[1:4] != [case, task, variant]:
            raise ServerDied(f"out-of-sequence reply {f[1:4]} to {[case, task, variant]}")
        return dict(zip(["status", "K", "wall_us", "cpu_us", "extra", "value"], f[4:10]))


def load_results(path):
    """The last record of each (system, case, task, variant)."""
    out = {}
    if path.exists():
        lines = path.read_text().splitlines()
        if lines:
            header = lines[0].split("\t")
            for line in lines[1:]:
                f = line.split("\t")
                if len(f) == len(header):
                    r = dict(zip(header, f))
                    out[(r["system"], r["case"], r["task"], r["variant"])] = r
    return out


def cmd_run(a):
    out = Path(a.output)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    cases = read_cases(a.cases)
    systems = a.systems.split(",")
    tasks = TASKS if a.tasks == "all" else a.tasks.split(",")
    unknown = set(tasks) - set(TASKS)
    if unknown:
        sys.exit(f"unknown tasks {sorted(unknown)}; choose from {TASKS}")
    respath = out / "results.tsv"
    done = {} if a.redo else load_results(respath)

    # Jobs, grouped by (system, task, variant, rank, kind) in increasing determinant;
    # dispatched round-robin over the groups (first the easiest case of every group,
    # then the second easiest, ...), so a group's timeouts are seen before its hard
    # cases are reached.
    groups = defaultdict(list)
    order = sorted(cases.values(), key=lambda c: (c["det"], c["id"]))
    for sysname in systems:
        for task in tasks:
            variants = ["default"] + (VARIANTS[sysname].get(task, []) if a.variants == "all" else [])
            for v in variants:
                for c in order:
                    groups[(sysname, task, v, c["rank"], c["kind"])].append(c["id"])
    fails = Counter()
    for (sysname, cid, task, v), r in done.items():
        if cid in cases and r["status"] in HARD_FAIL:
            fails[(sysname, task, v, cases[cid]["rank"], cases[cid]["kind"])] += 1
    queues = {s: [] for s in systems}
    for g, ids in groups.items():
        for pos, cid in enumerate(ids):
            if (g[0], cid, g[1], g[2]) not in done:
                queues[g[0]].append((pos, g, cid))
    for s in systems:
        queues[s].sort(key=lambda x: (x[0], x[1]))
        queues[s].reverse()                       # pop() from the end
    total = {s: len(queues[s]) for s in systems}
    print(f"{len(cases)} cases; jobs to run: " + ", ".join(f"{s} {total[s]}" for s in systems)
          + (f" ({len(done)} already in {respath})" if done else ""), flush=True)
    if a.dry_run or not any(total.values()):
        return

    lock = threading.Lock()
    new_file = not respath.exists() or a.redo
    # Append in the existing file's columns (older files have no job_s).
    columns = RESULT_COLUMNS if new_file else respath.read_text().split("\n", 1)[0].split("\t")
    R = open(respath, "w" if new_file else "a")
    if new_file:
        R.write("\t".join(RESULT_COLUMNS) + "\n")
    finished = Counter()
    status_counts = defaultdict(Counter)
    meta = {"started": time.strftime("%F %T"), "host": socket.gethostname(), "nproc": os.cpu_count(),
            "loadavg_start": os.getloadavg(), "args": vars(a), "servers": {}, "starts": {}}
    errors = []

    def record(sysname, cid, task, v, res):
        row = {"system": sysname, "case": cid, "task": task, "variant": v, "K": "0",
               "wall_us": "-", "cpu_us": "-", "extra": "-", "value": "-",
               "finished": time.strftime("%F %T"), "job_s": "-"}
        row.update(res)
        with lock:
            R.write("\t".join(str(row.get(c, "-")) for c in columns) + "\n")
            R.flush()
            finished[sysname] += 1
            status_counts[sysname][row["status"]] += 1

    stop = threading.Event()      # set on interrupt: take no more jobs, start no servers

    def next_job(sysname):
        with lock:
            if stop.is_set():
                return None
            q = queues[sysname]
            while q:
                pos, g, cid = q.pop()
                if a.max_timeouts and fails[g] >= a.max_timeouts:
                    skipped.append((sysname, cid, g[1], g[2]))
                    continue
                return g, cid
            return None

    skipped = []

    servers = []

    def worker(sysname, slot):
        srv = Server(sysname, slot, a, out)
        with lock:
            servers.append(srv)
        try:
            while True:
                job = next_job(sysname)
                if job is None:
                    break
                g, cid = job
                _, task, v, _, _ = g
                for attempt in range(2):      # a server that dies gets one restart per job
                    try:
                        if srv.proc is None:
                            with lock:
                                if stop.is_set():
                                    return
                                srv.launch()
                            srv.wait_ready()
                            with lock:
                                meta["servers"].setdefault(sysname, srv.ready)
                                meta["starts"][sysname] = meta["starts"].get(sysname, 0) + 1
                        t0 = time.time()
                        res = srv.request(cid, task, v)
                        res["job_s"] = f"{time.time() - t0:.3f}"
                        break
                    except ServerDied as e:
                        srv.kill()
                        if stop.is_set():
                            return                # killed by us: leave the job for a rerun
                        res = {"status": "crash", "value": str(e)[:300].replace("\t", " ")}
                        if "did not start" in str(e):
                            with lock:
                                errors.append(str(e))
                            raise
                if stop.is_set() and res["status"] in HARD_FAIL:
                    return
                if res["status"] in HARD_FAIL:
                    with lock:
                        fails[g] += 1
                record(sysname, cid, task, v, res)
        except ServerDied:
            pass
        finally:
            srv.kill()

    def progress():
        while any(t.is_alive() for t in threads):
            time.sleep(a.progress)
            with lock:
                parts = []
                for s in systems:
                    c = status_counts[s]
                    parts.append(f"{s} {finished[s]}/{total[s]} (" +
                                 " ".join(f"{k} {c[k]}" for k in sorted(c)) + ")")
            print(time.strftime("%T ") + "; ".join(parts) + f"; load {os.getloadavg()[0]:.0f}", flush=True)

    def interrupt(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupt)
    signal.signal(signal.SIGHUP, interrupt)
    threads = [threading.Thread(target=worker, args=(s, k), daemon=True)
               for s in systems for k in range(a.jobs)]
    for t in threads:
        t.start()
    threading.Thread(target=progress, daemon=True).start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        print("interrupted; finished jobs are saved, rerun to resume", flush=True)
    finally:
        # The servers run in their own sessions, so nothing else stops them.  With
        # stop set (under the lock) no worker launches another one.
        with lock:
            stop.set()
        for srv in list(servers):
            srv.kill()
        for t in threads:
            t.join(timeout=10)
    for (sysname, cid, task, v) in skipped:
        record(sysname, cid, task, v, {"status": "skipped"})
    R.close()
    meta["finished"] = time.strftime("%F %T")
    meta["loadavg_end"] = os.getloadavg()
    with open(out / "meta.jsonl", "a") as M:
        M.write(json.dumps(meta, default=str) + "\n")
    for e in errors:
        print("ERROR: " + e, file=sys.stderr)
    print("done: " + "; ".join(f"{s} " + " ".join(f"{k} {c}" for k, c in sorted(status_counts[s].items()))
                               for s in systems), flush=True)
    if errors:
        sys.exit(1)


# ---------------------------------------------------------------------------
# report

def per_call(r):
    """Seconds per call: the fastest block (the least disturbed by load)."""
    if r is None or r["status"] != "ok" or r["wall_us"] == "-":
        return None
    return min(int(x) for x in r["wall_us"].split(",")) / 1e6


def rough(r):
    """Whether a measurement's timed blocks stayed under 50 ms (K at its cap), so
    that Magma's 10 ms clock makes it imprecise."""
    t = per_call(r)
    return t is not None and int(r["K"]) * t < 0.05


def fmt_t(t):
    if t is None:
        return "-"
    if t < 1e-3:
        return f"{t * 1e6:.0f} µs"
    if t < 1:
        return f"{t * 1e3:.1f} ms"
    return f"{t:.2f} s"


def gmean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def short(v, n=60):
    return v if len(v) <= n else v[:n] + "…"


def cmd_report(a):
    out = Path(a.results)
    cases = read_cases(a.cases) if a.cases else read_cases(json.loads(
        (out / "meta.jsonl").read_text().splitlines()[-1])["args"]["cases"])
    res = load_results(out / "results.tsv")
    metas = [json.loads(l) for l in (out / "meta.jsonl").read_text().splitlines()] \
        if (out / "meta.jsonl").exists() else []
    L = []
    L.append("# Magma vs Oscar benchmark\n")
    for m in metas:
        servers = "; ".join(f"{s}: {d['version']} (load {d['load_s']} s, warm-up {d['warmup_s']} s)"
                            for s, d in m.get("servers", {}).items())
        L.append(f"- run {m['started']} – {m.get('finished', '?')} on {m['host']} ({m['nproc']} CPUs), "
                 f"-j {m['args']['jobs']} per system, timeout {m['args']['timeout']} s, load average "
                 f"{m['loadavg_start'][0]:.0f} → {m['loadavg_end'][0]:.0f}. {servers}"
                 if 'loadavg_end' in m else f"(unfinished). {servers}")
    L.append("")
    L.append("Times are per call (the fastest timed block; see bench_magma.m), ratio = Oscar time / "
             "Magma time over the cases where both succeeded (< 1: Oscar faster).  "
             "Default variants: Magma = the pipeline's code path, Oscar = its documented function.  "
             "`fails` counts timeout/oom/crash/error/skipped; `n/a` = unsupported.\n")

    mismatches = []
    summary = []
    tasks = [t for t in TASKS if any(k[2] == t for k in res)]
    for task in tasks:
        rows = defaultdict(list)
        for cid, c in cases.items():
            rm = res.get(("magma", cid, task, "default"))
            ro = res.get(("oscar", cid, task, "default"))
            if rm is None and ro is None:
                continue
            rows[(c["rank"], c["kind"])].append((cid, c, rm, ro))
        unsupported = lambda r: r is None or r["status"] == "unsupported"
        rows = {key: items for key, items in rows.items()      # drop what neither system supports
                if not all(unsupported(rm) and unsupported(ro) for _, _, rm, ro in items)}
        L.append(f"## {task}\n")
        if not rows:
            L.append("No cases" + (": it needs a second class of the genus, from fill output "
                                   "(lattice_basic_data) in the sampled tree.\n" if task == "isometric_other" else ".\n"))
            continue
        L.append("| rank | kind | cases | both ok | Magma only | Oscar only | Magma fails | Oscar fails "
                 "| n/a (M/O) | median Magma | median Oscar | ratio (geo. mean) | ratio range | mismatches |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        allr = []
        nrough = 0
        for (rank, kind) in sorted(rows):
            items = rows[(rank, kind)]
            both, monly, oonly = [], 0, 0
            mf = of = mna = ona = mism = 0
            for cid, c, rm, ro in items:
                sm = rm["status"] if rm else "missing"
                so = ro["status"] if ro else "missing"
                mna += sm == "unsupported"
                ona += so == "unsupported"
                mf += sm not in ("ok", "unsupported", "missing")
                of += so not in ("ok", "unsupported", "missing")
                if sm == "ok" and so == "ok":
                    both.append((per_call(rm), per_call(ro)))
                    nrough += rough(rm) or rough(ro)
                    if rm["value"] != ro["value"]:
                        mism += 1
                        mismatches.append((task, cid, rm["value"], ro["value"]))
                elif sm == "ok":
                    monly += 1
                elif so == "ok":
                    oonly += 1
            ratios = [o / m for m, o in both if m and o]
            allr += ratios
            g = gmean(ratios)
            rr = f"{min(ratios):.2g}–{max(ratios):.2g}" if ratios else "-"
            L.append(f"| {rank} | {kind} | {len(items)} | {len(both)} | {monly} | {oonly} | {mf} | {of} "
                     f"| {mna}/{ona} | {fmt_t(statistics.median([m for m, _ in both]) if both else None)} "
                     f"| {fmt_t(statistics.median([o for _, o in both]) if both else None)} "
                     f"| {f'{g:.3g}' if g else '-'} | {rr} | {mism} |")
            summary.append((task, rank, kind, len(items), len(both), monly, oonly, mf, of, g, mism))
        g = gmean(allr)
        L.append(f"\nAll ranks: ratio geometric mean {f'{g:.3g}' if g else '-'} over {len(allr)} cases."
                 + (f"  {nrough} of them have a timed block under 50 ms (calls too fast for the block "
                    f"cap), so their times are rough." if nrough else "") + "\n")

    # Variants against the same system's default.
    vlines = []
    for sysname, per_task in VARIANTS.items():
        for task, vs in per_task.items():
            for v in vs:
                ratios, mism, n, ok = [], 0, 0, 0
                for cid in cases:
                    rv = res.get((sysname, cid, task, v))
                    rd = res.get((sysname, cid, task, "default"))
                    if rv is None:
                        continue
                    n += 1
                    ok += rv["status"] == "ok"
                    if rv["status"] == "ok" and rd and rd["status"] == "ok":
                        tv, td = per_call(rv), per_call(rd)
                        if tv and td:
                            ratios.append(tv / td)
                        if rv["value"] != rd["value"]:
                            mism += 1
                            mismatches.append((f"{task} ({sysname} {v} vs default)", cid, rd["value"], rv["value"]))
                if n:
                    g = gmean(ratios)
                    vlines.append(f"| {sysname} | {task} | {v} | {n} | {ok} | {f'{g:.3g}' if g else '-'} | {mism} |")
    if vlines:
        L.append("## Variants\n")
        L.append("Each variant against its own system's default, over the cases where both succeeded.\n")
        L.append("| system | task | variant | cases | ok | time ratio variant/default | value mismatches |")
        L.append("|---|---|---|---|---|---|---|")
        L += vlines
        L.append("")

    # Values against what the census stores.
    exp_lines = []
    for task, col in (("class_number", "class_number"), ("mass", "mass")):
        for sysname in SYSTEMS:
            n = bad = 0
            for cid, c in cases.items():
                r = res.get((sysname, cid, task, "default"))
                if r and r["status"] == "ok" and c[col] != "-":
                    n += 1
                    stored = c[col][:-2] if c[col].endswith("/1") else c[col]   # older samples wrote n/1
                    if r["value"] != stored:
                        bad += 1
                        mismatches.append((f"{task} ({sysname} vs census)", cid, c[col], r["value"]))
            if n:
                exp_lines.append(f"| {task} | {sysname} | {n} | {bad} |")
    if exp_lines:
        L.append("## Against the census\n")
        L.append("| task | system | cases with a stored value | disagreements |")
        L.append("|---|---|---|---|")
        L += exp_lines
        L.append("")

    if mismatches:
        L.append(f"## Value mismatches ({len(mismatches)}; all in mismatches.tsv)\n")
        L.append("| task | case | Magma / reference | Oscar / variant |")
        L.append("|---|---|---|---|")
        for t, cid, x, y in mismatches[:a.max_listed]:
            L.append(f"| {t} | {cid} | `{short(x)}` | `{short(y)}` |")
        L.append("")
    with open(out / "mismatches.tsv", "w") as F:
        F.write("task\tcase\tfirst\tsecond\n")
        for t, cid, x, y in mismatches:
            F.write(f"{t}\t{cid}\t{x}\t{y}\n")

    # Where the server time went (job_s: each job as bench.py saw it, including the
    # timed blocks' untimed preparation and the timeouts).
    spent = defaultdict(lambda: [0, 0.0, 0])
    for (sysname, cid, task, v), r in res.items():
        js = r.get("job_s", "-")
        if js not in ("-", None, ""):
            e = spent[(sysname, task, v)]
            e[0] += 1
            e[1] += float(js)
            e[2] += r["status"] in HARD_FAIL
    if spent:
        starts = Counter()
        for m in metas:
            starts.update(m.get("starts", {}))
        L.append("## Where the server time went\n")
        L.append("From `job_s` (whole jobs, timeouts included).  Server starts: "
                 + ", ".join(f"{s} {k}" for s, k in sorted(starts.items()))
                 + " (each one also costs the start-up and warm-up listed at the top).\n")
        L.append("| system | task | variant | jobs | server hours | of which timeouts (count) |")
        L.append("|---|---|---|---|---|---|")
        for (sysname, task, v), (n, t, nto) in sorted(spent.items(), key=lambda x: -x[1][1]):
            L.append(f"| {sysname} | {task} | {v} | {n} | {t / 3600:.2f} | {nto} |")
        tot = defaultdict(float)
        for (sysname, _, _), (_, t, _) in spent.items():
            tot[sysname] += t
        L.append("\nTotal: " + ", ".join(f"{s} {t / 3600:.1f} server hours" for s, t in sorted(tot.items())) + ".\n")

    errs = Counter()
    example = {}
    for (sysname, cid, task, v), r in res.items():
        if r["status"] in ("error", "crash"):
            key = (sysname, task, v, r["status"], r["value"][:120])
            errs[key] += 1
            example.setdefault(key, cid)
    if errs:
        L.append("## Errors and crashes\n")
        L.append("| system | task | variant | status | count | message | example case |")
        L.append("|---|---|---|---|---|---|---|")
        for key, k in errs.most_common(a.max_listed):
            sysname, task, v, st, msg = key
            L.append(f"| {sysname} | {task} | {v} | {st} | {k} | `{msg}` | {example[key]} |")
        L.append("")
    with open(out / "summary.tsv", "w") as F:
        F.write("task\trank\tkind\tcases\tboth_ok\tmagma_only\toscar_only\tmagma_fails\toscar_fails"
                "\tratio_gmean\tmismatches\n")
        for row in summary:
            F.write("\t".join("-" if x is None else (f"{x:.4g}" if isinstance(x, float) else str(x))
                              for x in row) + "\n")
    (out / "report.md").write_text("\n".join(L) + "\n")
    print(f"wrote {out / 'report.md'}, {out / 'summary.tsv'} and {out / 'mismatches.tsv'} "
          f"({len(mismatches)} mismatches)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sample", help="draw a stratified sample of genera into a cases file")
    s.add_argument("--data-root", default=".", help="census tree (genera_basic/ etc.), default .")
    s.add_argument("--ranks", default="2-16", help="rank range, e.g. 2-16 (default)")
    s.add_argument("--kinds", default="def,indef", help="def (positive definite) and/or indef")
    s.add_argument("--per-stratum", type=int, default=3)
    s.add_argument("--min-det", type=int, default=1)
    s.add_argument("--max-det", type=int, default=0, help="0 = no bound")
    s.add_argument("--theta-vectors", type=float, default=1e5,
                   help="aim the theta precision at about this many vectors (default 1e5)")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("-o", "--output", default="bench_cases.tsv")

    r = sub.add_parser("run", help="run the tasks in Magma and Oscar")
    r.add_argument("cases")
    r.add_argument("-o", "--output", default="bench_results")
    r.add_argument("--systems", default="magma,oscar")
    r.add_argument("--tasks", default="all", help="comma-separated, default all: " + ",".join(TASKS))
    r.add_argument("--variants", default="default", choices=["default", "all"])
    r.add_argument("-j", "--jobs", type=int, default=4, help="servers per system (default 4)")
    r.add_argument("--timeout", type=int, default=600, help="seconds per job (default 600)")
    r.add_argument("--min-time", type=float, default=0.5,
                   help="repeat timed blocks until this many seconds (default 0.5)")
    r.add_argument("--max-reps", type=int, default=3, help="at most this many timed blocks (default 3)")
    r.add_argument("--max-timeouts", type=int, default=2,
                   help="skip the rest of a (system, task, variant, rank, kind) group after this many "
                        "timeouts; 0 = never (default 2)")
    r.add_argument("--memlimit", type=int, default=16, help="GB per server, 0 = none (default 16)")
    r.add_argument("--startup-timeout", type=int, default=1800)
    r.add_argument("--magma", default="magma")
    r.add_argument("--julia", default="julia")
    r.add_argument("--progress", type=int, default=60, help="seconds between progress lines")
    r.add_argument("--redo", action="store_true", help="start results.tsv afresh")
    r.add_argument("--dry-run", action="store_true", help="only count the jobs")

    t = sub.add_parser("report", help="summarise a results directory")
    t.add_argument("results")
    t.add_argument("--cases", help="cases file (default: the one the last run used)")
    t.add_argument("--max-listed", type=int, default=50)

    a = p.parse_args()
    {"sample": cmd_sample, "run": cmd_run, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
