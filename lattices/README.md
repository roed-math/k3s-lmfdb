# Lattice data pipeline

This folder generates the data behind the LMFDB's integral-lattice section: the
`lat_genera` and `lat_lattices` tables (see [`../schemas/`](../schemas) for the
column-by-column schemas and the label format), and — via the Oscar stage under
[`oscar/`](oscar) — the `lat_inclusions` table of embeddings between lattices.

## Pipeline overview

The driver is [`run_all.py`](run_all.py), which runs under `sage -python` and
orchestrates the stages below.  All stages write flat files under this
directory, sharded as `<folder>/<rank>/<nplus>/<identifier>` (see `LabelPath` in
[`connect_genus.m`](connect_genus.m)); each per-stage `.format` file lists the
pipe-separated columns of the corresponding records.

1. **List genera** (Sage, [`genus.py`](genus.py)): enumerate genus symbols for
   every signature `(nplus, nminus)` with `nminus <= rank/2` and `|det|` up to
   `--max-disc-ratio / rank` (plus the K3 extension to `D/(22-rank)` for
   `12 <= rank <= 20`, `nminus in {1,2}`, disabled by `--nokthree`).  Writes
   `genera_basic/`, including the LMFDB label and one representative Gram
   matrix per genus.
2. **Enumerate genera with Oscar** ([`oscar/pool.jl`](oscar/pool.jl) `--task
   enumerate`, on a pool of persistent Oscar workers): indefinite genera of rank
   >= 3 (Magma's spinor kernel overcounts some of them) and definite genera of rank
   <= 8 (Oscar is faster there; Magma from rank 9).  Writes one record per genus to
   `oscar_genus/`: the classes (indefinite ones in label order), with automorphism
   group orders and root systems for definite genera.  `--no-oscar-enum` skips it;
   `--oscar-jobs/-timeout/-memlimit/-batch`, `--oscar-definite-max-rank`,
   `--oscar-indefinite-min-rank` control it.
3. **Fill genera** (Magma, [`fill_genus.m`](fill_genus.m) via
   [`run_fill_genus.m`](run_fill_genus.m) under GNU parallel): compute all
   classes in each genus — from the Oscar record when there is an `ok` one,
   otherwise by default with the orbit-method enumeration of
   [`genus_reps.m`](genus_reps.m) (`--no-orth` falls back to the p-neighbour
   walk) — then per-lattice invariants (automorphism group, theta series,
   hashes, `is_universal`/`is_even_universal`, ...).  Writes
   `lattice_basic_data/`, `lattice_hashes/`, `genera_advanced/`.  Guarded by
   `--enum-masslimit/-timeout/-timelimit/-sizelimit/-adjlimit/-per-lattice-timeout`
   so pathological genera are bounded rather than hanging.  The classes of a
   definite genus are sorted by automorphism group order (larger first), minimum
   (larger first), root system (`RootString`/`CompareRootStrings` in
   [`orth_dec.mag`](orth_dec.mag), larger first), theta series, dual theta
   series and Gram matrix (`cmp_lat`), and that order is the label tiebreaker.
   A key whose computation timed out sorts after every computed one, so the
   labels are fixed by the run that enumerates the genus.
   **Festi-Veniani index** (Julia/Oscar, [`oscar/fix_festi_veniani.jl`](oscar/fix_festi_veniani.jl),
   right after fill): fill records the order of the image of O(L) on the
   discriminant group in `festi_veniani/` and leaves `festi_veniani_index` null;
   this computes |O(q_L)| per genus (Magma has no orthogonal groups of torsion
   quadratic forms) and writes [O(q_L) : image].  Chunks of `--fv-chunk` genera run
   under GNU parallel (`--fv-jobs`, `--fv-timelimit`); `--skip-festi-veniani` skips
   it.  Data filled before 2026-11 stored #Aut(L#/L) / |image| instead (the
   automorphism group of the abelian group, a different number whenever O(q_L) is
   smaller); convert it once with
   `julia --project=oscar oscar/fix_festi_veniani.jl --legacy .`, which first saves
   the image orders to `festi_veniani/`, so a rerun is harmless.  Keep
   `festi_veniani/` with the data.
   **Theta-Hecke pass** (Magma, [`theta_hecke.m`](theta_hecke.m) via
   [`run_theta_hecke.m`](run_theta_hecke.m), right after fill):
   `is_theta_distinguished` says whether the theta series tell the classes of a
   definite genus apart.  Fill decides it when the theta series differ early or
   agree to the Sturm precision (`ThetaSturmPrecision`, extending them when that
   is cheap); the rest it leaves null and records in `theta_hecke/`.  The pass
   decides those from the coefficients fill already has, using Eichler's
   commutation relation: summing theta series over the p-neighbours of a class is
   the Hecke operator T(p^2) + p + ... + p^(n-3) applied to its theta series.  So
   the relations among the theta series are stable under the adjacency matrix;
   a pair whose difference falls outside the largest stable space of relations
   that hold to the known precision has different theta series (no modular forms
   needed), and the remaining relations lie in the Hecke eigenspaces of the cusp
   space matching the adjacency eigenvalues, which a handful of coefficients
   usually determines -- then they hold exactly.  Genera sharing a cusp space run
   in one job (`--theta-hecke-chunk`, `--theta-hecke-jobs`), which builds the
   space (modular symbols; a q-expansion basis in half-integral weight) once;
   `--theta-hecke-timeout` bounds each adjacency matrix and space.  Verdicts go to
   `theta_hecke_done/` and into `genera_advanced`; eigenspace data per cusp space
   is kept in `theta_hecke_spaces/` (it depends only on the space, so `--clean`
   leaves it).  `--skip-theta-hecke` skips the pass.
4. **Tensor products** (Magma, [`run_tensor.m`](run_tensor.m)).
5. **Names** (Magma, [`run_basic_names.m`](run_basic_names.m)): atomic lattice
   names, merged across workers.
6. **Connect** (Magma, [`connect_genus.m`](connect_genus.m) via
   [`run_connect_genus.m`](run_connect_genus.m), two passes): everything that
   needs data across genera — duals, orthogonal decompositions, geometric
   invariants, names, and the universality/regularity columns (see below).
   Writes `lattice_advanced_data/`.  Label lookups `FindLabel` cannot decide --
   in an indefinite genus of rank >= 3 with several classes (an isometry test
   Magma lacks), or where hashing or Magma's isometry test overran -- are
   recorded in `label_lookups/` with their candidate classes, and the Oscar
   pool's `resolve` task ([`oscar/LMFDBLattices/src/resolve.jl`](oscar/LMFDBLattices/src/resolve.jl))
   fills them in with `is_isometric` between the two connect passes.  This covers
   `dual_label`, `primitive_scaling`, `even_sublattice`, `norm1_complement`,
   `root_complement` and the orthogonal factors (once every summand is labelled,
   the resolver writes the factor list and multiplicities, and pass 2 derives the
   decomposable fields and name from them); lookups from tensor/names stay null.
7. **Embeddings** (Julia/Oscar, the Oscar worker pool's `embed` task,
   [`oscar/LMFDBLattices/src/embed.jl`](oscar/LMFDBLattices/src/embed.jl);
   [`oscar/driver.jl`](oscar/driver.jl) runs one job by hand; design in
   [`oscar/INTEGRATION_PLAN.md`](oscar/INTEGRATION_PLAN.md)):
   prime-index overlattices, corank-1 primitive embeddings (rank n → n+1) and
   primitive-embedding counts into the K3 lattice, from the fill outputs only —
   no database connection.  Writes raw data under `inclusions/` and assembles
   `lat_inclusions.tbl` via [`oscar/write_table.jl`](oscar/write_table.jl).
   Runs last (it is the slowest optional stage); controlled by the `--embed-*`
   flags (task list, det cap, class-number cap, worker count, per-job timeout
   and memory cap —
   the corank-1 task is quadratic in the number of genera per rank, so start
   with a small `--embed-max-disc-ratio`).  Each task's domain genera are split
   into det-sorted chunks of `--embed-chunk`, one job each; a finished job is
   committed atomically, so a timeout loses one chunk, and rerunning without
   `--clean` (e.g. `--skip-list-genera --skip-enumerate-genera`) computes only
   the jobs that did not finish.  Each committed job also holds `units.jsonl`,
   one line per domain genus (task, classes, seconds, rows, status), merged into
   `embed_units.jsonl` after the stage for calibration.  Volume guard: a domain
   genus producing more than `--embed-max-rows-per-domain` rows (default 10^4)
   has its rows dropped and is reported as overflowed (by `run_all.py` and
   `write_table.jl`); its job's manifest lists it with the unloadable genera, so
   the universality certificates stop below its det.  Skipped with a warning if
   the Julia/Oscar environment is not set up.
8. **Refine** (Python, [`oscar/refine_universality.py`](oscar/refine_universality.py),
   inline in `run_all.py`): a bottom-up bitset pass over the embedding poset —
   each rank-n lattice gets a distinct power of two, masks are OR-ed up the
   corank-1 and overlattice edges, and then all-ones tests db-relative
   n-universality while mask-equals-genus-OR tests n-regularity.  Raise-only
   into the heuristic (`*_proved = false`) columns of `lattice_advanced_data`;
   mask deficits against the genus OR are written to `regularity_candidates.txt`
   for a Magma `RepresentsLattice` confirmation pass.

A small end-to-end run:

```sh
cd lattices
PATH="/path/to/magma/dir:$PATH" sage -python run_all.py -R 4 -D 120 -j 8 \
    --skip-embeddings --skip-tensor
```

Per-stage GNU parallel job logs (`fill.joblog`, `connect1.joblog`, ...) record
failures; per-label errors are printed as `ERROR: <label>: ...` on stdout, so a
retry list is `grep '^ERROR:' output | cut -d: -f2`.  `--clean` wipes the output
directories first so a rerun cannot silently blend with stale data.

## Installation

The pipeline needs Magma, SageMath, GNU parallel, and (for the embedding stage)
Julia with a development version of Oscar.  PARI/GP is used only through Sage's
bundled PARI (e.g. Denis Simon's `qfsolve.gp`, shipped with Sage, loaded by
[`lasVegas.py`](lasVegas.py)), so no standalone `gp` installation is required.

### Git submodules

After cloning, fetch the Magma dependencies
([FiniteGroups](FiniteGroups), [Parallel.magma](Parallel.magma),
[lattools](lattools)):

```sh
git submodule update --init
```

`Parallel.magma` points at `roed-math/Parallel.magma` (branch `fast-hardkill`):
the pinned `TimeoutCall` HardKill watchdog is not in `edgarcosta/Parallel.magma`.
In a clone made before that change, run `git submodule sync` first.

### Magma

Version 2.29 or later is required (`create_label.m` uses intrinsics introduced
in 2.29; minor releases matter — e.g. V2.29-4's `GenusRepresentatives` had a
definite-lattice bug that the orbit-method path avoids).  The `magma` binary
must be on `PATH` for the GNU parallel stages; on macOS the default install
location is not on `PATH`, so run the pipeline with
`PATH="/Applications/Magma:$PATH"`.

Magma must also carry the three required patches in [`magma_patches/`](magma_patches):
without `neighbors-2.30-1.patch`, `GenusRepresentatives` returns wrong classes for
indefinite binary lattices (lattices from another genus, missing classes); without
`mass-2.30-1.patch`, `Mass` is wrong or fails for many lattices of rank 2 mod 4 from
rank 14 on, which breaks the mass-verified enumerations; without
`diagonal-2.30-1.patch`, the spinor characters of indefinite genera of rank >= 3 are
wrong, so their class numbers and class representatives are.  Apply each from Magma's
root directory with `patch -p0 < .../magma_patches/<name>.patch` (try `--dry-run`
first on other versions); `magma -b check_magma.m` (from `lattices/`) prints
`MAGMA_OK` once all three are in, and `run_all.py` runs that check before the fill stage.

### SageMath

Sage 10.x provides the Python that `run_all.py` runs under (`sage -python`),
the genus-symbol enumeration, and the bundled PARI.  Install from
[sagemath.org](https://www.sagemath.org) or your package manager and make sure
`sage` is on `PATH`.

### GNU parallel

`brew install parallel` (macOS) or your package manager.

### Julia + Oscar (enumeration and embedding stages)

Until Oscar 1.9 is released, the embedding code needs the development branch
(the primitive-embedding refactor of
[oscar-system/Oscar.jl#5949](https://github.com/oscar-system/Oscar.jl/pull/5949)).
The project file in [`oscar/`](oscar) pins a known-good state:

```sh
brew install juliaup            # or curl -fsSL https://install.julialang.org | sh
juliaup add release             # Julia >= 1.10; 1.12.x tested
cd lattices/oscar
julia --project=. -e 'using Pkg; Pkg.instantiate()'
# smoke test (first run precompiles, 10-30 minutes):
julia --project=. -e 'using Oscar;
    @assert isdefined(Oscar, :_primitive_embeddings_in_unimodular_safe);
    println("Oscar OK ", Oscar.VERSION_NUMBER)'
```

Workers should run with `JULIA_NUM_THREADS=1`; budget a few GB of RSS each.

The pipeline's own Oscar code is the package [`oscar/LMFDBLattices`](oscar/LMFDBLattices)
(already in the project).  Its precompile workload runs every task once on small
lattices, so the compiled code is cached: the first definite `representatives`
call in a fresh process went from 110 s of compilation to 0.9 s (mimir, load ~160).
That precompilation happens the first time the package is loaded on a machine and
after any change to it (~5 min at low load); `run_all.py`'s Oscar preflight
triggers it.  The worker pool (`LMFDBLattices/src/pool.jl`, Distributed.jl) pays the
remaining ~30 s load once per worker, kills and replaces a worker whose genus
overruns `--timeout` or `--memlimit`, and resumes from existing records.  Test it
with `julia --project=oscar oscar/test_pool.jl`.

## Tests

```sh
cd lattices
tests/run.sh                    # the fast tests/test_*.m suites
tests/run.sh --slow             # also test_brandhorst and test_l12a_plus (~2.5 h)
tests/run.sh test_universality  # one suite, streaming its output
python3 tests/test_schema_consistency.py
python3 oscar/test_refine_universality.py
julia --project=oscar oscar/test_fix_festi_veniani.jl
julia --project=oscar oscar/test_hecke_fixes.jl
```

CI (`.github/workflows/tests.yml`) runs one job per Magma suite once the
`MAGMA_URL`/`MAGMA_USERNAME`/`MAGMA_PASSWORD` repository secrets are set.

A suite that calls anything forking via `TimeoutCall` must wrap its body in
`procedure() ... end procedure();`, as the `run_*.m` drivers do: a forked
child's exit moves the parent's read position in the script file, and the
parent then misparses the rest of it.  Note also that `TimeoutCall`'s
child-side `Alarm` cannot interrupt a single long C-level call; pass
`HardKill := true` (a watchdog that SIGKILLs the child, ~10 ms per call) for
anything that may enter one — see `Parallel.magma/TimeoutCall.m`.  Fill and
connect run every per-lattice computation that can run unbounded this way;
nothing can interrupt the batch process itself.  A child discards whatever Magma
initialises lazily on first use, so `WarmLatticeCaches` (in
[`fill_genus.m`](fill_genus.m)) initialises it once in the batch process.

## Rerunning and provenance

- `run_all.py --retry-missing` refills only the genera listed in
  `fill_missing_<rank>.txt` (errors and crashes left by an earlier run), one per
  job, with the guards given on its command line, and retries failed Oscar
  enumeration records.
- Every run appends to `run_provenance.jsonl` in `lattices/`: a `start` record
  (commit, whether the tree was dirty, submodule pins, host, arguments, Sage),
  `magma` and `oscar` records with the versions the checks found, and an `end`
  record with each stage's wall-clock time and the exit status.  Keep it with the
  data.

## Calibrating the fill stage

The fill guards (`--enum-*`) and the batch weights (`genus_work` in
`run_all.py`) were tuned at C=256.  [`calibrate_fill.py`](calibrate_fill.py)
measures them against a listed census.  Use a separate checkout (every stage
writes under `lattices/`), on the machine the census will run on, with load
comparable to production:

```sh
cd lattices
# 1. stratified sample: N genera per (rank, nminus, floor(log2 det)), copied into
#    ./genera_basic, plus the census population of every stratum
python3 calibrate_fill.py sample --census /path/to/census --per-stratum 4 --ranks 1-32
# 2. fill the sample with per-genus stats (here: the production guards)
sage -python run_all.py --skip-list-genera --clean --fill-stats -j 32 \
    --skip-embeddings --skip-tensor --skip-names --skip-connect
# 3. per-stratum outcomes, time quantiles, class numbers, slowest lattice, and
#    census CPU-hours extrapolated from the stratum populations
python3 calibrate_fill.py report --costs-out genus_costs.tsv > calibration_report.txt
# 4. later runs batch by the measured costs (expensive genera alone, longest first)
sage -python run_all.py ... --genus-costs genus_costs.tsv
```

To price connect as well, run it on the filled sample before the report:
`sage -python run_all.py --skip-list-genera --skip-enumerate-genera --connect-stats
-j 32 --skip-embeddings --skip-tensor --skip-names`.  The report then adds a connect
table.  Lookups into genera outside the sample find nothing and return quickly, so
connect times measured this way are a lower bound.

`--clean` here leaves the sampled `genera_basic/` alone (it is the input when
`--skip-list-genera` is given).  Each `fill_stats/` record gives the outcome
(`ok`, or the guard that decided it: `mass_skip`, `enum_fail`, `size_skip`,
`time_discard`), so one pilot with the production guards shows what coverage
they give; a second pass with generous guards on the strata where they trip
measures what those strata actually cost (a tripped guard is a censored
observation: its time is a lower bound).  Sampled genera start with a cold orth
cache (their rank+1 parents are not filled first, as they are in production),
so enumeration times are upper bounds.

## Benchmarking Magma against Oscar

[`benchmark/`](benchmark) times Magma and Oscar on the quantities the pipeline
stores that both can compute, and checks that they agree:

| task | Magma (default = pipeline code path) | Oscar |
| --- | --- | --- |
| `class_number` | `GenusReps` (stock: `GenusRepresentatives`) | `genus_representatives` (also `enumerate_definite_genus` `:orbit`/`:random`) |
| `mass` | `Mass` | `mass` |
| `aut_order` | `AutomorphismGroupFaster` on the `LLLGram`, as in fill (stock: `AutomorphismGroup`) | `isometry_group` (also `automorphism_group_order`) |
| `aut_label` | `IdentifyGroup` | `small_group_identification` |
| `minimum`, `kissing`, `successive_minima` | `Minimum`, `KissingNumber`, `SuccessiveMinima` | `minimum`, `kissing_number`, `successive_minima` |
| `root_system` | `RootString` | `root_lattice_recognition` |
| `theta` | `ThetaSeries` | counted from `short_vectors_iterator` (Oscar has no theta series) |
| `orth_decomp` | `OrthogonalDecompositionFaster` (stock: `OrthogonalDecomposition`) | `irreducible_components` |
| `isometric_same`, `isometric_other` | `IsIsometric`, `IsIsometricBinary` | `is_isometric` |
| `disc_image` | image of O(L) in Aut(L#/L), as fill's `aut_data` | `image_in_Oq` |

Only Magma has canonical forms, Hecke (adjacency) matrices, the hashes,
universality, Voronoi cells and t-designs, so those are not compared.  Most tasks
need a definite lattice in both systems; for indefinite lattices the comparison is
`class_number`, the isometry tests (Magma: rank 2 only) and `disc_image` (Oscar
only).  [`oscar/invariants.jl`](oscar/invariants.jl) has the Oscar computations
in the pipeline's value formats (root system strings, theta series), for reuse.

```sh
cd lattices
# 1. a stratified sample: N genera per (rank, definite/indefinite, floor(log2 det))
python3 benchmark/bench.py sample --data-root /path/to/census --ranks 2-16 \
    --per-stratum 3 -o bench_cases.tsv
# 2. every task on every case, -j servers per system, Magma and Oscar side by side
python3 benchmark/bench.py run bench_cases.tsv -o bench_results -j 8 --timeout 600
# 3. per task and rank: success counts, median times, Oscar/Magma time ratio, and
#    every value on which the two disagree
python3 benchmark/bench.py report bench_results     # -> bench_results/report.md
```

Point `--data-root` at a tree with fill output (`lattice_basic_data/`,
`genera_advanced/`) when there is one: the sample then uses the stored class .1 and
.2 (a non-isometric pair for `isometric_other`) and the report checks class numbers
against the census.  `run` appends to `bench_results/results.tsv` as jobs finish,
so rerunning it resumes; `--variants all` adds the alternative implementations,
`--tasks` picks some.  Within each (system, task, rank, definite/indefinite) group
the cases run by increasing determinant, and after `--max-timeouts` (default 2) the
rest of the group is recorded as `skipped`.

How it measures (details in the headers of
[`bench_magma.m`](benchmark/bench_magma.m) and
[`bench_oscar.jl`](benchmark/bench_oscar.jl)): each system runs as a long-lived
single-threaded server that `bench.py` feeds jobs over a pipe, warmed up first,
so start-up and Julia compilation are not timed (`report.md` lists them; Oscar's
warm-up took about 2.5 minutes on mimir at load ~30, and 10 minutes at load ~150).  Fast calls are timed in blocks of up to
10000 (1000 at rank 2), since Magma's clock has a 10 ms resolution, and every
call gets its own Gram matrix U G U^t: while a lattice is
referenced, Magma's `LatticeWithGram` returns that same object for an equal Gram
matrix, so repeated calls would reuse its cached minimum and short vectors.  Magma
bounds each job with a HardKill `TimeoutCall` (`--memlimit` via `SetMemoryLimit`);
Julia cannot interrupt a computation, so `bench.py` kills and restarts the Oscar
server on a timeout or past the memory limit (each restart costs that start-up
again).  The Oscar server runs `GC.gc()` before each timed block, so a block does
not pay for garbage left by earlier jobs.
`AutomorphismGroupFaster` is very slow on unreduced bases (see `CLAUDE.md`), so
the default Magma variant reduces first, as fill does; the stock intrinsics and
Oscar reduce internally.  `report.md` also shows where the server time went, per
system and task, from each job's total time (`job_s` in `results.tsv`).
Measured this way, the pipeline's `class_number` includes `GenusReps`' own
`TimeoutCall` fork (several ms), which dominates for small genera; the `stock`
variant shows Magma without it.  The Magma server runs in an empty directory, so
the orbit method starts with a cold `orth_cache/`, as in `calibrate_fill.py`.

## Universality and regularity columns

[`universality.m`](universality.m) computes `is_universal`,
`is_even_universal` (fill stage) and `universality`, `even_universality`,
`is_regular`, `regularity` with their `*_proved` companions (connect stage).
Proved values are theorems, from the known criterion sets — the 15-theorem, the
doubled 290-theorem, Kim–Kim–Oh's 2-universality set, Oh's `{I8, E8}` for
8-universality, Ko's values for `I_m` — met by explicitly non-embedded
witnesses; unproved values are heuristics from finite witness pools, to be
refined (and sometimes proved) by the closure pass over `lat_inclusions`
([`oscar/refine_universality.py`](oscar/refine_universality.py)).  Sources and
the soundness policy are documented at the top of `universality.m` and in
[`oscar/INTEGRATION_PLAN.md`](oscar/INTEGRATION_PLAN.md) §10.
