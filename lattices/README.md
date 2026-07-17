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
2. **Enumerate genera** (Magma, [`fill_genus.m`](fill_genus.m) via
   [`run_fill_genus.m`](run_fill_genus.m) under GNU parallel): compute all
   classes in each genus — by default with the orbit-method enumeration of
   [`genus_reps.m`](genus_reps.m) (`--no-orth` falls back to the p-neighbour
   walk) — then per-lattice invariants (automorphism group, theta series,
   hashes, `is_universal`/`is_even_universal`, ...).  Writes
   `lattice_basic_data/`, `lattice_hashes/`, `genera_advanced/`.  Guarded by
   `--enum-masslimit/-timeout/-timelimit/-sizelimit/-adjlimit` so pathological
   genera are bounded rather than hanging.
3. **Tensor products** (Magma, [`run_tensor.m`](run_tensor.m)).
4. **Names** (Magma, [`run_basic_names.m`](run_basic_names.m)): atomic lattice
   names, merged across workers.
5. **Connect** (Magma, [`connect_genus.m`](connect_genus.m) via
   [`run_connect_genus.m`](run_connect_genus.m), two passes): everything that
   needs data across genera — duals, orthogonal decompositions, geometric
   invariants, names, and the universality/regularity columns (see below).
   Writes `lattice_advanced_data/`.
6. **Embeddings** (Julia/Oscar, [`oscar/driver.jl`](oscar/driver.jl) under GNU
   parallel; design in [`oscar/INTEGRATION_PLAN.md`](oscar/INTEGRATION_PLAN.md)):
   prime-index overlattices, corank-1 primitive embeddings (rank n → n+1) and
   primitive-embedding counts into the K3 lattice, from the fill outputs only —
   no database connection.  Writes raw data under `inclusions/` and assembles
   `lat_inclusions.tbl` via [`oscar/write_table.jl`](oscar/write_table.jl).
   Runs last (it is the slowest optional stage); controlled by the `--embed-*`
   flags (task list, det cap, class-number cap, worker count, per-job timeout —
   the corank-1 task is quadratic in the number of genera per rank, so start
   with a small `--embed-max-disc-ratio`).  Skipped with a warning if the
   Julia/Oscar environment is not set up.
7. **Refine** (Python, [`oscar/refine_universality.py`](oscar/refine_universality.py),
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

### Magma

Version 2.29 or later is required (`create_label.m` uses intrinsics introduced
in 2.29; minor releases matter — e.g. V2.29-4's `GenusRepresentatives` had a
definite-lattice bug that the orbit-method path avoids).  The `magma` binary
must be on `PATH` for the GNU parallel stages; on macOS the default install
location is not on `PATH`, so run the pipeline with
`PATH="/Applications/Magma:$PATH"`.

### SageMath

Sage 10.x provides the Python that `run_all.py` runs under (`sage -python`),
the genus-symbol enumeration, and the bundled PARI.  Install from
[sagemath.org](https://www.sagemath.org) or your package manager and make sure
`sage` is on `PATH`.

### GNU parallel

`brew install parallel` (macOS) or your package manager.

### Julia + Oscar (embedding stage)

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

## Tests

```sh
cd lattices
magma -b run_tests.m                      # everything under tests/
magma -b target:=Universality run_tests.m # one file (substring match)
```

Some tests fork via `TimeoutCall` and print a harmless post-success
`bad syntax`.  Note that `TimeoutCall`'s child-side `Alarm` cannot interrupt a
single long C-level call; pass `HardKill := true` (a watchdog that SIGKILLs the
child) for anything that may enter one — see `Parallel.magma/TimeoutCall.m`.

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
