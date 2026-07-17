# Plan: integrating Stevell's Oscar embedding code into the lattice pipeline

Status: proposal (2026-07-16), written after reading `code/*.jl`, `run_all.py`, the
fill/connect Magma stages, and `schemas/lat_inclusions.md`.  Working repo:
`~/claude/k3s-lmfdb2` (the `k3s-lmfdb` clone is reserved for the genus-enumeration PR).

## 1. What the code computes

`code/lmfdb_lattice_enum.jl` produces three kinds of rows for `lat_inclusions`, all
keyed by LMFDB *lattice* labels and all attached to the **domain** (smaller) lattice:

| computation | input | output row | notes |
| --- | --- | --- | --- |
| **overlat** (`_prime_index_overlattices_lmfdb`) | all genera `Gs` of one rank `n`, det ≤ D | `M ⊂ L`, same rank, prime index `p`, `primitive=F`, `quotient_invs={p}`, `quotient_gens` = orbit reps of glue vectors | `det L = det M / p²`, so codomain is automatically within the same det bound — no extra enumeration needed |
| **k3** (`_embeddings_in_K3_lattice`) | even genera with `nminus ∈ {1,2}` | `M ↪ Λ_K3(-1)` (`II(19,3)`), `primitive=T`, `multiplicity` = # primitive embeddings (count only) | this is the headline feature for the K3 section (Néron–Severi / transcendental lattices) |
| **corank1** (`_primitive_embeddings_lmfdb`) | genera `G2s` of rank `n` (det ≤ D/n) × genera `G1s` of rank `n+1` (det ≤ D/(n+1)) | `M ↪ L` primitive, `quotient_rank=1`, `orthogonal_complement` = coordinates of the rank-1 complement (definite codomains only) | by far the most expensive: quadratic in # genera per rank pair |

`code/writing_file.jl` walks the raw `data/rank{n}/{label}/` files and emits
`lat_inclusions.tbl` in (approximately) psycodict `copy_from` format.

Two pieces of good news:

- **The det-bound conventions already match.** Stevell's corank1 convention
  (rank n up to D/n, rank n+1 up to D/(n+1)) is exactly `run_all.py`'s `Dbound`,
  and the `12 ≤ r ≤ 20, nminus ∈ {1,2}` det extension to `D/(22-r)` (the
  `--nokthree` flag) exists precisely to feed the k3 computation.  The hook at
  `run_all.py:183` is where this stage goes.
- **LMFDB.jl is not needed.**  Stevell queries the live LMFDB only to get
  (gram, label) pairs per genus and to filter out non-enumerated genera.  We have
  all of that locally: `lattice_basic_data/<r>/<nplus>/<lattice_label>` (fields per
  `lat_basic.format`, incl. `gram`, `label`) and `genera_advanced/.../<genus_label>`
  (`class_number`, first field).  A local loader replaces the DB connection and
  sidesteps the "dev LMFDB.jl + manual patch" problem entirely.

## 2. Environment

> **Status 2026-07-16: DONE.**  juliaup + Julia 1.12.6 (arm64) installed via brew;
> `lattices/oscar/Project.toml` + `Manifest.toml` hold Oscar `master` (1.9.0-DEV,
> post-#5949).  Smoke test passes: `_primitive_embeddings_in_unimodular_safe`,
> `k3_lattice`, `ZZLatGluingCtx` all present.  Run with
> `julia --project=lattices/oscar`.

Local Julia was **1.5.3 with Oscar 0.5.1 (from ~2020) — unusable**.  Setup step:

```sh
brew install juliaup && juliaup add release        # ≥ 1.10 required; release is fine
cd lattices/oscar
julia --project=. -e 'using Pkg; Pkg.add(name="Oscar", rev="master"); Pkg.instantiate()'
julia --project=. -e 'using Oscar; @show Oscar.VERSION_NUMBER;
                      @assert isdefined(Oscar, :_primitive_embeddings_in_unimodular_safe)'
```

- Pin to a specific `master` SHA (post-PR [#5949](https://github.com/oscar-system/Oscar.jl/pull/5949),
  merged 2026-07-15) and **commit `Project.toml` + `Manifest.toml`** so every worker and
  any future machine resolves identical versions.  Switch to the registry version when
  Oscar 1.9 ships (ask Stevell for the expected date).
- First `using Oscar` precompile is 10–30 min; do it once *before* fanning out workers
  (a cold parallel start stampedes the depot).  Workers run with `JULIA_NUM_THREADS=1`;
  budget ~2–6 GB RSS each (Oscar+GAP), so `-j 10`–`12` on the 16-core/128 GB local box.
- lovelace is the production target once calibrated (per-session permission required).

## 3. Components to build (all under `lattices/oscar/`)

Keep Stevell's files **vendored verbatim** in `code/` and put our glue next to them,
so his future updates drop in cleanly.  Patches to his behavior live in the wrapper,
and are listed in §7 so we can upstream them.

1. **`load_local.jl`** — replaces LMFDB.jl.
   `load_genera(data_root; rank, det_range, even_only, nminus_in, classlimit)`:
   - walk `genera_basic/<r>/<nplus>/` for cheap prefilters (det, `is_even`,
     `discriminant_group_invs`, signature — all columns of `genera_basic.format`);
   - a genus qualifies only if it is **completely enumerated**: `genera_advanced`
     exists and #files in `lattice_basic_data` for it equals `class_number` > 0.
     (Stevell's functions silently assume `:representatives` is the *full* class list;
     feeding partial lists produces wrong/missing rows, so this gate is correctness,
     not just hygiene.)
   - build `integer_lattice(gram)` per class, `set_attribute!(:lmfdb_label)`,
     group into `ZZGenus` via `genus(first)`, `set_attribute!(:representatives)`.
   - lazy: parse grams only for genera passing filters; the corank1 codomain pool is
     loaded once per job.

2. **`driver.jl`** — CLI entry point run by GNU parallel, one invocation per job:
   `julia --project driver.jl --task {overlat|k3|corank1} --rank n --data-root .. --out-root ..
    --domain-labels-file F [--codomain-det-lo a --codomain-det-hi b] [--classlimit c]`.
   It `include`s `code/lmfdb_lattice_enum.jl` and wraps it with:
   - **unit = one domain genus** (for corank1: domain genus × codomain det-window).
     Output goes to `inclusions/<r>/<nplus>/<domain_lattice_label>/…` (LabelPath
     sharding — a flat `rank12/` with 10⁵ subdirs would hurt), written to a temp file
     and atomically renamed, then a genus-level `.done.<task>[.<window>]` marker.
     Killed jobs (timeout, OOM) lose only the in-flight unit; **resume = rerun, driver
     skips marked units**.  No duplicate rows by construction.
   - the det-window markers make **incremental widening** free: raising the corank1
     codomain cap later only computes the new windows.
   - log per-unit wall time + row count to a `units.jsonl` ledger (feeds calibration),
     and record units *skipped* by a limit with the reason (provenance: the schema
     promises no completeness, but we should know exactly what was attempted).

3. **`write_table.jl`** — port of `writing_file.jl`: walk `inclusions/`, emit
   `lat_inclusions.tbl` in real psycodict format (3 header lines, then one row per
   line — Stevell's version writes a blank line after every row), with the k3 label
   computed from our own data rather than hard-coded, plus a built-in lint pass (§6).

4. **`setup.sh`** + a preflight check in the run_all stage that fails fast with a
   pointer to it if `julia --project -e 'using Oscar'` doesn't work.

## 4. run_all.py wiring

- Generalize `parallel()` to take an arbitrary command prefix (it currently hardcodes
  `magma -b`), and pass GNU parallel `--timeout` for this stage.
- Move the embeddings stage **after** connect: it's the slowest and most optional
  stage, and nothing downstream depends on it — a crash or cap there shouldn't delay
  the core tables.
- `build_embedding_inputs()`: after enumeration, walk `genera_basic`/`genera_advanced`
  (so unit counts and det histograms are *exact*, no estimation needed), apply the
  per-task cutoffs, and batch units into jobs by predicted cost (proxy until
  calibrated: Σ class_number² per genus, det-sorted) — same philosophy as
  `--batch-mass` for fill.  Jobs land in `embed_jobs.txt`, labels via per-job files
  (label lists can exceed ARG_MAX).  Keep jobs ≥ a few minutes of work: each pays
  ~20–40 s of Oscar startup.

### New flags (defaults are placeholders until calibration, §5)

| flag | default | meaning |
| --- | --- | --- |
| `--embed-tasks` | `k3,overlat,corank1` | comma list; `--skip-embeddings` still skips all |
| `--embed-jobs` | `= --jobs` | worker count for this stage (memory-bound, see §2) |
| `--embed-timelimit` | 3600 | GNU parallel `--timeout` per job (units inside survive) |
| `--embed-classlimit` | 100 | skip genera with class_number above this (domain or codomain) |
| `--embed-k3-min-rank` | 14 | k3 cost explodes as rank drops (complement rank 22−r is what Oscar must enumerate); extend 14→12 in later waves |
| `--embed-k3-max-disc-ratio` | 4096 | det cap for k3 domains (relative to 22−r, like Dbound) |
| `--embed-overlat-max-disc-ratio` | `= max-disc-ratio` | overlat is genus-local and cheap-ish; also `--embed-overlat-max-prank` (default 6) to skip huge discriminant p-ranks (prefilterable from `discriminant_group_invs`) |
| `--embed-corank1-max-disc-ratio` | 512 | applies Dbound-style to both sides; the pair count is the volume driver |
| `--embed-corank1-max-rank` | 16 | highest domain rank for corank1 |
| `--embed-calibrate` | off | sampling mode, §5 |
| `--embed-plan` | off | dry run: print units/predicted core-hours/predicted rows per stratum and exit |

## 5. Tuning cutoffs (the actual hard part)

Cost drivers differ per task, so cutoffs are per-task:

- **k3**: dominated by the orthogonal-complement genus (rank 22−r, det up to the
  domain's): Oscar effectively enumerates complement classes, and Kneser-style
  enumeration cost explodes with rank.  r=20 (complement rank 2) is cheap; r=12
  (complement rank 10) is brutal.  ⇒ run top-down in rank, cap det harder as rank
  decreases.  Output volume is trivial (one count per domain lattice).
- **overlat**: per-genus; cost ∝ class_number × (aut-group / discriminant-group
  orbit machinery: `image_in_Oq`, `double_cosets` in perm groups of size |O(q_p)|).
  Prefilter on discriminant p-rank and classlimit; otherwise run wide.
- **corank1**: #pairs ≈ (#genera rank n) × (#genera rank n+1) after signature/parity
  filters, each pair a nontrivial gluing computation plus per-result isometry
  recognition against the class list (Plesken–Souvignier per candidate for definite
  codomains).  This is where both compute *and* table volume can run away
  (`orthogonal_complement` arrays = multiplicity × rank numbers per row).

**Calibration mode** (`--embed-calibrate`): stratify qualifying units by
(task, rank, det decade), sample K (~5) units per stratum, run with a short timeout,
fit per-stratum cost `a·det^b` and rows/unit from the `units.jsonl` ledger, then emit
`embed_calibration.json` + a projection table (core-hours and table rows per stratum
at the *current* flags).  `--embed-plan` reuses the fit to answer "what would this
run cost" before committing.  Given a target budget (e.g. `--embed-budget-hours`),
choose caps greedily by value/cost: k3 first (feature value, tiny output), overlat
second, corank1 by increasing rank and det.  This mirrors the existing
`--enum-timelimit` TODO ("calibrate based on how much time we want to spend") — same
ledger can eventually serve both.

**Volume guards**: the ledger records rows/unit as it goes; a per-domain cap
(`--embed-max-rows-per-domain`, default 10⁴) marks a unit `.overflow` instead of
writing an unbounded payload, and `write_table.jl` reports total rows + overflow list
before anything is uploaded.  Never truncate silently — decide per case (or store
counts with empty vector arrays; needs a schema decision, see §8).

Two structural notes that make cutoffs safer than they look:

- Everything is **restartable and widenable**: caps only bound *this wave*; a later
  wave with bigger caps computes exactly the delta (det-window markers).
- Skipped ≠ silent: every unit excluded by a cap is in the ledger with a reason, so
  we can publish precise coverage statements later even though `lat_inclusions`
  itself promises no completeness.

## 6. Validation

1. **Slice diff vs Stevell**: run his exact `script_1.jl` slice (rank 1 det ≤ 500,
   rank 2 det ≤ 250) through our loader/driver and diff row-for-row against output he
   provides from his setup (label-mapping caveat in §8).
2. **Independent Magma cross-checks** (we have Magma locally, and it's a genuinely
   independent implementation):
   - corank1, rank 1 ↪ rank 2 (small det): primitive vectors of norm m up to Aut(L),
     enumerated directly; compare multiplicity and complement coordinates.
   - overlat (small det): enumerate index-p integral overlattices M ⊆ L ⊆ M♯ by
     brute force over order-p subgroups of M♯/M, dedup by Aut(M) then isometry;
     compare (codomain label, multiplicity).
3. **Literature spot checks for k3**: Morrison's uniqueness (ρ ≤ 10 ⇒ multiplicity 1
   for NS-type lattices), Nikulin's 2-elementary classification; existence must match
   Nikulin's criterion.
4. **Lint pass in `write_table.jl`**: every label exists in the generated lattice
   set; rank arithmetic consistent (`domain_rank + quotient_rank ≤ codomain_rank`,
   overlat det ratio a perfect square = p², parity of domain vs codomain for even
   overlattices); array shapes (`len(orthogonal_complement) = multiplicity ×
   codomain_rank` for definite corank1, `quotient_gens` length matches multiplicity ×
   domain rank); no duplicate (domain, codomain, quotient_invs) keys.

## 7. Issues found in the zip (to patch in the wrapper, and to upstream)

1. `writing_file.jl`: `i = parse(Int, last(_s))` takes the **last character** of
   `"rank12"` → rank 12 parses as 2.  Breaks for every rank ≥ 10.
2. `writing_file.jl`: `k3_label = "22.16.1.0.1"` is hard-coded and doesn't match our
   labeling (`r.nplus.det...` gives `22.19.1.0.x` for II(19,3); 16 is the signature
   *difference*, presumably the label scheme of whatever is on beta today).  Compute
   it from our own labeler.
3. `writing_file.jl` writes a blank line after every data row — not valid
   psycodict/COPY input.
4. `lmfdb_lattice_enum.jl`: `m = parse(Int, first(x[1]))` (first character of the
   codomain label) in the prim_emb filename — same ≥ rank 10 problem; the filename
   doesn't need the rank at all.
5. ~~The `ok` flag returned by the `_safe` functions is discarded.~~  Resolved by
   reading the Oscar source: the `_safe` functions end with
   `return !isempty(results), results`, so `ok` just means "found any" and
   discarding it is fine.  **But** probing known-true embeddings per stratum
   exposed a real Oscar bug behind the fast path: the unimodular route silently
   loses even sublattices of ODD unimodular lattices
   (`primitive_embeddings(I2, <2>)` returns `(false, [])` although `Z(1,1)` is
   such an embedding), while the generic machinery answers the same pair
   correctly — even→even, odd→odd and even→odd coprime/generic all check out.
   `driver.jl`'s dispatcher therefore routes even-domain → odd-unimodular pairs
   through `_primitive_embeddings_generic_safe`; report the fast-path bug
   upstream (Oscar issue) and to Stevell.
6. Appends open/close the output file once per row; fine at his scale, wasteful at
   ours — buffered per-unit writes in the wrapper.
6b. The `Vector{ZZGenus}` wrapper `_primitive_embeddings_lmfdb(G1s, G2s, path)`
   crashes: it `append!`s the per-genus method's return value, which is `nothing`
   (results go through `save_data_lmfdb`).  The driver loops the per-genus method
   directly instead.
7. `label_top` in the overlat loop is undefined if no class matches by isometry
   (would crash with `UndefVarError`); make it a hard error with context.
8. Multiplicity-0 k3 rows: `_embeddings_in_K3_lattice` writes a row even when the
   count is 0.  A 0-row is a *non-embeddability assertion* — valuable but a schema
   semantics question (and a volume question): confirm intent.

## 8. Questions for the call with Stevell

1. **Labels**: his data uses labels queried from (beta) LMFDB — which scheme?  Our
   pipeline regenerates everything under the new `r.nplus.det.j.x` scheme, so labels
   must come from our local files (already the plan), and any of his precomputed
   `data/` dumps would need a label map to be reused.  Which (rank, det) ranges has
   he already completed — worth importing for cross-validation even if we recompute?
2. Semantics of `ok=false` from the `_safe` variants (issue 5 above), and of
   multiplicity-0 rows (issue 8).
3. Equivalence relation behind `multiplicity` (up to O(codomain)? O×O?) — needed for
   the knowl, and for the Magma cross-check to count the same thing.
4. Observed costs on his cluster: wall times per (rank, det) stratum, worst
   offenders, memory; what job granularity worked (he mentioned splitting by rank and
   det ranges — what sizes?).
5. Roadmap: composite-index overlattices / corank ≥ 2 primitive embeddings (the
   schema supports general `quotient_invs`/`quotient_rank`) — do we reserve the
   incremental-widening design for those too?
6. Indefinite codomains get empty `orthogonal_complement` (no isometry transport) —
   acceptable long-term, or plan to fill later?
7. `images` column: schema says probably delete — confirm and drop from the writer.
8. Oscar 1.9 release timing (when can we unpin from master), and whether we should
   upstream the wrapper fixes (§7) to his repo/zip.

## 9. Phasing

| phase | content | exit criterion |
| --- | --- | --- |
| 0. env + smoke (~½ day, mostly precompile) — **env half done 2026-07-16 (see §2)** | juliaup, pinned project, tiny local dataset via `run_all.py -R 6 -D 300` (Magma is local), run all three tasks on it | rows produced end-to-end, `lat_inclusions.tbl` lints clean |
| 1. correctness | loader + driver + writer, §6.1–6.2 validation on the script_1 slice + ~100 random Magma cross-checks | zero discrepancies, or each one explained |
| 2. scale hardening | atomic units, ledgers, resume, timeout/OOM behavior, overflow guards | kill -9 mid-run + rerun = identical table |
| 3. calibration | `--embed-calibrate` on the production listing, projection report, pick caps to budget | signed-off core-hours + row-count projection |
| 4. production + upload | run (local or lovelace — needs authorization), `write_table.jl`, upload with the other lat_* tables from the *same* generation run, coverage note | `lat_inclusions` live on beta, coverage documented |

## 10. Derived representation statistics: universality / regularity columns

`lat_advanced.format` has four placeholder columns (`connect_genus.m:803-809`), and
`lat_basic.format` has two booleans already partially computed in `fill_genus.m`.
Definitions (schema `lat_lattices.md`, definite lattices only; "represents" = any
isometric embedding, not necessarily primitive):

- `is_universal` / `is_even_universal` (bool): represents all positive (even) integers.
- `universality` / `even_universality` (smallint): largest n such that every (even)
  positive definite rank-n lattice embeds.
- `is_regular` (bool) / `regularity` (smallint): rank-1 / largest-n version of
  "represents everything its genus represents".

### The three-layer architecture

**Layer 1 — rigorous, per-lattice, criterion sets (Magma, fill + connect).**
The literature gives finite criterion sets for exactly n = 1, 2, 8 (existence for all
n is Kim–Kim–Oh 2005, but explicit sets are open for n = 3–7 and for all even-n ≥ 2):

| claim | criterion | source | status |
| --- | --- | --- | --- |
| 1-universal | represents S₁ = {1,2,3,5,6,7,10,14,15} | Conway–Schneeberger, proof Bhargava (Contemp. Math. 272) | published; integer-matrix forms = our lattices |
| all even integers | represents 2·S₂₉₀, S₂₉₀ = {1,2,3,5,6,7,10,13,14,15,17,19,21,22,23,26,29,30,31,34,35,37,42,58,93,110,145,203,290} | Bhargava–Hanke 290-theorem via halving the even sublattice | **preprint only (2011, still unpublished 2026)** — flag in knowl |
| 2-universal | represents S₂ = {⟨1,1⟩, ⟨2,3⟩, ⟨3,3⟩, [2,1,2], [2,1,3], [2,1,4]} ([a,b,c] = Gram [[a,b],[b,c]]; [2,1,2] = A₂) | Kim–Kim–Oh (Contemp. Math. 249, 1999); unique minimal set: Kominers arXiv:0708.4340 | published |
| 8-universal | represents S₈ = {I₈, E₈} | B.-K. Oh, PAMS 128 (1999); uniqueness Kominers arXiv:0807.2099 | published |

Useful facts: k-universal follows from n-universal for k < n; **any even sublattice
of L absorbs every embedded even lattice**, so even-universality of L equals that of
its even sublattice (this is why 2·S₂₉₀ tested on L itself is correct even for odd
L); n-universal requires rank ≥ n+3 for n ≤ 5 (Oh, PAMS 128), giving cheap
short-circuits; **no even-2-criterion set exists** (Ji–Kim–Oh, J. Korean Math. Soc.
58 (2021) classify rank-5 even-2-universal candidates with 20 of 55 still unproven).

Value logic for the smallint columns (paired with the booleans
`universality_proved` / `even_universality_proved` / `regularity_proved`):
maintain `lower` = largest n *proven* via the S₁ → S₂ → S₈ ladder, and `upper` =
smallest k with an explicitly non-represented rank-(k+1) witness (pool: I_k,
⟨1,…,1,7⟩ — the complement of any I_{k−1} ↪ I_m is unimodular hence I_{m−k+1},
and 7 is not a sum of three squares — A_k, D_k, E₆/E₇/E₈, A₁^k, [2,1,4]-type even
binaries).  When lower = upper the value is stored with `proved = true` (a
theorem).  Otherwise the column holds the best *heuristic* — the largest k whose
entire witness pool provably embeds — with `proved = false`, and the Layer-3
closure pass refines it against the actual embedding data.  Null means not even a
heuristic survived the budgets.  Curated override: L ≅ I_m (detected by norm-1
vectors spanning) has universality min(m−3, 5) for m ≥ 4, proved — Ko 1937
(I_{n+3} is n-universal for n ≤ 5) plus the classical fact that norm-2 vectors of
Zᵐ generate only A/D root systems, so E₆ never embeds and the value caps at 5.

`RepresentsLattice(L, M)`: exact embedding test by depth-first search — vectors of
prescribed norm via `ShortVectors`, subsequent basis vectors with prescribed inner
products via integer affine solve + `CloseVectors` on the constraint kernel, with a
node cap (cap hit ⇒ undecided, never a wrong answer).

**Layer 2 — rigorous, per-genus (Magma, connect).**  `is_regular`/`regularity`:
class number 1 ⇒ regular with regularity = rank (genus-representation = some class
represents, so h = 1 makes every rank trivially regular).  Otherwise scan m up to
the common stored theta precision: a value some classmate represents but L misses is
an unconditional counterexample ⇒ `<false, 0>`.  No counterexample ⇒ `\N`: a finite
scan cannot *prove* regularity (deciding it is genuinely hard — Watson/913 regular
ternaries, 14 still GRH-conditional per Lemke Oliver), so positive claims beyond
h = 1 wait for Layer 3 + curation.

**Layer 3 — db-relative, post-Oscar closure pass (the "embedding data" layer).**
Once `lat_inclusions` exists: every embedding M ↪ L factors as finite-index
overlattice steps (prime by prime) followed by a chain of corank-1 primitive
embeddings, so the transitive closure of the two edge types is exactly
db-representability — *except* that chain intermediates can have determinants
beyond the enumeration window, so closure-absence is not non-representability.
Consequences: closure gives sound "represents" facts and per-(L, n) *verified
coverage bounds* ("L represents every rank-n lattice with |det| ≤ B_n(L)"), not
sound negatives.  The pass (`refine_universality.py`) therefore: (a) surfaces
regularity **candidate counterexamples**: a genus-mate-reachable M not reachable
into L.  The genus-representation half is rigorous (the closure path is a witness
embedding), but "L misses M" is only closure-relative — each candidate is finished
off by a `RepresentsLattice(L, M)` exhaustion in Magma, whose decided-false output
makes the counterexample a theorem (that combination is how regularity values
beyond h = 1 get decided); (b) upgrades the heuristic values in the columns —
db-verified coverage raises `universality`/`even_universality`/`regularity` with
`*_proved = false`, and matching witness upper bounds flip `*_proved` to true.
(Resolved 2026-07-17: heuristics go in the value columns; the `*_proved` booleans
record which values are theorems.)

### Smoke-run findings (2026-07-16, `run_all.py -R 4 -D 120` locally)

Two **pre-existing** blockers surfaced while validating the wiring end-to-end
(neither is caused by, nor blocks, the universality code — its intrinsics are
unit-tested directly and the connect wiring was verified by round-tripping the
stored theta and aut_group strings):

1. **Definite enumeration is broken on local Magma V2.29-4**: every definite
   genus (all 628 in the smoke) dies in Magma's own
   `package/Lattice/Lat/neighbors.m:911` with `Invalid attribute 'IsNatural'`
   (the workaround at `fill_genus.m:293` doesn't cover the path Adam's
   GenusRepresentatives takes).  This is exactly the territory of the
   genus-enumeration PR (jvoight/k3s-lmfdb#71); once that lands, rerun the smoke
   to see the universality/regularity columns on real definite records.
2. **Connect crashes on non-unimodular indefinite lattices**: `D := Dual(L)`
   rebuilds via `Lattice(B, M)`, which requires positive definiteness — 690/692
   indefinite labels failed.  Fix: dual from the inverse Gram with
   `CheckPositive := false` (flagged as a separate task).

### Consistency + open questions

- `is_universal` ⟺ `universality ≥ 1` and `is_even_universal` ⟺
  `even_universality ≥ 1` must hold between lat_basic (fill) and lat_advanced
  (connect): both call the same intrinsics in `universality.m`; add a lint.
- Existing bug (fixed as part of this work): `fill_genus.m` concluded
  `is_even_universal := true` from theta support on 2..30 — valid only for
  *doubly even* lattices (where halving stays integer-matrix and the 15-theorem
  applies); in general the doubled 290-set up to 580 is required (A₂ already shows
  2..30 can't suffice: its halved form is integer-valued).  Positive conclusions
  also need theta precision 580, which the stored series (~150) lacks — extended
  via the budgeted `ThetaSeriesIncremental`, leaving `\N` when the budget runs out.
- Call questions: (i) is relying on the unpublished Bhargava–Hanke 290 proof
  acceptable for `is_even_universal` (knowl footnote)?  (ii) ~~Layer-3 storage
  semantics~~ resolved via the `*_proved` columns.  (iii) should `universality`
  for even lattices be 0 (they miss ⟨1⟩) or `\N`-by-convention?  Currently 0,
  which is what the schema wording implies.
