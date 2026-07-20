# Sciex ZT Scanning-DIA — Session Handoff

**Purpose:** self-contained resume guide so another Claude (or the same account on a
different machine) can pick up this branch **without access to the local `~/.claude`
memory** — this document is the source of truth. Last updated end of the 2026-07 session.

---

## 0. TL;DR / current state

- **Branch:** `feat/sciex-zt-scanning` (off `develop`). Pushed to `origin`
  (`github.com/nwamsley1/Pioneer.jl`). Was developed in worktree `pio-sciex-zt`.
- **HEAD at handoff:** `0fc545743` — `perf(sciex-zt): cut collapse time 78% / alloc 89%`.
- **Current best result** (3-file `arrow_3rep`, k=6, MBR off): **79,056 precursor rows /
  12,829 protein-group rows.**
- **Everything is committed and byte-identity-verified.** The work below is DONE; the
  BACKLOG (§5) is analyzed but NOT implemented.

## 1. What "ZT scanning DIA" is (context for a fresh instance)

Sciex ZenoTOF 8600 **Zeno-Trap scanning DIA**: instead of fixed wide isolation windows,
Q1 scans continuously and ProteoWizard rebins it into contiguous ~1 m/z bins. A precursor
is therefore spread across **2k+1 consecutive Q1 bins within one cycle** (measured bin step
**S ≈ 1.022 m/z**; effective transmission width **W ≈ 7 m/z**, NOT the filename's "5"; the
empirical weight profile decays to baseline by **±6 bins ⇒ k=6**). One cycle ≈ 496 MS2
windows, so a precursor's ±k bins are consecutive scan indices.

**Model:** a "meta-scan" = the center bin (where precursor m/z falls in its ±S/2 window)
± k neighbor bins **in the same cycle**. After deconvolution, `collapse_to_metascans`
(`src/Routines/SearchDIA/SearchMethods/MainSearch/metascan_collapse.jl`) reduces the raw
per-(precursor,scan) PSMs to **one meta-PSM per (precursor, cycle)** — keeps the center
row's features + adds the 2k+1 raw weights (`metascan_w_*`) + shape features.

Env flag **`PIONEER_ZT_METASCAN_K=6`** turns the whole ZT path on (0/unset = conventional
DIA, byte-identical). Sample used = 3-proteome HYE (human/yeast/ecoli). Library:
`/Users/nathanwamsley/Data/SPEC_LIBS/3P_May7-2026_bitvec_10.poin`.

**Public data provenance:** PRIDE **PXD070049** (ZenoTOF 8600 LFQ benchmark); preprints
bioRxiv `2026.01.29.702266` (LFQ benchmark, the HYE Condition_A files) and
`2026.03.12.711261` (ZT-Scan-DIA method paper).

## 2. Commit arc this session (all on top of the earlier ZT support)

| Commit | What |
|---|---|
| `cf30cc4bd` | Refactor: run develop's chromatogram features **post-collapse** on the clean meta trace; dropped the bespoke "elution" impl (after collapse, develop's per-precursor aggregation IS the elution computation — same names/math). Auto-fixed `n_scans` (was a 13×-cycle artifact) + made chromatogram features nearly free (0.7s on 3.2M vs 35M). |
| `f641960d2` | Diag: split collapse vs chromatogram timing. |
| `ba322a4a1` | Docs: Part-3 efficiency plan + the rejected "A1" experiment. |
| `0104d129d` | **`zt_tri_pcor`** — mean-centered Pearson(weight profile, triangle). Centering fixes the flat-decoy compression of the uncentered `zt_tri_cosine`, and the weight-profile-vs-triangle signal was never in the 2nd-pass model → top-3 there. **+2.9% prec / +3.0% PG (3-file).** |
| `87b84259b` | Promote `zt_tri_cosine` + `zt_entropy` to the 2nd-pass model (`ADVANCED_FEATURE_SET`); pruned the other 7 profile descriptors (noise). **+1.3%/+0.5%** on top. |
| `9fa9272bf` | Diag: `PIONEER_DUMP_METASCANS` (meta-PSM dump) + `frag*_int` added to `PIONEER_DUMP_RAW_PSMS`, enabling the correctness check. |
| `0fc545743` | **Perf: collapse 30.9s→6.6s (−78%), 34.5GB→3.9GB (−89%)**, byte-identical. See §3. |

Cumulative feature gain over the shape-only baseline: **+4.2% prec / +3.5% PG.**
Docs on branch: `zt_part3_efficiency_plan.md`, `zt_shape_feature_experiments_plan.md`,
`zt_chromatogram_feature_redesign.tex`.

### Tested-and-rejected (don't redo)
- **A1** (drop `frag_corr_best_shape` + its consensus loop): −1.6% PG for only 2.3s. The
  consensus loop is cheap; the shape cost is the irreducible fragment gather. Reverted.
- **Summed-fragment aggregate** (`frag_sum_corr_{weight,tri}_shape`): 3-file neutral
  (−0.1%/+0.04%). Shape already covered by `frag_corr_strength_shape` + `zt_tri_pcor`.
  Reverted.

## 3. The collapse perf fix (how, so you understand the byte-identity guarantee)

Profiling (`PIONEER_DIAG_ALLOC` + `PIONEER_DIAG_PSM_FOOTPRINT`) showed `metascan_collapse`
allocated **34.5 GB/file** — as much as the whole deconvolution — for a 7GB→1.2GB
reduction. Three type-stability/copy fixes in `metascan_collapse.jl`, all byte-identical
(3-file result stayed EXACTLY 79,056/12,829):
1. **m/z accessors** (`getMz`→`Arrow.Primitive`, `getCenterMzs`→`Vector{Union{Missing,Float32}}`)
   boxed on 35M indexed accesses in the per-row guard → materialize concrete `Vector{Float32}` once.
2. **`fcols` typed as concrete `NTuple{8,Vector{Float32}}`** (was `Union{Nothing,NTuple}`,
   boxing `fcols[b]` on ~680M inner-gather iters) — the DOMINANT fix.
3. **Sort only the 11 read columns** via `_permute_f32`; pull meta rows from the original
   `psms` via `perm[center_rows]` instead of copying the full ~60-col table.

## 4. How to run a ZT search

Standard invocation (14-core machine; **never run two searches in parallel**):
```
env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV -u CONDA_PROMPT_MODIFIER \
  PIONEER_ZT_METASCAN_K=6 \
  julia --startup-file=no --project=. --threads 10 --gcthreads 5,1 \
  -e 'using Pioneer; SearchDIA("<config.json>")'
```
The `env -u CONDA_*` prefix avoids a conda/OpenMP abort (OMP Error #13). MBR is auto-
skipped for single-file runs; multi-file works. The 3-file config used all session was
`lwn_baseline/cfg_3rep_new.json` (points at `arrow_3rep`, MBR off, q=0.01).

## 5. BACKLOG (analyzed, NOT implemented — priority order)

### 5a. Dead post-deconv column removal — GENERAL win (all datasets), audit COMPLETE
The post-deconv PSM table is **35M rows × 61 cols ≈ 7 GB** and is the peak-RSS driver on
EVERY DIA search (not ZT-specific). Audit (6 parallel agents + verification, 0 production
read-sites confirmed) of `MainSearchScoredPSM` (`src/Routines/SearchDIA/PSMs/ScoredPSMs.jl`):

- **DEAD — remove field + computation:** `max_matched_residual`, `max_unmatched_residual`
  (both computed in `spectralDistanceMetrics.jl` ~L152-158,167-168; the scoring LOOP stays —
  needed for gof/hellinger/manhattan), `matched_ratio` (inline `ScoredPSMs.jl:245-246`,
  field `:138`; only the optional `PIONEER_DUMP_RAW_PSMS` dump reads it via `intersect` →
  safe), `best_rank` (field `:115`, assign `:227`; removing it likely makes
  `fusedMatch.jl`'s UnscoredPSM.best_rank compute dead too — cascade).
- **DEAD as stored FIELD, keep the LOCAL compute:** `total_ions_iso` (field `:76,197`; local
  at `:190` feeds `poisson`), `b_count` (field `:73,194`; local at `:189` feeds `total_ions`).
- **ALIVE — do NOT remove:** `error` — `features.jl:61,97` reads it to compute `err_norm`
  (a PRESCORE feature). Optional later opt: fold error→err_norm into `Score!`, drop `:error`.
- **NOT candidates:** the 24 frag-intensity cols (`frag*_int`, `fitted_frag*`, `shadow_frag*`)
  are ALIVE (feed `smoothed_2d_shadow_hellinger`).
- **Also mirror removals in `TuningScoredPSM`** (same file, ~L270-292, 335-357).
- Est. ~0.3–0.5 GB/file + dead compute. **Verify byte-identical IDs on ZT (79,056/12,829)
  AND a non-ZT dataset (Astral/Exploris) after removal.**
- Lead: `:tic` (listed with `:matched_ratio` in the `features.jl:1608` "computed-but-
  unwired" comment) — not a struct field; check separately.

### 5b. Batched / streaming collapse — the PEAK-MEMORY prize
Peak RSS (~27 GB) is dominated by the 35M×61 table (collapse churn is already fixed).
**Key property: the collapse is LOCAL** — each meta-PSM needs only ±k=±6 scans of one
precursor within one cycle (≤13 consecutive scan indices). So it can run **per contiguous-
scan batch with ±k overlap** + boundary center-ownership dedup, collapsing BEFORE
`library_search`'s `result = vcat(fetched...)` (`LibrarySearch.jl:227`) → vcat 3.2M meta-
PSMs instead of 35M raw rows. Peak ~7 GB → ~1.2 GB. Real refactor of the deconv→collapse
boundary; the boundary-overlap correctness is the risk; gate ZT-only; verify byte-identical.
Current partitioning is round-robin (needs contiguous batches). Also naturally parallel.

### 5c. Shadow-spectrum empirical spectral library (idea from this session)
Build an **empirical, interference-removed** `.poin` library from main-search results:
best-scoring PSM per precursor across runs → its **shadow spectrum** (the deconvolved-
observed, interference-subtracted per-fragment intensities) → library entry.
- **Shadow = `shadow_frag*_int`** in `SpectralScoresMainSearch` (`spectralDistanceMetrics.jl`):
  `shadow_peak = fitted_peak − r` where `r = Σ_all fitted − observed`, i.e. `observed − Σ_(other
  precursors) fitted`, clamped ≥0. (NOT `frag*_int`, which is raw matched;
  `frag*_smoothed_intensity` on disk is smoothed `frag*_int`, also not shadow.)
- **BLOCKER:** shadow cols are dropped at `scoring.jl:1206` (`select!(... Not([FITTED...,
  SHADOW...]))`) after feeding `smoothed_2d_shadow_hellinger`, so they never reach
  `temp_data/main_search_psms`. Step 1 = preserve them (a per-file sidecar like the existing
  `*.prec_prob.sidecar.arrow` pattern).
- Pieces: best-PSM-per-precursor aggregation (scores in the prec_prob sidecars / use
  `passing_psms`); map shadow rank→fragment m/z+annotation via the original library's ranked
  fragments (deconv `rank_at` = library order); write `.poin` reusing `BuildSpecLib`
  (`src/Routines/BuildSpecLib.jl`). Bonus: empirical iRT from the best PSM.
- Design calls: only top-8 fragments; apex-scan vs peak-integrated shadow (single scan is
  noisy); decoys (empirical only exists for targets → regenerate); coverage = detectable
  subset; per-precursor normalization.

### 5d. NCE-criterion A/B (idea from this session)
NCE tuning (`ParameterTuningSearch.jl:404-437`) sweeps NCE 21→40 (20 pts), deconvolves at
each, and picks per precursor the NCE with highest **`gof`** (deconvolution residual fit).
`gof`/shadow both come from residual `r`, so it implicitly picks the NCE where library@NCE
best explains the observed. Could try a more **direct** criterion — `:fitted_hellinger`
(fitted-vs-shadow, already computed on those PSMs) or a normalized shadow-vs-library cosine —
instead of `:gof`. Robust to candidate incompleteness for *relative* NCE ranking (same
candidate set across the grid); more sensitive to it for absolute shadow quality.

### 5e. Low-value leftovers
Trapezoid/empirical shape template + spectral-angle (likely marginal — profile is already
triangle-like). See `zt_shape_feature_experiments_plan.md`.

## 6. Env-flag reference (ZT + diagnostics)

| Flag | Effect |
|---|---|
| `PIONEER_ZT_METASCAN_K=6` | Enable the ZT meta-scan path (0/unset = conventional DIA, byte-identical) |
| `PIONEER_ZT_PROFILE_FEATURES` | Default on; appends the `zt_*` profile features to the main-search LGBM |
| `PIONEER_ZT_SHAPE_EXP` | Env-gate for shape-feature experiments (`ZT_SHAPE_EXP_FEATURES` vector, currently empty) |
| `PIONEER_DIAG_ALLOC` | Print per-`@alloc_bucket` allocation (+ fine-grained `collapse.*` prints) |
| `PIONEER_DIAG_PSM_FOOTPRINT` | Print post-deconv table rows×cols×bytes |
| `PIONEER_DUMP_RAW_PSMS=<dir>` | Dump raw pre-collapse PSMs (incl. `frag*_int`) |
| `PIONEER_DUMP_METASCANS=<dir>` | Dump the collapsed meta-PSM table (metascan_w_* + all ZT feats + target) |

## 7. Correctness / byte-identity harness

The session's correctness check recomputed every ZT shape/profile feature from raw inputs
and matched Pioneer **bit-exactly**. To reproduce: run with `PIONEER_ZT_METASCAN_K=6
PIONEER_DUMP_RAW_PSMS=<d> PIONEER_DUMP_METASCANS=<d>`, then independently recompute (a)
weight-profile feats from `metascan_w_*`, (b) fragment shape feats by reconstructing the
per-cycle profiles from the raw `frag*_int` (group by precursor, window scan_idx ∈ [c−6,c+6]).
The general byte-identity gate for any perf/removal change: 3-file ZT result must stay
**79,056 / 12,829** (LGBM is deterministic → identical features ⇒ identical count), and a
non-ZT dataset must be unchanged.

## 8. Data situation (IMPORTANT on a new machine)

- The local `arrow_3rep/` (converted 3-rep data) and all result dirs were **deleted to free
  space**. Paths were under `/Users/nathanwamsley/Data/SciexZT/` (won't exist / will differ
  on another machine).
- **Source of truth = RIS mzMLs:** `/Volumes/d.goldfarb/Active/RIS_Goldfarb_Lab/NTW/
  LFQ_ZenoTOF8600_ZTScanDIA_5Da_5min_50ng_Condition_A_REP{1,2,3}.mzML` (7.5 GB each). Mount
  RIS on the new machine to access.
- **Regenerate arrows:** Pioneer `convertMzML` (v0.7.5-dev) on the 3 mzMLs → an output dir,
  ~5 min/file, deterministic. Then point a config's `ms_data` at that dir.
- **For pure code work** (5a dead columns, 5b batched collapse) a non-ZT dataset is enough
  for the byte-identity check; ZT data is only needed for ZT-specific validation.

## 9. Notes for the resuming instance

- Commit only when asked; commit messages end with the project's `Co-Authored-By` trailer.
- ZT changes must stay byte-identical when `PIONEER_ZT_METASCAN_K` is unset (non-ZT path)
  and additive for existing features — the `hasproperty` filter in both LGBMs
  (`scoring.jl:124`, `pass1_oom.jl:78`) makes absent columns safe.
- North star (from the design notes): eventually UNIFY the ZT and conventional paths by
  parameterizing on S (center-m/z step) + W (window); conventional DIA = S=0 ⇒ k=0 where
  every ZT stage degenerates to the current path.
