# MBRRecoverySearch — Design + Progress

**Branch**: `feature/mbr-recovery-search`
**Base**: `develop` @ `1eb5c1886`

A new SearchMethod that targets B-bucket cells — Pioneer's biggest source of
missed IDs vs. DIA-NN. For every donor pid (q ≤ 0.01 in ≥1 file) missing from
a receiver file's `main_search_psms`, runs a targeted PMM deconvolution in the
file's calibrated RT window around the donor's projected RT, computes the 7
existing MBR_*_true features inline, and emits a per-file recovery sidecar.

## Motivation

From per-file attribution on HUMAN PEP < 0.01 precursors (set A = 58,134
precursors × 23 Olsen Exploris files, 1,337,082 cells):

  Stage where Pioneer-missed cell is lost   Cells at default 0.03 bitvec  %
  Bitvec rejected (no PSM)                  138,646                       61.5
  Scored, PEP > 0.9 (MainSearch filter)      77,041                       34.2
  Survived MainSearch, lost in scoring/FTR    9,752                        4.3

Loosening bitvec from 0.03 → 0.01 recovers most of the B bucket
(138,646 → 68,727) at +40% runtime, but ~half of the recovered cells trip
the PEP > 0.9 filter (B shifts to A, not to "both"). This new pathway targets
the B-bucket directly with bounded cost (~870K cells × few hundred µs each
≈ 100-200 s on 8 threads, +15-25% runtime) and bypasses both the bitvec gate
AND the PEP > 0.9 gate by using donor-driven evidence as the seed.

## Architectural pivot during implementation

**Original plan** (V1 of MBR_RECOVERY_PLAN.md): split PrecursorScoringSearch
into PrecursorPrescoreSearch (pass-1 LGBM + qval) and PrecursorScoringSearch
(MBR + FTR + final FDR), insert MBRRecoverySearch between them.

**Revised plan** (now in code): MBRRecoverySearch runs AFTER
PrecursorScoringSearch. Two reasons:

1. The post-PrecursorScoringSearch file
   `temp_data/main_search_psms/{file}.arrow` already has 117 columns
   including `trace_prob`, `weight`, `log2_intensity_explained`, `irt_obs`,
   `irt_pred`, `rt`, `log_by_ratio_m0`, `frag{1..8}_smoothed_intensity`,
   `target`, `cv_fold` — everything we need to identify donors and shape
   `_MBRDonorEntry` payloads.

2. The `q_value` column in those files is unpopulated (all zeros) at this
   stage. We recompute per-file qvals from `trace_prob + target` via
   `get_qvalues!` ourselves — minimal overhead.

This avoids the riskier refactor of splitting PrecursorScoringSearch. The
recovery FTR (matching recovery rows with counterfactual partners and
training a paired LightGBM) becomes a separate V1.5 follow-up rather than
a built-in part of V1.

## Pipeline position

```
... ParameterTuning → QuadTuning → BitVecCalibration → MainSearch
    → PrecursorScoringSearch         (existing — unchanged)
    → MBRRecoverySearch              (NEW — V1: emit recovery sidecar)
    → IntegrateChromatogramSearch    (will read recovery sidecar in V1.5)
    → ProteinInference / ProteinScoring / MaxLFQ
```

## Module layout (current state)

```
src/Routines/SearchDIA/SearchMethods/MBRRecoverySearch/
├── PLAN.md                  This file
├── types.jl                 [DONE] RecoveredSeedRow, params, results, constants.
                                    Results now caches donor_pool / donor_pid_set /
                                    fragment_keys + initialized flag (built once,
                                    survives reset_results!).
├── MBRRecoverySearch.jl     [V1.0] process_file!-driven per-cell extraction.
                                    ensure_recovery_initialized! builds the donor
                                    pool lazily on first file; extract_file_recovery_seeds
                                    threads over missing pids (per-thread scratch);
                                    summarize_results! is the final log.
├── donor_pool.jl            [DONE] build_recovery_donor_pool + helpers
├── mbr_features.jl          [DONE] compute_seed_mbr_features (7 MBR_*_true)
├── sidecar.jl               [DONE] write/read recovery_seed_sidecar.arrow
└── extraction.jl            [V1.1] extract_max_weight_in_window! — iRT-space
                                    window + primary weight + scan_idx ✓;
                                    secondary features now REAL: the walk uses
                                    FusedStandard (not FusedRTIndexed) so
                                    apply_main_scoring! fills the donor's
                                    b/y-ion + per-rank frag intensities, from
                                    which log2_intensity_explained, log_by_ratio,
                                    and smoothed_frag_sqrt are computed at the
                                    best scan (Step 1 done).
```

Registered in `src/importScripts.jl` as an ordered-include block; excluded
from the SearchMethods walkdir loader.

## What's done (checkpoint after V1.0 wiring)

- [x] Module skeleton, types, constants — compiles + instantiates
- [x] `list_post_scoring_psm_files`, `build_recovery_donor_pool`,
       `main_search_pid_set`, `_read_smoothed_frag_sqrt`
- [x] `write_recovery_seed_sidecar`, `read_recovery_seed_sidecar`,
       `recovery_seed_sidecar_path`
- [x] `ReceiverExtraction` struct + `compute_seed_mbr_features` (formulas
       mirror PrecursorScoringSearch/mbr_streaming.jl `_compute_mbr_inner!`
       exactly; verified against synthetic data — see commit message of
       previous push for sample output)
- [x] `extract_max_weight_in_window!` — reuses
       `collect_rt_window_precursors!` → `prepare_scan_peaks!` →
       `run_fused!(FusedRTIndexed, ...)` → `solve_deconvolution!` chain.
       Donor pid is manually prepended to the per-scan candidate list since
       it's NOT in `precursors_passing` (it's the missing pid we're recovering).
       **iRT-space window** (see decision note below): scans are visited when
       `abs(rt_irt_model(scan_rt) − donor.irt_obs) ≤ irt_tol`.
- [x] **V1.0 per-cell loop wired into `process_file!`** (Step 2):
       (1) `ensure_recovery_initialized!` builds the donor pool + fragment keys
           once, lazily, cached on results (survives `reset_results!`)
       (2) per file: map ms_file_idx → `main_search_psms/{name}.arrow`, compute
           missing-pid set, load rt_index from disk, thread extraction over
           missing pids using `getSearchData(ctx)[thread_id]` scratch
       (3) compute 7 MBR_*_true features inline, build RecoveredSeedRow, write
           sidecar; update counters (n_donor_pids, n_cells_attempted/emitted)
       (4) Pioneer precompiles + loads clean with all changes.

       NOTE: compilation ≠ correctness. The loop has not yet been *run* against
       real data — it needs upstream `main_search_psms/*.arrow` + `rt_indices/`
       that only exist mid-pipeline, so its first real exercise is Step 4 (e2e),
       which also requires Step 3 (registration — not yet in the pipeline).

## What's next

In order:

1. **[DONE] Finish extraction.jl secondary features**
   Approach taken (differs from the original Hs-column idea): ion type (b/y)
   is NOT stored in `Hs`, only `rank`. Rather than reconstruct ion types, the
   window walk now runs `run_fused!` with **`FusedStandard`** (was
   `FusedRTIndexed`) backed by `getMainUnscoredPsms(scratch)`, so
   `apply_main_scoring!` populates the donor column's `MainUnscoredPSM`
   (b_int, y_int, frag1..8_int) for free. Quant weights are unchanged
   (identical column set + matched fragments; the inline mz/iRT prec-checks
   are no-ops here — `irt_tol=Inf`, and every enumerated candidate already
   passes the quad mz window). At each best-weight scan we compute, mirroring
   MainSearch's `Score!` (ScoredPSMs.jl:199, 242-243) exactly:
   - `log2_intensity_explained = log2((b_int+y_int) / Σ scan intensity)`
   - `log_by_ratio (m0)        = log(b_int+1) − log(y_int+1)`
   - `smoothed_frag_sqrt[8]     = sqrt(L1-normalized frag{1..8}_int)`
     (single-scan: MainSearch's radius-1 chromatogram smoothing collapses to
     the raw per-rank value when a scan has no flanking cycles)
   Also fixed a latent skeleton bug: `n_frag_isotopes` was `UInt8(1)` but
   `run_fused!` wants `::Int64` (would have MethodError'd on first call).
   Pioneer compiles clean with the change.

2. **Per-cell extraction loop in performSearch!** (~half day, depends on
   how clean the SearchData scratch lifecycle access is from MBRRecoverySearch)
   - For each (donor pid × file) cell:
     - Pull thread-local `SearchDataStructures` scratch from search_context
     - Read per-file rt_index from `temp_data/rt_indices/{file}.arrow`
       (already written by PrecursorScoringSearch's `build_rt_indices!`)
     - Compute `center_rt = inverse(rt_irt_model)(donor.irt_obs)`
     - Get `rt_tol` from `getRtTolerance(ctx, ms_file_idx)` (per-file
       calibrated; floor at `params.rt_tol_floor_min`)
     - Call `extract_max_weight_in_window!`
     - If positive weight: call `compute_seed_mbr_features` and push
       RecoveredSeedRow into per-file vector
   - Parallelize per file (each has its own SearchData scratch) or per
     pid within a file (per-thread scratch).

3. **Wire into pipeline** (~15 min)
   - Add `("MBR Recovery", MBRRecoverySearch())` to the `searches` array
     in `src/Routines/SearchDIA.jl`, after `("Precursor Scoring", ...)`

4. **End-to-end test on 23 Olsen Exploris files** (~half hour, plus a
   ~15-min search run)
   - Same config we've been using:
     `/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/searches/olsen_exploris_3P_bitvec10_2026-06-29/config.json`
   - Validate per the metrics in the "Validation plan" section below.

5. **V1.5 — Recovery FTR** (deferred)
   - Pair recovery rows with cross-fold counterfactual partners using
     existing `build_counterfactual_partner_pools` machinery
   - Train a small paired LightGBM on (MBR_*_true vs MBR_*_false) for
     recovery rows only
   - Apply FTR qval threshold to filter recovery rows
   - Update `IntegrateChromatogramSearch` to ingest recovered passing rows

## Data structures (live in `types.jl`)

`RecoveredSeedRow` mirrors `_MBRDonorEntry`'s payload schema plus the 7
precomputed MBR_*_true features plus markers (`mbr_recovered_seed = true`,
`trace_prob = 0`).

`MBRRecoverySearchParameters` carries: donor_qval_threshold (default 0.01),
max_donors_per_pid (default 2), rt_tol_floor_min (default 0.05), isotope_err_bounds
((1,0) matching MainSearch), min_fraction_transmitted (from config).

`MBRRecoverySearchResults` is a mutable counter: n_donor_pids,
n_cells_attempted, n_cells_emitted, elapsed_sec.

## Validation plan

After V1 implementation finishes:

- **Headline**: Pioneer 0.03 + MBRRecoverySearch per-file precursor IDs
  vs DIA-NN per file. Current 0.03 baseline is 92.8 % of DIA-NN. Target: ≥ 96 %.
- **Stage attribution rerun**: B-bucket count drops from 138,646 to ≤ 80,000.
  DIA-only-B count drops from 51,679 to ≤ 25,000.
- **Cross-tool ground-truth check**: of recovered cells whose pid maps to
  a DIA-NN ID in the same file (canonical (seq, charge, mods) key match),
  recovery rate ≥ 60 %.
- **FDR sanity**: at recovery qval ≤ 0.01 (V1.5), entrapment/species-mismatch
  rate on recovered rows ≤ 1 %.
- **Runtime**: total search time ≤ 17 min (current 14.7 + ≤ 2 min for new stage).

## Validation artifacts that already exist

In `/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/searches/olsen_exploris_3P_bitvec10_2026-06-29/`:

- `FINDINGS.md` — full analysis log: bitvec sweep + mass-tol sweep + stage
  attribution + PG-precursor distribution + per-file Pioneer-vs-DIA-NN tables
- `compare_pioneer_diann.jl`, `compare_human_pep01.jl`, `stage_attribution.jl`,
  `bitvec_sweep_modkey.jl`, `pg_precursor_dist.jl`,
  `pg_per_file_peptide_support.jl`, `diann_compare_sweep.jl` — Julia
  analysis scripts that produced the numbers above

Recovery-validation scripts will be added there as the V1 implementation
finishes.

## How a future Claude session picks this up

1. Read this file first to understand the design + current status.
2. Read the existing code:
   `src/Routines/SearchDIA/SearchMethods/MBRRecoverySearch/{types,donor_pool,MBRRecoverySearch}.jl`
3. Read the integration targets (the files listed under step 3 above) to
   plan the extraction.jl implementation.
4. Existing MBR machinery to reuse:
   - `_MBRDonorEntry` struct: `PrecursorScoringSearch/mbr_streaming.jl:36`
   - `_insert_sorted_donor_entry!`: `mbr_streaming.jl:301`
   - `_mbr_hellinger`: `mbr_streaming.jl:131`
   - `MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT`: `mbr_streaming.jl:55`
   - 7-feature list: `mbr_ftr.jl:67-78`
   - Counterfactual pool builder: `mbr_pairing.jl:152`
5. Existing chromatogram extraction utilities to reuse (for extraction.jl):
   - `collect_rt_window_precursors!`: `IntegrateChromatogramsSearch/utils.jl:577`
   - `get_rt_tol`: same file
   - PMM solver: `src/utils/ML/spectralLinearRegression.jl`
   - Unified scan loop dispatch: `src/Routines/SearchDIA/process_scans.jl`
6. The user prefers concrete code snippets + concise prose over pseudo-code.
   Re-read `FINDINGS.md` in the search workspace for the data context that
   motivated this work.

## Reasoning notes for the future-self / reviewer

- We use `trace_prob` not `q_value` from the post-scoring file because
  `q_value` is unpopulated there. Compute per-file qval via
  `get_qvalues!(trace_prob, target, qvals; doSort=true)`.
- Decoys naturally come with qval ≤ 0.01 (1% FDR ≈ ~1% decoys at that
  threshold). Existing MBR's FTR pairing logic handles the imbalance via
  cross-class counterfactual matching, so we don't separately enumerate
  decoy donors.
- `trace_prob = 0` for recovery rows is a deliberate signal to a retrained
  FTR LightGBM. V1.5 will add `mbr_recovered_seed` as an explicit boolean
  feature to disambiguate.
- We chose PMM over Huber for extraction (matches MainSearch's weight
  semantics + faster) per user direction. PMM iters/tol use the production
  constants `DECONV_MAX_ITER` (1000) / `DECONV_CONVERGENCE_TOL` (0.01).
- **Window (REVISED — supersedes the original RT-minutes plan): iRT space.**
  The original design opened the window at `center_rt = inverse(rt_irt_model)(donor.irt_obs)`
  ± `getRtTolerance` (RT minutes). That depends on an iRT→RT inverse model that
  **does not exist** — `irt_rt_map` / `setIrtRtMap!` is dead code (never called
  anywhere in `src`), and no `inverse(::RtConversionModel)` is defined. Instead
  we window in iRT space: a scan is visited when
  `abs(rt_irt_model(scan_rt) − donor.irt_obs) ≤ irt_tol`, where `irt_tol =
  getIrtErrors(ctx)[ms_file_idx]`. This uses only the forward model (always
  available) and matches how the fragment-index search matches scans to
  precursors across files (`_precompute_scan_properties`,
  PartitionedFragmentIndex/search.jl:497). The per-scan *competitor* enumeration
  still runs in RT space via `get_rt_tol(rt_binned_tol, scan_rt)` (with an
  irt_tol/slope fallback when no rt_binned_tol exists), unchanged from
  IntegrateChromatogramSearch's `build_chromatograms`.
- **Architecture (REVISED): per-file work lives in `process_file!`, not
  `summarize_results!`.** The framework's `execute_search` loop already
  iterates `(ms_file_idx, spectra)` and calls `process_file!` per file, so that
  hook is the only place the receiver's `spectra` + per-file calibrated models
  are available. The donor pool is file-independent, so it is built once
  (lazily, on the first `process_file!`) and cached on the mutable results
  struct, which persists because `reset_results!` is a no-op. The old plan to
  do everything in `summarize_results!` could not get `spectra` and is why V0
  stalled.
- Cell enumeration: skip cells where donor pid is already in receiver's
  `main_search_psms`. Already-IDed pids stay available as competitors in
  the per-scan design matrix via `rt_index` lookup.
