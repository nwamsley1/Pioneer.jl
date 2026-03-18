# SecondPassSearch Pipeline

## Overview

SecondPassSearch is stage 6 of the 9-stage SearchDIA pipeline. It performs exhaustive precursor identification using calibrated parameters from earlier stages, applies a multi-model prescore cascade to filter low-quality precursors, and produces fully-featured PSM Arrow files for downstream LightGBM scoring.

The pipeline has three phases:
1. **Per-file search + prescore** (`process_file!`) — Huber deconvolution, per-file prescore filter, re-deconvolution
2. **Per-file feature computation** (`process_search_results!`) — full feature pipeline, fold-split Arrow output
3. **Cross-file global aggregation** (`summarize_results!`) — aggregate prescores globally, filter or re-deconvolve

---

## Phase 1: Per-File Search + Prescore (`process_file!`)

**File**: `SecondPassSearch.jl`, lines 329–378

For each MS file, `process_file!` does three things:

### 1.1 Load Fragment Index Matches

Loads pre-computed fragment index matches from disk (`load_fragment_index_matches`). This gives:
- `scan_to_prec_idx`: per-scan vector of ranges into `precursors_passed`
- `precursors_passed`: flat `Vector{UInt32}` of precursor IDs

These were built by earlier pipeline stages (FirstPassSearch/FragmentIndexSearch) to record which precursors have fragment evidence in which scans.

### 1.2 Initial Deconvolution (`perform_second_pass_search`)

**File**: `utils.jl`, lines 135–176

Partitions scans across threads, calls `process_scans_fragindex!()` in parallel. Each thread:
- For each scan, builds a sparse design matrix with columns = precursors present in that scan
- Solves the Huber-loss deconvolution to get per-precursor weights (abundance estimates)
- Records per-PSM features: weight, gof, spectral_contrast, fitted_manhattan_distance, scribe, y_count, etc.

Returns a raw PSM DataFrame with one row per (precursor, scan) observation.

### 1.3 Per-File Prescore Filter (`iterative_prescore_filter!`)

**File**: `utils.jl`, lines 1774–1969

A single-round prescore cascade that eliminates junk precursors before the expensive feature pipeline runs. Steps:

#### Step 1: Add columns + quality filter
- `add_second_search_columns!()` — adds RT, charge, target/decoy status
- `get_isotopes_captured!()` — determines isotope traces in each scan's quad isolation window
- Filters: `precursor_fraction_transmitted >= min_fraction_transmitted` and `weight > 1e-6`

#### Step 2: Add ML features
- `ensure_ms1_stub_columns!()` — placeholder MS1 columns (sentinels, since real MS1 join hasn't happened)
- `add_features!()` — irt_error, irt_diff, TIC, amino acid composition, sequence_length, etc.
- `add_precursor_ms2_features!()` — precursor-in-MS2 feature
- `sanitize_prescore_features!()` — replaces Inf/NaN with 0.0

Uses `ITERATIVE_PRESCORE_FEATURES` (10 features):
```
fitted_manhattan_distance, max_matched_residual, gof, max_unmatched_residual,
poisson, irt_error, y_count, scribe, err_norm, spectral_contrast
```

#### Step 3: Two-pass probit on ALL PSMs
- **Pass 1**: Seeds training from scribe q-value <= 0.01 targets + all decoys. Trains probit, applies scores.
- **Pass 2**: Retrains on probit q <= 0.01 targets + all decoys. Applies updated scores. (Bootstrap from rough to refined.)

#### Step 4: Best PSM per precursor + LightGBM
- Selects one PSM per precursor by highest `probit_score`
- Trains LightGBM (200 iterations, max_depth=5, 31 leaves) on `LGBM_RECOVERY_FEATURES` (~25 features, superset of probit features)
- Computes q-values on LightGBM predictions

#### Step 5: Save per-file scores
Writes `{precursor_idx, lgbm_prob, target}` to `temp_data/prescore_scores/{file_name}.arrow` — consumed later by `aggregate_prescore_globally!`.

#### Step 6: Collect surviving pairs
- All precursors with `q_main <= prescore_qvalue_threshold` (default 0.10) pass
- Collects all `(precursor_idx, scan_idx)` tuples from the original PSMs for passing precursors

#### Step 7: Re-run deconvolution (`rerun_search_with_filter`)

**File**: `utils.jl`, lines 1667–1709

Rebuilds the fragment index mapping keeping only the surviving `(precursor_idx, scan_idx)` pairs. Re-runs `perform_second_pass_search()` with the filtered mapping.

**Key insight**: by removing non-passing precursor×scan combinations from the design matrix, surviving precursors face less competition per scan, yielding cleaner Huber deconvolution weights.

Returns fresh raw PSMs (no feature columns yet).

---

## Phase 2: Per-File Feature Computation (`process_search_results!`)

**File**: `SecondPassSearch.jl`, lines 560–631

Calls `compute_psm_features!()` then writes fold-split Arrow files.

### 2.1 `compute_psm_features!()` — The Core Pipeline

**File**: `SecondPassSearch.jl`, lines 418–558

Transforms raw PSMs into fully-featured rows for ScoringSearch's LightGBM (~50 features). Steps:

1. **`add_second_search_columns!()`** — RT, charge, target/decoy
2. **`get_isotopes_captured!()`** — isotope traces per scan's quad window
3. **Filter by fraction transmitted** — removes PSMs below `min_fraction_transmitted`
4. **Filter zero-weight** — removes PSMs with `weight <= 0.0`
5. **`train_and_apply_prescore!()`** — lightweight probit for apex scan selection (uses `PRESCORE_FEATURES`: fitted_spectral_contrast, err_norm, log2_intensity_explained, gof, scribe, weight). Marks `best_scan = true` for the probit-selected apex per precursor group.
6. **`init_summary_columns!()` + `get_summary_scores!()`** — computes summary stats in a window around the apex (±2 scans): max_gof, max_matched_ratio, max_fitted_manhattan_distance, max_fitted_spectral_contrast, max_scribe, y_ions_sum, max_y_ions, num_scans, smoothness
7. **Filter to apex** — `filter!(x -> x.best_scan, psms)` → one row per precursor×isotope_trace
8. **MS1 join** — `parseMs1Psms()` + leftjoin with `_ms1` suffix, sentinel fill, ms1_ms2_rt_diff in iRT space, ms1_features_missing flag
9. **Column ordering** — `get_expected_column_order()` enforces consistent Arrow schemas across files
10. **`add_features!()`** — irt_error, irt_diff, TIC, amino acid composition, sequence_length, etc.
11. **`add_precursor_ms2_features!()`** — precursor-in-MS2 feature
12. **`initialize_prob_group_features!()`** — placeholder columns for ScoringSearch: trace_prob, prec_prob, global_prob, cv_fold, etc.

### 2.2 Fold-Split Arrow Output

After `compute_psm_features!` returns:
- Splits the DataFrame by `cv_fold` (0 or 1)
- Writes `{base_name}_fold0.arrow` and `{base_name}_fold1.arrow` to `temp_data/second_pass_psms/`
- Stores the base path (without fold suffix) via `setSecondPassPsms!()`

---

## Phase 3: Cross-File Global Aggregation (`summarize_results!`)

**File**: `SecondPassSearch.jl`, lines 641–662

Called once after all files are processed.

### 3.1 `aggregate_prescore_globally!()`

**File**: `utils.jl`, lines 2273–2377

Aggregates per-file LightGBM prescore probabilities into a global score per precursor:

1. **Load per-file scores** — reads each `prescore_scores/{file_name}.arrow`, builds `Dictionary{UInt32, Vector{Float32}}` mapping each precursor to its per-file probabilities
2. **Log-odds aggregation** — for each precursor:
   - `top_n = floor(sqrt(n_files))` (e.g., 2 for 4 files)
   - Sorts the precursor's probabilities descending, takes top-N
   - Converts each to log-odds: `log(p / (1-p))`
   - Averages log-odds, converts back to probability
3. **Global q-values** — `get_qvalues!()` on aggregated probabilities
4. **Build passing set** — all precursors with `q <= GLOBAL_PRESCORE_QVALUE_THRESHOLD` (0.15)

### 3.2 Branch: Filter vs Re-deconvolve

Controlled by `params.global_rerun_deconvolution` (default `false`).

#### Option A: Simple Row Filter (`filter_arrow_files_to_passing!`)

**File**: `utils.jl`, lines 2385–2412

- Reads each fold-split Arrow file
- `filter!(row -> row.precursor_idx in passing_precs, tbl)`
- Rewrites in-place
- **Preserves existing weights/features** from the per-file prescore pass

#### Option B: Re-Deconvolution (`rerun_globally_filtered!`)

**File**: `utils.jl`, lines 2422–2489

For each file:
1. **Load spectra** via `getMSData(msdr, ms_file_idx)`
2. **`rerun_search_with_precursor_filter()`** — rebuilds fragment index keeping only globally-passing precursors (all scans), re-runs Huber deconvolution
3. **`compute_psm_features!()`** — full feature pipeline on fresh PSMs
4. **Overwrite** fold-split Arrow files

**Key difference from per-file prescore re-run**: `rerun_search_with_precursor_filter()` filters by precursor identity only (keeps ALL scans for passing precursors), vs `rerun_search_with_filter()` which filters by specific (precursor, scan) pairs.

| Aspect | Per-file (`rerun_search_with_filter`) | Global (`rerun_search_with_precursor_filter`) |
|--------|---------------------------------------|-----------------------------------------------|
| Filter key | `Set{Tuple{UInt32, UInt32}}` | `Set{UInt32}` |
| Granularity | Per-scan: precursor kept in some scans, excluded from others | Per-precursor: all scans kept or none |
| Used by | `iterative_prescore_filter!` | `rerun_globally_filtered!` |
| Purpose | Remove per-scan junk before feature computation | Remove globally-rejected precursors across all files |

---

## Downstream: How ScoringSearch Consumes the Arrow Files

**File**: `ScoringSearch/ScoringSearch.jl`

1. `get_valid_fold_file_paths()` collects all existing `_fold0.arrow` and `_fold1.arrow` paths
2. Two loading paths:
   - **OOM path**: `ArrowFilePSMContainer(file_paths)` streams files, samples for training, applies predictions file-by-file
   - **In-memory path**: `load_psms_for_lightgbm(folder)` reads all fold files into one DataFrame
3. LightGBM trains on ~50 features using cross-validation folds (fold 0 trains model applied to fold 1, and vice versa)
4. After scoring, fold files are merged back per-file: `{base_name}.arrow`
5. Merged files feed into protein inference, chromatogram integration, and MaxLFQ

---

## Complete Data Flow Diagram

```
Per file (process_file! + process_search_results!):
  ┌─────────────────────────────────────────────────────────────────┐
  │ load_fragment_index_matches                                     │
  │   → (scan_to_prec_idx, precursors_passed)                      │
  │                                                                 │
  │ perform_second_pass_search (Huber deconvolution)                │
  │   → raw PSMs                                                    │
  │                                                                 │
  │ iterative_prescore_filter!                                      │
  │   ├── add columns + quality filter                              │
  │   ├── add ML features (with MS1 stubs)                          │
  │   ├── two-pass probit on all PSMs                               │
  │   ├── best PSM per precursor (by probit_score)                  │
  │   ├── LightGBM on best PSMs → per-file scores                  │
  │   ├── save scores to prescore_scores/{file}.arrow               │
  │   ├── collect surviving (prec, scan) pairs at q ≤ 0.10         │
  │   └── rerun_search_with_filter → fresh raw PSMs                 │
  │                                                                 │
  │ compute_psm_features!                                           │
  │   ├── columns, isotopes, quality filter                         │
  │   ├── probit prescore for apex selection                        │
  │   ├── summary scores (±2 scan window)                           │
  │   ├── MS1 join                                                  │
  │   ├── add_features, add_precursor_ms2_features                  │
  │   └── initialize probability columns                            │
  │                                                                 │
  │ Write fold-split Arrow:                                         │
  │   {name}_fold0.arrow, {name}_fold1.arrow                        │
  └─────────────────────────────────────────────────────────────────┘

After all files (summarize_results!):
  ┌─────────────────────────────────────────────────────────────────┐
  │ aggregate_prescore_globally!                                    │
  │   ├── load prescore_scores/{file}.arrow for all files           │
  │   ├── log-odds combine top-sqrt(N) per precursor                │
  │   ├── global q-values                                           │
  │   └── threshold at q ≤ 0.15 → passing_precs Set{UInt32}       │
  │                                                                 │
  │ ┌─── if global_rerun_deconvolution ───┐                         │
  │ │ rerun_globally_filtered!            │                         │
  │ │   for each file:                    │                         │
  │ │     load spectra                    │                         │
  │ │     rerun_search_with_precursor_    │                         │
  │ │       filter (all scans for passing │                         │
  │ │       precursors)                   │                         │
  │ │     compute_psm_features!           │                         │
  │ │     overwrite fold Arrow files      │                         │
  │ └─────────────────────────────────────┘                         │
  │ ┌─── else ────────────────────────────┐                         │
  │ │ filter_arrow_files_to_passing!      │                         │
  │ │   read fold Arrows, drop rows not   │                         │
  │ │   in passing set, rewrite in-place  │                         │
  │ └─────────────────────────────────────┘                         │
  └─────────────────────────────────────────────────────────────────┘

Downstream (ScoringSearch):
  ┌─────────────────────────────────────────────────────────────────┐
  │ Load fold-split Arrow files                                     │
  │   → LightGBM training (fold 0 trains → scores fold 1, etc.)    │
  │   → merge folds per file                                        │
  │   → protein inference, chromatogram integration, MaxLFQ         │
  └─────────────────────────────────────────────────────────────────┘
```
