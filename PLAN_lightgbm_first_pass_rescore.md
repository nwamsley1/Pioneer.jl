# Plan: LightGBM Rescoring of First-Pass Search PSMs

## Context

Currently, first-pass PSMs are scored with a probit model per-file using ~15 spectral/RT features. These per-file probabilities feed into `get_best_precursors_accross_runs()` which computes global probabilities via log-odds averaging. The probit model is a linear classifier — a LightGBM model trained on pooled cross-file data should achieve better target/decoy discrimination, yielding higher-quality global probabilities and ultimately better precursor selection for second-pass search.

**Approach**: After all files have been probit-scored, pool PSMs across files, sample up to 500K for LightGBM training with **cross-validation** (train on fold 0, predict fold 1; train on fold 1, predict fold 0), rescore ALL PSMs, then proceed to global probability calculation with the improved scores. Print diagnostic comparisons at multiple FDR thresholds to quantify the improvement.

## Files to Modify

| File | Purpose |
|------|---------|
| `src/.../FirstPassSearch/FirstPassSearch.jl` | Save feature columns to Arrow; call rescore step |
| `src/.../FirstPassSearch/lightgbm_rescore.jl` | **NEW** — LightGBM CV rescore function with diagnostics |
| `src/.../FirstPassSearch/getBestPrecursorsAccrossRuns.jl` | No changes (reads updated Arrow files as-is) |

All paths relative to `src/Routines/SearchDIA/SearchMethods/`.

## Changes

### 1. Preserve feature columns in Arrow files (`FirstPassSearch.jl`)

**`score_psms!()` (~line 359)**: Currently does `select!(psms, [:ms_file_idx, :score, ...])` which drops all feature columns. Change to also keep the probit feature columns and `:target`.

**`process_search_results!()` (~line 491-495)**: Currently writes only `[:ms_file_idx, :scan_idx, :precursor_idx, :rt, :irt_predicted, :q_value, :score, :prob, :scan_count, :PEP]`. Expand to also write the feature columns needed for LightGBM plus `:target`.

Feature columns to preserve (drop `:intercept` since LightGBM doesn't need it):
```julia
const FIRST_PASS_LGBM_FEATURES = [
    :spectral_contrast, :city_block, :entropy_score, :scribe,
    :charge2, :poisson, :irt_error, :missed_cleavage, :Mox,
    :TIC, :y_count, :err_norm, :spectrum_peak_count
]
```

Note: `:percent_theoretical_ignored` is conditionally removed when all zeros (line 323-325), so handle its absence gracefully.

### 2. New file: `lightgbm_rescore.jl`

Create `src/.../FirstPassSearch/lightgbm_rescore.jl` with a single public function:

```julia
const FIRST_PASS_LGBM_SAMPLE_SIZE = 500_000
const FIRST_PASS_LGBM_TRAIN_QVAL = 0.01f0

function rescore_first_pass_with_lightgbm!(
    search_context::SearchContext,
    fdr_scale_factor::Float32
)
```

**Algorithm**:

1. **Read & pool**: Collect paths from `getFirstPassPsms(getMSData(search_context))`, read all per-file Arrow files into a single DataFrame

2. **Add cv_fold column** from library:
   ```julia
   precursors = getPrecursors(getSpecLib(search_context))
   psms[!, :cv_fold] = [getCvFold(precursors, pid) for pid in psms.precursor_idx]
   ```
   (`getCvFold` defined in `LibraryIon.jl:712`)

3. **Sample**: Random sample of `min(FIRST_PASS_LGBM_SAMPLE_SIZE, nrow(all_psms))` rows, maintaining both folds in the sample

4. **Define features**: `FIRST_PASS_LGBM_FEATURES` filtered to columns actually present (handles missing `percent_theoretical_ignored`)

5. **Cross-validated training & prediction** (2-fold CV):
   ```
   For each fold in unique(cv_fold):
     train_data = sample rows where cv_fold != fold AND ((target & q_value ≤ 0.01) | !target)
     Train LightGBM on train_data
     Predict on ALL psms where cv_fold == fold (not just the sample)
   ```
   This ensures no PSM is ever scored by a model that trained on it.

6. **Guard**: if < 100 targets or < 100 decoys in either fold's training set, skip rescoring (log warning, return early)

7. **Diagnostic comparison**: Before overwriting scores, compute global probabilities BOTH ways and log comparison (see Section 3 below)

8. **Update scores**: Replace `prob` column with LightGBM predictions, recompute `q_value` via `get_qvalues!()` and `PEP` via `get_PEP!()`

9. **Write back**: Write rescored PSMs back to per-file Arrow files (grouped by `ms_file_idx`)

**Existing utilities to reuse**:
- `build_lightgbm_classifier()` — `src/utils/ML/lightgbm_utils.jl:53`
- `fit_lightgbm_model()` — `src/utils/ML/lightgbm_utils.jl:116`
- `lightgbm_predict()` / `predict()` — `src/utils/ML/lightgbm_utils.jl:134,161`
- `feature_matrix()` — `src/utils/ML/lightgbm_utils.jl:18`
- `get_qvalues!()` — `src/utils/ML/fdrUtilities.jl`
- `get_PEP!()` — `src/utils/ML/fdrUtilities.jl`
- `getCvFold()` — `src/structs/LibraryIon.jl:712`
- `_logodds_combine()` — `src/.../FirstPassSearch/getBestPrecursorsAccrossRuns.jl`

**LightGBM hyperparams** (conservative, similar to SimpleLightGBM in `model_config.jl`):
```julia
build_lightgbm_classifier(
    num_iterations = 100,
    max_depth = 6,
    num_leaves = 31,
    learning_rate = 0.05,
    feature_fraction = 0.5,
    bagging_fraction = 0.5,
    bagging_freq = 1,
    min_data_in_leaf = 500,
    min_gain_to_split = 1.0
)
```

### 3. Diagnostic comparison (inside `rescore_first_pass_with_lightgbm!`)

After LightGBM predictions are computed but BEFORE overwriting the Arrow files, compute global probabilities using both the original probit scores and the new LightGBM scores, then log a side-by-side comparison.

**Implementation**: Extract a helper function that mirrors the global probability logic from `get_best_precursors_accross_runs()` (lines 298-326):

```julia
function compute_global_fdr_diagnostics(
    psms::DataFrame,             # pooled PSMs with :precursor_idx, :prob, :target
    prec_is_decoy::AbstractVector{Bool},
    n_valid_files::Int,
    fdr_scale_factor::Float32,
    label::String                # "Probit" or "LightGBM"
)
```

This function:
1. Groups PSMs by precursor, picks best `prob` per precursor per file
2. Combines per-file probs via `_logodds_combine()` (same as `getBestPrecursorsAccrossRuns.jl:307`)
3. Computes global PEP via `get_PEP!()` and global q-values via `get_qvalues!()`
4. Logs at multiple FDR thresholds:

```
[ Info: === First-Pass Rescoring Diagnostics (Probit) ===
[ Info: Global PEP (log-odds): 898374 unique precursors (480858 targets, 417516 decoys) from 23 files
[ Info:   1% FDR: 37324 | 5% FDR: 50341 | 10% FDR: 58920 | 15% FDR: 63100 | 20% FDR: 66800
[ Info: === First-Pass Rescoring Diagnostics (LightGBM) ===
[ Info: Global PEP (log-odds): 898374 unique precursors (480858 targets, 417516 decoys) from 23 files
[ Info:   1% FDR: 41200 | 5% FDR: 55100 | 10% FDR: 63400 | 15% FDR: 67800 | 20% FDR: 71500
```

Call this function twice:
1. First with the original probit `prob` column (before overwriting)
2. Then with the LightGBM `prob` column (after prediction)

Then proceed with the LightGBM scores for the actual pipeline.

Note: `get_best_precursors_accross_runs()` will ALSO print its own diagnostics when it runs later — those will reflect the LightGBM scores since the Arrow files will have been updated.

### 4. Call rescore in `summarize_results!()` (`FirstPassSearch.jl`, ~line 520)

Insert the LightGBM rescore call **after** `map_retention_times!()` (line 564) and **before** `get_best_precursors_accross_runs!()` (line 566):

```julia
# Map retention times
map_retention_times!(search_context, params, all_psms_paths)

# LightGBM rescoring of first-pass PSMs (cross-validated)
fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
rescore_first_pass_with_lightgbm!(search_context, fdr_scale_factor)

# Process precursors (now uses LightGBM-rescored probabilities)
precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)
```

### 5. Include the new file

No changes needed. `importScripts.jl:234-244` uses `walkdir` to auto-include all `.jl` files in `SearchMethods/` subdirectories (except ParameterTuningSearch and ScoringSearch which are loaded explicitly). The new file will be picked up automatically.

## Data Flow

```
Per-file loop:
  search → add_psm_columns! → score_psms! (probit) → select_best_psms!
  → write Arrow (now includes feature columns + target)

summarize_results!:
  map_retention_times!
  → rescore_first_pass_with_lightgbm!   ← NEW
    ├─ read all Arrow files → pool into DataFrame
    ├─ add cv_fold from library
    ├─ sample up to 500K PSMs for training
    ├─ CV loop: train fold 0 → predict fold 1, train fold 1 → predict fold 0
    ├─ DIAGNOSTIC: compute global probs with probit scores → log at 1/5/10/15/20% FDR
    ├─ DIAGNOSTIC: compute global probs with LightGBM scores → log at 1/5/10/15/20% FDR
    ├─ update prob/q_value/PEP on all PSMs with LightGBM scores
    └─ write back to per-file Arrow files
  → get_best_precursors_accross_runs!    (reads updated Arrow files, prints its own diagnostics)
  → setPrecursorDict!
```

## Verification

1. `SearchDIA("./data/ecoli_test/ecoli_test_params.json")` — check:
   - Log shows Probit vs LightGBM diagnostic comparison at 5 FDR thresholds
   - No errors; pipeline completes end-to-end
   - Feature columns present in temp Arrow files
2. Create test JSON for Olsen Exploris regression test
3. Compare Probit vs LightGBM target counts at each FDR threshold in logs
