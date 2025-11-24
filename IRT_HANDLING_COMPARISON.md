# Comprehensive iRT Handling Comparison: `refine_library_irt` vs `develop`

**Date:** 2025-11-24
**Branch:** `refine_library_irt` (ahead of origin by 1 commit)
**Baseline:** `develop` branch
**Commit:** cd3fb5a0 - feat: Use library iRT for targeting, refined iRT for features only

## Executive Summary

The `refine_library_irt` branch introduces a **dual iRT system** that separates targeting from feature calculation:

- **Library iRT** (`irt_predicted`, `irt`): Used for **targeting** in chromatographic windows (unchanged)
- **Refined iRT** (`refined_irt_pred`, `refined_irt_obs`): Used for **feature calculation** in ML models (new)

This separation ensures robust targeting (library-based) while providing improved RT error features (refinement-based) for PSM discrimination in SecondPassSearch and ScoringSearch.

---

## Key Architectural Changes

### 1. New Infrastructure Files

#### `irt_refinement_utils.jl` (NEW in `refine_library_irt`)
**Location:** `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl`

**Purpose:** Implements per-file iRT refinement using amino acid composition

**Key Functions:**
```julia
# Creates 21-feature matrix (20 AA counts + library_irt) → error prediction
prepare_features_dataframe(sequences, library_irt, irt_errors)

# Trains linear regression model: error ~ irt_predicted + count_A + ... + count_V
fit_irt_refinement_model(sequences, irt_predicted, irt_observed; ms_file_idx, min_psms=20, train_fraction=0.67)

# Adds :refined_irt column to PSMs using streaming batch processing
add_refined_irt_column!(psms_path, refinement_model, search_context; batch_size=100_000)
```

**Algorithm:**
1. Filter sequences with sufficient PSMs (min 20)
2. Train/validation split (67/33)
3. Fit OLS: `error = intercept + β_irt * library_irt + Σ(β_aa * count_aa)`
4. Validate on held-out set
5. If MAE improves, retrain on full dataset; otherwise return `nothing`

**Does NOT exist in `develop` branch**

---

### 2. RT Conversion Models

#### `refine_library_irt` Branch
```julia
# RT → Library iRT (for targeting, unchanged)
library_irt = getRtIrtModel(search_context, ms_file_idx)(scan_rt)

# RT → Refined iRT (for features, NEW)
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(scan_rt)
```

#### `develop` Branch
```julia
# RT → iRT (single model for both targeting and features)
irt = getRtIrtModel(search_context, ms_file_idx)(scan_rt)
```

**Key Difference:**
- `refine_library_irt`: **Separate models** for targeting (`getRtIrtModel`) and features (`getRtToRefinedIrtModel`)
- `develop`: **Single model** (`getRtIrtModel`) for both purposes

---

### 3. Column Naming Convention

| Purpose                     | `develop` Branch    | `refine_library_irt` Branch       |
|-----------------------------|---------------------|-----------------------------------|
| **Library prediction**      | `:irt_predicted`    | `:irt_predicted` (unchanged)      |
| **Observed (for targeting)**| `:irt`              | `:irt` (library-based, unchanged) |
| **Refined prediction**      | N/A                 | `:refined_irt_pred` (NEW)         |
| **Observed (for features)** | `:irt` (same)       | `:refined_irt_obs` (NEW)          |

**Critical Insight:**
- In `develop`: `:irt` serves dual purpose (targeting + features)
- In `refine_library_irt`:
  - `:irt` = library-based observed iRT (for targeting)
  - `:refined_irt_obs` = observed iRT in refined space (for features)

---

## File-by-File Analysis

### FirstPassSearch: Initial RT Calibration

#### `FirstPassSearch/utils.jl`

**`develop` Branch (lines 95-110):**
```julia
# Line 98: Comment indicates library iRT usage
irt[i] = rt_irt(rt[i])  # No refinement comment

# Columns created
psms[!,:irt_predicted] = irt_pred  # Library iRT
psms[!,:rt] = rt
psms[!,:irt] = irt  # Observed iRT (library-based)
```

**`refine_library_irt` Branch (lines 96-122):**
```julia
# Line 98-99: Explicit comment about refinement timing
# Note: Using library iRT here since refinement happens after FirstPassSearch completes
irt[i] = rt_irt(rt[i])

# Columns created (same as develop)
psms[!,:irt_predicted] = irt_pred
psms[!,:rt] = rt
psms[!,:irt] = irt  # Still library-based at this stage
```

**Change:** Only documentation clarification; no functional changes in FirstPassSearch

**RT Model Fitting (lines 368-561 in `refine_library_irt`):**

The `map_retention_times!` function in `refine_library_irt` implements a 6-step pipeline:

1. **Step 1:** Fit RT → library_iRT spline (stored in `rt_irt_map`)
2. **Step 2:** Fit library_iRT → RT inverse spline (stored in `irt_rt_map`)
3. **Step 3:** Train iRT refinement model using `fit_irt_refinement_model`
4. **Step 4:** Fit RT → refined_iRT spline (stored in `rt_to_refined_irt_map`) *if refinement succeeds*
5. **Step 5:** Fit refined_iRT → RT inverse spline (stored in `refined_irt_to_rt_map`)
6. **Step 6:** Add `:refined_irt` column to PSMs file using `add_refined_irt_column!`

**Fallback Behavior:**
- If refinement fails (MAE doesn't improve), stores `IdentityModel()` for refined splines
- `:refined_irt` column still created (copies `:irt_predicted`)

---

### SecondPassSearch: Feature Extraction for ML

#### `SecondPassSearch/SecondPassSearch.jl`

**`develop` Branch (lines 490-553):**
```julia
# Line 495: Uses library iRT model for summary scores
getRtIrtModel(search_context, ms_file_idx)

# Line 550-553: MS1-MS2 RT diff using library iRT model
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                      Float32(-1),
                      abs.(rt_to_irt_model.(psms[!,:rt]) .- rt_to_irt_model.(psms[!,:rt_ms1]))))
```

**`refine_library_irt` Branch (lines 495-571):**
```julia
# Line 495: Uses REFINED iRT model for summary scores
getRtToRefinedIrtModel(search_context, ms_file_idx)

# Line 549-553: MS1-MS2 RT diff using REFINED iRT model
rt_to_refined_irt_model = getRtToRefinedIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                      Float32(-1),
                      abs.(rt_to_refined_irt_model.(psms[!,:rt]) .- rt_to_refined_irt_model.(psms[!,:rt_ms1]))))

# Lines 564-573: add_features! call passes BOTH models
add_features!(
    psms,
    search_context,
    getTICs(spectra),
    getMzArrays(spectra),
    ms_file_idx,
    getRtIrtModel(search_context, ms_file_idx),           # Library iRT for :irt_obs
    getRtToRefinedIrtModel(search_context, ms_file_idx),  # Refined iRT for :refined_irt_obs
    getPrecursorDict(search_context)
)
```

**Key Changes:**
1. **Summary Score Calculation** (line 495): Switches from `getRtIrtModel` → `getRtToRefinedIrtModel`
2. **MS1-MS2 RT Difference** (lines 549-553): Calculates in **refined iRT space** instead of library iRT space
3. **Dual Model Usage in add_features!** (lines 564-573): Passes **both** RT models:
   - `rt_to_irt_interp`: Library iRT model → calculates `:irt_obs` (line 911 in utils.jl)
   - `rt_to_refined_irt_interp`: Refined iRT model → calculates `:refined_irt_obs` (line 910)
4. **Feature Space:** All downstream RT-based features use refined iRT coordinates

**Targeting Unchanged:**
- Chromatographic window selection still uses **library iRT** (`getRtIrtModel` in `utils.jl`)
- Only **feature calculation** uses refined iRT

---

### ScoringSearch: ML Model Features

#### `model_config.jl`

**`develop` Branch (lines 44-105):**
```julia
const ADVANCED_FEATURE_SET = [
    # ... other features ...
    :irt_pred_qbin,    # Line 51: Quantile-binned library iRT
    :irt_error,        # Line 53: abs(irt_obs - irt_pred)
    :irt_diff,         # Line 54: Alternative error metric
    :ms1_ms2_rt_diff,  # Line 85: MS1-MS2 RT diff in iRT space
    # ... MS1 features ...
]
```

**`refine_library_irt` Branch (lines 44-105):**
```julia
const ADVANCED_FEATURE_SET = [
    # ... other features ...
    :refined_irt_pred_qbin,    # Line 51: Quantile-binned REFINED iRT
    :irt_error,                # Line 53: abs(refined_irt_obs - refined_irt_pred)
    :irt_diff,                 # Line 54: Alternative error metric
    :ms1_ms2_rt_diff,          # Line 85: MS1-MS2 RT diff in REFINED iRT space
    # ... MS1 features ...
]
```

**Feature Calculation (`percolatorSortOf.jl`):**

**`develop` Branch (lines 109-128):**
```julia
# Line 111: Uses library iRT for binning
psms[!, :irt_bin_idx] = getIrtBins(psms.irt_pred)

# Line 117: Error calculated using library iRT
irt_residual(psms, idx) = Float32(psms.irt_pred[idx] - psms.irt_obs[idx])
```

**`refine_library_irt` Branch (lines 111-128):**
```julia
# Line 111: Uses REFINED iRT for binning
psms[!, :irt_bin_idx] = getIrtBins(psms.refined_irt_pred)

# Line 117: Error calculated using REFINED iRT
irt_residual(psms, idx) = Float32(psms.refined_irt_pred[idx] - psms.refined_irt_obs[idx])
```

**Critical Changes:**
1. **iRT Binning:** Precursor pairing now uses **refined iRT bins** instead of library iRT bins
2. **Error Features:** `:irt_error` and `:irt_diff` computed using refined predictions/observations
3. **Quantile Binning:** `:refined_irt_pred_qbin` replaces `:irt_pred_qbin`

---

## Feature Flow Comparison

### `develop` Branch: Single iRT System

```
Library iRT ────────┐
                    │
                    ├──→ Targeting (RT windows)
                    │
                    └──→ Features (irt_error, irt_diff, ms1_ms2_rt_diff)
```

### `refine_library_irt` Branch: Dual iRT System

```
Library iRT ────────────────────────→ Targeting (RT windows, unchanged)

Observed RT ───┐
               ├──→ Refinement Model ──→ Refined iRT ──→ Features
Peptide Seq ───┘                         (irt_error, irt_diff,
                                          ms1_ms2_rt_diff)
```

---

## Correctness Verification

### ✅ Consistent Usage of `refined_irt_obs` and `refined_irt_pred`

**SecondPassSearch:**
- ✅ Uses `getRtToRefinedIrtModel` for feature extraction (line 495, 550, 570)
- ✅ MS1-MS2 RT difference calculated in refined iRT space (lines 549-553)

**ScoringSearch:**
- ✅ `ADVANCED_FEATURE_SET` uses `:refined_irt_pred_qbin` (model_config.jl:51)
- ✅ `irt_residual` function uses `refined_irt_pred` and `refined_irt_obs` (percolatorSortOf.jl:117)
- ✅ iRT binning uses `refined_irt_pred` (percolatorSortOf.jl:111)
- ✅ Quantile binning applied to `:refined_irt_pred` (score_psms.jl:103)

### ✅ Library iRT Preserved for Targeting

**SecondPassSearch/utils.jl:**
- ✅ `perform_second_pass_search` uses `getRtIrtModel` for targeting (unchanged)
- ✅ Chromatographic window selection based on library iRT predictions

**IntegrateChromatogramsSearch/utils.jl:**
- ✅ Uses `getRtIrtModel` for RT window extraction (unchanged)

**HuberTuningSearch/utils.jl:**
- ✅ Uses `getRtIrtModel` for parameter optimization (unchanged)

---

## Column Lifecycle Summary

| Stage            | Column Created           | Model Used                   | Purpose                     |
|------------------|--------------------------|------------------------------|-----------------------------|
| **FirstPassSearch** | `:irt_predicted`         | Library spline              | Library iRT values          |
|                  | `:irt`                   | Library RT→iRT spline       | Observed library iRT        |
|                  | `:refined_irt` (NEW)     | Refinement model (optional) | Refined library predictions |
| **SecondPassSearch** | `:ms1_ms2_rt_diff`       | `getRtToRefinedIrtModel` (NEW) | MS1-MS2 RT difference in refined space |
| **ScoringSearch** | `:refined_irt_obs` (derived) | RT→refined iRT spline   | Observed refined iRT        |
|                  | `:refined_irt_pred` (derived) | From `:refined_irt`    | Refined iRT predictions     |
|                  | `:irt_error`             | abs(refined_irt_obs - refined_irt_pred) | RT error feature    |
|                  | `:irt_bin_idx`           | getIrtBins(refined_irt_pred) | Precursor pairing bins     |

---

## Testing Recommendations

### 1. Feature Consistency Tests
```julia
@testset "Refined iRT Features" begin
    # Verify columns exist
    @test hasproperty(psms, :refined_irt_pred)
    @test hasproperty(psms, :refined_irt_obs)

    # Verify error calculation
    @test all(psms.irt_error .≈ abs.(psms.refined_irt_obs .- psms.refined_irt_pred))

    # Verify binning uses refined iRT
    bins = getIrtBins(psms.refined_irt_pred)
    @test bins == psms.irt_bin_idx
end
```

### 2. Targeting Consistency Tests
```julia
@testset "Library iRT Targeting" begin
    # Verify targeting unchanged
    library_irt = getRtIrtModel(context, file_idx)(scan_rt)

    # Verify chromatogram windows use library iRT
    @test window_center == library_irt
end
```

### 3. Refinement Model Tests
```julia
@testset "iRT Refinement" begin
    # Test model training
    model = fit_irt_refinement_model(sequences, irt_pred, irt_obs)

    # If model trained successfully
    if !isnothing(model)
        @test model.use_refinement == true
        @test model.mae_refined < model.mae_original
    end

    # Test fallback behavior
    psms = add_refined_irt_column!(psms_path, nothing, context)
    @test psms.refined_irt ≈ psms.irt_predicted  # Should copy when no refinement
end
```

### 4. Integration Tests
```julia
@testset "Full Pipeline" begin
    # Run SecondPassSearch
    execute_search(SecondPassSearch(), context, params)

    # Verify refined iRT used for features
    psms = load_psms(context)
    @test :refined_irt_pred in propertynames(psms)
    @test :ms1_ms2_rt_diff in propertynames(psms)

    # Run ScoringSearch
    execute_search(ScoringSearch(), context, params)

    # Verify ML features use refined iRT
    @test :refined_irt_pred_qbin in propertynames(psms)
end
```

---

## Migration Checklist

### ✅ Completed in `refine_library_irt` Branch

- [x] Create `irt_refinement_utils.jl` with refinement infrastructure
- [x] Add `getRtToRefinedIrtModel` accessor to SearchContext
- [x] Update SecondPassSearch to use refined iRT for features
- [x] Update ScoringSearch model_config.jl feature sets
- [x] Update percolatorSortOf.jl error calculations
- [x] Add `:refined_irt` column creation in FirstPassSearch
- [x] Preserve library iRT for targeting in all search methods

### ⚠️ Potential Future Enhancements

- [ ] Add diagnostic plots comparing library vs refined iRT predictions
- [ ] Report refinement model R² and MAE in search logs
- [ ] Implement cross-file refinement model sharing (if beneficial)
- [ ] Add refinement model serialization for reproducibility

---

## Conclusion

The `refine_library_irt` branch successfully implements a dual iRT system that:

1. **Preserves Targeting Robustness:** Library iRT continues to be used for all chromatographic window selection and targeting
2. **Improves Feature Quality:** Refined iRT provides more accurate RT error features for ML-based PSM discrimination
3. **Maintains Consistency:** All feature calculations consistently use `refined_irt_pred` and `refined_irt_obs`
4. **Handles Edge Cases:** Falls back gracefully when refinement doesn't improve MAE

The implementation is **correct and ready for merging** after validation testing confirms that:
- Targeting behavior is unchanged (library iRT)
- ML features correctly use refined iRT
- Fallback behavior works when refinement is disabled

**Key Success Metric:** PSM identification should improve or remain stable, with better RT error discrimination in regions where peptide composition affects retention behavior.
