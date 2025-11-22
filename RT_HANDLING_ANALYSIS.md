# Comprehensive Retention Time Handling Analysis
## refine_library_irt Branch vs develop Branch

**Date**: 2025-01-21 (Analysis) | 2025-01-21 (Bug Fix)
**Branch**: `refine_library_irt`
**Analysis Goal**: Identify why iRT refinement (50% MAE reduction) is not producing more identifications

---

## Executive Summary

### Critical Bug Fixed ✅

**Location**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl` lines 108 & 111

**Status**: **FIXED** in commit `37fd89e0` (2025-01-21)

**Impact**: RT indices were being built using **incorrect refined iRT values** (most recent PSM instead of best PSM), preventing identification improvements despite successful iRT refinement.

**Severity**: HIGH - This bug existed in both `refine_library_irt` and `develop` branches, affecting all Pioneer searches.

**Fix**: Two one-line changes implemented (details in Section 1.3)

### New Critical Bug Discovered 🔴

**Location**: `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl` (makeRTIndices function)

**Status**: **IDENTIFIED** (2025-01-21) - awaiting fix

**Discovery**: After fixing the first bug, testing revealed that iRT refinement doesn't improve IDs despite 50% MAE reduction. Investigation revealed a second critical bug in RT index imputation logic.

**Problem**: When imputing refined iRT for precursors not observed in a file, the code uses `best_refined_irt` from a **different** file (File A) that has a **different** refinement model, instead of calculating file-specific refined iRT using the **current** file's (File B) refinement model.

**Impact**:
- Systematic cross-file model mismatch in RT indices
- Prevents refinement from improving identifications
- Can make refined iRT **worse** than library iRT (explaining why refinement hurts IDs)

**Fix**: Use file-specific refinement model to calculate refined iRT for imputed precursors (details in Section 7)

**Priority**: **CRITICAL** - This is likely the primary reason iRT refinement isn't helping

---

## Table of Contents

1. [Critical Bug: getBestPrecursorsAccrossRuns.jl](#1-critical-bug-getbestprecursorsaccrossrunsjl)
2. [What Changed from Develop Branch](#2-what-changed-from-develop-branch)
3. [Complete RT Handling Flow](#3-complete-rt-handling-flow)
4. [Additional Potential Issues](#4-additional-potential-issues)
5. [Implementation Correctness Review](#5-implementation-correctness-review)
6. [Recommended Fix Plan](#6-recommended-fix-plan)
7. [Critical Imputation Bug Discovery](#7-critical-imputation-bug-discovery-2025-01-21-)

---

## 1. Critical Bug: getBestPrecursorsAccrossRuns.jl

### 1.1 Bug Description ✅ FIXED

**Status**: Fixed in commit `37fd89e0` on 2025-01-21

The `readPSMs!()` function in `getBestPrecursorsAccrossRuns.jl` (lines 74-129) processes PSMs from multiple files to identify the best match for each precursor across all runs. The bug occurred when updating an existing precursor entry.

**Correct behavior** (lines 97-102):
```julia
# Update best match if current is better
if (best_prob < prob)
    best_prob = prob                 # Updates local variable ✅
    best_refined_irt = refined_irt   # Updates local variable ✅
    best_scan_idx = scan_idx         # Updates local variable ✅
    best_ms_file_idx = best_ms_file_idx
end
```

**Previous incorrect behavior** (lines 107-115, BEFORE FIX):
```julia
# Store updated values back to dictionary
prec_to_best_prob[precursor_idx] = (
    best_prob = prob,                      # ❌ WRONG! Uses current PSM, not best
    best_ms_file_idx = best_ms_file_idx,   # ✅ Correct
    best_scan_idx = best_scan_idx,         # ✅ Correct
    best_refined_irt = refined_irt,        # ❌ WRONG! Uses current PSM, not best
    mean_refined_irt = mean_refined_irt,   # ✅ Correct
    var_refined_irt = var_refined_irt,     # ✅ Correct
    n = n,                                 # ✅ Correct
    mz = mz)                               # ✅ Correct
```

**What was wrong**:
- Lines 97-102 updated local variables `best_prob` and `best_refined_irt` when a better PSM was found
- Lines 107-115 should have stored these **updated local variables** back to the dictionary
- Instead, they stored the **current loop iteration values** (`prob` and `refined_irt`), discarding the best match logic

### 1.2 Impact Analysis

#### Immediate Impact
1. **RT indices built with wrong refined iRT values**
   - Uses refined iRT from most recently processed PSM file
   - Not the refined iRT from highest probability match
   - Different for each precursor depending on file processing order

2. **Inconsistent cross-file behavior**
   - Same precursor gets different "best" refined iRT depending on which file was processed last
   - Cross-run variance calculations are corrupted

3. **Explains lack of ID improvement**
   - iRT refinement model works correctly (50% MAE reduction demonstrated)
   - But RT indices use suboptimal refined iRT values for precursor lookup
   - RT windows centered on wrong refined iRT → miss valid matches

#### Evidence This Is The Problem
1. **Refinement shows 50% MAE improvement** - model is working
2. **No increase in identifications** - refined values not being used correctly
3. **Bug affects both prob and refined_irt** - exactly the values used for RT indexing
4. **Bug present in develop branch too** - long-standing issue (uses `best_irt = irt` on line 104)

### 1.3 Fix Implemented ✅

**Status**: Fixed in commit `37fd89e0` on 2025-01-21

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`

**Line 108** - Fixed best_prob:
```julia
# BEFORE (WRONG):
best_prob = prob,

# AFTER (CORRECT):
best_prob = best_prob,
```

**Line 111** - Fixed best_refined_irt:
```julia
# BEFORE (WRONG):
best_refined_irt = refined_irt,

# AFTER (CORRECT):
best_refined_irt = best_refined_irt,
```

**Commit details**:
- Commit: `37fd89e0`
- Message: "fix: Store best PSM values instead of current values in getBestPrecursorsAccrossRuns"
- Changes: 2 lines (1 file)
- Branch: `refine_library_irt`

### 1.4 Validation Testing (PENDING)

To validate the fix works correctly, verify:
1. `best_refined_irt` matches the refined iRT from the highest probability PSM
2. `best_prob` matches the highest probability across all files
3. Cross-run variance (`var_refined_irt`) is calculated correctly
4. Identification rates improve in SecondPassSearch vs FirstPassSearch

---

## 2. What Changed from Develop Branch

### 2.1 New Files

**`src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl`**
- `fit_irt_refinement_model()`: Trains linear regression on AA composition + library iRT
- `add_refined_irt_column!()`: Applies refinement model to PSM files
- 284 lines total

### 2.2 Modified Files

#### Core RT Conversion
**`src/structs/RetentionTimeConversionModel.jl`**
- Added `IrtRefinementModel` struct (lines 88-121)
- Callable interface: `refined_irt = model(sequence, library_irt)`
- Algorithm: Predicts error, returns `library_irt - predicted_error`

#### FirstPassSearch Changes
**`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`**
- Expanded `map_retention_times!()` from 2 steps to 6 steps (lines 329-561)
  - **Steps 1-2**: Library iRT alignment (unchanged)
  - **Steps 3-6**: NEW - Refinement training + refined iRT splines
- Modified `create_rt_indices!()` to use `best_refined_irt` (line 658)
- Modified `get_irt_errs()` to use `var_refined_irt` for tolerance (lines 744-753)

**`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`**
- Renamed all iRT fields: `best_irt` → `best_refined_irt`, `mean_irt` → `mean_refined_irt`, `var_irt` → `var_refined_irt`
- Changed function parameter: `rt_irt` → `rt_to_refined_irt`
- Updated all references throughout (60+ locations)

#### RT Index Building
**`src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl`**
- Changed `makeRTIndices()` parameter: `rt_to_irt_splines` → `rt_to_refined_irt_splines`
- Updated imputation logic to use refined iRT values
- Updated comments to clarify refined iRT usage

#### Post-FirstPassSearch Searches
**`src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`**
- Lines 181-183, 386-388: Use `getRtToRefinedIrtModel()` for RT window calculations
- Lines 834-971: Feature calculation using refined iRT
  - `:refined_irt_obs`: Observed refined iRT from scan RT
  - `:refined_irt_pred`: Predicted refined iRT using refinement model + library iRT
  - `:irt_error`: abs(refined_irt_obs - refined_irt_pred)
  - `:irt_diff`: abs(refined_irt_obs - best_refined_irt)
- Lines 949-950: DataFrame columns renamed to `:refined_irt_pred`/`:refined_irt_obs`

**`src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl`**
- Lines 225-227: Use `getRtToRefinedIrtModel()` for RT window calculations

**`src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`**
- Lines 254-256, 491-493: Use `getRtToRefinedIrtModel()` for RT window calculations

#### SearchContext Management
**`src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`**
- Added getter functions (lines 520-550):
  - `getRtToRefinedIrtModel(s, index)` - Returns RT → refined_iRT model with fallback to library iRT
  - `getRefinedIrtToRtModel(s, index)` - Returns refined_iRT → RT model with fallback
  - `getIrtRefinementModel(s, index)` - Returns refinement model for file
- Note: `getPredIrt()` unchanged - correctly returns library iRT from internal field

#### Scoring Pipeline
**`src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`**
- Lines 83, 103: Changed features_to_bin from `:irt_pred` → `:refined_irt_pred`

**`src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl`**
- Line 51: Comment updated to `:refined_irt_pred_qbin`
- Line 52: Commented feature updated to `:refined_irt_pred`
- Lines 353-354: Documentation example updated

**`src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`**
- Line 512: Column list updated from `:irt_obs` → `:refined_irt_obs`

**`src/Routines/SearchDIA/CommonSearchUtils/normalizeQuant.jl`**
- Lines 28-33, 42, 98: All `:irt_obs` → `:refined_irt_obs` references (7 locations)

**`src/utils/ML/percolatorSortOf.jl`**
- Line 111: `psms.irt_pred` → `psms.refined_irt_pred`
- Line 117: Both `psms.irt_pred` and `psms.irt_obs` → refined versions
- Line 124: `psms.irt_pred` → `psms.refined_irt_pred`

### 2.3 New SearchContext Fields

Three new dictionaries added to SearchContext:
1. `rt_to_refined_irt_map::Dict{Int, RtConversionModel}` - RT → refined_iRT splines per file
2. `refined_irt_to_rt_map::Dict{Int, RtConversionModel}` - refined_iRT → RT inverse splines per file
3. `irt_refinement_models::Dict{Int, Union{IrtRefinementModel, Nothing}}` - Refinement models per file

### 2.4 Summary of Changes

**Total files modified**: 13
**Total lines changed**: ~350 additions, ~150 deletions
**New functionality**: File-specific iRT refinement using amino acid composition
**Architecture**: Maintains library iRT alongside refined iRT, never overwrites library values

---

## 3. Complete RT Handling Flow

### 3.1 FirstPassSearch: RT Alignment Pipeline

**Function**: `map_retention_times!()` in `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`
**Lines**: 329-561

#### Step 1: Fit RT → library_iRT spline (lines 372-390)
```julia
# Uses high-confidence PSMs (prob > threshold)
rt_to_library_irt, valid_rt_library, valid_library_irt, _ = fit_irt_model(
    psms_df;
    lambda_penalty = 0.1,
    ransac_threshold = 1000,
    min_psms = 10,
    spline_degree = 3,
    max_knots = 7,
    outlier_threshold = 5.0
)
# Stores in rt_irt_map[ms_file_idx]
```

#### Step 2: Fit library_iRT → RT inverse spline (lines 392-408)
```julia
# Uses validated RT and library_iRT from Step 1
library_irt_to_rt, _, _, _ = fit_irt_model(
    inverse_df;
    # same parameters...
)
# Stores in irt_rt_map[ms_file_idx]
```

#### Step 3: Train iRT refinement model (lines 410-430)
```julia
# Calculate observed_irt for best PSMs using library spline
observed_irt = [rt_to_library_irt(rt) for rt in psms[:rt][best_hits]]

# Get sequences for best PSMs
best_sequences = [sequences[idx] for idx in psms[:precursor_idx][best_hits]]

# Train refinement model
refinement_model = fit_irt_refinement_model(
    best_sequences,
    psms[:irt_predicted][best_hits],  # Library iRT
    observed_irt,                      # Observed iRT from RT
    ms_file_idx = ms_file_idx,
    min_psms = 20,
    train_fraction = 0.67
)
# Stores in irt_refinement_models[ms_file_idx]
```

**Refinement Model Algorithm** (irt_refinement_utils.jl lines 106-220):
- **Features**: 20 amino acid counts + library_irt (21 features total)
- **Response variable**: error = library_irt - observed_irt
- **Algorithm**: Ordinary least squares linear regression
- **Train/Val split**: 67% train, 33% validation
- **Validation**: Only enables refinement if validation MAE < baseline MAE
- **Returns**: IrtRefinementModel with callable interface

#### Step 4: Fit RT → refined_iRT spline (lines 432-458)
```julia
# Only if refinement succeeds
if !isnothing(refinement_model) && refinement_model.use_refinement
    # Apply refinement model to get refined_irt for best PSMs
    refined_irt_values = [refinement_model(seq, lib_irt)
                         for (seq, lib_irt) in zip(best_sequences,
                                                    psms[:irt_predicted][best_hits])]

    # Fit RT → refined_iRT spline
    rt_to_refined_irt, valid_rt_refined, valid_refined_irt, _ = fit_irt_model(
        refined_psms_df;
        # same parameters...
    )
    # Stores in rt_to_refined_irt_map[ms_file_idx] (NEW FIELD)
end
```

#### Step 5: Fit refined_iRT → RT inverse spline (lines 460-476)
```julia
refined_irt_to_rt, _, _, _ = fit_irt_model(
    refined_irt_to_rt_df;
    # same parameters...
)
# Stores in refined_irt_to_rt_map[ms_file_idx] (NEW FIELD)
```

#### Step 6: Add :refined_irt column to PSMs file (lines 484-486)
```julia
add_refined_irt_column!(psms_path, refinement_model, search_context)
# Applies refinement model to all PSMs
# Fallback: copies :irt_predicted if no refinement
```

### 3.2 RT Index Building

**Function**: `makeRTIndices()` in `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl`
**Lines**: 70-113

**Called from**: `create_rt_indices!()` in FirstPassSearch/utils.jl line 689

```julia
# For each precursor in the library
for (prec_id, irt_mz) in pairs(prec_to_irt)
    if haskey(prec_set, prec_id)
        _irt_, prob = prec_set[prec_id]
        # High-confidence PSMs: use empirical refined iRT
        if (prob >= min_prob)
            irts[i] = rt_to_refined_irt(scan_rt)  # Empirical refined iRT ✅
            continue
        end
    end
    # Low-confidence or missing: impute from best across runs
    irts[i] = irt_mz.irt  # From prec_to_irt[prec_id].irt (best_refined_irt) ❌ BUG!
end

# Build RT index using refined iRT values
rt_index = buildRtIndex(irts, mzs, prec_ids; bin_width=0.1f0)
```

**Bug impact**: `prec_to_irt[prec_id].irt` contains `best_refined_irt` from getBestPrecursorsAccrossRuns, which has wrong values due to the bug!

### 3.3 iRT Tolerance Calculation

**Function**: `get_irt_errs()` in `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`
**Lines**: 721-763

```julia
# Component 1: Peak width in refined iRT space
for (file_idx, fwhm_stats) in pairs(fwhms)
    median_fwhm = fwhm_stats.median_fwhm
    mad_fwhm = fwhm_stats.mad_fwhm
    peak_width = median_fwhm + fwhm_nstd * mad_fwhm
end

# Component 2: Cross-run variance in refined iRT space
variance_ = collect(skipmissing(
    map(x -> (x[:n] > 2) ? sqrt(x[:var_refined_irt]/(x[:n] - 1)) : missing,
        prec_to_irt)
))
irt_std = median(variance_)

# Combined tolerance
irt_errs[file_idx] = peak_width + irt_nstd * irt_std
```

**Bug impact**: `var_refined_irt` is calculated from wrong refined iRT values, leading to incorrect cross-run variance estimates!

### 3.4 Post-FirstPassSearch: RT Window Calculations

**Used in**: SecondPassSearch, HuberTuningSearch, IntegrateChromatogramsSearch

**Pattern** (consistent across all methods):
```julia
# 1. Convert scan RT to refined iRT space
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(scan_rt)

# 2. Get tolerance for this file
irt_tol = getIrtErrors(search_context)[ms_file_idx]

# 3. Calculate RT window bounds in refined iRT space
irt_start = searchsortedfirst(rt_index.rt_bins, refined_irt - irt_tol, ...)
irt_stop = searchsortedlast(rt_index.rt_bins, refined_irt + irt_tol, ...)

# 4. Iterate through RT bins in window
for rt_bin_idx in irt_start:irt_stop
    # Process precursors in this refined iRT bin
end
```

**Locations**:
- SecondPassSearch/utils.jl: Lines 181-183 (MS2 path), Lines 386-388 (MS1 path)
- HuberTuningSearch/utils.jl: Line 226
- IntegrateChromatogramsSearch/utils.jl: Lines 254, 491

### 3.5 iRT-Based Feature Calculation

**Function**: `add_features!()` in `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`
**Lines**: 829-971

```julia
# Get refinement model for this file (before parallel loop)
refinement_model = getIrtRefinementModel(search_context, ms_file_idx)

# Inside parallel loop for each PSM
for i in 1:N
    # Observed refined iRT from scan RT
    refined_irt_obs[i] = rt_to_refined_irt_interp(rt[i])

    # Predicted refined iRT using refinement model + library iRT
    library_irt = getPredIrt(search_context, prec_idx)
    refined_irt_pred[i] = if !isnothing(refinement_model) && refinement_model.use_refinement
        refinement_model(precursor_sequence[prec_idx], library_irt)
    else
        library_irt  # Fallback to library iRT
    end

    # Feature: Difference from best refined iRT across runs
    irt_diff[i] = abs(refined_irt_obs[i] - prec_id_to_irt[prec_idx].best_refined_irt)

    # Feature: Prediction error
    irt_error[i] = abs(refined_irt_obs[i] - refined_irt_pred[i])

    # MS1-level refined iRT difference
    if !ms1_missing[i]
        ms1_refined_irt_obs = rt_to_refined_irt_interp(ms1_rt[i])
        ms1_irt_diff[i] = abs(ms1_refined_irt_obs - refined_irt_pred[i])
    end
end

# Store as DataFrame columns
psms[!,:refined_irt_obs] = refined_irt_obs
psms[!,:refined_irt_pred] = refined_irt_pred
psms[!,:irt_error] = irt_error
psms[!,:irt_diff] = irt_diff
```

**Features using refined iRT**:
- `:refined_irt_obs` - Observed refined iRT from scan RT
- `:refined_irt_pred` - Predicted refined iRT using model
- `:irt_error` - abs(refined_irt_obs - refined_irt_pred)
- `:irt_diff` - abs(refined_irt_obs - best_refined_irt) ❌ Uses buggy best_refined_irt!
- `:ms1_irt_diff` - MS1-level refined iRT difference

### 3.6 Scoring Model Feature Usage

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl`
**Lines**: 44-100

**ADVANCED_FEATURE_SET includes**:
- `:refined_irt_pred_qbin` (line 51) - Quantile-binned predicted refined iRT
- `:irt_error` (line 53) - Refined iRT prediction error
- `:irt_diff` (line 54) - Difference from best refined iRT ❌ Uses buggy value!
- `:ms1_ms2_rt_diff` (line 86) - MS1/MS2 RT difference (related)

**Quantile binning** (score_psms.jl lines 83, 103):
```julia
features_to_bin = [:prec_mz, :refined_irt_pred, :weight, :tic]
add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)
# Creates :refined_irt_pred_qbin column
```

---

## 4. Additional Potential Issues

### 4.1 Silent Fallback to Library iRT

**Location**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl` lines 520-529

```julia
function getRtToRefinedIrtModel(s::SearchContext, index::I) where {I<:Integer}
    if haskey(s.rt_to_refined_irt_map, index)
        return s.rt_to_refined_irt_map[index]
    elseif haskey(s.rt_irt_map, index)
        @debug "Refined iRT model not found for file $index, falling back to library iRT model"
        return s.rt_irt_map[index]
    else
        return IdentityModel()
    end
end
```

**Issue**: Uses `@debug` level for fallback warning
- Post-FirstPassSearch searches expect refined iRT models
- Fallback to library iRT would cause incorrect window calculations
- Should be `@warn` or `@user_warn` level for visibility

**Recommendation**: Change to `@user_warn` level

### 4.2 Potentially Narrower RT Windows

**Hypothesis**: Refinement reduces variance too much, making RT windows too narrow

**Evidence**:
- Refinement reduces MAE by ~50% (reported in logs)
- Lower error → lower cross-run variance
- Formula: `irt_tol = fwhm + nstd * sqrt(var_refined_irt)`
- Lower variance → narrower windows

**Potential impact**:
- Narrower windows might miss valid matches with slightly different RT behavior
- Could explain why more accurate iRT doesn't lead to more IDs

**Mitigation**:
- MAE validation check ensures refinement actually helps on validation set
- But validation set is same file, not cross-file behavior

**Recommendations**:
1. Log RT window widths before/after refinement for comparison
2. Track identification rates: FirstPass vs SecondPass vs HuberTuning
3. Check if RT window misses are increasing in SecondPassSearch
4. Consider adaptive tolerance scaling based on MAE improvement

### 4.3 No Quality Metrics Logged

**Missing observability**:
1. RT window width comparison (library iRT vs refined iRT space)
2. Identification rates per search method
3. Refinement success rate per file
4. Cross-run variance before/after refinement
5. FWHM distributions in library vs refined iRT space

**Recommendation**: Add comprehensive logging in FirstPassSearch after RT alignment:
```julia
@user_info "RT alignment statistics:"
@user_info "  Library iRT MAE: $library_mae"
@user_info "  Refined iRT MAE: $refined_mae ($(100*(1-refined_mae/library_mae))% improvement)"
@user_info "  RT window width (library): $library_window"
@user_info "  RT window width (refined): $refined_window"
@user_info "  Cross-run variance (library): $library_variance"
@user_info "  Cross-run variance (refined): $refined_variance"
```

### 4.4 Overfitting Risk

**Current approach**:
- 67/33 train/val split on same file
- Validates MAE improvement on held-out PSMs from **same file**
- Does not validate cross-file generalization

**Potential issue**:
- Model might overfit to specific file's RT behavior
- Could perform worse on other files in the experiment

**Current mitigation**:
- Each file gets its own refinement model (file-specific)
- No cross-file model sharing, so overfitting to one file doesn't affect others

**Recommendation**: Monitor per-file refinement success rates

---

## 5. Implementation Correctness Review

### 5.1 ✅ Correctly Implemented

1. **iRT refinement model architecture**
   - Sound algorithm: Linear regression on AA composition + library iRT
   - 21 features: 20 AA counts + library_irt
   - Predicts error = library_irt - observed_irt
   - Returns refined_irt = library_irt - predicted_error
   - Validation ensures improvement before enabling

2. **All post-FirstPassSearch searches**
   - Correctly use `getRtToRefinedIrtModel()` for RT → refined_iRT conversion
   - Consistent pattern across SecondPassSearch, HuberTuningSearch, IntegrateChromatogramsSearch
   - No searches still using library iRT models after FirstPassSearch

3. **RT window calculations**
   - Consistently use refined iRT space throughout pipeline
   - Tolerance calculated in refined iRT space
   - Window bounds calculated in refined iRT space
   - RT index built using refined iRT values

4. **Scoring features**
   - All iRT-based features use refined iRT values
   - `:refined_irt_pred` computed correctly using refinement model
   - `:refined_irt_obs` computed correctly from RT conversion
   - `:irt_error` and `:irt_diff` use correct refined iRT values
   - Column naming is consistent throughout pipeline

5. **Column naming consistency**
   - All DataFrame columns consistently use `refined_irt_pred` and `refined_irt_obs`
   - Updated across 5 files in commit e5608e3a
   - No references to old `:irt_pred` or `:irt_obs` column names remain

6. **Fallback mechanism**
   - `getRtToRefinedIrtModel()` correctly falls back to library iRT if refined unavailable
   - `IdentityModel()` returned if no models available
   - Ensures pipeline doesn't crash if refinement fails

7. **Library iRT preservation**
   - `getPredIrt()` correctly returns library iRT from internal `s.irt_obs` field
   - Never overwrites library iRT values
   - Maintains separate `rt_irt_map` (library) and `rt_to_refined_irt_map` (refined)

### 5.2 ✅ Critical Bug (FIXED)

**getBestPrecursorsAccrossRuns.jl** lines 108, 111
- **Status**: Fixed in commit `37fd89e0` on 2025-01-21
- **Issue**: Was storing current PSM values instead of best PSM values
- **Affected**: Both `best_prob` and `best_refined_irt`
- **Impact**: RT indices were being built with wrong refined iRT values

### 5.3 ⚠️ Areas Needing Improvement

1. **Fallback warning level** - Should be `@user_warn` instead of `@debug`
2. **Quality metrics logging** - Need RT window width comparisons, success rates
3. **Overfitting monitoring** - No cross-file validation metrics
4. **RT window width adaptation** - No adaptive scaling based on refinement quality

---

## 6. Recommended Fix Plan

### Phase 1: Fix Critical Bug ✅ COMPLETE

**Priority**: CRITICAL
**Status**: **COMPLETED** on 2025-01-21 (commit `37fd89e0`)
**Time taken**: ~30 minutes
**Testing**: Integration testing pending

#### 1.1 Code Fix ✅ COMPLETE
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`
**Commit**: `37fd89e0`

**Line 108**:
```julia
# BEFORE:
best_prob = prob,

# AFTER:
best_prob = best_prob,
```

**Line 111**:
```julia
# BEFORE:
best_refined_irt = refined_irt,

# AFTER:
best_refined_irt = best_refined_irt,
```

#### 1.2 Validation Test (RECOMMENDED)
Add test to verify correct behavior:
```julia
# Pseudo-code test
prec_to_best = Dictionary{UInt32, NamedTuple}()
# Add PSM with prob=0.9, refined_irt=100
# Add PSM for same precursor with prob=0.5, refined_irt=200
# Assert: prec_to_best[prec_id].best_prob == 0.9
# Assert: prec_to_best[prec_id].best_refined_irt == 100
```

#### 1.3 Integration Test (PENDING)
Run on small dataset to verify:
1. Identification rates improve in SecondPassSearch
2. Check that `best_refined_irt` values are consistent across files
3. Confirm RT window widths are reasonable

### Phase 2: Improve Observability (HIGH PRIORITY)

**Priority**: HIGH
**Estimated time**: 2-3 hours
**Testing**: Verify log output on test dataset

#### 2.1 Upgrade Fallback Warning
**File**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl` line 524

```julia
# BEFORE:
@debug "Refined iRT model not found for file $index, falling back to library iRT model"

# AFTER:
@user_warn "Refined iRT model not found for file $index, falling back to library iRT model. This may reduce identification accuracy."
```

#### 2.2 Add RT Alignment Statistics Logging
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

Add after Step 6 (line ~490):
```julia
if !isnothing(refinement_model) && refinement_model.use_refinement
    @user_info "  RT Alignment Statistics:"
    @user_info "    Refinement MAE improvement: $(refinement_model.mae_improvement)%"
    @user_info "    Training PSMs: $(refinement_model.n_train)"
    @user_info "    Validation PSMs: $(refinement_model.n_val)"
    # Add window width comparison if available
end
```

#### 2.3 Add RT Window Width Logging
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

In `get_irt_errs()` function (after line ~760):
```julia
@user_info "RT Window Statistics for file $file_idx:"
@user_info "  Peak width (FWHM): $peak_width"
@user_info "  Cross-run std: $irt_std"
@user_info "  Total tolerance: $irt_tol"
```

#### 2.4 Track Identification Rates
Add PSM count logging after each search method:
```julia
@user_info "PSM counts for file $file_idx:"
@user_info "  FirstPassSearch: $first_pass_count"
@user_info "  SecondPassSearch: $second_pass_count"
@user_info "  Improvement: $(100*(second_pass_count/first_pass_count - 1))%"
```

### Phase 3: Optional Enhancements (LOWER PRIORITY)

**Priority**: MEDIUM
**Estimated time**: 4-6 hours
**Testing**: Benchmark on multiple datasets

#### 3.1 Adaptive Tolerance Scaling
Adjust `irt_tol` based on refinement quality:
```julia
# If refinement reduced MAE significantly, can use tighter tolerance
refinement_factor = min(refined_mae / library_mae, 1.0)
adaptive_irt_tol = irt_tol * (0.7 + 0.3 * refinement_factor)
```

#### 3.2 Add Diagnostic Plots
**File**: New file `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/rt_alignment_plots.jl`

Generate plots comparing:
1. FWHM distributions: library iRT vs refined iRT
2. Cross-run variance: library iRT vs refined iRT
3. RT window widths: library vs refined
4. MAE distributions: before vs after refinement

#### 3.3 Cross-File Validation
Monitor per-file refinement success:
```julia
successful_refinements = count(files where use_refinement == true)
total_files = length(files)
@user_info "Refinement success rate: $successful_refinements / $total_files"
```

#### 3.4 Add Overfitting Checks
Log validation metrics per file:
```julia
@user_info "Refinement validation for file $file_idx:"
@user_info "  Train MAE: $train_mae"
@user_info "  Val MAE: $val_mae"
@user_info "  Generalization gap: $(val_mae - train_mae)"
```

### Phase 4: Documentation Updates (ONGOING)

**Priority**: MEDIUM
**Estimated time**: 1-2 hours

#### 4.1 Update CLAUDE.md Files
Add section on iRT refinement to:
- `src/Routines/SearchDIA/CLAUDE.md`
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/CLAUDE.md`

#### 4.2 Add Inline Documentation
Expand comments in key functions:
- `map_retention_times!()` - Explain all 6 steps
- `fit_irt_refinement_model()` - Document algorithm and validation
- `getRtToRefinedIrtModel()` - Explain fallback behavior

#### 4.3 Create User Documentation
Add section to user manual explaining:
- What is iRT refinement
- When it's beneficial
- How to interpret refinement metrics
- Troubleshooting refinement failures

---

## 7. Critical Imputation Bug Discovery (2025-01-21) 🔴

### 7.1 The Discovery

After fixing the critical bug in `getBestPrecursorsAccrossRuns.jl` (commit `37fd89e0`), testing revealed:

- ✅ **Bug fix works**: Significant increase in identifications on both branches
- ✅ **Same improvement on develop**: Bug fix on `develop` branch (library iRT only) gives **equal** ID boost
- ❌ **Refinement doesn't help**: Same number of IDs on `develop` vs `refine_library_irt` branch
- ⚠️ **Refinement might hurt**: Slightly **fewer** IDs with refinement compared to library iRT alone

**Key Insight**: The refinement model is working (50% MAE reduction validated), but something in the pipeline is preventing this improvement from translating into better identifications.

### 7.2 Root Cause Identified: Wrong Refinement Model for Imputation

**Location**: `makeRTIndices()` in `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl` lines 88-105

**Current imputation logic** (INCORRECT):
```julia
for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
    prec_ids[i] = prec_id
    irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}

    # High-confidence empirical matches
    if haskey(prec_set, prec_id)
        _irt_, prob = prec_set[prec_id]
        if (prob >= min_prob)
            irts[i] = _irt_  # Empirical refined iRT from THIS file ✅
            continue
        end
    end

    # IMPUTATION: Precursor not observed (or low confidence) in this file
    # ❌ BUG: Uses best_refined_irt from prec_to_irt dictionary
    irts[i] = irt  # irt = best_refined_irt from DIFFERENT file!
    mzs[i] = mz
end
```

Where `prec_to_irt` comes from (FirstPassSearch/utils.jl line 658):
```julia
# Maps each precursor to its best_refined_irt across ALL files
prec_to_irt = map(x -> (irt=x[:best_refined_irt], mz=x[:mz]), precursor_dict)
```

**The problem**:

1. `best_refined_irt` is the observed refined iRT from **File A** (whichever file had the best probability match for this precursor)
2. That refined iRT was calculated using **File A's refinement model**: `refinement_model_A(sequence, library_irt)`
3. We're now imputing it for **File B**, which has **different RT behavior** and a **different refinement model**: `refinement_model_B`
4. **File A's refined iRT doesn't account for File B's RT characteristics!**

### 7.3 Concrete Example

```julia
# Precursor 12345, sequence "PEPTIDE", library_irt = 45.0

# File A (best match, prob=0.95):
refinement_model_A(sequence="PEPTIDE", library_irt=45.0) = 50.0
# Observed refined_irt_A = 50.0
# This becomes best_refined_irt = 50.0

# File B (not observed or low prob):
# Current (WRONG): Use best_refined_irt = 50.0 from File A
irts[i] = 50.0  # ❌ Wrong! Uses File A's model for File B!

# Correct: Calculate using File B's model
refinement_model_B(sequence="PEPTIDE", library_irt=45.0) = 48.0  # Different!
irts[i] = 48.0  # ✅ Correct! File-specific prediction
```

**Why this causes problems**:
- File A and File B have different RT calibration, drift, or instrument characteristics
- Refinement models learn file-specific AA composition effects
- Using File A's refined iRT for File B creates **systematic errors**
- These errors can be **worse than using library iRT** (which is at least consistent across files)

### 7.4 Correct Solution ✅

**What we should do**:
```julia
# In makeRTIndices(), when building RT index for a specific file:
function makeRTIndices(
    temp_folder::String,
    psms_paths::Vector{String},
    prec_to_irt::Dictionary,  # Contains library iRT now, not best_refined_irt
    rt_to_refined_irt_splines::Any,
    refinement_models::Dict{Int, IrtRefinementModel},  # NEW: Pass refinement models
    precursors::PrecursorLibrary;  # NEW: Need access to sequences
    min_prob::AbstractFloat = 0.5
)
    for (key, psms_path) in enumerate(psms_paths)
        rt_to_refined_irt = rt_to_refined_irt_splines[key]
        refinement_model = refinement_models[key]  # Get THIS file's model

        for (i, (prec_id, irt_mz)) in enumerate(pairs(prec_to_irt)))
            # High-confidence: use empirical refined iRT
            if haskey(prec_set, prec_id) && prob >= min_prob
                irts[i] = rt_to_refined_irt(scan_rt)  # Empirical ✅
                continue
            end

            # IMPUTATION: Calculate file-specific refined iRT prediction
            library_irt = irt_mz.irt  # Now contains library iRT
            sequence = precursors[prec_id].sequence

            # Use THIS file's refinement model
            if !isnothing(refinement_model) && refinement_model.use_refinement
                irts[i] = refinement_model(sequence, library_irt)  # ✅ File-specific!
            else
                irts[i] = library_irt  # Fallback to library iRT
            end
            mzs[i] = mz
        end
    end
end
```

**Why this works**:
1. **File-specific**: Uses File B's refinement model for File B's RT index
2. **Consistent**: Same refinement approach as empirical matches
3. **Available**: All inputs already exist (library iRT, refinement model, sequence)
4. **Better**: Eliminates cross-file model mismatch

### 7.5 Implementation Changes Required

**Three files need modification**:

1. **FirstPassSearch/utils.jl** line 658:
```julia
# BEFORE:
prec_to_irt = map(x -> (irt=x[:best_refined_irt], mz=x[:mz]), precursor_dict)

# AFTER:
# Pass library iRT instead of best_refined_irt
prec_to_irt = map(x -> (irt=x[:library_irt], mz=x[:mz]), precursor_dict)
# OR keep both:
prec_to_irt = map(x -> (library_irt=x[:library_irt], best_refined_irt=x[:best_refined_irt], mz=x[:mz]), precursor_dict)
```

2. **FirstPassSearch/utils.jl** line 689 (create_rt_indices!):
```julia
# BEFORE:
rt_index_paths = makeRTIndices(
    rt_indices_folder,
    valid_psm_paths,
    prec_to_irt,
    valid_rt_models,
    min_prob=params.max_prob_to_impute
)

# AFTER:
rt_index_paths = makeRTIndices(
    rt_indices_folder,
    valid_psm_paths,
    prec_to_irt,  # Now contains library_irt
    valid_rt_models,
    getRefinementModels(search_context),  # NEW: Pass refinement models
    getPrecursors(search_context),  # NEW: Pass precursor library
    min_prob=params.max_prob_to_impute
)
```

3. **buildRTIndex.jl** lines 69-113:
```julia
# Update function signature and implementation as shown in Section 7.4
```

### 7.6 Expected Impact

This fix should:

1. ✅ **Eliminate cross-file model mismatch**: Each file uses its own refinement model
2. ✅ **Provide file-specific refined iRT for ALL precursors**: Not just empirical matches
3. ✅ **Allow refinement to improve IDs as expected**: Remove systematic imputation errors
4. ✅ **Quick win**: Simple fix, high probability of success

**Priority**: **CRITICAL - Test this BEFORE other hypotheses**

**Estimated time**: 1-2 hours (implementation + testing)

**Success criteria**:
- Refined branch shows 10-30% more IDs than library branch
- Improvement scales with refinement MAE improvement
- Cross-run variance actually decreases (not increases)

### 7.7 Why This Bug Wasn't Obvious

This bug was subtle because:

1. The code "works" - doesn't crash, produces refined iRT values
2. Refinement MAE validation shows improvement (within same file)
3. Bug only manifests when comparing cross-file imputation quality
4. Effect depends on how different refinement models are across files
5. If files have similar RT behavior, bug impact is small

**User question that revealed it**: "Why are we imputing from best_refined_irt when we have the refinement model for this file?"

This question correctly identified that we should be **calculating** file-specific refined iRT, not **borrowing** another file's refined iRT.

---

## Appendix A: Key Files Reference

### Core RT Conversion
- `src/structs/RetentionTimeConversionModel.jl` - RT/iRT conversion model types
- `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl` - Refinement model training

### FirstPassSearch
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl` - RT alignment pipeline
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl` - ❌ Bug location

### RT Index Building
- `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl` - RT index creation

### Post-FirstPassSearch Searches
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl` - Feature calculation
- `src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl` - RT window calculations
- `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl` - RT window calculations

### Scoring Pipeline
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl` - Feature binning
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl` - Feature sets
- `src/utils/ML/percolatorSortOf.jl` - ML training

---

## Appendix B: Commit History

### This Branch (refine_library_irt)
```
37fd89e0 fix: Store best PSM values instead of current values in getBestPrecursorsAccrossRuns ✅ CRITICAL BUG FIX
6511aaa5 docs: Add comprehensive RT handling analysis for refine_library_irt branch
e5608e3a fix: Update DataFrame column references from :irt_pred/:irt_obs to :refined_irt_pred/:refined_irt_obs
de0138ed refactor: Implement refined iRT model usage throughout search pipeline
b451cfd4 docs: Add plan for fixing rt_to_irt vs rt_to_refined_irt usage
df444eba fix: Complete field name updates for iRT refinement (best_refined_irt, mean_refined_irt, var_refined_irt)
f13fcb44 CHECKPOINT: Before iRT refinement implementation
```

### Related Commits (from develop)
```
1e56ac5c fix: Implement file-specific iRT refinement (CRITICAL BUG FIX)
4cc092f9 Rename all prob/probs variables to trace_prob/trace_probs
f5dfbaae speed up getBestPrecursorsAccrossRuns
```

---

## Appendix C: Testing Checklist

### After Bug Fix
- [ ] Verify `best_prob` stores highest probability across files
- [ ] Verify `best_refined_irt` matches refined iRT from best PSM
- [ ] Check cross-run variance calculations are reasonable
- [ ] Run full pipeline on test dataset (e.g., ecoli_test)
- [ ] Compare identification rates: FirstPass vs SecondPass
- [ ] Verify RT window widths are reasonable
- [ ] Check for any regressions in existing functionality

### After Observability Improvements
- [ ] Verify refinement statistics are logged
- [ ] Check RT window width logging is clear
- [ ] Confirm fallback warnings appear when expected
- [ ] Review identification rate comparisons

### Before Merge to Develop
- [ ] All tests pass
- [ ] No performance regressions
- [ ] Documentation updated
- [ ] Code review completed
- [ ] Benchmark results show improvement

---

## Status Update

**Date**: 2025-01-21 (Latest update)
**Current Status**: First bug fixed, second critical bug discovered

### Completed
- ✅ First critical bug analysis (Section 1)
- ✅ First bug fix implementation (commit `37fd89e0`)
- ✅ Integration testing revealed refinement doesn't improve IDs
- ✅ Second critical bug identified (Section 7 - RT index imputation)
- ✅ Documentation updates (this file + IRT_REFINEMENT_INVESTIGATION_PLAN.md)
- ✅ Code pushed to remote branch `refine_library_irt`

### Current Findings
- ✅ **First bug fix works**: Significant ID improvement on both branches
- ❌ **Refinement doesn't help**: Same IDs on develop (library iRT) vs refine_library_irt (refined iRT)
- 🔴 **New bug discovered**: RT index imputation uses wrong refinement model (cross-file mismatch)

### Pending (HIGH PRIORITY)
- 🔴 **CRITICAL**: Fix RT index imputation bug (Section 7.5)
  - Estimated time: 1-2 hours
  - Expected impact: Allow refinement to improve IDs as expected
  - Priority: Test this BEFORE other hypotheses
- ⏳ Integration testing after imputation fix
- ⏳ Validation that refinement improves IDs (10-30% expected)
- ⏳ Phase 2 observability improvements (if needed)

### Root Cause Analysis
**Why refinement doesn't improve IDs despite 50% MAE reduction:**

The imputation bug causes systematic cross-file model mismatch:
1. Precursor observed in File A with prob=0.95 → refined_irt_A = model_A(seq, lib_irt) = 50.0
2. Precursor NOT observed in File B → **uses refined_irt_A = 50.0** ❌ WRONG!
3. **Should calculate**: refined_irt_B = model_B(seq, lib_irt) = 48.0 ✅ CORRECT

**Impact**: File B's RT index uses File A's refined iRT (wrong model), creating systematic errors that can be **worse than library iRT** (explaining why refinement hurts IDs).

### Next Steps (Updated)
1. **IMMEDIATE**: Implement imputation fix (Section 7.5)
   - Modify `makeRTIndices()` to use file-specific refinement models
   - Pass refinement models and precursor library to function
   - Calculate refined iRT predictions instead of borrowing from other files
2. Run full pipeline on test dataset
3. Validate that refinement now improves IDs (expect 10-30% improvement)
4. If successful, consider Phase 2 logging for monitoring
5. Create pull request to merge into `develop` branch

### Expected Outcome (After Imputation Fix)
With both bugs fixed:
- RT indices use correct refined iRT values (best PSM, not most recent) ✅
- RT indices use file-specific refined iRT (File B's model for File B, not File A's) ✅
- Refinement improvements (50% MAE reduction) should translate to 10-30% more IDs
- Improvement should scale with refinement MAE quality per file

---

**End of Analysis**
