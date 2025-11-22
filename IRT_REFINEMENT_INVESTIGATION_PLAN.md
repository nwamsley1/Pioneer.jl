# iRT Refinement Investigation Plan
## Why Refinement Doesn't Improve Identifications Despite 50% MAE Reduction

**Date**: 2025-01-21
**Branch**: `refine_library_irt`
**Status**: Critical bug fixed, but refinement not providing expected benefits

---

## Executive Summary

### The Problem

After fixing the critical bug in `getBestPrecursorsAccrossRuns.jl` (commit `37fd89e0`):

- ✅ **Bug fix works**: Significant increase in identifications
- ✅ **Same improvement on develop**: Bug fix on develop branch (library iRT only) gives equal ID boost
- ❌ **iRT refinement not helping**: Despite 50% MAE reduction, refinement doesn't translate to more IDs
- ⚠️ **Refinement might hurt**: Slightly fewer IDs with refinement compared to library iRT alone

### Key Question

**If iRT refinement reduces MAE by 50%, why doesn't it improve identifications?**

The refinement model is working (validated MAE improvement), but something in the pipeline is preventing this improvement from translating into better results.

---

## Possible Root Causes

### Hypothesis 1: RT Windows Too Narrow ⚠️ HIGH LIKELIHOOD

**Theory**: Refinement reduces variance too much, making RT windows too narrow and missing valid matches.

**Evidence from RT_HANDLING_ANALYSIS.md**:
- RT tolerance formula: `irt_tol = peak_width + irt_nstd * sqrt(var_refined_irt)`
- If refinement reduces `var_refined_irt` by 50%, tolerance could drop proportionally
- Narrower windows → miss valid PSMs that have slightly different RT behavior

**Why this could happen**:
- Refinement model trained on same file it's applied to (67/33 split)
- Reduces variance within training data, but not necessarily cross-run variance
- Could be overfitting to specific RT behavior in validation set

**Test**: Compare RT window widths between library and refined iRT
- Expected: Refined windows significantly narrower (>30% reduction)
- Impact: Missing 10-20% of valid matches

---

### Hypothesis 2: Refinement Model Overfitting ⚠️ MEDIUM LIKELIHOOD

**Theory**: Model overfits to training data and doesn't generalize well.

**Evidence**:
- 67/33 train/val split on **same file**
- Validation MAE improves, but only on held-out PSMs from same file
- No cross-file validation or cross-run generalization test

**Why this could happen**:
- Model learns file-specific RT quirks rather than general AA composition effects
- Validation set has similar RT distribution to training set (same file)
- Different scans in same file might have different RT calibration drift

**Test**: Check if refinement helps or hurts cross-file imputation
- Compare refined iRT variance across files vs library iRT variance
- Expected: Refined variance could be higher (worse) if overfitting

---

### Hypothesis 3: Feature Quality Issues ⚠️ MEDIUM LIKELIHOOD

**Theory**: Refined iRT features are noisier than library iRT features, hurting ML scoring.

**Evidence**:
- ML scoring uses `:refined_irt_error`, `:irt_diff`, `:refined_irt_pred_qbin`
- If refined predictions are less consistent, these features become less informative
- LightGBM might learn to downweight or ignore refined iRT features

**Why this could happen**:
- Refinement adds model uncertainty on top of RT measurement uncertainty
- `irt_diff` uses `best_refined_irt` which varies by file (file-specific model)
- Inconsistent refinement across files → noisy features

**Test**: Compare feature AUC (target/decoy separation)
- Calculate AUC for `library_irt_error` vs `refined_irt_error`
- Expected: Refined features have lower AUC (worse discrimination)

---

### Hypothesis 4: RT Index Imputation Problems ⚠️ LOW-MEDIUM LIKELIHOOD

**Theory**: Imputed refined iRT values are worse than imputed library iRT values.

**Evidence**:
- Low-confidence PSMs get imputed `best_refined_irt` from cross-run best match
- `best_refined_irt` is file-specific (different refinement model per file)
- Imputing file A's refined iRT for use in file B might not work well

**Why this could happen**:
- Refinement models differ across files (different AA composition effects)
- File A's best_refined_irt doesn't account for file B's RT behavior
- Library iRT is more consistent across files (same for all files)

**Test**: Compare imputation accuracy
- For imputed precursors, check RT prediction error
- Expected: Refined imputation has higher error than library imputation

---

## Investigation Plan

### Phase 1: Data Collection & Diagnostic Logging 🔍

**Priority**: HIGH
**Time estimate**: 2-3 hours implementation + 30 min analysis
**Goal**: Understand what's happening during refinement and search

#### 1.1 RT Window Width Comparison

**Where**: `get_irt_errs()` in `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl` (lines 721-763)

**Add logging**:
```julia
# After calculating irt_errs for each file
for (file_idx, irt_err) in pairs(irt_errs)
    # Get library iRT statistics for comparison
    library_variance = calculate_library_irt_variance(prec_to_irt, file_idx)
    refined_variance = calculate_refined_irt_variance(prec_to_irt, file_idx)

    # Calculate what library iRT tolerance would be
    library_irt_tol = peak_width + irt_nstd * sqrt(library_variance)
    refined_irt_tol = irt_err

    @user_info "RT Window Statistics (File $file_idx):"
    @user_info "  Peak width (FWHM): $peak_width"
    @user_info "  Library iRT cross-run std: $(sqrt(library_variance))"
    @user_info "  Refined iRT cross-run std: $(sqrt(refined_variance))"
    @user_info "  Library iRT tolerance: $library_irt_tol"
    @user_info "  Refined iRT tolerance: $refined_irt_tol"
    @user_info "  Window width ratio (refined/library): $(refined_irt_tol/library_irt_tol)"

    if (refined_irt_tol/library_irt_tol) < 0.7
        @user_warn "  ⚠️  Refined windows are >30% narrower - may miss valid matches"
    end
end
```

**Analysis**:
- If ratio consistently < 0.7: **Hypothesis 1 confirmed** (narrow windows)
- If ratio ≈ 1.0: Windows are similar, problem is elsewhere
- If ratio > 1.0: Refinement is widening windows (unexpected, investigate)

#### 1.2 Refinement Model Quality Metrics

**Where**: `fit_irt_refinement_model()` in `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl` (lines 106-220)

**Add logging** (after line ~210):
```julia
@user_info "iRT Refinement Model Training (File $ms_file_idx):"
@user_info "  Training set: $n_train PSMs"
@user_info "  Validation set: $n_val PSMs"
@user_info "  Training MAE: $(round(train_mae, digits=3))"
@user_info "  Validation MAE: $(round(val_mae, digits=3))"
@user_info "  Baseline MAE: $(round(baseline_mae, digits=3))"
@user_info "  MAE improvement: $(round(100*(1 - val_mae/baseline_mae), digits=1))%"
@user_info "  Model coefficients range: [$(minimum(model.coefficients)), $(maximum(model.coefficients))]"
@user_info "  Use refinement: $use_refinement"

if use_refinement
    generalization_gap = val_mae - train_mae
    if generalization_gap > 0.1
        @user_warn "  ⚠️  Large generalization gap ($generalization_gap) - possible overfitting"
    end
end
```

**Analysis**:
- Large train/val gap: **Hypothesis 2 confirmed** (overfitting)
- Consistent improvement across files: Model is working as expected
- Some files fail validation: Check why (not enough PSMs? poor RT alignment?)

#### 1.3 Identification Rates Per Search Method

**Where**: After each search method in `SearchDIA.jl`

**Add logging**:
```julia
# After FirstPassSearch
first_pass_counts = count_psms_per_file(first_pass_results)
@user_info "FirstPassSearch Identifications:"
for (file_idx, count) in first_pass_counts
    @user_info "  File $file_idx: $count PSMs"
end

# After SecondPassSearch
second_pass_counts = count_psms_per_file(second_pass_results)
@user_info "SecondPassSearch Identifications:"
for (file_idx, count) in enumerate(second_pass_counts)
    improvement = count - first_pass_counts[file_idx]
    pct_improvement = 100 * improvement / first_pass_counts[file_idx]
    @user_info "  File $file_idx: $count PSMs (+$improvement, +$(round(pct_improvement, digits=1))%)"
end

# After ScoringSearch
scoring_counts = count_psms_per_file(scoring_results)
@user_info "After Scoring (FDR filtered):"
for (file_idx, count) in enumerate(scoring_counts)
    @user_info "  File $file_idx: $count PSMs"
end
```

**Analysis**:
- Compare SecondPass improvement: refine branch vs develop branch
- If develop shows bigger improvement: Refinement is hurting SecondPass
- If both similar: Problem is elsewhere (feature quality, scoring, etc.)

#### 1.4 Feature Distribution Analysis

**Where**: `add_features!()` in `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl` (lines 829-971)

**Add logging** (after line ~970):
```julia
# Calculate feature statistics
@user_info "iRT Feature Statistics (File $ms_file_idx):"
@user_info "  refined_irt_error: mean=$(round(mean(irt_error), digits=3)), std=$(round(std(irt_error), digits=3))"
@user_info "  irt_diff: mean=$(round(mean(irt_diff), digits=3)), std=$(round(std(irt_diff), digits=3))"
@user_info "  ms1_irt_diff: mean=$(round(mean(filter(!isnan, ms1_irt_diff)), digits=3)), std=$(round(std(filter(!isnan, ms1_irt_diff)), digits=3))"

# Calculate feature quality (AUC for target/decoy separation)
target_irt_error = irt_error[target]
decoy_irt_error = irt_error[.!target]
auc_irt_error = calculate_auc(target_irt_error, decoy_irt_error)
@user_info "  refined_irt_error AUC (target/decoy separation): $(round(auc_irt_error, digits=3))"

# Compare to what library iRT error would be (if available)
if haskey(search_context.rt_irt_map, ms_file_idx)
    library_irt_error = calculate_library_irt_error(...)
    auc_library = calculate_auc(library_irt_error[target], library_irt_error[.!target])
    @user_info "  library_irt_error AUC: $(round(auc_library, digits=3))"
    @user_info "  AUC improvement: $(round(100*(auc_irt_error - auc_library)/auc_library, digits=1))%"
end
```

**Analysis**:
- If refined AUC < library AUC: **Hypothesis 3 confirmed** (feature quality)
- If refined AUC > library AUC: Features are better, problem is elsewhere
- High std in refined features: More noise, could hurt ML training

---

### Phase 2: Controlled Comparison Experiments 🧪

**Priority**: CRITICAL
**Time estimate**: 3-4 hours
**Goal**: Direct comparison of library vs refined iRT on same data

#### 2.1 Compare RT Window Precursor Coverage

**Method**: Instrument RT window selection code

**Where**: SecondPassSearch/utils.jl lines 181-183, 386-388

**Implementation**:
```julia
# Before RT window calculation
total_precursors_in_index = length(rt_index.precursors)

# During RT window iteration
precursors_in_window = count_precursors_in_range(rt_index, irt_start, irt_stop)
precursor_coverage = precursors_in_window / total_precursors_in_index

# Log every N scans
if scan_idx % 100 == 0
    @debug "Scan $scan_idx: $precursors_in_window/$total_precursors_in_index precursors in RT window ($(round(100*precursor_coverage, digits=1))%)"
end
```

**Analysis**:
- Average coverage: refined vs library
- If refined coverage < library by >20%: Narrow windows are the problem
- If coverage similar: Window size is not the issue

#### 2.2 PSM Overlap Analysis

**Method**: Compare PSMs found in both branches

**Implementation**:
```julia
# Load PSMs from both branches
refined_psms = load_psms("refine_branch/second_pass/")
library_psms = load_psms("develop_branch/second_pass/")

# Find overlap
refined_ids = Set(zip(refined_psms.precursor_idx, refined_psms.scan_idx))
library_ids = Set(zip(library_psms.precursor_idx, library_psms.scan_idx))

overlap = length(intersect(refined_ids, library_ids))
refined_only = length(setdiff(refined_ids, library_ids))
library_only = length(setdiff(library_ids, refined_ids))

@user_info "PSM Overlap Analysis:"
@user_info "  Refined branch PSMs: $(length(refined_ids))"
@user_info "  Library branch PSMs: $(length(library_ids))"
@user_info "  Shared PSMs: $overlap ($(round(100*overlap/length(library_ids), digits=1))%)"
@user_info "  Refined-only PSMs: $refined_only"
@user_info "  Library-only PSMs: $library_only"
```

**Analysis**:
- Large library_only count: Refinement is missing valid PSMs (precision problem)
- Large refined_only count: Refinement finds different PSMs (investigate quality)
- Similar totals but different PSMs: Fundamental difference in RT targeting

#### 2.3 Feature Quality Comparison

**Method**: Compare discriminative power of features

**Implementation**:
```julia
# For both refined and library PSMs, calculate feature AUCs
features_to_compare = [
    :irt_error,
    :irt_diff,
    :spectral_contrast,
    :hyperscore
]

for feature in features_to_compare
    refined_auc = calculate_feature_auc(refined_psms, feature)
    library_auc = calculate_feature_auc(library_psms, feature)

    @user_info "Feature $feature AUC:"
    @user_info "  Refined: $(round(refined_auc, digits=3))"
    @user_info "  Library: $(round(library_auc, digits=3))"
    @user_info "  Difference: $(round(refined_auc - library_auc, digits=3))"
end
```

**Analysis**:
- If refined irt_error AUC worse: Features are noisier
- If other features similar: Problem is iRT-specific
- If all features worse in refined: Something fundamentally wrong

---

### Phase 3: Hypothesis Testing 🔬

**Priority**: MEDIUM
**Time estimate**: 1-2 hours per test
**Goal**: Confirm root cause with targeted experiments

#### 3.1 Test Hypothesis 1: Narrow Windows

**Method**: Artificially widen refined iRT windows

**Implementation** in `get_irt_errs()`:
```julia
# Add scaling factor to refined iRT tolerance
if use_refined_irt
    # Current calculation
    refined_irt_tol = peak_width + irt_nstd * irt_std

    # Test with wider windows
    window_scale_factor = 1.5  # 50% wider
    refined_irt_tol = refined_irt_tol * window_scale_factor

    @user_info "  Applied window scaling factor: $window_scale_factor"
    @user_info "  Adjusted refined iRT tolerance: $refined_irt_tol"
end
```

**Expected results**:
- If narrow windows are the problem: IDs should increase by 10-30%
- If no change: Problem is elsewhere
- If IDs decrease: Windows were already optimal

**Conclusion**:
- IDs improve → Implement adaptive tolerance scaling (see Phase 4.1)
- No change → Test next hypothesis

#### 3.2 Test Hypothesis 2: Overfitting / Imputation

**Method**: Use library iRT for imputation instead of refined iRT

**Implementation** in `makeRTIndices()` (buildRTIndex.jl lines 88-105):
```julia
for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
    prec_ids[i] = prec_id
    irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}

    if haskey(prec_set, prec_id)
        _irt_, prob = prec_set[prec_id]
        if (prob >= min_prob)
            # High confidence: use empirical refined iRT
            irts[i] = _irt_
            continue
        end
    end

    # LOW CONFIDENCE OR MISSING: Test using library iRT instead of refined
    # Original: irts[i] = irt  # refined iRT from best match
    # Test: Use library iRT instead
    library_irt = get_library_irt(precursors, prec_id)
    irts[i] = library_irt  # TESTING: Use library iRT for imputation
    mzs[i] = mz
end
```

**Expected results**:
- If imputation is the problem: IDs should improve
- If no change: Imputation is fine, problem is elsewhere

**Conclusion**:
- IDs improve → Use library iRT for imputation (hybrid approach)
- No change → Test next hypothesis

#### 3.3 Test Hypothesis 3: Feature Quality

**Method**: Use library iRT for features, refined iRT for RT windows

**Implementation** in `add_features!()`:
```julia
# Calculate BOTH library and refined iRT features
rt_to_library_irt = getRtIrtModel(search_context, ms_file_idx)
rt_to_refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)

for i in 1:N
    # Refined iRT for RT windows (keep current usage)
    refined_irt_obs[i] = rt_to_refined_irt(rt[i])
    refined_irt_pred[i] = refinement_model(sequence, library_irt)

    # TESTING: Calculate library iRT features for ML
    library_irt_obs = rt_to_library_irt(rt[i])
    library_irt_pred = getPredIrt(search_context, prec_idx)

    # Use library features instead of refined
    irt_error[i] = abs(library_irt_obs - library_irt_pred)  # TESTING
    irt_diff[i] = abs(library_irt_obs - best_library_irt)   # TESTING
end

# Store library features for ML
psms[!,:irt_error] = irt_error  # Now using library, not refined
psms[!,:irt_diff] = irt_diff    # Now using library, not refined
```

**Expected results**:
- If feature quality is the problem: IDs should improve
- If no change: Features are fine, problem is RT windows

**Conclusion**:
- IDs improve → Implement hybrid approach (see Phase 4.3)
- No change → Problem is likely RT window width (Hypothesis 1)

---

### Phase 4: Root Cause Solutions 🛠️

**Priority**: HIGH (after confirming root cause)
**Time estimate**: 2-4 hours per solution
**Goal**: Implement proper fix based on confirmed hypothesis

#### 4.1 Solution for Hypothesis 1: Adaptive Tolerance Scaling

**Problem**: Refined iRT windows too narrow, missing valid matches

**Solution**: Scale RT tolerance based on refinement quality

**Implementation**:
```julia
function get_irt_errs(
    fwhms::Dictionary,
    prec_to_irt::Dictionary,
    params::FirstPassSearchParameters,
    refinement_quality::Dictionary{Int64, Float32}  # NEW: MAE improvement per file
)
    irt_errs = Dictionary{Int64, Float32}()

    for (file_idx, fwhm_stats) in pairs(fwhms)
        # Base tolerance calculation
        peak_width = fwhm_stats.median_fwhm + params.fwhm_nstd * fwhm_stats.mad_fwhm
        refined_irt_std = calculate_cross_run_std(prec_to_irt, file_idx)
        base_tolerance = peak_width + params.irt_nstd * refined_irt_std

        # Adaptive scaling based on refinement quality
        if haskey(refinement_quality, file_idx)
            mae_improvement = refinement_quality[file_idx]  # e.g., 0.5 for 50% improvement

            # Scale factor: worse refinement → wider windows
            # Perfect refinement (100% improvement) → 1.0x (no scaling)
            # No improvement (0%) → 1.5x (50% wider)
            # Moderate improvement (50%) → 1.25x (25% wider)
            scale_factor = 1.0 + 0.5 * (1.0 - mae_improvement)

            adjusted_tolerance = base_tolerance * scale_factor

            @user_info "  File $file_idx adaptive tolerance:"
            @user_info "    MAE improvement: $(round(100*mae_improvement, digits=1))%"
            @user_info "    Scale factor: $(round(scale_factor, digits=2))x"
            @user_info "    Base tolerance: $(round(base_tolerance, digits=3))"
            @user_info "    Adjusted tolerance: $(round(adjusted_tolerance, digits=3))"

            irt_errs[file_idx] = adjusted_tolerance
        else
            irt_errs[file_idx] = base_tolerance
        end
    end

    return irt_errs
end
```

**Benefits**:
- Automatically widens windows when refinement is aggressive
- Maintains narrow windows when refinement is highly accurate
- File-specific adaptation

**Risks**:
- Could allow more false positives if too aggressive
- Need to tune scaling parameters

#### 4.2 Solution for Hypothesis 2: Improved Validation Strategy

**Problem**: Refinement model overfits to training file

**Solution 1**: Cross-file validation
```julia
function fit_irt_refinement_model_cross_validated(
    all_files_data::Vector{FileData},
    current_file_idx::Int
)
    # Train on current file
    train_data = all_files_data[current_file_idx]

    # Validate on OTHER files (cross-file generalization)
    val_files = [i for i in eachindex(all_files_data) if i != current_file_idx]
    val_data = vcat([all_files_data[i] for i in val_files]...)

    model = train_model(train_data)

    # Validate on different files
    cross_file_mae = evaluate_model(model, val_data)
    same_file_mae = evaluate_model(model, train_data.validation_set)

    @user_info "Cross-validation results:"
    @user_info "  Same-file MAE: $same_file_mae"
    @user_info "  Cross-file MAE: $cross_file_mae"
    @user_info "  Generalization gap: $(cross_file_mae - same_file_mae)"

    # Only use if cross-file performance is good
    use_refinement = cross_file_mae < baseline_mae * 0.9

    return model, use_refinement
end
```

**Solution 2**: Regularization
```julia
# Add L2 regularization to prevent overfitting
function fit_with_regularization(X, y; lambda=0.1)
    # Ridge regression instead of OLS
    β = (X'X + lambda * I) \ (X'y)
    return β
end
```

**Solution 3**: Stricter improvement threshold
```julia
# Only use refinement if improvement is substantial
use_refinement = (val_mae < baseline_mae * 0.8)  # Require >20% improvement
```

#### 4.3 Solution for Hypothesis 3: Hybrid Approach

**Problem**: Refined iRT features are noisier than library iRT features

**Solution**: Use refined iRT for RT windows, library iRT for ML features

**Implementation**:
```julia
# In SecondPassSearch - RT window calculation
# USE REFINED iRT for targeting (more accurate)
refined_irt = getRtToRefinedIrtModel(search_context, ms_file_idx)(scan_rt)
irt_start = find_rt_window_start(refined_irt, irt_tol)
irt_stop = find_rt_window_stop(refined_irt, irt_tol)

# In SecondPassSearch - Feature calculation
# USE LIBRARY iRT for features (more stable)
library_irt_obs = getRtIrtModel(search_context, ms_file_idx)(rt[i])
library_irt_pred = getPredIrt(search_context, prec_idx)

# Calculate features with library iRT
irt_error[i] = abs(library_irt_obs - library_irt_pred)
irt_diff[i] = abs(library_irt_obs - best_library_irt)

# Store in DataFrame with clear naming
psms[!,:library_irt_error] = irt_error  # For ML scoring
psms[!,:refined_irt] = refined_irt_obs  # For reference only
```

**Benefits**:
- Gets improved RT targeting from refinement
- Avoids noisy features in ML scoring
- Maintains backward compatibility

**Risks**:
- Adds complexity (two iRT systems)
- Need to update model_config.jl feature sets

#### 4.4 Solution for Hypothesis 4: Hybrid Imputation

**Problem**: Imputed refined iRT values don't generalize across files

**Solution**: Use library iRT for imputation, refined iRT for empirical matches

**Implementation** in `makeRTIndices()`:
```julia
for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
    if haskey(prec_set, prec_id)
        _irt_, prob = prec_set[prec_id]
        if (prob >= min_prob)
            # HIGH CONFIDENCE: Use empirical refined iRT
            irts[i] = rt_to_refined_irt(scan_rt)
            continue
        end
    end

    # LOW CONFIDENCE OR MISSING: Use library iRT for better generalization
    library_irt = precursors[prec_id].irt_predicted
    irts[i] = library_irt
    mzs[i] = mz
end
```

**Benefits**:
- Empirical refined iRT when available (most accurate)
- Library iRT for imputation (more consistent across files)
- Simple implementation

---

## Quick Comparison Experiment

**Before investing in full implementation, run this quick test:**

### Setup
1. Run SAME dataset on both branches:
   - `develop` branch (library iRT, with bug fix)
   - `refine_library_irt` branch (refined iRT, with bug fix)

2. Compare outputs at each stage:
   - FirstPassSearch PSMs
   - SecondPassSearch PSMs
   - Final identifications

### Key Metrics to Compare

```julia
# 1. RT Window Width
develop_rt_tolerance = get_average_tolerance(develop_irt_errs)
refine_rt_tolerance = get_average_tolerance(refine_irt_errs)
width_ratio = refine_rt_tolerance / develop_rt_tolerance

# 2. Precursor Coverage
develop_coverage = count_precursors_in_windows(develop_psms) / total_precursors
refine_coverage = count_precursors_in_windows(refine_psms) / total_precursors
coverage_ratio = refine_coverage / develop_coverage

# 3. PSM Overlap
shared = count_shared_psms(develop_psms, refine_psms)
develop_only = count_unique_psms(develop_psms, refine_psms)
refine_only = count_unique_psms(refine_psms, develop_psms)

# 4. Feature Quality
develop_irt_auc = calculate_auc(develop_psms.irt_error, develop_psms.target)
refine_irt_auc = calculate_auc(refine_psms.refined_irt_error, refine_psms.target)
auc_ratio = refine_irt_auc / develop_irt_auc

# Report
@user_info "Quick Comparison Results:"
@user_info "  RT window width ratio: $width_ratio"
@user_info "  Precursor coverage ratio: $coverage_ratio"
@user_info "  PSM overlap: $shared shared, $develop_only develop-only, $refine_only refine-only"
@user_info "  iRT feature AUC ratio: $auc_ratio"
```

### Interpretation

| Metric | Observation | Likely Hypothesis |
|--------|-------------|-------------------|
| width_ratio < 0.7 | Refined windows 30%+ narrower | **Hypothesis 1** (narrow windows) |
| coverage_ratio < 0.8 | Covering 20%+ fewer precursors | **Hypothesis 1** (narrow windows) |
| develop_only >> refine_only | Missing many valid PSMs | **Hypothesis 1** (narrow windows) |
| auc_ratio < 0.9 | Features 10%+ worse | **Hypothesis 3** (feature quality) |
| Large generalization gap | Train/val difference significant | **Hypothesis 2** (overfitting) |

**Decision tree**:
- If `width_ratio < 0.7` AND `coverage_ratio < 0.8`: → Implement Solution 4.1 (adaptive tolerance)
- If `auc_ratio < 0.9`: → Implement Solution 4.3 (hybrid approach)
- If generalization gap large: → Implement Solution 4.2 (better validation)
- If multiple issues: → Implement hybrid solution combining fixes

---

## Recommended Immediate Next Steps

### Step 1: Update RT_HANDLING_ANALYSIS.md ⏰ 30 minutes
- Document new finding: bug fix works, refinement doesn't help
- Add this investigation plan as appendix or reference

### Step 2: Quick Comparison Experiment ⏰ 1 hour
- Run same dataset on both branches
- Collect the 4 key metrics above
- Identify most likely hypothesis

### Step 3: Implement Phase 1 Logging ⏰ 2-3 hours
- Add logging to all 4 sections
- Run test search with logging enabled
- Review logs to confirm hypothesis

### Step 4: Test Most Likely Hypothesis ⏰ 1-2 hours
- Implement one of the Phase 3 tests
- Run test search
- Compare results

### Step 5: Implement Proper Solution ⏰ 2-4 hours
- Based on confirmed hypothesis
- Implement appropriate Phase 4 solution
- Test on full dataset

### Step 6: Documentation & PR ⏰ 1-2 hours
- Document findings and solution
- Update analysis documents
- Create PR with fix

**Total estimated time**: 1-2 days for complete investigation and fix

---

## Success Criteria

**Minimum success**: Understand why refinement doesn't help
- Identify root cause with confidence
- Document findings for future work

**Good success**: Fix identified issue
- Refinement provides 5-15% ID improvement
- Improvement scales with refinement quality

**Excellent success**: Refinement provides expected benefits
- Refinement provides 20-40% ID improvement
- Validates iRT refinement approach
- Opens door for further RT modeling improvements

---

**End of Investigation Plan**
