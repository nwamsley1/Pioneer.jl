# Plan: Complete Removal of Quantile Binning (qbin) Features

## Overview
This plan outlines the complete removal of the quantile binning feature system from Pioneer.jl. The most recent commit (0899293c) already removed some qbin functionality, but this plan addresses all remaining code and parameters.

## Background
Quantile binning was originally implemented to discretize continuous features (prec_mz, irt_pred, weight, tic) into bins for machine learning models. After evaluation, the decision was made to remove this feature entirely and use only the original continuous features.

## Scope of Changes

### 1. Core Feature Generation Code
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl`

#### Remove Function (Lines 329-436)
- **Function**: `add_quantile_binned_features!(df::DataFrame, features::Vector{Symbol}, n_bins::Int=100)`
- **Purpose**: Creates quantile-binned versions of features with `_qbin` suffix
- **Action**: Delete entire function including docstring

#### Clean Up Feature Set Definitions
- **ADVANCED_FEATURE_SET** (Lines 44-105):
  - Already has qbin features commented out (lines 47, 51, 77, 80)
  - **Action**: Remove commented-out lines entirely

- **REDUCED_FEATURE_SET** (Lines 107-131):
  - Currently USES qbin features: `:prec_mz_qbin`, `:irt_pred_qbin`, `:weight_qbin`, `:tic_qbin` (lines 109, 111, 121)
  - **Action**: Replace with continuous versions: `:prec_mz`, `:irt_pred`, `:weight`, `:tic`

### 2. Feature Generation Calls
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

#### Remove Binning Calls (3 locations)
1. **Lines 83-86** - Out-of-memory path (DISABLED):
```julia
# Add quantile-binned features before training
features_to_bin = [:prec_mz, :irt_pred, :weight, :tic]
@user_info "Creating quantile-binned features with $n_quantile_bins bins: $(join(string.(features_to_bin), ", "))"
add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)
```
**Action**: Delete these 4 lines

2. **Lines 103-106** - In-memory path (ACTIVE):
```julia
# Add quantile-binned features before training
features_to_bin = [:prec_mz, :irt_pred, :weight, :tic]
@user_info "Creating quantile-binned features with $n_quantile_bins bins: $(join(string.(features_to_bin), ", "))"
add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)
```
**Action**: Delete these 4 lines

3. **Line 69** - Function parameter:
```julia
n_quantile_bins::Int64,
```
**Action**: Remove this parameter from `score_precursor_isotope_traces()` function signature

### 3. Diagnostic Logging
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

#### Remove qbin Diagnostics (3 locations)
1. **Lines 364-371** - LightGBM training diagnostics:
```julia
# Diagnostic: Report which quantile-binned features are being used
qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
if !isempty(qbin_features)
    @user_info "LightGBM using $(length(qbin_features)) quantile-binned features: ..."
    ...
end
```
**Action**: Delete entire diagnostic block

2. **Lines 411-420** - Probit regression diagnostics:
```julia
# Diagnostic: Report which quantile-binned features are being used
qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
if !isempty(qbin_features)
    @user_info "Probit using $(length(qbin_features)) quantile-binned features: ..."
    ...
end
```
**Action**: Delete entire diagnostic block

3. **Lines 808-815** - OOM diagnostics (in commented-out function):
```julia
# Diagnostic: Report which quantile-binned features are being used
qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
if !isempty(qbin_features)
    @user_info "OOM LightGBM using $(length(qbin_features)) quantile-binned features: ..."
    ...
end
```
**Action**: Delete entire diagnostic block

### 4. Parameter Removal
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

#### Remove from ScoringSearchParameters Struct
- **Line 67**: Field definition `n_quantile_bins::Int64`
- **Line 109**: Constructor initialization `Int64(ml_params.n_quantile_bins)`
- **Line 397**: Parameter passing in `score_precursor_isotope_traces()` call

**Action**: Remove all three references

**File**: `assets/example_config/defaultSearchParams.json`

#### Remove Parameter Definition
- **Line 148**: `"n_quantile_bins": 25,`

**Action**: Delete this line (note: may need to handle trailing comma on previous line)

**File**: `src/Routines/SearchDIA/ParseInputs/paramsChecks.jl`

**Action**:
- Search for validation logic for `n_quantile_bins`
- Remove any validation checks for this parameter

### 5. Test Configuration Files
**Files**: All `*params*.json` files found by glob pattern

Potentially affected files:
- `test/test_config/ecoli_test_params.json`
- `data/precompile/search_*.json`
- Various test scenario parameter files

**Action**:
- Search each file for `n_quantile_bins`
- Remove parameter if present
- Verify JSON remains valid after removal

## Implementation Order

1. **Phase 1 - Remove Parameter Passing**
   - Remove parameter from `score_psms.jl` function signature
   - Remove parameter from `ScoringSearch.jl` struct and calls
   - Remove from `defaultSearchParams.json`

2. **Phase 2 - Remove Feature Generation**
   - Remove calls to `add_quantile_binned_features!()`
   - Remove diagnostic logging for qbin features

3. **Phase 3 - Remove Core Function**
   - Delete `add_quantile_binned_features!()` function
   - Clean up feature set definitions (remove commented qbin features)
   - Replace qbin features with continuous versions in REDUCED_FEATURE_SET

4. **Phase 4 - Clean Up Test Files**
   - Remove parameter from all test configuration files
   - Verify parameter validation doesn't break

## Testing Strategy

### Minimal Integration Test
Run the E. coli test to ensure basic functionality:
```bash
julia --project=. -e 'using Pioneer; SearchDIA("test/test_config/ecoli_test_params.json")'
```

### Validation Checks
1. **Compile Check**: Verify all files load without errors
2. **Feature Sets**: Verify models use continuous features, not binned
3. **Parameter Parsing**: Ensure removing parameter doesn't break config loading
4. **Model Training**: Verify LightGBM and Probit models train successfully

### Expected Behavior Changes
- Model training should use continuous features (`:prec_mz`, `:irt_pred`, `:weight`, `:tic`)
- No `_qbin` columns should be created in PSM DataFrames
- No diagnostic messages about quantile binning
- Parameter validation should not check for `n_quantile_bins`

## Rollback Plan
If issues arise:
1. Git branch: `fix-mbr-column-naming` has all original code
2. Commit hash: 0899293c already removed some qbin functionality
3. Can selectively revert individual files if needed

## Notes
- The out-of-memory processing path (lines 77-98 in score_psms.jl) is already DISABLED (hardcoded `if false`), so removing qbin calls there is safe
- REDUCED_FEATURE_SET currently depends on qbin features - must replace with continuous versions
- ADVANCED_FEATURE_SET already has qbin features commented out - good indicator of intent

## Success Criteria
✅ No references to `qbin` in active code (only in comments if desired)
✅ No references to `n_quantile_bins` in code or config files
✅ `add_quantile_binned_features!()` function completely removed
✅ Feature sets use only continuous features
✅ All tests pass
✅ No breaking changes to user-facing API
