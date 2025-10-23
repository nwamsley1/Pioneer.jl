# Plan: Fix apply_mbr_filter to use trace_prob and MBR_boosted_trace_prob

## Problem
The `apply_mbr_filter` function and its callees use `:prob` in several places, but should use `:trace_prob` for finalized predictions. Additionally, the conditional check should verify column existence rather than using `params.match_between_runs`.

## Files to Modify

### 1. ScoringSearch.jl (Line 310-312)
**Current Code:**
```julia
if params.match_between_runs
    # Apply MBR filter to MBR_boosted_trace_prob column (modifies in place)
    apply_mbr_filter!(merged_df, params)
```

**Issue:** Should check if `MBR_boosted_trace_prob` column exists instead of checking `params.match_between_runs`.

**Fix:**
```julia
if hasproperty(merged_df, :MBR_boosted_trace_prob)
    # Apply MBR filter to MBR_boosted_trace_prob column (modifies in place)
    apply_mbr_filter!(merged_df, params)
```

### 2. scoring_interface.jl - Multiple Functions

#### Line 161: `train_and_evaluate(ThresholdFilter)`
**Current Code:**
```julia
if isempty(candidate_data) || !hasproperty(candidate_data, :prob)
    return nothing
end

# candidate_labels represents bad transfer flags
τ = get_ftr_threshold(
    candidate_data.prob,
    candidate_labels,
    params.max_MBR_false_transfer_rate
)

# Handle edge case where threshold is infinite (no valid threshold found)
if isinf(τ)
    n_passing = 0
else
    n_passing = sum(candidate_data.prob .>= τ)
end

return FilterResult("Threshold", candidate_data.prob, τ, n_passing)
```

**Issue:** References `.prob` instead of `.trace_prob`.

**Fix:**
```julia
if isempty(candidate_data) || !hasproperty(candidate_data, :trace_prob)
    return nothing
end

# candidate_labels represents bad transfer flags
τ = get_ftr_threshold(
    candidate_data.trace_prob,
    candidate_labels,
    params.max_MBR_false_transfer_rate
)

# Handle edge case where threshold is infinite (no valid threshold found)
if isinf(τ)
    n_passing = 0
else
    n_passing = sum(candidate_data.trace_prob .>= τ)
end

return FilterResult("Threshold", candidate_data.trace_prob, τ, n_passing)
```

#### Line 426: `select_mbr_features()`
**Current Code:**
```julia
candidate_features = [
    :prob,
    :irt_error, :ms1_ms2_rt_diff, :MBR_max_pair_prob, :MBR_best_irt_diff,
    :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio
]
```

**Issue:** Uses `:prob` instead of `:trace_prob`.

**Fix:**
```julia
candidate_features = [
    :trace_prob,
    :irt_error, :ms1_ms2_rt_diff, :MBR_max_pair_prob, :MBR_best_irt_diff,
    :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio
]
```

#### Line 583: `add_best_trace_indicator()`
**Current Code:**
```julia
# Group-based operation for combined traces
transform!(groupby(df, :precursor_idx),
          :prob => (p -> begin
              best_idx = argmax(p)
              result = falses(length(p))
              result[best_idx] = true
              result
          end) => :best_trace)
```

**Issue:** Uses `:prob` instead of `:trace_prob`.

**Fix:**
```julia
# Group-based operation for combined traces
transform!(groupby(df, :precursor_idx),
          :trace_prob => (p -> begin
              best_idx = argmax(p)
              result = falses(length(p))
              result[best_idx] = true
              result
          end) => :best_trace)
```

## Changes Summary

| File | Line(s) | Change |
|------|---------|--------|
| ScoringSearch.jl | 310 | Change `if params.match_between_runs` to `if hasproperty(merged_df, :MBR_boosted_trace_prob)` |
| scoring_interface.jl | 161 | Change `!hasproperty(candidate_data, :prob)` to `!hasproperty(candidate_data, :trace_prob)` |
| scoring_interface.jl | 167 | Change `candidate_data.prob` to `candidate_data.trace_prob` |
| scoring_interface.jl | 176 | Change `candidate_data.prob` to `candidate_data.trace_prob` |
| scoring_interface.jl | 179 | Change `candidate_data.prob` to `candidate_data.trace_prob` |
| scoring_interface.jl | 426 | Change `:prob` to `:trace_prob` in candidate_features |
| scoring_interface.jl | 583 | Change `:prob` to `:trace_prob` in transform! |

## Rationale

1. **Column Existence Check**: Using `hasproperty()` is more robust than checking parameter flags, as it directly verifies the data state.

2. **Consistent Column Usage**: The code writes final predictions to `:trace_prob` and `:MBR_boosted_trace_prob`, so all downstream operations should reference these columns, not the temporary `:prob` column.

3. **MBR Filter Scope**: `apply_mbr_filter` operates on finalized predictions after training completes, so it should use `:trace_prob` as input features.

4. **Best Trace Selection**: The best trace indicator should be based on final trace probabilities (`:trace_prob`), not temporary working probabilities.

## Testing

After implementing these changes:
1. Run with MBR enabled to verify filtering works correctly
2. Run without MBR to verify non-MBR path still works
3. Verify that output is sorted by `MBR_boosted_prec_prob` (descending)
4. Check that `MBR_candidate` and `MBR_transfer_q_value` columns are correctly populated
