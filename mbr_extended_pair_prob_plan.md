# MBR Extended Pair Probability Features Implementation Plan

## Overview

This plan outlines the implementation of additional MBR (Match Between Runs) pair probability statistics alongside the existing `:MBR_max_pair_prob` feature. The goal is to add `:MBR_mean_pair_prob`, `:MBR_median_pair_prob`, and `:MBR_min_pair_prob` features to provide more comprehensive statistical information about matching precursor pairs across runs.

## Current Implementation Analysis

### Existing MBR_max_pair_prob Structure

The current implementation in `src/utils/ML/percolatorSortOf.jl` tracks only the top 2 best matching pairs per precursor:

```julia
@NamedTuple{
    best_prob_1::Float32,        # Highest probability match
    best_prob_2::Float32,        # Second highest probability match
    best_log2_weights_1::Vector{Float32},
    best_log2_weights_2::Vector{Float32},
    # ... other features for top 2 matches
}
```

Key locations where `MBR_max_pair_prob` is used:
- **Line 799**: `psms_subset.MBR_max_pair_prob[i] = scores.best_prob_1`
- **Line 807**: `psms_subset.MBR_max_pair_prob[i] = scores.best_prob_2`
- **Line 1064**: `sub_psms.MBR_max_pair_prob[i] = sub_psms.prob[best_idx]`

### Current Limitations

1. **Limited Statistics**: Only tracks the maximum (best) probability
2. **Two-Match Limit**: Only stores top 2 matches, insufficient for statistical measures
3. **Binary Choice**: Uses either `best_prob_1` or `best_prob_2` based on file index matching

## Implementation Strategy

### Phase 1: Extend Data Structures

#### 1.1 Modify Precursor Score Storage (Memory-Efficient Approach)

**File**: `src/utils/ML/percolatorSortOf.jl` around line 911

**Current**:
```julia
@NamedTuple{best_prob_1::Float32, best_prob_2::Float32, ...}
```

**Proposed (Memory-Efficient)**:
```julia
@NamedTuple{
    best_prob_1::Float32,        # Existing: Highest probability
    best_prob_2::Float32,        # Existing: Second highest probability
    worst_prob_1::Float32,       # NEW: Lowest probability
    worst_prob_2::Float32,       # NEW: Second lowest probability
    mean_prob::Float32,          # NEW: Running mean of all probabilities
    count_pairs::Int32,          # NEW: Count of valid matching pairs
    best_log2_weights_1::Vector{Float32},
    best_log2_weights_2::Vector{Float32},
    # ... existing fields
}
```

#### 1.2 Update DataFrame Initialization

**File**: `src/utils/ML/percolatorSortOf.jl` around line 1083

**Add new columns**:
```julia
if match_between_runs
    psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
    psms[!, :MBR_mean_pair_prob]            = zeros(Float32, n)  # NEW
    psms[!, :MBR_min_pair_prob]             = zeros(Float32, n)  # NEW
    psms[!, :MBR_num_pairs]                 = zeros(Int32, n)    # NEW: Count of valid pairs
    # ... existing columns
end
```

### Phase 2: Data Collection Enhancement

#### 2.1 Modify Score Collection Logic (Memory-Efficient)

**File**: `src/utils/ML/percolatorSortOf.jl` around line 694-740

**Current Logic**: Only tracks top 2 probabilities
**Proposed Logic**: Track running statistics (best 2, worst 2, mean, count)

```julia
# Helper function for updating statistics
function update_pair_statistics(current_stats, new_prob::Float32)
    # Update count and running mean
    new_count = current_stats.count_pairs + 1
    new_mean = (current_stats.mean_prob * current_stats.count_pairs + new_prob) / new_count

    # Update best probabilities (existing logic enhanced)
    new_best_1, new_best_2 = if new_prob > current_stats.best_prob_1
        (new_prob, current_stats.best_prob_1)
    elseif new_prob > current_stats.best_prob_2
        (current_stats.best_prob_1, new_prob)
    else
        (current_stats.best_prob_1, current_stats.best_prob_2)
    end

    # Update worst probabilities (new logic)
    new_worst_1, new_worst_2 = if new_count == 1
        (new_prob, zero(Float32))  # First probability
    elseif new_count == 2
        (min(current_stats.worst_prob_1, new_prob), max(current_stats.worst_prob_1, new_prob))
    elseif new_prob < current_stats.worst_prob_1
        (new_prob, current_stats.worst_prob_1)  # New minimum
    elseif new_prob < current_stats.worst_prob_2
        (current_stats.worst_prob_1, new_prob)  # New second minimum
    else
        (current_stats.worst_prob_1, current_stats.worst_prob_2)  # No change
    end

    return merge(current_stats, (
        best_prob_1 = new_best_1,
        best_prob_2 = new_best_2,
        worst_prob_1 = new_worst_1,
        worst_prob_2 = new_worst_2,
        mean_prob = new_mean,
        count_pairs = new_count
        # ... preserve existing fields
    ))
end

# In getBestScorePerPrec! function
if haskey(prec_to_best_score_new, key)
    scores = prec_to_best_score_new[key]
    new_scores = update_pair_statistics(scores, prob)
    prec_to_best_score_new[key] = new_scores
else
    # Initialize new entry
    insert!(prec_to_best_score_new, key, (
        best_prob_1 = prob,
        best_prob_2 = zero(Float32),
        worst_prob_1 = prob,         # NEW
        worst_prob_2 = zero(Float32), # NEW
        mean_prob = prob,            # NEW
        count_pairs = Int32(1),      # NEW
        # ... existing fields
    ))
end
```

#### 2.2 Statistical Computation Functions (Memory-Efficient)

**Statistics are computed incrementally during collection - no additional functions needed!**

The running statistics approach eliminates the need for separate statistical computation functions since all values are maintained incrementally:

- **Mean**: Computed using running average formula: `new_mean = (old_mean * old_count + new_value) / new_count`
- **Min**: Tracked as `worst_prob_1`
- **Max**: Tracked as `best_prob_1` (existing)
- **Count**: Tracked as `count_pairs`

This approach provides **O(1) time and space complexity** per probability update.

### Phase 3: Feature Assignment

#### 3.1 Update Feature Assignment Logic

**File**: `src/utils/ML/percolatorSortOf.jl` around line 799-813

**Current**:
```julia
psms_subset.MBR_max_pair_prob[i] = scores.best_prob_1  # or best_prob_2
```

**Proposed (Memory-Efficient)**:
```julia
# Existing max assignment
psms_subset.MBR_max_pair_prob[i] = scores.best_prob_1  # or best_prob_2

# NEW: Direct assignment from precomputed statistics
psms_subset.MBR_mean_pair_prob[i] = scores.mean_prob
psms_subset.MBR_min_pair_prob[i] = scores.worst_prob_1
psms_subset.MBR_num_pairs[i] = scores.count_pairs
```

#### 3.2 Handle Missing Data Cases

**File**: `src/utils/ML/percolatorSortOf.jl` around line 813, 1053

**Update missing data assignment**:
```julia
# When no valid matches exist
psms_subset.MBR_max_pair_prob[i]        = -1.0f0
psms_subset.MBR_mean_pair_prob[i]       = -1.0f0  # NEW
psms_subset.MBR_min_pair_prob[i]        = -1.0f0  # NEW
psms_subset.MBR_num_pairs[i]            = 0        # NEW
psms_subset.MBR_is_missing[i]           = true
```

### Phase 4: Feature Selection Integration

#### 4.1 Update Feature Selection

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl` around line 384

**Add new features to candidate list**:
```julia
candidate_features = [
    :prob, :irt_error, :irt_pred, :Mox, :tic, :rt_diff,
    :MBR_max_pair_prob, :MBR_best_irt_diff,
    :MBR_mean_pair_prob,      # NEW
    :MBR_min_pair_prob,       # NEW
    :MBR_num_pairs,           # NEW
    :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio, :MBR_num_runs
]
```

### Phase 5: Testing and Validation

#### 5.1 Unit Tests

**File**: `test/Routines/SearchDIA/SearchMethods/ScoringSearch/test_scoring_interface.jl`

**Add test cases**:
```julia
@testset "MBR Extended Pair Probability Features" begin
    # Test feature selection includes new features
    df = DataFrame(
        prob = [0.9, 0.8, 0.7, 0.6],
        MBR_max_pair_prob = [0.8, 0.9, 0.7, 0.6],
        MBR_mean_pair_prob = [0.7, 0.8, 0.6, 0.5],    # NEW
        MBR_median_pair_prob = [0.75, 0.85, 0.65, 0.55], # NEW
        MBR_min_pair_prob = [0.6, 0.7, 0.5, 0.4],     # NEW
        MBR_num_pairs = [3, 4, 2, 2],                  # NEW
        some_other_feature = [1, 2, 3, 4]
    )

    features = select_mbr_features(df)

    @test :MBR_mean_pair_prob in features
    @test :MBR_median_pair_prob in features
    @test :MBR_min_pair_prob in features
    @test :MBR_num_pairs in features
end

@testset "Pair Probability Statistics" begin
    # Test statistical computation function
    probs = [0.9f0, 0.8f0, 0.7f0, 0.6f0]
    mean_p, median_p, min_p, count_p = compute_pair_prob_statistics(probs)

    @test mean_p ≈ 0.75f0
    @test median_p ≈ 0.75f0
    @test min_p ≈ 0.6f0
    @test count_p ≈ 4.0f0

    # Test edge cases
    empty_mean, empty_median, empty_min, empty_count = compute_pair_prob_statistics(Float32[])
    @test empty_mean == -1.0f0
    @test empty_median == -1.0f0
    @test empty_min == -1.0f0
    @test empty_count == 0.0f0
end
```

#### 5.2 Integration Tests

**Update existing test cases** to include new features:
- Update mock data generation to include new columns
- Verify that new features are properly populated during MBR workflows
- Test that missing data cases are handled correctly

## Implementation Considerations

### Memory Usage

**Impact**: Adding `all_pair_probs::Vector{Float32}` will increase memory usage
**Mitigation**:
- Monitor memory usage with profiling during testing
- Consider limiting the maximum number of stored probabilities if memory becomes an issue
- Option to use `all_pair_probs::Union{Nothing, Vector{Float32}}` to save memory when disabled

### Backward Compatibility

**Approach**: All changes are additive, maintaining existing behavior
- Keep existing `MBR_max_pair_prob` logic unchanged
- New features have default values (-1.0f0) for missing data, consistent with existing patterns
- Feature selection remains optional based on column availability

### Performance Considerations

**Optimizations**:
- Compute statistics only when needed (non-empty probability vectors)
- Use efficient sorting algorithms for median computation
- Consider pre-allocating vectors to avoid repeated allocations

### Data Validation

**Quality Checks**:
- Ensure `MBR_max_pair_prob ≥ MBR_mean_pair_prob ≥ MBR_min_pair_prob` (when all are valid)
- Verify `MBR_num_pairs` matches the actual count of valid probabilities
- Add assertions to catch data inconsistencies during development

## File Modification Summary

| File | Lines | Modification Type | Description |
|------|-------|------------------|-------------|
| `src/utils/ML/percolatorSortOf.jl` | 911-925 | Struct Extension | Add `all_pair_probs` field to score storage |
| `src/utils/ML/percolatorSortOf.jl` | 1083-1089 | Column Addition | Initialize new DataFrame columns |
| `src/utils/ML/percolatorSortOf.jl` | 694-740 | Logic Enhancement | Collect all pair probabilities |
| `src/utils/ML/percolatorSortOf.jl` | 799-813 | Feature Assignment | Compute and assign new statistics |
| `src/utils/ML/percolatorSortOf.jl` | 1053 | Missing Data | Handle missing values for new features |
| `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl` | 384 | Feature Selection | Add new features to candidate list |
| `test/Routines/SearchDIA/SearchMethods/ScoringSearch/test_scoring_interface.jl` | Various | Test Addition | Add comprehensive test coverage |

## Implementation Timeline

1. **Phase 1-2** (Data Structures & Collection): ~2-3 days
2. **Phase 3** (Feature Assignment): ~1 day
3. **Phase 4** (Integration): ~0.5 days
4. **Phase 5** (Testing): ~1-2 days

**Total Estimated Effort**: 4.5-6.5 days

## Success Criteria

- [ ] All new features (`:MBR_mean_pair_prob`, `:MBR_median_pair_prob`, `:MBR_min_pair_prob`, `:MBR_num_pairs`) are available in MBR workflows
- [ ] Existing `:MBR_max_pair_prob` behavior remains unchanged
- [ ] New features are properly included in ML model training when available
- [ ] Comprehensive test coverage with edge case handling
- [ ] No significant performance regression in MBR processing
- [ ] Memory usage remains within acceptable bounds

This implementation will provide richer statistical information about cross-run precursor matching, potentially improving the accuracy of MBR filtering models by capturing the full distribution of matching quality rather than just the maximum value.