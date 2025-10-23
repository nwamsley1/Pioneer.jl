# Plan: Replace Global Q-value Spline with Dictionary Lookup

## Problem Statement

Currently, global q-values are calculated using a spline interpolation approach:
1. Merge all PSMs sorted by `global_prob`
2. Calculate q-values using `get_qvalue_spline()` which creates an interpolation function
3. Map q-values back to PSMs by interpolating from their `global_prob` values

However, `global_prob` is calculated at the **precursor level** (one value per precursor_idx across all runs), so each precursor should have exactly **one** global q-value. Using spline interpolation is unnecessary and potentially introduces imprecision.

## Proposed Solution

Replace spline interpolation with a direct dictionary mapping:
1. After calculating `global_prob`, reduce to one row per precursor
2. Calculate q-values on this reduced dataset
3. Create `Dict{UInt32, Float32}` mapping `precursor_idx -> global_qval`
4. Replace `add_interpolated_column` with dictionary lookup operation

## Benefits

1. **Exactness**: No interpolation error - each precursor gets its exact q-value
2. **Efficiency**: Dictionary lookup is O(1) vs spline evaluation overhead
3. **Clarity**: Makes it explicit that global q-values are per-precursor, not per-PSM
4. **Memory**: Dictionary is more compact than spline representation
5. **Correctness**: Ensures all PSMs from same precursor get identical global q-value

## Current Implementation

### Step 5: Calculate global_prob (ScoringSearch.jl, lines 373-412)
```julia
# Calculate global_prob by grouping by precursor_idx
transform!(groupby(merged_df, :precursor_idx),
           :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
```
**Result**: Each precursor_idx has one `global_prob` value, but PSMs still have multiple rows per precursor (one per run/file).

### Step 6: Merge by global_prob (ScoringSearch.jl, lines 414-419)
```julia
stream_sorted_merge(filtered_refs, results.merged_quant_path, :MBR_boosted_global_prob, :target;
                   batch_size=10_000_000, reverse=[true,true])
```
**Result**: All PSMs merged, sorted by global_prob. Multiple rows per precursor.

### Step 7: Calculate global q-values via spline (ScoringSearch.jl, lines 421-426)
```julia
results.precursor_global_qval_interp[] = get_precursor_global_qval_spline(results.merged_quant_path, params, search_context)
```
**Result**: Spline interpolation function stored in `results.precursor_global_qval_interp[]`.

### Step 10: Add global q-values via interpolation (ScoringSearch.jl, lines 457, 466)
```julia
add_interpolated_column(:MBR_boosted_global_qval, :MBR_boosted_global_prob, results.precursor_global_qval_interp[])
add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[])
```
**Result**: Q-values interpolated from global_prob for each PSM row.

## Proposed Implementation

### New Function: get_precursor_global_qval_dict (ScoringSearch.jl)

```julia
"""
    get_precursor_global_qval_dict(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    -> Dict{UInt32, Float32}

Calculate global q-values at the precursor level and return as dictionary mapping.

# Process
1. Load merged PSM data
2. Reduce to one row per precursor (group by precursor_idx, take first/max global_prob)
3. Calculate q-values on reduced dataset
4. Return Dict{precursor_idx => global_qval}
"""
function get_precursor_global_qval_dict(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    # Determine which score column to use
    score_col = params.match_between_runs ? :MBR_boosted_global_prob : :global_prob

    # Load merged PSMs
    merged_table = Arrow.Table(merged_path)
    df = DataFrame(merged_table)

    # Reduce to one row per precursor
    # Group by precursor_idx and keep the row with maximum global_prob (they should all be the same)
    precursor_df = combine(groupby(df, :precursor_idx)) do group
        # All rows for same precursor should have same global_prob
        # Just take the first one along with its target label
        (global_prob = first(group[!, score_col]),
         target = first(group.target))
    end

    # Sort by global_prob descending
    sort!(precursor_df, :global_prob, rev=true)

    # Calculate q-values using standard FDR calculation
    n = nrow(precursor_df)
    qvals = Vector{Float32}(undef, n)

    # Use library FDR scale factor
    fdr_scale = getLibraryFdrScaleFactor(search_context)

    # Calculate q-values (same logic as get_qvalues! but for precursor level)
    get_qvalues!(precursor_df.global_prob, precursor_df.target, qvals; fdr_scale_factor=fdr_scale)

    # Create dictionary mapping precursor_idx -> global_qval
    qval_dict = Dict{UInt32, Float32}()
    for i in 1:n
        qval_dict[precursor_df.precursor_idx[i]] = qvals[i]
    end

    return qval_dict
end
```

### New Pipeline Operation: add_dict_column (PipelineOperations.jl)

```julia
"""
    add_dict_column(new_col::Symbol, key_col::Symbol, lookup_dict::Dict{K,V}) where {K,V}

Add a new column by looking up values in a dictionary based on a key column.

# Arguments
- `new_col`: Name of the new column to create
- `key_col`: Name of the column containing keys for dictionary lookup
- `lookup_dict`: Dictionary mapping keys to values

# Example
```julia
pipeline = TransformPipeline() |>
    add_dict_column(:global_qval, :precursor_idx, precursor_qval_dict)
```
"""
function add_dict_column(new_col::Symbol, key_col::Symbol, lookup_dict::Dict{K,V}) where {K,V}
    return AddDictColumn(new_col, key_col, lookup_dict)
end

struct AddDictColumn{K,V} <: PipelineOp
    new_col::Symbol
    key_col::Symbol
    lookup_dict::Dict{K,V}
end

function apply_op(op::AddDictColumn, df::DataFrame)
    # Create new column by looking up each key
    df[!, op.new_col] = [get(op.lookup_dict, key, missing) for key in df[!, op.key_col]]
    return df
end

function describe_op(op::AddDictColumn)
    return "add_dict_column($(op.new_col) from $(op.key_col))"
end
```

### Modified Step 7 (ScoringSearch.jl, lines 421-426)

**Before:**
```julia
# Step 7: Calculate global precursor q-values
step7_time = @elapsed begin
    results.precursor_global_qval_interp[] = get_precursor_global_qval_spline(results.merged_quant_path, params, search_context)
end
```

**After:**
```julia
# Step 7: Calculate global precursor q-values
step7_time = @elapsed begin
    results.precursor_global_qval_dict[] = get_precursor_global_qval_dict(results.merged_quant_path, params, search_context)
end
```

### Modified Step 10 (ScoringSearch.jl, lines 453-473)

**Before:**
```julia
if params.match_between_runs
    qvalue_filter_pipeline = TransformPipeline() |>
        add_interpolated_column(:MBR_boosted_global_qval, :MBR_boosted_global_prob, results.precursor_global_qval_interp[]) |>
        add_interpolated_column(:MBR_boosted_qval, :MBR_boosted_prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :MBR_boosted_prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([...])
else
    qvalue_filter_pipeline = TransformPipeline() |>
        add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
        add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([...])
end
```

**After:**
```julia
if params.match_between_runs
    qvalue_filter_pipeline = TransformPipeline() |>
        add_dict_column(:MBR_boosted_global_qval, :precursor_idx, results.precursor_global_qval_dict[]) |>
        add_interpolated_column(:MBR_boosted_qval, :MBR_boosted_prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :MBR_boosted_prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([...])
else
    qvalue_filter_pipeline = TransformPipeline() |>
        add_dict_column(:global_qval, :precursor_idx, results.precursor_global_qval_dict[]) |>
        add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([...])
end
```

### Modified ScoringSearchResults struct (ScoringSearch.jl)

**Before:**
```julia
mutable struct ScoringSearchResults <: SearchResults
    merged_quant_path::String
    precursor_qval_interp::Ref{Union{Nothing, Any}}
    precursor_pep_interp::Ref{Union{Nothing, Any}}
    precursor_global_qval_interp::Ref{Union{Nothing, Any}}
end
```

**After:**
```julia
mutable struct ScoringSearchResults <: SearchResults
    merged_quant_path::String
    precursor_qval_interp::Ref{Union{Nothing, Any}}
    precursor_pep_interp::Ref{Union{Nothing, Any}}
    precursor_global_qval_dict::Ref{Dict{UInt32, Float32}}
end
```

### Modified Step 11: Re-calculate q-values (ScoringSearch.jl, lines 487-595)

The recalculation logic also needs to use dictionary instead of spline:

**Lines 519-542 (MBR global qval recalculation):**

**Before:**
```julia
# Calculate global q-value spline using MBR-boosted global scores
global_qval_interp = get_qvalue_spline(
    results.merged_quant_path, :MBR_boosted_global_prob, true;
    min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
    fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
)

# Recalculate MBR_boosted_global_qval column from MBR-boosted scores
recalculate_global_qvalue_pipeline = TransformPipeline() |>
    add_interpolated_column(:MBR_boosted_global_qval, :MBR_boosted_global_prob, global_qval_interp)
```

**After:**
```julia
# Calculate global q-value dictionary using MBR-boosted global scores
global_qval_dict = get_precursor_global_qval_dict(results.merged_quant_path, params, search_context)

# Recalculate MBR_boosted_global_qval column via dictionary lookup
recalculate_global_qvalue_pipeline = TransformPipeline() |>
    add_dict_column(:MBR_boosted_global_qval, :precursor_idx, global_qval_dict)
```

**Similar changes for non-MBR path (lines 564-582).**

## Files to Modify

### 1. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`
- Add `get_precursor_global_qval_dict()` function
- Modify `ScoringSearchResults` struct to use `precursor_global_qval_dict` instead of `precursor_global_qval_interp`
- Update Step 7 to call new function
- Update Step 10 to use `add_dict_column` instead of `add_interpolated_column` for global q-values
- Update Step 11 recalculation logic to use dictionary

### 2. `src/utils/FileOperations/pipeline/PipelineOperations.jl`
- Add `add_dict_column()` function
- Add `AddDictColumn` struct
- Implement `apply_op()` for `AddDictColumn`
- Implement `describe_op()` for `AddDictColumn`
- Export `add_dict_column`

### 3. Optional: Remove obsolete function
- `get_precursor_global_qval_spline()` can be marked as deprecated or removed

## Implementation Steps

1. **Add pipeline operation** (PipelineOperations.jl):
   - Implement `add_dict_column()` and supporting structs
   - Test with simple example

2. **Add dictionary calculation** (ScoringSearch.jl):
   - Implement `get_precursor_global_qval_dict()`
   - Test that it produces correct q-values

3. **Modify results struct** (ScoringSearch.jl):
   - Change `precursor_global_qval_interp` to `precursor_global_qval_dict`
   - Update initialization

4. **Update Step 7** (ScoringSearch.jl):
   - Call new dictionary function instead of spline function

5. **Update Step 10** (ScoringSearch.jl):
   - Replace `add_interpolated_column` with `add_dict_column` for global q-values
   - Keep interpolation for prec_prob-based columns (qval, pep)

6. **Update Step 11** (ScoringSearch.jl):
   - Replace spline with dictionary in recalculation logic

7. **Testing**:
   - Run full pipeline on test dataset
   - Verify global q-values are identical across PSMs from same precursor
   - Compare results with old spline method (should be nearly identical)

## Testing Strategy

### Unit Test
```julia
@testset "Global Q-value Dictionary" begin
    # Create test data with multiple PSMs per precursor
    df = DataFrame(
        precursor_idx = [1, 1, 1, 2, 2, 3],
        global_prob = [0.9, 0.9, 0.9, 0.7, 0.7, 0.5],
        target = [true, true, true, false, false, true]
    )

    # Calculate q-values
    qval_dict = get_precursor_global_qval_dict(df, ...)

    # Verify one entry per precursor
    @test length(qval_dict) == 3
    @test haskey(qval_dict, 1)
    @test haskey(qval_dict, 2)
    @test haskey(qval_dict, 3)
end
```

### Integration Test
```julia
# Run full pipeline
SearchDIA("test_params.json")

# Load results
psms = Arrow.Table("output/passing_psms/file1.arrow") |> DataFrame

# Verify all PSMs from same precursor have same global_qval
for (prec_idx, group) in pairs(groupby(psms, :precursor_idx))
    @test allequal(group.global_qval)
end
```

## Backward Compatibility

This change is **not** backward compatible with existing code that expects spline interpolation objects. However:
- The change is internal to ScoringSearch
- No external API changes
- Only affects ScoringSearchResults struct
- Old result objects cannot be reused (not typically done anyway)

## Performance Considerations

**Expected improvements:**
- Dictionary lookup: O(1) vs spline evaluation: O(log n) or O(1) with interpolation overhead
- Memory: Dict is more compact than spline representation
- Clarity: Code is simpler and more explicit

**No performance concerns:**
- Dictionary size is number of precursors (typically 10K-100K)
- Lookup is very fast
- No regression expected

## Edge Cases to Handle

1. **Missing precursor_idx**: Dictionary lookup returns `missing`
2. **No targets/decoys**: Q-value calculation may fail - handle gracefully
3. **Empty datasets**: Return empty dictionary
4. **Duplicate precursor entries**: Should not occur after reduction, but verify

## Alternative Considered

**Keep spline but reduce to precursor level first**: This would still use interpolation but on reduced data. Rejected because:
- Dictionary is simpler and more direct
- No benefit to interpolation when we have exact values
- Dictionary makes the per-precursor nature explicit
