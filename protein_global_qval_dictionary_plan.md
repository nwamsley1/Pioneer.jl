# Plan: Replace Protein Group Global Q-value Spline with Dictionary Lookup

## Problem Statement

Currently, protein group global q-values are calculated using spline interpolation:
1. Calculate `global_pg_score` by taking max `pg_score` across all files for each protein
2. Merge all protein groups sorted by `global_pg_score`
3. Calculate q-values using `get_qvalue_spline()` which creates an interpolation function
4. Map q-values back to protein groups by interpolating from their `global_pg_score` values

However, `global_pg_score` is calculated at the **protein group level** (one value per (protein_name, target, entrap_id) tuple across all runs), so each protein group should have exactly **one** global q-value. Using spline interpolation is unnecessary and potentially introduces imprecision.

## Proposed Solution

Replace spline interpolation with a direct dictionary mapping:
1. After calculating `global_pg_score`, reduce to one row per protein group
2. Calculate q-values on this reduced dataset
3. Create dictionary mapping `(protein_name, target, entrap_id) -> global_pg_qval`
4. Replace `add_interpolated_column` with dictionary lookup operation

## Benefits

1. **Exactness**: No interpolation error - each protein group gets its exact q-value
2. **Efficiency**: Dictionary lookup is O(1) vs spline evaluation overhead
3. **Clarity**: Makes it explicit that global q-values are per-protein-group
4. **Correctness**: Ensures all instances of same protein group get identical global q-value
5. **Consistency**: Matches the approach used for precursor global q-values

## Current Implementation

### Step 16: Calculate global_pg_score (ScoringSearch.jl, lines 686-689)
```julia
acc_to_max_pg_score = calculate_and_add_global_scores!(pg_refs)
```
**Result**: Each protein group (protein_name, target, entrap_id) has one `global_pg_score`, but protein group tables still have multiple rows per protein (one per file).

### Step 17-18: Sort and merge by global_pg_score (lines 691-701)
```julia
sort_file_by_keys!(pg_refs, :global_pg_score, :target; reverse=[true, true])
stream_sorted_merge(pg_refs, sorted_pg_scores_path, :global_pg_score, :target;
                   batch_size=1000000, reverse=[true,true])
```
**Result**: All protein groups merged, sorted by global_pg_score. Multiple rows per protein group.

### Step 19: Calculate global q-values via spline (lines 703-706)
```julia
search_context.global_pg_score_to_qval[] = get_protein_global_qval_spline(sorted_pg_scores_path, params)
```
**Result**: Spline interpolation function stored in `search_context.global_pg_score_to_qval[]`.

### Step 23: Add global q-values via interpolation (line 728)
```julia
add_interpolated_column(:global_pg_qval, :global_pg_score, search_context.global_pg_score_to_qval[])
```
**Result**: Q-values interpolated from global_pg_score for each protein group row.

## Proposed Implementation

### Dictionary Key Type

Since protein groups are identified by (protein_name, target, entrap_id), we need a composite key:

```julia
# Define key type for protein groups
struct ProteinGroupKey
    protein_name::String
    target::Bool
    entrap_id::UInt8
end

# Make it hashable for dictionary use
Base.hash(pg::ProteinGroupKey, h::UInt) = hash((pg.protein_name, pg.target, pg.entrap_id), h)
Base.isequal(a::ProteinGroupKey, b::ProteinGroupKey) =
    a.protein_name == b.protein_name && a.target == b.target && a.entrap_id == b.entrap_id
```

### New Function: get_protein_global_qval_dict (ScoringSearch.jl)

```julia
"""
    get_protein_global_qval_dict(merged_path::String, params::ScoringSearchParameters)
    -> Dict{ProteinGroupKey, Float32}

Calculate global q-values at the protein group level and return as dictionary mapping.

Since global_pg_score is calculated per protein group (one value per (protein_name, target, entrap_id)
across all files), each protein group should have exactly one global q-value. This function:
1. Loads merged protein group data
2. Reduces to one row per protein group
3. Calculates q-values on the reduced dataset
4. Returns Dict{ProteinGroupKey => global_pg_qval}

This approach is more accurate than spline interpolation since we have exact values for each protein group.
"""
function get_protein_global_qval_dict(merged_path::String, params::ScoringSearchParameters)
    # Load merged protein groups
    merged_table = Arrow.Table(merged_path)
    df = DataFrame(merged_table)

    # Reduce to one row per protein group
    # All rows for the same protein group should have the same global_pg_score
    protein_group_df = combine(groupby(df, [:protein_name, :target, :entrap_id])) do group
        (global_pg_score = first(group.global_pg_score),)
    end

    # Sort by global_pg_score descending for q-value calculation
    sort!(protein_group_df, :global_pg_score, rev=true)

    # Calculate q-values using standard FDR calculation
    n = nrow(protein_group_df)
    qvals = Vector{Float32}(undef, n)

    # Calculate q-values
    get_qvalues!(protein_group_df.global_pg_score, protein_group_df.target, qvals)

    # Create dictionary mapping protein group key -> global_pg_qval
    qval_dict = Dict{ProteinGroupKey, Float32}()
    for i in 1:n
        key = ProteinGroupKey(
            protein_group_df.protein_name[i],
            protein_group_df.target[i],
            protein_group_df.entrap_id[i]
        )
        qval_dict[key] = qvals[i]
    end

    return qval_dict
end
```

### Modified Pipeline Operation: add_dict_column_composite_key

We need a variant of `add_dict_column` that can handle composite keys:

```julia
"""
    add_dict_column_composite_key(new_col::Symbol, key_cols::Vector{Symbol}, lookup_dict::Dict{ProteinGroupKey,V}) where {V}

Add a new column by looking up values in a dictionary using multiple key columns to form a composite key.

# Arguments
- `new_col`: Name of the new column to create
- `key_cols`: Vector of column names to use for creating composite keys (e.g., [:protein_name, :target, :entrap_id])
- `lookup_dict`: Dictionary mapping ProteinGroupKey to values

# Example
```julia
pipeline = TransformPipeline() |>
    add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], pg_qval_dict)
```
"""
function add_dict_column_composite_key(new_col::Symbol, key_cols::Vector{Symbol}, lookup_dict::Dict{ProteinGroupKey,V}) where {V}
    desc = "add_dict_column_composite_key($new_col from $(join(key_cols, ", ")))"
    op = function(df)
        # Pre-allocate result vector
        result = Vector{Union{V, Missing}}(undef, nrow(df))

        # Look up each composite key in dictionary
        for i in 1:nrow(df)
            key = ProteinGroupKey(
                df.protein_name[i],
                df.target[i],
                df.entrap_id[i]
            )
            result[i] = get(lookup_dict, key, missing)
        end

        df[!, new_col] = result
        return df
    end
    return desc => op
end
```

**Alternative (simpler):** Just use a tuple as the key instead of a custom struct:

```julia
# Use tuple keys: Dict{Tuple{String, Bool, UInt8}, Float32}

function get_protein_global_qval_dict(merged_path::String, params::ScoringSearchParameters)
    # ... load and reduce data ...

    # Create dictionary using tuple keys
    qval_dict = Dict{Tuple{String, Bool, UInt8}, Float32}()
    for i in 1:n
        key = (protein_group_df.protein_name[i],
               protein_group_df.target[i],
               protein_group_df.entrap_id[i])
        qval_dict[key] = qvals[i]
    end
    return qval_dict
end

# Then in pipeline operation:
function add_dict_column_composite_key(new_col::Symbol, key_cols::Vector{Symbol},
                                      lookup_dict::Dict{Tuple{String,Bool,UInt8},V}) where {V}
    desc = "add_dict_column_composite_key($new_col from $(join(key_cols, ", ")))"
    op = function(df)
        result = Vector{Union{V, Missing}}(undef, nrow(df))
        for i in 1:nrow(df)
            key = (df.protein_name[i], df.target[i], df.entrap_id[i])
            result[i] = get(lookup_dict, key, missing)
        end
        df[!, new_col] = result
        return df
    end
    return desc => op
end
```

**Recommendation:** Use the tuple approach for simplicity. It's more idiomatic Julia and doesn't require defining a new type.

### Modified Step 19 (ScoringSearch.jl, lines 703-706)

**Before:**
```julia
# Step 19: Calculate global protein q-values
step19_time = @elapsed begin
    search_context.global_pg_score_to_qval[] = get_protein_global_qval_spline(sorted_pg_scores_path, params)
end
```

**After:**
```julia
# Step 19: Calculate global protein q-values
step19_time = @elapsed begin
    search_context.global_pg_score_to_qval_dict[] = get_protein_global_qval_dict(sorted_pg_scores_path, params)
end
```

### Modified Step 23 (ScoringSearch.jl, lines 725-736)

**Before:**
```julia
# Step 23: Add q-values and passing flags to protein groups
step23_time = @elapsed begin
    protein_qval_pipeline = TransformPipeline() |>
        add_interpolated_column(:global_pg_qval, :global_pg_score, search_context.global_pg_score_to_qval[]) |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
        add_interpolated_column(:pg_pep, :pg_score, search_context.pg_score_to_pep[]) |>
        filter_by_multiple_thresholds([
            (:global_pg_qval, params.q_value_threshold),
            (:pg_qval, params.q_value_threshold)
        ])
    apply_pipeline!(pg_refs, protein_qval_pipeline)
end
```

**After:**
```julia
# Step 23: Add q-values and passing flags to protein groups
step23_time = @elapsed begin
    protein_qval_pipeline = TransformPipeline() |>
        add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], search_context.global_pg_score_to_qval_dict[]) |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
        add_interpolated_column(:pg_pep, :pg_score, search_context.pg_score_to_pep[]) |>
        filter_by_multiple_thresholds([
            (:global_pg_qval, params.q_value_threshold),
            (:pg_qval, params.q_value_threshold)
        ])
    apply_pipeline!(pg_refs, protein_qval_pipeline)
end
```

### Modified SearchContext fields

Need to check where `global_pg_score_to_qval` is defined in SearchContext and change its type from interpolation to dictionary.

**Before:**
```julia
global_pg_score_to_qval::Ref{Any}  # Interpolation object
```

**After:**
```julia
global_pg_score_to_qval_dict::Ref{Dict{Tuple{String,Bool,UInt8}, Float32}}  # Dictionary
```

### Update Step 24: update_psms_with_probit_scores_refs (lines 753-761)

This function receives `global_pg_score_to_qval` - need to check what it does with it and update accordingly.

**Current:**
```julia
update_psms_with_probit_scores_refs(
    paired_files,
    acc_to_max_pg_score,
    search_context.pg_score_to_qval[],
    search_context.global_pg_score_to_qval[]
)
```

Need to investigate this function to see how it uses the global q-value mapping.

## Files to Modify

### 1. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`
- Add `get_protein_global_qval_dict()` function
- Update Step 19 to call new dictionary function
- Update Step 23 to use `add_dict_column_composite_key`
- Update Step 24 if needed

### 2. `src/utils/FileOperations/pipeline/PipelineOperations.jl`
- Add `add_dict_column_composite_key()` function for composite key lookups
- Export the function

### 3. `src/Routines/SearchDIA/SearchContext.jl` (or wherever SearchContext is defined)
- Change `global_pg_score_to_qval` field from `Ref{Any}` to `Ref{Dict{Tuple{String,Bool,UInt8}, Float32}}`
- Update initialization

### 4. Check `update_psms_with_probit_scores_refs` function
- Verify how it uses `global_pg_score_to_qval`
- Update if needed to work with dictionary instead of interpolation

### 5. Optional: Remove obsolete function
- `get_protein_global_qval_spline()` can be marked as deprecated or removed

## Implementation Steps

1. **Add composite key pipeline operation** (PipelineOperations.jl):
   - Implement `add_dict_column_composite_key()` using tuple keys
   - Test with simple example

2. **Add dictionary calculation** (ScoringSearch.jl):
   - Implement `get_protein_global_qval_dict()`
   - Test that it produces correct q-values

3. **Modify SearchContext**:
   - Change `global_pg_score_to_qval` → `global_pg_score_to_qval_dict`
   - Update type and initialization

4. **Update Step 19** (ScoringSearch.jl):
   - Call dictionary function instead of spline function

5. **Update Step 23** (ScoringSearch.jl):
   - Replace `add_interpolated_column` with `add_dict_column_composite_key` for global q-values
   - Keep interpolation for pg_score-based columns

6. **Update Step 24** (if needed):
   - Check and update `update_psms_with_probit_scores_refs` if necessary

7. **Testing**:
   - Run full pipeline on test dataset
   - Verify global q-values are identical across all instances of same protein group
   - Compare results with old spline method (should be nearly identical)

## Testing Strategy

### Unit Test
```julia
@testset "Protein Group Global Q-value Dictionary" begin
    # Create test data with multiple instances per protein group
    df = DataFrame(
        protein_name = ["P1", "P1", "P1", "P2", "P2"],
        target = [true, true, true, false, false],
        entrap_id = [UInt8(0), UInt8(0), UInt8(0), UInt8(0), UInt8(0)],
        global_pg_score = [0.9, 0.9, 0.9, 0.7, 0.7]
    )

    # Calculate q-values
    qval_dict = get_protein_global_qval_dict(df, params)

    # Verify one entry per protein group
    @test length(qval_dict) == 2
    @test haskey(qval_dict, ("P1", true, UInt8(0)))
    @test haskey(qval_dict, ("P2", false, UInt8(0)))
end
```

### Integration Test
```julia
# Run full pipeline
SearchDIA("test_params.json")

# Load results
pg = Arrow.Table("output/passing_proteins/file1.arrow") |> DataFrame

# Verify all instances of same protein group have same global_pg_qval
for (pg_key, group) in pairs(groupby(pg, [:protein_name, :target, :entrap_id]))
    @test allequal(group.global_pg_qval)
end
```

## Edge Cases to Handle

1. **Missing protein groups**: Dictionary lookup returns `missing`
2. **No targets/decoys**: Q-value calculation may fail - handle gracefully
3. **Empty datasets**: Return empty dictionary
4. **Protein name encoding**: Ensure string comparisons work correctly with special characters

## Performance Considerations

**Expected improvements:**
- Dictionary lookup: O(1) vs spline evaluation overhead
- Memory: Dict is more compact than spline representation
- Clarity: Code is simpler and more explicit

**No performance concerns:**
- Dictionary size is number of unique protein groups (typically 1K-10K)
- Lookup is very fast
- No regression expected

## Backward Compatibility

This change is **not** backward compatible with existing code that expects spline interpolation objects. However:
- The change is internal to ScoringSearch
- SearchContext field name changes from `global_pg_score_to_qval` to `global_pg_score_to_qval_dict`
- Old SearchContext objects cannot be reused (not typically done anyway)

## Notes

- This mirrors the approach used for precursor global q-values
- Maintains consistency in how global q-values are handled across the pipeline
- Protein groups are more complex than precursors due to composite keys, but tuple keys handle this elegantly
