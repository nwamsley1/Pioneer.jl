# Filter PSMs by Q-value Pipeline API Plan

## Status: ✅ COMPLETED (2025-01)

## Overview

This document outlines the plan to replace the monolithic `filter_psms_by_qvalue` function with a composable Pipeline API, following the successful pattern established for `sort_and_filter_quant_tables_refs`.

## Problem Statement

Current issues with `filter_psms_by_qvalue`:
- Combines multiple operations (add columns, filter, write) in one function
- Hard-coded column names and filtering logic
- Not reusable for similar filtering workflows
- Operations are hidden from the calling code
- Difficult to test individual components

## Solution: Composable Pipeline Operations

### New Pipeline Operations Needed

1. **`add_interpolated_column(new_col, source_col, interpolator)`**
   - Adds a new column by applying an interpolation function to an existing column
   - Example: `add_interpolated_column(:global_qval, :global_prob, global_qval_interp)`

2. **`filter_by_threshold(col, threshold; comparison = :<=)`**
   - Filter rows where column meets threshold condition
   - Example: `filter_by_threshold(:global_qval, 0.01)`

3. **`filter_by_multiple_thresholds(conditions)`**
   - Apply multiple threshold filters with AND logic
   - Example: `filter_by_multiple_thresholds([(:global_qval, 0.01), (:qval, 0.01)])`

4. **`write_to_folder(output_folder; preserve_basename = true)`**
   - Terminal operation that writes filtered results to new location
   - Preserves original filename structure

### Implementation Example

Replace current usage:
```julia
passing_refs = filter_psms_by_qvalue(
    filtered_refs,
    passing_psms_folder,
    getPrecursors(getSpecLib(search_context)),
    results.precursor_global_qval_interp[],
    results.precursor_qval_interp[],
    params.q_value_threshold
)
```

With explicit pipeline:
```julia
# Build q-value filtering pipeline
qvalue_filter_pipeline = TransformPipeline() |>
    add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
    add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
    filter_by_multiple_thresholds([
        (:global_qval, params.q_value_threshold),
        (:qval, params.q_value_threshold)
    ])

# Apply pipeline and write to new location
passing_refs = PSMFileReference[]
for ref in filtered_refs
    if exists(ref)
        # Apply transformations
        apply_pipeline!(ref, qvalue_filter_pipeline)
        
        # Write to new location
        output_path = joinpath(passing_psms_folder, basename(file_path(ref)))
        new_ref = write_transformed(ref, output_path)
        push!(passing_refs, new_ref)
    end
end
```

## Implementation Phases

### Phase 1: Add New Operations to FileOperations.jl
- Implement `add_interpolated_column`
- Implement `filter_by_threshold` and `filter_by_multiple_thresholds`
- Add tests for each operation

### Phase 2: Add Write Operations
- Implement `write_transformed` function that creates new FileReference
- Ensure proper metadata tracking (sort state, schema)
- Handle file path management

### Phase 3: Update ScoringSearch.jl
- Replace `filter_psms_by_qvalue` call with explicit pipeline
- Ensure passing_refs are properly collected
- Test with real data

### Phase 4: Deprecate Old Function
- Add deprecation warning to `filter_psms_by_qvalue`
- Provide migration guide in warning message
- Keep function for backward compatibility

## Benefits

1. **Transparency**: Q-value calculation and filtering logic is explicit
2. **Flexibility**: Easy to adjust thresholds or add additional filters
3. **Reusability**: Operations can be used in other contexts
4. **Testability**: Each operation can be unit tested
5. **Maintainability**: Clear what transformations are applied

## Alternative Design Considerations

### Option 1: Combined Filter-and-Write Operation
```julia
# Single operation that filters and writes
filtered_refs = filter_and_write_by_qvalue(
    filtered_refs,
    passing_psms_folder,
    qvalue_filter_pipeline
)
```

### Option 2: Pipeline with Terminal Operation
```julia
# Pipeline includes write as terminal operation
qvalue_pipeline = TransformPipeline() |>
    add_interpolated_column(...) |>
    filter_by_threshold(...) |>
    write_to_folder(passing_psms_folder)

# Returns new references
passing_refs = apply_pipeline_batch(filtered_refs, qvalue_pipeline)
```

### Recommendation: Option 2
- Keeps pipeline composable
- Write operation is explicit
- Batch processing is clear

## Testing Strategy

1. Unit tests for interpolation operations
2. Tests for threshold filtering with various conditions
3. Tests for write operations and reference creation
4. Integration test comparing old vs new implementation
5. Performance benchmarks

## Future Extensions

- Support for OR logic in filters
- Parallel pipeline application
- Streaming write operations for large files
- Custom interpolation functions

## Migration Path

1. Implement new operations alongside existing code
2. Create comparison tests to ensure identical results
3. Update ScoringSearch to use new API
4. Monitor performance in production
5. Deprecate old function after validation period

## Example: Other Use Cases

```julia
# Filter by peptide length
length_filter_pipeline = TransformPipeline() |>
    add_column(:peptide_length, row -> length(row.sequence)) |>
    filter_by_threshold(:peptide_length, 7, comparison = :>=) |>
    filter_by_threshold(:peptide_length, 50, comparison = :<=)

# Complex scoring filters
score_filter_pipeline = TransformPipeline() |>
    add_column(:combined_score, row -> row.score1 * 0.7 + row.score2 * 0.3) |>
    filter_by_threshold(:combined_score, 0.95) |>
    sort_by([:combined_score], rev=[true])
```

This approach continues the pattern of making file operations explicit, composable, and testable throughout the SearchDIA pipeline.