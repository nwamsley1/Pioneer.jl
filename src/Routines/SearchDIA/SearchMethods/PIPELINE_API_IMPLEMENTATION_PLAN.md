# Pipeline API Implementation Plan

## Status: ✅ COMPLETED (2025-01)

## Overview

This document outlines the implementation of a composable Pipeline API for FileReference operations. The goal is to replace monolithic functions like `sort_and_filter_quant_tables_refs` with clear, testable, and efficient pipeline operations.

## Problem Statement

Current issues with `sort_and_filter_quant_tables_refs`:
- Performs 6+ operations in one opaque function
- Hard to understand what transformations are applied
- Difficult to test individual operations
- Not reusable for other similar workflows
- Operations are hidden from the calling code

## Solution: Pipeline API

A composable API that:
1. Makes each operation explicit and named
2. Executes all operations in a single file pass
3. Properly tracks file metadata (sort state, schema)
4. Is testable at the operation level
5. Can be reused across different search methods

## Design

### Core Types

```julia
struct TransformPipeline
    operations::Vector{Pair{String, Function}}  # description => operation
    post_actions::Vector{Function}  # Actions after transform (e.g., mark_sorted!)
end

struct PipelineOperation
    operation::Pair{String, Function}
    post_action::Function
end
```

### Usage Example

```julia
# Clear, composable pipeline
pipeline = TransformPipeline() |>
    add_best_trace_indicator(isotope_type, best_traces) |>
    rename_column(:prob, :trace_prob) |>
    select_columns(necessary_cols) |>
    filter_rows(row -> row.best_trace; desc="keep_best_traces") |>
    remove_columns(:best_trace) |>
    sort_by([:global_prob, :target], rev=[true, true])

# Apply to files
for ref in file_refs
    apply_pipeline!(ref, pipeline)
end
```

## Implementation Phases

### Phase 1: Core Pipeline Infrastructure
- Add TransformPipeline type to FileOperations.jl
- Implement |> operator for chaining
- Create apply_pipeline! function
- Handle PipelineOperation with post-actions

### Phase 2: Operation Builders
Basic operations:
- `add_column(name, compute_fn)` - Add computed column
- `rename_column(old, new)` - Rename column
- `select_columns(cols)` - Keep only specified columns
- `remove_columns(cols...)` - Remove specified columns
- `filter_rows(predicate; desc)` - Filter rows by condition
- `sort_by(cols; rev)` - Sort with automatic state tracking

### Phase 3: Sort Integration
- Override `Base.sort!` for FileReference types
- Ensure sort state is always tracked
- sort_by operation includes post-action for mark_sorted!

### Phase 4: Specialized Operations
Scoring-specific:
- `add_best_trace_indicator(isotope_type, best_traces)`
- `get_quant_necessary_columns()` - Standard column set

### Phase 5: Update ScoringSearch.jl
Replace:
```julia
sort_and_filter_quant_tables_refs(
    second_pass_refs, 
    params.isotope_tracetype,
    :global_prob,
    best_traces 
)
```

With:
```julia
# Define columns we need
necessary_cols = get_quant_necessary_columns()

# Build explicit pipeline
quant_pipeline = TransformPipeline() |>
    add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
    rename_column(:prob, :trace_prob) |>
    select_columns(vcat(necessary_cols, :best_trace)) |>
    filter_rows(row -> row.best_trace; desc="keep_best_traces") |>
    remove_columns(:best_trace) |>
    sort_by([:global_prob, :target], rev=[true, true])

# Apply to all files
for ref in second_pass_refs
    if exists(ref)
        apply_pipeline!(ref, quant_pipeline)
    end
end
```

### Phase 6: Cleanup
- Deprecate sort_and_filter_quant_tables_refs
- Add documentation
- Update tests

## Benefits

1. **Clarity**: Each operation is explicit in the calling code
2. **Testability**: Individual operations can be unit tested
3. **Efficiency**: Single pass through file regardless of operation count
4. **Reusability**: Operations can be composed differently for other uses
5. **Maintainability**: Easy to modify or extend pipelines
6. **Safety**: Type system ensures proper operation composition

## Migration Strategy

1. Implement alongside existing code
2. Test thoroughly with real data
3. Update ScoringSearch to use new API
4. Deprecate old function with helpful error message
5. Update other search methods to use pipelines where beneficial

## Future Extensions

- Parallel pipeline execution for multiple files
- Pipeline optimization (operation reordering)
- Pipeline visualization/debugging tools
- Lazy evaluation for even better performance

## Code Organization

```
FileOperations.jl
├── Core Pipeline Types
├── Pipeline Application Logic
├── Basic Operation Builders
└── Utility Functions

scoring_interface.jl
├── Scoring-Specific Operations
├── Column Set Definitions
└── Deprecated Functions

ScoringSearch.jl
└── Direct pipeline usage (no wrapper functions)
```

## Testing Strategy

1. Unit tests for each operation builder
2. Tests for pipeline composition
3. Tests for sort state tracking
4. Integration tests with real Arrow files
5. Performance benchmarks vs. current implementation

## Example: Other Use Cases

```julia
# In IntegrateChromatogramsSearch
integration_pipeline = TransformPipeline() |>
    add_column(:peak_area, compute_peak_area) |>
    add_column(:peak_quality, assess_quality) |>
    filter_rows(row -> row.peak_quality > 0.8) |>
    sort_by([:ms_file_idx, :rt])

# In MaxLFQSearch
normalization_pipeline = TransformPipeline() |>
    add_column(:normalized_intensity, normalize) |>
    filter_rows(row -> !ismissing(row.normalized_intensity)) |>
    sort_by([:protein_group, :precursor_idx])
```

This approach makes file transformations explicit, testable, and reusable across the entire SearchDIA pipeline.