# FileOperations and ScoringSearch Interaction Summary

## Overview

This document summarizes how FileOperations.jl and ScoringSearch interact, documenting the abstraction layers and data flow patterns established during the refactoring.

## Architecture Layers

### 1. FileReference Abstraction Layer
**Purpose**: Type-safe metadata about files without direct file access

```julia
# Abstract base type
abstract type FileReference end

# Concrete implementations
PSMFileReference      # For peptide-spectrum match files
ProteinGroupFileReference  # For protein group files

# Common interface
file_path(ref)    # Get file path
exists(ref)       # Check existence  
row_count(ref)    # Get row count
is_sorted_by(ref, cols...)  # Check sort state
```

### 2. FileOperations Layer
**Purpose**: All file I/O operations, enforcing safety and efficiency

Key functions used by ScoringSearch:
- `stream_transform()` - Transform file with streaming, write to new location
- `transform_and_write!()` - Transform file in-place
- `merge_psm_files()` - Merge sorted PSM files by key
- `merge_protein_groups_by_score()` - Specialized protein group merge
- `apply_protein_inference()` - Wrapper for protein inference algorithm

### 3. ScoringSearch Interface Layer  
**Purpose**: High-level operations specific to scoring workflow

Key functions in `scoring_interface.jl`:
- `filter_psms_by_qvalue()` - Creates new filtered files
- `sort_and_filter_quant_tables_refs()` - In-place filtering/sorting
- `add_global_scores_to_psms_refs()` - Adds scores in-place
- `perform_probit_regression_protein()` - ML protein scoring

## Data Flow in ScoringSearch

### Step 1: PSM Scoring
```julia
# Score PSMs using XGBoost (modifies files in second_pass_folder)
score_precursor_isotope_traces(
    second_pass_folder,
    file_paths,
    precursors,
    ...
)
```

### Step 2: Best Trace Selection
```julia
# Find best isotope traces
best_traces = get_best_traces(second_pass_psms_paths)

# Create references for second pass PSMs
second_pass_refs = [PSMFileReference(path) for path in second_pass_psms_paths]
```

### Step 3: Filter and Sort
```julia
# In-place modification of second pass files
filtered_refs = sort_and_filter_quant_tables_refs(
    second_pass_refs,
    isotope_tracetype,
    :global_prob,
    best_traces
)
# Returns same references after modifying files
```

### Step 4: Merge for Scoring
```julia
# Merge PSM files sorted by probability
merge_psm_files(filtered_refs, merged_path, :global_prob)
```

### Step 5: Q-value Filtering
```julia
# Creates NEW files in passing_psms_folder
passing_refs = filter_psms_by_qvalue(
    filtered_refs,
    passing_psms_folder,
    precursors,
    qval_interp,
    threshold
)
# Returns references to new filtered files
```

### Step 6: Protein Inference
```julia
# Creates protein group files
psm_to_pg_path = perform_protein_inference(
    passing_psm_paths,
    precursors,
    output_folder,
    entrapment_groups,
    target_proteins,
    decoy_proteins
)
```

### Step 7: Protein Scoring
```julia
# Create paired references for internal use
paired_files = PairedSearchFiles[]
for (idx, ref) in enumerate(passing_refs)
    if haskey(psm_to_pg_path, file_path(ref))
        pg_path = psm_to_pg_path[file_path(ref)]
        push!(paired_files, PairedSearchFiles(ref, pg_path, idx))
    end
end

# Use for protein regression/scoring
scoring_refs = ScoringSearchResultRefs(paired_files)
```

## Key Design Patterns

### 1. Reference Lifecycle
- Created when needed, not stored long-term
- Used for type safety and metadata tracking
- Disposed after operations complete

### 2. In-Place vs New Files
- **In-place**: Sorting, adding columns, filtering best traces
- **New files**: Q-value filtering, protein inference
- Choice based on efficiency and workflow requirements

### 3. Streaming Operations
- Never load full datasets into memory
- Process in configurable batches
- Maintain sort order during operations

### 4. Function Parameter Pattern
```julia
# Avoid naming conflicts by passing functions
apply_protein_inference(psm_ref, output_path, precursors, 
                       protein_inference_fn=getProteinGroupsDict)
```

## Interface Benefits

### Type Safety
- Can't accidentally pass protein file where PSM expected
- Compiler catches type mismatches
- Clear distinction between file types

### Efficiency
- Single-pass operations where possible
- Streaming prevents memory issues
- Batched processing for large files

### Maintainability  
- File operations centralized in FileOperations.jl
- Search logic separated from I/O concerns
- Easy to test with mock references

### Flexibility
- Generic operations work on abstract FileReference
- Specialized operations for specific types
- Easy to add new file types

## Common Usage Patterns

### Pattern 1: Transform and Track
```julia
ref = PSMFileReference(path)
transform_and_write!(ref) do df
    # Modify dataframe
    df.new_col = compute_values(df)
    return df
end
mark_sorted!(ref, :new_col)
```

### Pattern 2: Create New from Old
```julia
input_ref = PSMFileReference(input_path)
output_ref = stream_transform(input_ref, output_path) do batch
    # Transform batch
    filter(row -> row.score > threshold, batch)
end
```

### Pattern 3: Merge Multiple Files
```julia
refs = [PSMFileReference(path) for path in paths]
merged_ref = merge_psm_files(refs, output_path, :score)
```

## Future Improvements

### Composable Pipeline Interface
User requested ability to chain operations:
```julia
pipeline = add_column(:best_trace) |>
           filter_rows(row -> row.best_trace) |>
           sort_by(:score, :target)
           
apply_pipeline!(refs, pipeline)
```

This would make complex operations clearer while maintaining single-pass efficiency.

## Summary

The FileOperations/ScoringSearch interaction demonstrates:
1. Clean separation between I/O and business logic
2. Type-safe file handling preventing common errors
3. Efficient streaming operations for large datasets
4. Flexible abstraction supporting future enhancements

The refactoring successfully encapsulates file operations while maintaining performance and adding safety.