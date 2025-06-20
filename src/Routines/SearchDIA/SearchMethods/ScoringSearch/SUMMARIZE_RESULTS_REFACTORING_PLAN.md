# Plan to Modify summarize_results! in ScoringSearch to Use New Abstractions

## Overview
This plan details how to update the `summarize_results!` function in ScoringSearch.jl starting from Step 3 to use the new file reference abstractions instead of direct file operations.

## Current File Operations Analysis

### Step 3: Process Quantification Results
- **Function**: `sort_and_filter_quant_tables`
- **Current**: Takes paths, loads/modifies/writes files directly
- **Files Modified**: Second pass PSM files (filtered and sorted)

### Step 4: Merge PSM Scores (Global)
- **Function**: `merge_sorted_psms_scores`
- **Current**: Takes paths, creates merged output file
- **Files Created**: Merged PSM file at `results.merged_quant_path`

### Step 6: Re-sort and Merge PSM Scores (Experiment-wide)
- **Function**: `sort_quant_tables` then `merge_sorted_psms_scores`
- **Current**: Modifies PSM files in place, then merges
- **Files Modified**: Individual PSM files re-sorted

### Step 8: Filter PSMs by Q-value
- **Function**: `get_psms_passing_qval`
- **Current**: Creates new filtered PSM files
- **Files Created**: Passing PSM files in `passing_psms_folder`

### Step 11: Protein Inference
- **Function**: `perform_protein_inference`
- **Current**: Updates PSM files, creates protein group files
- **Files Modified**: PSM files (adds protein assignments)
- **Files Created**: Protein group files

### Step 12: Protein Probit Regression
- **Function**: `perform_protein_probit_regression`
- **Current**: Updates protein group files with probit scores
- **Files Modified**: Protein group files

### Steps 13-19: Protein Group Processing
- Various functions that sort, merge, and update protein groups
- Already partially using references (Step 18)

## Implementation Plan

### Phase 1: Create Reference-Based Wrappers

#### 1.1 Add to scoring_interface.jl:
```julia
"""
    process_and_filter_psms(psm_refs::Vector{PSMFileReference}, 
                          output_dir::String,
                          isotope_tracetype::IsotopeTraceType,
                          prob_col::Symbol,
                          best_traces::Set) -> Vector{PSMFileReference}
                          
Filter and sort PSM files based on best traces.
"""
function process_and_filter_psms(psm_refs::Vector{PSMFileReference}, 
                               output_dir::String,
                               isotope_tracetype::IsotopeTraceType,
                               prob_col::Symbol,
                               best_traces::Set)
    # Implementation using FileOperations
    filtered_refs = PSMFileReference[]
    
    for (idx, ref) in enumerate(psm_refs)
        output_path = joinpath(output_dir, "filtered_$(idx).arrow")
        
        # Use streaming filter
        filtered_ref = stream_filter(ref, output_path) do batch
            # Filter logic for best traces
            filter_best_traces(batch, best_traces, isotope_tracetype)
        end
        
        # Sort the filtered file
        sort_file_by_keys!(filtered_ref, prob_col; reverse=true)
        
        push!(filtered_refs, filtered_ref)
    end
    
    return filtered_refs
end

"""
    merge_psm_files(psm_refs::Vector{PSMFileReference},
                   output_path::String,
                   sort_col::Symbol) -> PSMFileReference
                   
Merge multiple PSM files sorted by specified column.
"""
function merge_psm_files(psm_refs::Vector{PSMFileReference},
                       output_path::String,
                       sort_col::Symbol)
    # Ensure all files are sorted
    for ref in psm_refs
        if !is_sorted_by(ref, sort_col)
            @warn "Not sorted by key ! $sort_col"
            sort_file_by_keys!(ref, sort_col; reverse=true)
        end
    end
    
    # Use existing merge function
    return merge_psm_scores(psm_refs, output_path, sort_col)
end
```

#### 1.2 Update utils.jl functions to accept references:
```julia
# Modified signature
function sort_and_filter_quant_tables_refs(
    psm_refs::Vector{PSMFileReference},
    merged_quant_path::String,
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Set
)
    # Use reference-based operations
    for ref in psm_refs
        # Filter and sort using FileOperations
        filtered_df = stream_filter(ref, ...) do batch
            # filtering logic
        end
        sort_file_by_keys!(ref, prob_col; reverse=true)
    end
end
```

### Phase 2: Update summarize_results! Implementation

#### 2.1 Create references early:
```julia
# After Step 2, create references for second pass PSMs
second_pass_refs = [PSMFileReference(path) for path in getSecondPassPsms(getMSData(search_context))]
```

#### 2.2 Update each step:

**Step 3:**
```julia
# OLD:
sort_and_filter_quant_tables(
    getSecondPassPsms(getMSData(search_context)),
    results.merged_quant_path,
    params.isotope_tracetype,
    :global_prob,
    best_traces
)

# NEW:
filtered_refs = process_and_filter_psms(
    second_pass_refs,
    temp_folder,
    params.isotope_tracetype,
    :global_prob,
    best_traces
)
```

**Step 4:**
```julia
# OLD:
merge_sorted_psms_scores(
    getSecondPassPsms(getMSData(search_context)),
    results.merged_quant_path,
    :global_prob
)

# NEW:
merged_ref = merge_psm_files(
    filtered_refs,
    results.merged_quant_path,
    :global_prob
)
```

**Step 6:**
```julia
# OLD:
sort_quant_tables(
    getSecondPassPsms(getMSData(search_context)),
    results.merged_quant_path,
    :prec_prob
)

# NEW:
for ref in filtered_refs
    sort_file_by_keys!(ref, :prec_prob; reverse=true)
end
```

**Step 8:**
```julia
# OLD:
get_psms_passing_qval(
    getPrecursors(getSpecLib(search_context)),
    getPassingPsms(getMSData(search_context)),
    passing_psms_folder,
    passing_psms_paths,
    ...
)

# NEW:
passing_refs = filter_psms_by_qvalue(
    filtered_refs,
    passing_psms_folder,
    getPrecursors(getSpecLib(search_context)),
    results.precursor_global_qval_interp[],
    results.precursor_qval_interp[],
    params.q_value_threshold
)
```

### Phase 3: Handle Paired Files Consistently

#### 3.1 Update protein inference to return references:
```julia
# In perform_protein_inference, return PairedSearchFiles directly
paired_files = perform_protein_inference_refs(
    passing_refs,
    passing_proteins_folder,
    getPrecursors(getSpecLib(search_context)),
    protein_to_possible_peptides,
    min_peptides = params.min_peptides
)
```

#### 3.2 Store references immediately:
```julia
# Step 11.5 (new): Store references right after protein inference
scoring_refs = ScoringSearchResultRefs(paired_files)
store_results!(search_context, ScoringSearch, scoring_refs)
```

#### 3.3 Use stored references for remaining steps:
```julia
# Step 12: Use references
pg_refs = get_protein_refs(scoring_refs)
perform_protein_probit_regression_refs(
    pg_refs,
    params.max_psms_in_memory,
    qc_folder
)

# Step 13: Use references
acc_to_max_pg_score = calculate_global_protein_scores_refs(pg_refs)
```

### Phase 4: Eliminate Direct File Operations

#### 4.1 Replace all writeArrow calls:
```julia
# OLD:
writeArrow(sorted_pg_scores_path, DataFrame())

# NEW:
write_arrow_file(sorted_pg_scores_path, DataFrame())
```

#### 4.2 Replace file removal with reference operations:
```julia
# OLD:
if isfile(sorted_pg_scores_path)
    rm(sorted_pg_scores_path)
end

# NEW:
if exists(merged_pg_ref)
    # Clear file through reference
    update_file_content!(merged_pg_ref, DataFrame())
end
```

## Benefits

1. **Type Safety**: Can't accidentally pass PSM path where protein group expected
2. **Metadata Tracking**: Sort state automatically maintained
3. **Consistency**: All file operations go through FileOperations
4. **Error Prevention**: Validation built into reference operations
5. **Future Flexibility**: Easy to change implementation later

## Migration Strategy

1. **Step 1**: Add new reference-based functions alongside existing ones
2. **Step 2**: Update summarize_results! to use references internally
3. **Step 3**: Deprecate path-based functions
4. **Step 4**: Remove old implementations

## Testing Plan

1. **Unit Tests**: Test each new reference-based function
2. **Integration Test**: Run full pipeline with references
3. **Performance Test**: Ensure no regression
4. **Validation Test**: Check file consistency throughout

## Key Principles

1. **Create references early**: Convert paths to references at the start
2. **Pass references**: Functions should accept references, not paths
3. **Return references**: Functions should return references for created files
4. **Use FileOperations**: All file I/O through the abstraction layer
5. **Maintain state**: References track sort state and metadata

## Example Implementation Pattern

```julia
# OLD: Path-based
function process_files(paths::Vector{String}, output_dir::String)
    for path in paths
        df = DataFrame(Arrow.Table(path))
        # process df
        Arrow.write(joinpath(output_dir, basename(path)), df)
    end
end

# NEW: Reference-based
function process_files(refs::Vector{PSMFileReference}, output_dir::String)
    output_refs = PSMFileReference[]
    for ref in refs
        output_path = joinpath(output_dir, basename(file_path(ref)))
        output_ref = stream_transform(ref, output_path) do batch
            # process batch
        end
        push!(output_refs, output_ref)
    end
    return output_refs
end
```

This plan provides a systematic approach to modernizing the file handling in summarize_results! while maintaining backward compatibility and improving robustness.