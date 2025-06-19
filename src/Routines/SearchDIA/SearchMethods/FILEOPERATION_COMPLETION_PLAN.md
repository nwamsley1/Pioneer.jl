# Plan to Complete FileOperations Integration

## Overview
Complete the migration of remaining `writeArrow` calls to use FileOperations abstractions for better encapsulation and automatic metadata management.

## Analysis of Remaining Work

### 1. OOM Version of perform_probit_analysis ✅ COMPLETED
Currently uses paths and direct file operations. Need to:
- ✅ Update function signature to accept `ProteinGroupFileReference` objects
- ✅ Use `row_count(ref)` instead of loading tables to count rows
- ✅ Fixed table loading to use `file_path(ref)`

### 2. Remaining writeArrow Calls in utils.jl
Found 11 writeArrow calls that need evaluation:

#### Category A: Should Use FileOperations (5 calls)
1. **sort_and_filter_quant_tables** (line 138) - ✅ Created sort_and_filter_quant_tables_refs
2. **add_global_scores_to_psms** (line 1126) - ✅ Created add_global_scores_to_psms_refs
3. **update_psms_with_probit_scores** (line 1269) - ✅ Created update_psms_with_probit_scores_refs
4. **write_protein_groups** (line 1845, 2027) - Creates new protein group files
5. **get_proteins_passing_qval** (line 2066) - ✅ Created get_proteins_passing_qval_refs

#### Category B: Keep as writeArrow (6 calls)
1. **merge_sorted_psms_scores** (lines 265, 339) - Complex streaming merge
2. **write_protein_groups_arrow** (line 586) - Low-level writer
3. **perform_probit_analysis_oom** (line 1484) - Already handled by apply_probit_scores!
4. **update_best_psms!** (line 164) - Temporary working files

## Implementation Status

### Phase 1: Update perform_probit_analysis_oom ✅ COMPLETED

#### 1.1 Update Function Signature
```julia
function perform_probit_analysis_oom(
    pg_refs::Vector{ProteinGroupFileReference},  # Changed from paths
    total_protein_groups::Int, 
    max_protein_groups_in_memory::Int, 
    qc_folder::String
)
```

#### 1.2 Use References for Sampling
```julia
# Use row_count instead of loading table
for ref in pg_refs
    if exists(ref)
        n_rows = row_count(ref)
        n_sample = ceil(Int, n_rows * sampling_ratio)
        # ... sampling logic ...
    end
end
```

#### 1.3 Apply Model Using References
```julia
# Replace direct writeArrow with:
apply_probit_scores!(pg_refs, β_fitted, feature_names)

# Count statistics after applying
for ref in pg_refs
    df = DataFrame(Arrow.Table(file_path(ref)))
    total_targets += sum(df.target)
    total_decoys += sum(.!df.target)
end
```

### Phase 2: Convert Key writeArrow Calls ✅ MOSTLY COMPLETED

#### 2.1 sort_and_filter_quant_tables ✅ COMPLETED
Created reference-based version:
```julia
function sort_and_filter_quant_tables_refs(
    psm_refs::Vector{PSMFileReference},
    isotope_trace_type::IsotopeTraceType,
    prob_col::Symbol,
    best_traces::Union{Nothing, Set}
)
    for ref in psm_refs
        if best_traces !== nothing
            # Use stream_filter and sort_file_by_keys!
            filtered_ref = stream_filter(ref, temp_path) do batch
                # filtering logic
            end
            mv(temp_path, file_path(ref), force=true)
        end
        sort_file_by_keys!(ref, prob_col, :target; reverse=true)
    end
end
```

#### 2.2 add_global_scores_to_psms ✅ COMPLETED
Created add_global_scores_to_psms_refs:
```julia
# Create PSM references
psm_refs = [PSMFileReference(path) for path in passing_psms_paths]

# Add global scores using references
for ref in psm_refs
    add_column_to_file!(ref, :global_pg_score) do batch
        # calculate global scores
    end
end
```

#### 2.3 update_psms_with_probit_scores ✅ COMPLETED
Created update_psms_with_probit_scores_refs:
```julia
function update_psms_with_probit_scores_refs(
    psm_refs::Vector{PSMFileReference},
    pg_refs::Vector{ProteinGroupFileReference},
    acc_to_max_pg_score::Dict,
    pg_score_to_qval,
    global_pg_score_to_qval
)
    # Use update_psms_with_scores or create new wrapper
    for (psm_ref, pg_ref) in zip(psm_refs, pg_refs)
        transform_and_write!(psm_ref) do df
            # Add all score columns
            # Return transformed df
        end
    end
end
```

### Phase 3: Handle Protein Group Writing

#### 3.1 Create write_protein_groups_ref
```julia
function write_protein_groups_ref(
    output_path::String,
    protein_groups_dict::Dict,
    features_dict::Dict
) -> ProteinGroupFileReference
    # Build DataFrame
    df = create_protein_groups_dataframe(protein_groups_dict, features_dict)
    
    # Write using writeArrow
    writeArrow(output_path, df)
    
    # Return reference
    return ProteinGroupFileReference(output_path)
end
```

#### 3.2 Update get_proteins_passing_qval ✅ COMPLETED
Created get_proteins_passing_qval_refs:
```julia
function get_proteins_passing_qval_refs(
    pg_refs::Vector{ProteinGroupFileReference},
    output_dir::String,
    global_pg_score_to_qval,
    pg_score_to_qval,
    q_value_threshold::Float32
)
    for ref in pg_refs
        output_path = joinpath(output_dir, basename(file_path(ref)))
        
        # Use stream_filter
        filtered_ref = stream_filter(ref, output_path) do row
            # Check q-value thresholds
            row.global_pg_qval <= q_value_threshold || 
            row.pg_qval <= q_value_threshold
        end
    end
end
```

## Benefits

1. **Consistency**: All file operations use same abstraction
2. **Metadata Management**: Automatic schema/row count updates
3. **Type Safety**: Can't mix PSM and protein references
4. **Memory Efficiency**: Better use of streaming operations
5. **Debugging**: Easier to track file operations

## Testing Strategy

1. **Unit Tests**: Test each converted function
2. **Integration Test**: Run full pipeline with references
3. **Performance Test**: Ensure no regression
4. **Memory Test**: Verify OOM handling still works

## Implementation Order

1. Update perform_probit_analysis_oom (simplest change)
2. Convert sort_and_filter_quant_tables
3. Update PSM scoring functions
4. Handle protein group writing
5. Update get_proteins_passing_qval

This approach maintains backward compatibility while progressively migrating to the reference-based system.