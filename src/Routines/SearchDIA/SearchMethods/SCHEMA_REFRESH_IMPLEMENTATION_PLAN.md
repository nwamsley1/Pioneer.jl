# Detailed Plan for Approach 3: Use FileOperations for All Schema Changes

## Overview
Ensure all file operations that modify schema go through the FileOperations abstraction layer, which automatically maintains FileReference metadata integrity.

## Core Principle
- Never use `writeArrow` or `Arrow.write` directly (except within FileOperations.jl)
- All file modifications must go through FileOperations.jl
- FileReference objects automatically stay synchronized
- Use `writeArrow` from utils/writeArrow.jl for Windows compatibility

## Implementation Plan

### Phase 1: Add Missing FileOperations Functions

#### 1.1 Add a general-purpose write function
```julia
# In FileOperations.jl
"""
    write_arrow_file(ref::FileReference, df::DataFrame) -> FileReference
    
Write a DataFrame to the file referenced by ref, updating all metadata.
Uses writeArrow from utils/writeArrow.jl to handle Windows file locking issues.
"""
function write_arrow_file(ref::FileReference, df::DataFrame)
    # Use writeArrow which handles Windows-specific issues
    writeArrow(file_path(ref), df)
    
    # Update reference metadata
    new_ref = create_reference(file_path(ref), typeof(ref))
    ref.schema = schema(new_ref)
    ref.row_count = row_count(new_ref)
    ref.sorted_by = ()  # Reset sort state as we don't know if df is sorted
    
    return ref
end
```

#### 1.2 Add a transform-and-write function
```julia
"""
    transform_and_write!(ref::FileReference, transform_fn::Function) -> FileReference
    
Load entire file, apply transformation, and write back.
For operations that need full dataset access (like sorting).
"""
function transform_and_write!(ref::FileReference, transform_fn::Function)
    validate_exists(ref)
    
    # Load entire file
    df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
    
    # Apply transformation
    transformed_df = transform_fn(df)
    
    # Write back using writeArrow for Windows compatibility
    return write_arrow_file(ref, transformed_df)
end
```

#### 1.3 Add column operations with sorting
```julia
"""
    add_column_and_sort!(ref::FileReference, col_name::Symbol, 
                        compute_fn::Function, sort_keys::Symbol...; 
                        reverse::Bool=false) -> FileReference
    
Add a column and immediately sort by specified keys.
"""
function add_column_and_sort!(ref::FileReference, col_name::Symbol,
                             compute_fn::Function, sort_keys::Symbol...;
                             reverse::Bool=false)
    # Add column
    add_column_to_file!(ref, col_name, compute_fn)
    
    # Sort by keys
    sort_file_by_keys!(ref, sort_keys...; reverse=reverse)
    
    return ref
end
```

#### 1.4 Update existing functions to use writeArrow
```julia
# Update add_column_to_file! to use writeArrow instead of Arrow.write
# Update sort_file_by_keys! to use writeArrow instead of Arrow.write
# Update all other functions in FileOperations.jl that write files
```

### Phase 2: Create Reference-Based Wrappers in scoring_interface.jl

#### 2.1 Wrap calculate_global_protein_scores
```julia
"""
    calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    
Calculate global protein scores and add them to files via references.
"""
function calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    # First pass: collect max scores
    acc_to_max_pg_score = Dict{Tuple{String,Bool,UInt8}, Float32}()
    
    for ref in pg_refs
        process_with_memory_limit(ref) do batch
            for row in eachrow(batch)
                key = (row.protein_name, row.target, row.entrap_id)
                old = get(acc_to_max_pg_score, key, -Inf32)
                acc_to_max_pg_score[key] = max(row.pg_score, old)
            end
        end
    end
    
    # Second pass: add global_pg_score column and sort
    for ref in pg_refs
        add_column_and_sort!(ref, :global_pg_score, 
            batch -> begin
                scores = Vector{Float32}(undef, nrow(batch))
                for i in 1:nrow(batch)
                    key = (batch.protein_name[i], batch.target[i], batch.entrap_id[i])
                    scores[i] = get(acc_to_max_pg_score, key, batch.pg_score[i])
                end
                scores
            end,
            :global_pg_score, :target;  # sort keys
            reverse=true
        )
    end
    
    return acc_to_max_pg_score
end
```

#### 2.2 Wrap perform_probit_analysis operations
```julia
"""
    apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference}, 
                        β_fitted::Vector{Float64}, feature_names::Vector{Symbol})
    
Apply probit regression scores to protein group files.
"""
function apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference},
                             β_fitted::Vector{Float64},
                             feature_names::Vector{Symbol})
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Calculate probit scores
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted)
            
            # Update scores
            df[!, :old_pg_score] = copy(df.pg_score)
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by new scores
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            return df
        end
    end
end
```

### Phase 3: Update All Direct File Operations

#### 3.1 Replace writeArrow in calculate_global_protein_scores
```julia
# OLD (lines 1104-1128):
for pg_path in passing_pg_paths
    df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
    # ... add global_pg_score ...
    sort!(df, [:global_pg_score, :target], rev = [true, true])
    writeArrow(pg_path, df)
end

# NEW:
pg_refs = [ProteinGroupFileReference(path) for path in passing_pg_paths]
calculate_and_add_global_scores!(pg_refs)
```

#### 3.2 Replace writeArrow in perform_probit_analysis
```julia
# OLD (lines 1580-1612):
for pg_path in pg_paths
    df = DataFrame(Tables.columntable(Arrow.Table(pg_path)))
    # ... calculate probit scores ...
    df[!, :pg_score] = Float32.(prob_scores)
    sort!(df, [:pg_score, :target], rev = [true, true])
    writeArrow(pg_path, df)
end

# NEW:
pg_refs = [ProteinGroupFileReference(path) for path in pg_paths]
apply_probit_scores!(pg_refs, β_fitted, feature_names)
```

#### 3.3 Update other writeArrow calls
Similar pattern for all other writeArrow calls:
1. Create FileReference for the file
2. Use appropriate FileOperations function
3. Let FileOperations handle metadata updates

### Phase 4: Update ScoringSearch.jl Integration

#### 4.1 Maintain references throughout pipeline
```julia
# After protein inference (Step 11)
pg_refs = [ref.protein_ref for ref in scoring_refs.paired_files]

# Step 12: Probit regression
perform_protein_probit_regression_refs(pg_refs, ...)

# Step 13: Global scores (using refs)
acc_to_max_pg_score = calculate_and_add_global_scores!(pg_refs)

# Step 17: Sort protein tables (already using refs)
for ref in pg_refs
    sort_file_by_keys!(ref, :global_pg_score; reverse=true)
end
```

#### 4.2 Update function signatures
Update functions to accept references instead of paths where appropriate.

### Phase 5: Add Safety Checks

#### 5.1 Add validation in FileOperations
```julia
# Validate schema changes are consistent
function validate_schema_change(old_schema::FileSchema, new_schema::FileSchema)
    # Check that no columns were accidentally removed
    removed = setdiff(old_schema.columns, new_schema.columns)
    if !isempty(removed)
        @warn "Columns removed during operation: $removed"
    end
end
```

## Benefits of This Approach

1. **Automatic Synchronization**: FileReference metadata always matches file state
2. **No Manual Updates**: No need to remember to call refresh_schema!
3. **Centralized Logic**: All file operations in one place
4. **Type Safety**: Can't forget to update references
5. **Future Flexibility**: Easy to change file format or add optimizations
6. **Windows Compatibility**: Using writeArrow ensures proper file handling on Windows

## Migration Strategy

1. **Phase 1**: Add new FileOperations functions (non-breaking)
2. **Phase 2**: Create wrapper functions in scoring_interface.jl
3. **Phase 3**: Update one function at a time to use new wrappers
4. **Phase 4**: Test each change incrementally
5. **Phase 5**: Remove or deprecate old functions

## Testing Plan

1. **Unit Tests**:
   - Test each new FileOperations function
   - Verify schema updates correctly
   - Check sort state preservation
   - Test on Windows systems

2. **Integration Tests**:
   - Run full ScoringSearch pipeline
   - Verify all columns present in output
   - Check sorting is maintained

3. **Performance Tests**:
   - Compare timing before/after changes
   - Monitor memory usage

## Risk Mitigation

- Keep old functions during transition
- Add extensive logging during migration
- Test on small datasets first
- Have rollback plan ready
- Test specifically on Windows systems

## Example Usage After Implementation

```julia
# All file operations go through references
pg_refs = [ProteinGroupFileReference(path) for path in pg_paths]

# Add column - reference automatically updated
add_column_to_file!(ref, :new_score, compute_fn)

# Transform entire file - reference automatically updated  
transform_and_write!(ref) do df
    # Any transformation
    sort!(df, :score)
    return df
end

# Sort - reference tracks sort state
sort_file_by_keys!(ref, :score, :target; reverse=true)

# Schema always current
@assert has_column(schema(ref), :new_score)
```

This approach ensures complete abstraction of file operations while maintaining reference integrity automatically and handling Windows-specific file system issues.