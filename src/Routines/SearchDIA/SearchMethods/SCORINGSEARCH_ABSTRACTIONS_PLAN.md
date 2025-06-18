# ScoringSearch Abstractions Integration Plan

## Overview
This document outlines the systematic plan to use the new file reference abstractions throughout the ScoringSearch pipeline.

## Current Status

### ✅ Completed
1. **FileOperations.jl enhancements**:
   - `stream_sorted_merge` now auto-sorts files with warnings
   - `sort_file_by_keys!` uses `Tables.columntable` for safe in-memory operations
   - Generic N-key merge support implemented

2. **ScoringSearch.jl updates**:
   - File references created and stored in SearchContext after protein inference
   - Reference-based merging used for protein groups

3. **Sort consistency fixes**:
   - `write_protein_groups_arrow` now sorts by single key (`:pg_score`)
   - `sort_protein_tables` uses reference-based sorting

## Phase 1: Core Abstractions (DONE)
- ✅ Auto-sorting in merge operations
- ✅ Safe file modification with Tables.columntable
- ✅ Reference-based operations in key functions

## Phase 2: Complete Reference Integration

### 2.1 Add FileOperations Wrappers
Create these functions in FileOperations.jl:

```julia
# Write operations through references
function write_arrow_file(path::String, df::DataFrame) -> FileReference
    Arrow.write(path, df)
    # Determine reference type based on content
    if haskey(names(df), :precursor_idx)
        return PSMFileReference(path)
    elseif haskey(names(df), :protein_name)
        return ProteinGroupFileReference(path)
    else
        error("Unknown file type")
    end
end

# Update file through reference
function update_file_content!(ref::FileReference, df::DataFrame)
    Arrow.write(file_path(ref), df)
    # Update metadata if needed
    ref.row_count = nrow(df)
    ref.sorted_by = ()  # Reset sort state
end
```

### 2.2 Replace Direct File Operations

#### In utils.jl:
- `writeArrow` calls → `write_arrow_file` or `update_file_content!`
- Direct `Arrow.Table` reads (when modifying) → through FileOperations
- Keep read-only `DataFrame(Arrow.Table(...))` as-is

#### In protein_inference_helpers.jl:
- `writeArrow` in line 366 → `update_file_content!`
- `Arrow.write` in line 445 → `write_arrow_file`

#### In score_psms.jl:
- Any `writeArrow` calls → appropriate FileOperations wrapper

### 2.3 Streaming Operations for Large Files

Add to FileOperations.jl:
```julia
# Process large files in chunks
function process_file_in_chunks(ref::FileReference, 
                               chunk_fn::Function;
                               chunk_size=100_000)
    tbl = Arrow.Table(file_path(ref))
    partitions = Tables.partitioner(tbl, chunk_size)
    
    results = []
    for partition in partitions
        chunk_result = chunk_fn(DataFrame(Tables.columntable(partition)))
        push!(results, chunk_result)
    end
    
    return results
end
```

## Phase 3: Enhanced Validation

### 3.1 Add Consistency Checks
```julia
# Validate PSM-protein group pairing
function validate_paired_files(paired::PairedSearchFiles)
    psm_proteins = get_unique_proteins(paired.psm_ref)
    pg_proteins = get_protein_names(paired.protein_ref)
    
    missing = setdiff(psm_proteins, pg_proteins)
    if !isempty(missing)
        @warn "PSMs reference missing proteins" n_missing=length(missing)
    end
end
```

### 3.2 Add to ScoringSearch workflow
- After protein inference: validate paired files
- Before merging: ensure files are properly paired
- After updates: verify consistency

## Phase 4: Performance Optimizations

### 4.1 Lazy Loading
- Keep memory-mapped reads for analysis
- Only load into memory when modifying

### 4.2 Batch Operations
- Group file operations to minimize I/O
- Use streaming for large merges

### 4.3 Parallel Safety
- Ensure FileReference updates are thread-safe
- Use file locking for concurrent access

## Implementation Priority

1. **High Priority** (prevents errors):
   - ✅ Auto-sorting in merges
   - ✅ Safe file modifications
   - Consistency validation

2. **Medium Priority** (improves design):
   - FileOperations wrappers for all writes
   - Replace direct Arrow.write calls
   - Add validation functions

3. **Low Priority** (nice to have):
   - Performance optimizations
   - Enhanced diagnostics
   - Parallel safety improvements

## Benefits

1. **Robustness**: Auto-sorting prevents merge failures
2. **Safety**: Centralized file operations prevent corruption
3. **Maintainability**: Single point of change for I/O
4. **Performance**: Memory-mapped reads, in-memory writes only when needed
5. **Debugging**: Better error messages and validation

## Next Steps

1. Add FileOperations wrappers (Phase 2.1)
2. Gradually replace direct file operations (Phase 2.2)
3. Add validation to catch issues early (Phase 3)
4. Monitor performance and optimize as needed (Phase 4)

## Notes

- Keep read-only operations memory-mapped for performance
- All modifications should go through FileOperations
- Maintain backward compatibility during transition
- Add tests for new FileOperations functions