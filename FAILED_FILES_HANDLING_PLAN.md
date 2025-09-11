# Comprehensive Plan for Handling Failed Files in Pioneer.jl Pipeline

## Problem Statement

The Pioneer.jl pipeline has 9 sequential search methods, where files can fail at any stage. When files fail:
1. They don't generate expected output files (Arrow files, intermediate results)
2. Downstream methods must handle missing files gracefully
3. Aggregate methods (MaxLFQ, ScoringSearch) must only process successful files
4. QC plots should only include successful files

Currently, the handling is inconsistent, leading to issues like:
- MaxLFQ missing the third successful file
- QC plots including failed files
- UndefRefError when accessing undefined array elements

## Critical Design Constraint: Index-Based File Tracking

**The file index (position in the original MS data array) is the primary key throughout the pipeline.**

When the pipeline starts:
1. MS data files are loaded into an array: `[file1, file2, file3, ..., file12]`
2. Each file gets an index: 1, 2, 3, ..., 12
3. This index is used EVERYWHERE:
   - `ms_file_idx` parameter in all search methods
   - Array indices in `ArrowTableReference` fields
   - Keys in dictionaries for models (RT, mass error, etc.)
   - File naming for intermediate Arrow files

**This index must be preserved even when files fail!**

## Existing Infrastructure

We already have robust infrastructure for tracking failed files:

### SearchContext Failed File Tracking
```julia
# In SearchContext
failed_files::Set{Int64}  # Stores indices of failed files
file_failure_reasons::Dict{Int64, String}  # Maps index to failure reason

# Helper functions
markFileFailed!(ctx::SearchContext, ms_file_idx::Int64, reason::String)
isFileFailed(ctx::SearchContext, ms_file_idx::Int64)::Bool
shouldSkipFile(ctx::SearchContext, ms_file_idx::Int64)::Bool
getValidFileIndices(ctx::SearchContext)::Vector{Int64}  # Returns indices of non-failed files
```

### ArrowTableReference Structure
```julia
struct ArrowTableReference{N}
    file_paths::NTuple{N, String}  # Original MS data paths
    file_id_to_name::NTuple{N, String}
    first_pass_psms::Vector{String}    # Length N, indexed by ms_file_idx
    second_pass_psms::Vector{String}   # Length N, indexed by ms_file_idx
    passing_psms::Vector{String}       # Length N, indexed by ms_file_idx
    passing_proteins::Vector{String}   # Length N, indexed by ms_file_idx
    rt_index_paths::Vector{String}     # Length N, indexed by ms_file_idx
    failed_search_indicator::Vector{Bool}  # Length N, indexed by ms_file_idx
end
```

**Key insight**: These arrays have length N (total files) and use `ms_file_idx` as the index.

## The Current Problem

When `writeArrow` is called in search methods:
```julia
writeArrow(getPassingPsms(getMSData(search_context))[ms_file_idx], passing_psms)
```

This should:
1. Write the Arrow file to disk
2. Store the path at index `ms_file_idx` in the array

But the array is initialized with `Vector{String}(undef, n)`, creating undefined references that cause issues when:
- Failed files never set their paths (remain undefined)
- Downstream methods try to iterate over all paths

## Proposed Solution

### 1. Fix Arrow Path Storage

#### Option A: Initialize with Empty Strings (Recommended)
```julia
# In ArrowTableReference constructor, change:
Vector{String}(undef, n) → fill("", n)
```

This ensures:
- No undefined references
- Empty string indicates "no file generated"
- Safe to check with `!isempty(path) && isfile(path)`

#### Option B: Track Path Setting
Ensure `writeArrow` properly sets the path:
```julia
function writeArrow(path_slot::Ref{String}, data::DataFrame, ms_file_idx::Int64, base_dir::String, file_name::String)
    full_path = joinpath(base_dir, file_name)
    Arrow.write(full_path, data)
    path_slot[] = full_path  # Store the path at the correct index
end
```

### 2. Standardized Path Retrieval Functions

Create utility functions that respect the index-based system:

```julia
"""
Get paths for valid files, maintaining index association.
Returns: Vector of (index, path) tuples
"""
function get_valid_indexed_paths(path_array::Vector{String}, search_context::SearchContext)
    valid_indices = getValidFileIndices(search_context)
    indexed_paths = Tuple{Int64, String}[]
    
    for idx in valid_indices
        if idx <= length(path_array) && !isempty(path_array[idx]) && isfile(path_array[idx])
            push!(indexed_paths, (idx, path_array[idx]))
        end
    end
    
    return indexed_paths
end

"""
Get paths for valid files as a simple vector.
"""
function get_valid_paths_only(path_array::Vector{String}, search_context::SearchContext)
    return [path for (_, path) in get_valid_indexed_paths(path_array, search_context)]
end

"""
Get file names for valid files, preserving index order.
"""
function get_valid_file_names(search_context::SearchContext)
    valid_indices = getValidFileIndices(search_context)
    all_names = getFileIdToName(getMSData(search_context))
    
    # Return names in index order
    return [all_names[idx] for idx in sort(valid_indices)]
end
```

### 3. Update MaxLFQSearch (Immediate Fix)

```julia
# In MaxLFQSearch summarize_results!

# Get all passing PSM paths array (length = total files)
all_passing_psm_paths = getPassingPsms(getMSData(search_context))

# Get only valid file paths with their indices
valid_indexed_paths = get_valid_indexed_paths(all_passing_psm_paths, search_context)

# Extract just the paths for PSMFileReference creation
existing_passing_psm_paths = [path for (_, path) in valid_indexed_paths]

# Get file names in the same order
successful_file_names = get_valid_file_names(search_context)

@user_info "MaxLFQ will process $(length(successful_file_names)) successful files: $(join(successful_file_names, ", "))"

# Create file references
psm_refs = [PSMFileReference(path) for path in existing_passing_psm_paths]

# Use successful_file_names for all downstream processing
```

### 4. File Naming Convention for Arrow Files

To maintain the index connection, Arrow files should be named using either:

#### Option A: Index-based naming
```julia
"psms_$(ms_file_idx).arrow"  # psms_1.arrow, psms_2.arrow, etc.
```

#### Option B: Name-based with index lookup
```julia
"$(file_name).arrow"  # Use getFileIdToName to map back to index
```

Currently using Option B, which works but requires name→index mapping.

### 5. ScoringSearch Updates

ScoringSearch writes many intermediate files. Ensure all use consistent naming:

```julia
# For each intermediate file type:
for ms_file_idx in getValidFileIndices(search_context)
    file_name = getFileIdToName(getMSData(search_context), ms_file_idx)
    output_path = joinpath(output_dir, "$(file_name)_scored.arrow")
    # Process and write file
    # Store path at correct index
    scored_paths[ms_file_idx] = output_path
end
```

### 6. Verification Functions

Add verification to ensure index consistency:

```julia
"""
Verify that all paths in an array correspond to their correct indices.
"""
function verify_path_indices(path_array::Vector{String}, search_context::SearchContext)
    for idx in getValidFileIndices(search_context)
        if !isempty(path_array[idx])
            expected_name = getFileIdToName(getMSData(search_context), idx)
            actual_name = basename(path_array[idx])
            if !contains(actual_name, expected_name)
                @warn "Index mismatch: idx=$idx, expected=$expected_name, actual=$actual_name"
            end
        end
    end
end
```

## Implementation Steps

### Immediate (Fix current issue):
1. Update MaxLFQSearch to use `getValidFileIndices()` properly
2. Ensure it captures all 3 successful files

### Short Term:
1. Change `Vector{String}(undef, n)` to `fill("", n)` in ArrowTableReference
2. Add utility functions to SearchMethods.jl
3. Update all aggregate methods to use these utilities

### Long Term:
1. Add index verification throughout pipeline
2. Standardize file naming conventions
3. Add comprehensive tests for failed file scenarios

## Benefits

1. **Index Preservation**: Original file indices maintained throughout
2. **No Undefined References**: Empty strings instead of undefined
3. **Clear Tracking**: Failed vs successful files clearly identified
4. **Robust Aggregation**: Methods correctly process only successful files
5. **Debugging**: Easy to trace which index corresponds to which file

## Key Insight

The index is the unchanging identifier that links:
- Original MS data file
- All intermediate Arrow files
- All models and parameters
- Final results

By preserving this index relationship and properly tracking which indices represent failed files, we can robustly handle partial pipeline failures.