# Failed File Tracking Investigation - Pioneer.jl

## Executive Summary

Pioneer has a comprehensive failed file tracking system in `SearchContext`, but **ScoringSearch has a critical bug** where it doesn't properly handle the mismatch between file indices when writing back filtered results. This causes empty groups to appear in supposedly "valid" files, leading to the `argmax` error on empty collections.

## Failed File Tracking System

### Core Components

1. **SearchContext Fields**:
   ```julia
   failed_files::Set{Int64}           # Set of failed file indices
   file_failure_reasons::Dict{Int64, String}  # Reasons for each failure
   ```

2. **Key Functions**:
   - `markFileFailed!(ctx, ms_file_idx, reason)` - Mark a file as failed
   - `isFileFailed(ctx, ms_file_idx)` - Check if file is failed
   - `shouldSkipFile(ctx, ms_file_idx)` - Check if file should be skipped
   - `getValidFileIndices(ctx)` - Get list of valid (non-failed) file indices
   - `getFailedFiles(ctx)` - Get set of all failed files

### How Search Methods Handle Failed Files

#### ✅ FirstPassSearch (Correct)
```julia
# In processChunk!
if shouldSkipFile(search_context, ms_file_idx)
    # Create empty results for failed file
    markFileFailed!(search_context, ms_file_idx, reason)
    return empty_results
end

# In summarize_results!
valid_indices = getValidFileIndices(search_context)
valid_psms_paths = [all_psms_paths[i] for i in valid_indices]
```

#### ✅ SecondPassSearch (Correct)
Similar pattern - checks `shouldSkipFile` and only processes valid files.

#### ✅ ParameterTuningSearch (Correct)
```julia
if shouldSkipFile(search_context, ms_file_idx)
    # Skip processing, return default parameters
    return default_mass_tolerance
end
```

#### ❌ ScoringSearch (BUGGY)
Has multiple issues with failed file handling.

## The Bug in ScoringSearch

### Issue 1: Index Mismatch When Writing Back Results

**Location**: `ScoringSearch.jl:333-336`

```julia
# Step 2: After MBR filtering, write back to files
for (idx, ref) in enumerate(second_pass_refs)  # ← BUG: idx is 1,2,3... not original file indices!
    sub_df = merged_df[merged_df.ms_file_idx .== idx, :]  # ← WRONG: ms_file_idx is original index
    write_arrow_file(ref, sub_df)
end
```

**The Problem**:
1. `second_pass_refs` contains only valid files (e.g., files 1,3,5 if file 2,4 failed)
2. `enumerate(second_pass_refs)` gives indices 1,2,3...
3. But `ms_file_idx` in the DataFrame uses original file indices (1,3,5)
4. This causes empty DataFrames to be written for some files!

**Example**:
- Original files: [1,2,3,4,5]
- Failed files: [2,4]
- Valid files: [1,3,5]
- `second_pass_refs` enumeration: [(1, ref1), (2, ref3), (3, ref5)]
- When `idx=2`, it looks for `ms_file_idx==2`, but file 2 failed! No data found!
- Empty DataFrame written to ref3 (which is actually file 3's reference)

### Issue 2: No Check for Empty Groups After Writing

After writing potentially empty DataFrames back to files, Step 3 proceeds without checking:

```julia
# Step 3: Find Best Isotope Traces
best_traces = get_best_traces(valid_second_pass_psms, ...)  # Uses same paths
# Step 4: Process with add_best_trace_indicator
apply_pipeline!(second_pass_refs, quant_processing_pipeline)  # Can have empty groups!
```

The `add_best_trace_indicator` function then calls `argmax` on empty groups, causing the error.

### Issue 3: Incomplete Filtering After MBR

MBR filtering can create situations where some precursors have no valid traces:
```julia
# MBR creates cross-run references
# Some precursors may pass MBR filtering but have no data in current run
# These create empty groups when grouped by precursor_idx
```

## Why This Wasn't Caught

1. **Silent Failure**: Empty DataFrames are valid Arrow files, so no error occurs during writing
2. **Delayed Error**: The error only manifests when `argmax` is called on empty groups
3. **Testing Gap**: Tests typically use all-valid or all-failed files, not mixed scenarios
4. **Assumption Violation**: Code assumes enumerate index matches ms_file_idx

## Data Flow Diagram

```
Valid Files: [1,3,5] (files 2,4 failed earlier)
                ↓
Step 1: Score with XGBoost
    second_pass_refs = [ref1, ref3, ref5]
                ↓
Step 2: Merge and compute probabilities
    merged_df has ms_file_idx ∈ {1,3,5}
                ↓
    ❌ Write back loop:
    for (idx, ref) in enumerate(second_pass_refs)
        idx=1 → looks for ms_file_idx==1 ✓ (finds data)
        idx=2 → looks for ms_file_idx==2 ✗ (no data! file 2 failed)
        idx=3 → looks for ms_file_idx==3 ✗ (wrong! should look for 5)
                ↓
Step 3: Some refs now have empty DataFrames
                ↓
Step 4: argmax fails on empty groups
```

## Correct Implementation Pattern

### Option 1: Use Original Indices
```julia
valid_indices = getValidFileIndices(search_context)
for (i, ms_file_idx) in enumerate(valid_indices)
    ref = second_pass_refs[i]
    sub_df = merged_df[merged_df.ms_file_idx .== ms_file_idx, :]
    write_arrow_file(ref, sub_df)
end
```

### Option 2: Create Index Mapping
```julia
# Create mapping from enumeration to actual file indices
file_idx_mapping = Dict(i => valid_indices[i] for i in 1:length(valid_indices))
for (enum_idx, ref) in enumerate(second_pass_refs)
    actual_file_idx = file_idx_mapping[enum_idx]
    sub_df = merged_df[merged_df.ms_file_idx .== actual_file_idx, :]
    write_arrow_file(ref, sub_df)
end
```

### Option 3: Check for Empty DataFrames
```julia
for (idx, ref) in enumerate(second_pass_refs)
    sub_df = merged_df[merged_df.ms_file_idx .== idx, :]
    if nrow(sub_df) == 0
        @warn "No data for file index $idx after filtering"
        # Either skip or handle appropriately
        continue
    end
    write_arrow_file(ref, sub_df)
end
```

## Other Vulnerable Locations

Similar patterns that might have the same issue:

1. **IntegrateChromatogramSearch** - Check how it maps indices
2. **MaxLFQSearch** - Uses `getValidFileIndices` correctly
3. **Any place using `enumerate` with filtered file lists**

## Recommendations

### Immediate Fix
1. Fix the index mismatch in ScoringSearch Step 2
2. Add empty group protection to `add_best_trace_indicator`
3. Add validation after writing DataFrames to check for empty files

### Long-term Improvements
1. **Standardize Index Handling**: Create helper functions for proper index mapping
2. **Add Invariant Checks**: Assert that written DataFrames are non-empty for valid files
3. **Enhance Testing**: Add tests with mixed valid/failed file scenarios
4. **Document Pattern**: Create clear documentation on handling filtered file lists

### Example Helper Function
```julia
"""
Write DataFrame subsets back to file references, properly handling index mapping.
"""
function write_subsets_to_refs(merged_df::DataFrame, refs::Vector{PSMFileReference},
                               valid_indices::Vector{Int64})
    for (i, ref) in enumerate(refs)
        actual_file_idx = valid_indices[i]
        sub_df = merged_df[merged_df.ms_file_idx .== actual_file_idx, :]
        if nrow(sub_df) == 0
            @warn "No data for file $actual_file_idx after filtering"
        end
        write_arrow_file(ref, sub_df)
    end
end
```

## Conclusion

Pioneer's failed file tracking system is well-designed and works correctly in most search methods. However, ScoringSearch has a critical bug where it mismatches indices when writing filtered data back to files. This causes empty DataFrames to be written for valid files, which then triggers the `argmax` error on empty collections.

The fix requires:
1. Correcting the index mapping in ScoringSearch Step 2
2. Adding defensive programming for empty groups
3. Improving test coverage for mixed valid/failed file scenarios