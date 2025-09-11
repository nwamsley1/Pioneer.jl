# File Handling Architecture in Pioneer.jl

## Overview

Pioneer.jl uses a robust file handling system designed to handle pipeline failures gracefully. The system distinguishes between **persistent file identifiers** (`ms_file_idx`) and **positional indices** used for iteration, allowing files to fail at any stage without breaking downstream analysis.

## Core Concepts

### ms_file_idx: Persistent File Identifier

Every MS data file in the pipeline is assigned a unique `ms_file_idx` based on its position in the original file list. This identifier:
- **Remains constant** throughout the entire pipeline
- **Is stored in PSM data** as a column `:ms_file_idx`
- **Maps to specific file paths** via `ArrowTableReference.file_paths[ms_file_idx]`
- **Maps to file names** via `ArrowTableReference.file_id_to_name[ms_file_idx]`

**Example:**
```
Original files: [file1.arrow, file2.arrow, file3.arrow, file4.arrow]
ms_file_idx:    [1,           2,           3,           4]

If file3.arrow fails:
- file1.arrow still has ms_file_idx = 1
- file2.arrow still has ms_file_idx = 2  
- file4.arrow still has ms_file_idx = 4
- PSM data contains ms_file_idx values [1, 2, 4] only
```

### Failed File Tracking

Failed files are tracked in:
- `SearchContext.failed_files::Set{Int64}` - Set of failed ms_file_idx values
- `SearchContext.file_failure_reasons::Dict{Int64, String}` - Failure reasons

### ArrowTableReference Structure

```julia
struct ArrowTableReference{N} <: MassSpecDataReference
    file_paths::NTuple{N, String}           # All original file paths
    file_id_to_name::NTuple{N, String}      # Parsed file names
    first_pass_psms::Vector{String}         # PSM output paths
    second_pass_psms::Vector{String}
    passing_psms::Vector{String}
    passing_proteins::Vector{String}
    rt_index_paths::Vector{String}
    failed_search_indicator::Vector{Bool}   # Per-file failure status
end
```

**Key Point:** All arrays in `ArrowTableReference` have **exactly N positions** corresponding to the original file count, ensuring `ms_file_idx` can always be used as a direct array index.

## File Iteration Patterns

### ✅ CORRECT: Use ms_file_idx consistently

```julia
# Main search iteration (SearchMethods.jl)
original_file_indices = collect(1:length(msdr))
for (enum_idx, spectra) in ProgressBar(enumerate(msdr))
    ms_file_idx = original_file_indices[enum_idx]  # Consistent identifier
    
    # Process file with ms_file_idx
    results = process_file!(search_method, ms_file_idx, spectra)
    
    # Store results using ms_file_idx
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, result_path)
end
```

### ✅ CORRECT: Filter by ms_file_idx in PSM data

```julia
# ScoringSearch file write-back
for ref in second_pass_refs
    ms_file_idx = extract_ms_file_idx(file_path(ref))
    if ms_file_idx !== nothing
        # Filter merged data by ms_file_idx 
        sub_df = merged_df[merged_df.ms_file_idx .== ms_file_idx, :]
        write_arrow_file(ref, sub_df)
    end
end
```

### ✅ CORRECT: Dict-based lookups for cross-file data

```julia
# RT conversion models (buildRTIndex.jl)
if haskey(psms, :ms_file_idx) && !isempty(psms[:ms_file_idx])
    file_idx = Int64(first(psms[:ms_file_idx]))
    if haskey(rt_to_irt_splines, file_idx)
        rt_to_irt = rt_to_irt_splines[file_idx]
    elseif haskey(rt_to_irt_splines, enumerate_key)
        rt_to_irt = rt_to_irt_splines[enumerate_key]  # Fallback
    end
end
```

### ✅ CORRECT: Array iteration with zip

```julia
# QC plotting with matched arrays
for (ms_table_path, short_fname) in zip(ms_table_paths, short_fnames)
    # Process matched pairs
    plot_data(ms_table_path, short_fname)
end
```

### ❌ INCORRECT: enumerate with array indexing

```julia
# BROKEN: assumes enumerate index matches file positions
for (enum_idx, file_path) in enumerate(file_paths)
    file_name = file_names[enum_idx]  # WRONG: may be misaligned
end
```

### ❌ INCORRECT: Using enumerate index as ms_file_idx

```julia
# BROKEN: enumerate index != ms_file_idx when files fail
for (file_idx, spectra) in enumerate(msdr)
    # file_idx is positional (1,2,3,4), not persistent identifier
    setResults(context, file_idx, results)  # WRONG
end
```

## Key Components and Functions

### SearchContext Accessors (All use ms_file_idx)

```julia
# Dict-based lookups (safe with failed files)
getIrtErrors(search_context)[ms_file_idx]          # Dict{Int64, Float32}
getNceModel(search_context)[ms_file_idx]           # Dict{Int64, NceModel}  
getMassErrorModel(search_context)[ms_file_idx]     # Dict{Int64, MassErrorModel}
getQuadTransmissionModel(search_context)[ms_file_idx] # Dict{Int64, QuadTransmissionModel}
```

### ArrowTableReference Accessors (All require full arrays)

```julia
# Direct array indexing (requires position ms_file_idx to exist)
getMSData(msdr, ms_file_idx)                    # msdr.file_paths[ms_file_idx]
getFileIdToName(ref, ms_file_idx)              # ref.file_id_to_name[ms_file_idx] 
getFirstPassPsms(ref, ms_file_idx)             # ref.first_pass_psms[ms_file_idx]
getSecondPassPsms(ref, ms_file_idx)            # ref.second_pass_psms[ms_file_idx]
getPassingPsms(ref, ms_file_idx)               # ref.passing_psms[ms_file_idx]

# Setters
setFirstPassPsms!(ref, ms_file_idx, path)      # ref.first_pass_psms[ms_file_idx] = path
setSecondPassPsms!(ref, ms_file_idx, path)     # ref.second_pass_psms[ms_file_idx] = path  
setPassingPsms!(ref, ms_file_idx, path)        # ref.passing_psms[ms_file_idx] = path
```

### Utility Functions

```julia
# Extract ms_file_idx from PSM files
function extract_ms_file_idx(file_path::String) -> Int64?
    table = Arrow.Table(file_path)
    return haskey(table, :ms_file_idx) ? Int64(first(table[:ms_file_idx])) : nothing
end

# Get successful files only
function get_valid_file_names_by_indices(search_context) -> Vector{String}
    all_names = collect(getFileIdToName(getMSData(search_context)))
    valid_indices = getValidFileIndices(search_context)
    return [all_names[i] for i in valid_indices]
end
```

## Locations Using File Indices

### Fixed Locations (Using Dict-based lookups)
- ✅ `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl:80-83` - RT spline lookups
- ✅ All SearchContext model accessors (mass error, RT conversion, etc.)
- ✅ `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl:151-158` - File write-back

### Fixed Locations (Using enumerate correctly)  
- ✅ `src/Routines/SearchDIA/SearchMethods/SearchMethods.jl:54-56` - Main file iteration
- ✅ `src/Routines/SearchDIA/WriteOutputs/qcPlots.jl:375` - TIC plotting with zip

### Critical Locations (Direct array access - requires full arrays)
- ⚠️ `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:346` - `getMSData(msdr, ms_file_idx)`
- ⚠️ `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:348-382` - All ArrowTableReference accessors

### MaxLFQ Integration
- ✅ `src/Routines/SearchDIA/SearchMethods/MaxLFQSearch/MaxLFQSearch.jl:224` - Uses all file names with `collect(getFileIdToName())`

## Testing Guidelines

### Integration Testing
Test with mixed successful/failed files in different orders:
```julia
test_files = [
    "successful1.arrow",     # ms_file_idx = 1
    "failed_dummy.arrow",    # ms_file_idx = 2 (will fail)  
    "successful2.arrow",     # ms_file_idx = 3
    "successful3.arrow",     # ms_file_idx = 4
    "failed_dummy2.arrow"    # ms_file_idx = 5 (will fail)
]
```

### Unit Testing
- Verify PSM data contains correct `ms_file_idx` values
- Test array accessor functions with all positions
- Verify failed file tracking updates correctly

### Common Test Scenarios
1. **All files succeed** - Should work identically to before
2. **Early files fail** - Later files should retain correct ms_file_idx
3. **Late files fail** - Early files should be unaffected  
4. **Mixed success/failure** - Pipeline should complete with correct indices
5. **Duplicate files** - Should handle gracefully without index conflicts

## Migration Checklist

When adding new file-related code:

1. ✅ **Use ms_file_idx consistently** - Don't use enumerate indices as file identifiers
2. ✅ **Use Dict lookups for cross-file data** - Avoid array indexing for dynamic data
3. ✅ **Filter PSM data by ms_file_idx** - Not by positional indices
4. ✅ **Test with failed files** - Ensure robustness to pipeline failures
5. ✅ **Use zip for parallel arrays** - Avoid enumerate with separate indexing

## Anti-Patterns to Avoid

1. **enumerate index as file identifier**
   ```julia
   for (file_idx, data) in enumerate(files)
       store_result(file_idx, result)  # WRONG: file_idx is positional
   ```

2. **Array indexing with dynamic size assumptions**
   ```julia
   file_names = get_successful_files()  # Length = N_successful
   for ms_file_idx in [1, 2, 5, 7]    # Original indices
       name = file_names[ms_file_idx]   # WRONG: may exceed bounds
   ```

3. **Mixing enumerate and ms_file_idx**
   ```julia
   for (i, file) in enumerate(files)
       ms_file_idx = get_file_idx(file)
       results[i] = process(file)        # WRONG: i != ms_file_idx
   ```

## Summary

The Pioneer.jl file handling system provides robust failure handling through:

1. **Persistent identifiers** (`ms_file_idx`) that survive file failures
2. **Separation of concerns** between file identification and iteration
3. **Full array structures** that maintain all positions for direct indexing
4. **Dict-based lookups** for dynamic cross-file data
5. **Explicit failed file tracking** for transparent error handling

This architecture ensures that successful files can complete analysis regardless of which other files fail or in what order they are processed.