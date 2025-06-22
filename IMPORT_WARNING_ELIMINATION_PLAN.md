# Import Warning Elimination Plan

## Problem Statement

The current `importScripts()` system in SearchDIA has several issues:
1. **Method Overwriting Warnings**: When the same function is defined in multiple files, Julia shows warnings
2. **Inconsistent Import Logic**: Mix of direct `include()` and undefined `safe_include!` 
3. **Code Errors**: Typos and scope issues in the current implementation
4. **No Conflict Detection**: System doesn't know which files would cause conflicts

## Root Cause Analysis

The current `safe_include!` pattern only prevents including the same **file** twice:
```julia
function safe_include!(files_loaded::Set{String}, file_path::String)
    if file_path ∉ files_loaded  # Only checks FILE paths
        include(file_path)        # Different files can still define same functions
        push!(files_loaded, file_path)
    end
end
```

This doesn't prevent **method overwriting warnings** when different files define the same function.

## Enhanced Solution: Method-Aware Safe Import

### Core Concept
Create a system that tracks which **methods/functions** have been defined, not just which files have been loaded.

### Implementation Strategy

#### 1. Function Name Extraction
Parse Julia files to extract function definitions without executing them:
```julia
function extract_function_names(file_path::String)::Set{Symbol}
    # Parse file content to find function definitions
    # Return set of function names that would be defined
end
```

#### 2. Enhanced Import Tracker
```julia
struct SafeImportTracker
    files_loaded::Set{String}
    methods_defined::Set{Symbol}  # Track function names
    method_sources::Dict{Symbol, String}  # Track which file defined each method
    conflicts_skipped::Vector{String}  # Track skipped files
end
```

#### 3. Conflict-Aware Import Function
```julia
function safe_include_with_method_check!(tracker::SafeImportTracker, file_path::String)
    # 1. Check if file already loaded
    if file_path in tracker.files_loaded
        return false
    end
    
    # 2. Extract potential function definitions
    potential_methods = extract_function_names(file_path)
    
    # 3. Check for conflicts with already-defined methods
    conflicts = intersect(potential_methods, tracker.methods_defined)
    if !isempty(conflicts)
        @warn "Skipping $file_path - would redefine: $conflicts"
        push!(tracker.conflicts_skipped, file_path)
        return false
    end
    
    # 4. Safe to include - no conflicts detected
    include(file_path)
    push!(tracker.files_loaded, file_path)
    union!(tracker.methods_defined, potential_methods)
    for method in potential_methods
        tracker.method_sources[method] = file_path
    end
    return true
end
```

### Benefits of This Approach

1. **Eliminates Method Overwriting Warnings**: System prevents conflicting method definitions
2. **Clear Conflict Reporting**: Shows exactly which files are skipped and why
3. **Better Debugging**: Tracks which file defined each method
4. **Automatic Conflict Resolution**: No manual file management needed
5. **Maintains Functionality**: SearchDIA continues to work with available methods

### Integration Plan

1. **Phase 1**: Implement function name extraction
2. **Phase 2**: Create enhanced import tracker
3. **Phase 3**: Replace current importScripts with enhanced version
4. **Phase 4**: Test with existing files and identify natural conflicts
5. **Phase 5**: Add file_utils.jl back and let the system automatically handle conflicts

### Alternative: Simple Conflict Detection

If function parsing proves complex, a simpler approach:
1. Try to include each file in a test module
2. If it causes method overwriting warnings, skip it
3. Report which files were skipped and why

This would still solve the core problem without requiring sophisticated parsing.