# Comprehensive Plan for Failed File Handling in SearchDIA Pipeline

## Problem Statement

The SearchDIA pipeline currently fails when individual MS files encounter errors (e.g., dummy/test files with incorrect schemas, corrupt data, etc.). These failures cause crashes in downstream search methods instead of gracefully handling the problematic files and continuing with valid data.

## Key Issues Identified

### 1. Type System Issues
- Functions expect specific concrete types (e.g., `SplineRtConversionModel`) instead of abstract types
- Fallback models (e.g., `IdentityModel`) don't match expected concrete types
- Need to fix function signatures to accept abstract base types

### 2. Missing Failed File Tracking
- No centralized mechanism to mark files as failed
- Failed files continue to be processed by downstream methods
- No way to exclude failed files from cross-run analyses

### 3. Inadequate Error Propagation
- Errors in individual files crash entire pipeline
- Missing graceful degradation for partial failures
- Downstream methods assume all files succeeded

## Comprehensive Solution Plan

### Phase 1: Fix Type System Issues

#### 1.1 Audit and Fix Function Signatures
**Target Functions:**
- `getBestPrecursorsAccrossRuns.jl:56` - `readPSMs!` function
- Any other functions expecting concrete RT conversion model types

**Actions:**
```julia
# Before (problematic):
function readPSMs!(args..., rt_model::SplineRtConversionModel, args...)

# After (fixed):
function readPSMs!(args..., rt_model::RtConversionModel, args...)
```

**Files to Check:**
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl`
- Any file using RT conversion models in function signatures
- Search for patterns: `::SplineRtConversionModel`, `::IdentityModel`

#### 1.2 Verify Type Hierarchy
**Ensure proper inheritance:**
```julia
abstract type RtConversionModel end
struct SplineRtConversionModel <: RtConversionModel
struct IdentityModel <: RtConversionModel
```

### Phase 2: Implement Failed File Tracking System

#### 2.1 Add Failed File State to SearchContext
**Location:** Core SearchContext structure
```julia
# Add to SearchContext
mutable struct SearchContext{...}
    # ... existing fields ...
    failed_files::Set{Int64}  # Track failed file indices
    file_failure_reasons::Dict{Int64, String}  # Track why files failed
end

# Helper functions
function markFileFailed!(ctx::SearchContext, ms_file_idx::Int64, reason::String)
function isFileFailed(ctx::SearchContext, ms_file_idx::Int64)::Bool
function getFailedFiles(ctx::SearchContext)::Set{Int64}
function getFailureReason(ctx::SearchContext, ms_file_idx::Int64)::String
```

#### 2.2 Integration Points for Failure Marking
**FirstPassSearch - Primary Detection Point:**
```julia
# In FirstPassSearch/FirstPassSearch.jl process_file! catch block
catch e
    file_name = getFileIdToName(getMSData(search_context), ms_file_idx)
    reason = "FirstPassSearch failed: $(typeof(e))"
    markFileFailed!(search_context, ms_file_idx, reason)
    @user_warn "Marked file $file_name as failed: $reason"
    # ... create empty results ...
end
```

**Additional Detection Points:**
- Parameter tuning search failures
- Data loading/parsing errors
- Severe quality issues (no usable scans)

### Phase 3: Implement File Skipping Logic

#### 3.1 Core Skipping Function
```julia
function shouldSkipFile(search_context::SearchContext, ms_file_idx::Int64)::Bool
    return isFileFailed(search_context, ms_file_idx)
end

function getValidFileIndices(search_context::SearchContext)::Vector{Int64}
    total_files = length(getMSData(search_context))
    return [i for i in 1:total_files if !isFileFailed(search_context, i)]
end
```

#### 3.2 Integration in Each Search Method
**Template for all search methods:**
```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    # Check if file should be skipped
    if shouldSkipFile(search_context, ms_file_idx)
        file_name = getFileIdToName(getMSData(search_context), ms_file_idx)
        @user_warn "Skipping $(getSearchMethodName()) for failed file: $file_name"
        return results  # Return early with unchanged results
    end
    
    # Normal processing...
    try
        # ... method-specific logic ...
    catch e
        # Method-specific failure handling
        file_name = getFileIdToName(getMSData(search_context), ms_file_idx)
        reason = "$(getSearchMethodName()) failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        @user_warn "File $file_name failed in $(getSearchMethodName()): $reason"
        # Create appropriate empty/fallback results
    end
end
```

### Phase 4: Handle Cross-Run Analyses

#### 4.1 Filter Failed Files from Cross-Run Functions
**Key functions to modify:**
- `getBestPrecursorsAccrossRuns()` - Filter out failed file PSMs
- `map_retention_times!()` - Skip failed files in RT mapping
- Any quantification across runs
- Protein inference across runs

**Example implementation:**
```julia
function getBestPrecursorsAccrossRuns(psms_paths, prec_mzs, rt_irt; max_q_val)
    # Filter out paths corresponding to failed files
    valid_indices = getValidFileIndices(search_context)
    valid_psms_paths = [psms_paths[i] for i in valid_indices]
    valid_rt_irt = Dict(i => rt_irt[i] for i in valid_indices if haskey(rt_irt, i))
    
    if isempty(valid_psms_paths)
        @user_warn "No valid files for cross-run precursor analysis"
        return empty_result()
    end
    
    # Continue with valid files only
    # ... existing logic ...
end
```

#### 4.2 Update RT Mapping Logic
```julia
function map_retention_times!(search_context, results, params)
    valid_files = getValidFileIndices(search_context)
    
    for ms_file_idx in valid_files  # Only process valid files
        # ... existing RT mapping logic ...
    end
    
    # Set identity models for failed files
    for failed_idx in getFailedFiles(search_context)
        setRtIrtMap!(search_context, IdentityModel(), failed_idx)
        setIrtRtMap!(search_context, IdentityModel(), failed_idx)
    end
end
```

### Phase 5: Output and Reporting

#### 5.1 Failed File Summary Report
```julia
function generateFailedFileReport(search_context::SearchContext, output_dir::String)
    failed_files = getFailedFiles(search_context)
    
    if isempty(failed_files)
        @user_info "All files processed successfully"
        return
    end
    
    report_path = joinpath(output_dir, "failed_files_report.txt")
    open(report_path, "w") do io
        println(io, "Failed Files Report")
        println(io, "==================")
        println(io, "Total files processed: $(length(getMSData(search_context)))")
        println(io, "Failed files: $(length(failed_files))")
        println(io, "Success rate: $(round((1 - length(failed_files)/length(getMSData(search_context)))*100, digits=1))%")
        println(io, "")
        
        for file_idx in failed_files
            file_name = getFileIdToName(getMSData(search_context), file_idx)
            reason = getFailureReason(search_context, file_idx)
            println(io, "File $file_idx: $file_name")
            println(io, "  Reason: $reason")
            println(io, "")
        end
    end
    
    @user_warn "$(length(failed_files)) files failed processing. See report: $report_path"
end
```

#### 5.2 Empty Output Files for Failed Files
```julia
function createEmptyOutputsForFailedFiles(search_context::SearchContext, method_name::String)
    """Create empty but valid output files for failed files to maintain output consistency"""
    for failed_idx in getFailedFiles(search_context)
        file_name = getFileIdToName(getMSData(search_context), failed_idx)
        # Create empty Arrow files with correct schema
        createEmptyMethodOutput(method_name, file_name, failed_idx)
    end
end
```

### Phase 6: Testing Strategy

#### 6.1 Unit Tests
- Test failed file tracking functions
- Test file skipping logic
- Test type system fixes

#### 6.2 Integration Tests
- Run pipeline with mix of valid and invalid files
- Verify failed files are properly excluded
- Confirm valid files process correctly

#### 6.3 Test Files
- Keep existing dummy files as negative test cases
- Add edge cases: partially corrupt files, wrong schema versions
- Test graceful degradation scenarios

### Phase 7: Implementation Priority

#### High Priority (Immediate)
1. **Fix type system issue in `getBestPrecursorsAccrossRuns.jl`**
   - Change `SplineRtConversionModel` to `RtConversionModel` in function signatures
   - This fixes the immediate crash

2. **Add basic failed file tracking to SearchContext**
   - Implement `markFileFailed!` and `isFileFailed` functions
   - Update FirstPassSearch to mark failed files

#### Medium Priority
3. **Implement file skipping in downstream methods**
   - Add skipping logic to each search method's `process_file!`
   - Update cross-run analysis functions

4. **Add comprehensive error handling**
   - Wrap remaining vulnerable functions in try-catch
   - Ensure all methods can handle failed files gracefully

#### Lower Priority
5. **Enhanced reporting and diagnostics**
   - Failed file reports
   - Empty output file generation
   - Comprehensive testing

### Phase 8: Long-term Improvements

#### 8.1 Configurable Failure Tolerance
```julia
# Add to parameters
struct PioneerParameters
    # ... existing fields ...
    max_failed_file_percentage::Float64  # e.g., 0.2 = allow up to 20% failures
    critical_file_failure_threshold::Int  # e.g., fail pipeline if >X critical files fail
end
```

#### 8.2 Advanced Recovery Mechanisms
- Attempt alternative processing for borderline files
- Partial processing for files with specific issues
- Automatic file format detection and conversion

#### 8.3 Performance Optimizations
- Early detection of problematic files
- Skip expensive operations on known-bad files
- Parallel validation of file integrity

## Implementation Checklist

### Immediate Actions (Fix Current Crash)
- [ ] Fix `getBestPrecursorsAccrossRuns.jl` function signature
- [ ] Verify `RtConversionModel` abstract type hierarchy
- [ ] Test with current dummy files

### Phase 1 Implementation  
- [ ] Add failed file tracking to SearchContext
- [ ] Implement helper functions for failure tracking
- [ ] Update FirstPassSearch with failure marking
- [ ] Add file skipping logic to cross-run functions

### Phase 2 Implementation
- [ ] Add skipping logic to all search methods
- [ ] Comprehensive error handling in vulnerable functions
- [ ] Testing with mixed valid/invalid file sets

### Phase 3 Implementation
- [ ] Failed file reporting
- [ ] Empty output file generation
- [ ] Documentation and user guidance

## Success Metrics

1. **Robustness:** Pipeline completes successfully with mix of valid and invalid files
2. **Informativeness:** Clear reporting of which files failed and why
3. **Consistency:** Output structure maintained even with failed files
4. **Performance:** Minimal overhead from failure tracking
5. **Usability:** Users can easily identify and address file issues

This plan provides a systematic approach to making the SearchDIA pipeline robust against individual file failures while maintaining full functionality for valid data.