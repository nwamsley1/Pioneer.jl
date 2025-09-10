# Search Method Error Handling Template

This template shows how to add comprehensive error handling to all search methods to support failed file tracking.

## Template for process_file! Functions

Every search method should implement this pattern in their `process_file!` function:

```julia
function process_file!(
    results::YourSearchResults,
    params::YourSearchParameters, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # Check if file should be skipped due to previous failure
    if shouldSkipFile(search_context, ms_file_idx)
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        @user_warn "Skipping $(getSearchMethodName()) for previously failed file: $file_name"
        return results  # Return early with unchanged results
    end
    
    # Normal processing wrapped in try-catch
    try
        # Your normal search method logic here...
        # Example:
        # - Data validation
        # - Search algorithm execution  
        # - Result processing
        
    catch e
        # Handle failures gracefully
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        
        reason = "$(getSearchMethodName()) failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        @user_warn "$(getSearchMethodName()) failed for MS data file: $file_name. Error type: $(typeof(e)). File marked as failed."
        
        # Create appropriate empty/fallback results for this method
        # This depends on what your method normally produces
        createFallbackResults!(results, ms_file_idx)
    end
    
    return results
end
```

## Method-Specific Fallback Creation

Each search method needs to implement `createFallbackResults!` appropriate to their output:

### Parameter Tuning Methods
```julia
function createFallbackResults!(results::ParameterTuningSearchResults, ms_file_idx::Int64)
    # Set conservative default parameters
    results.fragment_tolerance[ms_file_idx] = 20.0f0  # Conservative 20 ppm
    results.mass_error_model[ms_file_idx] = MassErrorModel(0.0f0, (20.0f0, 20.0f0))
end
```

### NCE Tuning Methods  
```julia
function createFallbackResults!(results::NceTuningSearchResults, ms_file_idx::Int64)
    # Set identity NCE model (no calibration)
    results.nce_models[ms_file_idx] = IdentityNceModel()
end
```

### Quad Tuning Methods
```julia
function createFallbackResults!(results::QuadTuningSearchResults, ms_file_idx::Int64)
    # Set default quad transmission model
    setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
end
```

### PSM Search Methods (FirstPass, SecondPass, etc.)
```julia
function createFallbackResults!(results::FirstPassSearchResults, ms_file_idx::Int64)
    # Create empty PSM DataFrame with proper schema
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        irt_predicted = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[],
        scan_count = UInt32[],
        fwhm = Float32[]
    )
    results.psms[] = empty_psms
end
```

## Cross-Run Analysis Functions

Functions that analyze data across multiple files need to filter out failed files:

```julia
function cross_run_analysis(search_context::SearchContext, ...)
    # Get only valid files
    valid_indices = getValidFileIndices(search_context)
    
    if isempty(valid_indices)
        @user_warn "No valid files for cross-run analysis"
        return empty_result()
    end
    
    # Filter input data to valid files only
    valid_data = filter_to_valid_files(input_data, valid_indices)
    
    # Proceed with analysis on valid files
    result = perform_analysis(valid_data, ...)
    
    return result
end
```

## Loop-Based Processing

When iterating over files, always check for failed status:

```julia
# Bad - processes all files
for ms_file_idx in 1:length(getMSData(search_context))
    process_file_data(ms_file_idx)
end

# Good - skips failed files  
for ms_file_idx in 1:length(getMSData(search_context))
    if shouldSkipFile(search_context, ms_file_idx)
        continue
    end
    process_file_data(ms_file_idx)
end

# Better - only iterate valid files
valid_files = getValidFileIndices(search_context)
for ms_file_idx in valid_files
    process_file_data(ms_file_idx)
end
```

## Error Message Best Practices

1. **Include file name** for user debugging
2. **Show error type** but not full error object (prevents infinite display loops)
3. **Use consistent messaging** across methods
4. **Log at appropriate level** (@user_warn for expected failures, @user_error for unexpected)

```julia
# Good error handling
file_name = try
    getFileIdToName(getMSData(search_context), ms_file_idx)
catch
    "file_$ms_file_idx"
end

reason = "$(getSearchMethodName()) failed: $(typeof(e))"
markFileFailed!(search_context, ms_file_idx, reason)
@user_warn "$(getSearchMethodName()) failed for MS data file: $file_name. Error type: $(typeof(e)). File marked as failed."

# Bad - can cause infinite loops
@user_error "Search failed: $e"  # Don't show full error object
```

## Files That Need Updates

Apply this template to the following search method files:

### Core Search Methods
- [ ] `ParameterTuningSearch/ParameterTuningSearch.jl`
- [ ] `NceTuningSearch/NceTuningSearch.jl` 
- [ ] `QuadTuningSearch/QuadTuningSearch.jl`
- [x] `FirstPassSearch/FirstPassSearch.jl` (already updated)
- [ ] `HuberTuningSearch/HuberTuningSearch.jl`
- [ ] `SecondPassSearch/SecondPassSearch.jl`
- [ ] `ScoringSearch/ScoringSearch.jl`
- [ ] `IntegrateChromatogramSearch/IntegrateChromatogramSearch.jl`
- [ ] `MaxLFQSearch/MaxLFQSearch.jl`

### Cross-Run Analysis Functions
- [x] `FirstPassSearch/getBestPrecursorsAccrossRuns.jl` (already updated)
- [x] `FirstPassSearch/utils.jl` - map_retention_times! (already updated)
- [ ] Any other functions that process data across multiple files

### Key Points for Each Method

1. **ParameterTuningSearch**: Set conservative mass tolerance fallbacks
2. **NceTuningSearch**: Use identity NCE models for failed files  
3. **QuadTuningSearch**: Use default quad transmission models
4. **HuberTuningSearch**: Use default Huber delta values
5. **SecondPassSearch**: Create empty PSM results with proper schema
6. **ScoringSearch**: Skip failed files in ML training, create empty protein groups
7. **IntegrateChromatogramSearch**: Skip chromatogram extraction for failed files
8. **MaxLFQSearch**: Exclude failed files from quantification matrix

## Testing Strategy

1. **Unit tests**: Test each method's error handling with problematic inputs
2. **Integration test**: Run full pipeline with mix of valid and invalid files
3. **Edge cases**: Test with all files failed, no files failed, partial failures
4. **Output validation**: Ensure output file structure is maintained even with failures

This template ensures consistent, robust error handling across all search methods while maintaining pipeline functionality even when individual files fail.