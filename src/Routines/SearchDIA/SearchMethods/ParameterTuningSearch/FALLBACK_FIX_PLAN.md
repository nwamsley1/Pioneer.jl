# Fix Plan for ParameterTuningSearch Fallback RT Model Issue

## Problem Analysis

When ParameterTuningSearch fails to collect enough PSMs and uses fallback parameters, it doesn't properly store RT models in SearchContext. This causes downstream methods (like FirstPassSearch) to fail with:
```
KeyError: key 1 not found
```
when trying to access RT models via `getRtTolFunc(search_context, ms_file_idx)`.

## Root Cause

In `process_file!`, when using fallback parameters (lines 341-357), the code:
1. Sets conservative mass error model in SearchContext ✓
2. Creates identity RT model in results ✓
3. **FAILS to store RT model in SearchContext** ✗

The missing step causes `search_context.irt_to_rt_model` dictionary to not have an entry for the file index.

## Current Fallback Code (Lines 341-357)

```julia
# Use conservative fallback parameters
@warn "Failed to converge mass error model after $n_attempts attempts for file $ms_file_idx. Using conservative defaults."

# Set conservative mass error model
fallback_mass_err = MassErrorModel(0.0f0, (50.0f0, 50.0f0))
results.mass_err_model[] = fallback_mass_err
setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)

# Set identity RT model
results.rt_to_irt_model[] = SplineRtConversionModel(identity)
results.irt = Float32[]
results.rt = Float32[]

# Track that we used fallback
push!(warnings, "FALLBACK: Used conservative parameters")
fallback_used = true
```

## Enhanced Solution with Cross-File Parameter Borrowing

### 1. Primary Fix: Store RT Model in SearchContext
```julia
# CRITICAL: Store RT model in SearchContext for downstream methods
setRtToIrtSpline!(search_context, ms_file_idx, identity_rt_model)
```

### 2. Enhanced Fallback: Borrow Parameters from Neighboring Files

Instead of using generic conservative parameters, attempt to borrow from successfully tuned files:

```julia
function get_fallback_parameters(search_context::SearchContext, ms_file_idx::Int, results::ParameterTuningSearchResults)
    # Try to find a successfully tuned file's parameters
    borrowed_from = nothing
    fallback_mass_err = nothing
    fallback_rt_model = nothing
    
    # Check previous files first
    for file_idx in (ms_file_idx-1):-1:1
        if hasRtToIrtModel(search_context, file_idx) && hasMassErrorModel(search_context, file_idx)
            borrowed_from = file_idx
            fallback_mass_err = getMassErrorModel(search_context, file_idx)
            fallback_rt_model = getRtToIrtModel(search_context, file_idx)
            break
        end
    end
    
    # If no previous file, check next files
    if borrowed_from === nothing
        for file_idx in (ms_file_idx+1):getNumFiles(search_context)
            if hasRtToIrtModel(search_context, file_idx) && hasMassErrorModel(search_context, file_idx)
                borrowed_from = file_idx
                fallback_mass_err = getMassErrorModel(search_context, file_idx)
                fallback_rt_model = getRtToIrtModel(search_context, file_idx)
                break
            end
        end
    end
    
    # If found parameters from another file
    if borrowed_from !== nothing
        @warn "Using parameters from file $borrowed_from for file $ms_file_idx due to insufficient PSMs"
        push!(results.warnings, "BORROWED: Using parameters from file $borrowed_from")
        return fallback_mass_err, fallback_rt_model, borrowed_from
    else
        # Use conservative defaults if no other files available
        @warn "No successfully tuned files available. Using conservative defaults for file $ms_file_idx"
        push!(results.warnings, "FALLBACK: Used conservative default parameters")
        return MassErrorModel(0.0f0, (50.0f0, 50.0f0)), SplineRtConversionModel(identity), nothing
    end
end

# In process_file! fallback section:
fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(search_context, ms_file_idx, results)

# Set models in both results and SearchContext
results.mass_err_model[] = fallback_mass_err
setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)

results.rt_to_irt_model[] = fallback_rt_model
setRtToIrtSpline!(search_context, ms_file_idx, fallback_rt_model)

# Track borrowing in results for reporting
if borrowed_from !== nothing
    results.borrowed_from = borrowed_from
end
```

### 3. Always Generate QC Plots

Even when insufficient data for model fitting, generate plots showing:
- The limited data points available
- The fallback/borrowed model overlay
- Clear indication of fallback/borrowed status

```julia
function generate_fallback_qc_plots(
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int,
    psms::DataFrame,
    borrowed_from::Union{Nothing, Int}
)
    # Generate RT alignment plot even with limited data
    rt_plot_path = joinpath(params.qc_plot_folder, "rt_alignment", "file_$(ms_file_idx)_rt_alignment.pdf")
    
    # Plot whatever PSMs we have
    if size(psms, 1) > 0
        plot_rt_alignment(
            psms.rt,
            psms.irt,
            results.rt_to_irt_model[],
            rt_plot_path,
            title = borrowed_from !== nothing ? 
                "RT Alignment (Borrowed from File $borrowed_from)" : 
                "RT Alignment (Fallback - Insufficient Data)",
            show_warning = true
        )
    else
        # Create empty plot with warning message
        create_empty_plot(
            rt_plot_path,
            "No PSMs Available for RT Alignment",
            "Using $(borrowed_from !== nothing ? "borrowed" : "default") parameters"
        )
    end
    
    # Similar for mass error plot
    mass_error_plot_path = joinpath(params.qc_plot_folder, "mass_error", "file_$(ms_file_idx)_mass_error.pdf")
    
    if length(results.ppm_errs[]) > 0
        plot_mass_error_distribution(
            results.ppm_errs[],
            results.mass_err_model[],
            mass_error_plot_path,
            title = borrowed_from !== nothing ?
                "Mass Error (Borrowed from File $borrowed_from)" :
                "Mass Error (Fallback - Insufficient Data)",
            show_warning = true
        )
    else
        create_empty_plot(
            mass_error_plot_path,
            "No Fragment Matches for Mass Error",
            "Using $(borrowed_from !== nothing ? "borrowed" : "default") tolerances"
        )
    end
end
```

## Implementation Steps

1. **Add helper functions**:
   - `get_fallback_parameters()` - Find and borrow parameters from other files
   - `generate_fallback_qc_plots()` - Create plots even with limited data
   - `hasRtToIrtModel()`, `hasMassErrorModel()` - Check if models exist

2. **Update fallback logic in process_file!**:
   - Try to borrow from neighboring files first
   - Always store models in SearchContext
   - Always generate QC plots

3. **Enhance reporting**:
   - Track which file parameters were borrowed from
   - Include in summary report
   - Show in QC plots

4. **Update cross-run learning**:
   - Don't include borrowed parameters in learning
   - Mark files that used borrowed parameters

## Testing Strategy

1. **Test parameter borrowing**:
   ```julia
   # Scenario 1: File 1 fails, File 2 succeeds
   # Expected: File 1 borrows from File 2
   
   # Scenario 2: File 2 fails, Files 1 & 3 succeed  
   # Expected: File 2 borrows from File 1 (previous preferred)
   
   # Scenario 3: All files fail
   # Expected: All use conservative defaults
   ```

2. **Test QC plot generation**:
   - Verify plots generated even with 0 PSMs
   - Check borrowed/fallback status shown in titles
   - Ensure warning indicators visible

3. **Test pipeline continuation**:
   - Confirm FirstPassSearch doesn't crash
   - Verify borrowed parameters work correctly
   - Check results quality with borrowed vs default

## Expected Behavior After Fix

### When Insufficient PSMs with Other Successful Files:
1. Warning: "Using parameters from file 2 for file 1 due to insufficient PSMs"
2. Parameters borrowed from nearest successful file
3. QC plots generated showing borrowed status
4. Pipeline continues without errors

### When All Files Have Insufficient PSMs:
1. Warning: "No successfully tuned files available. Using conservative defaults"
2. Conservative defaults used (±50 ppm, identity RT)
3. QC plots generated showing fallback status
4. Pipeline continues without errors

## Enhanced Reporting

### Parameter Tuning Report Should Include:
```
File 1: BORROWED from File 2
  - Mass tolerance: ±15.3 ppm (borrowed)
  - RT model: Spline (borrowed)
  - PSMs collected: 210/1000
  - Warnings: Insufficient PSMs after 3 attempts

File 2: SUCCESS
  - Mass tolerance: ±15.3 ppm
  - RT model: Spline (5 knots)
  - PSMs collected: 2500/1000
```

### Cross-Run Parameter Report:
```
Successfully tuned: 1/2 files
Borrowed parameters: 1 file
Fallback defaults: 0 files

Parameter borrowing:
  File 1 <- File 2
```

## Code Locations

1. **Main fallback logic**: 
   - File: `ParameterTuningSearch.jl`
   - Lines: ~341-357

2. **New helper functions**:
   - Add to: `ParameterTuningSearch.jl` or `utils.jl`
   - Functions: `get_fallback_parameters()`, `generate_fallback_qc_plots()`

3. **QC plot generation**:
   - File: `qc_plots.jl` (may need creation)
   - Functions: `plot_rt_alignment()`, `plot_mass_error_distribution()`, `create_empty_plot()`

4. **Reporting enhancements**:
   - File: `summarize.jl`
   - Update: `generate_parameter_report()`, `generate_cross_run_report()`

## Risk Assessment

- **Low Risk**: Parameter borrowing uses proven models from same experiment
- **Medium Complexity**: Requires tracking borrowed relationships
- **High Value**: Significantly improves robustness for challenging data
- **User Transparency**: Clear warnings and plot annotations

## Benefits

1. **Better parameter estimates** than generic defaults
2. **Improved consistency** across files in same experiment  
3. **Always generates QC plots** for diagnostics
4. **Clear tracking** of parameter sources
5. **Graceful degradation** when data quality varies