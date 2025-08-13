# Fix Plan for generate_summary_report Error

## Problem
The function `generate_summary_report` is called in `summarize_results!` but is not defined anywhere. This causes an UndefVarError when the function tries to generate a summary report after merging PDFs.

## Investigation Needed
1. Check if `generate_summary_report` exists elsewhere in the codebase
2. Determine what this function should do
3. Either find the function or implement it or remove the call

## Analysis

### Current Usage
In `summarize_results!` at line 789:
```julia
# Generate summary report
generate_summary_report(results, search_context, qc_plots_folder)
```

### Expected Functionality
Based on the context and parameters, this function likely should:
1. Create a summary text/PDF report of the parameter tuning results
2. Include statistics about convergence, parameters found, warnings, etc.
3. Save it to the qc_plots_folder

## Solution Options

### Option 1: Remove the Call (Quick Fix)
Simply comment out or remove the call if the summary report is not essential:
```julia
# Generate summary report (TODO: implement if needed)
# generate_summary_report(results, search_context, qc_plots_folder)
```

### Option 2: Implement Basic Summary Report
Create a simple function that writes a summary to a text file:
```julia
function generate_summary_report(results::ParameterTuningSearchResults, 
                                search_context::SearchContext, 
                                output_folder::String)
    summary_path = joinpath(output_folder, "parameter_tuning_summary.txt")
    
    open(summary_path, "w") do io
        println(io, "Parameter Tuning Summary Report")
        println(io, "=" ^ 50)
        println(io, "Generated: $(Dates.now())")
        println(io, "")
        
        # Get diagnostics
        diagnostics = getDiagnostics(results)
        file_statuses = values(diagnostics.file_statuses)
        
        # Overall statistics
        println(io, "Overall Statistics:")
        println(io, "  Total files: $(length(file_statuses))")
        println(io, "  Converged: $(sum(s.converged for s in file_statuses))")
        println(io, "  Used fallback: $(sum(s.used_fallback for s in file_statuses))")
        println(io, "")
        
        # File-by-file results
        println(io, "File-by-File Results:")
        println(io, "-" ^ 50)
        for status in sort(collect(file_statuses), by=s->s.file_idx)
            println(io, "File: $(status.file_name) (index: $(status.file_idx))")
            println(io, "  Converged: $(status.converged)")
            println(io, "  Iterations: $(status.n_iterations)")
            println(io, "  PSMs: $(status.psm_count)")
            println(io, "  Mass offset: $(round(status.mass_offset, digits=2)) ppm")
            println(io, "  Mass tolerance: ($(round(status.mass_tolerance[1], digits=1)), $(round(status.mass_tolerance[2], digits=1))) ppm")
            if !isempty(status.warnings)
                println(io, "  Warnings: $(join(status.warnings, "; "))")
            end
            println(io, "")
        end
        
        # Parameter history statistics
        history = getParameterHistory(results)
        if length(history.file_parameters) > 0
            println(io, "Parameter Statistics:")
            println(io, "-" ^ 50)
            
            converged_params = [p for p in values(history.file_parameters) if p.converged]
            if !isempty(converged_params)
                offsets = [p.mass_offset for p in converged_params]
                left_tols = [p.mass_tolerance[1] for p in converged_params]
                right_tols = [p.mass_tolerance[2] for p in converged_params]
                
                println(io, "Mass Offset (converged files):")
                println(io, "  Median: $(round(median(offsets), digits=2)) ppm")
                println(io, "  Range: [$(round(minimum(offsets), digits=2)), $(round(maximum(offsets), digits=2))] ppm")
                println(io, "")
                
                println(io, "Mass Tolerance (converged files):")
                println(io, "  Left: median=$(round(median(left_tols), digits=1)), range=[$(round(minimum(left_tols), digits=1)), $(round(maximum(left_tols), digits=1))] ppm")
                println(io, "  Right: median=$(round(median(right_tols), digits=1)), range=[$(round(minimum(right_tols), digits=1)), $(round(maximum(right_tols), digits=1))] ppm")
            end
        end
    end
    
    @info "Generated parameter tuning summary report: $summary_path"
end
```

### Option 3: Check for Existing Implementation
Search the codebase for similar report generation functions that might be used as a template or might be the intended function with a different name.

## Recommended Solution

**Go with Option 2** - Implement a basic text summary report. This provides:
1. Immediate fix for the error
2. Useful diagnostic information
3. Can be enhanced later if needed
4. Minimal complexity

## Implementation Steps

1. Add the `generate_summary_report` function to ParameterTuningSearch.jl or utils.jl
2. Import `Dates` if not already imported
3. Import `Statistics` for median calculation if needed
4. Test with a parameter tuning run

## Alternative Quick Fix

If time is critical, use Option 1 and comment out the call:
```julia
# TODO: Implement summary report generation
# generate_summary_report(results, search_context, qc_plots_folder)
```

This eliminates the error immediately while preserving the intent for future implementation.