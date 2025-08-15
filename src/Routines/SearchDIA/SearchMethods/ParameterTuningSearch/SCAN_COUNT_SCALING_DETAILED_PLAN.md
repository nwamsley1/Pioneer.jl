# Detailed Plan: Simplified Scan Count Scaling Strategy

## Status: DETAILED IMPLEMENTATION PLAN

## Overview

Restructure the parameter tuning process to use a cleaner scan count scaling approach:
- Run all 3 phases with a fixed scan count
- If no convergence, scale up the scan count and repeat all phases
- Continue until reaching maximum scan count
- Maintain per-iteration mass tolerance scaling within each phase

## Current vs Proposed Architecture

### Current Architecture
```
foreach iteration:
    scan_count = min(initial * 2^(iteration-1), max_scans)
    foreach phase (if needed):
        foreach iteration_in_phase:
            mass_tolerance *= scale_factor
            try convergence with current scans and tolerance
```

### Proposed Architecture
```
scan_count = initial_scan_count
while scan_count <= max_scan_count:
    foreach phase (1, 2, 3):
        foreach iteration_in_phase:
            mass_tolerance *= scale_factor  # Still scales within phase
            try convergence with fixed scan_count
    if converged: 
        return success
    scan_count = min(scan_count * scan_scale_factor, max_scan_count)
```

## Detailed Implementation

### 1. Update IterationState Structure

```julia
mutable struct IterationState
    # Phase tracking (unchanged)
    current_phase::Int64
    current_iteration_in_phase::Int64
    total_iterations::Int64
    phase_bias_shifts::Vector{Float32}
    converged::Bool
    collection_tolerance::Float32
    
    # NEW: Scan scaling tracking
    scan_attempt::Int64              # Which attempt (1, 2, 3, ...)
    current_scan_count::Int64        # Scan count for current attempt
    max_scan_count_reached::Bool     # Flag when we hit the maximum
    
    function IterationState()
        new(1, 0, 0, Float32[], false, 0.0f0, 1, 0, false)
    end
end
```

### 2. Add Scan Scale Factor to Parameters

```julia
struct IterationSettings
    init_mass_tol_ppm::Float32
    mass_tolerance_scale_factor::Float32
    iterations_per_phase::Int64
    scan_scale_factor::Float32  # NEW: Factor to scale scan count between attempts
    
    function IterationSettings(
        init_tol::Float32 = 20.0f0,
        mass_scale_factor::Float32 = 2.0f0,
        iterations_per_phase::Int64 = 3,
        scan_scale_factor::Float32 = 2.0f0  # NEW parameter with default
    )
        @assert init_tol > 0.0f0 "Initial tolerance must be positive"
        @assert mass_scale_factor > 1.0f0 "Mass scale factor must be greater than 1"
        @assert iterations_per_phase > 0 "Iterations per phase must be positive"
        @assert scan_scale_factor > 1.0f0 "Scan scale factor must be greater than 1"
        
        new(init_tol, mass_scale_factor, iterations_per_phase, scan_scale_factor)
    end
end
```

### 3. Parameter Extraction Update

In `ParameterTuningSearchParameters` constructor:
```julia
# Require iteration_settings to be present
iteration_settings = tuning_params.iteration_settings
IterationSettings(
    Float32(iter.init_mass_tol_ppm),           # Required
    Float32(iter.mass_tolerance_scale_factor), # Required
    Int64(iter.iterations_per_phase),          # Required
    Float32(iter.scan_scale_factor)            # Required
)
```

### 4. New Core Functions

#### 4.1 run_single_phase_attempt
```julia
"""
    run_single_phase_attempt(phase::Int64, scan_count::Int64, iteration_state, 
                            results, params, search_context, ms_file_idx, spectra)

Run a single phase (3 iterations) with fixed scan count.
Returns true if converged, false otherwise.
"""
function run_single_phase_attempt(
    phase::Int64,
    scan_count::Int64,
    iteration_state::IterationState,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    settings = getIterationSettings(params)
    
    # Set up phase
    iteration_state.current_phase = phase
    iteration_state.current_iteration_in_phase = 0
    
    # Apply phase bias shift
    reset_for_new_phase!(search_context, ms_file_idx, params, phase, iteration_state)
    
    # Create filtered spectra with fixed scan count
    filtered_spectra = FilteredMassSpecData(
        spectra,
        max_scans = scan_count,
        topn = something(getTopNPeaks(params), 200),
        target_ms_order = UInt8(2)
    )
    
    # Run iterations within phase (with mass tolerance scaling)
    for iter in 1:settings.iterations_per_phase
        iteration_state.current_iteration_in_phase = iter
        iteration_state.total_iterations += 1
        
        @info "Scan Attempt $(iteration_state.scan_attempt), " *
              "Phase $phase, Iteration $iter " *
              "(Total: $(iteration_state.total_iterations)) " *
              "with $scan_count scans"
        
        # Try Strategy 1
        psms, mass_err_model, ppm_errs = execute_strategy(
            1, filtered_spectra, spectra, 
            search_context, params, ms_file_idx, iteration_state
        )
        
        if check_and_store_convergence!(
            results, search_context, params, ms_file_idx,
            psms, mass_err_model, ppm_errs, 
            "Strategy 1 (Attempt $(iteration_state.scan_attempt), Phase $phase)",
            iteration_state, filtered_spectra, spectra
        )
            return true
        end
        
        # Expand mass tolerance for next iteration (except first)
        if iter > 1
            expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                  settings.mass_tolerance_scale_factor)
        end
        
        # Try Strategy 2
        psms, mass_err_model, ppm_errs = execute_strategy(
            2, filtered_spectra, spectra,
            search_context, params, ms_file_idx, iteration_state
        )
        
        if check_and_store_convergence!(
            results, search_context, params, ms_file_idx,
            psms, mass_err_model, ppm_errs, 
            "Strategy 2 (Attempt $(iteration_state.scan_attempt), Phase $phase)",
            iteration_state, filtered_spectra, spectra
        )
            return true
        end
        
        # Expand mass tolerance for next iteration
        if iter < settings.iterations_per_phase
            expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                  settings.mass_tolerance_scale_factor)
        end
    end
    
    return false
end
```

#### 4.2 run_all_phases_with_scan_count
```julia
"""
    run_all_phases_with_scan_count(scan_count, iteration_state, results, 
                                   params, search_context, ms_file_idx, spectra)

Run all 3 phases with a fixed scan count.
Returns true if any phase achieves convergence.
"""
function run_all_phases_with_scan_count(
    scan_count::Int64,
    iteration_state::IterationState,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    @info "="^60
    @info "Starting scan attempt $(iteration_state.scan_attempt) with $scan_count scans"
    @info "="^60
    
    # Update iteration state
    iteration_state.current_scan_count = scan_count
    
    # Run all 3 phases
    for phase in 1:3
        @info "-"^40
        @info "Beginning Phase $phase of Attempt $(iteration_state.scan_attempt)"
        @info "-"^40
        
        converged = run_single_phase_attempt(
            phase, scan_count, iteration_state,
            results, params, search_context, ms_file_idx, spectra
        )
        
        if converged
            @info "✓ Converged in Phase $phase of Attempt $(iteration_state.scan_attempt)!"
            return true
        end
        
        @info "Phase $phase completed without convergence"
    end
    
    @info "All 3 phases completed without convergence for $scan_count scans"
    return false
end
```

### 5. Redesigned process_file!

```julia
function process_file!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    
    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    warnings = String[]
    final_psm_count = 0
    
    # Initialize iteration state
    iteration_state = IterationState()
    settings = getIterationSettings(params)
    
    # Get scan count parameters
    scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    scan_scale_factor = settings.scan_scale_factor
    
    try
        @info "Processing file: $parsed_fname (index: $ms_file_idx)"
        @info "Scan scaling strategy: initial=$scan_count, max=$max_scans, scale_factor=$scan_scale_factor"
        @info "Mass tolerance scaling: init=$(settings.init_mass_tol_ppm) ppm, " *
              "scale_factor=$(settings.mass_tolerance_scale_factor), " *
              "iterations_per_phase=$(settings.iterations_per_phase)"
        
        # Initialize models
        initialize_models!(search_context, ms_file_idx, params)
        
        # Main scan scaling loop
        attempt_count = 0
        while scan_count <= max_scans
            attempt_count += 1
            iteration_state.scan_attempt = attempt_count
            
            # Run all phases with current scan count
            converged = run_all_phases_with_scan_count(
                scan_count, iteration_state, results,
                params, search_context, ms_file_idx, spectra
            )
            
            if converged
                iteration_state.converged = true
                @info "SUCCESS: Converged after $attempt_count attempt(s) " *
                      "with $scan_count scans"
                break
            end
            
            # Check if we're at max scans
            if scan_count >= max_scans
                @warn "Reached maximum scan count ($max_scans) without convergence"
                iteration_state.max_scan_count_reached = true
                break
            end
            
            # Scale up scan count for next attempt
            next_scan_count = Int64(ceil(scan_count * scan_scale_factor))
            
            # Cap at maximum
            if next_scan_count > max_scans
                if scan_count < max_scans
                    # Do one final attempt at exactly max_scans
                    scan_count = max_scans
                    @info "Scaling to maximum scan count for final attempt"
                else
                    # Already tried at max, give up
                    break
                end
            else
                scan_count = next_scan_count
                @info "Scaling scan count to $scan_count for next attempt"
            end
        end
        
        # Get final PSM count if converged
        if converged
            final_psm_count = size(results.rt, 1)  # Or however we track this
        end
        
        @info "Completed $(iteration_state.total_iterations) total iterations " *
              "across $(iteration_state.scan_attempt) attempt(s). " *
              "Converged: $converged"
        
    catch e
        @warn "Parameter tuning failed for file $parsed_fname with error" exception=e
        push!(warnings, "ERROR: $(string(e))")
    end
    
    # Store results and handle fallback if needed
    if !converged
        @warn "Failed to converge for file $ms_file_idx after " *
              "$(iteration_state.scan_attempt) attempts"
        
        # Apply fallback strategy
        fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
            search_context, params, ms_file_idx, warnings
        )
        
        setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)
        setRtConversionModel!(search_context, ms_file_idx, fallback_rt_model)
        
        # Update results with fallback
        results.mass_err_model[] = fallback_mass_err
        results.rt_to_irt_model[] = fallback_rt_model
    end
    
    # Record tuning status
    store_final_results!(
        results, search_context, params, ms_file_idx,
        converged, iteration_state.scan_attempt, final_psm_count, warnings
    )
    
    # Add to diagnostics
    record_file_status!(
        results.diagnostics, ms_file_idx, parsed_fname,
        converged, !converged, warnings, iteration_state
    )
end
```

### 6. Remove Old Functions

Functions to remove or deprecate:
- `add_scans_for_iteration!` - No longer needed
- Old iteration logic in current `process_file!`

### 7. Update Logging

Add clear section separators for better readability:
```julia
# In run_all_phases_with_scan_count
@info "="^60
@info "SCAN ATTEMPT $attempt_number: $scan_count scans"
@info "="^60

# Between phases
@info "-"^40
@info "Phase $phase_number"
@info "-"^40
```

## Configuration

### JSON Structure
```json
{
    "parameter_tuning": {
        "search_settings": {
            "initial_scan_count": 500,
            "max_parameter_tuning_scans": 8000,
            "iteration_settings": {
                "init_mass_tol_ppm": 20.0,
                "mass_tolerance_scale_factor": 1.5,
                "iterations_per_phase": 3,
                "scan_scale_factor": 2.0  // NEW: configurable scan scaling
            }
        }
    }
}
```

### Parameter Validation
Add to `paramsChecks.jl`:
```julia
# All iteration_settings fields are now required
check_param(iter_settings, "init_mass_tol_ppm", Real)
check_param(iter_settings, "mass_tolerance_scale_factor", Real)
check_param(iter_settings, "iterations_per_phase", Integer)
check_param(iter_settings, "scan_scale_factor", Real)

# Validate ranges
if iter_settings["mass_tolerance_scale_factor"] <= 1.0
    error("mass_tolerance_scale_factor must be greater than 1.0")
end
if iter_settings["scan_scale_factor"] <= 1.0
    error("scan_scale_factor must be greater than 1.0")
end
```

## Testing Strategy
1. Unit test new functions with mock data
2. Integration test with ecoli dataset
3. Compare convergence rates old vs new approach
4. Verify scan count progression follows expected pattern

## Example Execution Flow

With configuration:
- `initial_scan_count`: 500
- `max_parameter_tuning_scans`: 8000  
- `scan_scale_factor`: 2.0
- `mass_tolerance_scale_factor`: 1.5
- `iterations_per_phase`: 3

```
Attempt 1 (500 scans):
  Phase 1: 
    Iter 1: 20 ppm, bias 0 → No convergence
    Iter 2: 30 ppm, bias adjusted → No convergence
    Iter 3: 45 ppm, bias adjusted → No convergence
  Phase 2:
    Iter 1: 20 ppm, bias +67.5 → No convergence
    Iter 2: 30 ppm, bias adjusted → No convergence
    Iter 3: 45 ppm, bias adjusted → No convergence
  Phase 3:
    Iter 1: 20 ppm, bias -67.5 → No convergence
    Iter 2: 30 ppm, bias adjusted → No convergence
    Iter 3: 45 ppm, bias adjusted → No convergence

Attempt 2 (1000 scans):
  Phase 1:
    Iter 1: 20 ppm, bias 0 → No convergence
    Iter 2: 30 ppm, bias adjusted → CONVERGED!
```

## Benefits

1. **Clearer Structure**: Scan scaling is separate from phase/iteration logic
2. **More Efficient**: Don't waste time with many low-scan iterations
3. **Better Logging**: Clear delineation between attempts
4. **Predictable Resource Usage**: Know exactly how many attempts will be made
5. **Maintains Flexibility**: Still have mass tolerance scaling within phases

## Implementation Order

1. **Phase 1**: Update data structures (IterationState, IterationSettings)
2. **Phase 2**: Implement new core functions (run_single_phase_attempt, run_all_phases_with_scan_count)
3. **Phase 3**: Rewrite process_file! to use new approach
4. **Phase 4**: Remove old functions and clean up
5. **Phase 5**: Update tests and documentation
6. **Phase 6**: Update JSON validation and examples

## Risk Mitigation

- **Risk**: Breaking existing workflows
  - **Mitigation**: Clear error messages when required parameters are missing
  - **Mitigation**: Provide updated example JSON files
  
- **Risk**: Different convergence behavior
  - **Mitigation**: Extensive testing on representative datasets
  
- **Risk**: Increased runtime for easy cases
  - **Mitigation**: Early convergence still exits immediately

## Summary

This plan restructures the parameter tuning to cleanly separate scan count scaling (between attempts) from mass tolerance scaling (within phases). The result is more predictable, easier to understand, and more efficient resource usage while maintaining all the flexibility of the current approach.