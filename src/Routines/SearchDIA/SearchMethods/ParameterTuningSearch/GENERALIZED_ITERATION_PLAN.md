# Plan: Generalized Iteration Strategy with Configurable Scaling and Reset

## Status: ✅ IMPLEMENTED (2025-01-14)

## Overview

This plan generalizes the parameter tuning iteration strategy to support:
1. **Configurable mass tolerance scaling** - User-defined factor for tolerance expansion (not just doubling)
2. **Phase-based reset mechanism** - After N iterations, reset to initial tolerance with bias shift
3. **Bias shift exploration** - Try positive and negative bias offsets to handle extreme mass errors

## Current Behavior

### Existing Implementation
- **Scan Growth**: Doubles each iteration (500→1000→2000→4000→8000)
- **Mass Tolerance**: Doubles when Strategy 2 is executed (hardcoded 2.0x factor)
- **Max Iterations**: Fixed at 5 iterations
- **No Reset**: Once max iterations reached, uses fallback parameters

### Problems with Current Approach
1. **Fixed scaling factor** - Always doubles, may be too aggressive or conservative
2. **No bias exploration** - If true bias is ±40 ppm but initial window is ±20 ppm, fails to find PSMs
3. **Single-phase only** - No mechanism to try different parameter spaces
4. **Limited adaptability** - Cannot adjust strategy based on data characteristics

## Proposed Solution

### New Configuration Parameters

Add to JSON under `parameter_tuning.search_settings`:
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 1.5,      // Factor to scale tolerance each iteration (default: 2.0)
    "iterations_per_phase": 3,                // Iterations before reset (default: 3)
    "max_phases": 3,                          // Maximum number of phases (default: 3)
    "bias_shift_strategy": "alternating",    // How to shift bias: "alternating", "positive_first", "negative_first"
    "bias_shift_magnitude": "max_tolerance"  // Amount to shift: "max_tolerance" or specific value
}
```

### Phase-Based Iteration Structure

```
Phase 1 (iterations 1-3): Initial parameters
  Iteration 1: Initial tolerance, zero bias
  Iteration 2: Scaled tolerance (1.5x), adjusted bias
  Iteration 3: Scaled tolerance (2.25x), adjusted bias
  
Phase 2 (iterations 4-6): Reset with positive bias shift
  Iteration 4: Initial tolerance, +max_tolerance bias
  Iteration 5: Scaled tolerance (1.5x), adjusted bias
  Iteration 6: Scaled tolerance (2.25x), adjusted bias
  
Phase 3 (iterations 7-9): Reset with negative bias shift
  Iteration 7: Initial tolerance, -max_tolerance bias
  Iteration 8: Scaled tolerance (1.5x), adjusted bias
  Iteration 9: Scaled tolerance (2.25x), adjusted bias
```

## Implementation Details

### 1. New Type Definitions

#### Add to `types.jl`:
```julia
"""
Configuration for iteration behavior in parameter tuning.
"""
struct IterationSettings
    mass_tolerance_scale_factor::Float32
    iterations_per_phase::Int64
    max_phases::Int64
    bias_shift_strategy::Symbol  # :alternating, :positive_first, :negative_first
    bias_shift_magnitude::Union{Float32, Symbol}  # Number or :max_tolerance
end

# Add to ParameterTuningSearchParameters
struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # ... existing fields ...
    iteration_settings::IterationSettings
end
```

### 2. Iteration State Tracking

#### Add iteration state structure:
```julia
"""
Tracks state across phases and iterations.
"""
mutable struct IterationState
    current_phase::Int64
    current_iteration_in_phase::Int64
    total_iterations::Int64
    phase_bias_shifts::Vector{Float32}  # Bias shift applied at start of each phase
    converged::Bool
end

function next_iteration!(state::IterationState, settings::IterationSettings)
    state.total_iterations += 1
    state.current_iteration_in_phase += 1
    
    if state.current_iteration_in_phase > settings.iterations_per_phase
        state.current_phase += 1
        state.current_iteration_in_phase = 1
    end
    
    return state.current_phase <= settings.max_phases
end
```

### 3. Modified Helper Functions

#### Update `expand_mass_tolerance!`:
```julia
function expand_mass_tolerance!(search_context, ms_file_idx, params, scale_factor::Float32)
    current_model = getMassErrorModel(search_context, ms_file_idx)
    current_left = getLeftTol(current_model)
    current_right = getRightTol(current_model)
    
    new_left = current_left * scale_factor
    new_right = current_right * scale_factor
    
    new_model = create_capped_mass_model(
        getMassOffset(current_model),
        new_left,
        new_right,
        getMaxTolerancePpm(params)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
    
    @info "Scaled tolerance by factor $(scale_factor): " *
          "from ($(round(current_left, digits=1)), $(round(current_right, digits=1))) " *
          "to ($(round(getLeftTol(new_model), digits=1)), $(round(getRightTol(new_model), digits=1))) ppm"
end
```

#### New function `reset_for_new_phase!`:
```julia
function reset_for_new_phase!(search_context, ms_file_idx, params, phase::Int64, iteration_state::IterationState)
    settings = params.iteration_settings
    initial_tolerance = getFragTolPpm(params)
    
    # Determine bias shift for this phase
    bias_shift = calculate_phase_bias_shift(phase, settings, params)
    push!(iteration_state.phase_bias_shifts, bias_shift)
    
    # Reset to initial tolerance with new bias
    new_model = create_capped_mass_model(
        bias_shift,
        initial_tolerance,
        initial_tolerance,
        getMaxTolerancePpm(params)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
    
    @info "Phase $phase: Reset to initial tolerance with bias shift $(round(bias_shift, digits=1)) ppm"
end

function calculate_phase_bias_shift(phase::Int64, settings::IterationSettings, params)::Float32
    if phase == 1
        return 0.0f0  # No shift in first phase
    end
    
    # Determine magnitude
    magnitude = if settings.bias_shift_magnitude == :max_tolerance
        getMaxTolerancePpm(params)
    else
        Float32(settings.bias_shift_magnitude)
    end
    
    # Determine direction based on strategy
    if settings.bias_shift_strategy == :alternating
        return phase % 2 == 0 ? magnitude : -magnitude
    elseif settings.bias_shift_strategy == :positive_first
        return phase == 2 ? magnitude : -magnitude
    else  # :negative_first
        return phase == 2 ? -magnitude : magnitude
    end
end
```

### 4. Updated Main Loop in `process_file!`

```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    warnings = String[]
    final_psm_count = 0
    
    # Initialize iteration state
    iteration_state = IterationState(1, 0, 0, Float32[], false)
    settings = params.iteration_settings
    
    try
        @info "Processing file: $parsed_fname (index: $ms_file_idx)"
        @info "Iteration settings: scale_factor=$(settings.mass_tolerance_scale_factor), " *
              "iterations_per_phase=$(settings.iterations_per_phase), max_phases=$(settings.max_phases)"
        
        # Initialize models
        initialize_models!(search_context, ms_file_idx, params)
        
        # Initialize filtered spectra
        filtered_spectra = initialize_filtered_spectra(spectra, params)
        
        # Main convergence loop
        while next_iteration!(iteration_state, settings)
            # Check if we're starting a new phase
            if iteration_state.current_iteration_in_phase == 1 && iteration_state.current_phase > 1
                # Reset for new phase
                reset_for_new_phase!(search_context, ms_file_idx, params, 
                                   iteration_state.current_phase, iteration_state)
                
                # Reset filtered spectra for new phase
                filtered_spectra = initialize_filtered_spectra(spectra, params)
            end
            
            # Add scans for this iteration
            if !add_scans_for_iteration!(filtered_spectra, params, 
                                        iteration_state.total_iterations - 1)
                @info "Reached maximum scan count"
                continue  # Move to next phase if we can't add more scans
            end
            
            @info "Phase $(iteration_state.current_phase), " *
                  "Iteration $(iteration_state.current_iteration_in_phase) " *
                  "(Total: $(iteration_state.total_iterations))"
            
            # Strategy 1: Try with current parameters
            psms, mass_err_model, ppm_errs = execute_strategy(
                1, filtered_spectra, spectra, 
                search_context, params, ms_file_idx, results
            )
            
            # Check convergence after Strategy 1
            if check_and_store_convergence!(
                results, search_context, params, ms_file_idx,
                psms, mass_err_model, ppm_errs, 
                "Strategy 1 (Phase $(iteration_state.current_phase))"
            )
                converged = true
                iteration_state.converged = true
                final_psm_count = psms !== nothing ? size(psms, 1) : 0
                break
            else
                # Strategy 2: Scale tolerance AND adjust bias
                # Calculate cumulative scale factor for this iteration within phase
                cumulative_scale = settings.mass_tolerance_scale_factor ^ 
                                  iteration_state.current_iteration_in_phase
                
                # Only expand if we're past the first iteration of the phase
                if iteration_state.current_iteration_in_phase > 1
                    expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                          settings.mass_tolerance_scale_factor)
                end
                
                # Now execute strategy 2 with the scaled tolerance
                psms, mass_err_model, ppm_errs = execute_modified_strategy_2(
                    filtered_spectra, spectra, search_context, 
                    params, ms_file_idx, results
                )
                
                # Check convergence after Strategy 2
                if check_and_store_convergence!(
                    results, search_context, params, ms_file_idx,
                    psms, mass_err_model, ppm_errs, 
                    "Strategy 2 (Phase $(iteration_state.current_phase), " *
                    "scale=$(round(cumulative_scale, digits=2))x)"
                )
                    converged = true
                    iteration_state.converged = true
                    final_psm_count = psms !== nothing ? size(psms, 1) : 0
                    break
                end
            end
            
            @info "Phase $(iteration_state.current_phase), " *
                  "Iteration $(iteration_state.current_iteration_in_phase) complete"
        end
        
        # Log phase summary
        @info "Completed $(iteration_state.total_iterations) total iterations across " *
              "$(iteration_state.current_phase) phases. Converged: $converged"
        
        if !isempty(iteration_state.phase_bias_shifts)
            @info "Phase bias shifts attempted: $(iteration_state.phase_bias_shifts)"
        end
        
        # Store results
        store_final_results!(results, search_context, params, ms_file_idx, 
                           converged, iteration_state.total_iterations, 
                           final_psm_count, warnings)
        
        # Handle non-convergence
        if !converged
            push!(warnings, "Failed after $(iteration_state.current_phase) phases")
            handle_non_convergence!(results, search_context, ms_file_idx, 
                                   iteration_state.total_iterations, warnings)
        end
        
    catch e
        handle_error!(results, search_context, ms_file_idx, e, parsed_fname)
    end
    
    return results
end
```

### 5. Modified Strategy 2

Since we're now using configurable scaling, we need a modified execute_strategy that doesn't hardcode the doubling:

```julia
function execute_modified_strategy_2(filtered_spectra, spectra, search_context, 
                                    params, ms_file_idx, results)
    @info "Strategy 2: Adjusting bias with current tolerance"
    
    # Collect PSMs with current (already scaled) tolerance to estimate bias
    psms_temp, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                       params, ms_file_idx, 
                                       "with current tolerance for bias estimation")
    
    # Fit mass error model to get bias estimate
    mass_err_temp, _, _ = fit_models_from_psms(psms_temp, spectra, search_context, 
                                               params, ms_file_idx)
    
    if mass_err_temp === nothing
        @info "No valid mass error model for bias adjustment, using current tolerance only"
        psms = psms_temp
    else
        # Adjust bias based on the fitted model
        adjust_mass_bias!(search_context, ms_file_idx, mass_err_temp, params)
        
        # Collect final PSMs with adjusted bias
        psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                       params, ms_file_idx, 
                                       "with current tolerance and adjusted bias")
    end
    
    # Fit models from collected PSMs
    mass_err_model, ppm_errs, psm_count = fit_models_from_psms(psms, spectra, 
                                                                search_context, 
                                                                params, ms_file_idx)
    
    return psms, mass_err_model, ppm_errs
end
```

## Benefits of This Approach

### 1. **Adaptive Exploration**
- Systematically explores different parameter spaces
- Handles cases where true bias is far from zero
- Configurable scaling allows fine-tuning for different instruments

### 2. **Better Coverage**
- Phase 1: Explores around zero bias
- Phase 2: Explores positive bias region
- Phase 3: Explores negative bias region
- Ensures all reasonable parameter spaces are tested

### 3. **Flexibility**
- Users can configure based on their instrument characteristics
- Scale factor can be conservative (1.2) or aggressive (3.0)
- Number of phases and iterations per phase are adjustable

### 4. **Graceful Degradation**
- If one phase fails, moves to the next
- Still falls back to borrowing parameters if all phases fail
- Maintains compatibility with existing fallback mechanisms

## Configuration Examples

### Conservative Settings (High-accuracy instruments)
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 1.2,
    "iterations_per_phase": 4,
    "max_phases": 2,
    "bias_shift_strategy": "alternating",
    "bias_shift_magnitude": 20.0
}
```

### Aggressive Settings (Lower-accuracy or unknown samples)
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 2.5,
    "iterations_per_phase": 2,
    "max_phases": 4,
    "bias_shift_strategy": "alternating",
    "bias_shift_magnitude": "max_tolerance"
}
```

### Default Settings (Backward compatible)
```json
"iteration_settings": {
    "mass_tolerance_scale_factor": 2.0,
    "iterations_per_phase": 3,
    "max_phases": 3,
    "bias_shift_strategy": "alternating",
    "bias_shift_magnitude": "max_tolerance"
}
```

## Implementation Steps

1. ✅ **Add configuration structures** to `types.jl`
2. ✅ **Update parameter extraction** in `ParameterTuningSearchParameters` constructor
3. ✅ **Implement iteration state tracking** structures and functions
4. ✅ **Modify helper functions** to accept scale factors
5. ✅ **Implement phase reset logic** with bias shifting
6. ✅ **Update main loop** in `process_file!`
7. ✅ **Update execute_strategy** to handle modified Strategy 2
8. ✅ **Add logging** for phase transitions and bias shifts
9. ✅ **Update documentation** and example configurations
10. 🔄 **Test with various parameter combinations**

## Implementation Summary

Successfully implemented the generalized iteration strategy with:
- `IterationSettings` struct for configuration
- `IterationState` for tracking progress across phases
- Phase-based reset mechanism with bias shifting
- Configurable tolerance scaling factor
- Comprehensive logging of phase transitions
- Full backward compatibility with default settings

## Testing Strategy

### Unit Tests
1. Test iteration state transitions
2. Test bias shift calculations for different strategies
3. Test tolerance scaling with different factors
4. Test phase reset behavior

### Integration Tests
1. Test with files that have zero bias (should converge in Phase 1)
2. Test with files that have +30 ppm bias (should converge in Phase 2)
3. Test with files that have -30 ppm bias (should converge in Phase 3)
4. Test with impossible parameters (should use fallback after all phases)

### Performance Tests
1. Compare convergence speed with different scale factors
2. Measure overhead of phase resets
3. Validate that early convergence stops iteration

## Backward Compatibility

To maintain backward compatibility:
1. Default values match current behavior (2x scaling, 3 iterations)
2. If `iteration_settings` is not present in JSON, use defaults
3. Old configuration files will continue to work unchanged
4. Phase structure is transparent to downstream methods

## Future Enhancements

1. **Adaptive scaling** - Adjust scale factor based on PSM yield
2. **Smart phase selection** - Skip phases based on initial results
3. **Parameter persistence** - Save successful parameters for similar files
4. **Diagnostic reporting** - Track which phase/iteration succeeded for QC
5. **Machine learning** - Predict optimal parameters based on file characteristics

## Summary

This plan provides a robust, configurable framework for parameter tuning that:
- Handles extreme mass biases through phase-based exploration
- Allows user control over convergence strategy
- Maintains backward compatibility
- Provides clear diagnostic information
- Sets foundation for future enhancements

The implementation is modular, testable, and integrates cleanly with the existing codebase while significantly improving the robustness of parameter tuning for challenging datasets.