# Simplified Convergence Strategy for ParameterTuningSearch

## Overview
Remove the confusing `execute_strategy` function and implement a clearer convergence loop directly in the main flow.

## Core Concepts

### Phase Structure (3 phases total)
1. **Phase 1**: Start with zero bias shift
2. **Phase 2**: Start with positive bias shift (+max_tolerance)
3. **Phase 3**: Start with negative bias shift (-max_tolerance)

### Within Each Phase
- Run up to `iterations_per_phase` attempts
- Each iteration expands the mass tolerance
- Each iteration tries to find convergence through bias adjustment

## Detailed Pseudo-Code

```pseudocode
function run_single_phase_attempt(
    phase, scan_count, iteration_state, results, params, search_context, ms_file_idx, spectra
)
    settings = getIterationSettings(params)
    
    # Set up phase
    iteration_state.current_phase = phase
    iteration_state.current_iteration_in_phase = 0
    
    # Apply phase-specific initial bias shift
    initial_bias_shift = get_phase_bias_shift(phase, settings)  # 0, +max_tol, -max_tol
    apply_initial_bias_shift!(search_context, ms_file_idx, initial_bias_shift)
    
    # Reset mass tolerance to initial value for this phase
    reset_mass_tolerance!(search_context, ms_file_idx, settings.init_mass_tol_ppm)
    
    # Create filtered spectra with fixed scan count for entire phase
    filtered_spectra = FilteredMassSpecData(
        spectra,
        max_scans = scan_count,
        topn = getTopNPeaks(params),
        target_ms_order = 2
    )
    
    # Run iterations within phase
    for iter in 1:settings.iterations_per_phase
        iteration_state.current_iteration_in_phase = iter
        iteration_state.total_iterations += 1
        
        log_iteration_start(phase, iter, scan_count)
        
        # Get current tolerance for this iteration
        current_tolerance = calculate_iteration_tolerance(
            settings.init_mass_tol_ppm,
            settings.mass_tolerance_scale_factor,
            iter
        )
        # Note: current_tolerance = init_tol * (scale_factor ^ (iter - 1))
        
        # Step 1: Collect PSMs at current tolerance (no bias adjustment yet)
        psms_initial = collect_psms(
            filtered_spectra, spectra, search_context, params, ms_file_idx
        )
        
        # Step 2: Fit models to get initial convergence check
        mass_err_model_initial, rt_model, ppm_errs = fit_models_from_psms(
            psms_initial, spectra, search_context, params, ms_file_idx
        )
        
        # Step 3: Check convergence with initial collection
        if check_convergence(mass_err_model_initial, params)
            # Test tolerance expansion post-convergence
            expanded_result = test_tolerance_expansion!(
                search_context, params, ms_file_idx,
                psms_initial, mass_err_model_initial, ppm_errs,
                current_tolerance, filtered_spectra, spectra
            )
            
            if expanded_result.improved
                store_results!(results, expanded_result)
                log_convergence_success("Initial collection (expanded)", phase, iter)
                return true
            else
                store_results!(results, psms_initial, mass_err_model_initial, ppm_errs)
                log_convergence_success("Initial collection", phase, iter)
                return true
            end
        end
        
        # Step 4: If not converged, expand tolerance for bias adjustment attempt
        if iter < settings.iterations_per_phase  # Don't expand on last iteration
            expand_mass_tolerance!(
                search_context, ms_file_idx, params,
                settings.mass_tolerance_scale_factor
            )
            # Now tolerance is: init_tol * (scale_factor ^ iter)
        end
        
        # Step 5: Adjust mass bias based on initial fitting
        if mass_err_model_initial !== nothing
            new_bias = getMassOffset(mass_err_model_initial)
            adjust_mass_bias!(search_context, ms_file_idx, new_bias)
            log_bias_adjustment(old_bias, new_bias)
        end
        
        # Step 6: Collect PSMs with adjusted bias (and expanded tolerance)
        psms_adjusted = collect_psms(
            filtered_spectra, spectra, search_context, params, ms_file_idx
        )
        
        # Step 7: Fit models again with bias-adjusted PSMs
        mass_err_model_adjusted, rt_model, ppm_errs = fit_models_from_psms(
            psms_adjusted, spectra, search_context, params, ms_file_idx
        )
        
        # Step 8: Check convergence with bias-adjusted collection
        if check_convergence(mass_err_model_adjusted, params)
            # Test tolerance expansion post-convergence
            expanded_result = test_tolerance_expansion!(
                search_context, params, ms_file_idx,
                psms_adjusted, mass_err_model_adjusted, ppm_errs,
                current_tolerance * settings.mass_tolerance_scale_factor,
                filtered_spectra, spectra
            )
            
            if expanded_result.improved
                store_results!(results, expanded_result)
                log_convergence_success("Bias-adjusted (expanded)", phase, iter)
                return true
            else
                store_results!(results, psms_adjusted, mass_err_model_adjusted, ppm_errs)
                log_convergence_success("Bias-adjusted", phase, iter)
                return true
            end
        end
        
        log_iteration_complete(phase, iter, "No convergence")
    end
    
    log_phase_complete(phase, "No convergence after all iterations")
    return false
end
```

## Helper Functions

```pseudocode
function get_phase_bias_shift(phase, settings)
    # Calculate the maximum tolerance that would be reached
    max_tolerance = settings.init_mass_tol_ppm * 
                   (settings.mass_tolerance_scale_factor ^ settings.iterations_per_phase)
    
    if phase == 1
        return 0.0  # No bias shift for phase 1
    elseif phase == 2
        return +max_tolerance  # Positive shift for phase 2
    elseif phase == 3
        return -max_tolerance  # Negative shift for phase 3
    end
end

function calculate_iteration_tolerance(init_tol, scale_factor, iteration)
    # Tolerance grows exponentially with iteration
    return init_tol * (scale_factor ^ (iteration - 1))
end

function check_convergence(mass_err_model, params)
    if mass_err_model === nothing
        return false
    end
    
    # Check if the fitted parameters are reasonable
    mass_offset = getMassOffset(mass_err_model)
    left_tol = getLeftTol(mass_err_model)
    right_tol = getRightTol(mass_err_model)
    
    # Convergence criteria:
    # 1. Mass offset is reasonably small (within initial tolerance)
    # 2. Fitted tolerances are reasonable (not too wide)
    
    init_tol = params.iteration_settings.init_mass_tol_ppm
    
    # Mass offset should be less than 1/4 of initial tolerance
    if abs(mass_offset) > (init_tol / 4)
        return false
    end
    
    # Average fitted tolerance should be less than initial tolerance
    avg_fitted_tol = (left_tol + right_tol) / 2
    if avg_fitted_tol > init_tol
        return false
    end
    
    return true
end
```

## Main Scan Scaling Loop

```pseudocode
function process_file!(results, params, search_context, ms_file_idx, spectra)
    iteration_state = IterationState()
    settings = getIterationSettings(params)
    
    scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    
    # Initialize models with defaults
    initialize_models!(search_context, ms_file_idx, params)
    
    # Main scan scaling loop
    attempt = 0
    while scan_count <= max_scans
        attempt += 1
        iteration_state.scan_attempt = attempt
        
        log_scan_attempt_start(attempt, scan_count)
        
        # Run all 3 phases with current scan count
        for phase in 1:3
            log_phase_start(phase, attempt)
            
            converged = run_single_phase_attempt(
                phase, scan_count, iteration_state,
                results, params, search_context, ms_file_idx, spectra
            )
            
            if converged
                log_success(attempt, phase, scan_count)
                return results
            end
        end
        
        # Scale up scan count for next attempt
        if scan_count >= max_scans
            log_max_scans_reached()
            break
        end
        
        next_scan_count = ceil(scan_count * settings.scan_scale_factor)
        scan_count = min(next_scan_count, max_scans)
        
        log_scan_scaling(scan_count)
    end
    
    # Apply fallback strategy if no convergence
    apply_fallback_strategy!(results, search_context, params, ms_file_idx)
    
    return results
end
```

## Key Changes from Current Implementation

1. **Remove `execute_strategy` function** - All logic is now inline and clear
2. **Explicit tolerance management** - Tolerance is calculated based on iteration number
3. **Clear bias adjustment flow** - Always try initial collection, then bias adjustment
4. **Simplified convergence checking** - Single function with clear criteria
5. **Phase initialization** - Each phase starts fresh with its specific bias shift
6. **Post-convergence expansion** - Always test if expanding tolerance improves results

## Iteration Flow Summary

For each phase and iteration:
1. **Initial Collection**: Collect PSMs at current tolerance
2. **Initial Fit**: Fit models and check convergence
3. **Expand Tolerance**: Increase tolerance for next attempts (if not converged)
4. **Adjust Bias**: Update bias based on fitted model
5. **Adjusted Collection**: Collect PSMs with new bias
6. **Adjusted Fit**: Fit models and check convergence
7. **Post-Convergence Test**: If converged, test with expanded tolerance

## Benefits

1. **Clarity**: Each step has a clear purpose
2. **No Hidden State**: All tolerance and bias changes are explicit
3. **Predictable**: Easy to trace through the logic
4. **Testable**: Each component can be tested independently
5. **Maintainable**: Clear separation of concerns

## Configuration Example

```json
{
    "iteration_settings": {
        "init_mass_tol_ppm": 20.0,
        "mass_tolerance_scale_factor": 1.5,
        "iterations_per_phase": 3,
        "scan_scale_factor": 2.0
    }
}
```

With these settings:
- Iteration 1: 20 ppm
- Iteration 2: 30 ppm (20 * 1.5)
- Iteration 3: 45 ppm (20 * 1.5^2)
- Phase 2 starts with +45 ppm bias
- Phase 3 starts with -45 ppm bias