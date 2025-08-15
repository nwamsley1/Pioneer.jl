# Simplified Convergence Strategy for ParameterTuningSearch (V2)

## Overview
Remove the confusing `execute_strategy` function and implement a clearer convergence loop directly in the main flow. Filtered spectra is collected once and reused/appended throughout all phases and scan scaling attempts.

## Core Concepts

### Filtered Spectra Management
- Create filtered spectra ONCE at the beginning with initial scan count
- Reuse the same filtered spectra across all 3 phases
- Only append more scans when scaling up between attempts
- Never recreate filtered spectra

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
function run_single_phase(
    phase, filtered_spectra, iteration_state, results, params, search_context, ms_file_idx, spectra
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
    
    # Run iterations within phase (using the SAME filtered_spectra)
    for iter in 1:settings.iterations_per_phase
        iteration_state.current_iteration_in_phase = iter
        iteration_state.total_iterations += 1
        
        log_iteration_start(phase, iter, length(filtered_spectra))
        
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

## Main Process File Function

```pseudocode
function process_file!(results, params, search_context, ms_file_idx, spectra)
    iteration_state = IterationState()
    settings = getIterationSettings(params)
    
    initial_scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    scan_scale_factor = settings.scan_scale_factor
    
    # Initialize models with defaults
    initialize_models!(search_context, ms_file_idx, params)
    
    # Create filtered spectra ONCE with initial scan count
    filtered_spectra = FilteredMassSpecData(
        spectra,
        max_scans = initial_scan_count,
        topn = getTopNPeaks(params),
        target_ms_order = 2
    )
    
    # Track current scan count
    current_scan_count = initial_scan_count
    attempt = 0
    converged = false
    
    # Main scan scaling loop
    while true
        attempt += 1
        iteration_state.scan_attempt = attempt
        
        log_scan_attempt_start(attempt, current_scan_count)
        
        # Run all 3 phases with current filtered_spectra
        for phase in 1:3
            log_phase_start(phase, attempt, length(filtered_spectra))
            
            converged = run_single_phase(
                phase, filtered_spectra, iteration_state,
                results, params, search_context, ms_file_idx, spectra
            )
            
            if converged
                log_success(attempt, phase, current_scan_count)
                return results
            end
        end
        
        # Check if we've reached or exceeded max scans
        if current_scan_count >= max_scans
            log_message("Reached maximum scan count ($max_scans), no convergence achieved")
            break
        end
        
        # Calculate next scan count
        next_scan_count = ceil(current_scan_count * scan_scale_factor)
        
        # Check if next iteration would exceed max
        if next_scan_count > max_scans
            # Do one final attempt with exactly max_scans
            additional_scans = max_scans - current_scan_count
            log_message("Adding final $additional_scans scans to reach maximum of $max_scans")
            
            append!(filtered_spectra, max_additional_scans = additional_scans)
            current_scan_count = max_scans
            
            # Run one more iteration with max scans
            attempt += 1
            iteration_state.scan_attempt = attempt
            
            for phase in 1:3
                log_phase_start(phase, attempt, length(filtered_spectra))
                
                converged = run_single_phase(
                    phase, filtered_spectra, iteration_state,
                    results, params, search_context, ms_file_idx, spectra
                )
                
                if converged
                    log_success(attempt, phase, current_scan_count)
                    return results
                end
            end
            
            # Final attempt failed, break out
            log_message("Final attempt with maximum scans failed to converge")
            break
            
        else
            # Normal scaling - append more scans
            additional_scans = next_scan_count - current_scan_count
            log_message("Appending $additional_scans scans (total will be $next_scan_count)")
            
            append!(filtered_spectra, max_additional_scans = additional_scans)
            current_scan_count = next_scan_count
        end
    end
    
    # Apply fallback strategy if no convergence
    if !converged
        log_warning("Failed to converge after $attempt attempts")
        apply_fallback_strategy!(results, search_context, params, ms_file_idx)
    end
    
    return results
end
```

## Key Differences from Previous Version

1. **Single FilteredMassSpecData Creation**
   - Created ONCE at the beginning with initial scan count
   - Reused across ALL phases and attempts
   - Only appended to when scaling up scan count

2. **Scan Scaling Logic**
   - After all 3 phases complete without convergence
   - Append additional scans to existing filtered_spectra
   - Special handling for final attempt at max_scans

3. **Memory Efficiency**
   - No repeated creation of filtered spectra
   - Incremental addition of scans
   - Maintains scan history throughout process

## Example Execution Flow

With settings:
- `initial_scan_count`: 500
- `max_parameter_tuning_scans`: 8000
- `scan_scale_factor`: 2.0
- `iterations_per_phase`: 3

```
Create FilteredMassSpecData with 500 scans

Attempt 1 (500 scans):
  Phase 1: Use existing 500 scans → No convergence
  Phase 2: Use existing 500 scans → No convergence  
  Phase 3: Use existing 500 scans → No convergence
  Append 500 more scans (total: 1000)

Attempt 2 (1000 scans):
  Phase 1: Use existing 1000 scans → No convergence
  Phase 2: Use existing 1000 scans → No convergence
  Phase 3: Use existing 1000 scans → No convergence
  Append 1000 more scans (total: 2000)

Attempt 3 (2000 scans):
  Phase 1: Use existing 2000 scans → No convergence
  Phase 2: Use existing 2000 scans → Converged!
  Success!
```

## Benefits of This Approach

1. **Efficiency**: No redundant spectra filtering
2. **Consistency**: Same spectra collection used across phases
3. **Progressive**: Builds on previous scan collections
4. **Clear Logic**: Explicit scan management
5. **Memory Friendly**: Single filtered spectra object

## Helper Functions

```pseudocode
function append!(filtered_spectra::FilteredMassSpecData; max_additional_scans::Int)
    # Append additional scans to existing filtered spectra
    # This adds to the existing collection, doesn't recreate
    current_count = length(filtered_spectra)
    target_count = current_count + max_additional_scans
    
    # Add scans from the underlying spectra object
    add_scans_from_source!(filtered_spectra, target_count)
    
    log_message("Filtered spectra now contains $(length(filtered_spectra)) scans")
end

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
```

## Summary

The key insight is that filtered spectra represents our working set of MS2 scans, and we should:
1. Create it once at the start
2. Reuse it across all phases within an attempt
3. Only append more scans when scaling up between attempts
4. Never recreate or reset it

This makes the algorithm more efficient and the logic clearer.