# Simplified Convergence Strategy for ParameterTuningSearch (V3)

## Overview
Clear convergence strategy with initial PSM collection outside the iteration loop, followed by expand/adjust cycles within the loop.

## Core Concepts

### Phase Structure (3 phases total)
1. **Phase 1**: Start with zero bias shift
2. **Phase 2**: Start with positive bias shift (+max_tolerance)
3. **Phase 3**: Start with negative bias shift (-max_tolerance)

### Within Each Phase
1. **Initial attempt** (before loop): Collect PSMs and check convergence at base tolerance
2. **Iteration loop** (up to `iterations_per_phase`): Expand tolerance → Adjust bias → Collect → Check convergence

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
    
    log_phase_start(phase, initial_bias_shift, settings.init_mass_tol_ppm)
    
    # ========================================
    # INITIAL ATTEMPT (before iteration loop)
    # ========================================
    
    # Collect PSMs at initial tolerance with phase bias
    psms_initial = collect_psms(
        filtered_spectra, spectra, search_context, params, ms_file_idx
    )
    
    # Fit models and check initial convergence
    mass_err_model, rt_model, ppm_errs = fit_models_from_psms(
        psms_initial, spectra, search_context, params, ms_file_idx
    )
    
    if check_convergence(mass_err_model, params)
        # Test tolerance expansion post-convergence
        expanded_result = test_tolerance_expansion!(
            search_context, params, ms_file_idx,
            psms_initial, mass_err_model, ppm_errs,
            settings.init_mass_tol_ppm, filtered_spectra, spectra
        )
        
        if expanded_result.improved
            store_results!(results, expanded_result)
            log_convergence_success("Initial attempt (expanded)", phase, 0)
        else
            store_results!(results, psms_initial, mass_err_model, ppm_errs)
            log_convergence_success("Initial attempt", phase, 0)
        end
        return true
    end
    
    log_message("Initial attempt did not converge, starting expand/adjust iterations")
    
    # ========================================
    # ITERATION LOOP (expand tolerance → adjust bias cycles)
    # ========================================
    
    for iter in 1:settings.iterations_per_phase
        iteration_state.current_iteration_in_phase = iter
        iteration_state.total_iterations += 1
        
        log_iteration_start(phase, iter)
        
        # Step 1: EXPAND TOLERANCE
        expand_mass_tolerance!(
            search_context, ms_file_idx, params,
            settings.mass_tolerance_scale_factor
        )
        current_tolerance = settings.init_mass_tol_ppm * 
                          (settings.mass_tolerance_scale_factor ^ iter)
        log_message("Expanded tolerance to $current_tolerance ppm")
        
        # Step 2: ADJUST BIAS (based on previous fitting)
        if mass_err_model !== nothing
            new_bias = getMassOffset(mass_err_model)
            adjust_mass_bias!(search_context, ms_file_idx, new_bias)
            log_message("Adjusted bias to $new_bias ppm")
        end
        
        # Step 3: COLLECT PSMs with expanded tolerance and adjusted bias
        psms_adjusted = collect_psms(
            filtered_spectra, spectra, search_context, params, ms_file_idx
        )
        
        # Step 4: FIT MODELS with new PSMs
        mass_err_model, rt_model, ppm_errs = fit_models_from_psms(
            psms_adjusted, spectra, search_context, params, ms_file_idx
        )
        
        # Step 5: CHECK CONVERGENCE
        if check_convergence(mass_err_model, params)
            # Test tolerance expansion post-convergence
            expanded_result = test_tolerance_expansion!(
                search_context, params, ms_file_idx,
                psms_adjusted, mass_err_model, ppm_errs,
                current_tolerance, filtered_spectra, spectra
            )
            
            if expanded_result.improved
                store_results!(results, expanded_result)
                log_convergence_success("Iteration $iter (expanded)", phase, iter)
            else
                store_results!(results, psms_adjusted, mass_err_model, ppm_errs)
                log_convergence_success("Iteration $iter", phase, iter)
            end
            return true
        end
        
        log_message("Iteration $iter did not converge")
    end
    
    log_phase_complete(phase, "No convergence after $(settings.iterations_per_phase) iterations")
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
        
        log_separator("=", 60)
        log_message("SCAN ATTEMPT $attempt: $current_scan_count scans")
        log_separator("=", 60)
        
        # Run all 3 phases with current filtered_spectra
        for phase in 1:3
            log_separator("-", 40)
            log_message("PHASE $phase (Attempt $attempt)")
            log_separator("-", 40)
            
            converged = run_single_phase(
                phase, filtered_spectra, iteration_state,
                results, params, search_context, ms_file_idx, spectra
            )
            
            if converged
                log_success("CONVERGED in Phase $phase, Attempt $attempt with $current_scan_count scans")
                return results
            end
        end
        
        # Check if we've reached or exceeded max scans
        if current_scan_count >= max_scans
            log_warning("Reached maximum scan count ($max_scans), no convergence achieved")
            break
        end
        
        # Calculate next scan count
        next_scan_count = Int(ceil(current_scan_count * scan_scale_factor))
        
        # Check if next iteration would exceed max
        if next_scan_count > max_scans
            if current_scan_count < max_scans
                # Do one final attempt with exactly max_scans
                additional_scans = max_scans - current_scan_count
                log_message("Adding final $additional_scans scans to reach maximum of $max_scans")
                
                append!(filtered_spectra, max_additional_scans = additional_scans)
                current_scan_count = max_scans
                # Loop will continue for one more attempt
            else
                # Already at max, break
                break
            end
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

## Phase Execution Flow Diagram

```
Phase N:
│
├── Apply phase bias shift (0, +max_tol, -max_tol)
├── Reset tolerance to init_mass_tol_ppm
│
├── [INITIAL ATTEMPT]
│   ├── Collect PSMs at base tolerance
│   ├── Fit models
│   ├── Check convergence
│   └── If converged → SUCCESS
│
└── [ITERATION LOOP] (i = 1 to iterations_per_phase)
    ├── Expand tolerance (× scale_factor)
    ├── Adjust bias (from previous model)
    ├── Collect PSMs
    ├── Fit models
    ├── Check convergence
    └── If converged → SUCCESS
```

## Example Execution with Specific Settings

Settings:
- `init_mass_tol_ppm`: 20.0
- `mass_tolerance_scale_factor`: 1.5
- `iterations_per_phase`: 3

### Phase 1 Execution:
```
Initial: Bias = 0, Tolerance = 20 ppm
  → Collect PSMs → Fit → Check convergence
  → Not converged

Iteration 1: Tolerance = 30 ppm (20 × 1.5)
  → Expand tolerance to 30
  → Adjust bias based on initial fit
  → Collect PSMs → Fit → Check convergence
  → Not converged

Iteration 2: Tolerance = 45 ppm (20 × 1.5²)
  → Expand tolerance to 45
  → Adjust bias based on iteration 1 fit
  → Collect PSMs → Fit → Check convergence
  → Not converged

Iteration 3: Tolerance = 67.5 ppm (20 × 1.5³)
  → Expand tolerance to 67.5
  → Adjust bias based on iteration 2 fit
  → Collect PSMs → Fit → Check convergence
  → Not converged
```

### Phase 2 Execution:
```
Initial: Bias = +67.5, Tolerance = 20 ppm (reset)
  → Collect PSMs → Fit → Check convergence
  → Not converged

Iteration 1: Tolerance = 30 ppm
  → Expand tolerance to 30
  → Adjust bias based on initial fit
  → Collect PSMs → Fit → Check convergence
  → CONVERGED!
```

## Key Points

1. **Initial attempt is special**: Each phase starts with a collection at base tolerance before any expansion
2. **Iteration loop is expand→adjust→collect**: Always expand first, then adjust bias, then collect
3. **Bias adjustment uses previous model**: The bias from the previous fit guides the next collection
4. **Tolerance resets per phase**: Each phase starts fresh at init_mass_tol_ppm
5. **Phase bias persists**: The initial phase bias shift remains throughout the phase

## Helper Functions

```pseudocode
function get_phase_bias_shift(phase, settings)
    # Calculate the maximum tolerance that would be reached
    max_tolerance = settings.init_mass_tol_ppm * 
                   (settings.mass_tolerance_scale_factor ^ settings.iterations_per_phase)
    
    if phase == 1
        return 0.0  # No bias shift
    elseif phase == 2
        return +max_tolerance  # Positive shift
    elseif phase == 3
        return -max_tolerance  # Negative shift
    end
end

function check_convergence(mass_err_model, params)
    if mass_err_model === nothing
        return false
    end
    
    mass_offset = getMassOffset(mass_err_model)
    left_tol = getLeftTol(mass_err_model)
    right_tol = getRightTol(mass_err_model)
    init_tol = params.iteration_settings.init_mass_tol_ppm
    
    # Convergence criteria:
    # 1. Mass offset should be less than 1/4 of initial tolerance
    if abs(mass_offset) > (init_tol / 4)
        return false
    end
    
    # 2. Average fitted tolerance should be less than initial tolerance
    avg_fitted_tol = (left_tol + right_tol) / 2
    if avg_fitted_tol > init_tol
        return false
    end
    
    return true
end

function test_tolerance_expansion!(...)
    # After convergence, test if expanding tolerance by 50% finds more PSMs
    # Returns expanded results if improvement, otherwise returns original
end
```

## Benefits of This Structure

1. **Clear separation**: Initial attempt vs iteration cycles
2. **Predictable flow**: Always expand→adjust→collect in iterations
3. **Efficient**: No redundant operations
4. **Debuggable**: Each step has clear purpose and logging
5. **Maintainable**: Simple logic flow without hidden complexity