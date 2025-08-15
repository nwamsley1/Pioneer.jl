# Simplified Convergence Strategy for ParameterTuningSearch (V4)

## Overview
Clear convergence strategy with initial PSM collection outside the iteration loop, followed by expand→collect→adjust→collect cycles within the loop.

## Core Concepts

### Phase Structure (3 phases total)
1. **Phase 1**: Start with zero bias shift
2. **Phase 2**: Start with positive bias shift (+max_tolerance)
3. **Phase 3**: Start with negative bias shift (-max_tolerance)

### Within Each Phase
1. **Initial attempt** (before loop): Collect PSMs and check convergence at base tolerance
2. **Iteration loop** (up to `iterations_per_phase`): 
   - Expand tolerance
   - Collect PSMs to determine bias
   - Adjust bias based on fitting
   - Collect PSMs again with adjusted bias
   - Check convergence

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
    # ITERATION LOOP (expand → collect → adjust → collect cycles)
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
        
        # Step 2: COLLECT PSMs with expanded tolerance (to determine bias)
        psms_for_bias = collect_psms(
            filtered_spectra, spectra, search_context, params, ms_file_idx
        )
        log_message("Collected $(nrow(psms_for_bias)) PSMs with expanded tolerance")
        
        # Step 3: FIT MODEL to determine bias adjustment
        mass_err_model_for_bias, _, _ = fit_models_from_psms(
            psms_for_bias, spectra, search_context, params, ms_file_idx
        )
        
        # Step 4: ADJUST BIAS based on fitted model
        if mass_err_model_for_bias !== nothing
            new_bias = getMassOffset(mass_err_model_for_bias)
            current_model = getMassErrorModel(search_context, ms_file_idx)
            old_bias = getMassOffset(current_model)
            
            # Update the bias while keeping the expanded tolerance
            adjusted_model = MassErrorModel(
                new_bias,
                (getLeftTol(current_model), getRightTol(current_model))
            )
            setMassErrorModel!(search_context, ms_file_idx, adjusted_model)
            
            log_message("Adjusted bias from $old_bias to $new_bias ppm")
        else
            log_message("Could not fit model for bias adjustment, keeping current bias")
        end
        
        # Step 5: COLLECT PSMs again with adjusted bias
        psms_adjusted = collect_psms(
            filtered_spectra, spectra, search_context, params, ms_file_idx
        )
        log_message("Collected $(nrow(psms_adjusted)) PSMs with adjusted bias")
        
        # Step 6: FIT FINAL MODELS with bias-adjusted PSMs
        mass_err_model, rt_model, ppm_errs = fit_models_from_psms(
            psms_adjusted, spectra, search_context, params, ms_file_idx
        )
        
        # Step 7: CHECK CONVERGENCE
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
    ├── Collect PSMs with expanded tolerance
    ├── Fit model to determine bias
    ├── Adjust bias based on fitting
    ├── Collect PSMs with adjusted bias
    ├── Fit final models
    ├── Check convergence
    └── If converged → SUCCESS
```

## Detailed Example with Specific Settings

Settings:
- `init_mass_tol_ppm`: 20.0
- `mass_tolerance_scale_factor`: 1.5
- `iterations_per_phase`: 3

### Phase 1 Execution:
```
Initial Attempt:
  Bias = 0, Tolerance = 20 ppm
  → Collect PSMs
  → Fit models
  → Check convergence
  → Not converged

Iteration 1:
  → Expand: Tolerance = 30 ppm (20 × 1.5)
  → Collect: PSMs with 30 ppm tolerance
  → Fit: Determine bias should be +2.5 ppm
  → Adjust: Set bias to +2.5 ppm (keeping 30 ppm tolerance)
  → Collect: PSMs with 30 ppm tolerance and +2.5 ppm bias
  → Fit: Final models
  → Check convergence
  → Not converged

Iteration 2:
  → Expand: Tolerance = 45 ppm (20 × 1.5²)
  → Collect: PSMs with 45 ppm tolerance
  → Fit: Determine bias should be +1.8 ppm
  → Adjust: Set bias to +1.8 ppm (keeping 45 ppm tolerance)
  → Collect: PSMs with 45 ppm tolerance and +1.8 ppm bias
  → Fit: Final models
  → Check convergence
  → Not converged

Iteration 3:
  → Expand: Tolerance = 67.5 ppm (20 × 1.5³)
  → Collect: PSMs with 67.5 ppm tolerance
  → Fit: Determine bias should be +0.5 ppm
  → Adjust: Set bias to +0.5 ppm (keeping 67.5 ppm tolerance)
  → Collect: PSMs with 67.5 ppm tolerance and +0.5 ppm bias
  → Fit: Final models
  → Check convergence
  → Not converged
```

### Phase 2 Execution:
```
Initial Attempt:
  Bias = +67.5 ppm (phase shift), Tolerance = 20 ppm (reset)
  → Collect PSMs
  → Fit models
  → Check convergence
  → Not converged

Iteration 1:
  → Expand: Tolerance = 30 ppm
  → Collect: PSMs with 30 ppm tolerance and +67.5 ppm initial bias
  → Fit: Determine bias should actually be +65.2 ppm
  → Adjust: Set bias to +65.2 ppm (keeping 30 ppm tolerance)
  → Collect: PSMs with 30 ppm tolerance and +65.2 ppm bias
  → Fit: Final models
  → Check convergence
  → CONVERGED! ✓
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

## Key Points

1. **Two PSM collections per iteration**: First to determine bias, second with adjusted bias
2. **Bias determination needs data**: Can't adjust bias without first collecting PSMs at expanded tolerance
3. **Each iteration has 7 clear steps**: Expand → Collect → Fit → Adjust → Collect → Fit → Check
4. **Initial attempt is simpler**: Just collect → fit → check (no expand/adjust)
5. **Phase bias persists**: The initial phase bias shift remains as the starting point

## Why Two Collections?

The two-collection approach is necessary because:
1. **First collection**: With expanded tolerance but current bias - determines what the bias SHOULD be
2. **Second collection**: With expanded tolerance AND corrected bias - the actual PSMs we evaluate

This ensures we're always adjusting bias based on actual data at the current tolerance level.

## Benefits of This Structure

1. **Data-driven bias adjustment**: Bias is determined from actual PSMs, not guessed
2. **Clear causality**: Each step follows logically from the previous
3. **No hidden state**: All parameter changes are explicit
4. **Efficient convergence**: Bias adjustment is informed by current tolerance
5. **Debuggable**: Can track exactly why each decision was made