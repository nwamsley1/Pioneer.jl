# Plan to Combine Strategies 2 and 3

## Current Structure
Currently we have three separate strategies:
1. **Strategy 1**: Collect PSMs with current parameters → CHECK CONVERGENCE
2. **Strategy 2**: Expand mass tolerance → NO CHECK (always followed by Strategy 3)
3. **Strategy 3**: Adjust mass bias → CHECK CONVERGENCE

## Problem with Current Structure
- Strategy 2 never checks convergence
- Strategy 2 is ALWAYS followed by Strategy 3
- This creates unnecessary complexity and redundant code

## Proposed New Structure
Two strategies instead of three:
1. **Strategy 1**: Collect PSMs with current parameters → CHECK CONVERGENCE
2. **Strategy 2**: Expand tolerance + Adjust bias → CHECK CONVERGENCE

## Implementation Changes

### 1. Modify `execute_strategy` Function

**Current Implementation:**
```julia
function execute_strategy(strategy_num::Int, filtered_spectra, spectra, search_context, 
                         params, ms_file_idx, results)
    
    if strategy_num == 1
        @info "Strategy 1: Collecting PSMs with current parameters"
        psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context, 
                                       params, ms_file_idx, "with current parameters")
    
    elseif strategy_num == 2
        @info "Strategy 2: Expanding mass tolerance"
        expand_mass_tolerance!(search_context, ms_file_idx, params)
        psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                       params, ms_file_idx, "with expanded tolerance")
    
    elseif strategy_num == 3
        @info "Strategy 3: Adjusting mass bias"
        # Need mass error model from previous attempt
        psms_temp, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                           params, ms_file_idx, "for bias estimation")
        mass_err_temp, _, _ = fit_models_from_psms(psms_temp, spectra, search_context, 
                                                   params, ms_file_idx)
        
        if mass_err_temp === nothing
            @info "Skipping Strategy 3 - no valid mass error model"
            return nothing, nothing, nothing
        end
        
        adjust_mass_bias!(search_context, ms_file_idx, mass_err_temp, params)
        psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                       params, ms_file_idx, "with adjusted bias")
    end
    
    # Fit models from collected PSMs
    mass_err_model, ppm_errs, psm_count = fit_models_from_psms(psms, spectra, 
                                                                search_context, params, ms_file_idx)
    
    return psms, mass_err_model, ppm_errs
end
```

**New Implementation:**
```julia
function execute_strategy(strategy_num::Int, filtered_spectra, spectra, search_context, 
                         params, ms_file_idx, results)
    
    if strategy_num == 1
        @info "Strategy 1: Collecting PSMs with current parameters"
        psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context, 
                                       params, ms_file_idx, "with current parameters")
    
    elseif strategy_num == 2
        @info "Strategy 2: Expanding mass tolerance and adjusting bias"
        
        # First expand mass tolerance
        expand_mass_tolerance!(search_context, ms_file_idx, params)
        
        # Collect PSMs with expanded tolerance to estimate bias
        psms_temp, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                           params, ms_file_idx, "with expanded tolerance for bias estimation")
        
        # Fit mass error model to get bias estimate
        mass_err_temp, _, _ = fit_models_from_psms(psms_temp, spectra, search_context, 
                                                   params, ms_file_idx)
        
        if mass_err_temp === nothing
            @info "No valid mass error model for bias adjustment, using expanded tolerance only"
            psms = psms_temp  # Use the PSMs we already collected
        else
            # Adjust bias based on the fitted model
            adjust_mass_bias!(search_context, ms_file_idx, mass_err_temp, params)
            
            # Collect final PSMs with expanded tolerance AND adjusted bias
            psms, _ = collect_and_log_psms(filtered_spectra, spectra, search_context,
                                           params, ms_file_idx, "with expanded tolerance and adjusted bias")
        end
    
    else
        error("Invalid strategy number: $strategy_num. Only strategies 1 and 2 are supported.")
    end
    
    # Fit models from collected PSMs
    mass_err_model, ppm_errs, psm_count = fit_models_from_psms(psms, spectra, 
                                                                search_context, params, ms_file_idx)
    
    return psms, mass_err_model, ppm_errs
end
```

### 2. Modify Main Loop in `process_file!`

**Current Implementation:**
```julia
# Strategy 1: Try with current parameters
psms, mass_err_model, ppm_errs = execute_strategy(
    1, filtered_spectra, spectra, 
    search_context, params, ms_file_idx, results
)

# Check convergence after Strategy 1
if check_and_store_convergence!(
    results, search_context, params, ms_file_idx,
    psms, mass_err_model, ppm_errs, "Strategy 1"
)
    converged = true
    final_psm_count = psms !== nothing ? size(psms, 1) : 0
else
    # Strategy 2: Expand tolerance (no convergence check)
    psms, mass_err_model, ppm_errs = execute_strategy(
        2, filtered_spectra, spectra,
        search_context, params, ms_file_idx, results
    )
    
    # Strategy 3: Adjust bias (always follows Strategy 2)
    psms, mass_err_model, ppm_errs = execute_strategy(
        3, filtered_spectra, spectra,
        search_context, params, ms_file_idx, results
    )
    
    # Check convergence after Strategy 3
    if check_and_store_convergence!(
        results, search_context, params, ms_file_idx,
        psms, mass_err_model, ppm_errs, "Strategy 3"
    )
        converged = true
        final_psm_count = psms !== nothing ? size(psms, 1) : 0
    end
end
```

**New Implementation:**
```julia
# Strategy 1: Try with current parameters
psms, mass_err_model, ppm_errs = execute_strategy(
    1, filtered_spectra, spectra, 
    search_context, params, ms_file_idx, results
)

# Check convergence after Strategy 1
if check_and_store_convergence!(
    results, search_context, params, ms_file_idx,
    psms, mass_err_model, ppm_errs, "Strategy 1"
)
    converged = true
    final_psm_count = psms !== nothing ? size(psms, 1) : 0
else
    # Strategy 2: Expand tolerance AND adjust bias (combined)
    psms, mass_err_model, ppm_errs = execute_strategy(
        2, filtered_spectra, spectra,
        search_context, params, ms_file_idx, results
    )
    
    # Check convergence after Strategy 2
    if check_and_store_convergence!(
        results, search_context, params, ms_file_idx,
        psms, mass_err_model, ppm_errs, "Strategy 2 (expanded tolerance + bias adjustment)"
    )
        converged = true
        final_psm_count = psms !== nothing ? size(psms, 1) : 0
    end
end
```

## Benefits of This Change

1. **Simpler Logic**: Only two strategies instead of three
2. **Clearer Intent**: Strategy 2 is clearly "expand and adjust"
3. **Less Code Duplication**: No need to handle strategy 3 separately
4. **More Efficient**: One less call to `execute_strategy`
5. **Maintains Same Behavior**: Still expands tolerance then adjusts bias

## Testing Considerations

1. Verify that the combined strategy produces the same results
2. Check that convergence checking still works correctly
3. Ensure logging is clear about what's happening
4. Test edge cases (no PSMs found, no valid mass error model, etc.)

## Implementation Steps

1. Modify `execute_strategy` to handle only strategies 1 and 2
2. Update the main loop in `process_file!` to call only strategies 1 and 2
3. Update logging messages to be clear about the combined operation
4. Test the changes with sample data
5. Update any documentation that references the three-strategy approach