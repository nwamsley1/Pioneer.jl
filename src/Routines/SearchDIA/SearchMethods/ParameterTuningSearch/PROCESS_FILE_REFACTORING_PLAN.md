# Process File Refactoring Plan

## Current Issues

### 1. Code Duplication
- **PSM Collection Pattern**: Repeated 3 times (lines 317, 381, 424)
  ```julia
  psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
  final_psm_count = size(psms, 1)
  @info "Collected $final_psm_count PSMs..."
  ```

- **Fragment Matching Pattern**: Repeated 3 times (lines 323, 387, 430)
  ```julia
  fragments = get_matched_fragments(spectra, psms, search_context, params, ms_file_idx)
  if length(fragments) > 0
      mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
  ```

- **Convergence Check & Model Storage**: Repeated 2 times (lines 329-348, 437-456)
  ```julia
  if check_convergence(psms, mass_err_model, ...)
      @info "Converged after..."
      setMassErrorModel!(search_context, ms_file_idx, mass_err_model)
      # Store ppm_errs for plotting (5 lines of duplicate code)
      # Store RT model (2 lines of duplicate code)
      converged = true
      break
  ```

- **Mass Error Model Setting**: Repeated pattern for setting/getting models

### 2. Unnecessary Variables
- `ppm_errs` - initialized to `nothing` and barely used
- `initial_mass_offset` - set but never meaningfully used (always 0)
- Multiple temporary tolerance variables that could be simplified

### 3. Control Flow Issues
- Deeply nested if-else blocks (up to 4 levels)
- Unclear separation between the 3 strategies
- Mixed concerns (model fitting, logging, state management)
- Inconsistent handling of edge cases (no PSMs, no fragments)

### 4. Poor Separation of Concerns
- Main loop mixes:
  - Scan management
  - PSM collection
  - Model fitting
  - Convergence checking
  - Result storage
  - Logging

## Proposed Refactoring

### 1. Extract Helper Functions

```julia
# Collect PSMs and log results
function collect_and_log_psms(filtered_spectra, spectra, search_context, params, ms_file_idx, context_msg::String)
    psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    psm_count = size(psms, 1)
    @info "Collected $psm_count PSMs $context_msg"
    return psms, psm_count
end

# Fit models from PSMs
function fit_models_from_psms(psms, spectra, search_context, params, ms_file_idx)
    psm_count = size(psms, 1)
    
    if psm_count == 0
        return nothing, nothing, 0
    end
    
    fragments = get_matched_fragments(spectra, psms, search_context, params, ms_file_idx)
    
    if length(fragments) == 0
        return nothing, nothing, psm_count
    end
    
    mass_err_model, ppm_errs = fit_mass_err_model(params, fragments)
    return mass_err_model, ppm_errs, psm_count
end

# Check convergence and store results if converged
function check_and_store_convergence!(results, search_context, params, ms_file_idx, 
                                      psms, mass_err_model, ppm_errs, strategy_name::String)
    if mass_err_model === nothing || psms === nothing
        return false
    end
    
    current_model = getMassErrorModel(search_context, ms_file_idx)
    
    if !check_convergence(psms, mass_err_model, current_model, ppm_errs, getMinPsms(params))
        return false
    end
    
    @info "Converged after $strategy_name"
    
    # Store mass error model
    setMassErrorModel!(search_context, ms_file_idx, mass_err_model)
    
    # Store ppm errors for plotting
    if ppm_errs !== nothing && length(ppm_errs) > 0
        resize!(results.ppm_errs, 0)
        append!(results.ppm_errs, ppm_errs)
        @info "Stored $(length(results.ppm_errs)) ppm errors for plotting"
    end
    
    # Store RT model
    rt_model_data = fit_irt_model(params, psms)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
    @info "Stored $(length(results.rt)) RT points and $(length(results.irt)) iRT points"
    
    return true
end

# Add scans using doubling strategy
function add_scans_for_iteration!(filtered_spectra, params, iteration::Int)
    if iteration == 0
        additional_scans = getInitialScanCount(params)
    else
        current_scans = length(filtered_spectra)
        target_scans = min(current_scans * 2, getMaxParameterTuningScans(params))
        additional_scans = target_scans - current_scans
        
        if additional_scans <= 0
            @info "Reached maximum scan count of $(getMaxParameterTuningScans(params))"
            return false  # Signal to stop iterations
        end
    end
    
    @info "Iteration $(iteration+1): Adding $additional_scans scans (total: $(length(filtered_spectra) + additional_scans))"
    append!(filtered_spectra; max_additional_scans = additional_scans)
    return true
end

# Execute a single convergence strategy
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

# Helper to expand mass tolerance
function expand_mass_tolerance!(search_context, ms_file_idx, params)
    current_model = getMassErrorModel(search_context, ms_file_idx)
    current_left = getLeftTol(current_model)
    current_right = getRightTol(current_model)
    
    new_left = current_left * 2.0f0
    new_right = current_right * 2.0f0
    
    new_model = create_capped_mass_model(
        getMassOffset(current_model),
        new_left,
        new_right,
        getMaxTolerancePpm(params)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
    
    actual_left = getLeftTol(new_model)
    actual_right = getRightTol(new_model)
    @info "Expanded tolerance from ($(round(current_left, digits=1)), $(round(current_right, digits=1))) " *
          "to ($(round(actual_left, digits=1)), $(round(actual_right, digits=1))) ppm"
end

# Helper to adjust mass bias
function adjust_mass_bias!(search_context, ms_file_idx, new_mass_err_model, params)
    current_model = getMassErrorModel(search_context, ms_file_idx)
    new_bias = getMassOffset(new_mass_err_model)
    old_bias = getMassOffset(current_model)
    
    @info "Adjusting mass bias from $(round(old_bias, digits=2)) to $(round(new_bias, digits=2)) ppm"
    
    new_model = create_capped_mass_model(
        new_bias,
        getLeftTol(current_model),
        getRightTol(current_model),
        getMaxTolerancePpm(params)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
end
```

### 2. Simplified Main Loop Structure

```julia
function process_file!(results, params, search_context, ms_file_idx, spectra)
    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    n_attempts = 0
    warnings = String[]
    final_psm_count = 0
    
    try
        @info "Processing file: $parsed_fname (index: $ms_file_idx)"
        
        # Initialize models
        initialize_models!(search_context, ms_file_idx, params)
        
        # Initialize filtered spectra
        filtered_spectra = initialize_filtered_spectra(spectra, params)
        
        # Main convergence loop
        for iteration in 0:4  # Max 5 iterations
            n_attempts = iteration + 1
            
            # Add scans for this iteration
            if !add_scans_for_iteration!(filtered_spectra, params, iteration)
                break  # Reached max scans
            end
            
            # Try three strategies in sequence
            for strategy in 1:3
                psms, mass_err_model, ppm_errs = execute_strategy(
                    strategy, filtered_spectra, spectra, 
                    search_context, params, ms_file_idx, results
                )
                
                # Check convergence
                if check_and_store_convergence!(
                    results, search_context, params, ms_file_idx,
                    psms, mass_err_model, ppm_errs, "Strategy $strategy"
                )
                    converged = true
                    final_psm_count = size(psms, 1)
                    break
                end
            end
            
            if converged
                break
            end
            
            @info "Iteration $n_attempts complete. Moving to next iteration..."
        end
        
        # Store results
        store_final_results!(results, search_context, params, ms_file_idx, 
                           converged, n_attempts, final_psm_count, warnings)
        
        # Handle non-convergence
        if !converged
            handle_non_convergence!(results, search_context, ms_file_idx, 
                                   n_attempts, warnings)
        end
        
    catch e
        handle_error!(results, search_context, ms_file_idx, e, parsed_fname)
    end
end
```

### 3. Benefits of Refactoring

#### Reduced Code Duplication
- Single function for PSM collection and logging
- Single function for model fitting
- Single function for convergence checking and storage
- Eliminated ~100 lines of duplicate code

#### Improved Clarity
- Clear separation of strategies
- Linear control flow
- Each function has single responsibility
- Consistent error handling

#### Better Maintainability
- Changes to PSM collection only need one update
- Strategies are isolated and testable
- Easier to add new strategies
- Clearer parameter flow

#### Reduced Complexity
- Maximum nesting reduced from 4 to 2 levels
- Each function is <30 lines
- Clear naming conventions
- Explicit return values

### 4. Implementation Steps

1. **Phase 1**: Extract helper functions
   - `collect_and_log_psms`
   - `fit_models_from_psms`
   - `check_and_store_convergence!`

2. **Phase 2**: Extract strategy functions
   - `add_scans_for_iteration!`
   - `execute_strategy`
   - Strategy-specific helpers

3. **Phase 3**: Refactor main loop
   - Simplify control flow
   - Remove unnecessary variables
   - Clean up initialization

4. **Phase 4**: Testing and validation
   - Ensure identical behavior
   - Test edge cases
   - Performance verification

### 5. Testing Checklist

- [ ] Verify convergence with good data
- [ ] Test handling of no PSMs
- [ ] Test handling of no fragments
- [ ] Test max scan limit
- [ ] Test all three strategies
- [ ] Verify fallback handling
- [ ] Check error handling
- [ ] Validate logging output

## Summary

This refactoring will:
- Reduce function size from ~300 lines to ~100 lines
- Eliminate ~100 lines of duplicate code
- Improve readability and maintainability
- Make the control flow clear and testable
- Preserve all existing functionality