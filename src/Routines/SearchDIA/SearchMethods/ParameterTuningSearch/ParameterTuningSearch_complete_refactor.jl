# Complete refactored version with all helper functions properly extracted

#==========================================================
# Helper Functions (extracted from old process_file!)
#==========================================================

"""
    collect_and_log_psms(filtered_spectra, spectra, search_context, params, ms_file_idx, context_msg::String)

Collect PSMs and log the count with context message.
"""
function collect_and_log_psms(filtered_spectra, spectra, search_context, params, ms_file_idx, context_msg::String)
    psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    psm_count = size(psms, 1)
    @info "Collected $psm_count PSMs $context_msg"
    return psms, psm_count
end

"""
    fit_models_from_psms(psms, spectra, search_context, params, ms_file_idx)

Fit mass error and RT models from PSMs.
"""
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

"""
    check_and_store_convergence!(results, search_context, params, ms_file_idx, psms, mass_err_model, ppm_errs, strategy_name)

Check convergence criteria and store results if converged.
"""
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

"""
    add_scans_for_iteration!(filtered_spectra, params, iteration::Int)

Add scans using doubling strategy for current iteration.
Returns false if maximum scan count reached.
"""
function add_scans_for_iteration!(filtered_spectra, params, iteration::Int)
    if iteration == 0
        additional_scans = getInitialScanCount(params)
    else
        current_scans = length(filtered_spectra)
        target_scans = min(current_scans * 2, getMaxParameterTuningScans(params))
        additional_scans = target_scans - current_scans
        
        if additional_scans <= 0
            @info "Reached maximum scan count of $(getMaxParameterTuningScans(params))"
            return false
        end
    end
    
    @info "Iteration $(iteration+1): Adding $additional_scans scans (total: $(length(filtered_spectra) + additional_scans))"
    append!(filtered_spectra; max_additional_scans = additional_scans)
    return true
end

"""
    expand_mass_tolerance!(search_context, ms_file_idx, params)

Double the current mass tolerance up to maximum.
"""
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

"""
    adjust_mass_bias!(search_context, ms_file_idx, new_mass_err_model, params)

Adjust mass bias while keeping current tolerances.
"""
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

"""
    execute_strategy(strategy_num::Int, filtered_spectra, spectra, search_context, params, ms_file_idx, results)

Execute one of three convergence strategies.
"""
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

"""
    initialize_models!(search_context, ms_file_idx, params)

Initialize mass error and quad transmission models for file.
"""
function initialize_models!(search_context, ms_file_idx, params)
    # Use fixed initial parameters from JSON configuration
    initial_bias = 0.0f0  # Always start with zero bias
    initial_tolerance = getFragTolPpm(params)
    max_tolerance = getMaxTolerancePpm(params)
    
    # Set initial mass error model
    setMassErrorModel!(search_context, ms_file_idx, create_capped_mass_model(
        initial_bias,
        initial_tolerance,
        initial_tolerance,
        max_tolerance
    ))
    
    # Set initial quad transmission model
    setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))
    
    @info "Initial mass error model for file $ms_file_idx: $(getMassErrorModel(search_context, ms_file_idx))"
end

"""
    initialize_filtered_spectra(spectra, params)

Create FilteredMassSpecData with initial settings.
"""
function initialize_filtered_spectra(spectra, params)
    topn_peaks = something(getTopNPeaks(params), 200)
    return FilteredMassSpecData(
        spectra,
        max_scans = 0,  # Start with zero scans
        topn = topn_peaks,
        target_ms_order = UInt8(2)  # Only MS2 for presearch
    )
end

"""
    store_final_results!(results, search_context, params, ms_file_idx, converged, n_attempts, final_psm_count, warnings)

Store tuning results and diagnostic information.
"""
function store_final_results!(results, search_context, params, ms_file_idx, 
                              converged, n_attempts, final_psm_count, warnings)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Log final status
    @info "Final data status for file $ms_file_idx:"
    @info "  - RT points: $(length(results.rt))"
    @info "  - iRT points: $(length(results.irt))" 
    @info "  - PPM errors: $(length(results.ppm_errs))"
    @info "  - Converged: $converged"
    
    # Record tuning results
    tuning_result = TuningResults(
        getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        (getLeftTol(getMassErrorModel(search_context, ms_file_idx)), 
         getRightTol(getMassErrorModel(search_context, ms_file_idx))),
        converged,
        final_psm_count,
        n_attempts,
        warnings
    )
    store_tuning_results!(results.parameter_history, ms_file_idx, tuning_result)
    
    # Record diagnostic status
    status = ParameterTuningStatus(
        ms_file_idx,
        parsed_fname,
        converged,
        !converged,  # used_fallback
        !converged ? "Failed to converge after $n_attempts attempts" : "",
        n_attempts,
        final_psm_count,
        getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        (getLeftTol(getMassErrorModel(search_context, ms_file_idx)), 
         getRightTol(getMassErrorModel(search_context, ms_file_idx))),
        warnings
    )
    record_tuning_status!(results.diagnostics, status)
    
    # Store final mass error model in results
    results.mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
end

"""
    handle_non_convergence!(results, search_context, ms_file_idx, n_attempts, warnings)

Handle case when parameter tuning fails to converge.
"""
function handle_non_convergence!(results, search_context, ms_file_idx, n_attempts, warnings)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    @warn "Failed to converge mass error model after $n_attempts attempts for file $ms_file_idx."
    
    # Try to borrow parameters from neighboring files or use defaults
    fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
        search_context, ms_file_idx, warnings
    )
    
    # Apply fallback/borrowed parameters
    results.mass_err_model[] = fallback_mass_err
    results.rt_to_irt_model[] = fallback_rt_model
    
    # Store models in SearchContext for downstream methods
    setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
    setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)
    
    # Set default IRT error for fallback/borrowed parameters
    if borrowed_from !== nothing && haskey(getIrtErrors(search_context), borrowed_from)
        getIrtErrors(search_context)[ms_file_idx] = getIrtErrors(search_context)[borrowed_from]
    else
        getIrtErrors(search_context)[ms_file_idx] = 2.0f0  # Conservative default
    end
    
    if borrowed_from !== nothing
        borrowed_fname = getParsedFileName(search_context, borrowed_from)
        @warn "PARAMETER_TUNING_BORROWED: File $parsed_fname borrowed parameters from $borrowed_fname"
    else
        @warn "PARAMETER_TUNING_FALLBACK: File $parsed_fname used fallback values (±50 ppm, identity RT)"
    end
end

"""
    handle_error!(results, search_context, ms_file_idx, e, parsed_fname)

Handle errors during parameter tuning.
"""
function handle_error!(results, search_context, ms_file_idx, e, parsed_fname)
    @warn "Parameter tuning failed for file $parsed_fname with error." exception=(e, catch_backtrace())
    
    error_warnings = String["Error: $(string(e))"]
    fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
        search_context, ms_file_idx, error_warnings
    )
    
    # Apply fallback parameters
    results.mass_err_model[] = fallback_mass_err
    results.rt_to_irt_model[] = fallback_rt_model
    setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
    setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)
    
    # Set default IRT error
    if borrowed_from !== nothing && haskey(getIrtErrors(search_context), borrowed_from)
        getIrtErrors(search_context)[ms_file_idx] = getIrtErrors(search_context)[borrowed_from]
    else
        getIrtErrors(search_context)[ms_file_idx] = 2.0f0
    end
    
    # Record error status
    mass_tols = (getLeftTol(results.mass_err_model[]), getRightTol(results.mass_err_model[]))
    status = ParameterTuningStatus(
        ms_file_idx,
        parsed_fname,
        false,  # not converged
        true,   # used fallback
        "Exception: $(typeof(e))",
        0,      # n_iterations
        0,      # psm_count
        getMassOffset(results.mass_err_model[]),
        mass_tols,
        error_warnings
    )
    record_tuning_status!(results.diagnostics, status)
    
    if borrowed_from !== nothing
        borrowed_fname = getParsedFileName(search_context, borrowed_from)
        @warn "PARAMETER_TUNING_ERROR: File $parsed_fname had error $(typeof(e)). Borrowed from $borrowed_fname."
    else
        @warn "PARAMETER_TUNING_ERROR: File $parsed_fname had error $(typeof(e)). Used fallback values."
    end
    
    # Store final model in results
    results.mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
end

#==========================================================
# Main Refactored process_file! Function
#==========================================================

"""
Process a single MS file to determine optimal mass error and RT parameters.
"""
function process_file!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    
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
                    final_psm_count = psms !== nothing ? size(psms, 1) : 0
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
    
    return results
end