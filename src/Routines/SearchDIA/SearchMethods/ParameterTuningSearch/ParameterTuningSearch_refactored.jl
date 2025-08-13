# New refactored process_file! function that uses the extracted helper functions

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