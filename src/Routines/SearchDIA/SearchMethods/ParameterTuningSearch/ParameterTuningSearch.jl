# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    ParameterTuningSearch

Search method for tuning mass error and retention time parameters.

This search:
1. Estimates mass error distributions for fragment matching
2. Fits retention time alignment between library and empirical data
3. Generates QC plots for parameter estimation quality
4. Stores models in SearchContext for use by other search methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 2),
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "frag_err_quantile" => 0.01,
        "min_samples" => 1000,
        "max_qval" => 0.01,
        "max_presearch_iters" => 10,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "sample_rate" => 0.1,
        "abreviate_precursor_calc" => false
    ),
    :irt_mapping_params => Dict(
        "n_sigma_tol" => 3
    )
)

# Execute search
results = execute_search(ParameterTuningSearch(), search_context, params)
```
"""
# Type definitions moved to types.jl

#==========================================================
Results Access Methods
==========================================================#
getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getRtToIrtModel(ptsr::ParameterTuningSearchResults) = ptsr.rt_to_irt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path
getDiagnostics(ptsr::ParameterTuningSearchResults) = ptsr.diagnostics
getParameterHistory(ptsr::ParameterTuningSearchResults) = ptsr.parameter_history

function set_rt_to_irt_model!(
    ptsr::ParameterTuningSearchResults, 
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64,
    model::Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32}
) where {P<:ParameterTuningSearchParameters}
    
    ptsr.rt_to_irt_model[] = model[1]
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    
    #parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
end

"""
    get_fallback_parameters(search_context, ms_file_idx, warnings)

Attempt to borrow parameters from successfully tuned neighboring files.
Falls back to conservative defaults if no successfully tuned files are available.

Returns: (mass_error_model, rt_model, borrowed_from_idx)
"""
function get_fallback_parameters(
    search_context::SearchContext, 
    ms_file_idx::Int64,
    warnings::Vector{String}
)
    # Try to find a successfully tuned file's parameters
    borrowed_from = nothing
    fallback_mass_err = nothing
    fallback_rt_model = nothing
    
    ms_data = getMSData(search_context)
    n_files = length(ms_data.file_paths)
    
    # Check previous files first (prefer neighboring files)
    for file_idx in (ms_file_idx-1):-1:1
        if haskey(search_context.rt_irt_map, file_idx) && 
           haskey(search_context.mass_error_model, file_idx)
            borrowed_from = file_idx
            fallback_mass_err = search_context.mass_error_model[file_idx]
            fallback_rt_model = search_context.rt_irt_map[file_idx]
            break
        end
    end
    
    # If no previous file, check subsequent files
    if borrowed_from === nothing
        for file_idx in (ms_file_idx+1):n_files
            if haskey(search_context.rt_irt_map, file_idx) && 
               haskey(search_context.mass_error_model, file_idx)
                borrowed_from = file_idx
                fallback_mass_err = search_context.mass_error_model[file_idx]
                fallback_rt_model = search_context.rt_irt_map[file_idx]
                break
            end
        end
    end
    
    # Return borrowed parameters or conservative defaults
    if borrowed_from !== nothing
        borrowed_fname = getParsedFileName(search_context, borrowed_from)
        @info "Borrowing parameters from file $borrowed_fname (index $borrowed_from) for file $ms_file_idx"
        push!(warnings, "BORROWED: Using parameters from file $borrowed_fname")
        return fallback_mass_err, fallback_rt_model, borrowed_from
    else
        # Use conservative defaults if no other files available
        @info "No successfully tuned files available. Using conservative defaults for file $ms_file_idx"
        push!(warnings, "FALLBACK: Used conservative default parameters (no borrowing available)")
        return MassErrorModel(0.0f0, (50.0f0, 50.0f0)), IdentityModel(), nothing
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::ParameterTuningSearch, params::Any) = ParameterTuningSearchParameters(params)

function init_search_results(::ParameterTuningSearchParameters, search_context::SearchContext)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)
    rt_alingment_plots = joinpath(qc_dir, "rt_alignment_plots")
    !isdir(rt_alingment_plots) && mkdir(rt_alingment_plots)
    mass_error_plots = joinpath(qc_dir, "mass_error_plots")
    !isdir(mass_error_plots) && mkdir(mass_error_plots)
    ms1_mass_error_plots = joinpath(qc_dir, "ms1_mass_error_plots")
    !isdir(ms1_mass_error_plots ) && mkdir(ms1_mass_error_plots )
    return ParameterTuningSearchResults(
        Base.Ref{MassErrorModel}(),
        Ref{RtConversionModel}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        qc_dir,
        ParameterTuningDiagnostics(),
        ParameterHistory()
    )
end

#==========================================================
Helper Functions for Refactored process_file!
==========================================================#

"""
Collect PSMs from filtered spectra using library search.
"""
function collect_psms(
    filtered_spectra::FilteredMassSpecData,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:ParameterTuningSearchParameters}
    
    # Perform library search on filtered data
    psms = library_search(filtered_spectra, search_context, params, ms_file_idx)
    
    if !iszero(size(psms, 1))
        # CRITICAL: Map filtered scan indices back to original
        # library_search returns scan_idx relative to filtered_spectra
        psms[!, :filtered_scan_idx] = psms[!, :scan_idx]
        psms[!, :scan_idx] = [
            getOriginalScanIndex(filtered_spectra, idx) 
            for idx in psms[!, :filtered_scan_idx]
        ]
        
        # Add metadata columns using ORIGINAL spectra
        precursors = getPrecursors(getSpecLib(search_context))
        add_tuning_search_columns!(
            psms,
            spectra,
            getIsDecoy(precursors),
            getIrt(precursors),
            getCharge(precursors),
            getRetentionTimes(spectra),
            getTICs(spectra)
        )
        
        # Score and filter PSMs
        filter_and_score_psms!(psms, params, search_context)
        
        # Clean up temporary column
        if "filtered_scan_idx" in names(psms)
            select!(psms, Not(:filtered_scan_idx))
        end
    end
    
    return psms
end

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
Main Refactored process_file! Function
==========================================================#

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

"""
Process search results and generate QC plots for a single MS file.
"""
function process_search_results!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    try
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    mass_error_folder = getMassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Always generate plots, even with limited or no data
    # This ensures we have diagnostic output for all files
    
    # Generate RT alignment plot
    rt_plot_path = joinpath(rt_alignment_folder, parsed_fname*".pdf")
    @info "Generating RT plot for $parsed_fname: RT points = $(length(results.rt)), iRT points = $(length(results.irt))"
    if length(results.rt) > 0
        generate_rt_plot(results, rt_plot_path, parsed_fname)
        @info "Generated normal RT plot"
    else
        # Create a diagnostic plot showing fallback/borrowed status
        generate_fallback_rt_plot(results, rt_plot_path, parsed_fname, search_context, ms_file_idx)
        @info "Generated fallback RT plot"
    end
    
    # Generate mass error plot
    mass_plot_path = joinpath(mass_error_folder, parsed_fname*".pdf")
    @info "Generating mass error plot for $parsed_fname: PPM errors = $(length(results.ppm_errs))"
    if length(results.ppm_errs) > 0
        generate_mass_error_plot(results, parsed_fname, mass_plot_path)
        @info "Generated normal mass error plot"
    else
        # Create a diagnostic plot showing fallback/borrowed status
        generate_fallback_mass_error_plot(results, mass_plot_path, parsed_fname, search_context, ms_file_idx)
        @info "Generated fallback mass error plot"
    end
    
    # Clear plotting data to save memory
    resize!(results.rt, 0)
    resize!(results.irt, 0)
    resize!(results.ppm_errs, 0)
    catch e
        @warn "Failed to generate plots for file $ms_file_idx" exception=(e, catch_backtrace())
    end
end

"""
Reset results state between files.
"""
function reset_results!(results::ParameterTuningSearchResults)
    # Clear data vectors
    resize!(results.irt, 0)
    resize!(results.rt, 0)
    resize!(results.ppm_errs, 0)
    
    # Models are per-file, so they get reset at the start of each file
    # No need to reset them here
end

"""
Summarize results across all MS files.
"""
function summarize_results!(results::ParameterTuningSearchResults, params::P, search_context::SearchContext) where {P<:ParameterTuningSearchParameters}
    # Combine individual plots into merged PDFs
    try
        rt_plots_folder = getRtAlignPlotFolder(search_context)
        mass_error_plots_folder = getMassErrPlotFolder(search_context)
        qc_plots_folder = getQcPlotsFolder(results)  # Fixed: use results instead of search_context
        
        # Merge RT alignment plots
        rt_merged_path = joinpath(qc_plots_folder, "rt_alignment_combined.pdf")
        rt_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(rt_plots_folder, join=true))
        if !isempty(rt_pdf_files)
            merge_pdfs(rt_pdf_files, rt_merged_path)
            @info "Merged $(length(rt_pdf_files)) RT alignment plots to $rt_merged_path"
        else
            @info "No RT alignment plots to merge"
        end
        
        # Merge mass error plots
        mass_merged_path = joinpath(qc_plots_folder, "mass_error_combined.pdf")
        mass_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(mass_error_plots_folder, join=true))
        if !isempty(mass_pdf_files)
            merge_pdfs(mass_pdf_files, mass_merged_path)
            @info "Merged $(length(mass_pdf_files)) mass error plots to $mass_merged_path"
        else
            @info "No mass error plots to merge"
        end
        
        # Generate summary report
        # TODO: Implement generate_summary_report if detailed report needed
        # For now, the diagnostic summary below provides the key information
        
    catch e
        @warn "Failed to merge QC plots" exception=(e, catch_backtrace())
    end
    
    # Log diagnostic summary
    diagnostics = getDiagnostics(results)
    # Fixed: use values() to iterate over dictionary values
    file_statuses = values(diagnostics.file_statuses)
    converged_count = sum(s.converged for s in file_statuses)
    fallback_count = sum(s.used_fallback for s in file_statuses)
    total_files = length(file_statuses)
    
    @info "Parameter Tuning Summary:"
    @info "  - Total files: $total_files"
    @info "  - Converged: $converged_count"
    @info "  - Used fallback: $fallback_count"
    
    # Log any warnings
    for status in file_statuses
        if !isempty(status.warnings)
            @warn "File $(status.file_name) warnings:" warnings=status.warnings
        end
    end
end