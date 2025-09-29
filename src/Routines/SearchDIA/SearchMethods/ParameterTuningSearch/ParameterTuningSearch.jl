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
    :isotope_err_bounds => (3, 1),
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
Iteration State Management Functions
==========================================================#

"""
    calculate_phase_bias_shift(phase::Int64, params)

Calculate the bias shift for a given phase.
Phase 1: 0, Phase 2: +max_tolerance, Phase 3: -max_tolerance
"""
function calculate_phase_bias_shift(phase::Int64, params)::Float32
    max_tol = getMaxTolerancePpm(params)
    
    if phase == 1
        return 0.0f0  # No shift in first phase
    elseif phase == 2
        return max_tol  # Positive shift in second phase
    elseif phase == 3
        return -max_tol  # Negative shift in third phase
    else
        error("Invalid phase: $phase. Only phases 1-3 are supported.")
    end
end

"""
    reset_for_new_phase!(search_context, ms_file_idx, params, phase, iteration_state)

Reset parameters for a new phase with appropriate bias shift.
"""
function reset_for_new_phase!(search_context, ms_file_idx, params, phase::Int64, iteration_state::IterationState)
    initial_tolerance = getFragTolPpm(params)
    
    # Calculate and store bias shift for this phase (simplified - no settings needed)
    bias_shift = calculate_phase_bias_shift(phase, params)
    push!(iteration_state.phase_bias_shifts, bias_shift)
    
    # Reset to initial tolerance with new bias
    new_model = MassErrorModel(
        bias_shift,
        (initial_tolerance, initial_tolerance)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
end

"""
    update_best_attempt!(iteration_state, psm_count, mass_err_model, rt_model_data, 
                        ppm_errs, phase, score, iteration, scan_count)

Update the best attempt tracking if the current attempt has more PSMs than the previous best.
This allows us to use the best parameters as fallback when convergence fails.
"""
function update_best_attempt!(
    iteration_state::IterationState,
    psm_count::Int64,
    mass_err_model::Union{Nothing, MassErrorModel},
    rt_model_data::Union{Nothing, Tuple},
    ppm_errs::Union{Nothing, Vector{Float64}},
    phase::Int64,
    score::UInt8,
    iteration::Int64,
    scan_count::Int64
)
    # Only update if we have more PSMs than the previous best and valid models
    if psm_count > iteration_state.best_psm_count && mass_err_model !== nothing
        iteration_state.best_psm_count = psm_count
        iteration_state.best_mass_error_model = mass_err_model
        iteration_state.best_rt_model = rt_model_data
        iteration_state.best_ppm_errs = ppm_errs
        iteration_state.best_phase = phase
        iteration_state.best_score = score
        iteration_state.best_iteration = iteration
        iteration_state.best_scan_count = scan_count
    end
end

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
    
    # CRITICAL: Store RT model in SearchContext for downstream methods
    setRtIrtMap!(search_context, model[1], ms_file_idx)
    
    #parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
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
        Plots.Plot[],  # rt_plots
        Plots.Plot[],  # mass_plots
        qc_dir,
        ParameterTuningDiagnostics(),
        ParameterHistory(),
        Ref{Union{Nothing, IterationState}}(nothing)  # current_iteration_state
    )
end

#==========================================================
Helper Functions for Refactored process_file!
==========================================================#

"""
Collect PSMs from filtered spectra using library search.
Returns tuple (psms::DataFrame, was_filtered::Bool) where was_filtered indicates if PSMs passed scoring/filtering.
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
    was_filtered = false
    
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
        n_filtered = filter_and_score_psms!(psms, params, search_context)
        
        # Note: filter_and_score_psms! now always returns the actual PSM count
        # An empty DataFrame will naturally have 0 PSMs
        
        # Clean up temporary column if DataFrame not empty
        if !isempty(psms) && "filtered_scan_idx" in names(psms)
            select!(psms, Not(:filtered_scan_idx))
        end
    end
    
    return psms, was_filtered
end

"""
    collect_and_log_psms(filtered_spectra, spectra, search_context, params, ms_file_idx, context_msg::String)

Collect PSMs and log the count with context message.
Returns (psms, psm_count) where psm_count is the number of filtered PSMs (0 if filtering failed).
"""
function collect_and_log_psms(filtered_spectra, spectra, search_context, params, ms_file_idx, context_msg::String)
    psms, was_filtered = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    psm_count = size(psms, 1)
    
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
    check_and_store_convergence!(results, search_context, params, ms_file_idx, 
                                 psms, mass_err_model, ppm_errs, strategy_name::String,
                                 iteration_state::IterationState, filtered_spectra, spectra)

Check convergence criteria and store results if converged.
Optionally tests tolerance expansion if enabled.
"""
function check_and_store_convergence!(results, search_context, params, ms_file_idx, 
                                      psms, mass_err_model, ppm_errs, strategy_name::String,
                                      iteration_state::IterationState, filtered_spectra, spectra)
    if mass_err_model === nothing || psms === nothing
        return false
    end
    
    current_model = getMassErrorModel(search_context, ms_file_idx)
    
    ppm_errs_count = ppm_errs !== nothing ? length(ppm_errs) : 0

    if !check_convergence(psms, mass_err_model, current_model, ppm_errs, getMinPsms(params))
        return false
    end
    
    # Test tolerance expansion (always enabled)
    final_psms = psms
    final_model = mass_err_model
    final_ppm_errs = ppm_errs
    was_expanded = false
    
    collection_tol = iteration_state.collection_tolerance
    
    if collection_tol > 0.0f0
        final_psms, final_model, final_ppm_errs, was_expanded = test_tolerance_expansion!(
            search_context, params, ms_file_idx,
            psms, mass_err_model, ppm_errs,
            collection_tol, filtered_spectra, spectra
        )
    else
        @user_warn "Invalid collection tolerance, skipping expansion test"
    end
    
    # Store final mass error model (original or expanded)
    setMassErrorModel!(search_context, ms_file_idx, final_model)
    
    # Store final ppm errors for plotting
    if final_ppm_errs !== nothing && length(final_ppm_errs) > 0
        resize!(results.ppm_errs, 0)
        append!(results.ppm_errs, final_ppm_errs)
    end
    
    # Store RT model from final PSMs (should always have filtered PSMs at convergence)
    if !isempty(final_psms)
        rt_model_data = fit_irt_model(params, final_psms)
        set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
    else
        @user_warn "No PSMs available for RT model at convergence - this should not happen"
    end
    
    return true
end

# DEPRECATED: add_scans_for_iteration! removed - use run_single_phase_attempt and run_all_phases_with_scan_count instead

"""
    expand_mass_tolerance!(search_context, ms_file_idx, params, scale_factor::Float32 = 2.0f0)

Scale the current mass tolerance by the given factor up to maximum.
"""
function expand_mass_tolerance!(search_context, ms_file_idx, params, scale_factor::Float32 = 2.0f0)
    current_model = getMassErrorModel(search_context, ms_file_idx)
    current_left = getLeftTol(current_model)
    current_right = getRightTol(current_model)
    
    new_left = current_left * scale_factor
    new_right = current_right * scale_factor
    
    new_model = MassErrorModel(
        getMassOffset(current_model),
        (new_left, new_right)
    )
    
    setMassErrorModel!(search_context, ms_file_idx, new_model)
end


"""
    execute_strategy(strategy_num::Int, filtered_spectra, spectra, search_context, params, ms_file_idx, iteration_state)

Execute one of two convergence strategies:
- Strategy 1: Try with current parameters
- Strategy 2: Expand mass tolerance AND adjust bias
"""
# execute_strategy function removed - logic now inline in run_single_phase

"""
    initialize_models!(search_context, ms_file_idx, params)

Initialize mass error and quad transmission models for file.
"""
function initialize_models!(search_context, ms_file_idx, params)
    # Use fixed initial parameters from JSON configuration
    initial_bias = 0.0f0  # Always start with zero bias
    initial_tolerance = getFragTolPpm(params)
    
    # Set initial mass error model
    setMassErrorModel!(search_context, ms_file_idx, MassErrorModel(
        initial_bias,
        (initial_tolerance, initial_tolerance)
    ))
    
    # Set initial quad transmission model
    setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))
end


"""
    store_final_results!(results, search_context, params, ms_file_idx, converged, n_attempts, final_psm_count, iteration_state)

Store tuning results and diagnostic information, including best attempt fallback status.
"""
function store_final_results!(results, search_context, params, ms_file_idx, 
                              converged, n_attempts, final_psm_count, iteration_state::IterationState)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Determine convergence type
    convergence_type = if converged
        "CONVERGED"
    elseif iteration_state.best_psm_count > 0
        "BEST_ATTEMPT ($(iteration_state.best_psm_count) PSMs)"
    else
        "FAILED"
    end
    
    # Record tuning results
    tuning_result = TuningResults(
        getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        (getLeftTol(getMassErrorModel(search_context, ms_file_idx)), 
         getRightTol(getMassErrorModel(search_context, ms_file_idx))),
        converged,
        final_psm_count,
        n_attempts
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
         getRightTol(getMassErrorModel(search_context, ms_file_idx)))
    )
    record_tuning_status!(results.diagnostics, status)
    
    # Store final mass error model in results
    results.mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
end



#==========================================================
New Scan Scaling Functions
==========================================================#

"""
    run_single_phase_attempt(phase::Int64, scan_count::Int64, iteration_state, 
                            results, params, search_context, ms_file_idx, spectra)

Run a single phase (with configured iterations) using fixed scan count.
Returns true if converged, false otherwise.
"""
function run_single_phase(
    phase::Int64,
    filtered_spectra::FilteredMassSpecData,
    iteration_state::IterationState,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    settings = getIterationSettings(params)
    
    # Set up phase
    iteration_state.current_phase = phase
    iteration_state.current_iteration_in_phase = 0
    
    # Get all score thresholds to try
    score_thresholds = getMinIndexSearchScores(params)
    psm_count = 0
    # NEW: Score threshold loop INSIDE each phase
    for (score_idx, min_score) in enumerate(score_thresholds)
        
        # Update current score threshold
        setCurrentMinScore!(params, min_score)
        
        # Apply phase bias shift and reset tolerance for this score
        reset_for_new_phase!(search_context, ms_file_idx, params, phase, iteration_state)
    
        # ========================================
        # INITIAL ATTEMPT (before iteration loop)
        # ========================================
        
        # Collect PSMs at initial tolerance with phase bias and current score
        psms_initial, _ = collect_and_log_psms(
            filtered_spectra, spectra, search_context, 
            params, ms_file_idx, "initial attempt with min_score=$min_score"
        )

        # Fit models and check initial convergence
        mass_err_model, ppm_errs, psm_count = fit_models_from_psms(
            psms_initial, spectra, search_context, params, ms_file_idx
        )
        
        # Track collection tolerance if we have a model
        if mass_err_model !== nothing
            current_model = getMassErrorModel(search_context, ms_file_idx)
            iteration_state.collection_tolerance = (getLeftTol(current_model) + getRightTol(current_model)) / 2.0f0
            
            # Track best attempt even if not converged
            if psm_count > 0 && !isempty(psms_initial)
                # Fit RT model for best attempt tracking (only if we have filtered PSMs)
                rt_model_data = fit_irt_model(params, psms_initial)
                
                # Update best attempt with filtered PSM count
                update_best_attempt!(
                    iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
                    phase, min_score, 0, iteration_state.current_scan_count
                )
            end
        end
        
        if mass_err_model !== nothing && check_and_store_convergence!(
            results, search_context, params, ms_file_idx,
            psms_initial, mass_err_model, ppm_errs,
            "Initial attempt (Phase $phase, min_score=$min_score)",
            iteration_state, filtered_spectra, spectra
        )
            return true
        end
        
        # ========================================
        # ITERATION LOOP (expand → collect → adjust → collect cycles)
        # ========================================
        
        for iter in 1:settings.iterations_per_phase
            iteration_state.current_iteration_in_phase = iter
            iteration_state.total_iterations += 1
            
            # Step 1: EXPAND TOLERANCE
            expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                  settings.mass_tolerance_scale_factor)
            current_tolerance = settings.init_mass_tol_ppm * 
                              (settings.mass_tolerance_scale_factor ^ iter)
            current_model = getMassErrorModel(search_context, ms_file_idx)

            # Step 2: COLLECT PSMs with expanded tolerance (to determine bias)
            psms_for_bias, _ = collect_and_log_psms(
                filtered_spectra, spectra, search_context,
                params, ms_file_idx, "with expanded tolerance and min_score=$min_score"
            )
            
            # Step 3: FIT MODEL to determine bias adjustment
            mass_err_for_bias, _, _ = fit_models_from_psms(
                psms_for_bias, spectra, search_context, params, ms_file_idx
            )
            
            # Step 4: ADJUST BIAS based on fitted model
            if mass_err_for_bias !== nothing
                new_bias = getMassOffset(mass_err_for_bias)
                current_model = getMassErrorModel(search_context, ms_file_idx)
                old_bias = getMassOffset(current_model)
                
                # Update the bias while keeping the expanded tolerance
                adjusted_model = MassErrorModel(
                    new_bias,
                    (getLeftTol(current_model), getRightTol(current_model))
                )
                setMassErrorModel!(search_context, ms_file_idx, adjusted_model)
            end
            
            # Step 5: COLLECT PSMs again with adjusted bias
            psms_adjusted, _ = collect_and_log_psms(
                filtered_spectra, spectra, search_context,
                params, ms_file_idx, "with adjusted bias and min_score=$min_score"
            )
            
            # Step 6: FIT FINAL MODELS with bias-adjusted PSMs
            mass_err_model, ppm_errs, psm_count = fit_models_from_psms(
                psms_adjusted, spectra, search_context, params, ms_file_idx
            )
            
            # Update collection tolerance tracking
            current_model = getMassErrorModel(search_context, ms_file_idx)
            iteration_state.collection_tolerance = (getLeftTol(current_model) + getRightTol(current_model)) / 2.0f0
            
            # Track best attempt after each iteration
            if mass_err_model !== nothing && psm_count > 0 && !isempty(psms_adjusted)
                # Fit RT model for best attempt tracking (only if we have filtered PSMs)
                rt_model_data = fit_irt_model(params, psms_adjusted)
                
                # Update best attempt with filtered PSM count
                update_best_attempt!(
                    iteration_state, psm_count, mass_err_model, rt_model_data, ppm_errs,
                    phase, min_score, iter, iteration_state.current_scan_count
                )
            end
            
            # Step 7: CHECK CONVERGENCE
            if check_and_store_convergence!(
                results, search_context, params, ms_file_idx,
                psms_adjusted, mass_err_model, ppm_errs,
                "Phase $phase, Iteration $iter with min_score=$min_score",
                iteration_state, filtered_spectra, spectra
            )
                return true
            end
        end
    end  # End of score threshold loop

    return false
end

"""
    run_all_phases_with_scan_count(scan_count, iteration_state, results, 
                                   params, search_context, ms_file_idx, spectra)

Run all 3 phases with a fixed scan count.
Returns true if any phase achieves convergence.
"""
function run_all_phases_with_scan_count(
    filtered_spectra::FilteredMassSpecData,
    iteration_state::IterationState,
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    
    # Update iteration state
    iteration_state.current_scan_count = length(filtered_spectra)
    # Run all 3 phases with the same filtered_spectra
    for phase in 1:3
        converged = run_single_phase(
            phase, filtered_spectra, iteration_state,
            results, params, search_context, ms_file_idx, spectra
        )
        
        if converged
            return true
        end
    end
    @user_warn "All 3 phases completed without convergence"
    return false
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
    
    # Check if file should be skipped due to previous failure
    if is_file_failed(search_context, ms_file_idx)
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        @user_warn "Skipping ParameterTuningSearch for previously failed file: $file_name"
        return results  # Return early with unchanged results
    end
    
    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    final_psm_count = 0
    
    # Initialize iteration state
    iteration_state = IterationState()
    settings = getIterationSettings(params)
    
    # Get scan count parameters
    scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    scan_scale_factor = settings.scan_scale_factor
    
    # Define filtered_spectra outside try block for use in fallback
    filtered_spectra = nothing
    
    try
        
        # Initialize models
        initialize_models!(search_context, ms_file_idx, params)
        
        # Create filtered spectra ONCE with initial scan count
        try
            filtered_spectra = FilteredMassSpecData(
                spectra,
                max_scans = scan_count,
                topn = something(getTopNPeaks(params), 200),
                target_ms_order = UInt8(2)
            )
            
            # Check if we have any usable scans
            if length(filtered_spectra) == 0
                file_name = try
                    getFileIdToName(getMSData(search_context), ms_file_idx)
                catch
                    "file_$ms_file_idx"
                end
                @user_warn "MS data file $file_name contains no usable MS2 scans for parameter tuning. Using conservative defaults."
                # Force fallback by setting converged = false and raising an internal exception flag
                iteration_state.failed_with_exception = true
                # Skip to fallback logic by jumping out of the try block
                throw(ErrorException("No usable MS2 scans"))
            end
        catch e
            throw(e)
            if isa(e, ErrorException) && e.msg == "No usable MS2 scans"
                # Re-throw our controlled exception to handle in outer catch
                rethrow(e)
            else
                # Handle other FilteredMassSpecData creation errors
                file_name = try
                    getFileIdToName(getMSData(search_context), ms_file_idx)
                catch
                    "file_$ms_file_idx"
                end
                @user_warn "Failed to create filtered spectra for MS data file: $file_name. Error type: $(typeof(e)). Using conservative defaults."
                # Log full stacktrace for diagnostics
                try
                    bt = catch_backtrace()
                    @user_error sprint(showerror, e, bt)
                catch
                end
                iteration_state.failed_with_exception = true
                throw(e)  # Re-throw to be caught by outer catch block
            end
        end
        
        # Track current scan count
        current_scan_count = scan_count
        attempt_count = 0
        
        # Main scan scaling loop
        while true
            attempt_count += 1
            iteration_state.scan_attempt = attempt_count
            
            # Run all phases with current filtered_spectra
            converged = run_all_phases_with_scan_count(
                filtered_spectra, iteration_state, results,
                params, search_context, ms_file_idx, spectra
            )
            
            if converged
                iteration_state.converged = true
                break
            end
            
            # Check if we've reached or exceeded max scans
            if current_scan_count >= max_scans
                iteration_state.max_scan_count_reached = true
                break
            end
            
            # Calculate next scan count
            next_scan_count = Int64(ceil(current_scan_count * scan_scale_factor))
            
            # Check if next iteration would exceed max
            if next_scan_count > max_scans
                if current_scan_count < max_scans
                    # Do one final attempt with exactly max_scans
                    additional_scans = max_scans - current_scan_count
                    append!(filtered_spectra; max_additional_scans = additional_scans)
                    current_scan_count = max_scans
                    # Loop will continue for one more attempt
                else
                    # Already at max, break
                    break
                end
            else
                # Normal scaling - append more scans
                additional_scans = next_scan_count - current_scan_count
                append!(filtered_spectra; max_additional_scans = additional_scans)
                current_scan_count = next_scan_count
            end
        end
        #@debug_l1 "manual set mass err model . "
        #setMassErrorModel!(search_context, ms_file_idx, 
        #MassErrorModel(getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        #(7.0f0, 7.0f0)))
        # Get final PSM count if converged
        if converged
            final_psm_count = size(results.rt, 1)  # Track from results
        end
        
        converged_score = converged ? getMinIndexSearchScore(params) : nothing
        if !converged
            @user_warn "Completed $(iteration_state.total_iterations) total iterations " *
                   "across $(iteration_state.scan_attempt) attempt(s). " *
                   "✗ Did not converge with any score threshold: $(getMinIndexSearchScores(params))"
        end
        
    catch e
        # Get the actual file name for clear error reporting
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        
        # Don't mark as failed - just continue with defaults
        @user_warn "Parameter tuning failed for MS data file: $file_name. Error type: $(typeof(e)). Using conservative default parameters to continue analysis."
        # Log full stacktrace for diagnostics
        try
            bt = catch_backtrace()
            @user_error sprint(showerror, e, bt)
        catch
        end
        # During debugging, rethrow to stop the pipeline and surface the root cause
        rethrow(e)
        
        # Set conservative defaults to allow pipeline to continue
        converged = false
        iteration_state.failed_with_exception = true
        
        # Set default IRT error to infinite tolerance (no IRT filtering)
        getIrtErrors(search_context)[ms_file_idx] = typemax(Float32)
        
        # Clear any partial state
        if filtered_spectra !== nothing
            # If we got this far, we had some data structure initialized
        end
    end
    
    # Store results and handle fallback if needed
    if !converged
        # Check if we have a best attempt to use
        if iteration_state.best_mass_error_model !== nothing
            left_tol = round(getLeftTol(iteration_state.best_mass_error_model), digits = 1)
            right_tol = round(getRightTol(iteration_state.best_mass_error_model), digits = 1)
            @user_warn "Failed to converge for file $ms_file_idx after $(iteration_state.scan_attempt) attempts. \n" *
                  "Using best iteration: Phase $(iteration_state.best_phase) \n" *
                  "(bias=$(round(calculate_phase_bias_shift(iteration_state.best_phase, params), digits=1)) ppm), \n" *
                  "Score threshold $(iteration_state.best_score), \n" *
                  "Iteration $(iteration_state.best_iteration). \n" *
                  "Yielded $(iteration_state.best_psm_count) PSMs with \n" *
                  "mass offset $(round(getMassOffset(iteration_state.best_mass_error_model), digits=1)) ppm and \n" *
                  "tolerance -$left_tol/+$right_tol ppm \n"
            
            # Apply best attempt models
            setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
            
            if iteration_state.best_rt_model !== nothing
                set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                    iteration_state.best_rt_model)
            else
                # If no RT model, use identity
                setRtIrtMap!(search_context, IdentityModel(), ms_file_idx)
                results.rt_to_irt_model[] = IdentityModel()
            end
            
            # Store best ppm errors for plotting
            if iteration_state.best_ppm_errs !== nothing
                resize!(results.ppm_errs, 0)
                append!(results.ppm_errs, iteration_state.best_ppm_errs)
            end
            
            # Test 1.5x tolerance expansion on best iteration (only if we have filtered_spectra)
            if filtered_spectra !== nothing
                @debug_l1 "Testing 1.5x tolerance expansion on best iteration parameters..."
                expanded_model = MassErrorModel(
                    getMassOffset(iteration_state.best_mass_error_model),
                    (getLeftTol(iteration_state.best_mass_error_model) * 1.5f0,
                     getRightTol(iteration_state.best_mass_error_model) * 1.5f0)
                )
                
                # Apply expanded model and collect PSMs
                setMassErrorModel!(search_context, ms_file_idx, expanded_model)
                expanded_psms, expanded_ppm_errs = collect_psms_with_model(
                    filtered_spectra, search_context, params, ms_file_idx, spectra
                )
                
                if size(expanded_psms, 1) > iteration_state.best_psm_count
                    # Refit model with expanded PSMs
                    if length(expanded_ppm_errs) > 0
                        # Calculate mass error from expanded PSMs
                        expanded_fragments = get_matched_fragments(
                            spectra, expanded_psms, search_context, params, ms_file_idx
                        )
                        if length(expanded_fragments) > 0
                            refitted_model, refitted_ppm_errs = fit_mass_err_model(params, expanded_fragments)
                            if refitted_model !== nothing
                                # Use expanded results
                                psm_increase = size(expanded_psms, 1) - iteration_state.best_psm_count
                                percent_increase = round(100 * psm_increase / iteration_state.best_psm_count, digits=1)
                                
                                @user_info "Tolerance expansion yielded $(size(expanded_psms, 1)) PSMs " *
                                          "(+$psm_increase, $percent_increase% increase)"
                                
                                iteration_state.best_mass_error_model = refitted_model
                                iteration_state.best_psm_count = size(expanded_psms, 1)
                                iteration_state.best_ppm_errs = refitted_ppm_errs
                                
                                # Update results with expanded data
                                results.mass_err_model[] = refitted_model
                                resize!(results.ppm_errs, 0)
                                append!(results.ppm_errs, refitted_ppm_errs)
                                setMassErrorModel!(search_context, ms_file_idx, refitted_model)
                                
                                @user_info "EXPANDED_BEST_ATTEMPT: Tolerance expansion increased PSMs by $percent_increase%"
                            else
                                # Revert to best iteration model if refitting failed
                                setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
                                results.mass_err_model[] = iteration_state.best_mass_error_model
                            end
                        else
                            # No fragments, revert to best iteration model
                            setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
                            results.mass_err_model[] = iteration_state.best_mass_error_model
                        end
                    else
                        # No PPM errors, revert to best iteration model
                        setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
                        results.mass_err_model[] = iteration_state.best_mass_error_model
                    end
                else
                    # Expansion didn't help, revert to best iteration model
                    setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)
                    results.mass_err_model[] = iteration_state.best_mass_error_model
                    @debug_l1 "Tolerance expansion did not significantly improve PSM count"
                end
            else
                # Could not test tolerance expansion - no filtered spectra available
                @user_warn "Cannot test tolerance expansion - filtered spectra not available"
            end
            
            # Build detailed warning message
            left_tol = round(getLeftTol(iteration_state.best_mass_error_model), digits = 1)
            right_tol = round(getRightTol(iteration_state.best_mass_error_model), digits = 1)
            final_psm_count = iteration_state.best_psm_count
            
        else
            # Determine why we failed to provide appropriate messaging
            if iteration_state.failed_with_exception
                file_name = try
                    getFileIdToName(getMSData(search_context), ms_file_idx)
                catch
                    "file_$ms_file_idx"
                end
                @user_warn "Processing failed with exception for MS data file: $file_name. Using conservative default parameters to continue analysis."
            else
                # No valid attempts found any PSMs - use conservative defaults or borrow
                @user_warn "Failed to converge for file $ms_file_idx after $(iteration_state.scan_attempt) attempts. " *
                          "No attempts found sufficient PSMs (best was $(iteration_state.best_psm_count)), using fallback strategy"
            end
            
            fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
                search_context, ms_file_idx
            )
            
            setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)
            setRtIrtMap!(search_context, fallback_rt_model, ms_file_idx)
            
            # Update results with fallback
            results.mass_err_model[] = fallback_mass_err
            results.rt_to_irt_model[] = fallback_rt_model
            
            if borrowed_from === nothing
                @user_info "CONSERVATIVE_FALLBACK: No valid attempts, used defaults"
            end
        end
    end
    
    # Record tuning status
    store_final_results!(
        results, search_context, params, ms_file_idx,
        converged, iteration_state.scan_attempt, final_psm_count, iteration_state
    )
    
    # Add to diagnostics
    record_file_status!(
        results.diagnostics, ms_file_idx, parsed_fname,
        converged, !converged, iteration_state
    )

    # Store iteration_state in results for use in process_search_results!
    results.current_iteration_state[] = iteration_state

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
        
        # Get iteration_state from results
        iteration_state = results.current_iteration_state[]
        # Note: No mass-error buffer is applied. Plots reflect the fitted model.
        
        # Always generate plots, even with limited or no data
        # This ensures we have diagnostic output for all files
        
        # Generate RT alignment plot (only store in memory, no individual files)
        if length(results.rt) > 0
            # If we have RT data, generate regular plot
            rt_plot = generate_rt_plot(results, parsed_fname)
            push!(results.rt_plots, rt_plot)  # Store for combined PDF
        elseif iteration_state !== nothing && iteration_state.best_rt_model !== nothing
            # Use best iteration data if available
            rt_plot = generate_best_iteration_rt_plot_in_memory(results, parsed_fname, iteration_state)
            push!(results.rt_plots, rt_plot)
        else
            # Create a diagnostic plot showing fallback/borrowed status
            fallback_plot = generate_fallback_rt_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(results.rt_plots, fallback_plot)  # Store for combined PDF
            end
        end
        
        # Generate mass error plot (only store in memory, no individual files)
        m = getMassErrorModel(search_context, ms_file_idx)
        #buffered_model = MassErrorModel(
        #    getMassOffset(m),
        #    (getLeftTol(m) + 1.0f0, getRightTol(m) + 1.0f0)
        #)
        buffered_model = m  # No additional buffer applied
        setMassErrorModel!(
            search_context,
            ms_file_idx,
            buffered_model)

        results.mass_err_model[] = buffered_model  # Update results with buffered model

        if length(results.ppm_errs) > 0
            # If we have mass error data, generate regular plot
            mass_plot = generate_mass_error_plot(results, parsed_fname)
            push!(results.mass_plots, mass_plot)  # Store for combined PDF
        elseif iteration_state !== nothing && iteration_state.best_ppm_errs !== nothing
            # Use best iteration data if available
            mass_plot = generate_best_iteration_mass_error_plot_in_memory(results, parsed_fname, iteration_state)
            push!(results.mass_plots, mass_plot)
        else
            # Create a diagnostic plot showing fallback/borrowed status
            fallback_plot = generate_fallback_mass_error_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(results.mass_plots, fallback_plot)  # Store for combined PDF
            end
        end
        
        # Update models in search context
        setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
        setRtIrtMap!(search_context, getRtToIrtModel(results), ms_file_idx)
        
        # Clear plotting data to save memory
        resize!(results.rt, 0)
        resize!(results.irt, 0)
        resize!(results.ppm_errs, 0)
        results.current_iteration_state[] = nothing  # Clear iteration state after use
    catch e
        # Plot-generation failures are non-fatal and should not mark the file as failed.
        # Keep downstream processing intact; log and continue.
        @user_warn "Failed to generate plots for file $ms_file_idx" exception=(e, catch_backtrace())
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

# Deprecated: apply_final_mass_error_buffer! removed – buffer is applied per-file before plotting

"""
Summarize results across all MS files.
"""
function summarize_results!(results::ParameterTuningSearchResults, params::P, search_context::SearchContext) where {P<:ParameterTuningSearchParameters}
    # Combine individual plots into merged PDFs using save_multipage_pdf
    
    try
        rt_plots_folder = getRtAlignPlotFolder(search_context)
        mass_error_plots_folder = getMassErrPlotFolder(search_context)
        
        # Create combined RT alignment PDF from collected plots
        if !isempty(results.rt_plots)
            rt_combined_path = joinpath(rt_plots_folder, "rt_alignment_plots.pdf")
            try
                if isfile(rt_combined_path)
                    safeRm(rt_combined_path, nothing)
                end
            catch e
                @user_warn "Could not clear existing RT plots file: $e"
            end
            save_multipage_pdf(results.rt_plots, rt_combined_path)
            empty!(results.rt_plots)  # Clear to free memory
        end
        
        # Create combined mass error PDF from collected plots
        if !isempty(results.mass_plots)
            mass_combined_path = joinpath(mass_error_plots_folder, "mass_error_plots.pdf")
            try
                if isfile(mass_combined_path)
                    safeRm(mass_combined_path, nothing)
                end
            catch e
                @user_warn "Could not clear existing mass error plots file: $e"
            end
            save_multipage_pdf(results.mass_plots, mass_combined_path)
            empty!(results.mass_plots)  # Clear to free memory
        end
        
        # Generate summary report
        # TODO: Implement generate_summary_report if detailed report needed
        # For now, the diagnostic summary below provides the key information
        
    catch e
        @user_warn "Failed to merge QC plots" exception=(e, catch_backtrace())
    end
    
    # Apply buffer to all mass error models AFTER plots are generated
    # This ensures plots show the actual fitted values, not buffered ones
    # Do not apply an additional buffer here; models were buffered once per file
    
    # Log diagnostic summary
    diagnostics = getDiagnostics(results)
    # Fixed: use values() to iterate over dictionary values
    file_statuses = values(diagnostics.file_statuses)
    converged_count = sum(s.converged for s in file_statuses)
    fallback_count = sum(s.used_fallback for s in file_statuses)
    total_files = length(file_statuses)
    
    @user_info "Parameter Tuning Summary:" * "\n" *
               "  - Total files processed: $total_files" * "\n" *
               "  - Converged: $converged_count" * "\n" *
               "  - Used fallback: $fallback_count"
    
end
