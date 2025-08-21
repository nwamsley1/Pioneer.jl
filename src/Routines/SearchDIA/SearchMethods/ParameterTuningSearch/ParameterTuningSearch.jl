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
    
    @debug_l1 "Phase $phase: Reset to initial tolerance (±$(round(initial_tolerance, digits=1)) ppm) " *
           "with bias shift $(round(bias_shift, digits=1)) ppm"
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
        
        @debug_l1 "    📊 New best attempt: $(psm_count) PSMs " *
              "(Phase $phase, Score $score, Iteration $iteration)"
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
        @debug_l1 "Borrowing parameters from file $borrowed_fname (index $borrowed_from) for file $ms_file_idx"
        push!(warnings, "BORROWED: Using parameters from file $borrowed_fname")
        return fallback_mass_err, fallback_rt_model, borrowed_from
    else
        # Use conservative defaults if no other files available
        @debug_l1 "No successfully tuned files available. Using conservative defaults for file $ms_file_idx"
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
        Plots.Plot[],  # rt_plots
        Plots.Plot[],  # mass_plots
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
        
        # Score and filter PSMs - check if filtering succeeded
        n_filtered = filter_and_score_psms!(psms, params, search_context)
        was_filtered = n_filtered != -1
        
        # If filtering failed, return empty DataFrame
        if !was_filtered
            psms = DataFrame()
        end
        
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
    
    if was_filtered
        @debug_l1 "Collected $psm_count filtered PSMs $context_msg"
    else
        @debug_l1 "Collected 0 filtered PSMs $context_msg (filtering failed)"
    end
    
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
    @user_warn "size(psms) = $(size(psms, 1)), " *
           "mass_err_model = $(mass_err_model), " *
           "ppm_errs = $(ppm_errs_count)" *
           " min_psms = $(getMinPsms(params))"
    if !check_convergence(psms, mass_err_model, current_model, ppm_errs, getMinPsms(params))
        return false
    end
    
    @debug_l1 "Converged after $strategy_name"
    
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
        
        if was_expanded
            @debug_l1 "Using expanded tolerance results for final model"
        end
    else
        @user_warn "Invalid collection tolerance, skipping expansion test"
    end
    
    # Store final mass error model (original or expanded)
    setMassErrorModel!(search_context, ms_file_idx, final_model)
    
    # Store final ppm errors for plotting
    if final_ppm_errs !== nothing && length(final_ppm_errs) > 0
        resize!(results.ppm_errs, 0)
        append!(results.ppm_errs, final_ppm_errs)
        @debug_l1 "Stored $(length(results.ppm_errs)) ppm errors for plotting"
    end
    
    # Store RT model from final PSMs (should always have filtered PSMs at convergence)
    if !isempty(final_psms)
        rt_model_data = fit_irt_model(params, final_psms)
        set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
        @debug_l1 "Stored $(length(results.rt)) RT points and $(length(results.irt)) iRT points"
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
    
    @debug_l1 "Scaled tolerance by factor $(scale_factor): " *
           "from ($(round(current_left, digits=1)), $(round(current_right, digits=1))) " *
           "to ($(round(new_left, digits=1)), $(round(new_right, digits=1))) ppm"
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
    
    @debug_l1 "Initial mass error model for file $ms_file_idx: $(getMassErrorModel(search_context, ms_file_idx))"
end


"""
    store_final_results!(results, search_context, params, ms_file_idx, converged, n_attempts, final_psm_count, warnings, iteration_state)

Store tuning results and diagnostic information, including best attempt fallback status.
"""
function store_final_results!(results, search_context, params, ms_file_idx, 
                              converged, n_attempts, final_psm_count, warnings, iteration_state::IterationState)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Determine convergence type
    convergence_type = if converged
        "CONVERGED"
    elseif iteration_state.best_psm_count > 0
        "BEST_ATTEMPT ($(iteration_state.best_psm_count) PSMs)"
    else
        "FAILED"
    end
    
    # Log final status
    @debug_l1 "Final parameter tuning status for file $ms_file_idx ($parsed_fname):"
    @debug_l1 "  - Status: $convergence_type"
    @debug_l1 "  - RT points: $(length(results.rt))"
    @debug_l1 "  - iRT points: $(length(results.irt))" 
    @debug_l1 "  - PPM errors: $(length(results.ppm_errs))"
    @debug_l1 "  - Final PSM count: $final_psm_count"
    
    if !converged && iteration_state.best_psm_count > 0
        @debug_l1 "  - Best attempt: Phase $(iteration_state.best_phase), " *
              "Score $(iteration_state.best_score), " *
              "Iteration $(iteration_state.best_iteration), " *
              "Scans: $(iteration_state.best_scan_count)"
    end
    
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
    
    @debug_l1 "="^50
    @debug_l1 "Beginning Phase $phase of Attempt $(iteration_state.scan_attempt)"
    @debug_l1 "Phase bias shift: $(calculate_phase_bias_shift(phase, params)) ppm"
    @debug_l1 "="^50
    
    # Get all score thresholds to try
    score_thresholds = getMinIndexSearchScores(params)
    
    # NEW: Score threshold loop INSIDE each phase
    for (score_idx, min_score) in enumerate(score_thresholds)
        @debug_l1 ""
        @debug_l1 "  Phase $phase - Trying min_score = $min_score (threshold $score_idx of $(length(score_thresholds)))"
        @debug_l1 "  " * "-"^40
        
        # Update current score threshold
        setCurrentMinScore!(params, min_score)
        
        # Apply phase bias shift and reset tolerance for this score
        reset_for_new_phase!(search_context, ms_file_idx, params, phase, iteration_state)
    
        # ========================================
        # INITIAL ATTEMPT (before iteration loop)
        # ========================================
        
        @debug_l1 "    Initial attempt at base tolerance ($(settings.init_mass_tol_ppm) ppm) with min_score=$min_score"
        
        # Collect PSMs at initial tolerance with phase bias and current score
        psms_initial, _ = collect_and_log_psms(
            filtered_spectra, spectra, search_context, 
            params, ms_file_idx, "initial attempt with min_score=$min_score"
        )
        
        # Fit models and check initial convergence
        mass_err_model, ppm_errs, psm_count = fit_models_from_psms(
            psms_initial, spectra, search_context, params, ms_file_idx
        )
        
        # Only log MAD if we have ppm_errs
        if ppm_errs !== nothing
            @user_warn "    mass_err_model $mass_err_model, mad(ppm_errs) $(mad(ppm_errs)), min_score=$min_score"
        else
            @user_warn "    mass_err_model $mass_err_model, no ppm_errs available (psm_count: $psm_count), min_score=$min_score"
        end
        
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
            @debug_l1 "  ✓ Converged in Phase $phase with min_score=$min_score (initial tolerance)"
            return true
        end
        
        @debug_l1 "    Initial attempt with min_score=$min_score did not converge, starting expand/adjust iterations"
        
        # ========================================
        # ITERATION LOOP (expand → collect → adjust → collect cycles)
        # ========================================
        
        for iter in 1:settings.iterations_per_phase
            iteration_state.current_iteration_in_phase = iter
            iteration_state.total_iterations += 1
            
            @debug_l1 ""
            @debug_l1 "    Phase $phase, Score $min_score, Iteration $iter (Total: $(iteration_state.total_iterations))"
            
            # Step 1: EXPAND TOLERANCE
            expand_mass_tolerance!(search_context, ms_file_idx, params, 
                                  settings.mass_tolerance_scale_factor)
            current_tolerance = settings.init_mass_tol_ppm * 
                              (settings.mass_tolerance_scale_factor ^ iter)
            current_model = getMassErrorModel(search_context, ms_file_idx)
            @debug_l1 "      Step 1: Expanded tolerance to $(round(current_tolerance, digits=1)) ppm while bias is $(getMassOffset(current_model)) ppm"
            
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
                
                @debug_l1 "      Step 4: Adjusted bias from $(round(old_bias, digits=2)) to $(round(new_bias, digits=2)) ppm"
            else
                @debug_l1 "      Step 4: Could not fit model for bias adjustment, keeping current bias"
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
                @debug_l1 "  ✓ Converged in Phase $phase with min_score=$min_score (iteration $iter)"
                return true
            end
            
            @debug_l1 "      Iteration $iter with min_score=$min_score did not converge"
        end
        
        @debug_l1 "  Phase $phase with min_score=$min_score did not converge"
    end  # End of score threshold loop
    
    @debug_l1 "Phase $phase completed without convergence (tried all score thresholds)"
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
    @debug_l1 "="^60
    @debug_l1 "Starting scan attempt $(iteration_state.scan_attempt) with $(length(filtered_spectra)) scans"
    @debug_l1 "="^60
    
    # Update iteration state
    iteration_state.current_scan_count = length(filtered_spectra)
    
    # Run all 3 phases with the same filtered_spectra
    for phase in 1:3
        converged = run_single_phase(
            phase, filtered_spectra, iteration_state,
            results, params, search_context, ms_file_idx, spectra
        )
        
        if converged
            @debug_l1 "✓ Converged in Phase $phase of Attempt $(iteration_state.scan_attempt)!"
            return true
        end
    end
    
    @debug_l1 "All 3 phases completed without convergence"
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
    
    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    warnings = String[]
    final_psm_count = 0
    
    # Initialize iteration state
    iteration_state = IterationState()
    settings = getIterationSettings(params)
    
    # Get scan count parameters
    scan_count = getInitialScanCount(params)
    max_scans = getMaxParameterTuningScans(params)
    scan_scale_factor = settings.scan_scale_factor
    
    try
         @debug_l1 "Processing file: $parsed_fname (index: $ms_file_idx)"
         @debug_l1 "Scan scaling strategy: initial=$scan_count, max=$max_scans, scale_factor=$scan_scale_factor"
         @debug_l1 "Mass tolerance scaling: init=$(settings.init_mass_tol_ppm) ppm, " *
               "scale_factor=$(settings.mass_tolerance_scale_factor), " *
               "iterations_per_phase=$(settings.iterations_per_phase)"
        
        # Initialize models
        initialize_models!(search_context, ms_file_idx, params)
        
        # Create filtered spectra ONCE with initial scan count
        filtered_spectra = FilteredMassSpecData(
            spectra,
            max_scans = scan_count,
            topn = something(getTopNPeaks(params), 200),
            target_ms_order = UInt8(2)
        )
        
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
                @debug_l1 "SUCCESS: Converged after $attempt_count attempt(s) " *
                #       "with $current_scan_count scans"
                break
            end
            
            # Check if we've reached or exceeded max scans
            if current_scan_count >= max_scans
                @user_warn "Reached maximum scan count ($max_scans) without convergence"
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
                    @debug_l1 "Adding final $additional_scans scans to reach maximum of $max_scans"
                    
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
                @debug_l1 "Appending $additional_scans scans (total will be $next_scan_count)"
                
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
        if converged
            @debug_l1 "Completed $(iteration_state.total_iterations) total iterations " *
                   "across $(iteration_state.scan_attempt) attempt(s). " *
                   "✓ CONVERGED with min_score=$converged_score"
            push!(warnings, "Converged with min_score=$converged_score")
        else
            @debug_l1 "Completed $(iteration_state.total_iterations) total iterations " *
                   "across $(iteration_state.scan_attempt) attempt(s). " *
                   "✗ Did not converge with any score threshold: $(getMinIndexSearchScores(params))"
        end
        
    catch e
        @user_warn "Parameter tuning failed for file $parsed_fname with error" exception=e
        push!(warnings, "ERROR: $(string(e))")
    end
    
    # Store results and handle fallback if needed
    if !converged
        @user_warn "Failed to converge for file $ms_file_idx after " *
              "$(iteration_state.scan_attempt) attempts"
        
        # Check if we have a best attempt to use
        if iteration_state.best_mass_error_model !== nothing
            @debug_l1 "Using best attempt as fallback: " *
                  "$(iteration_state.best_psm_count) PSMs from " *
                  "Phase $(iteration_state.best_phase), " *
                  "Score $(iteration_state.best_score), " *
                  "Iteration $(iteration_state.best_iteration)"
            
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
            
            # Update results with best attempt
            results.mass_err_model[] = iteration_state.best_mass_error_model
            
            # Update warnings and final PSM count
            push!(warnings, "BEST_ATTEMPT_FALLBACK: Used parameters from best attempt " *
                           "($(iteration_state.best_psm_count) PSMs, " *
                           "Phase $(iteration_state.best_phase), " *
                           "Score $(iteration_state.best_score))")
            
            final_psm_count = iteration_state.best_psm_count
            
        else
            # No valid attempts at all - use conservative defaults or borrow
            @user_warn "No valid attempts found, using fallback strategy"
            
            fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
                search_context, ms_file_idx, warnings
            )
            
            setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)
            setRtIrtMap!(search_context, fallback_rt_model, ms_file_idx)
            
            # Update results with fallback
            results.mass_err_model[] = fallback_mass_err
            results.rt_to_irt_model[] = fallback_rt_model
            
            if borrowed_from === nothing
                push!(warnings, "CONSERVATIVE_FALLBACK: No valid attempts, used defaults")
            end
        end
    end
    
    # Record tuning status
    store_final_results!(
        results, search_context, params, ms_file_idx,
        converged, iteration_state.scan_attempt, final_psm_count, warnings, iteration_state
    )
    
    # Add to diagnostics
    record_file_status!(
        results.diagnostics, ms_file_idx, parsed_fname,
        converged, !converged, warnings, iteration_state
    )

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
        
        # Generate RT alignment plot (only store in memory, no individual files)
        @debug_l1 "Generating RT plot for $parsed_fname: RT points = $(length(results.rt)), iRT points = $(length(results.irt))"
        if length(results.rt) > 0
            rt_plot = generate_rt_plot(results, parsed_fname)
            push!(results.rt_plots, rt_plot)  # Store for combined PDF
            @debug_l1 "Generated normal RT plot"
        else
            # Create a diagnostic plot showing fallback/borrowed status
            fallback_plot = generate_fallback_rt_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(results.rt_plots, fallback_plot)  # Store for combined PDF
            end
            @debug_l1 "Generated fallback RT plot"
        end
        
        # Generate mass error plot (only store in memory, no individual files)
        @debug_l1 "Generating mass error plot for $parsed_fname: PPM errors = $(length(results.ppm_errs))"
        if length(results.ppm_errs) > 0
            mass_plot = generate_mass_error_plot(results, parsed_fname)
            push!(results.mass_plots, mass_plot)  # Store for combined PDF
            @debug_l1 "Generated normal mass error plot"
        else
            # Create a diagnostic plot showing fallback/borrowed status
            fallback_plot = generate_fallback_mass_error_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(results.mass_plots, fallback_plot)  # Store for combined PDF
            end
            @debug_l1 "Generated fallback mass error plot"
        end
        
        # Update models in search context
        setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
        setRtIrtMap!(search_context, getRtToIrtModel(results), ms_file_idx)
        
        # Clear plotting data to save memory
        resize!(results.rt, 0)
        resize!(results.irt, 0)
        resize!(results.ppm_errs, 0)
    catch e
        @user_warn "Failed to generate plots for file $ms_file_idx" exception=(e, catch_backtrace())
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
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
    # Combine individual plots into merged PDFs using save_multipage_pdf
    @debug_l1 "Writing combined QC plots..."
    
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
            @debug_l1 "Created combined RT alignment PDF with $(length(results.rt_plots)) plots"
            empty!(results.rt_plots)  # Clear to free memory
        else
            @debug_l1 "No RT alignment plots to combine"
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
            @debug_l1 "Created combined mass error PDF with $(length(results.mass_plots)) plots"
            empty!(results.mass_plots)  # Clear to free memory
        else
            @debug_l1 "No mass error plots to combine"
        end
        
        # Generate summary report
        # TODO: Implement generate_summary_report if detailed report needed
        # For now, the diagnostic summary below provides the key information
        
    catch e
        @user_warn "Failed to merge QC plots" exception=(e, catch_backtrace())
    end
    
    # Log diagnostic summary
    diagnostics = getDiagnostics(results)
    # Fixed: use values() to iterate over dictionary values
    file_statuses = values(diagnostics.file_statuses)
    converged_count = sum(s.converged for s in file_statuses)
    fallback_count = sum(s.used_fallback for s in file_statuses)
    total_files = length(file_statuses)
    
    @debug_l1 "Parameter Tuning Summary:"
    @debug_l1 "  - Total files: $total_files"
    @debug_l1 "  - Converged: $converged_count"
    @debug_l1 "  - Used fallback: $fallback_count"
    
    # Log any warnings
    for status in file_statuses
        if !isempty(status.warnings)
            @user_warn "File $(status.file_name) warnings:" warnings=status.warnings
        end
    end
end