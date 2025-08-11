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

"""
Main file processing method for parameter tuning search.
"""
function process_file!(
    results::ParameterTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    
    # Initialize diagnostics for this file
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    warnings = String[]
    final_psm_count = 0

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

    # Get parsed filename for diagnostics
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    try
        # Get initial parameters from cross-run learning
        initial_params = get_initial_parameters(results.parameter_history, params)
        
        # Initialize tracking variables
        warnings = String[]
        converged = false
        n_attempts = 0
        min_psms_for_fitting = 1000
        final_psm_count = 0
        
        # Count total MS2 scans for logging
        total_ms2_count = 0
        for i in 1:length(spectra)
            if getMsOrder(spectra, i) == UInt8(2)
                total_ms2_count += 1
            end
        end
        
        # Use bias estimate from cross-run learning or default
        initial_bias = initial_params.bias_estimate
        
        # Store initial parameters for resetting at each iteration
        initial_mass_offset = initial_bias
        initial_tolerance = initial_params.initial_tolerance
        max_tolerance = getMaxTolerancePpm(params)
        
        # Set initial mass error model
        results.mass_err_model[] = create_capped_mass_model(
            initial_mass_offset,
            initial_tolerance,
            initial_tolerance,
            max_tolerance
        )
        
        ppm_errs = nothing
        
        # Declare psms and final_psm_count outside loop for scope access
        psms = nothing
        final_psm_count = 0

        setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
        setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))
        n, N = 0, 5
        
        # Initialize filtered_spectra with zero scans
        filtered_spectra = FilteredMassSpecData(
            spectra,
            max_scans = 0,  # Start with zero scans
            topn = 200,     # Top 200 peaks per scan
            target_ms_order = UInt8(2)  # Only MS2 for presearch
        ) 
        @info "Initial mass error model getMassErrorModel for ms_file_idx $ms_file_idx: $(getMassErrorModel(search_context, ms_file_idx))"
        while n < N
            # Reset to initial parameters at start of each iteration (except first)
            if n > 0
                @info "Resetting to initial mass error parameters for iteration $(n+1)"
                results.mass_err_model[] = create_capped_mass_model(
                    initial_mass_offset,
                    initial_tolerance,
                    initial_tolerance,
                    max_tolerance
                )
                setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            end
            
            #1) First collect/sample more psms 
            # Calculate additional scans to add
            if n == 0
                # First iteration: add initial scan count
                additional_scans = params.initial_scan_count
            else
                # Subsequent iterations: expand up to expanded_scan_count limit
                current_scans = length(filtered_spectra)
                additional_scans = min(params.expanded_scan_count - current_scans, 2500)
            end
            
            @info "Iteration $(n+1): Adding $additional_scans scans (currently $(length(filtered_spectra)) scans)"
            append!(filtered_spectra; max_additional_scans = additional_scans)
            @info "Strategy 1: Collecting PSMs with current parameters"
            psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
            final_psm_count = size(psms, 1)
            @info "Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans"
            
            fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            if check_convergence(psms, mass_err_model, results.mass_err_model[])
                @info "Converged after scan expansion"
                results.mass_err_model[] = mass_err_model
                
                # Store ppm_errs for plotting
                if ppm_errs_new !== nothing && length(ppm_errs_new) > 0
                    resize!(results.ppm_errs, 0)  # Clear any old data
                    append!(results.ppm_errs, ppm_errs_new .+ getMassOffset(mass_err_model))
                    @info "Stored $(length(results.ppm_errs)) ppm errors for plotting"
                end
                
                set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                fit_irt_model(params, psms))
                @info "Stored $(length(results.rt)) RT points and $(length(results.irt)) iRT points for plotting"
                
                converged = true
                break
            end

            #2) Expand mass tolerance
            @info "Strategy 2: Expanding mass tolerance" 
            current_left_tol = getLeftTol(results.mass_err_model[])
            current_right_tol = getRightTol(results.mass_err_model[])
            new_left_tol = current_left_tol * 1.5f0
            new_right_tol = current_right_tol * 1.5f0
            
            results.mass_err_model[] = create_capped_mass_model(
                getMassOffset(results.mass_err_model[]),
                new_left_tol,
                new_right_tol,
                max_tolerance
            )
            
            # Log actual values after capping
            actual_left = getLeftTol(results.mass_err_model[])
            actual_right = getRightTol(results.mass_err_model[])
            if actual_left < new_left_tol || actual_right < new_right_tol
                @info "Expanded tolerance from ($(round(current_left_tol, digits=1)), $(round(current_right_tol, digits=1))) to ($(round(actual_left, digits=1)), $(round(actual_right, digits=1))) ppm (capped at $max_tolerance)"
            else
                @info "Expanded tolerance from ($(round(current_left_tol, digits=1)), $(round(current_right_tol, digits=1))) to ($(round(actual_left, digits=1)), $(round(actual_right, digits=1))) ppm"
            end
            
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
            final_psm_count = size(psms, 1)
            @info "Collected $final_psm_count PSMs with expanded tolerance"
            
            fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)

            #3) Adjust mass bias
            @info "Strategy 3: Adjusting mass bias" 
            new_bias = getMassOffset(mass_err_model)
            old_bias = getMassOffset(results.mass_err_model[])
            @info "Adjusting mass bias from $(round(old_bias, digits=2)) to $(round(new_bias, digits=2)) ppm"
            
            # Keep current tolerances but update bias
            current_left_tol = getLeftTol(results.mass_err_model[])
            current_right_tol = getRightTol(results.mass_err_model[])
            
            results.mass_err_model[] = create_capped_mass_model(
                new_bias,
                current_left_tol,
                current_right_tol,
                max_tolerance
            )
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
            final_psm_count = size(psms, 1)
            @info "Collected $final_psm_count PSMs with adjusted bias"
            
            fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            
            if check_convergence(psms, mass_err_model, results.mass_err_model[])
                @info "Converged after bias adjustment"
                results.mass_err_model[] = mass_err_model
                
                # Store ppm_errs for plotting
                if ppm_errs_new !== nothing && length(ppm_errs_new) > 0
                    resize!(results.ppm_errs, 0)  # Clear any old data
                    append!(results.ppm_errs, ppm_errs_new .+ getMassOffset(mass_err_model))
                    @info "Stored $(length(results.ppm_errs)) ppm errors for plotting"
                end
                
                set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                fit_irt_model(params, psms))
                @info "Stored $(length(results.rt)) RT points and $(length(results.irt)) iRT points for plotting"
                
                converged = true
                break
            end
            
            n += 1
            @info "Iteration $(n) complete. Moving to next iteration..."
        end

        # Main convergence loop
        #=
        while n_attempts < 5
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))

            # Collect PSMs through iterations
            psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
            final_psm_count = size(psms, 1)
            psm_improvement = final_psm_count - prev_psm_count
            psm_improvement_pct = prev_psm_count > 0 ? (psm_improvement / prev_psm_count * 100) : 0
            
            @info "Iteration $(n_attempts + 1): Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans ($(round(psm_improvement_pct, digits=1))% change)"
            
            # Check if we have enough PSMs to proceed with full fitting
            if final_psm_count < min_psms_for_fitting
                @warn "Insufficient PSMs: $final_psm_count < $min_psms_for_fitting"
                push!(warnings, "Iteration $(n_attempts + 1): $final_psm_count PSMs")
                
                # Strategy 1: First attempt with very few PSMs - try scan expansion
                if n_attempts == 0 && length(filtered_spectra) < params.expanded_scan_count
                    additional_scans = params.expanded_scan_count - length(filtered_spectra)
                    n_added = append!(filtered_spectra; max_additional_scans = additional_scans)
                    @info "Strategy: Expanding scan sampling - added $n_added scans (now $(length(filtered_spectra)) total)"
                    push!(warnings, "Expanded scan sampling to $(length(filtered_spectra)) scans")
                    prev_psm_count = final_psm_count
                    n_attempts += 1
                    continue
                end
                
                # Strategy 2: Have some PSMs (≥50) - try bias correction if we haven't just done it
                if final_psm_count >= min_psms_for_bias && !bias_corrected_last
                    @info "Strategy: Attempting mass bias correction with $final_psm_count PSMs"
                    
                    fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
                    
                    if length(fragments) > 0
                        temp_mass_err_model, _ = fit_mass_err_model(params, fragments)
                        new_bias = getMassOffset(temp_mass_err_model)
                        
                        # Keep current tolerance but update bias
                        current_tol = max(getLeftTol(results.mass_err_model[]), 
                                         getRightTol(results.mass_err_model[]))
                        results.mass_err_model[] = MassErrorModel(new_bias, (getLeftTol(results.mass_err_model[]), getRightTol(results.mass_err_model[])))
                        
                        @info "Corrected mass bias: $(round(new_bias, digits=2)) ppm (was $(round(getMassOffset(results.mass_err_model[]), digits=2)) ppm)"
                        push!(warnings, "Bias correction: $(round(new_bias, digits=2)) ppm")
                        bias_corrected_last = true
                    else
                        @warn "Could not calculate bias - no fragment matches"
                        bias_corrected_last = false
                    end
                    
                # Strategy 3: Significant improvement (>20%) after tolerance increase - try bias correction
                elseif psm_improvement_pct > 20 && final_psm_count >= min_psms_for_bias && !bias_corrected_last
                    @info "Strategy: PSMs improved by $(round(psm_improvement_pct, digits=1))% - attempting bias correction"
                    
                    psms_for_bias = first(psms, min(final_psm_count, 200))
                    fragments = get_matched_fragments(spectra, psms_for_bias, results, search_context, params, ms_file_idx)
                    
                    if length(fragments) > 0
                        temp_mass_err_model, _ = fit_mass_err_model(params, fragments)
                        new_bias = getMassOffset(temp_mass_err_model)
                        
                        current_tol = max(getLeftTol(results.mass_err_model[]), 
                                         getRightTol(results.mass_err_model[]))
                        results.mass_err_model[] = MassErrorModel(new_bias, (current_tol, current_tol))
                        
                        @info "Corrected mass bias after improvement: $(round(new_bias, digits=2)) ppm"
                        push!(warnings, "Bias correction after improvement: $(round(new_bias, digits=2)) ppm")
                        bias_corrected_last = true
                    else
                        bias_corrected_last = false
                    end
                    
                # Strategy 4: Too few PSMs or bias didn't help - expand tolerance
                else
                    current_tol = max(getLeftTol(results.mass_err_model[]), 
                                     getRightTol(results.mass_err_model[]))
                    new_tol = min(current_tol * 1.5f0, 50.0f0)
                    
                    if new_tol > current_tol
                        # Keep current bias if we have corrected it
                        current_bias = getMassOffset(results.mass_err_model[])
                        results.mass_err_model[] = MassErrorModel(current_bias, (new_tol, new_tol))
                        
                        @info "Strategy: Expanding mass tolerance from $(round(current_tol, digits=1)) to $(round(new_tol, digits=1)) ppm"
                        push!(warnings, "Expanded tolerance to $(round(new_tol, digits=1)) ppm")
                        bias_corrected_last = false
                    else
                        @warn "Already at maximum tolerance (50 ppm) - cannot expand further"
                        push!(warnings, "At maximum tolerance")
                        # No more strategies available
                        if n_attempts >= 2  # Give at least 3 attempts total
                            break
                        end
                    end
                end
                
                prev_psm_count = final_psm_count
                n_attempts += 1
                continue
            else
                break
                converged = true
            end
            n_attempts += 1
        end
        =#
            # Fit RT alignment model
            #set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
            #                    fit_irt_model(params, psms))

            # Get fragments and fit mass error model
            #fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            #mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            #ppm_errs = ppm_errs_new  # Update for later use
            
            # Check convergence conditions
            # Get the current tolerance being used for search
            #if abs(getMassOffset(mass_err_model))>max(abs(getLeftTol(mass_err_model)), abs(getRightTol(mass_err_model)))/4
            #    prev_mass_err = getMassOffset(mass_err_model)
            #    results.mass_err_model[] = MassErrorModel(getMassOffset(mass_err_model), (getLeftTol(mass_err_model), getRightTol(mass_err_model)))
            #else
            #    results.mass_err_model[] = MassErrorModel(getMassOffset(mass_err_model) + prev_mass_err, (getLeftTol(mass_err_model), getRightTol(mass_err_model)))
            #    converged = true
            #    break
            #end
            #n_attempts += 1
        #end
        
        # Log final data status before storing results
        @info "Final data status for file $ms_file_idx:"
        @info "  - RT points: $(length(results.rt))"
        @info "  - iRT points: $(length(results.irt))" 
        @info "  - PPM errors: $(length(results.ppm_errs))"
        @info "  - Converged: $converged"
        
        # Record tuning results
        tuning_result = TuningResults(
            getMassOffset(results.mass_err_model[]),
            (getLeftTol(results.mass_err_model[]), getRightTol(results.mass_err_model[])),
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
            getMassOffset(results.mass_err_model[]),
            (getLeftTol(results.mass_err_model[]), getRightTol(results.mass_err_model[])),
            warnings
        )
        record_tuning_status!(results.diagnostics, status)
        
        if !converged
            @warn "Failed to converge mass error model after $n_attempts attempts for file $ms_file_idx."
            
            # Try to borrow parameters from neighboring files or use defaults
            fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
                search_context, ms_file_idx, warnings
            )
            
            # Apply fallback/borrowed parameters
            results.mass_err_model[] = fallback_mass_err
            results.rt_to_irt_model[] = fallback_rt_model
            
            # CRITICAL: Store models in SearchContext for downstream methods
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)
            
            # Set default IRT error for fallback/borrowed parameters
            if borrowed_from !== nothing && haskey(getIrtErrors(search_context), borrowed_from)
                # Use borrowed IRT error
                getIrtErrors(search_context)[ms_file_idx] = getIrtErrors(search_context)[borrowed_from]
            else
                # Use conservative default (2 minutes)
                getIrtErrors(search_context)[ms_file_idx] = 2.0f0
            end
            
            if borrowed_from !== nothing
                borrowed_fname = getParsedFileName(search_context, borrowed_from)
                @warn "PARAMETER_TUNING_BORROWED: File $parsed_fname borrowed parameters from $borrowed_fname"
            else
                @warn "PARAMETER_TUNING_FALLBACK: File $parsed_fname used fallback values (±50 ppm, identity RT)"
            end
        else
            # Normal case - append collected errors
            if ppm_errs !== nothing
                append!(results.ppm_errs, ppm_errs)
            end
        end
        println("\n")
    catch e
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        @warn "Parameter tuning failed for file $parsed_fname with error." exception=(e, catch_backtrace())
        
        # Try to borrow parameters from neighboring files or use defaults
        error_warnings = String["Error: $(string(e))"]
        fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
            search_context, ms_file_idx, error_warnings
        )
        
        # Apply fallback/borrowed parameters
        results.mass_err_model[] = fallback_mass_err
        results.rt_to_irt_model[] = fallback_rt_model
        
        # CRITICAL: Store models in SearchContext for downstream methods
        setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
        setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)
        
        # Set default IRT error for exception case
        if borrowed_from !== nothing && haskey(getIrtErrors(search_context), borrowed_from)
            # Use borrowed IRT error
            getIrtErrors(search_context)[ms_file_idx] = getIrtErrors(search_context)[borrowed_from]
        else
            # Use conservative default (2 minutes)
            getIrtErrors(search_context)[ms_file_idx] = 2.0f0
        end
        
        # Record error in diagnostics
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
    
    # Update models in search context
    setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
    
    setRtIrtMap!(search_context, getRtToIrtModel(results), ms_file_idx)
    catch
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        nothing
    end
end

function reset_results!(ptsr::ParameterTuningSearchResults)
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    resize!(ptsr.ppm_errs, 0)
end


"""
Summarize results across all files and merge QC plots.
"""
function summarize_results!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:ParameterTuningSearchParameters}
    
    @info "Parameter Tuning Search Summary"
    @info "================================"
    
    # Generate diagnostics report
    diagnostics = results.diagnostics
    ms_data = getMSData(search_context)
    n_files = length(ms_data.file_paths)
    
    @info "Total files processed: $n_files"
    @info "Successfully tuned: $(diagnostics.n_successful)"
    @info "Used fallback parameters: $(diagnostics.n_fallback)"
    @info "Failed completely: $(diagnostics.n_failed)"
    
    # Generate parameter tuning report
    out_dir = getDataOutDir(search_context)
    report_path = joinpath(out_dir, "parameter_tuning_report.txt")
    generate_parameter_tuning_report(diagnostics, report_path)
    @info "Parameter tuning report saved to: $report_path"
    
    # Generate cross-run statistics report if applicable
    if results.parameter_history.global_stats.n_successful_files > 0
        cross_run_report_path = joinpath(out_dir, "cross_run_parameter_report.txt")
        generate_cross_run_report(results.parameter_history, cross_run_report_path)
        @info "Cross-run parameter report saved to: $cross_run_report_path"
    end
    
    # Warn if many files used fallback
    if diagnostics.n_fallback > 0
        @warn "$(diagnostics.n_fallback) file(s) used fallback parameters. Check parameter_tuning_report.txt for details."
    end
    
    @info "Merging QC plots..."
    
    # Merge RT alignment plots
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    output_path = joinpath(rt_alignment_folder, "rt_alignment_plots.pdf")
    try
        if isfile(output_path)
            safeRm(output_path, nothing)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    rt_plots = [joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) 
    if endswith(x, ".pdf")]
    
    if !isempty(rt_plots)
        merge_pdfs(rt_plots, 
                    output_path, 
                  cleanup=true)
    end
    
    # Merge mass error plots
    mass_error_folder = getMassErrPlotFolder(search_context)
    output_path = joinpath(mass_error_folder, "mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    mass_plots = [joinpath(mass_error_folder, x) for x in readdir(mass_error_folder) 
                    if endswith(x, ".pdf")]

    if !isempty(mass_plots)
        merge_pdfs(mass_plots, 
                  output_path, 
                  cleanup=true)
    end

    @info "QC plot merging complete"
end


