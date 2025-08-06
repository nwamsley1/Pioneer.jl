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
        
        # Create filtered/sampled data structure with initial scan count
        filtered_spectra = FilteredMassSpecData(
            spectra,
            max_scans = params.initial_scan_count,  # Default: 2500
            topn = params.topn_peaks,
            target_ms_order = UInt8(2)  # Only MS2 for presearch
        )
        
        # Count total MS2 scans for logging
        total_ms2_count = 0
        for i in 1:length(spectra)
            if getMsOrder(spectra, i) == UInt8(2)
                total_ms2_count += 1
            end
        end
        
        @info "Initial sampling: $(length(filtered_spectra)) MS2 scans from $total_ms2_count total for file $ms_file_idx"
        
        # Use bias estimate from cross-run learning or default
        initial_bias = initial_params.bias_estimate
        
        # Set initial mass error model
        results.mass_err_model[] = MassErrorModel(
            initial_bias,
            (initial_params.initial_tolerance, initial_params.initial_tolerance)
        )
        
        prev_mass_err = 0.0f0
        init_mass_tol = initial_params.initial_tolerance
        ppm_errs = nothing
        
        # Main convergence loop
        while n_attempts < 5
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))

            # Collect PSMs through iterations
            psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
            final_psm_count = size(psms, 1)
            @info "Iteration $(n_attempts + 1): Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans"
            
            # Check if we have enough PSMs to proceed
            if final_psm_count < min_psms_for_fitting
                @warn "Insufficient PSMs: $final_psm_count < $min_psms_for_fitting"
                push!(warnings, "Insufficient PSMs: $final_psm_count < $min_psms_for_fitting")
                
                # First retry: Expand to more scans
                if n_attempts == 0 && length(filtered_spectra) < params.expanded_scan_count
                    additional_scans = params.expanded_scan_count - length(filtered_spectra)
                    n_added = append!(filtered_spectra, max_additional_scans = additional_scans)
                    @info "Expanded sampling: added $n_added scans (now $(length(filtered_spectra)) total)"
                    push!(warnings, "Expanded scan sampling to $(length(filtered_spectra)) scans")
                    n_attempts += 1
                    continue
                end
                
                # Second retry: Increase mass tolerance
                if n_attempts == 1
                    current_tol = getLeftTol(results.mass_err_model[]) 
                    new_tol = min(current_tol * 1.5f0, 50.0f0)
                    results.mass_err_model[] = MassErrorModel(
                        getMassOffset(results.mass_err_model[]),
                        (new_tol, new_tol)
                    )
                    @info "Expanding mass tolerance from $current_tol to $new_tol ppm"
                    push!(warnings, "Expanded tolerance to $new_tol ppm due to low PSM count")
                    n_attempts += 1
                    continue
                end
                
                # Give up after two retries
                break
            end
            
            # Fit RT alignment model
            set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                fit_irt_model(params, psms))

            # Get fragments and fit mass error model
            fragments = get_matched_fragments(spectra, psms, results, search_context, params, ms_file_idx)
            mass_err_model, ppm_errs_new = fit_mass_err_model(params, fragments)
            ppm_errs = ppm_errs_new  # Update for boundary checking
            
            # Check boundary sampling adequacy
            boundary_result = check_boundary_sampling(ppm_errs_new, mass_err_model)
            if !boundary_result.adequate_sampling
                @warn boundary_result.diagnostic_message
                push!(warnings, boundary_result.diagnostic_message)
                # Expand tolerance based on boundary analysis
                expanded_model = expand_tolerance(mass_err_model, boundary_result.suggested_expansion_factor)
                results.mass_err_model[] = expanded_model
                init_mass_tol *= boundary_result.suggested_expansion_factor
                n_attempts += 1
                continue
            end
            
            # Check convergence conditions
            if abs(getMassOffset(mass_err_model)) > (getFragTolPpm(params)/4)
                prev_mass_err = getMassOffset(mass_err_model)
                @warn "Mass error offset is too large: $(getMassOffset(mass_err_model)). Retrying with updated estimate."
                push!(warnings, "Large mass offset: $(getMassOffset(mass_err_model)) ppm")
                results.mass_err_model[] = MassErrorModel(
                    getMassOffset(mass_err_model), 
                    (getFragTolPpm(params), getFragTolPpm(params))
                )
            elseif (getLeftTol(mass_err_model) + getRightTol(mass_err_model)) > 1.8*(init_mass_tol) 
                init_mass_tol *= 1.5
                results.mass_err_model[] = MassErrorModel(
                    0.0f0, 
                    (Float32(init_mass_tol), Float32(init_mass_tol))
                )
                @warn "Mass error model is too wide: $(getLeftTol(mass_err_model) + getRightTol(mass_err_model)). Retrying with wider initial mass tolerance: $init_mass_tol ppm"
                push!(warnings, "Wide tolerance: $(getLeftTol(mass_err_model) + getRightTol(mass_err_model)) ppm")
            else
                results.mass_err_model[] = MassErrorModel(
                    getMassOffset(mass_err_model) + prev_mass_err, 
                    (getLeftTol(mass_err_model), getRightTol(mass_err_model))
                )
                converged = true
                break
            end
            n_attempts += 1
        end
        
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
            @warn "Failed to converge mass error model after $n_attempts attempts for file $ms_file_idx. Using conservative defaults." 
            # Apply fallback parameters
            results.mass_err_model[] = MassErrorModel(
                0.0f0,  # No bias assumption
                (50.0f0, 50.0f0)  # ±50 ppm tolerance
            )
            results.rt_to_irt_model[] = IdentityModel()  # No RT conversion
            
            @warn "PARAMETER_TUNING_FALLBACK: File $parsed_fname used fallback values (±50 ppm, identity RT)"
        else
            # Normal case - append collected errors
            if ppm_errs !== nothing
                append!(results.ppm_errs, ppm_errs)
            end
        end
        
    catch e
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        @warn "Parameter tuning failed for file $parsed_fname. Using conservative defaults." exception=(e, catch_backtrace())
        
        # Apply fallback parameters
        results.mass_err_model[] = MassErrorModel(
            0.0f0,  # No bias assumption
            (50.0f0, 50.0f0)  # ±50 ppm tolerance
        )
        results.rt_to_irt_model[] = IdentityModel()  # No RT conversion
        
        # Record error in diagnostics
        status = ParameterTuningStatus(
            ms_file_idx,
            parsed_fname,
            false,  # not converged
            true,   # used fallback
            "Exception: $(typeof(e))",
            0,      # n_iterations
            0,      # psm_count
            0.0f0,  # mass_offset
            (50.0f0, 50.0f0),  # tolerance
            ["Error: $(string(e))"]
        )
        record_tuning_status!(results.diagnostics, status)
        
        @warn "PARAMETER_TUNING_ERROR: File $parsed_fname had error $(typeof(e)). Used fallback values."
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
    
    # Generate and save RT alignment plot
    rt_plot_path = joinpath(rt_alignment_folder, parsed_fname*".pdf")
    generate_rt_plot(results, rt_plot_path, parsed_fname)
    
    # Generate and save mass error plot
    mass_plot_path = joinpath(mass_error_folder, parsed_fname*".pdf")
    generate_mass_error_plot(results, parsed_fname, mass_plot_path)
    
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


