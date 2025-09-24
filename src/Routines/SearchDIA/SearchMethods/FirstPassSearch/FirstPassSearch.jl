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
    FirstPassSearch

Initial search method to identify PSMs and establish retention time calibration.

This search:
1. Performs initial PSM identification with learned scoring
2. Calculates retention time indices and FWHM statistics
3. Maps between library and empirical retention times
4. Generates RT calibration curves for subsequent searches

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (1, 0),
    :first_search_params => Dict(
        "n_train_rounds_probit" => 10,
        "max_iter_probit" => 100,
        "max_q_value_probit_rescore" => 0.01,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "max_precursors_passing" => 5000,
        "abreviate_precursor_calc" => false
    ),
    :summarize_first_search_params => Dict(
        "min_inference_points" => 1000,
        "max_q_val_for_irt" => 0.01,
        "max_precursors" => 10000,
        "max_irt_bin_size" => 0.1,
        "max_prob_to_impute" => 0.99
    ),
    :irt_mapping_params => Dict(
        "min_prob" => 0.9
    )
)

# Execute search
results = execute_search(FirstPassSearch(), search_context, params)
```
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#


"""
Results container for first pass search.
Holds FWHM statistics, PSM file paths, and model updates.
"""
struct FirstPassSearchResults <: SearchResults
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32,mad_fwhm::Float32}}
    psms::Base.Ref{DataFrame}
    ms1_mass_err_model::Base.Ref{<:MassErrorModel}
    ms1_ppm_errs::Vector{Float32}
    ms1_mass_plots::Vector{Plots.Plot}
    qc_plots_folder_path::String
end

"""
Parameters for first pass search.
Configures PSM identification, scoring, and RT calibration.
"""
struct FirstPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    frag_tol_ppm::Float32
    ms1_tol_ppm::Float32
    frag_err_quantile::Float32
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    match_between_runs::Bool
    relative_improvement_threshold::Float32
    
    # Scoring parameters
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
    max_PEP::Float32
    
    # RT parameters
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    prec_estimation::P

    function FirstPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        first_params = params.first_search
        frag_params = first_params.fragment_settings
        score_params = first_params.scoring_settings
        rt_params = params.rt_alignment
        irt_mapping_params = first_params.irt_mapping
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        # Determine precursor estimation strategy
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            0.0f0,  # No transmission threshold for first pass
            0.0f0,  # No fragment tolerance for first pass
            Float32(params.parameter_tuning.iteration_settings.ms1_tol_ppm),  # MS1 tolerance from config
            Float32(params.parameter_tuning.search_settings.frag_err_quantile),
            # Handle min_score as either single value or array (use first value if array)
            begin
                min_score_raw = frag_params.min_score
                if min_score_raw isa Vector
                    UInt8(first(min_score_raw))
                else
                    UInt8(min_score_raw)
                end
            end,
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1), # max_best_rank
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            global_params.match_between_runs,
            Float32(frag_params.relative_improvement_threshold),
            
            Int64(score_params.n_train_rounds),
            Int64(score_params.max_iterations),
            Float32(score_params.max_q_value_probit_rescore),
            Float32(score_params.max_PEP),
            
            Int64(1000), # Default min_inference_points
            Float32(rt_params.min_probability),
            Float32(rt_params.min_probability),
            Float32(0.1), # Default max_irt_bin_size
            Float32(irt_mapping_params.max_prob_to_impute_irt),  # Default max_prob_to_impute
            Float32(irt_mapping_params.fwhm_nstd),   # Default fwhm_nstd
            Float32(irt_mapping_params.irt_nstd),   # Default irt_nstd
            prec_estimation
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = FirstPassSearchParameters(params)
getMs1MassErrorModel(ptsr::FirstPassSearchResults) = ptsr.ms1_mass_err_model[]
getMs1TolPpm(params::FirstPassSearchParameters) = params.ms1_tol_ppm

function init_search_results(
    ::FirstPassSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
    !isdir(temp_folder) && mkdir(temp_folder)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    ms1_mass_error_plots = joinpath(qc_dir, "ms1_mass_error_plots")
    !isdir(ms1_mass_error_plots ) && mkdir(ms1_mass_error_plots )
    return FirstPassSearchResults(
        Dictionary{Int64, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}(),
        Base.Ref{DataFrame}(),
        Base.Ref{MassErrorModel}(),
        Vector{Float32}(),
        Plots.Plot[],
        qc_dir
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the first pass search.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:FirstPassSearchParameters}

    """
    Perform library search with current parameters.
    """
    function perform_library_search(
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        
        return library_search(spectra, search_context, params, ms_file_idx)
    end

    """
    Process PSMs from library search.
    """
    function process_psms!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)

        """
        Select best PSMs based on criteria.
        """
        function select_best_psms!(
            psms::DataFrame,
            precursor_mzs::AbstractVector,
            params::FirstPassSearchParameters,
            search_context::SearchContext
        )
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
            get_best_psms!(
                psms,
                precursor_mzs,
                max_PEP=params.max_PEP,
                fdr_scale_factor=fdr_scale_factor
            )
        end

        rt_model = getRtIrtModel(search_context, ms_file_idx)
        # Add columns
        add_psm_columns!(psms, spectra, search_context, rt_model, ms_file_idx)
        
        # Score PSMs
        score_psms!(psms, params, search_context)
        # Get best PSMs
        select_best_psms!(
            psms,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            params,
            search_context
        )
        return psms
    end

    """
    Add necessary columns to PSM DataFrame.
    """
    function add_psm_columns!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        rt_model::RtConversionModel,
        ms_file_idx::Int64)
        add_main_search_columns!(
            psms,
            getModel(rt_model),
            getStructuralMods(getPrecursors(getSpecLib(search_context))),
            getMissedCleavages(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getIrt(getPrecursors(getSpecLib(search_context))),
            getCharge(getPrecursors(getSpecLib(search_context))),
            getRetentionTimes(spectra),
            getTICs(spectra),
            getMzArrays(spectra)
        )
        # Calculate RT values
        psms[!, :irt_observed] = rt_model.(psms[!, :rt])
        psms[!, :irt_error] = Float16.(abs.(psms[!, :irt_observed] .- psms[!, :irt_predicted]))
        psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)
        psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    end

    """
    Score PSMs using probit model.
    """
    function score_psms!(
        psms::DataFrame,
        params::FirstPassSearchParameters,
        search_context::SearchContext)
        column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe, :percent_theoretical_ignored,
            :charge2, :poisson, :irt_error, 
            :missed_cleavage, 
            :Mox,
            #:charge, Only works with charge 2 if at least 3 charge states presence. otherwise singular error
            #:b_count, might be good for non-tryptic enzymes
            :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
        ]

        # Avoid singular error if no peaks were ignored
        if maximum(psms.percent_theoretical_ignored) == 0
            deleteat!(column_names, findfirst(==(:percent_theoretical_ignored), column_names))
        end


        # Select scoring columns
        select!(psms, vcat(column_names, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target]))
        # Score PSMs
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        try
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        catch
            column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe,
            :charge2, :poisson, :irt_error, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
            ]
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        end
        # Process scores
       
        select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target])
        get_probs!(psms, psms[!,:score])
        sort!(psms, :rt)
    end

    try
        # Get models and update fragment lookup table
        psms = perform_library_search(spectra, search_context, params, ms_file_idx)
        results.psms[] = process_psms!(psms, spectra, search_context, params, ms_file_idx)

        temp_psms = results.psms[] 
        temp_psms = temp_psms[temp_psms[!,:q_value].<=0.001,:]
        most_intense = sortperm(temp_psms[!,:log2_summed_intensity], rev = true)
        ms1_errs = vcat(
            mass_error_search(
                spectra,
                temp_psms[most_intense[1:(min(3000, length(most_intense)))],:scan_idx],
                temp_psms[most_intense[1:(min(3000, length(most_intense)))],:precursor_idx],
                UInt32(ms_file_idx),
                getSpecLib(search_context),
                getSearchData(search_context),
                MassErrorModel(
                0.0f0,
                (getMs1TolPpm(params), getMs1TolPpm(params))  # Use MS1 tolerance from JSON config
                ),
                params,
                MS1CHROM()
            )...
        )
        if length(ms1_errs) > 1
            mad_dev = mad(ms1_errs; normalize=true)
            med_errs = median(ms1_errs)
            low_bound, high_bound = med_errs - mad_dev*7, med_errs + mad_dev*7
            filter!(x->(low_bound<x)&(high_bound>x), ms1_errs)
            ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
            results.ms1_mass_err_model[] = ms1_mass_err_model
            append!(results.ms1_ppm_errs, ms1_ppm_errs)
            select!(results.psms[] , Not(:log2_summed_intensity))
        else
            #ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
            #Default to MS2 pattern 
            results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
            append!(results.ms1_ppm_errs, Float32[])
            select!(results.psms[] , Not(:log2_summed_intensity))
        end
    catch e
        # Get file name for debugging
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end

        reason = "FirstPassSearch failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        @user_warn "First pass search failed for MS data file: $file_name. Error type: $(typeof(e)). Creating empty results to continue pipeline."
        
        # Create an empty but properly structured DataFrame to avoid downstream errors
        empty_psms = DataFrame(
            ms_file_idx = UInt32[],
            scan_idx = UInt32[], 
            precursor_idx = UInt32[],
            rt = Float32[],
            irt_predicted = Float32[],
            q_value = Float32[],
            score = Float32[], 
            prob = Float32[],
            scan_count = UInt32[],
            fwhm = Float32[]  # Add this to prevent missing column error
        )
        results.psms[] = empty_psms
        
        # Set default mass error model
        results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
        #rethrow(e)
    end

    return results
end

"""
Initial file processing complete, no additional processing needed.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:FirstPassSearchParameters}
    psms = results.psms[]
    fwhms = skipmissing(psms[!, :fwhm])
    fwhm_points = count(!ismissing, fwhms)
    if fwhm_points >= 1#params.min_inference_points
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = median(fwhms),
            mad_fwhm = mad(fwhms, normalize=true)))
    else
        @user_warn "Insuficient fwhm_points to estimate for $ms_file_idx"
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = 0.2f0,
            mad_fwhm = 0.2f0))
    end
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    temp_path = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms", parsed_fname * ".arrow")
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
            :irt_predicted, :q_value, :score, :prob, :scan_count])
    )
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, temp_path)

    #####
    #MS1 mass tolerance 
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    # Generate mass error plot
    push!(results.ms1_mass_plots, generate_ms1_mass_error_plot(results, parsed_fname))
    # Update models in search context
    setMs1MassErrorModel!(search_context, ms_file_idx, getMs1MassErrorModel(results))
end

"""
No cleanup needed between files.
"""
function reset_results!(results::FirstPassSearchResults)
    empty!(results.psms[])
    resize!(results.ms1_ppm_errs, 0)
    return nothing
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}
    
    """
    Process precursors and calculate iRT errors.
    """
    function get_best_precursors_accross_runs!(
        search_context::SearchContext,
        results::FirstPassSearchResults,
        params::FirstPassSearchParameters
    )

        # Filter out failed files
        valid_indices = get_valid_file_indices(search_context)
        all_psms_paths = getFirstPassPsms(getMSData(search_context))
        valid_psms_paths = [all_psms_paths[i] for i in valid_indices]
        
        # Create RT-IRT map for valid files only
        all_rt_irt = getRtIrtModel(search_context)
        valid_rt_irt = Dict{Int64, RtConversionModel}(i => all_rt_irt[i] for i in valid_indices if haskey(all_rt_irt, i))
        
        if isempty(valid_psms_paths)
            @user_warn "No valid files for cross-run precursor analysis"
            return Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}}()
        end
        
        # Get best precursors from valid files only
        return get_best_precursors_accross_runs(
            valid_psms_paths,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            valid_rt_irt,
            max_q_val=params.max_q_val_for_irt
        )
    end
    # Map retention times
    map_retention_times!(search_context, results, params)
    # Process precursors
    precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)

    if false==true#params.match_between_runs==true
        #######
        #Each target has a corresponding decoy and vice versa
        #Add the complement targets/decoys to the precursor dict 
        #if the `sibling_peptide_scores` parameter is set to true
        #In the target/decoy scoring (see SearchMethods/ScoringSearch)
        #the maximum score for each target/decoy pair is shared accross runs
        #in an iterative training scheme. 
        precursors = getPrecursors(getSpecLib(search_context))
        i = 1
        for (pid, val) in pairs(precursor_dict)
            i += 1
            setPredIrt!(search_context, pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
            partner_pid = getPartnerPrecursorIdx(precursors)[pid]
            if ismissing(partner_pid)
                continue
            end

            # If the partner needs to be added, then give it the irt of the currently identified precursor
            # Otherwise if the partner was ID'ed, it should keep its original predicted iRT
            if !haskey(precursor_dict, partner_pid)
                insert!(precursor_dict, partner_pid, val)
                setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
            else
                setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[partner_pid])
            end
            
        end
    else
        for (pid, val) in pairs(precursor_dict)
            setPredIrt!(search_context, pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
        end
    end

    setPrecursorDict!(search_context, precursor_dict)
    # Calculate RT indices
    create_rt_indices!(search_context, results, precursor_dict, params)
    
    # Merge mass error plots
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    output_path = joinpath(ms1_mass_error_folder, "ms1_mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @user_warn "Could not clear existing file: $e"
    end

    if !isempty(results.ms1_mass_plots)
        save_multipage_pdf(results.ms1_mass_plots, output_path)
        empty!(results.ms1_mass_plots)
    end
end

