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
    global_pep_threshold::Float32
    # RT parameters
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    plot_rt_alignment::Bool
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
            Float32(score_params.global_pep_threshold),

            Int64(1000), # Default min_inference_points
            Float32(rt_params.min_probability),
            Float32(rt_params.min_probability),
            Float32(0.1), # Default max_irt_bin_size
            Float32(irt_mapping_params.max_prob_to_impute_irt),  # Default max_prob_to_impute
            Float32(irt_mapping_params.fwhm_nstd),   # Default fwhm_nstd
            Float32(irt_mapping_params.irt_nstd),   # Default irt_nstd
            Bool(hasproperty(irt_mapping_params, :plot_rt_alignment) ? irt_mapping_params.plot_rt_alignment : false),
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
    frag_match_folder = joinpath(getDataOutDir(search_context), "temp_data", "fragment_index_matches")
    !isdir(frag_match_folder) && mkdir(frag_match_folder)
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
Bypass mode: only runs fragment index search, skips probit scoring and MS1 mass error.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:FirstPassSearchParameters}

    try
        # Run fragment index search only (no getPSMS, no probit scoring)
        scan_to_prec_idx, precursors_passed = fragment_index_search_only(
            spectra, search_context, params, ms_file_idx
        )

        # Write fragment index matches to Arrow
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        output_path = joinpath(
            getDataOutDir(search_context), "temp_data", "fragment_index_matches",
            parsed_fname * ".arrow"
        )
        write_fragment_index_matches(scan_to_prec_idx, precursors_passed, output_path)
        setFragmentIndexMatches!(getMSData(search_context), ms_file_idx, output_path)

        # Create empty PSMs DataFrame (needed by process_search_results!)
        results.psms[] = DataFrame(
            ms_file_idx = UInt32[], scan_idx = UInt32[], precursor_idx = UInt32[],
            rt = Float32[], irt_predicted = Float32[], q_value = Float32[],
            score = Float32[], prob = Float32[], scan_count = UInt32[],
            fwhm = Union{Missing, Float32}[], PEP = Float16[]
        )

        # Use default MS1 mass error model (no PSMs to estimate from)
        results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)

    catch e
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end

        reason = "FirstPassSearch failed: $e"
        markFileFailed!(search_context, ms_file_idx, reason)
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "First pass search failed for MS data file: $file_name. Error: $e."

        try
            bt = catch_backtrace()
            @user_error "Full error details:\n" * sprint(showerror, e, bt)
        catch
        end

        results.psms[] = DataFrame(
            ms_file_idx = UInt32[], scan_idx = UInt32[], precursor_idx = UInt32[],
            rt = Float32[], irt_predicted = Float32[], q_value = Float32[],
            score = Float32[], prob = Float32[], scan_count = UInt32[],
            fwhm = Union{Missing, Float32}[], PEP = Float16[]
        )
        results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
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
    # Bypass mode: set default FWHM values (no PSMs to measure)
    insert!(results.fwhms, ms_file_idx, (
        median_fwhm = 0.5f0,
        mad_fwhm = 0.2f0))

    # Update MS1 mass error model with default
    setMs1MassErrorModel!(search_context, ms_file_idx, getMassErrorModel(search_context, ms_file_idx))
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
Bypass mode: build precursor dict from fragment index match files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}

    # Build precursor dict from all fragment index match files
    valid_indices = get_valid_file_indices(search_context)
    precursors = getPrecursors(getSpecLib(search_context))
    prec_irts = getIrt(precursors)
    prec_mzs = getMz(precursors)

    # Collect all unique precursor IDs across all fragment index match files
    all_prec_ids = Set{UInt32}()
    for idx in valid_indices
        match_path = getFragmentIndexMatches(getMSData(search_context), idx)
        isempty(match_path) && continue
        !isfile(match_path) && continue
        tbl = Arrow.Table(match_path)
        union!(all_prec_ids, Set{UInt32}(tbl[:precursor_idx]))
    end

    @user_info "Bypass FirstPassSearch: $(length(all_prec_ids)) unique precursors from fragment index"

    # Build minimal precursor dict with library iRT values
    precursor_dict = Dictionary{UInt32, @NamedTuple{
        best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32,
        best_irt::Float32, mean_irt::Union{Missing, Float32},
        var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32
    }}()

    for pid in all_prec_ids
        insert!(precursor_dict, pid, (
            best_prob = 1.0f0,
            best_ms_file_idx = UInt32(0),
            best_scan_idx = UInt32(0),
            best_irt = prec_irts[pid],
            mean_irt = missing,
            var_irt = missing,
            n = missing,
            mz = prec_mzs[pid]
        ))
        setPredIrt!(search_context, pid, prec_irts[pid])
    end

    setPrecursorDict!(search_context, precursor_dict)

    # Skip: map_retention_times!, create_rt_indices!, MS1 mass error plots
    # Second pass will use fragment index matches directly.
end

