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
    FragmentIndexSearch

Fragment index matching search. Runs fragment_index_search_only() per file,
writes scan→precursor mappings to Arrow, and builds a minimal precursor dict
for downstream searches.

This is the first real search step: it identifies which precursors have
fragment matches in each scan, without deconvolution or scoring.
"""
struct FragmentIndexSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for fragment index search.
"""
struct FragmentIndexSearchResults <: SearchResults
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32,mad_fwhm::Float32}}
    psms::Base.Ref{DataFrame}
    ms1_mass_err_model::Base.Ref{<:MassErrorModel}
    ms1_ppm_errs::Vector{Float32}
    ms1_mass_plots::Vector{Plots.Plot}
    qc_plots_folder_path::String
end

"""
Parameters for fragment index search.
Only needs isotope error bounds, minimum score, and spec order.
"""
struct FirstPassSearchParameters <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    spec_order::Set{Int64}

    function FirstPassSearchParameters(params::PioneerParameters)
        isotope_bounds = params.global_settings.isotope_settings.err_bounds_first_pass
        frag_idx_params = params.fragment_index_search
        new(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            UInt8(frag_idx_params.min_score),
            Set{Int64}([2])
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FragmentIndexSearch, params::Any) = FirstPassSearchParameters(params)
getMs1MassErrorModel(r::FragmentIndexSearchResults) = r.ms1_mass_err_model[]

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
    !isdir(ms1_mass_error_plots) && mkdir(ms1_mass_error_plots)
    return FragmentIndexSearchResults(
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
Process a single MS file: run fragment index search only, write matches to Arrow.
"""
function process_file!(
    results::FragmentIndexSearchResults,
    params::FirstPassSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)

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

        # Create empty PSMs DataFrame (needed by interface)
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

        reason = "FragmentIndexSearch failed: $e"
        markFileFailed!(search_context, ms_file_idx, reason)
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "Fragment index search failed for MS data file: $file_name. Error: $e."

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
No additional per-file processing needed for fragment index search.
"""
function process_search_results!(
    results::FragmentIndexSearchResults,
    params::FirstPassSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
)
    # Set default FWHM values (no PSMs to measure)
    insert!(results.fwhms, ms_file_idx, (
        median_fwhm = 0.5f0,
        mad_fwhm = 0.2f0))

    # Update MS1 mass error model with default
    setMs1MassErrorModel!(search_context, ms_file_idx, getMassErrorModel(search_context, ms_file_idx))
end

"""
No cleanup needed between files.
"""
function reset_results!(results::FragmentIndexSearchResults)
    empty!(results.psms[])
    resize!(results.ms1_ppm_errs, 0)
    return nothing
end

"""
Summarize results across all files.
Build precursor dict from fragment index match files.
"""
function summarize_results!(
    results::FragmentIndexSearchResults,
    params::FirstPassSearchParameters,
    search_context::SearchContext
)

    r(t) = round(t; digits=2)
    t_total_start = time()

    # Collect unique precursor IDs from all fragment index match files
    valid_indices = get_valid_file_indices(search_context)

    # Collect all unique precursor IDs across all fragment index match files
    t1_start = time()
    all_prec_ids = Set{UInt32}()
    for idx in valid_indices
        match_path = getFragmentIndexMatches(getMSData(search_context), idx)
        isempty(match_path) && continue
        !isfile(match_path) && continue
        tbl = Arrow.Table(match_path)
        union!(all_prec_ids, Set{UInt32}(tbl[:precursor_idx]))
    end
    t1 = time() - t1_start

    @user_info "FragmentIndexSearch: $(length(all_prec_ids)) unique precursors from fragment index"

    t_total = time() - t_total_start
    @info "FragmentIndexSearch summarize: arrow_read=$(r(t1))s, total=$(r(t_total))s"
end
