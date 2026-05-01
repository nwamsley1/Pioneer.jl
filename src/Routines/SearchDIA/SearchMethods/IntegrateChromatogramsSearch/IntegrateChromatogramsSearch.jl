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
    IntegrateChromatogramSearch

Search method for analyzing chromatograms to get quantitative information.

This search:
1. Uses precursor and trace information from previous searches
2. Builds chromatograms for each precursor
3. Integrates areas for quantification
4. Incorporates isotope pattern information
"""
struct IntegrateChromatogramSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for chromatogram integration search.
"""
struct IntegrateChromatogramSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}  # Chromatogram data per file
end

"""
Parameters for chromatogram integration search.
"""
struct IntegrateChromatogramSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}

    # Chromatogram parameters
    wh_smoothing_strength::Float32

    # Deconvolution parameters (MS2)
    lambda::Float32
    reg_type::RegularizationType
    deconvolution_solver::DeconvolutionSolver
    max_iter_outer::Int64
    max_diff::Float32

    # Analysis strategies
    isotope_tracetype::I
    prec_estimation::P

    function IntegrateChromatogramSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        search_params = params.search

        # Hardcoded — formerly read from global.isotope_settings (removed from
        # the public schema; values shipped at defaults in every config).
        isotope_trace_type = SeperateTraces()
        min_fraction_transmitted = 0.25f0
        prec_estimation = PartialPrecCapture()

        # max_rank hardcoded to typemax(UInt8); the library already enforces
        # a per-precursor cap at build time.
        max_frag_rank = UInt8(255)

        # n_isotopes lives at search.n_isotopes; falls back to the legacy
        # nested location for old configs (see _resolve_n_isotopes).
        n_isotopes_val = _resolve_n_isotopes(search_params)

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(1), UInt8(0)),  # isotope err_bounds
            Float32(min_fraction_transmitted),
            Int64(n_isotopes_val),
            max_frag_rank,
            Set{Int64}([2]),

            Float32(1e-6),    # wh_smoothing_strength

            Float32(0.0),     # lambda (no regularization)
            NoNorm(),         # reg_type
            PoissonMMSolver(),  # hardcoded: PMM is the production solver
            DECONV_MAX_ITER,          # max_iter_outer
            DECONV_CONVERGENCE_TOL,   # max_diff

            isotope_trace_type,
            prec_estimation
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::IntegrateChromatogramSearch, params::Any) = IntegrateChromatogramSearchParameters(params)

function init_search_results(::IntegrateChromatogramSearchParameters, search_context::SearchContext)
    return IntegrateChromatogramSearchResults(
        Ref(DataFrame())
    )
end

"""
Process a single file for chromatogram integration.
"""
function process_file!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData) where {P<:IntegrateChromatogramSearchParameters}

    # Check if required files exist (e.g. upstream step skipped this file)
    rt_index_path = getRtIndex(getMSData(search_context), ms_file_idx)
    passing_psms_path = getPassingPsms(getMSData(search_context), ms_file_idx)

    if isempty(rt_index_path) || isempty(passing_psms_path)
        file_name = getFileIdToName(getMSData(search_context), ms_file_idx)
        @debug_l2 "Skipping IntegrateChromatogramSearch for file $file_name - missing required files from previous steps"
        return results
    end

    # Build retention time index for efficient precursor lookup
    rt_index = buildRtIndex(
        DataFrame(Arrow.Table(rt_index_path)),
        bin_rt_size = 0.1)

    # Load PSMs that passed previous filtering steps. Decoys are intentionally
    # kept here — ProteinInferenceSearch and ProteinScoringSearch need them
    # for protein-level FDR / PEP calibration. Final decoy suppression for
    # output happens later in MaxLFQSearch when output.write_decoys=false.
    passing_psms = DataFrame(Tables.columntable(Arrow.Table(passing_psms_path)))

    # If there are no PSMs to integrate (e.g. sparse / empty file), skip
    # chromatogram extraction entirely. Downstream steps treat an empty
    # results table as a no-op file.
    if nrow(passing_psms) == 0
        results.psms[] = passing_psms
        return results
    end

    # Initialize columns to store integration results
    # peak_area: Integrated area of chromatographic peak
    # new_best_scan: Updated apex scan after refinement
    passing_psms[!, :peak_area] = zeros(Float32, nrow(passing_psms))
    passing_psms[!, :new_best_scan] = zeros(UInt32, nrow(passing_psms))
    passing_psms[!, :points_integrated] = zeros(UInt32, nrow(passing_psms))
    # Extract chromatograms for all passing PSMs
    chromatograms = extract_chromatograms(
        spectra,
        passing_psms,
        rt_index,
        search_context,
        params,
        ms_file_idx,
        MS2CHROM(),
    )
    # Save unsorted chromatograms for sorting benchmarks (first file only)
    if ms_file_idx == 1
        out_dir = getDataOutDir(search_context)
        bench_df = copy(chromatograms)
        bench_df[!, :rt_milliminutes] = round.(UInt32, bench_df.rt .* 1000)
        Arrow.write(joinpath(out_dir, "unsorted_chroms_ms2.arrow"), bench_df)
    end
    # MS1 chromatogram extraction is currently unwired; the MS1
    # build_chromatograms body is block-commented in utils.jl pending a
    # fused port. The ms1_quant knob has been removed from the public
    # config schema.
    #Arrow.write(joinpath(out_dir, "test_chroms_ms1.arrow"), ms1_chromatograms)
    #jldsave("/Users/nathanwamsley/Desktop/test_chroms_ms1.jld2"; ms1_chromatograms)
    if seperateTraces(params.isotope_tracetype)
        get_isotopes_captured!(
            chromatograms,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            chromatograms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )
    end

    if seperateTraces(params.isotope_tracetype)
        fast_df_sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])
    else
        fast_df_sort!(chromatograms, [:precursor_idx, :rt])
    end

    # CombineTraces: compute fraction_transmitted AFTER sort (order-independent per-row calc)
    if !seperateTraces(params.isotope_tracetype)
        get_isotopes_captured!(
            chromatograms,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            chromatograms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra);
            compute_isotope_set=false
        )
    end

    # Integrate chromatographic peaks for each precursor (skip if no chromatograms extracted)
    if nrow(chromatograms) > 0
        integrate_precursors(
            chromatograms,
            params.min_fraction_transmitted,
            passing_psms[!, :precursor_idx],
            passing_psms[!, :scan_idx],
            passing_psms[!, :peak_area],
            passing_psms[!, :new_best_scan],
            passing_psms[!, :points_integrated],
            λ = params.wh_smoothing_strength
        )
    end
    # MS1 integration disabled — see extract_chromatograms call above.
    # Clear chromatograms to free memory
    chromatograms = nothing

    # Store processed PSMs in results
    results.psms[] = passing_psms

    return results
end

function process_search_results!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:IntegrateChromatogramSearchParameters}

    passing_psms = results.psms[]

    # Skip processing if no PSMs (e.g. file had no rt_index or passing PSMs)
    if nrow(passing_psms) == 0 || ncol(passing_psms) == 0
        @debug_l2 "No PSMs to process for file $ms_file_idx in IntegrateChromatogramSearch results"
        return nothing
    end

    parsed_fname = getFileIdToName(getMSData(search_context), ms_file_idx)
    # Process final PSMs
    process_final_psms!(
        passing_psms,
        search_context,
        parsed_fname,
        ms_file_idx
    )
    # Save results
    writeArrow(getPassingPsms(getMSData(search_context))[ms_file_idx], passing_psms)
    return nothing
end

function reset_results!(results::IntegrateChromatogramSearchResults)
    results.psms[] = DataFrame()
    GC.gc()
    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    ::P,
    ::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    return nothing
end
