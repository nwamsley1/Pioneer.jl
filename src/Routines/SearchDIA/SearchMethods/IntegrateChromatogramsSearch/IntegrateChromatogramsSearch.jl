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
1. Loads the best PSMs selected by `MainSearch`
2. Builds MS2 chromatograms around each precursor's refined RT
3. Annotates chromatogram rows with precursor transmission for weighted smoothing
4. Integrates each peak with WH smoothing, second-derivative bounds, baseline
   subtraction, and trapezoidal area calculation

When separate isotope traces are requested, chromatogram rows also get an
`:isotopes_captured` grouping label so integration can choose one trace per
precursor.
"""
struct IntegrateChromatogramSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for chromatogram integration search.
"""
struct IntegrateChromatogramSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}  # PSM rows for one file after integration
end

function _resolve_chromatogram_trace_type(
    params::PioneerParameters,
    min_fraction_transmitted::Float32,
)
    trace_mode = params.optimization.chromatogram_integration.trace_mode

    if trace_mode == "combined"
        return CombineTraces(min_fraction_transmitted)
    elseif trace_mode == "separate"
        return SeperateTraces()
    end

    throw(ArgumentError(
        "Invalid optimization.chromatogram_integration.trace_mode=$(repr(trace_mode)); " *
        "expected \"combined\" or \"separate\"",
    ))
end

function default_chromatogram_integration_huber_solver()
    return HuberSolver(
        300.0f0,      # delta (hardcoded until tuning is reintroduced)
        0.0f0,        # lambda (no regularization)
        Int64(50),    # max_iter_newton
        Int64(100),   # max_iter_bisection
        10.0f0,       # accuracy_newton
        10.0f0,       # accuracy_bisection
        NoNorm(),
    )
end

function _resolve_chromatogram_deconvolution_solver(params::PioneerParameters)
    solver_name = get(
        params.optimization.chromatogram_integration,
        :deconvolution_solver,
        "huber",
    )

    if solver_name == "huber"
        return default_chromatogram_integration_huber_solver()
    elseif solver_name == "pmm"
        return PoissonMMSolver()
    end

    throw(ArgumentError(
        "Invalid optimization.chromatogram_integration.deconvolution_solver=$(repr(solver_name)); " *
        "expected \"huber\" or \"pmm\".",
    ))
end

with_chromatogram_huber_delta(solver::HuberSolver, delta::Float32) = with_huber_delta(solver, delta)
with_chromatogram_huber_delta(solver::DeconvolutionSolver, ::Float32) = solver

function calibrated_chromatogram_deconvolution_solver(
    search_context::SearchContext,
    solver::DeconvolutionSolver,
)
    return with_chromatogram_huber_delta(solver, getHuberDelta(search_context))
end

function chromatogram_integration_solver_label(solver::HuberSolver)
    return "HuberSolver(delta=$(Float64(solver.delta)))"
end

function chromatogram_integration_solver_label(::PoissonMMSolver)
    return "PoissonMMSolver"
end

function chromatogram_integration_solver_label(solver::DeconvolutionSolver)
    return string(nameof(typeof(solver)))
end

"""
Parameters for chromatogram integration search.

This stage configures trace mode and the deconvolution solver. Peak bounds are
selected inside `integrate_chrom`.
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
    match_between_runs::Bool
    q_value_threshold::Float32
    pep_bin_size::Int64

    function IntegrateChromatogramSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        search_params = params.search

        min_fraction_transmitted = 0.25f0
        isotope_trace_type = _resolve_chromatogram_trace_type(
            params,
            Float32(min_fraction_transmitted),
        )
        prec_estimation = PartialPrecCapture()

        # max_rank hardcoded to typemax(UInt8); the library already enforces
        # a per-precursor cap at build time.
        max_frag_rank = UInt8(255)

        # n_isotopes lives at search.n_isotopes; falls back to the legacy
        # nested location for old configs (see _resolve_n_isotopes).
        n_isotopes_val = _resolve_n_isotopes(search_params)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings
        mbr = hasproperty(global_params, :match_between_runs) ?
              Bool(global_params.match_between_runs) : true

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(1), UInt8(0)),  # isotope err_bounds
            Float32(min_fraction_transmitted),
            Int64(n_isotopes_val),
            max_frag_rank,
            Set{Int64}([2]),

            Float32(1e-6),    # wh_smoothing_strength

            Float32(0.0),     # lambda (no regularization)
            NoNorm(),         # reg_type
            _resolve_chromatogram_deconvolution_solver(params),
            DECONV_MAX_ITER,          # max_iter_outer
            DECONV_CONVERGENCE_TOL,   # max_diff

            isotope_trace_type,
            prec_estimation,
            mbr,
            _resolve_q_value_threshold(global_params),
            Int64(ml_params.pep_bin_size),
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::IntegrateChromatogramSearch, params::Any) = IntegrateChromatogramSearchParameters(params)

function write_intermediate_chromatogram_debug_plots(
    ::IntegrateChromatogramSearchParameters,
)
    return DEBUG_CONSOLE_LEVEL[] >= 1
end

function init_search_results(::IntegrateChromatogramSearchParameters, search_context::SearchContext)
    return IntegrateChromatogramSearchResults(
        Ref(DataFrame())
    )
end

function _filter_mbr_sidecars_by_rows!(main_path::String, selected_rows::Vector{Int})
    for suffix in (PASS1_SIDECAR_SUFFIX, MBR_SIDECAR_SUFFIX)
        side_path = main_path * suffix
        isfile(side_path) || continue
        side = DataFrame(Tables.columntable(Arrow.Table(side_path)))
        isempty(selected_rows) ? writeArrow(side_path, side[Int[], :]) :
            writeArrow(side_path, side[selected_rows, :])
    end
    return nothing
end

"""
Process a single file for chromatogram integration.

The integration input is the passing PSM table from `MainSearch`; no main-search
intervals are passed to the bound picker. Chromatograms are extracted,
transmission-annotated, sorted for the configured trace mode, then integrated.
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

    if !seperateTraces(params.isotope_tracetype)
        selected_rows = select_combined_trace_seed_rows_by_score(passing_psms)
        if params.match_between_runs
            _filter_mbr_sidecars_by_rows!(passing_psms_path, selected_rows)
        end
        passing_psms = passing_psms[selected_rows, :]
    end

    # Initialize columns written by chromatogram integration.
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
    # MS1 chromatogram extraction is currently unwired; the MS1
    # build_chromatograms body is block-commented in utils.jl pending a
    # fused port. The ms1_quant knob has been removed from the public
    # config schema.
    #Arrow.write(joinpath(out_dir, "test_chroms_ms1.arrow"), ms1_chromatograms)
    #jldsave("/Users/nathanwamsley/Desktop/test_chroms_ms1.jld2"; ms1_chromatograms)
    if nrow(chromatograms) > 0
        # WH smoothing uses precursor transmission as both a correction factor
        # and an observation weight. Separate-trace mode also uses isotope
        # labels for chromatogram grouping.
        get_isotopes_captured!(
            chromatograms,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            chromatograms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra),
            compute_isotope_set = compute_chromatogram_isotope_sets(params.isotope_tracetype),
        )
    end
    sort_chromatograms_for_integration!(chromatograms, params.isotope_tracetype)

    # Integrate chromatographic peaks for each precursor (skip if no chromatograms extracted)
    if nrow(chromatograms) > 0
        psm_isotopes_captured = nothing
        if seperateTraces(params.isotope_tracetype)
            # In separate-trace mode, choose the transmitted isotope trace used
            # for quantification and mirror that choice onto the PSM/output
            # columns. Combined mode keeps the upstream PSM isotope metadata.
            selected_quant_trace = select_quant_trace_by_transmission(chromatograms)
            apply_quant_trace_selection!(
                passing_psms,
                selected_quant_trace,
            )
            psm_isotopes_captured = passing_psms[!, :isotopes_captured]
        end
        integrate_precursors(
            chromatograms,
            params.isotope_tracetype,
            params.min_fraction_transmitted,
            passing_psms[!, :precursor_idx],
            passing_psms[!, :scan_idx],
            passing_psms[!, :peak_area],
            passing_psms[!, :new_best_scan],
            passing_psms[!, :points_integrated],
            isotopes_captured = psm_isotopes_captured,
            λ = params.wh_smoothing_strength,
        )
        add_mbr_integrated_spectra_to_psms!(
            passing_psms,
            chromatograms,
            getRtIrtModel(search_context, ms_file_idx),
            bitvec_rank_table = getBitVecExcessRanks(search_context, ms_file_idx),
        )
        if params.match_between_runs && isfile(passing_psms_path * PASS1_SIDECAR_SUFFIX)
            selected_rows = Int[
                row for row in axes(passing_psms, 1)
                if isfinite(Float32(passing_psms[row, :peak_area])) &&
                   Float32(passing_psms[row, :peak_area]) > 0.0f0 &&
                   isfinite(Float32(passing_psms[row, MBR_INTEGRATED_APEX_IRT_COLUMN]))
            ]
            if length(selected_rows) != nrow(passing_psms)
                _filter_mbr_sidecars_by_rows!(passing_psms_path, selected_rows)
                passing_psms = passing_psms[selected_rows, :]
            end
        end
        if write_intermediate_chromatogram_debug_plots(params)
            debug_write_target_chromatogram_plots(
                chromatograms,
                passing_psms,
                params.min_fraction_transmitted,
                params.wh_smoothing_strength,
                getDataOutDir(search_context),
                ms_file_idx,
                getFileIdToName(getMSData(search_context), ms_file_idx),
            )
        end
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
    output_path = getPassingPsms(getMSData(search_context))[ms_file_idx]
    if params.match_between_runs && isfile(output_path * PASS1_SIDECAR_SUFFIX)
        writeArrow(output_path, passing_psms)
        return nothing
    end

    # Process final PSMs
    process_final_psms!(
        passing_psms,
        search_context,
        parsed_fname,
        ms_file_idx
    )
    # Save results
    writeArrow(output_path, passing_psms)
    return nothing
end

function reset_results!(results::IntegrateChromatogramSearchResults)
    results.psms[] = DataFrame()
    GC.gc()
    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    passing_paths = String[
        path for path in getPassingPsms(getMSData(search_context))
        if !isempty(path) && isfile(path)
    ]
    if params.match_between_runs &&
       any(path -> isfile(path * PASS1_SIDECAR_SUFFIX), passing_paths)
        summary = finalize_mbr_after_chromatogram_integration!(
            passing_paths,
            getPrecursors(getSpecLib(search_context)),
            getFragmentLookupTable(getSpecLib(search_context));
            q_value_threshold = params.q_value_threshold,
            donor_q_threshold = MBR_DONOR_Q_THRESHOLD,
            min_pep_points_per_bin = params.pep_bin_size,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context),
            merged_path = joinpath(getDataOutDir(search_context), "merged_quant.arrow"),
        )
        @debug_l1 "Post-integration MBR completed: files=$(summary.n_files), " *
                  "candidates=$(summary.n_candidates), recovered=$(summary.n_recovered)"

        for (ms_file_idx, path) in enumerate(getPassingPsms(getMSData(search_context)))
            (isempty(path) || !isfile(path)) && continue
            psms = DataFrame(Tables.columntable(Arrow.Table(path)))
            if nrow(psms) == 0 || ncol(psms) == 0
                writeArrow(path, psms)
                continue
            end
            process_final_psms!(
                psms,
                search_context,
                getFileIdToName(getMSData(search_context), ms_file_idx),
                ms_file_idx,
            )
            writeArrow(path, psms)
        end
    end
    return nothing
end
