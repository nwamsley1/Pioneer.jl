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
    boundary_selection_method::String
    learned_boundary_max_train_groups::Int64
    learned_boundary_min_positive::Int64
    learned_boundary_min_negative::Int64
    learned_boundary_min_crossrun_refs::Int64
    learned_boundary_max_isotope_trace_groups_per_file::Int64

    function IntegrateChromatogramSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        search_params = params.search
        chrom_params = params.optimization.chromatogram_integration

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
        boundary_selection_method = String(get(
            chrom_params,
            :boundary_selection_method,
            "learned",
        ))
        if !(boundary_selection_method in ("second_derivative", "learned"))
            throw(ArgumentError(
                "Invalid optimization.chromatogram_integration.boundary_selection_method=$(repr(boundary_selection_method)); " *
                "expected \"second_derivative\" or \"learned\".",
            ))
        end
        learned_boundary_max_train_groups = Int64(get(
            chrom_params,
            :learned_boundary_max_train_groups,
            20_000,
        ))
        learned_boundary_min_positive = Int64(get(
            chrom_params,
            :learned_boundary_min_positive,
            25,
        ))
        learned_boundary_min_negative = Int64(get(
            chrom_params,
            :learned_boundary_min_negative,
            25,
        ))
        learned_boundary_min_crossrun_refs = Int64(get(
            chrom_params,
            :learned_boundary_min_crossrun_refs,
            1,
        ))
        learned_boundary_max_isotope_trace_groups_per_file = Int64(get(
            chrom_params,
            :learned_boundary_max_isotope_trace_groups_per_file,
            2_000,
        ))

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
            boundary_selection_method,
            learned_boundary_max_train_groups,
            learned_boundary_min_positive,
            learned_boundary_min_negative,
            learned_boundary_min_crossrun_refs,
            learned_boundary_max_isotope_trace_groups_per_file,
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::IntegrateChromatogramSearch, params::Any) = IntegrateChromatogramSearchParameters(params)

function write_intermediate_chromatogram_debug_plots(
    params::IntegrateChromatogramSearchParameters,
)
    return params.boundary_selection_method != "learned"
end

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

    @debug_l1 "Chromatogram integration solver for file $ms_file_idx: " *
              "$(chromatogram_integration_solver_label(calibrated_chromatogram_deconvolution_solver(search_context, params.deconvolution_solver)))"

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
    else
        # Combined-trace quantification still integrates all scan points together,
        # but we also keep isotope-capture states so the learned boundary model
        # can train on isotope-specific side traces.
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
    sort_chromatograms_for_integration!(chromatograms, params.isotope_tracetype)

    # Integrate chromatographic peaks for each precursor (skip if no chromatograms extracted)
    if nrow(chromatograms) > 0
        psm_isotopes_captured = nothing
        selected_quant_trace = nothing
        if seperateTraces(params.isotope_tracetype)
            selected_quant_trace = select_quant_trace_by_transmission(chromatograms)
            apply_quant_trace_selection!(
                passing_psms,
                selected_quant_trace,
            )
            psm_isotopes_captured = passing_psms[!, :isotopes_captured]
        end
        boundary_candidate_data = params.boundary_selection_method == "learned" ?
            Vector{Any}(nothing, nrow(passing_psms)) :
            nothing
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
            mainsearch_1pct_start_scan = hasproperty(passing_psms, :mainsearch_1pct_start_scan) ?
                passing_psms[!, :mainsearch_1pct_start_scan] :
                nothing,
            mainsearch_1pct_stop_scan = hasproperty(passing_psms, :mainsearch_1pct_stop_scan) ?
                passing_psms[!, :mainsearch_1pct_stop_scan] :
                nothing,
            rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx),
            boundary_candidate_data = boundary_candidate_data,
        )
        if boundary_candidate_data !== nothing
            extra_candidate_rows = if hasproperty(chromatograms, :isotopes_captured)
                collect_isotope_boundary_candidate_rows(
                    chromatograms,
                    passing_psms,
                    selected_quant_trace,
                    search_context,
                    ms_file_idx,
                    params.min_fraction_transmitted,
                    params.wh_smoothing_strength,
                    max_groups = params.learned_boundary_max_isotope_trace_groups_per_file,
                    rng = Random.MersenneTwister(1776 + ms_file_idx),
                )
            else
                NamedTuple[]
            end
            write_boundary_candidate_table!(
                boundary_candidate_data,
                passing_psms,
                search_context,
                ms_file_idx,
                extra_candidate_rows,
            )
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

function debug_write_selected_boundary_chromatogram_plots(
    selected_candidates::AbstractDataFrame,
    params::P,
    search_context::SearchContext,
) where {P<:IntegrateChromatogramSearchParameters}
    DEBUG_CONSOLE_LEVEL[] >= 1 || return nothing
    target_precursor_idxs = DEBUG_CHROM_TARGET_PRECURSOR_IDXS[]
    isempty(target_precursor_idxs) && return nothing
    nrow(selected_candidates) == 0 && return nothing

    selected_targets = selected_candidates[
        in.(UInt32.(selected_candidates.precursor_idx), Ref(target_precursor_idxs)),
        :,
    ]
    nrow(selected_targets) == 0 && return nothing

    ms_data = getMSData(search_context)
    for group in groupby(selected_targets, :ms_file_idx)
        ms_file_idx = Int(first(group.ms_file_idx))
        rt_index_path = getRtIndex(ms_data, ms_file_idx)
        passing_psms_path = getPassingPsms(ms_data, ms_file_idx)
        if isempty(rt_index_path) || isempty(passing_psms_path) ||
           !isfile(rt_index_path) || !isfile(passing_psms_path)
            debug_l1("chromatogram_debug_plot ms_file_idx=$(ms_file_idx) " *
                     "skipped=true reason=missing_final_debug_inputs")
            continue
        end

        target_group_precursors = Set(UInt32.(group.precursor_idx))
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(passing_psms_path)))
        psm_mask = in.(UInt32.(passing_psms.precursor_idx), Ref(target_group_precursors))
        passing_psms = passing_psms[psm_mask, :]
        nrow(passing_psms) == 0 && continue

        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(rt_index_path)),
            bin_rt_size = 0.1,
        )
        spectra = getMSData(ms_data, ms_file_idx)
        chromatograms = extract_chromatograms(
            spectra,
            passing_psms,
            rt_index,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM(),
        )
        nrow(chromatograms) == 0 && continue

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
        )
        sort_chromatograms_for_integration!(chromatograms, params.isotope_tracetype)

        debug_write_target_chromatogram_plots(
            chromatograms,
            passing_psms,
            params.min_fraction_transmitted,
            params.wh_smoothing_strength,
            getDataOutDir(search_context),
            ms_file_idx,
            getFileIdToName(ms_data, ms_file_idx);
            selected_boundary_candidates = group,
        )
    end

    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    params.boundary_selection_method == "learned" || return nothing

    candidates = load_boundary_candidate_tables(search_context)
    nrow(candidates) == 0 && return nothing

    training_candidates = hasproperty(candidates, :target) ?
        candidates[candidates.target .== true, :] :
        copy(candidates)
    if nrow(training_candidates) == 0
        @user_warn "Learned chromatogram boundary model skipped: no target candidate rows available; keeping second-derivative bounds."
        return nothing
    end

    models = train_boundary_candidate_models_by_cv_fold(
        training_candidates;
        max_groups = params.learned_boundary_max_train_groups,
        min_positive = params.learned_boundary_min_positive,
        min_negative = params.learned_boundary_min_negative,
        min_crossrun_refs = params.learned_boundary_min_crossrun_refs,
        rng = Random.MersenneTwister(1776),
    )

    if isempty(models) || all(isnothing, values(models))
        @user_warn "Learned chromatogram boundary models skipped: insufficient labeled candidate diversity; keeping second-derivative bounds."
        return nothing
    end

    log_boundary_cv_model_feature_importances(models)
    selected = select_boundary_candidate_rows_crossfold(candidates, models)
    log_boundary_candidate_category_tally(selected)
    if DEBUG_CONSOLE_LEVEL[] >= 1
        debug_candidates = copy(candidates)
        label_boundary_candidate_targets!(
            debug_candidates;
            min_crossrun_refs = params.learned_boundary_min_crossrun_refs,
        )
        score_boundary_candidates_crossfold!(debug_candidates, models)
        log_boundary_candidate_debug(debug_candidates, selected)
    end
    updated = apply_selected_boundary_candidates!(selected, search_context)
    debug_write_selected_boundary_chromatogram_plots(selected, params, search_context)
    @user_info "Learned chromatogram boundary model selected bounds for $updated precursor traces."
    return nothing
end
