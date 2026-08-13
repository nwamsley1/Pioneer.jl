struct HuberTuningSearch <: TuningMethod end

struct HuberTuningSearchResults <: SearchResults
    huber_delta::Base.Ref{Float32}
    tuning_psms::Vector{DataFrame}
end

struct HuberTuningSearchParameters{C<:IntegrateChromatogramSearchParameters} <: FragmentIndexSearchParameters
    chromatogram_params::C
    enabled::Bool
    base_solver::HuberSolver
    delta_grid::Vector{Float32}
    min_pct_diff::Float32
    max_psms_for_huber::Int64

    function HuberTuningSearchParameters(params::PioneerParameters)
        chromatogram_params = IntegrateChromatogramSearchParameters(params)
        enabled = chromatogram_params.deconvolution_solver isa HuberSolver
        base_solver = enabled ?
            chromatogram_params.deconvolution_solver::HuberSolver :
            default_chromatogram_integration_huber_solver()

        new{typeof(chromatogram_params)}(
            chromatogram_params,
            enabled,
            base_solver,
            chromatogram_huber_delta_grid(base_solver.delta),
            DEFAULT_HUBER_MIN_PCT_DIFF,
            Int64(DEFAULT_MAX_PSMS_FOR_HUBER),
        )
    end
end

get_parameters(::HuberTuningSearch, params::Any) = HuberTuningSearchParameters(params)

function init_search_results(::HuberTuningSearchParameters, ::SearchContext)
    return HuberTuningSearchResults(
        Ref(300.0f0),
        DataFrame[],
    )
end

function process_file!(
    results::HuberTuningSearchResults,
    params::HuberTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
)
    params.enabled || return results

    rt_index_path = getRtIndex(getMSData(search_context), ms_file_idx)
    passing_psms_path = getPassingPsms(getMSData(search_context), ms_file_idx)
    if isempty(rt_index_path) || isempty(passing_psms_path)
        return results
    end

    passing_psms = DataFrame(Tables.columntable(Arrow.Table(passing_psms_path)))
    isempty(passing_psms) && return results

    calibration_psms = select_huber_calibration_psms(
        passing_psms,
        params.max_psms_for_huber,
    )
    isempty(calibration_psms) && return results

    rt_index = buildRtIndex(
        DataFrame(Arrow.Table(rt_index_path)),
        bin_rt_size = 0.1,
    )

    tuning_psms = perform_huber_calibration_search(
        spectra,
        calibration_psms,
        rt_index,
        search_context,
        params,
        ms_file_idx,
    )
    nrow(tuning_psms) > 0 && push!(results.tuning_psms, tuning_psms)

    return results
end

function process_search_results!(
    ::HuberTuningSearchResults,
    ::HuberTuningSearchParameters,
    ::SearchContext,
    ::Int64,
    ::MassSpecData,
)
    return nothing
end

function reset_results!(::HuberTuningSearchResults)
    return nothing
end

function summarize_results!(
    results::HuberTuningSearchResults,
    params::HuberTuningSearchParameters,
    search_context::SearchContext,
)
    params.enabled || return nothing

    fallback_delta = params.base_solver.delta
    if isempty(results.tuning_psms)
        results.huber_delta[] = fallback_delta
        setHuberDelta!(search_context, fallback_delta)
        @user_warn "No Huber calibration observations found; using default delta $(fallback_delta)"
        return nothing
    end

    try
        all_psms = vcat(results.tuning_psms...)
        optimal_delta = estimate_optimal_huber_delta(
            all_psms,
            params.delta_grid,
            params.min_pct_diff,
        )
        results.huber_delta[] = optimal_delta
        setHuberDelta!(search_context, optimal_delta)
        @debug_l1 "Global Huber delta calibration selected delta=$(optimal_delta) from $(nrow(all_psms)) observations"
    catch e
        results.huber_delta[] = fallback_delta
        setHuberDelta!(search_context, fallback_delta)
        @user_warn "Failed to determine optimal Huber delta, using default delta $(fallback_delta)" exception=e
    end

    return nothing
end
