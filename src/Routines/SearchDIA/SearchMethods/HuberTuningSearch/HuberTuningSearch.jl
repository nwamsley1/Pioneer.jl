struct HuberTuningSearch <: TuningMethod end

struct HuberTuningSearchResults <: SearchResults
    huber_delta::Base.Ref{Float32}
    huber_histogram::Dict{Int, Int}
    n_observations::Base.RefValue{Int}
    n_selected::Base.RefValue{Int}
    n_curves::Base.RefValue{Int}
end

HuberTuningSearchResults(delta, histogram, observations) =
    HuberTuningSearchResults(delta, histogram, observations, Ref(0), Ref(0))

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
        Dict{Int, Int}(),
        Ref(0),
    )
end

function execute_search(::HuberTuningSearch, search_context::SearchContext, params::PioneerParameters)
    tuning_params = HuberTuningSearchParameters(params)
    tuning_params.enabled || return nothing
    Random.seed!(1844)
    results = init_search_results(tuning_params, search_context)
    winners = search_context.huber_calibration_winners
    files = huber_winner_files(winners)
    @debug_l1 "Global Huber selection: winners=$(count(w -> w.file_idx != 0, winners)) files=$(length(files))"
    try
        for file_idx in ProgressBar(files)
            idx = Int64(file_idx)
            spectra = getMSData(getMSData(search_context), idx)
            process_file!(results, tuning_params, search_context, idx, spectra)
        end
        summarize_results!(results, tuning_params, search_context)
    finally
        search_context.huber_calibration_winners = HuberCalibrationWinner[]
    end
    return nothing
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
        global_huber_psms(passing_psms, search_context.huber_calibration_winners, ms_file_idx),
        params.max_psms_for_huber,
    )
    isempty(calibration_psms) && return results
    results.n_selected[] += nrow(calibration_psms)

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
        ms_file_idx;
        competing_psms = passing_psms,
    )
    results.n_curves[] += nrow(tuning_psms) ÷ length(params.delta_grid)
    accumulate_huber_histogram!(
        results.huber_histogram, tuning_psms, params.delta_grid, params.min_pct_diff,
    )
    results.n_observations[] += nrow(tuning_psms)
    if ms_file_idx % 100 == 0
        @debug_l1 "Huber calibration summary: file_idx=$ms_file_idx selected_psms=$(results.n_selected[]) evaluated_curves=$(results.n_curves[]) delta_evaluations=$(results.n_observations[]) accepted_curves=$(sum(values(results.huber_histogram); init=0)) retained_bins=$(length(results.huber_histogram))"
    end

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
    if results.n_observations[] == 0
        results.huber_delta[] = fallback_delta
        setHuberDelta!(search_context, fallback_delta)
        @user_warn "No Huber calibration observations found; using default delta $(fallback_delta)"
        return nothing
    end

    try
        optimal_delta = estimate_optimal_huber_delta(results.huber_histogram)
        results.huber_delta[] = optimal_delta
        setHuberDelta!(search_context, optimal_delta)
        @debug_l1 "Global Huber delta calibration selected delta=$(optimal_delta) from selected_psms=$(results.n_selected[]) evaluated_curves=$(results.n_curves[]) delta_evaluations=$(results.n_observations[]) accepted_curves=$(sum(values(results.huber_histogram); init=0)); retained_bins=$(length(results.huber_histogram))"
    catch e
        results.huber_delta[] = fallback_delta
        setHuberDelta!(search_context, fallback_delta)
        @user_warn "Failed to determine optimal Huber delta, using default delta $(fallback_delta)" exception=e
    end

    return nothing
end
