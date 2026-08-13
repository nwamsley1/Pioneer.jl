const DEFAULT_HUBER_DELTA_EXP = 1.5f0
const DEFAULT_HUBER_DELTA_ITERS = 15
const DEFAULT_MAX_PSMS_FOR_HUBER = 5000
const DEFAULT_HUBER_MIN_PCT_DIFF = 10.0f0

function chromatogram_huber_delta_grid(
    delta0::Float32;
    delta_exp::Float32 = DEFAULT_HUBER_DELTA_EXP,
    delta_iters::Int64 = DEFAULT_HUBER_DELTA_ITERS,
)
    return Float32[delta0 * (delta_exp^i) for i in -4:(delta_iters + 6)]
end

function empty_huber_tuning_results()
    return DataFrame(
        ms_file_idx = UInt32[],
        precursor_idx = UInt32[],
        scan_idx = UInt32[],
        huber_delta = Float32[],
        weight = Float32[],
    )
end

function select_huber_calibration_psms(psms::DataFrame, max_psms::Int64)
    max_psms <= 0 && return psms[1:0, :]
    nrow(psms) <= max_psms && return psms

    psm_counts = combine(
        groupby(psms, :scan_idx, sort = false),
        nrow => :n_psms,
    )
    sort!(psm_counts, :n_psms, rev = true)

    scans_to_keep = Set{UInt32}()
    total_psms = 0
    for row in eachrow(psm_counts)
        push!(scans_to_keep, UInt32(row.scan_idx))
        total_psms += row.n_psms
        total_psms >= max_psms && break
    end

    return filter(row -> UInt32(row.scan_idx) in scans_to_keep, psms)
end

function huber_scan_precursor_mapping(psms::DataFrame)
    scan_to_prec = Dict{UInt32, Vector{UInt32}}()
    for row in eachrow(psms)
        scan_idx = UInt32(row.scan_idx)
        if !haskey(scan_to_prec, scan_idx)
            scan_to_prec[scan_idx] = UInt32[]
        end
        push!(scan_to_prec[scan_idx], UInt32(row.precursor_idx))
    end
    return scan_to_prec
end

function huber_precursor_rt_map(psms::DataFrame)
    precursor_rt_map = Dict{UInt32, Float32}()
    sizehint!(precursor_rt_map, nrow(psms))
    for row in eachrow(psms)
        precursor_rt_map[UInt32(row.precursor_idx)] = Float32(row.rt)
    end
    return precursor_rt_map
end

function perform_huber_calibration_search(
    spectra::MassSpecData,
    calibration_psms::DataFrame,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64,
)
    isempty(calibration_psms) && return empty_huber_tuning_results()

    thread_tasks = partition_scans(spectra, Threads.nthreads(), ms_order_select = 2)
    scan_to_prec = huber_scan_precursor_mapping(calibration_psms)
    precursor_rt_map = huber_precursor_rt_map(calibration_psms)
    precursor_set = Set(calibration_psms[!, :precursor_idx])  # shared read-only across threads (was N copies)

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            return process_huber_calibration_scans!(
                last(thread_task),
                precursor_set,
                precursor_rt_map,
                scan_to_prec,
                rt_index,
                spectra,
                search_context,
                search_data,
                params,
                ms_file_idx,
            )
        end
    end

    thread_results = fetch.(tasks)
    isempty(thread_results) && return empty_huber_tuning_results()
    return vcat(thread_results...)
end

function process_huber_calibration_scans!(
    scan_range::AbstractVector{Int64},
    precursors_passing::Set{UInt32},
    precursor_rt_map::Dict{UInt32, Float32},
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    rt_index::retentionTimeIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64,
)
    chrom_params = params.chromatogram_params
    tuning_results = empty_huber_tuning_results()

    Hs = getHsFused(search_data)
    weights = getTempWeights(search_data)
    colnorm2 = getColNorm2(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    fused_scratch = getFusedScratch(search_data)
    corr_mz = getScanCorrectedMz(search_data)
    obs_low = getScanObsLow(search_data)
    obs_high = getScanObsHigh(search_data)
    isotopes_buf = getIsotopes(search_data)
    prec_trans_buf = getPrecursorTransmission(search_data)
    id_to_col = getIdToCol(search_data)
    unscored_psms = getTuningUnscoredPsms(search_data)

    rt_bin_start, rt_bin_stop = 1, 1
    prec_mz = zero(Float32)
    precs_temp = getPrecIds(search_data)
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    has_rt_tol = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol = has_rt_tol ? getRtTolerance(search_context, ms_file_idx) : nothing
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    nce_model = getNceModel(search_context, ms_file_idx)
    mass_error_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)
    prec_irt_arr = getIrt(precursors)
    frag_lookup = getFragmentLookupTable(spec_lib)
    kind = FusedRTIndexed(chrom_params.prec_estimation, UInt8(chrom_params.max_frag_rank))
    huber_solvers = HuberSolver[
        with_chromatogram_huber_delta(params.base_solver, delta)
        for delta in params.delta_grid
    ]

    for scan_idx in scan_range
        ((scan_idx < 1) | (scan_idx > length(spectra))) && continue
        scan_key = UInt32(scan_idx)
        haskey(scan_to_prec, scan_key) || continue
        getMsOrder(spectra, scan_idx) ∉ chrom_params.spec_order && continue

        rt = getRetentionTime(spectra, scan_idx)
        if rt_binned_tol !== nothing
            rt_tol_local = get_rt_tol(rt_binned_tol, Float32(rt))
        else
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end

        rt_bin_start_new = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local, lt = (r, x) -> r.lb < x) - 1, 1)
        rt_bin_stop_new = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local, lt = (x, r) -> r.ub > x) + 1, length(rt_index.rt_bins))

        prec_mz_new = getCenterMz(spectra, scan_idx)
        if (rt_bin_start_new != rt_bin_start) || (rt_bin_stop_new != rt_bin_stop) || (prec_mz_new != prec_mz)
            rt_bin_start = rt_bin_start_new
            rt_bin_stop = rt_bin_stop_new
            prec_mz = prec_mz_new
        end

        quad_func = getQuadTransmissionFunction(
            quad_model,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx),
        )

        prec_temp_size = collect_rt_window_precursors!(
            precs_temp, rt_index, rt_bin_start, rt_bin_stop,
            precursors_passing, prec_mz_arr, prec_charge_arr, prec_sulfur_arr,
            getIsoSplines(search_data), quad_func, prec_trans_buf,
            chrom_params.isotope_err_bounds, chrom_params.min_fraction_transmitted,
            precursor_rt_map, Float32(rt), rt_binned_tol, rt_tol_local,
        )

        if prec_temp_size == 0
            reset!(id_to_col)
            reset!(Hs)
            continue
        end

        scan_mz = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(
            corr_mz, obs_low, obs_high,
            mass_error_model, scan_mz, scan_int, Float32(rt),
        )

        nmatches, nmisses = run_fused!(
            kind,
            Hs, unscored_psms, id_to_col, fused_scratch,
            corr_mz, obs_low, obs_high, peak_mz_len,
            isotopes_buf, prec_trans_buf,
            frag_lookup, nce_model,
            precs_temp, 1:prec_temp_size,
            prec_mz_arr, prec_charge_arr, prec_sulfur_arr, prec_irt_arr,
            getIsoSplines(search_data), quad_func, mass_error_model,
            scan_int, 0f0, Float32(Inf),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx)),
            chrom_params.n_frag_isotopes,
            chrom_params.isotope_err_bounds,
        )

        if nmatches > 2
            n_active_cols = n_active(id_to_col)
            if n_active_cols > length(weights)
                new_entries = n_active_cols - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
            end

            target_precursors = scan_to_prec[scan_key]
            for (delta, solver) in zip(params.delta_grid, huber_solvers)
                initialize_weights!(id_to_col, weights, precursor_weights)
                solve_deconvolution!(
                    solver,
                    Hs, residuals, weights, colnorm2,
                    getMu(search_data), getObserved(search_data),
                    chrom_params.max_iter_outer, chrom_params.max_diff,
                )
                update_precursor_weights!(id_to_col, weights, precursor_weights)

                for precursor_idx in target_precursors
                    col = id_to_col[precursor_idx]
                    iszero(col) && continue
                    push!(tuning_results, (
                        UInt32(ms_file_idx),
                        precursor_idx,
                        scan_key,
                        delta,
                        weights[col],
                    ))
                end
            end
        end

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    return tuning_results
end

function estimate_optimal_huber_delta(
    psms::DataFrame,
    delta_grid::Vector{Float32},
    min_pct_diff::Float32,
)
    delta_col = hasproperty(psms, :huber_delta) ? :huber_delta : Symbol("huber_δ")
    group_cols = hasproperty(psms, :ms_file_idx) ?
        [:ms_file_idx, :precursor_idx, :scan_idx] :
        [:precursor_idx, :scan_idx]

    grouped_psms = groupby(psms, group_cols)
    curves = combine(grouped_psms) do sdf
        process_huber_curve(sdf[!, :weight], sdf[!, delta_col])
    end

    filter!(row -> row.n == length(delta_grid), curves)
    filter!(row -> row.wdiff > (min_pct_diff / 100), curves)
    filter!(row -> !ismissing(row.huber50), curves)
    isempty(curves) && throw(ArgumentError("No complete Huber calibration curves passed filtering"))

    curves[!, :huber50] = ceil.(Int, curves[!, :huber50])
    huber_hist = combine(groupby(curves, :huber50), nrow)
    sort!(huber_hist, :huber50)
    huber_hist[!, :prob] = huber_hist[!, :nrow] ./ sum(huber_hist[!, :nrow])
    huber_hist[!, :cum_prob] = Float32.(cumsum(huber_hist[!, :prob]))

    return get_median_huber_delta(
        huber_hist[!, :cum_prob],
        huber_hist[!, :huber50],
    )
end

function process_huber_curve(
    weights::AbstractVector{Float32},
    huber_deltas::AbstractVector{Float32},
)
    order = sortperm(huber_deltas)
    sorted_weights = @view weights[order]
    sorted_deltas = @view huber_deltas[order]
    min_w, max_w = minimum(sorted_weights), maximum(sorted_weights)
    huber50 = missing
    w50 = min_w + (max_w - min_w) / 2

    if length(sorted_weights) > 1
        for i in 1:(length(sorted_weights) - 1)
            if (w50 >= sorted_weights[i]) && (w50 <= sorted_weights[i + 1])
                huber50 = sorted_deltas[i] + (sorted_deltas[i + 1] - sorted_deltas[i]) / 2
            elseif (w50 <= sorted_weights[i]) && (w50 >= sorted_weights[i + 1])
                huber50 = sorted_deltas[i] + (sorted_deltas[i + 1] - sorted_deltas[i]) / 2
            end
        end
    end

    return (
        min = min_w,
        max = max_w,
        n = length(sorted_weights),
        huber50 = huber50,
        w50 = w50,
        wdiff = iszero(min_w) ? Float32(Inf) : (max_w - min_w) / min_w,
    )
end

function get_median_huber_delta(
    cum_prob::Vector{Float32},
    delta::AbstractVector{<:Integer};
    cum_prob_threshold::Float32 = 0.5f0,
)
    n = length(cum_prob)
    n == 1 && return Float32(first(delta))

    for i in 1:(n - 1)
        if cum_prob[i + 1] >= cum_prob_threshold
            x1, x2 = cum_prob[i], cum_prob[i + 1]
            y1, y2 = delta[i], delta[i + 1]
            slope = (y2 - y1) / (x2 - x1)
            midpoint = (x2 + x1) / 2
            return Float32(y1 + slope * (midpoint - x1))
        end
    end

    @user_warn "Could not estimate Huber delta"
    return Float32(first(delta))
end
