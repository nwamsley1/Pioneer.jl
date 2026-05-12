# MS1 feature restoration for the MainSearch -> PrecursorScoringSearch handoff.
#
# The historical SecondPassSearch wrote precursor-isotope evidence into
# `*_ms1` columns before the final LightGBM model. MainSearch now owns that
# handoff, so this file provides a lightweight MS1 isotope matcher and the
# schema join needed by the final scorer.

const MS1_BASE_SCHEMA = [
    (:m0, Bool, false),
    (:n_iso, UInt8, UInt8(0)),
    (:big_iso, UInt8, UInt8(0)),
    (:m0_error, Float16, Float16(-1)),
    (:error, Float16, Float16(-1)),
    (:spectral_contrast, Float16, Float16(-1)),
    (:fitted_spectral_contrast, Float16, Float16(-1)),
    (:gof, Float16, Float16(-1)),
    (:max_matched_residual, Float16, Float16(-1)),
    (:max_unmatched_residual, Float16, Float16(-1)),
    (:fitted_manhattan_distance, Float16, Float16(-1)),
    (:matched_ratio, Float16, Float16(-1)),
    (:weight, Float32, Float32(-1)),
    (:ms_file_idx, UInt32, UInt32(0)),
    (:scan_idx, UInt32, UInt32(0)),
    (:rt, Float32, Float32(-1)),
    (:rt_max_intensity, Float32, Float32(-1)),
    (:rt_diff_max_intensity, Float32, Float32(-1)),
    (:pair_idx, UInt32, UInt32(0)),
]

const MS1_COMPUTED_COLUMNS = [:ms1_ms2_rt_diff, :ms1_features_missing]

@inline _ms1_col(name::Symbol) = Symbol(String(name) * "_ms1")

@inline function _ms1_float16(x)
    xf = Float32(x)
    return Float16(isfinite(xf) ? xf : 0.0f0)
end

function _empty_ms1_feature_frame()
    df = DataFrame(precursor_idx = UInt32[])
    for (name, typ, _) in MS1_BASE_SCHEMA
        name === :precursor_idx && continue
        df[!, name] = typ[]
    end
    return df
end

function _coerce_ms1_base_schema!(ms1_psms::DataFrame)
    for (name, typ, sentinel) in MS1_BASE_SCHEMA
        if hasproperty(ms1_psms, name)
            ms1_psms[!, name] = typ.(coalesce.(ms1_psms[!, name], sentinel))
        else
            ms1_psms[!, name] = fill(sentinel, nrow(ms1_psms))
        end
        disallowmissing!(ms1_psms, name)
    end
    return ms1_psms
end

function _ms1_expected_column_order(psms::DataFrame)
    all_cols = propertynames(psms)
    ms1_cols_ordered = [_ms1_col(name) for (name, _, _) in MS1_BASE_SCHEMA]
    ms1_related = Set(vcat(ms1_cols_ordered, MS1_COMPUTED_COLUMNS))
    non_ms1_cols = [c for c in all_cols if !(c in ms1_related)]
    return vcat(non_ms1_cols, ms1_cols_ordered, MS1_COMPUTED_COLUMNS)
end

function _join_ms1_features!(
    psms::DataFrame,
    ms1_psms::DataFrame,
    rt_to_irt_model,
)
    if nrow(psms) == 0
        for (name, typ, _) in MS1_BASE_SCHEMA
            psms[!, _ms1_col(name)] = typ[]
        end
        psms[!, :ms1_ms2_rt_diff] = Float32[]
        psms[!, :ms1_features_missing] = Bool[]
        return psms
    end

    if nrow(ms1_psms) > 0
        _coerce_ms1_base_schema!(ms1_psms)
        ms1_psms = combine(groupby(ms1_psms, :precursor_idx), names(ms1_psms, Not(:precursor_idx)) .=> first .=> names(ms1_psms, Not(:precursor_idx)))
        psms_tmp = leftjoin(
            psms,
            ms1_psms,
            on = :precursor_idx,
            makeunique = true,
            renamecols = "" => "_ms1",
        )
        miss_mask = ismissing.(psms_tmp[!, _ms1_col(first(MS1_BASE_SCHEMA)[1])])
        empty!(psms)
        append!(psms, psms_tmp; cols=:union)
    else
        miss_mask = trues(nrow(psms))
        for (name, typ, sentinel) in MS1_BASE_SCHEMA
            psms[!, _ms1_col(name)] = fill(sentinel, nrow(psms))
        end
    end

    for (name, typ, sentinel) in MS1_BASE_SCHEMA
        col = _ms1_col(name)
        if hasproperty(psms, col)
            psms[!, col] = typ.(coalesce.(psms[!, col], sentinel))
        else
            psms[!, col] = fill(sentinel, nrow(psms))
        end
        disallowmissing!(psms, col)
    end

    rt_ms1 = psms[!, :rt_ms1]
    rt_ms2 = psms[!, :rt]
    ms1_ms2_rt_diff = Vector{Float32}(undef, nrow(psms))
    @inbounds for i in eachindex(rt_ms1)
        if rt_ms1[i] == Float32(-1)
            ms1_ms2_rt_diff[i] = Float32(-1)
        else
            ms1_ms2_rt_diff[i] = abs(Float32(rt_to_irt_model(rt_ms2[i])) - Float32(rt_to_irt_model(rt_ms1[i])))
        end
    end
    psms[!, :ms1_ms2_rt_diff] = ms1_ms2_rt_diff
    psms[!, :ms1_features_missing] = miss_mask

    select!(psms, _ms1_expected_column_order(psms))
    return psms
end

function _score_ms1_candidate_scan(
    isotopes::AbstractVector{<:Isotope{Float32}},
    scan_mz::AbstractArray{Union{Missing, Float32}},
    scan_int::AbstractArray{Union{Missing, Float32}},
    mem::AbstractMassErrorModel,
    corrected::Vector{Float32},
    obs_low::Vector{Float32},
    obs_high::Vector{Float32},
)
    n_peaks = prepare_scan_peaks!(corrected, obs_low, obs_high, mem, scan_mz, scan_int)
    return _score_ms1_candidate_scan_prepared(
        isotopes, scan_int, mem, corrected, obs_low, obs_high, n_peaks
    )
end

function _score_ms1_candidate_scan_prepared(
    isotopes::AbstractVector{<:Isotope{Float32}},
    scan_int::AbstractArray{Union{Missing, Float32}},
    mem::AbstractMassErrorModel,
    corrected::Vector{Float32},
    obs_low::Vector{Float32},
    obs_high::Vector{Float32},
    n_peaks::Int,
)
    n_peaks == 0 && return _missing_ms1_scan_score()

    sum_pred2 = 0.0f0
    sum_pred_obs = 0.0f0
    sum_obs2 = 0.0f0
    sum_obs = 0.0f0
    matched_pred = 0.0f0
    unmatched_pred = 0.0f0
    ppm_abs_sum = 0.0f0
    m0_error = Float32(-1)
    n_iso = UInt8(0)
    big_iso = UInt8(0)
    m0 = false

    preds = Vector{Float32}(undef, length(isotopes))
    obs = Vector{Float32}(undef, length(isotopes))
    matched = falses(length(isotopes))

    @inbounds for (i, iso) in enumerate(isotopes)
        pred = max(Float32(getIntensity(iso)), 0.0f0)
        iso_mz = Float32(getMZ(iso))
        half_width = conservative_half_width(mem, iso_mz)
        low, high = match_window(iso_mz, half_width)
        start_idx = bsearch_hybrid(corrected, low, 1, n_peaks)
        peak_idx, _, matched_mz = start_idx <= n_peaks ?
            scan_for_nearest_in_window(corrected, obs_low, obs_high, start_idx, n_peaks, iso_mz, high) :
            (0, typemax(Float32), 0.0f0)

        obs_i = peak_idx > 0 ? Float32(coalesce(scan_int[peak_idx], 0.0f0)) : 0.0f0
        preds[i] = pred
        obs[i] = obs_i
        matched[i] = peak_idx > 0

        sum_pred2 += pred * pred
        sum_obs2 += obs_i * obs_i
        sum_obs += obs_i

        if peak_idx > 0
            iso_idx = getIsoIdx(iso)
            n_iso += one(UInt8)
            big_iso = max(big_iso, iso_idx)
            matched_pred += pred
            sum_pred_obs += pred * obs_i
            ppm_err = compute_ppm_err(iso_mz, matched_mz)
            ppm_abs = abs(ppm_err)
            ppm_abs_sum += ppm_abs
            if iso_idx == one(UInt8)
                m0 = true
                m0_error = ppm_abs
            end
        else
            unmatched_pred += pred
        end
    end

    weight = sum_pred2 > 0.0f0 ? sum_pred_obs / sum_pred2 : 0.0f0
    if !(m0 && n_iso >= UInt8(2) && weight >= 1.0f-6)
        return _missing_ms1_scan_score()
    end

    sum_fitted = 0.0f0
    sum_fitted2 = 0.0f0
    fitted_dot_obs = 0.0f0
    manhattan = 0.0f0
    max_matched_resid = 0.0f0
    max_unmatched_fitted = 0.0f0
    matched_fitted_sum = 0.0f0

    @inbounds for i in eachindex(preds)
        fitted = max(weight * preds[i], 0.0f0)
        resid = abs(fitted - obs[i])
        sum_fitted += fitted
        sum_fitted2 += fitted * fitted
        fitted_dot_obs += fitted * obs[i]
        manhattan += resid
        if matched[i]
            max_matched_resid = max(max_matched_resid, resid)
            matched_fitted_sum += fitted
        else
            max_unmatched_fitted = max(max_unmatched_fitted, fitted)
        end
    end

    eps_score = 1.0f-10
    gof = sum_fitted > 0.0f0 ? -log2(manhattan / sum_fitted + eps_score) : 0.0f0
    max_matched_residual = matched_fitted_sum > 0.0f0 ?
        -log2(max_matched_resid / matched_fitted_sum + eps_score) : 0.0f0
    max_unmatched_residual = sum_fitted > 0.0f0 && max_unmatched_fitted > 0.0f0 ?
        -log2(max_unmatched_fitted / sum_fitted + eps_score) : 0.0f0
    fitted_manhattan_distance = sum_obs > 0.0f0 ?
        -log2(manhattan / sum_obs + eps_score) : 0.0f0
    spectral_denom = sqrt(sum_pred2 * sum_obs2)
    spectral_contrast = spectral_denom > 0.0f0 ? sum_pred_obs / spectral_denom : 0.0f0
    fitted_spectral_denom = sqrt(sum_fitted2 * sum_obs2)
    fitted_spectral_contrast = fitted_spectral_denom > 0.0f0 ? fitted_dot_obs / fitted_spectral_denom : 0.0f0
    matched_ratio = matched_pred > 0.0f0 && unmatched_pred > 0.0f0 ?
        min(log2(matched_pred / unmatched_pred), 10.0f0) : 0.0f0

    return (
        has_features = true,
        m0 = m0,
        n_iso = n_iso,
        big_iso = big_iso,
        m0_error = _ms1_float16(m0_error),
        error = _ms1_float16(log2(ppm_abs_sum + 1.0f-6)),
        spectral_contrast = _ms1_float16(spectral_contrast),
        fitted_spectral_contrast = _ms1_float16(fitted_spectral_contrast),
        gof = _ms1_float16(gof),
        max_matched_residual = _ms1_float16(max_matched_residual),
        max_unmatched_residual = _ms1_float16(max_unmatched_residual),
        fitted_manhattan_distance = _ms1_float16(fitted_manhattan_distance),
        matched_ratio = _ms1_float16(matched_ratio),
        weight = Float32(isfinite(weight) ? weight : 0.0f0),
    )
end

function _missing_ms1_scan_score()
    return (
        has_features = false,
        m0 = false,
        n_iso = UInt8(0),
        big_iso = UInt8(0),
        m0_error = Float16(-1),
        error = Float16(-1),
        spectral_contrast = Float16(-1),
        fitted_spectral_contrast = Float16(-1),
        gof = Float16(-1),
        max_matched_residual = Float16(-1),
        max_unmatched_residual = Float16(-1),
        fitted_manhattan_distance = Float16(-1),
        matched_ratio = Float16(-1),
        weight = Float32(-1),
    )
end

function _candidate_rt_index(psms::DataFrame, prec_mz_arr)
    candidate = unique(psms[:, [:precursor_idx, :rt]])
    sort!(candidate, :rt)
    prec_ids = UInt32.(candidate[!, :precursor_idx])
    rts = Float32.(candidate[!, :rt])
    mzs = Float32[prec_mz_arr[pid] for pid in prec_ids]
    return buildRtIndex(rts, mzs, prec_ids, 0.1f0)
end

function _build_ms1_isotope_dict(candidate_ids, prec_mz_arr, prec_charge_arr, prec_sulfur_arr, iso_splines)
    isotopes_dict = Dict{UInt32, Vector{Isotope{Float32}}}()
    for pid in candidate_ids
        charge = prec_charge_arr[pid]
        mono_mz = Float32(prec_mz_arr[pid])
        precursor_mass = (mono_mz * charge) - charge
        sulfur_idx = min(Int64(prec_sulfur_arr[pid]), 5)
        isotopes = Vector{Isotope{Float32}}(undef, 5)
        iso_mz = mono_mz
        @inbounds for iso_idx in 1:5
            intensity = max(Float32(iso_splines(sulfur_idx, iso_idx - 1, precursor_mass)), 0.0f0)
            isotopes[iso_idx] = Isotope(iso_mz, intensity, UInt8(iso_idx), UInt32(pid))
            iso_mz += Float32(NEUTRON / charge)
        end
        isotopes_dict[UInt32(pid)] = isotopes
    end
    return isotopes_dict
end

@inline function _ms1_local_rt_tol(rt_to_irt_model, irt_tol::Float32, rt::Float32)
    if !isfinite(irt_tol)
        return Float32(Inf)
    end
    h = 0.1f0
    local_slope = abs((Float32(rt_to_irt_model(rt + h)) - Float32(rt_to_irt_model(rt - h))) / (2.0f0 * h))
    return irt_tol / max(local_slope, 0.01f0)
end

function _mainsearch_ms1_mass_error_model(search_context::SearchContext, ms_file_idx::Int64)
    if haskey(search_context.ms1_mass_error_model, ms_file_idx)
        return search_context.ms1_mass_error_model[ms_file_idx]
    end
    return getMassErrorModel(search_context, ms_file_idx)
end

function _collect_ms1_feature_psms(
    psms::DataFrame,
    spectra::MassSpecData,
    search_context::SearchContext,
    ms_file_idx::Int64,
)
    nrow(psms) == 0 && return _empty_ms1_feature_frame()

    precursors = getPrecursors(getSpecLib(search_context))
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)
    prec_pair_idxs = getPairIdx(precursors)
    candidate_ids = unique(UInt32.(psms[!, :precursor_idx]))
    candidate_set = Set(candidate_ids)
    ms2_rt_lookup = Dict{UInt32, Float32}(UInt32(row.precursor_idx) => Float32(row.rt) for row in eachrow(psms))

    isotopes_dict = _build_ms1_isotope_dict(
        candidate_ids,
        prec_mz_arr,
        prec_charge_arr,
        prec_sulfur_arr,
        getIsoSplines(getSearchData(search_context)[1]),
    )

    rt_index = _candidate_rt_index(psms, prec_mz_arr)
    isempty(rt_index.rt_bins) && return _empty_ms1_feature_frame()

    mem = _mainsearch_ms1_mass_error_model(search_context, ms_file_idx)
    rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = haskey(getIrtErrors(search_context), ms_file_idx) ?
        getIrtErrors(search_context)[ms_file_idx] : Float32(Inf)

    corrected = Float32[]
    obs_low = Float32[]
    obs_high = Float32[]
    best_by_prec = Dict{UInt32, NamedTuple}()
    max_weight_rt = Dict{UInt32, Tuple{Float32, Float32}}()
    best_rt_diff = Dict{UInt32, Float32}()

    for scan_idx in 1:length(spectra)
        getMsOrder(spectra, scan_idx) == 1 || continue
        rt = Float32(getRetentionTime(spectra, scan_idx))
        rt_tol = _ms1_local_rt_tol(rt_to_irt_model, irt_tol, rt)
        rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
        rt_start <= rt_stop || continue

        scan_mz = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        n_peaks = prepare_scan_peaks!(corrected, obs_low, obs_high, mem, scan_mz, scan_int)
        n_peaks == 0 && continue

        for rt_bin_idx in rt_start:rt_stop
            for (pid, _) in rt_index.rt_bins[rt_bin_idx].prec
                pid in candidate_set || continue
                score = _score_ms1_candidate_scan_prepared(
                    isotopes_dict[pid], scan_int, mem,
                    corrected, obs_low, obs_high, n_peaks,
                )
                score.has_features || continue

                ms2_rt = ms2_rt_lookup[pid]
                rt_diff = abs(rt - ms2_rt)
                current_max = get(max_weight_rt, pid, (-Float32(Inf), Float32(-1)))
                if score.weight > first(current_max)
                    max_weight_rt[pid] = (score.weight, rt)
                end

                if !haskey(best_by_prec, pid) || rt_diff < best_rt_diff[pid]
                    pair_idx = extract_pair_idx(prec_pair_idxs, pid)
                    best_by_prec[pid] = (
                        precursor_idx = pid,
                        m0 = score.m0,
                        n_iso = score.n_iso,
                        big_iso = score.big_iso,
                        m0_error = score.m0_error,
                        error = score.error,
                        spectral_contrast = score.spectral_contrast,
                        fitted_spectral_contrast = score.fitted_spectral_contrast,
                        gof = score.gof,
                        max_matched_residual = score.max_matched_residual,
                        max_unmatched_residual = score.max_unmatched_residual,
                        fitted_manhattan_distance = score.fitted_manhattan_distance,
                        matched_ratio = score.matched_ratio,
                        weight = score.weight,
                        ms_file_idx = UInt32(ms_file_idx),
                        scan_idx = UInt32(scan_idx),
                        rt = rt,
                        rt_max_intensity = get(max_weight_rt, pid, (score.weight, rt))[2],
                        rt_diff_max_intensity = abs(get(max_weight_rt, pid, (score.weight, rt))[2] - ms2_rt),
                        pair_idx = pair_idx,
                    )
                    best_rt_diff[pid] = rt_diff
                end
            end
        end
    end

    if isempty(best_by_prec)
        return _empty_ms1_feature_frame()
    end

    rows = collect(values(best_by_prec))
    raw_ms1 = DataFrame(rows)
    for row_idx in 1:nrow(raw_ms1)
        pid = raw_ms1[row_idx, :precursor_idx]
        if haskey(max_weight_rt, pid)
            _, rt_max = max_weight_rt[pid]
            raw_ms1[row_idx, :rt_max_intensity] = rt_max
            raw_ms1[row_idx, :rt_diff_max_intensity] = abs(rt_max - ms2_rt_lookup[pid])
        end
    end

    return _expand_ms1_features_to_candidate_pairs(raw_ms1, candidate_ids, prec_pair_idxs)
end

function _expand_ms1_features_to_candidate_pairs(ms1_psms::DataFrame, candidate_ids, prec_pair_idxs)
    nrow(ms1_psms) == 0 && return ms1_psms

    best_by_pair = Dict{UInt32, Int}()
    for i in 1:nrow(ms1_psms)
        pid = ms1_psms[i, :precursor_idx]
        pair_idx = extract_pair_idx(prec_pair_idxs, pid)
        pair_key = iszero(pair_idx) ? pid : pair_idx
        if !haskey(best_by_pair, pair_key) || ms1_psms[i, :weight] > ms1_psms[best_by_pair[pair_key], :weight]
            best_by_pair[pair_key] = i
        end
    end

    rows = NamedTuple[]
    for pid in candidate_ids
        pair_idx = extract_pair_idx(prec_pair_idxs, pid)
        pair_key = iszero(pair_idx) ? pid : pair_idx
        haskey(best_by_pair, pair_key) || continue
        src = best_by_pair[pair_key]
        push!(rows, (
            precursor_idx = UInt32(pid),
            m0 = ms1_psms[src, :m0],
            n_iso = ms1_psms[src, :n_iso],
            big_iso = ms1_psms[src, :big_iso],
            m0_error = ms1_psms[src, :m0_error],
            error = ms1_psms[src, :error],
            spectral_contrast = ms1_psms[src, :spectral_contrast],
            fitted_spectral_contrast = ms1_psms[src, :fitted_spectral_contrast],
            gof = ms1_psms[src, :gof],
            max_matched_residual = ms1_psms[src, :max_matched_residual],
            max_unmatched_residual = ms1_psms[src, :max_unmatched_residual],
            fitted_manhattan_distance = ms1_psms[src, :fitted_manhattan_distance],
            matched_ratio = ms1_psms[src, :matched_ratio],
            weight = ms1_psms[src, :weight],
            ms_file_idx = ms1_psms[src, :ms_file_idx],
            scan_idx = ms1_psms[src, :scan_idx],
            rt = ms1_psms[src, :rt],
            rt_max_intensity = ms1_psms[src, :rt_max_intensity],
            rt_diff_max_intensity = ms1_psms[src, :rt_diff_max_intensity],
            pair_idx = pair_idx,
        ))
    end

    return isempty(rows) ? _empty_ms1_feature_frame() : DataFrame(rows)
end

function add_ms1_features!(
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
)
    ms1_psms = _collect_ms1_feature_psms(psms, spectra, search_context, ms_file_idx)
    _join_ms1_features!(psms, ms1_psms, getRtIrtModel(search_context, ms_file_idx))
    return psms
end
