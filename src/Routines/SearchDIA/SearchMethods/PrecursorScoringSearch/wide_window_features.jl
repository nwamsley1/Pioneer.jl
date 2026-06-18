const WIDE_WINDOW_CYCLE_MARGIN = 3
const WIDE_WINDOW_CORR_THRESHOLD = 0.7f0

@inline function _wide_window_pcor(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n < 2 && return 0f0
    mean_x = 0f0
    mean_y = 0f0
    @inbounds for i in 1:n
        mean_x += Float32(x[i])
        mean_y += Float32(y[i])
    end
    inv_n = 1f0 / Float32(n)
    mean_x *= inv_n
    mean_y *= inv_n

    sum_xx = 0f0
    sum_yy = 0f0
    sum_xy = 0f0
    @inbounds for i in 1:n
        dx = Float32(x[i]) - mean_x
        dy = Float32(y[i]) - mean_y
        sum_xx += dx * dx
        sum_yy += dy * dy
        sum_xy += dx * dy
    end
    denom = sqrt(sum_xx * sum_yy)
    return denom > 0f0 ? Float32(sum_xy / denom) : 0f0
end

@inline function _wide_window_zero_feature_values()
    return (
        wide_ms1_m0_candidate_fraction = 0f0,
        wide_frag_candidate_fraction = 0f0,
        wide_ms1_frag_sum_corr = 0f0,
        wide_frag_corr_mean = 0f0,
        wide_n_correlated_fragments = UInt8(0),
        wide_n_correlated_fragments_bitvec_rank = UInt16(0),
        wide_frag_corr_strength = 0f0,
        wide_frag_corr_effective_n = 0f0,
        wide_frag_corr_best_m0 = 0f0,
        wide_signal_support = 0f0,
    )
end

function _wide_window_feature_values(
    ms1_m0::AbstractVector,
    fragments::AbstractMatrix,
    candidate_mask::AbstractVector{Bool};
    bitvec_rank_table = nothing,
)
    n_scans = length(ms1_m0)
    n_scans == 0 && return _wide_window_zero_feature_values()

    n_frags = size(fragments, 2)
    frag_sum = zeros(Float32, n_scans)
    ms1_total = 0f0
    ms1_candidate = 0f0
    frag_total = 0f0
    frag_candidate = 0f0
    signal_scans = 0

    @inbounds for scan_i in 1:n_scans
        m0 = max(Float32(ms1_m0[scan_i]), 0f0)
        summed_frag = 0f0
        for frag_i in 1:n_frags
            summed_frag += max(Float32(fragments[scan_i, frag_i]), 0f0)
        end
        frag_sum[scan_i] = summed_frag
        ms1_total += m0
        frag_total += summed_frag
        if candidate_mask[scan_i]
            ms1_candidate += m0
            frag_candidate += summed_frag
        end
        (m0 > 0f0 || summed_frag > 0f0) && (signal_scans += 1)
    end

    has_signal = falses(n_frags)
    @inbounds for frag_i in 1:n_frags
        for scan_i in 1:n_scans
            if Float32(fragments[scan_i, frag_i]) > 0f0
                has_signal[frag_i] = true
                break
            end
        end
    end

    pair_corr = zeros(Float32, n_frags, n_frags)
    pair_sum = 0f0
    pair_count = 0
    @inbounds for frag_a in 1:n_frags
        has_signal[frag_a] || continue
        for frag_b in (frag_a + 1):n_frags
            has_signal[frag_b] || continue
            corr = _wide_window_pcor(view(fragments, :, frag_a), view(fragments, :, frag_b))
            pair_corr[frag_a, frag_b] = corr
            pair_corr[frag_b, frag_a] = corr
            pair_sum += corr
            pair_count += 1
        end
    end

    n_correlated = UInt8(0)
    correlated_mask = UInt16(0)
    corr_strength = 0f0
    corr_sumsq = 0f0
    other_sum = Vector{Float32}(undef, n_scans)
    @inbounds for frag_i in 1:n_frags
        has_signal[frag_i] || continue
        for scan_i in 1:n_scans
            other_sum[scan_i] = frag_sum[scan_i] - max(Float32(fragments[scan_i, frag_i]), 0f0)
        end
        corr = _wide_window_pcor(view(fragments, :, frag_i), other_sum)
        positive_corr = min(max(corr, 0f0), 1f0)
        corr_strength += positive_corr
        corr_sumsq += positive_corr * positive_corr
        if corr > WIDE_WINDOW_CORR_THRESHOLD
            n_correlated += UInt8(1)
            correlated_mask |= UInt16(1) << (frag_i - 1)
        end
    end

    best_frag = 0
    best_consensus = typemin(Float32)
    @inbounds for frag_i in 1:n_frags
        has_signal[frag_i] || continue
        consensus = 0f0
        n_pairs = 0
        for other_frag_i in 1:n_frags
            (other_frag_i == frag_i || !has_signal[other_frag_i]) && continue
            consensus += pair_corr[frag_i, other_frag_i]
            n_pairs += 1
        end
        avg_consensus = n_pairs > 0 ? Float32(consensus / n_pairs) : typemin(Float32)
        if avg_consensus > best_consensus
            best_consensus = avg_consensus
            best_frag = frag_i
        end
    end

    return (
        wide_ms1_m0_candidate_fraction = ms1_total > 0f0 ? Float32(ms1_candidate / ms1_total) : 0f0,
        wide_frag_candidate_fraction = frag_total > 0f0 ? Float32(frag_candidate / frag_total) : 0f0,
        wide_ms1_frag_sum_corr = _wide_window_pcor(ms1_m0, frag_sum),
        wide_frag_corr_mean = pair_count > 0 ? Float32(pair_sum / pair_count) : 0f0,
        wide_n_correlated_fragments = n_correlated,
        wide_n_correlated_fragments_bitvec_rank = _bitvec_pattern_rank(bitvec_rank_table, correlated_mask),
        wide_frag_corr_strength = corr_strength,
        wide_frag_corr_effective_n = corr_sumsq > 0f0 ? Float32((corr_strength * corr_strength) / corr_sumsq) : 0f0,
        wide_frag_corr_best_m0 = best_frag > 0 ? _wide_window_pcor(view(fragments, :, best_frag), ms1_m0) : 0f0,
        wide_signal_support = Float32(signal_scans / n_scans),
    )
end

@inline function _wide_peak_intensity(
    mz::AbstractVector{Float32},
    intens,
    target::Float32,
    ppm_tol::Float32,
)
    n_peaks = length(mz)
    n_peaks == 0 && return 0f0
    tol = target * ppm_tol * 1f-6
    low_mz = target - tol
    high_mz = target + tol
    start_idx = bsearch_hybrid(mz, low_mz, 1, n_peaks)
    best_diff = Inf32
    best_intensity = 0f0
    @inbounds for peak_i in start_idx:n_peaks
        peak_mz = mz[peak_i]
        peak_mz > high_mz && break
        intensity = intens[peak_i]
        ismissing(intensity) && continue
        diff = abs(peak_mz - target)
        if diff < best_diff
            best_diff = diff
            best_intensity = Float32(intensity)
        end
    end
    return best_intensity
end

@inline function _wide_peak_intensity(mz, intens, target::Float32, ppm_tol::Float32)
    n_peaks = length(mz)
    n_peaks == 0 && return 0f0
    tol = target * ppm_tol * 1f-6
    low_mz = target - tol
    high_mz = target + tol
    best_diff = Inf32
    best_intensity = 0f0
    @inbounds for peak_i in 1:n_peaks
        peak_mz = mz[peak_i]
        ismissing(peak_mz) && continue
        peak_mz_f32 = Float32(peak_mz)
        peak_mz_f32 < low_mz && continue
        peak_mz_f32 > high_mz && break
        intensity = intens[peak_i]
        ismissing(intensity) && continue
        diff = abs(peak_mz_f32 - target)
        if diff < best_diff
            best_diff = diff
            best_intensity = Float32(intensity)
        end
    end
    return best_intensity
end

@inline function _wide_fragment_peak_intensity(
    scan_corrected_mz::Vector{Float32},
    scan_obs_low::Vector{Float32},
    scan_obs_high::Vector{Float32},
    scan_int,
    n_peaks::Int,
    target::Float32,
    mem::AbstractMassErrorModel,
)
    n_peaks == 0 && return 0f0
    half_width = conservative_half_width(mem, target)
    conservative_low, conservative_high = match_window(target, half_width)
    start_idx = bsearch_hybrid(scan_corrected_mz, conservative_low, 1, n_peaks)
    best_peak, _, _ = scan_for_nearest_in_window(
        scan_corrected_mz,
        scan_obs_low,
        scan_obs_high,
        start_idx,
        n_peaks,
        target,
        conservative_high,
    )
    best_peak == 0 && return 0f0
    intensity = scan_int[best_peak]
    return ismissing(intensity) ? 0f0 : Float32(intensity)
end

@inline function _wide_fragment_mono_peak_intensity(
    scan_corrected_mz::Vector{Float32},
    scan_obs_low::Vector{Float32},
    scan_obs_high::Vector{Float32},
    scan_int,
    n_peaks::Int,
    frag::AltimeterFragment,
    mem::AbstractMassErrorModel,
)
    return _wide_fragment_peak_intensity(
        scan_corrected_mz,
        scan_obs_low,
        scan_obs_high,
        scan_int,
        n_peaks,
        Float32(getMz(frag)),
        mem,
    )
end

@inline function _wide_window_key(center_mz, isolation_width)
    (ismissing(center_mz) || ismissing(isolation_width)) && return (Int32(0), Int32(0))
    return (
        Int32(round(Int, Float32(center_mz) * 100f0)),
        Int32(round(Int, Float32(isolation_width) * 100f0)),
    )
end

function _wide_window_index_spectra(spectra)
    n_scans = length(spectra)
    window_scans = Dict{Tuple{Int32, Int32}, Vector{Int32}}()
    scan_to_window_key = fill((Int32(0), Int32(0)), n_scans)
    scan_to_window_pos = zeros(Int32, n_scans)
    ms1_scan_idxs = Int32[]
    ms1_scan_rts = Float32[]

    for scan_idx in 1:n_scans
        ms_order = getMsOrder(spectra, scan_idx)
        if ms_order == UInt8(1)
            push!(ms1_scan_idxs, Int32(scan_idx))
            push!(ms1_scan_rts, Float32(getRetentionTime(spectra, scan_idx)))
        elseif ms_order == UInt8(2)
            key = _wide_window_key(
                getCenterMz(spectra, scan_idx),
                getIsolationWidthMz(spectra, scan_idx),
            )
            key == (Int32(0), Int32(0)) && continue
            scans = get!(window_scans, key, Int32[])
            push!(scans, Int32(scan_idx))
            scan_to_window_key[scan_idx] = key
            scan_to_window_pos[scan_idx] = Int32(length(scans))
        end
    end

    scan_to_ms1 = zeros(Int32, n_scans)
    n_ms1 = length(ms1_scan_idxs)
    if n_ms1 > 0
        @inbounds for scan_idx in 1:n_scans
            scan_rt = Float32(getRetentionTime(spectra, scan_idx))
            pos = searchsortedfirst(ms1_scan_rts, scan_rt)
            scan_to_ms1[scan_idx] = if pos == 1
                ms1_scan_idxs[1]
            elseif pos > n_ms1
                ms1_scan_idxs[end]
            else
                after_delta = abs(ms1_scan_rts[pos] - scan_rt)
                before_delta = abs(ms1_scan_rts[pos - 1] - scan_rt)
                before_delta <= after_delta ? ms1_scan_idxs[pos - 1] : ms1_scan_idxs[pos]
            end
        end
    end

    return (
        window_scans = window_scans,
        scan_to_window_key = scan_to_window_key,
        scan_to_window_pos = scan_to_window_pos,
        scan_to_ms1 = scan_to_ms1,
    )
end

function _wide_expanded_window_scans(candidate_scans::Vector{Int32}, scan_index)
    bounds = Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}}()
    @inbounds for scan in candidate_scans
        pos = scan_index.scan_to_window_pos[Int(scan)]
        pos == 0 && continue
        key = scan_index.scan_to_window_key[Int(scan)]
        if haskey(bounds, key)
            low_pos, high_pos = bounds[key]
            bounds[key] = (min(low_pos, pos), max(high_pos, pos))
        else
            bounds[key] = (pos, pos)
        end
    end

    expanded = Int32[]
    for (key, (low_pos, high_pos)) in bounds
        scans = scan_index.window_scans[key]
        first_pos = max(1, Int(low_pos) - WIDE_WINDOW_CYCLE_MARGIN)
        last_pos = min(length(scans), Int(high_pos) + WIDE_WINDOW_CYCLE_MARGIN)
        for pos in first_pos:last_pos
            push!(expanded, scans[pos])
        end
    end
    sort!(expanded)
    unique!(expanded)
    return expanded
end

function _wide_scans_between_bounds(scan_index, scan_min::Integer, scan_max::Integer)
    key = scan_index.scan_to_window_key[Int(scan_min)]
    key == (Int32(0), Int32(0)) && return Int32[]
    scan_index.scan_to_window_key[Int(scan_max)] == key || return Int32[]
    low_pos = scan_index.scan_to_window_pos[Int(scan_min)]
    high_pos = scan_index.scan_to_window_pos[Int(scan_max)]
    low_pos == 0 && return Int32[]
    high_pos == 0 && return Int32[]
    if low_pos > high_pos
        low_pos, high_pos = high_pos, low_pos
    end
    scans = scan_index.window_scans[key]
    return Int32[scans[pos] for pos in Int(low_pos):Int(high_pos)]
end

function _wide_candidate_scans_from_core(tbl, group_start::Int, group_stop::Int, scan_index)
    candidate_scans = Int32[]
    @inbounds for row in group_start:group_stop
        append!(
            candidate_scans,
            _wide_scans_between_bounds(
                scan_index,
                tbl.wide_core_scan_min[row],
                tbl.wide_core_scan_max[row],
            ),
        )
    end
    if isempty(candidate_scans)
        @inbounds for row in group_start:group_stop
            push!(candidate_scans, Int32(tbl.scan_idx[row]))
        end
    end
    sort!(candidate_scans)
    unique!(candidate_scans)
    return candidate_scans
end

function _wide_top_fragment_idxs(
    frag_lookup::LibraryFragmentLookup,
    frag_list,
    nce_model,
    pid::UInt32,
    prec_charge::UInt8,
    prec_mz::Float32,
)
    frag_idxs = zeros(UInt64, 8)
    pred_intensities = fill(-Inf32, 8)
    spline_data = getSplineData(frag_lookup, nce_model, prec_charge, prec_mz)
    frag_range = getPrecFragRange(frag_lookup, pid)
    @inbounds for frag_idx in frag_range
        frag = frag_list[Int(frag_idx)]
        isIso(frag) && continue
        rank = Int(getRank(frag))
        1 <= rank <= 8 || continue
        pred_intensity = max(Float32(getIntensity(frag, spline_data)), 0f0)
        if pred_intensity > pred_intensities[rank]
            pred_intensities[rank] = pred_intensity
            frag_idxs[rank] = UInt64(frag_idx)
        end
    end
    return frag_idxs
end

function _wide_collect_feature_values(
    spectra,
    scan_index,
    expanded_scans::Vector{Int32},
    candidate_scans::Vector{Int32},
    prec_mz::Float32,
    fragment_idxs::Vector{UInt64},
    frag_list,
    ms1_ppm_tol::Float32,
    ms1_ppm_offset::Float32,
    frag_mem::AbstractMassErrorModel;
    bitvec_rank_table = nothing,
)
    n_scans = length(expanded_scans)
    n_scans == 0 && return _wide_window_zero_feature_values()

    ms1_m0 = zeros(Float32, n_scans)
    fragments = zeros(Float32, n_scans, 8)
    candidate_mask = falses(n_scans)
    ms1_target = prec_mz * (1f0 + ms1_ppm_offset * 1f-6)
    scan_corrected_mz = Float32[]
    scan_obs_low = Float32[]
    scan_obs_high = Float32[]
    candidate_pos = 1
    n_candidate_scans = length(candidate_scans)

    @inbounds for (expanded_pos, scan_i32) in pairs(expanded_scans)
        while candidate_pos <= n_candidate_scans && candidate_scans[candidate_pos] < scan_i32
            candidate_pos += 1
        end
        candidate_mask[expanded_pos] = candidate_pos <= n_candidate_scans &&
            candidate_scans[candidate_pos] == scan_i32

        scan_idx = Int(scan_i32)
        ms1_idx = scan_index.scan_to_ms1[scan_idx]
        if ms1_idx > 0
            ms1_m0[expanded_pos] = _wide_peak_intensity(
                getMzArray(spectra, Int(ms1_idx)),
                getIntensityArray(spectra, Int(ms1_idx)),
                ms1_target,
                ms1_ppm_tol,
            )
        end

        mz = getMzArray(spectra, scan_idx)
        intensities = getIntensityArray(spectra, scan_idx)
        scan_rt = Float32(getRetentionTime(spectra, scan_idx))
        n_peaks = prepare_scan_peaks!(
            scan_corrected_mz,
            scan_obs_low,
            scan_obs_high,
            frag_mem,
            mz,
            intensities,
            scan_rt,
        )

        for rank in 1:8
            frag_idx = fragment_idxs[rank]
            frag_idx > 0 || continue
            fragments[expanded_pos, rank] = _wide_fragment_mono_peak_intensity(
                scan_corrected_mz,
                scan_obs_low,
                scan_obs_high,
                intensities,
                n_peaks,
                frag_list[Int(frag_idx)],
                frag_mem,
            )
        end
    end

    return _wide_window_feature_values(
        ms1_m0,
        fragments,
        candidate_mask;
        bitvec_rank_table = bitvec_rank_table,
    )
end

function _wide_scatter_features!(columns, group_start::Int, group_stop::Int, features)
    @inbounds for row in group_start:group_stop
        columns.ms1_candidate[row] = features.wide_ms1_m0_candidate_fraction
        columns.frag_candidate[row] = features.wide_frag_candidate_fraction
        columns.ms1_frag_corr[row] = features.wide_ms1_frag_sum_corr
        columns.frag_corr_mean[row] = features.wide_frag_corr_mean
        columns.n_correlated[row] = features.wide_n_correlated_fragments
        columns.n_correlated_rank[row] = features.wide_n_correlated_fragments_bitvec_rank
        columns.corr_strength[row] = features.wide_frag_corr_strength
        columns.corr_effective_n[row] = features.wide_frag_corr_effective_n
        columns.best_m0[row] = features.wide_frag_corr_best_m0
        columns.signal_support[row] = features.wide_signal_support
    end
    return nothing
end

function add_wide_window_features_to_fold_file!(
    fpath::String,
    spectra,
    search_context::SearchContext,
    ms_file_idx::Integer,
    precursors::LibraryPrecursors,
    frag_lookup::LibraryFragmentLookup,
    nce_model,
)
    tbl = DataFrame(Tables.columntable(Arrow.Table(fpath)); copycols = true)
    n = nrow(tbl)

    tbl[!, :wide_ms1_m0_candidate_fraction] = zeros(Float32, n)
    tbl[!, :wide_frag_candidate_fraction] = zeros(Float32, n)
    tbl[!, :wide_ms1_frag_sum_corr] = zeros(Float32, n)
    tbl[!, :wide_frag_corr_mean] = zeros(Float32, n)
    tbl[!, :wide_n_correlated_fragments] = zeros(UInt8, n)
    tbl[!, :wide_n_correlated_fragments_bitvec_rank] = zeros(UInt16, n)
    tbl[!, :wide_frag_corr_strength] = zeros(Float32, n)
    tbl[!, :wide_frag_corr_effective_n] = zeros(Float32, n)
    tbl[!, :wide_frag_corr_best_m0] = zeros(Float32, n)
    tbl[!, :wide_signal_support] = zeros(Float32, n)

    if n == 0
        writeArrow(fpath, tbl)
        return n
    end

    scan_index = _wide_window_index_spectra(spectra)
    ms1_mem = getMs1MassErrorModel(search_context, ms_file_idx)
    ms1_ppm_tol = Float32(max(10f0, getLeftTol(ms1_mem), getRightTol(ms1_mem)))
    ms1_ppm_offset = Float32(getMassOffset(ms1_mem))
    frag_mem = getMassErrorModel(search_context, ms_file_idx)
    bitvec_rank_table = getBitVecExcessRanks(search_context, Int64(ms_file_idx))
    prec_mzs = getMz(precursors)
    prec_charges = getCharge(precursors)
    frag_list = getFragments(frag_lookup)
    columns = (
        ms1_candidate = tbl[!, :wide_ms1_m0_candidate_fraction],
        frag_candidate = tbl[!, :wide_frag_candidate_fraction],
        ms1_frag_corr = tbl[!, :wide_ms1_frag_sum_corr],
        frag_corr_mean = tbl[!, :wide_frag_corr_mean],
        n_correlated = tbl[!, :wide_n_correlated_fragments],
        n_correlated_rank = tbl[!, :wide_n_correlated_fragments_bitvec_rank],
        corr_strength = tbl[!, :wide_frag_corr_strength],
        corr_effective_n = tbl[!, :wide_frag_corr_effective_n],
        best_m0 = tbl[!, :wide_frag_corr_best_m0],
        signal_support = tbl[!, :wide_signal_support],
    )

    row = 1
    @inbounds while row <= n
        pid = UInt32(tbl.precursor_idx[row])
        group_start = row
        while row <= n && UInt32(tbl.precursor_idx[row]) == pid
            row += 1
        end
        group_stop = row - 1

        candidate_scans = _wide_candidate_scans_from_core(tbl, group_start, group_stop, scan_index)
        expanded_scans = _wide_expanded_window_scans(candidate_scans, scan_index)
        if isempty(expanded_scans)
            _wide_scatter_features!(columns, group_start, group_stop, _wide_window_zero_feature_values())
            continue
        end

        prec_mz = Float32(prec_mzs[pid])
        prec_charge = UInt8(prec_charges[pid])
        fragment_idxs = _wide_top_fragment_idxs(
            frag_lookup,
            frag_list,
            nce_model,
            pid,
            prec_charge,
            prec_mz,
        )
        features = _wide_collect_feature_values(
            spectra,
            scan_index,
            expanded_scans,
            candidate_scans,
            prec_mz,
            fragment_idxs,
            frag_list,
            ms1_ppm_tol,
            ms1_ppm_offset,
            frag_mem;
            bitvec_rank_table = bitvec_rank_table,
        )
        _wide_scatter_features!(columns, group_start, group_stop, features)
    end

    writeArrow(fpath, tbl)
    return n
end

function add_wide_window_features_to_fold_files!(
    search_context::SearchContext,
    valid_file_indices::Vector{Int},
)
    isempty(valid_file_indices) && return nothing

    ms_ref = getMSData(search_context)
    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    frag_lookup = getFragmentLookupTable(spec_lib)
    n_files = 0
    n_rows = 0

    for ms_file_idx in valid_file_indices
        spectra = getMSData(ms_ref, ms_file_idx)
        nce_model = getNceModel(search_context, ms_file_idx)
        for fold in UInt8[0, 1]
            fpath = getSecondPassPsmsFold(ms_ref, ms_file_idx, fold)
            (isempty(fpath) || !isfile(fpath)) && continue
            n_rows += add_wide_window_features_to_fold_file!(
                fpath,
                spectra,
                search_context,
                ms_file_idx,
                precursors,
                frag_lookup,
                nce_model,
            )
            n_files += 1
        end
    end

    @debug_l1 "Wide-window cross-run features: added $(length(WIDE_WINDOW_FEATURES)) features " *
              "to $n_files fold files ($n_rows rows)"
    return nothing
end
