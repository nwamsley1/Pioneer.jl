const WIDE_WINDOW_CYCLE_MARGIN = 3
const FLANKING_SCAN_CENTRIC_BATCH_SIZE = 50000
const FLANKING_WINDOW_TEMP_COLUMNS = [
    :flanking_core_scan_min,
    :flanking_core_scan_max,
    :flanking_core_ms1_m0_signal,
    :flanking_core_frag_signal,
]

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

@inline function _wide_candidate_signal_fraction(core_signal::Float32, flank_signal::Float32)
    total_signal = core_signal + flank_signal
    return total_signal > 0f0 ? Float32(core_signal / total_signal) : 0f0
end

@inline function _wide_window_zero_feature_values(
    core_ms1_m0_signal::Float32 = 0f0,
    core_frag_signal::Float32 = 0f0,
)
    core_ms1_m0_signal = max(core_ms1_m0_signal, 0f0)
    core_frag_signal = max(core_frag_signal, 0f0)
    return (
        flanking_ms1_m0_candidate_fraction = _wide_candidate_signal_fraction(core_ms1_m0_signal, 0f0),
        flanking_frag_candidate_fraction = _wide_candidate_signal_fraction(core_frag_signal, 0f0),
        flanking_ms1_frag_sum_corr = 0f0,
        flanking_frag_corr_mean = 0f0,
        flanking_frag_corr_strength = 0f0,
        flanking_frag_corr_effective_n = 0f0,
        flanking_frag_corr_best_m0 = 0f0,
        flanking_signal_support = 0f0,
        frag_apex_gt2x_flank_bitvec_rank = UInt16(0),
    )
end

function _wide_flank_feature_values(
    core_ms1_m0_signal::Float32,
    core_frag_signal::Float32,
    core_fragments::NTuple{8, Float32},
    flank_ms1_m0::AbstractVector,
    flank_fragments::AbstractMatrix;
    bitvec_rank_table = nothing,
)
    core_ms1_m0_signal = max(core_ms1_m0_signal, 0f0)
    core_frag_signal = max(core_frag_signal, 0f0)
    n_scans = length(flank_ms1_m0)
    n_scans == 0 && return _wide_window_zero_feature_values(core_ms1_m0_signal, core_frag_signal)

    n_frags = size(flank_fragments, 2)
    frag_sum = zeros(Float32, n_scans)
    ms1_flank_signal = 0f0
    frag_flank_signal = 0f0
    signal_scans = 0

    @inbounds for scan_i in 1:n_scans
        m0 = max(Float32(flank_ms1_m0[scan_i]), 0f0)
        summed_frag = 0f0
        for frag_i in 1:n_frags
            summed_frag += max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        frag_sum[scan_i] = summed_frag
        ms1_flank_signal += m0
        frag_flank_signal += summed_frag
        (m0 > 0f0 || summed_frag > 0f0) && (signal_scans += 1)
    end

    has_signal = falses(n_frags)
    @inbounds for frag_i in 1:n_frags
        for scan_i in 1:n_scans
            if Float32(flank_fragments[scan_i, frag_i]) > 0f0
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
            corr = _wide_window_pcor(view(flank_fragments, :, frag_a), view(flank_fragments, :, frag_b))
            pair_corr[frag_a, frag_b] = corr
            pair_corr[frag_b, frag_a] = corr
            pair_sum += corr
            pair_count += 1
        end
    end

    apex_mask = UInt16(0)
    corr_strength = 0f0
    corr_sumsq = 0f0
    other_sum = Vector{Float32}(undef, n_scans)
    @inbounds for frag_i in 1:n_frags
        flank_avg = 0f0
        for scan_i in 1:n_scans
            flank_avg += max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        flank_avg /= Float32(n_scans)
        core_apex = max(core_fragments[frag_i], 0f0)
        if core_apex > 0f0 && core_apex >= 2f0 * flank_avg
            apex_mask |= UInt16(1) << (frag_i - 1)
        end

        has_signal[frag_i] || continue
        for scan_i in 1:n_scans
            other_sum[scan_i] = frag_sum[scan_i] - max(Float32(flank_fragments[scan_i, frag_i]), 0f0)
        end
        corr = _wide_window_pcor(view(flank_fragments, :, frag_i), other_sum)
        positive_corr = min(max(corr, 0f0), 1f0)
        corr_strength += positive_corr
        corr_sumsq += positive_corr * positive_corr
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
        flanking_ms1_m0_candidate_fraction = _wide_candidate_signal_fraction(core_ms1_m0_signal, ms1_flank_signal),
        flanking_frag_candidate_fraction = _wide_candidate_signal_fraction(core_frag_signal, frag_flank_signal),
        flanking_ms1_frag_sum_corr = _wide_window_pcor(flank_ms1_m0, frag_sum),
        flanking_frag_corr_mean = pair_count > 0 ? Float32(pair_sum / pair_count) : 0f0,
        flanking_frag_corr_strength = corr_strength,
        flanking_frag_corr_effective_n = corr_sumsq > 0f0 ? Float32((corr_strength * corr_strength) / corr_sumsq) : 0f0,
        flanking_frag_corr_best_m0 = best_frag > 0 ? _wide_window_pcor(view(flank_fragments, :, best_frag), flank_ms1_m0) : 0f0,
        flanking_signal_support = Float32(signal_scans / n_scans),
        frag_apex_gt2x_flank_bitvec_rank = _bitvec_pattern_rank(bitvec_rank_table, apex_mask),
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

function _wide_append_flanking_window_scans!(
    flank_scans::Vector{Int32},
    scan_index,
    scan_min::Integer,
    scan_max::Integer,
)
    key = scan_index.scan_to_window_key[Int(scan_min)]
    key == (Int32(0), Int32(0)) && return false
    scan_index.scan_to_window_key[Int(scan_max)] == key || return false
    low_pos = scan_index.scan_to_window_pos[Int(scan_min)]
    high_pos = scan_index.scan_to_window_pos[Int(scan_max)]
    low_pos == 0 && return false
    high_pos == 0 && return false
    if low_pos > high_pos
        low_pos, high_pos = high_pos, low_pos
    end

    scans = scan_index.window_scans[key]
    for pos in max(1, Int(low_pos) - WIDE_WINDOW_CYCLE_MARGIN):(Int(low_pos) - 1)
        push!(flank_scans, scans[pos])
    end
    for pos in (Int(high_pos) + 1):min(length(scans), Int(high_pos) + WIDE_WINDOW_CYCLE_MARGIN)
        push!(flank_scans, scans[pos])
    end
    return true
end

function _wide_update_flanking_window_bounds!(
    bounds::Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}},
    scan_index,
    scan_min::Integer,
    scan_max::Integer,
)
    key = scan_index.scan_to_window_key[Int(scan_min)]
    key == (Int32(0), Int32(0)) && return false
    scan_index.scan_to_window_key[Int(scan_max)] == key || return false
    low_pos = scan_index.scan_to_window_pos[Int(scan_min)]
    high_pos = scan_index.scan_to_window_pos[Int(scan_max)]
    low_pos == 0 && return false
    high_pos == 0 && return false
    if low_pos > high_pos
        low_pos, high_pos = high_pos, low_pos
    end
    if haskey(bounds, key)
        prev_low, prev_high = bounds[key]
        bounds[key] = (min(prev_low, low_pos), max(prev_high, high_pos))
    else
        bounds[key] = (low_pos, high_pos)
    end
    return true
end

function _wide_append_flanking_window_bounds!(
    flank_scans::Vector{Int32},
    scan_index,
    bounds::Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}},
)
    for (key, (low_pos, high_pos)) in bounds
        scans = scan_index.window_scans[key]
        for pos in max(1, Int(low_pos) - WIDE_WINDOW_CYCLE_MARGIN):(Int(low_pos) - 1)
            push!(flank_scans, scans[pos])
        end
        for pos in (Int(high_pos) + 1):min(length(scans), Int(high_pos) + WIDE_WINDOW_CYCLE_MARGIN)
            push!(flank_scans, scans[pos])
        end
    end
    sort!(flank_scans)
    unique!(flank_scans)
    return flank_scans
end

function _wide_flank_window_scans_from_core!(
    flank_scans::Vector{Int32},
    scan_index,
    scan_min::Integer,
    scan_max::Integer,
)
    empty!(flank_scans)
    _wide_append_flanking_window_scans!(flank_scans, scan_index, scan_min, scan_max)
    return flank_scans
end

function _wide_group_flank_window_scans_from_core!(
    flank_scans::Vector{Int32},
    scan_index,
    core_scan_min::AbstractVector{UInt32},
    core_scan_max::AbstractVector{UInt32},
    fallback_scan_idxs::AbstractVector{UInt32},
    group_start::Int,
    group_stop::Int,
)
    empty!(flank_scans)
    bounds = Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}}()
    valid_core_bounds = false
    @inbounds for row in group_start:group_stop
        valid_core_bounds |= _wide_update_flanking_window_bounds!(
            bounds,
            scan_index,
            core_scan_min[row],
            core_scan_max[row],
        )
    end
    if !valid_core_bounds
        @inbounds for row in group_start:group_stop
            _wide_update_flanking_window_bounds!(
                bounds,
                scan_index,
                fallback_scan_idxs[row],
                fallback_scan_idxs[row],
            )
        end
    end
    return _wide_append_flanking_window_bounds!(flank_scans, scan_index, bounds)
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

struct FlankingWindowGroupBuffer
    group_start::Int
    group_stop::Int
    core_ms1_m0_signal::Float32
    core_frag_signal::Float32
    core_fragments::NTuple{8, Float32}
    ms1_target::Float32
    # Top-8 fragments precomputed once at group-build time and sorted by m/z so
    # the per-scan fill can thread a monotonic bsearch cursor (Opt 1). Stored as
    # fixed NTuples (inline, no per-group heap allocation); only entries
    # 1:n_frags are valid. Entry s maps to original intensity-rank
    # `sorted_ranks[s]` (the column to write).
    n_frags::Int
    sorted_ranks::NTuple{8, Int}
    sorted_targets::NTuple{8, Float32}
    sorted_lows::NTuple{8, Float32}
    sorted_highs::NTuple{8, Float32}
    ms1_m0::Vector{Float32}
    fragments::Matrix{Float32}
end

"""
    _wide_group_fragment_windows(frag_idxs, frag_list, mem)

From the rank-indexed `frag_idxs` (0 = absent), build the match windows for the
valid fragments sorted by ascending target m/z. `mem` supplies the (m/z-only,
scan-invariant) half-width, so `(target, low, high)` are computed once per group
instead of once per (scan x fragment). Returns `(n, ranks, targets, lows, highs)`
as fixed `NTuple{8}`s (only entries `1:n` are valid); `ranks[s]` is the original
intensity rank. `@inline` + non-escaping `MVector` stack scratch + in-place
insertion sort => no per-group heap allocation at the (inlined) call site.
"""
@inline function _wide_group_fragment_windows(
    frag_idxs::Vector{UInt64},
    frag_list,
    mem::AbstractMassErrorModel,
)
    ranks = zero(MVector{8, Int})
    targets = zero(MVector{8, Float32})
    lows = zero(MVector{8, Float32})
    highs = zero(MVector{8, Float32})
    n = 0
    @inbounds for rank in 1:8
        fi = frag_idxs[rank]
        fi > 0 || continue
        n += 1
        ranks[n] = rank
        targets[n] = Float32(getMz(frag_list[Int(fi)]))
    end
    # Insertion sort the first n entries (n <= 8) by ascending target m/z,
    # carrying the original rank alongside. Stable; allocation-free.
    @inbounds for i in 2:n
        tv = targets[i]
        rv = ranks[i]
        j = i - 1
        while j >= 1 && targets[j] > tv
            targets[j + 1] = targets[j]
            ranks[j + 1] = ranks[j]
            j -= 1
        end
        targets[j + 1] = tv
        ranks[j + 1] = rv
    end
    @inbounds for s in 1:n
        hw = conservative_half_width(mem, targets[s])
        lo, hi = match_window(targets[s], hw)
        lows[s] = lo
        highs[s] = hi
    end
    return n, Tuple(ranks), Tuple(targets), Tuple(lows), Tuple(highs)
end

function _wide_add_group_ms2_work!(
    ms2_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    group_idx::Int,
    flank_scans::Vector{Int32},
)
    @inbounds for (flank_pos, scan_i32) in pairs(flank_scans)
        work_items = get!(() -> Tuple{Int, Int}[], ms2_work, scan_i32)
        push!(work_items, (group_idx, flank_pos))
    end
    return ms2_work
end

function _wide_add_group_ms1_work!(
    ms1_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    scan_index,
    group_idx::Int,
    flank_scans::Vector{Int32},
)
    @inbounds for (flank_pos, scan_i32) in pairs(flank_scans)
        ms1_i32 = scan_index.scan_to_ms1[Int(scan_i32)]
        ms1_i32 > 0 || continue
        work_items = get!(() -> Tuple{Int, Int}[], ms1_work, ms1_i32)
        push!(work_items, (group_idx, flank_pos))
    end
    return ms1_work
end

function _wide_fill_scan_centric_ms1_m0!(
    groups::Vector{FlankingWindowGroupBuffer},
    ms1_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    spectra,
    ms1_ppm_tol::Float32,
)
    ms1_keys = collect(keys(ms1_work))
    Threads.@threads for key_idx in eachindex(ms1_keys)
        ms1_i32 = ms1_keys[key_idx]
        work_items = ms1_work[ms1_i32]
        mz = getMzArray(spectra, Int(ms1_i32))
        intensities = getIntensityArray(spectra, Int(ms1_i32))

        @inbounds for (group_idx, flank_pos) in work_items
            group = groups[group_idx]
            group.ms1_m0[flank_pos] = _wide_peak_intensity(
                mz,
                intensities,
                group.ms1_target,
                ms1_ppm_tol,
            )
        end
    end
    return nothing
end

function _wide_fill_scan_centric_fragments!(
    groups::Vector{FlankingWindowGroupBuffer},
    ms2_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    spectra,
    frag_list,
    frag_mem::AbstractMassErrorModel,
)
    scan_keys = collect(keys(ms2_work))
    max_thread_id = Threads.maxthreadid()
    scan_corrected_mz = [Float32[] for _ in 1:max_thread_id]
    scan_obs_low = [Float32[] for _ in 1:max_thread_id]
    scan_obs_high = [Float32[] for _ in 1:max_thread_id]

    Threads.@threads for key_idx in eachindex(scan_keys)
        scan_i32 = scan_keys[key_idx]
        work_items = ms2_work[scan_i32]
        scan_idx = Int(scan_i32)
        mz = getMzArray(spectra, scan_idx)
        intensities = getIntensityArray(spectra, scan_idx)
        scan_rt = Float32(getRetentionTime(spectra, scan_idx))
        thread_idx = Threads.threadid()
        n_peaks = prepare_scan_peaks!(
            scan_corrected_mz[thread_idx],
            scan_obs_low[thread_idx],
            scan_obs_high[thread_idx],
            frag_mem,
            mz,
            intensities,
            scan_rt,
        )

        @inbounds for (group_idx, flank_pos) in work_items
            group = groups[group_idx]
            n_frags = group.n_frags
            sorted_ranks = group.sorted_ranks
            sorted_targets = group.sorted_targets
            sorted_lows = group.sorted_lows
            sorted_highs = group.sorted_highs
            fragments = group.fragments
            # Fragments are sorted by ascending m/z, so each successive bsearch
            # can start from the previous fragment's first-peak index: the search
            # range only ever shrinks, and the result is identical to searching
            # [1, n_peaks] every time (low_{s+1} >= low_s).
            start_idx = 1
            for s in 1:n_frags
                start_idx = bsearch_hybrid(
                    scan_corrected_mz[thread_idx], sorted_lows[s], start_idx, n_peaks,
                )
                best_peak, _, _ = scan_for_nearest_in_window(
                    scan_corrected_mz[thread_idx],
                    scan_obs_low[thread_idx],
                    scan_obs_high[thread_idx],
                    start_idx,
                    n_peaks,
                    sorted_targets[s],
                    sorted_highs[s],
                )
                if best_peak != 0
                    intensity = intensities[best_peak]
                    fragments[flank_pos, sorted_ranks[s]] =
                        ismissing(intensity) ? 0f0 : Float32(intensity)
                end
            end
        end
    end
    return nothing
end

function _wide_reduce_scan_centric_group_batch!(columns, groups, bitvec_rank_table)
    Threads.@threads for group_idx in eachindex(groups)
        group = groups[group_idx]
        features = _wide_flank_feature_values(
            group.core_ms1_m0_signal,
            group.core_frag_signal,
            group.core_fragments,
            group.ms1_m0,
            group.fragments;
            bitvec_rank_table = bitvec_rank_table,
        )
        _wide_scatter_features!(columns, group.group_start, group.group_stop, features)
    end
    return nothing
end

function _wide_flush_scan_centric_group_batch!(
    columns,
    groups::Vector{FlankingWindowGroupBuffer},
    ms1_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    ms2_work::Dict{Int32, Vector{Tuple{Int, Int}}},
    spectra,
    frag_list,
    frag_mem::AbstractMassErrorModel,
    ms1_ppm_tol::Float32,
    bitvec_rank_table,
)
    _wide_fill_scan_centric_ms1_m0!(
        groups,
        ms1_work,
        spectra,
        ms1_ppm_tol,
    )
    _wide_fill_scan_centric_fragments!(
        groups,
        ms2_work,
        spectra,
        frag_list,
        frag_mem,
    )

    _wide_reduce_scan_centric_group_batch!(columns, groups, bitvec_rank_table)
    empty!(groups)
    empty!(ms1_work)
    empty!(ms2_work)
    return nothing
end

function _wide_scatter_features!(columns, group_start::Int, group_stop::Int, features)
    @inbounds for row in group_start:group_stop
        columns.ms1_candidate[row] = features.flanking_ms1_m0_candidate_fraction
        columns.frag_candidate[row] = features.flanking_frag_candidate_fraction
        columns.ms1_frag_corr[row] = features.flanking_ms1_frag_sum_corr
        columns.frag_corr_mean[row] = features.flanking_frag_corr_mean
        columns.corr_strength[row] = features.flanking_frag_corr_strength
        columns.corr_effective_n[row] = features.flanking_frag_corr_effective_n
        columns.best_m0[row] = features.flanking_frag_corr_best_m0
        columns.signal_support[row] = features.flanking_signal_support
        columns.apex_gt2x_rank[row] = features.frag_apex_gt2x_flank_bitvec_rank
    end
    return nothing
end

function add_wide_window_features_to_table!(
    tbl::DataFrame,
    spectra,
    search_context,
    ms_file_idx::Integer,
    precursors,
    frag_lookup,
    nce_model,
)
    n = nrow(tbl)

    flanking_ms1_m0_candidate_fraction = zeros(Float32, n)
    flanking_frag_candidate_fraction = zeros(Float32, n)
    flanking_ms1_frag_sum_corr = zeros(Float32, n)
    flanking_frag_corr_mean = zeros(Float32, n)
    flanking_frag_corr_strength = zeros(Float32, n)
    flanking_frag_corr_effective_n = zeros(Float32, n)
    flanking_frag_corr_best_m0 = zeros(Float32, n)
    flanking_signal_support = zeros(Float32, n)
    frag_apex_gt2x_flank_bitvec_rank = zeros(UInt16, n)

    tbl[!, :flanking_ms1_m0_candidate_fraction] = flanking_ms1_m0_candidate_fraction
    tbl[!, :flanking_frag_candidate_fraction] = flanking_frag_candidate_fraction
    tbl[!, :flanking_ms1_frag_sum_corr] = flanking_ms1_frag_sum_corr
    tbl[!, :flanking_frag_corr_mean] = flanking_frag_corr_mean
    tbl[!, :flanking_frag_corr_strength] = flanking_frag_corr_strength
    tbl[!, :flanking_frag_corr_effective_n] = flanking_frag_corr_effective_n
    tbl[!, :flanking_frag_corr_best_m0] = flanking_frag_corr_best_m0
    tbl[!, :flanking_signal_support] = flanking_signal_support
    tbl[!, :frag_apex_gt2x_flank_bitvec_rank] = frag_apex_gt2x_flank_bitvec_rank

    if n == 0
        select!(tbl, Not(FLANKING_WINDOW_TEMP_COLUMNS))
        return n
    end

    precursor_idxs = tbl[!, :precursor_idx]::Vector{UInt32}
    scan_idxs = tbl[!, :scan_idx]::Vector{UInt32}
    flanking_core_scan_min = tbl[!, :flanking_core_scan_min]::Vector{UInt32}
    flanking_core_scan_max = tbl[!, :flanking_core_scan_max]::Vector{UInt32}
    flanking_core_ms1_m0_signal = tbl[!, :flanking_core_ms1_m0_signal]::Vector{Float32}
    flanking_core_frag_signal = tbl[!, :flanking_core_frag_signal]::Vector{Float32}
    core_fragment_cols = Tuple(tbl[!, col]::Vector{Float32} for col in SMOOTHED_FRAGMENT_INTENSITY_COLUMNS)
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
        ms1_candidate = flanking_ms1_m0_candidate_fraction,
        frag_candidate = flanking_frag_candidate_fraction,
        ms1_frag_corr = flanking_ms1_frag_sum_corr,
        frag_corr_mean = flanking_frag_corr_mean,
        corr_strength = flanking_frag_corr_strength,
        corr_effective_n = flanking_frag_corr_effective_n,
        best_m0 = flanking_frag_corr_best_m0,
        signal_support = flanking_signal_support,
        apex_gt2x_rank = frag_apex_gt2x_flank_bitvec_rank,
    )
    flank_scans = Int32[]
    groups = FlankingWindowGroupBuffer[]
    ms1_work = Dict{Int32, Vector{Tuple{Int, Int}}}()
    ms2_work = Dict{Int32, Vector{Tuple{Int, Int}}}()

    row = 1
    @inbounds while row <= n
        pid = precursor_idxs[row]
        group_start = row
        while row <= n && precursor_idxs[row] == pid
            row += 1
        end
        group_stop = row - 1

        _wide_group_flank_window_scans_from_core!(
            flank_scans,
            scan_index,
            flanking_core_scan_min,
            flanking_core_scan_max,
            scan_idxs,
            group_start,
            group_stop,
        )

        core_ms1_m0_signal = flanking_core_ms1_m0_signal[group_start]
        core_frag_signal = flanking_core_frag_signal[group_start]
        if isempty(flank_scans)
            _wide_scatter_features!(
                columns,
                group_start,
                group_stop,
                _wide_window_zero_feature_values(core_ms1_m0_signal, core_frag_signal),
            )
            continue
        end

        prec_mz = Float32(prec_mzs[pid])
        prec_charge = UInt8(prec_charges[pid])
        core_fragments = ntuple(rank -> Float32(core_fragment_cols[rank][group_start]), 8)
        fragment_idxs = _wide_top_fragment_idxs(
            frag_lookup,
            frag_list,
            nce_model,
            pid,
            prec_charge,
            prec_mz,
        )
        n_frags, sorted_ranks, sorted_targets, sorted_lows, sorted_highs =
            _wide_group_fragment_windows(fragment_idxs, frag_list, frag_mem)
        ms1_m0 = zeros(Float32, length(flank_scans))
        ms1_target = prec_mz * (1f0 + ms1_ppm_offset * 1f-6)
        fragments = zeros(Float32, length(flank_scans), 8)
        push!(
            groups,
            FlankingWindowGroupBuffer(
                group_start,
                group_stop,
                core_ms1_m0_signal,
                core_frag_signal,
                core_fragments,
                ms1_target,
                n_frags,
                sorted_ranks,
                sorted_targets,
                sorted_lows,
                sorted_highs,
                ms1_m0,
                fragments,
            ),
        )
        _wide_add_group_ms1_work!(ms1_work, scan_index, length(groups), flank_scans)
        _wide_add_group_ms2_work!(ms2_work, length(groups), flank_scans)

        if length(groups) >= FLANKING_SCAN_CENTRIC_BATCH_SIZE
            _wide_flush_scan_centric_group_batch!(
                columns,
                groups,
                ms1_work,
                ms2_work,
                spectra,
                frag_list,
                frag_mem,
                ms1_ppm_tol,
                bitvec_rank_table,
            )
        end
    end

    _wide_flush_scan_centric_group_batch!(
        columns,
        groups,
        ms1_work,
        ms2_work,
        spectra,
        frag_list,
        frag_mem,
        ms1_ppm_tol,
        bitvec_rank_table,
    )

    select!(tbl, Not(FLANKING_WINDOW_TEMP_COLUMNS))
    return n
end
