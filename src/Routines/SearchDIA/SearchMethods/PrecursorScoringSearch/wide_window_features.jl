# Raw wide-window chromatogram features for the experiment-wide
# PrecursorScoringSearch LightGBM.
#
# Existing fragment/MS1 chromatogram features are computed only over PSM rows
# that survived the fragment-index search. These features inspect the raw
# same-window MS2 scans immediately outside that candidate range, plus nearest
# MS1 scans, so the model can learn whether signal truly disappears beyond the
# fragment-index support.

const WIDE_WINDOW_CYCLE_MARGIN = 6
const WIDE_WINDOW_FRAGMENT_PPM_TOL = 20f0
const WIDE_WINDOW_CORR_THRESHOLD = 0.7f0

@inline function _wide_window_pcor(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(ArgumentError("correlation vectors must have equal length"))
    n < 2 && return 0f0
    mx = 0f0
    my = 0f0
    @inbounds for i in 1:n
        mx += Float32(x[i])
        my += Float32(y[i])
    end
    inv_n = 1f0 / Float32(n)
    mx *= inv_n
    my *= inv_n
    sx = 0f0
    sy = 0f0
    sxy = 0f0
    @inbounds for i in 1:n
        dx = Float32(x[i]) - mx
        dy = Float32(y[i]) - my
        sx += dx * dx
        sy += dy * dy
        sxy += dx * dy
    end
    denom = sqrt(sx * sy)
    return denom > 0f0 ? Float32(sxy / denom) : 0f0
end

@inline function _wide_window_zero_feature_values()
    return (
        wide_ms1_m0_candidate_fraction = 0f0,
        wide_frag_candidate_fraction = 0f0,
        wide_ms1_frag_sum_corr = 0f0,
        wide_frag_corr_mean = 0f0,
        wide_n_correlated_fragments = UInt8(0),
        wide_frag_corr_best_m0 = 0f0,
        wide_signal_support = 0f0,
    )
end

function _wide_window_feature_values(
    ms1_m0::AbstractVector,
    fragments::AbstractMatrix,
    candidate_mask::AbstractVector{Bool},
)
    n_scans = length(ms1_m0)
    n_scans == length(candidate_mask) ||
        throw(ArgumentError("candidate_mask length must match ms1_m0"))
    size(fragments, 1) == n_scans ||
        throw(ArgumentError("fragment matrix row count must match ms1_m0"))
    n_scans == 0 && return _wide_window_zero_feature_values()

    n_frags = size(fragments, 2)
    frag_sum = zeros(Float32, n_scans)

    ms1_total = 0f0
    ms1_candidate = 0f0
    frag_total = 0f0
    frag_candidate = 0f0
    signal_scans = 0

    @inbounds for i in 1:n_scans
        m0 = max(Float32(ms1_m0[i]), 0f0)
        fs = 0f0
        for r in 1:n_frags
            fs += max(Float32(fragments[i, r]), 0f0)
        end
        frag_sum[i] = fs
        ms1_total += m0
        frag_total += fs
        if candidate_mask[i]
            ms1_candidate += m0
            frag_candidate += fs
        end
        (m0 > 0f0 || fs > 0f0) && (signal_scans += 1)
    end

    ms1_candidate_fraction = ms1_total > 0f0 ? Float32(ms1_candidate / ms1_total) : 0f0
    frag_candidate_fraction = frag_total > 0f0 ? Float32(frag_candidate / frag_total) : 0f0
    ms1_frag_sum_corr = _wide_window_pcor(ms1_m0, frag_sum)
    signal_support = Float32(signal_scans / n_scans)

    has_signal = falses(n_frags)
    @inbounds for r in 1:n_frags
        for i in 1:n_scans
            if Float32(fragments[i, r]) > 0f0
                has_signal[r] = true
                break
            end
        end
    end

    pair_corr = zeros(Float32, n_frags, n_frags)
    pair_sum = 0f0
    pair_count = 0
    @inbounds for r1 in 1:n_frags
        has_signal[r1] || continue
        for r2 in (r1 + 1):n_frags
            has_signal[r2] || continue
            c = _wide_window_pcor(view(fragments, :, r1), view(fragments, :, r2))
            pair_corr[r1, r2] = c
            pair_corr[r2, r1] = c
            pair_sum += c
            pair_count += 1
        end
    end
    frag_corr_mean = pair_count > 0 ? Float32(pair_sum / pair_count) : 0f0

    n_correlated = UInt8(0)
    other_sum = Vector{Float32}(undef, n_scans)
    @inbounds for r in 1:n_frags
        has_signal[r] || continue
        for i in 1:n_scans
            other_sum[i] = frag_sum[i] - max(Float32(fragments[i, r]), 0f0)
        end
        _wide_window_pcor(view(fragments, :, r), other_sum) > WIDE_WINDOW_CORR_THRESHOLD &&
            (n_correlated += UInt8(1))
    end

    best_r = 0
    best_consensus = typemin(Float32)
    @inbounds for r in 1:n_frags
        has_signal[r] || continue
        consensus = 0f0
        npairs = 0
        for r2 in 1:n_frags
            (r2 == r || !has_signal[r2]) && continue
            consensus += pair_corr[r, r2]
            npairs += 1
        end
        avg = npairs > 0 ? Float32(consensus / npairs) : typemin(Float32)
        if avg > best_consensus
            best_consensus = avg
            best_r = r
        end
    end
    best_m0_corr = best_r > 0 ? _wide_window_pcor(view(fragments, :, best_r), ms1_m0) : 0f0

    return (
        wide_ms1_m0_candidate_fraction = ms1_candidate_fraction,
        wide_frag_candidate_fraction = frag_candidate_fraction,
        wide_ms1_frag_sum_corr = ms1_frag_sum_corr,
        wide_frag_corr_mean = frag_corr_mean,
        wide_n_correlated_fragments = n_correlated,
        wide_frag_corr_best_m0 = best_m0_corr,
        wide_signal_support = signal_support,
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
    lo_mz = target - tol
    hi_mz = target + tol
    i0 = bsearch_hybrid(mz, lo_mz, 1, n_peaks)
    best_diff = Inf32
    best_int = 0f0
    @inbounds for j in i0:n_peaks
        m = mz[j]
        m > hi_mz && break
        val = intens[j]
        ismissing(val) && continue
        diff = abs(m - target)
        if diff < best_diff
            best_diff = diff
            best_int = Float32(val)
        end
    end
    return best_int
end

@inline function _wide_peak_intensity(mz, intens, target::Float32, ppm_tol::Float32)
    n_peaks = length(mz)
    n_peaks == 0 && return 0f0
    tol = target * ppm_tol * 1f-6
    lo_mz = target - tol
    hi_mz = target + tol
    best_diff = Inf32
    best_int = 0f0
    @inbounds for j in 1:n_peaks
        m = mz[j]
        ismissing(m) && continue
        mf = Float32(m)
        mf < lo_mz && continue
        mf > hi_mz && break
        val = intens[j]
        ismissing(val) && continue
        diff = abs(mf - target)
        if diff < best_diff
            best_diff = diff
            best_int = Float32(val)
        end
    end
    return best_int
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
        order = getMsOrder(spectra, scan_idx)
        if order == UInt8(1)
            push!(ms1_scan_idxs, Int32(scan_idx))
            push!(ms1_scan_rts, Float32(getRetentionTime(spectra, scan_idx)))
        elseif order == UInt8(2)
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
                d_after = abs(ms1_scan_rts[pos] - scan_rt)
                d_before = abs(ms1_scan_rts[pos - 1] - scan_rt)
                d_before <= d_after ? ms1_scan_idxs[pos - 1] : ms1_scan_idxs[pos]
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
    n_scans = length(scan_index.scan_to_window_pos)
    @inbounds for scan in candidate_scans
        s = Int(scan)
        (s < 1 || s > n_scans) && continue
        pos = scan_index.scan_to_window_pos[s]
        pos == 0 && continue
        key = scan_index.scan_to_window_key[s]
        if haskey(bounds, key)
            lo, hi = bounds[key]
            bounds[key] = (min(lo, pos), max(hi, pos))
        else
            bounds[key] = (pos, pos)
        end
    end

    expanded = Int32[]
    for (key, (lo, hi)) in bounds
        scans = scan_index.window_scans[key]
        first_pos = max(1, Int(lo) - WIDE_WINDOW_CYCLE_MARGIN)
        last_pos = min(length(scans), Int(hi) + WIDE_WINDOW_CYCLE_MARGIN)
        for pos in first_pos:last_pos
            push!(expanded, scans[pos])
        end
    end
    sort!(expanded)
    unique!(expanded)
    return expanded
end

function _wide_top_fragment_mzs(
    frag_lookup::LibraryFragmentLookup,
    frag_list,
    nce_model,
    pid::UInt32,
    prec_charge::UInt8,
    prec_mz::Float32,
)
    mzs = zeros(Float32, 6)
    pred_ints = fill(-Inf32, 6)
    spline_data = getSplineData(frag_lookup, nce_model, prec_charge, prec_mz)
    frag_range = getPrecFragRange(frag_lookup, pid)
    @inbounds for frag_idx in frag_range
        frag = frag_list[Int(frag_idx)]
        isIso(frag) && continue
        rank = Int(getRank(frag))
        1 <= rank <= 6 || continue
        pred_int = max(Float32(getIntensity(frag, spline_data)), 0f0)
        if pred_int > pred_ints[rank]
            pred_ints[rank] = pred_int
            mzs[rank] = Float32(getMz(frag))
        end
    end
    return mzs
end

function _wide_collect_feature_values(
    spectra,
    scan_index,
    expanded_scans::Vector{Int32},
    candidate_scans::Vector{Int32},
    prec_mz::Float32,
    fragment_mzs::Vector{Float32},
    ms1_ppm_tol::Float32,
    ms1_ppm_offset::Float32,
)
    n = length(expanded_scans)
    n == 0 && return _wide_window_zero_feature_values()

    candidate_set = Set(candidate_scans)
    ms1_m0 = zeros(Float32, n)
    fragments = zeros(Float32, n, 6)
    candidate_mask = falses(n)
    ms1_target = prec_mz * (1f0 + ms1_ppm_offset * 1f-6)

    @inbounds for (i, scan_i32) in pairs(expanded_scans)
        scan_idx = Int(scan_i32)
        candidate_mask[i] = scan_i32 in candidate_set

        ms1_idx = scan_index.scan_to_ms1[scan_idx]
        if ms1_idx > 0
            ms1_m0[i] = _wide_peak_intensity(
                getMzArray(spectra, Int(ms1_idx)),
                getIntensityArray(spectra, Int(ms1_idx)),
                ms1_target,
                ms1_ppm_tol,
            )
        end

        mz = getMzArray(spectra, scan_idx)
        intens = getIntensityArray(spectra, scan_idx)
        for r in 1:6
            target = fragment_mzs[r]
            target > 0f0 || continue
            fragments[i, r] = _wide_peak_intensity(mz, intens, target, WIDE_WINDOW_FRAGMENT_PPM_TOL)
        end
    end

    return _wide_window_feature_values(ms1_m0, fragments, candidate_mask)
end

function _wide_precursor_groups(precursor_idx)
    perm = sortperm(precursor_idx)
    starts = Int[]
    ends = Int[]
    isempty(perm) && return perm, starts, ends
    start_i = 1
    prev = UInt32(precursor_idx[perm[1]])
    @inbounds for i in 2:length(perm)
        pid = UInt32(precursor_idx[perm[i]])
        if pid != prev
            push!(starts, start_i)
            push!(ends, i - 1)
            start_i = i
            prev = pid
        end
    end
    push!(starts, start_i)
    push!(ends, length(perm))
    return perm, starts, ends
end

function _wide_scatter_features!(columns, rows::Vector{Int}, features)
    @inbounds for row in rows
        columns.ms1_candidate[row] = features.wide_ms1_m0_candidate_fraction
        columns.frag_candidate[row] = features.wide_frag_candidate_fraction
        columns.ms1_frag_corr[row] = features.wide_ms1_frag_sum_corr
        columns.frag_corr_mean[row] = features.wide_frag_corr_mean
        columns.n_correlated[row] = features.wide_n_correlated_fragments
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
    tbl[!, :wide_frag_corr_best_m0] = zeros(Float32, n)
    tbl[!, :wide_signal_support] = zeros(Float32, n)

    if n == 0 || !all(c -> hasproperty(tbl, c), (:precursor_idx, :scan_idx))
        writeArrow(fpath, tbl)
        return n
    end

    scan_index = _wide_window_index_spectra(spectra)
    if isempty(scan_index.window_scans)
        writeArrow(fpath, tbl)
        return n
    end

    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    ms1_ppm_tol = Float32(max(10f0, getLeftTol(mem), getRightTol(mem)))
    ms1_ppm_offset = Float32(getMassOffset(mem))
    prec_mzs = getMz(precursors)
    prec_charges = getCharge(precursors)
    frag_list = getFragments(frag_lookup)

    perm, starts, ends = _wide_precursor_groups(tbl.precursor_idx)
    columns = (
        ms1_candidate = tbl[!, :wide_ms1_m0_candidate_fraction],
        frag_candidate = tbl[!, :wide_frag_candidate_fraction],
        ms1_frag_corr = tbl[!, :wide_ms1_frag_sum_corr],
        frag_corr_mean = tbl[!, :wide_frag_corr_mean],
        n_correlated = tbl[!, :wide_n_correlated_fragments],
        best_m0 = tbl[!, :wide_frag_corr_best_m0],
        signal_support = tbl[!, :wide_signal_support],
    )

    n_precursors = length(starts)
    Threads.@threads :static for group_idx in 1:n_precursors
        i_start = starts[group_idx]
        i_end = ends[group_idx]
        n_rows = i_end - i_start + 1
        rows = Vector{Int}(undef, n_rows)
        candidate_scans = Vector{Int32}(undef, n_rows)
        @inbounds for k in 1:n_rows
            row = perm[i_start + k - 1]
            rows[k] = row
            candidate_scans[k] = Int32(tbl.scan_idx[row])
        end

        pid = UInt32(tbl.precursor_idx[rows[1]])
        if Int(pid) > length(prec_mzs)
            _wide_scatter_features!(columns, rows, _wide_window_zero_feature_values())
            continue
        end

        expanded_scans = _wide_expanded_window_scans(candidate_scans, scan_index)
        if isempty(expanded_scans)
            _wide_scatter_features!(columns, rows, _wide_window_zero_feature_values())
            continue
        end

        prec_mz = Float32(prec_mzs[pid])
        prec_charge = UInt8(prec_charges[pid])
        fragment_mzs = _wide_top_fragment_mzs(
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
            fragment_mzs,
            ms1_ppm_tol,
            ms1_ppm_offset,
        )
        _wide_scatter_features!(columns, rows, features)
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
    t = @elapsed begin
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
    end

    @debug_l1 "Wide-window cross-run features: added $(length(WIDE_WINDOW_FEATURES)) features " *
              "to $n_files fold files ($n_rows rows) in $(round(t, digits = 2))s"
    return nothing
end
