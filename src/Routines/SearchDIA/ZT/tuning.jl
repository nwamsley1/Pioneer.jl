# Scanning-quad (ZT) hooks in ParameterTuningSearch.

"""Q1 bins per RT bin in the scan priority order for scanning-quad files."""
const TUNING_ZT_Q1_BINS = 15

"""
    tuning_scan_priority(geom, spectra) -> Vector{Int32}

MS2 scan order for parameter tuning: develop's TIC-first order, or on a scanning-quad file the
Q1-stratified order (`get_ms2_scan_priority_order_q1`).
"""
tuning_scan_priority(::Nothing, spectra::MassSpecData) = get_ms2_scan_priority_order(spectra)
tuning_scan_priority(::ZTGeometry, spectra::MassSpecData) =
    get_ms2_scan_priority_order_q1(spectra, TUNING_ZT_Q1_BINS)

"""
    zt_keeps_quad_model(geom) -> Bool

Whether parameter tuning must keep the file's installed transmission model rather than reset it
to a square box of the stated isolation width. On a scanning quad the stated width is only the
Q1 step, and the model installed from the geometry spans the whole meta-scan; resetting it would
leave MainSearch deconvolving with a ~±0.5 Da box, giving every expanded neighbour-bin candidate
zero transmission.
"""
zt_keeps_quad_model(::Nothing) = false
zt_keeps_quad_model(::ZTGeometry) = true

"""
    get_ms2_scan_priority_order_q1(spectra::MassSpecData, n_q1_bins::Int) -> Vector{Int32}

Scanning-quad variant of `get_ms2_scan_priority_order`. Same outer round-robin over
`SCAN_PRIORITY_N_RT_BINS` RT bins, but each RT bin keeps its own rotation over `n_q1_bins` Q1
(window-centre m/z) bins and a visit returns the highest-TIC undrawn scan of the NEXT Q1 bin in
that rotation. No Q1 bin is revisited until every other Q1 bin of that RT bin has given a scan.

Why: on a swept quad the highest-TIC 1-Da windows of successive cycles are the same windows,
the ones holding the dominant precursor, so TIC-first within an RT bin resamples a handful of
precursors (measured on a 54-min nano file: 3 peptides supplied 613 of 1,311 tuning PSMs).
Exactly sized throughout; no push!.
"""
function get_ms2_scan_priority_order_q1(spectra::MassSpecData, n_q1_bins::Int)
    n_rt = SCAN_PRIORITY_N_RT_BINS
    rt_bins, _, _, _ = compute_rt_bins(spectra, n_rt)
    ms_orders = getMsOrders(spectra)
    cmz = getCenterMzs(spectra)
    tics = getTICs(spectra)
    ms2 = Int32[i for i in 1:length(spectra) if ms_orders[i] == 2 && !ismissing(cmz[i])]
    n = length(ms2)
    n == 0 && return Int32[]
    mz_lo = minimum(Float32(cmz[i]) for i in ms2); mz_hi = maximum(Float32(cmz[i]) for i in ms2)
    mz_w = max(mz_hi - mz_lo, eps(Float32)) / n_q1_bins
    q1_of(i) = clamp(ceil(Int, (Float32(cmz[i]) - mz_lo) / mz_w), 1, n_q1_bins)
    cell(i) = (Int(rt_bins[i]) - 1) * n_q1_bins + q1_of(i)      # 1-based, RT-major
    ncell = n_rt * n_q1_bins
    # counting sort into cells
    counts = zeros(Int32, ncell)
    @inbounds for i in ms2; counts[cell(i)] += 1; end
    starts = Vector{Int32}(undef, ncell + 1); starts[1] = 1
    @inbounds for c in 1:ncell; starts[c + 1] = starts[c] + counts[c]; end
    sorted = Vector{Int32}(undef, n)
    fill!(counts, 0)
    @inbounds for i in ms2
        c = cell(i); sorted[starts[c] + counts[c]] = i; counts[c] += 1
    end
    @inbounds for c in 1:ncell
        lo, hi = starts[c], starts[c + 1] - 1
        lo < hi && sort!(view(sorted, lo:hi), by = idx -> tics[idx], rev = true)
    end
    # emit: outer pass over RT bins, per-RT-bin Q1 cursor
    pos = copy(starts)                       # next undrawn scan per cell
    q1_cursor = ones(Int, n_rt)
    out = Vector{Int32}(undef, n)
    w = 1
    while w <= n
        for r in 1:n_rt
            base = (r - 1) * n_q1_bins
            for _ in 1:n_q1_bins                     # find the next non-empty Q1 bin of this RT bin
                q = q1_cursor[r]
                q1_cursor[r] = q == n_q1_bins ? 1 : q + 1
                c = base + q
                if pos[c] < starts[c + 1]
                    out[w] = sorted[pos[c]]; pos[c] += 1; w += 1
                    break
                end
            end
            w > n && break
        end
    end
    return out
end
