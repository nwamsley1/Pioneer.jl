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
    ZTGeometry

Q1 bin lattice for a scanning ("swept") quadrupole DIA acquisition.

The recorded `isolationWidthMz` is only the Q1 *step*: the physical quadrupole window is several
Da wide and swept across m/z, and each recorded bin accumulates while the quad advances one step.
One precursor's fragments therefore appear in many adjacent bins with a smooth weight profile,
which is what `metascan_k` spans.

Fields:
- `bin_step`      median adjacent-`centerMz` difference within a cycle (Da). The true Q1 step.
- `nominal_width` median recorded `isolationWidthMz` (Da). Used only for the tiling check.
- `bins_per_ramp` median MS2 scans per centerMz ramp (a "cycle" per `compute_cycle_idxs`;
                  MS1 scans are not delimiters, they inherit the current cycle index).
- `metascan_k`    expansion half-width in bins. A config value, identical for every file.
- `transmission_fwhm` FWHM (Da) of the effective transmission profile. A PHYSICAL property of
                  the acquisition, deliberately independent of `metascan_k`: `k` decides how much
                  of the profile we sample, not what shape it has.

`bin_step` is deliberately taken from `centerMz` differences rather than `isolationWidthMz`:
the recorded width is dithered 50/50 between two quantized values (e.g. 1.0200 / 1.0241) whose
*average* is the true step, so a median over that column is degenerate and lands ~2e-3 Da low.
"""
struct ZTGeometry
    bin_step::Float32
    nominal_width::Float32
    bins_per_ramp::Int32
    metascan_k::Int32
    transmission_fwhm::Float32
    # Fitted transmission half-base (Da) from QuadTuningSearch; 0 until fitted. When > 0 the
    # collapse template is a triangle of this width instead of the Gaussian above.
    template_h::Float32
end

"""Expansion half-width used before quad tuning when `acquisition.metascan_k` is absent."""
const ZT_METASCAN_K_DEFAULT = 6

"""
    zt_metascan_k_is_derived(params) -> Bool

True when the config leaves `acquisition.metascan_k` unset (or 0), so QuadTuningSearch may
replace the provisional default with the value implied by the fitted transmission profile.
"""
function zt_metascan_k_is_derived(params)
    acq = params.acquisition
    return !(hasproperty(acq, :metascan_k) && Int(acq.metascan_k) > 0)
end

"""
Fraction of the transmission half-base `h` used by the per-meta-scan triangle regression: the
fit sees only |Δ| <= `ZT_FIT_CORE_FRACTION * h` and extrapolates to the zero crossing. The
empirical medians sit above a triangle in the core and below it in the tails (the profile is
closer to Gaussian), so the tails would drag the slope; but too small a core sees only the
convex apex and the slope comes out too shallow. Expressed as a fraction of h, not in Da or
bins, so every method fits the same part of its own profile.

Measured (A_REP1, both methods, k configured = 6):
- 5 Da:  0.62 h (= the validated 4.0 Da window) -> h 6.48, k_implied 6, R² 0.92;
         0.50 h -> h 6.65, k_implied 7, R² 0.89, -78 precursors.
- 10 Da: 0.31 h (a fixed 4 Da) -> h 13.59, k_implied 7, R² 0.905, IQR 3.3 Da;
         0.50 h -> h 12.75, k_implied 6, R² 0.937.
0.62 keeps ~4 bins per side on both methods. Before the fit exists, `h` is taken as the
expansion span `metascan_k * bin_step` (the geometry's own estimate).
"""
const ZT_FIT_CORE_FRACTION = 0.62f0

"""Fit half-width in Da before any fit exists: half the expansion span."""
zt_fit_limit_da(g::ZTGeometry) = ZT_FIT_CORE_FRACTION * Float32(g.metascan_k) * g.bin_step

"""Fit half-width in Da once `h` is known."""
zt_fit_limit_da(h::Real) = ZT_FIT_CORE_FRACTION * Float32(h)

"""
Target candidates (precursor x scan pairs after expansion) per main-search chunk on a
scanning-quad file. Exact, not predicted: the fragment index + expansion run once over the whole
file before chunking. Rows are 0.13-0.26 of candidates on the ZT files seen, so 60M candidates
is ~10-15M raw rows per chunk. nano15 (634M candidates, 164M rows in one pass) swapped a 48 GB
machine. `PIONEER_ZT_CHUNK_CANDIDATES` overrides it for tuning.
"""
const ZT_CHUNK_CANDIDATES_DEFAULT = 30_000_000   # 60M -> 30M (2026-09-17): chunk plateau 26.9 -> 22.5 GB on EV1109, no time cost
zt_chunk_candidates() = something(tryparse(Int, get(ENV, "PIONEER_ZT_CHUNK_CANDIDATES", "")),
                                  ZT_CHUNK_CANDIDATES_DEFAULT)

"""
    zt_cycle_scan_ranges(spectra) -> Vector{UnitRange{Int}}

Contiguous MS2 scan ranges, one per acquisition cycle, in scan order.
"""
function zt_cycle_scan_ranges(spectra::MassSpecData)
    cycles = getCycleIdxs(spectra)
    n = length(spectra)
    out = UnitRange{Int}[]
    i = 1
    while i <= n
        if getMsOrder(spectra, i) != 2
            i += 1; continue
        end
        c = cycles[i]; j = i
        while j + 1 <= n && getMsOrder(spectra, j + 1) == 2 && cycles[j + 1] == c
            j += 1
        end
        push!(out, i:j)
        i = j + 1
    end
    return out
end

"""
    zt_candidate_count(cycle_ranges, scan_to_prec_idx) -> Int

Exact number of (precursor, scan) candidates in the given cycles, from the per-scan ranges the
fragment index + expansion produced.
"""
function zt_candidate_count(cycle_ranges::Vector{UnitRange{Int}},
                            scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}})
    n = 0
    @inbounds for r in cycle_ranges, si in r
        rng = scan_to_prec_idx[si]
        ismissing(rng) || (n += length(rng))
    end
    return n
end

"""
    zt_cycle_chunks_by_candidates(ranges, scan_to_prec_idx, target) -> Vector{Vector{UnitRange{Int}}}

Group consecutive cycles into chunks of >= `target` candidates. Every chunk boundary is a cycle
boundary, so no meta-scan is split: candidate expansion and the collapse are both confined to
one cycle. Dense elution regions get many small chunks, empty regions one large chunk.
"""
function zt_cycle_chunks_by_candidates(ranges::Vector{UnitRange{Int}},
                                       scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                                       target::Int)
    chunks = Vector{Vector{UnitRange{Int}}}()
    cur = UnitRange{Int}[]; acc = 0
    for r in ranges
        push!(cur, r); acc += zt_candidate_count([r], scan_to_prec_idx)
        if acc >= target
            push!(chunks, cur); cur = UnitRange{Int}[]; acc = 0
        end
    end
    isempty(cur) || push!(chunks, cur)
    return chunks
end

"""
    zt_thread_tasks(cycle_ranges, k, n_threads) -> Vector{Tuple{Int, Vector{Int}}}

Deal one chunk's cycles to threads so that at any moment every thread is working the SAME m/z
segment of the ramp, on different cycles. Segment width is one meta-scan (2k+1 bins). Thread t
owns cycles t, t+T, t+2T, ... and walks them segment-major: segment 1 of each of its cycles, then
segment 2, and so on. All threads therefore touch the same precursors' library entries at the
same time (cache locality), and each thread processes a full meta-scan width consecutively (what
a warm-started solver and per-precursor template reuse will need). Per-thread vectors are exactly
sized up front; nothing grows.
"""
function zt_thread_tasks(cycle_ranges::Vector{UnitRange{Int}}, k::Int, n_threads::Int)
    C = length(cycle_ranges)
    W = 2k + 1
    T = min(n_threads, C)
    counts = zeros(Int, T)
    @inbounds for c in 1:C
        counts[mod1(c, T)] += length(cycle_ranges[c])
    end
    tasks = [(t, Vector{Int}(undef, counts[t])) for t in 1:T]
    pos = ones(Int, T)
    L = maximum(length, cycle_ranges)
    S = cld(L, W)
    @inbounds for s in 1:S, c in 1:C
        r = cycle_ranges[c]
        lo = first(r) + (s - 1) * W
        hi = min(lo + W - 1, last(r))
        lo > hi && continue
        t = mod1(c, T); v = last(tasks[t]); p = pos[t]
        for si in lo:hi
            v[p] = si; p += 1
        end
        pos[t] = p
    end
    return tasks
end

"""Number of complete cycles sampled by `detect_zt_geometry`."""
const ZT_GEOM_SAMPLE_CYCLES = 8

"""
Default effective transmission FWHM (Da). Measured on the reference ZT 5Da data by two
independent routes — isotope-ratio integration (6.30 Da) and a Razo fit with its width bound
relaxed (6.84-6.97 Da) — agreeing with a weight-profile fit from the prior effort (~6.7 Da).
Override per acquisition with `acquisition.transmission_fwhm_mz`.
"""
const ZT_TRANSMISSION_FWHM_DEFAULT = 6.5f0

"""
Overhang (Da) added to the recorded half-width for the ZT fragment-index candidacy box. Zero, so
candidacy is exactly the precursor's own Q1 bin (+/-w/2 ~ +/-0.51 Da) — the narrow half of the
two-width quad model. Isotope slack is supplied separately by `isotope_err_bounds`, which widens
the low side by ~0.5 Da to catch precursors whose M+1 sits in the bin.
"""
const ZT_FRAG_OVERHANG = 0.0f0

"""
    zt_frag_overhang() -> Float32

Candidacy overhang. Zero: candidacy is exactly the precursor's own Q1 bin.

Measured on A_REP1 — widening to 0.5 (a +/-1.01 Da box) costs +26% emissions
(30.9M -> 39.1M) for +1.9% centers and +71 IDs (25,141 -> 25,212). `filter_to_center_bin!`
discards everything outside +/-w/2 regardless, so those extra emissions are computed and thrown
away; the only channel by which they help is a marginally different LUT calibration. Not worth it.
"""
zt_frag_overhang() = ZT_FRAG_OVERHANG

"""
    zt_candidacy_tol() -> Float32

Wide-emit candidacy half-width in Da (`PIONEER_ZT_CANDIDACY_TOL`); 0 disables it.

Wide-emit widens the fragment-index box so a precursor gets an emission CHANCE in every bin of
its meta-scan, then `map_any_hit_to_center!` re-anchors each emission to the precursor's own bin.
It survives if it cleared the bitvec in ANY bin, rather than needing its center bin to clear.
"""
zt_candidacy_tol() =
    something(tryparse(Float32, get(ENV, "PIONEER_ZT_CANDIDACY_TOL", "")), 0.0f0)

"""
    zt_candidacy_overhang(g::ZTGeometry) -> Float32

Overhang for the fragment-index candidacy box: the wide-emit half-width when one is set, else
the precursor's own bin. Single source of truth, so BitVecCalibration's LUT keeps mirroring
candidacy whichever mode is active — a mismatch there cost 2,008 IDs when measured.
"""
function zt_candidacy_overhang(g::ZTGeometry)
    tol = zt_candidacy_tol()
    hw = g.nominal_width / 2
    return tol > hw ? tol - hw : zt_frag_overhang()
end

"""
Extra margin (Da) by which the deconvolution box exceeds the meta-scan expansion span, so the
box never clips a bin that `expand_to_metascans!` deliberately included.
"""
const ZT_DECONV_MARGIN = 0.5f0

"""
    zt_deconv_overhang(g::ZTGeometry) -> Float32

Overhang for the deconvolution box: enough to span `+/-k` bins plus a margin. With
`SquareQuadModel` the resulting window is `centerMz +/- (w/2 + overhang)`, so this makes the box
strictly wider than the expansion span and keeps `metascan_k` the single control.
"""
zt_deconv_overhang(g::ZTGeometry) =
    Float32(g.metascan_k) * g.bin_step + ZT_DECONV_MARGIN

"""
    detect_zt_geometry(spectra::MassSpecData, metascan_k::Integer) -> Union{Nothing, ZTGeometry}

Measure the Q1 bin lattice from the first `ZT_GEOM_SAMPLE_CYCLES` complete cycles.

Sampling a handful of cycles is sufficient: within a single file the lattice is constant to
within Float32 granularity (ramp start std 0.0, step std ~3e-8 over 818 cycles on the reference
ZT data). Returns `nothing` when the file carries no usable MS2 isolation metadata.
"""
function detect_zt_geometry(spectra::MassSpecData, metascan_k::Integer,
                            transmission_fwhm::Real = ZT_TRANSMISSION_FWHM_DEFAULT)
    widths       = Float32[]
    spacings     = Float32[]
    cycle_counts = Int32[]

    unset       = typemax(UInt32)
    prev_cycle  = unset
    prev_center = NaN32
    n_in_cycle  = 0

    @inbounds for si in 1:length(spectra)
        getMsOrder(spectra, si) == 2 || continue
        w = getIsolationWidthMz(spectra, si); ismissing(w) && continue
        c = getCenterMz(spectra, si);         ismissing(c) && continue
        cyc = getCycleIdx(spectra, si)

        if prev_cycle == unset
            prev_cycle = cyc
        elseif cyc != prev_cycle
            push!(cycle_counts, Int32(n_in_cycle))
            length(cycle_counts) >= ZT_GEOM_SAMPLE_CYCLES && break
            prev_cycle  = cyc
            n_in_cycle  = 0
            prev_center = NaN32
        end

        n_in_cycle += 1
        push!(widths, Float32(w))
        # abs() so either sweep direction works; prev_center is reset at every cycle boundary,
        # so this never spans a ramp wrap.
        isnan(prev_center) || push!(spacings, abs(Float32(c) - prev_center))
        prev_center = Float32(c)
    end

    (isempty(widths) || isempty(spacings)) && return nothing

    return ZTGeometry(
        median(spacings),
        median(widths),
        isempty(cycle_counts) ? Int32(n_in_cycle) : Int32(round(median(cycle_counts))),
        Int32(metascan_k),
        Float32(transmission_fwhm),
        0f0,
    )
end

"""
    zt_tiles_contiguously(g::ZTGeometry; rtol = 0.02) -> Bool

True when the Q1 windows tile edge-to-edge (`nominal_width ≈ bin_step`) — the signature of a
reconstructed sweep. Static-quad DIA overlaps or leaves gaps and misses by tens of percent.

The 2% tolerance absorbs the width dithering (0.2% on the reference ZT data).
"""
zt_tiles_contiguously(g::ZTGeometry; rtol::Float64 = 0.02) =
    g.bin_step > 0 && abs(g.nominal_width - g.bin_step) <= rtol * g.bin_step

"""
    zt_span_mz(g::ZTGeometry) -> Float32

Full m/z width spanned by a `±metascan_k` expansion.
"""
zt_span_mz(g::ZTGeometry) = Float32(2 * g.metascan_k + 1) * g.bin_step


"""
    zt_with_metascan_k(g::ZTGeometry, k::Integer) -> ZTGeometry

Same measured lattice, different expansion half-width. Used when `metascan_k` is DERIVED from
the fitted transmission profile (`k_implied = round(h / bin_step)`) after quad tuning.
"""
zt_with_metascan_k(g::ZTGeometry, k::Integer) =
    ZTGeometry(g.bin_step, g.nominal_width, g.bins_per_ramp, Int32(k), g.transmission_fwhm,
               g.template_h)

"""
    zt_with_template_h(g::ZTGeometry, h::Real) -> ZTGeometry

Same geometry with the fitted transmission half-base installed, so `zt_transmission_template`
returns the measured triangle rather than the configured Gaussian.
"""
zt_with_template_h(g::ZTGeometry, h::Real) =
    ZTGeometry(g.bin_step, g.nominal_width, g.bins_per_ramp, g.metascan_k, g.transmission_fwhm,
               Float32(h))

"""
    zt_transmission_template(g::ZTGeometry, k::Int) -> (Vector{Float32}, Float32)

Expected weight profile across the `2k+1` sampled bins, and its Euclidean norm.

The shape comes from the acquisition's transmission FWHM, NOT from `k`: the profile is a fixed
physical property and `k` only decides how much of it we sample. A `k`-dependent template (the
earlier `1 - |j|/(k+1)` triangle) makes the same physical bin claim different transmission at
different `k` — 0.25 at j=3 with k=3 versus 0.57 with k=6, against a true value near 0.52 — so
`zt_tri_cosine` was measured against a different yardstick at each `k`.
"""
function zt_transmission_template(g::ZTGeometry, k::Int)
    t = if g.template_h > 0f0
        # Measured triangle (QuadTuningSearch): linear to zero at +/-h.
        Float32[max(0f0, 1f0 - abs(Float32(j) * g.bin_step) / g.template_h) for j in -k:k]
    else
        σ = g.transmission_fwhm / 2.3548f0             # FWHM -> Gaussian sigma, in Da
        Float32[exp(-(Float32(j) * g.bin_step)^2 / (2f0 * σ * σ)) for j in -k:k]
    end
    return t, sqrt(sum(x -> x * x, t))
end

