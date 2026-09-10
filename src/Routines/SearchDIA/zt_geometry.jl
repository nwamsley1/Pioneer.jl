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

`bin_step` is deliberately taken from `centerMz` differences rather than `isolationWidthMz`:
the recorded width is dithered 50/50 between two quantized values (e.g. 1.0200 / 1.0241) whose
*average* is the true step, so a median over that column is degenerate and lands ~2e-3 Da low.
"""
struct ZTGeometry
    bin_step::Float32
    nominal_width::Float32
    bins_per_ramp::Int32
    metascan_k::Int32
end

"""Number of complete cycles sampled by `detect_zt_geometry`."""
const ZT_GEOM_SAMPLE_CYCLES = 8

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
function detect_zt_geometry(spectra::MassSpecData, metascan_k::Integer)
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
