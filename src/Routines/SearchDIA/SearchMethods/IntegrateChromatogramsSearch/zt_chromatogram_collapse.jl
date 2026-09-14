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

# ============================================================================
# Scanning-quad (ZT) chromatogram collapse.
#
# IntegrateChromatogramsSearch re-extracts chromatograms from the raw spectra rather than
# reusing MainSearch's collapsed table, so on a scanning quad each precursor's trace carries one
# point per Q1 bin — ~2k+1 points per cycle, spanning only a few tens of milliseconds. The
# integrator assumes rt is monotonic in scan_idx and that positions are row indices, so it would
# treat those as 2k+1 separate time samples when they are ONE time point measured 2k+1 ways.
#
# Collapsing to one point per (precursor, cycle) before smoothing restores the one-point-per-
# cycle cadence the integrator expects.
# ============================================================================

"""
    MetascanReduction

How a precursor's metascan bins within one cycle are reduced to a single chromatogram point.
The choice is a genuine tradeoff, so it dispatches rather than being hardcoded:

- [`SumMetascan`](@ref)      — assumption-free, the most direct measure of total transmitted signal
- [`TriangleMetascan`](@ref) — down-weights outer bins, where transmission is low and interference
                               is proportionally worse

A third option, weighting by the *measured* transmission profile, is the maximum-likelihood
estimator for abundance given a known profile shape. It is not implemented yet because it needs
the per-file transmission curve threaded through; add it as another subtype when it is.
"""
abstract type MetascanReduction end

"""
Sum the bin intensities. With `intensity(bin) ~ abundance * T(offset) + noise` this estimates
`abundance * sum(T)`, where `sum(T)` depends only on the precursor's sub-bin m/z offset and is
therefore constant for a given precursor across runs — so ratios are preserved. Noise
accumulates over all bins.
"""
struct SumMetascan <: MetascanReduction end

"""
Triangle-weighted mean of the bin intensities, matching the template the shape features use.
Trades a little signal for less noise: outer bins contribute least, which is where transmission
is lowest and interference proportionally largest.
"""
struct TriangleMetascan <: MetascanReduction end

"""
Triangle-weighted SUM over only the innermost `±half_width` bins, with the triangle still shaped
by the full `k` — so the weight at the window edge is `1 - half_width/(k+1)` (0.714 at ±2 for
k=6), NOT zero. Keeps the high-transmission core and discards the outer bins entirely, where
transmission is lowest (~8% at ±6 on a 6.3 Da FWHM) and interference is proportionally worst.

A sum rather than a normalized mean, so it stays on the same footing as `SumMetascan`: the weight
total depends only on which bins are present, which is constant per precursor across runs, so
ratios are preserved.
"""
struct TriangleWindowMetascan <: MetascanReduction
    half_width::Int
end

"""
    reduce_metascan(reduction, intensities, offsets, k) -> Float32

Reduce one cycle's metascan bins to a single value. `offsets[i]` is the bin's scan-index offset
from the center bin (0 at center), and `intensities[i]` its deconvolved intensity.
"""
@inline function reduce_metascan(::SumMetascan, intensities::AbstractVector{Float32},
                                 offsets::AbstractVector{Int}, k::Int)
    s = zero(Float32)
    @inbounds for v in intensities
        s += v
    end
    return s
end

@inline function reduce_metascan(::TriangleMetascan, intensities::AbstractVector{Float32},
                                 offsets::AbstractVector{Int}, k::Int)
    kf = Float32(k + 1)
    num = zero(Float32); den = zero(Float32)
    @inbounds for i in eachindex(intensities)
        t = max(0f0, 1f0 - abs(Float32(offsets[i])) / kf)
        num += t * intensities[i]
        den += t
    end
    return den > 0f0 ? num / den : 0f0
end

@inline function reduce_metascan(r::TriangleWindowMetascan, intensities::AbstractVector{Float32},
                                 offsets::AbstractVector{Int}, k::Int)
    kf = Float32(k + 1)
    hw = r.half_width
    s = zero(Float32)
    @inbounds for i in eachindex(intensities)
        o = offsets[i]
        if abs(o) <= hw
            s += max(0f0, 1f0 - abs(Float32(o)) / kf) * intensities[i]
        end
    end
    return s
end

"""
    collapse_chromatograms_to_metascans(chromatograms, spectra, precursors, k; reduction) -> DataFrame

Collapse each precursor's metascan bins within a cycle into ONE chromatogram point, keeping the
center bin's row (so `rt`, `scan_idx` and every other column stay consistent and `rt` stays
monotonic in `scan_idx`) and replacing its `:intensity` with the reduction over the cycle's bins.

A row is the cycle's center iff its scan's `centerMz` is the closest of that cycle's bins to the
precursor m/z. Mirrors `collapse_to_metascans` in MainSearch, which keeps center rows the same
way. Returns `chromatograms` unchanged when `k <= 0` or the table is empty.
"""
function collapse_chromatograms_to_metascans(
    chromatograms::DataFrame,
    spectra::MassSpecData,
    precursors,
    k::Int;
    reduction::MetascanReduction = SumMetascan(),
)
    n = nrow(chromatograms)
    (n == 0 || k <= 0) && return chromatograms

    # Function barrier: these accessors box on every index and are read once per row.
    prec_mz::Vector{Float32} = Vector{Float32}(getMz(precursors))
    cmzs::Vector{Float32} = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
    cycles::Vector{UInt32} = Vector{UInt32}(getCycleIdxs(spectra))

    pid0 = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    scn0 = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    int0 = chromatograms[!, :intensity]::AbstractVector{Float32}

    # Sort by (precursor, scan). Scan index increases with cycle, so a precursor's cycles come
    # out contiguous and in order — no need for a three-level key.
    sortkeys = Vector{UInt64}(undef, n)
    @inbounds for i in 1:n
        sortkeys[i] = (UInt64(pid0[i]) << 32) | UInt64(scn0[i])
    end
    perm = sortperm(sortkeys; alg = QuickSort)
    sortkeys = UInt64[]

    pid = pid0[perm]; scn = scn0[perm]; ints = Float32.(int0[perm])

    center_rows = Int[];      sizehint!(center_rows, n ÷ (2k + 1) + 1)
    reduced     = Float32[];  sizehint!(reduced,     n ÷ (2k + 1) + 1)
    grp_int = Float32[]; grp_off = Int[]

    i = 1
    @inbounds while i <= n
        # one (precursor, cycle) group
        j0 = i
        p = pid[i]; c0 = cycles[scn[i]]
        while i <= n && pid[i] == p && cycles[scn[i]] == c0
            i += 1
        end
        lo, hi = j0, i - 1

        # center = the bin whose centerMz is nearest this precursor's m/z
        pm = prec_mz[p]
        best = lo; bestd = abs(pm - cmzs[scn[lo]])
        for r in (lo + 1):hi
            d = abs(pm - cmzs[scn[r]])
            if d < bestd
                bestd = d; best = r
            end
        end

        empty!(grp_int); empty!(grp_off)
        cs = Int(scn[best])
        for r in lo:hi
            push!(grp_int, ints[r])
            push!(grp_off, Int(scn[r]) - cs)
        end

        push!(center_rows, best)
        push!(reduced, reduce_metascan(reduction, grp_int, grp_off, k))
    end

    meta = chromatograms[perm[center_rows], :]
    meta[!, :intensity] = reduced
    return meta
end
