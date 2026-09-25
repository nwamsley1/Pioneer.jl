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

- [`SumMetascan`](@ref)            — assumption-free total transmitted signal; the fallback
- [`MatchedFilterMetascan`](@ref)  — least-squares abundance under the per-file fitted triangle;
                                     the default when a triangle was fitted
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
Matched filter with the MEASURED transmission profile: bin weights `T_j = max(0, 1 - |j·S| / h)`
from the per-file fitted triangle (`h` in Da, `S` the bin step), estimate
`a = Σ T_j·w_j / Σ T_j²` — the least-squares abundance given `w_j ≈ a·T_j + noise`. Unlike the
sums it is invariant to WHICH bins are present (a missing bin drops out of both sums), and it
weights bins by their expected signal, so low-transmission outer bins contribute little. Scaled
by `Σ T_j` over the full ±k so the result is on the same footing as `SumMetascan` (total
transmitted signal), which keeps the downstream intensity scale unchanged.
"""
struct MatchedFilterMetascan <: MetascanReduction
    h_over_step::Float32          # fitted h / bin_step
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

@inline function reduce_metascan(r::MatchedFilterMetascan, intensities::AbstractVector{Float32},
                                 offsets::AbstractVector{Int}, k::Int)
    num = zero(Float32); den = zero(Float32); tsum = zero(Float32)
    @inbounds for i in eachindex(intensities)
        t = max(0f0, 1f0 - abs(Float32(offsets[i])) / r.h_over_step)
        num += t * intensities[i]
        den += t * t
    end
    @inbounds for j in -k:k
        tsum += max(0f0, 1f0 - abs(Float32(j)) / r.h_over_step)
    end
    return den > 0f0 ? (num / den) * tsum : 0f0
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

# --- IntegrateChromatogramsSearch hooks. The `::Nothing` methods (not a ZT file) are no-ops.

"""
    zt_collapse_chromatograms(geom, chromatograms, spectra, search_context) -> DataFrame

Collapse each precursor's meta-scan bins within a cycle into one chromatogram point before
smoothing and integration. Extraction selects a precursor in every scan its transmission window
covers, so on ZT the raw trace carries ~2k+1 points per cycle spanning a few tens of ms, which
the integrator would read as that many separate time samples.

The reduction is a matched filter on the per-file fitted triangle; vs a plain sum (3A+3B 5 Da,
2026-09-17) per-precursor log2(A/B) MAD 0.18/0.21/0.31 -> 0.17/0.18/0.27 (H/Y/E), CVs -0.5 pt,
medians unchanged. A plain sum when no triangle was fitted.
"""
zt_collapse_chromatograms(::Nothing, chromatograms::DataFrame, spectra, search_context) = chromatograms

function zt_collapse_chromatograms(g::ZTGeometry, chromatograms::DataFrame, spectra::MassSpecData,
                                   search_context::SearchContext)
    (g.metascan_k > 0 && nrow(chromatograms) > 0) || return chromatograms
    n_pre = nrow(chromatograms)
    red = g.template_h > 0f0 ? MatchedFilterMetascan(g.template_h / g.bin_step) : SumMetascan()
    out = collapse_chromatograms_to_metascans(chromatograms, spectra,
                                              getPrecursors(getSpecLib(search_context)),
                                              Int(g.metascan_k); reduction = red)
    @user_info "ZT chromatogram collapse (k=$(g.metascan_k), $(typeof(red))): " *
               "$n_pre -> $(nrow(out)) points"
    return out
end

"""
    zt_reset_transmission!(geom, chromatograms)

A collapsed meta-scan point already integrates the whole transmission window, so a per-point
transmission correction in WH smoothing would double-count it: set it to 1.
`get_isotopes_captured!` still runs first, because it also produces `isotopes_captured`, which
SeperateTraces mode groups on.
"""
zt_reset_transmission!(::Nothing, chromatograms::DataFrame) = nothing
function zt_reset_transmission!(g::ZTGeometry, chromatograms::DataFrame)
    g.metascan_k > 0 && (chromatograms[!, :precursor_fraction_transmitted] .= one(Float32))
    return nothing
end
