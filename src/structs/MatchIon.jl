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

abstract type MatchIon{T<:AbstractFloat} <: Ion{T} end

"""
    MzSortEntry

8-byte sidecar for sorting fragment ions by m/z without moving full 24B DetailedFrag structs.
The hot comparison loop in matchPeaks! reads only the 4B mz field; the idx field resolves
to the full template only on match/miss writes (infrequent relative to comparisons).
"""
struct MzSortEntry
    mz::Float32     # 4B — sort key, used directly by matchPeaks! comparisons
    idx::UInt32     # 4B — index into ion_templates array
end               # 8B, no padding

"""
    ion_match_lt(a, b) -> Bool

Compare two MatchIon instances by (peak_ind, prec_id) without tuple allocation.
Used as `lt` argument to `sort!` in hot loops (~40k calls per file).
"""
@inline ion_match_lt(a, b) = (a.peak_ind < b.peak_ind) || (a.peak_ind == b.peak_ind && a.prec_id < b.prec_id)

function reset!(fms::Vector{M}, last_non_empty::Int64) where {M<:MatchIon{Float32}}
    for i in range(1, last_non_empty)
        fms[i] = M()
    end
end

"""
    FragmentMatch{T} <: MatchIon{T}

Match between a fragment ion and a mass spectrum peak. 32 bytes for Float32.
"""
struct FragmentMatch{T<:AbstractFloat} <: MatchIon{T}
    predicted_intensity::T    # 4B  (offset 0)
    intensity::T              # 4B  (offset 4)
    theoretical_mz::T         # 4B  (offset 8)
    match_mz::T               # 4B  (offset 12)
    peak_ind::UInt32          # 4B  (offset 16) — was Int64
    prec_id::UInt32           # 4B  (offset 20)
    frag_index::UInt8         # 1B  (offset 24)
    frag_charge::UInt8        # 1B  (offset 25)
    ion_type::UInt8           # 1B  (offset 26)
    is_isotope::Bool          # 1B  (offset 27)
    predicted_rank::UInt8     # 1B  (offset 28)
end                           # 29B → padded to 32B

FragmentMatch{Float64}() = FragmentMatch(Float64(0), Float64(0), Float64(0), Float64(0), zero(UInt32), zero(UInt32), UInt8(0), UInt8(0), zero(UInt8), false, zero(UInt8))
FragmentMatch{Float32}() = FragmentMatch(Float32(0), Float32(0), Float32(0), Float32(0), zero(UInt32), zero(UInt32), UInt8(0), UInt8(0), zero(UInt8), false, zero(UInt8))


getMZ(f::FragmentMatch) = f.theoretical_mz
getMatchMZ(f::FragmentMatch) = f.match_mz
getMatchMz(f::FragmentMatch) = f.match_mz
getPredictedIntensity(f::FragmentMatch) = f.predicted_intensity
getIntensity(f::FragmentMatch) = f.intensity

getPeakInd(f::FragmentMatch) = f.peak_ind
getFragInd(f::FragmentMatch) = f.frag_index
getCharge(f::FragmentMatch) = f.frag_charge
getIonType(f::FragmentMatch) = f.ion_type

getPrecID(f::FragmentMatch) = f.prec_id
getRank(f::FragmentMatch) = f.predicted_rank
isIsotope(f::FragmentMatch) = f.is_isotope
getIsoIdx(f::FragmentMatch) = UInt8(f.is_isotope)

"""
    UnmatchedIon

Lightweight struct for unmatched (missed) fragment/precursor ions.
Only stores the three fields actually read by buildDesignMatrix!.
"""
struct UnmatchedIon
    prec_id::UInt32              # 4B
    predicted_intensity::Float32  # 4B
    is_isotope::Bool             # 1B
end                              # 9B → padded to 12B

UnmatchedIon() = UnmatchedIon(zero(UInt32), zero(Float32), false)
getPrecID(f::UnmatchedIon) = f.prec_id
getPredictedIntensity(f::UnmatchedIon) = f.predicted_intensity
getIsoIdx(f::UnmatchedIon) = UInt8(f.is_isotope)

struct PrecursorMatch{T<:AbstractFloat} <: MatchIon{T}
    predicted_intensity::T
    intensity::T
    theoretical_mz::T
    observed_mz::T
    iso_idx::UInt8
    peak_ind::Int64
    prec_id::UInt32
end

getPrecID(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.prec_id
getIntensity(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.intensity
getPeakInd(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.peak_ind
getPredictedIntensity(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.predicted_intensity
getIsoIdx(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.iso_idx
getMZ(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.theoretical_mz
getMatchMZ(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.observed_mz
getMatchMz(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.observed_mz
PrecursorMatch{Float32}() = PrecursorMatch(zero(Float32), 
zero(Float32), 
zero(Float32), 
zero(Float32),
zero(UInt8), zero(Int64), zero(UInt32))
