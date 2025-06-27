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

function reset!(fms::Vector{M}, last_non_empty::Int64) where {M<:MatchIon{Float32}}
    for i in range(1, last_non_empty)
        fms[i] = M()
    end
end

"""
    FragmentMatch

Type that represents a match between a fragment ion and a mass spectrum peak

### Fields

- transition::Transition -- Represents a fragment ion
- intensity::Float32 -- Intensity of the matching peak
- match_mz::Float32 -- M/Z of the matching empirical peak. 
    NOT the transition mass and may differ from getMZ(transition) by some error
- count::UInt8 -- Number of matches (may want to count the number of matches if the spectrum is researched at 
    different mass offsets, such as in a cross correlation score)
- peak_int::Int64 -- Index of the matching peak in the mass spectrum. 

### Examples

- `FragmentMatch(transition::Transition, intensity::Float32, mass::Float32, count::UInt8, peak_ind::Int64) --
    default internal constructor
- `FragmentMatch()` -- constructor for null/empty precursor

### GetterMethods

- getMZ(f::FragmentMatch) = getMZ(f.transition)
- getLow(f::FragmentMatch) = getLow(f.transition)
- getHigh(f::FragmentMatch) = getHigh(f.transition)
- getPrecID(f::FragmentMatch) = getPrecID(f.transition)
- getCharge(f::FragmentMatch) = getCharge(f.transition)
- getIsotope(f::FragmentMatch) = getIsotope(f.transition)
- getIonType(f::FragmentMatch) = getIonType(f.transition)
- getInd(f::FragmentMatch) = getInd(f.transition)
"""
struct FragmentMatch{T<:AbstractFloat} <: MatchIon{T}
    predicted_intensity::T
    intensity::T
    theoretical_mz::T
    match_mz::T
    peak_ind::Int64
    frag_index::UInt8
    frag_charge::UInt8
    frag_isotope::UInt8
    ion_type::UInt8
    is_isotope::Bool
    prec_id::UInt32
    count::UInt8
    scan_idx::UInt32
    ms_file_idx::UInt32
    predicted_rank::UInt8
end

FragmentMatch{Float64}() = FragmentMatch(Float64(0), Float64(0), Float64(0), Float64(0), 0, UInt8(0), UInt8(0), UInt8(0),'y', false, UInt32(0), UInt8(0), UInt32(0), UInt32(0), zero(UInt8))
FragmentMatch{Float32}() = FragmentMatch(Float32(0), Float32(0), Float32(0), Float32(0), 0, UInt8(0), UInt8(0), UInt8(0), zero(UInt8), false, UInt32(0), UInt8(0), UInt32(0), UInt32(0), zero(UInt8))


getMZ(f::FragmentMatch) = f.theoretical_mz
getMatchMZ(f::FragmentMatch) = f.match_mz
getPredictedIntensity(f::FragmentMatch) = f.predicted_intensity
getIntensity(f::FragmentMatch) = f.intensity

getPeakInd(f::FragmentMatch) = f.peak_ind
getFragInd(f::FragmentMatch) = f.frag_index
getCharge(f::FragmentMatch) = f.frag_charge
getIsotope(f::FragmentMatch) = f.frag_isotope
getIonType(f::FragmentMatch) = f.ion_type

getPrecID(f::FragmentMatch) = f.prec_id
getCount(f::FragmentMatch) = f.count
getScanID(f::FragmentMatch) = f.scan_idx
getMSFileID(f::FragmentMatch) = f.ms_file_idx
getRank(f::FragmentMatch) = f.predicted_rank
isIsotope(f::FragmentMatch) = f.is_isotope
getIsoIdx(f::FragmentMatch) = UInt8(f.is_isotope)

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
getMatchMz(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.observed_mz
PrecursorMatch{Float32}() = PrecursorMatch(zero(Float32), 
zero(Float32), 
zero(Float32), 
zero(Float32),
zero(UInt8), zero(Int64), zero(UInt32))
