using Dictionaries
abstract type Match end
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
struct FragmentMatch{T<:AbstractFloat} <: Match
    predicted_intensity::T
    intensity::T
    theoretical_mz::T
    match_mz::T
    peak_ind::Int64
    frag_index::UInt8
    frag_charge::UInt8
    frag_isotope::UInt8
    ion_type::Char
    is_isotope::Bool
    prec_id::UInt32
    count::UInt8
    scan_idx::UInt32
    ms_file_idx::UInt32
    predicted_rank::UInt8
end

FragmentMatch{Float64}() = FragmentMatch(Float64(0), Float64(0), Float64(0), Float64(0), 0, UInt8(0), UInt8(0), UInt8(0),'y', false, UInt32(0), UInt8(0), UInt32(0), UInt32(0), zero(UInt8))
FragmentMatch{Float32}() = FragmentMatch(Float32(0), Float32(0), Float32(0), Float32(0), 0, UInt8(0), UInt8(0), UInt8(0),'y', false, UInt32(0), UInt8(0), UInt32(0), UInt32(0), zero(UInt8))


getFragMZ(f::FragmentMatch) = f.theoretical_mz
getMatchMZ(f::FragmentMatch) = f.match_mz
getPredictedIntenisty(f::FragmentMatch) = f.predicted_intensity
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

struct PrecursorMatch{T<:AbstractFloat} <: Match
    predicted_intensity::T
    intensity::T
    peak_ind::Int64
    prec_id::UInt32
end

getPrecID(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.prec_id
getIntensity(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.intensity
getPeakInd(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.peak_ind
getPredictedIntenisty(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = pm.predicted_intensity
PrecursorMatch{Float32}() = PrecursorMatch(zero(Float32), zero(Float32), zero(Int64), zero(UInt32))
getRank(pm::PrecursorMatch{T}) where {T<:AbstractFloat} = one(UInt8)
getFragInd(::PrecursorMatch{Float32}) = Inf
getIonType(::PrecursorMatch{Float32}) = 'y'

abstract type FragmentIndexType end

getFragMZ(f::FragmentIndexType) = f.frag_mz
getPrecID(f::FragmentIndexType) = f.prec_id
getPrecCharge(f::FragmentIndexType) = f.prec_charge

import Base.<
import Base.>

<(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) < x
<(x::T, y::FragmentIndexType) where {T<:Real} = <(y, x)
>(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) > x
>(x::T, y::FragmentIndexType) where {T<:Real} = >(y, x)

struct FragmentIon{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    prec_id::UInt32
    prec_mz::T
    prec_charge::UInt8
end

getPrecMZ(f::FragmentIon) = f.prec_mz

struct LibraryFragment{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    frag_charge::UInt8
    is_y_ion::Bool
    is_isotope::Bool
    ion_position::UInt8
    ion_index::UInt8
    intensity::Float32
    prec_charge::UInt8
    prec_id::UInt32
    rank::UInt8
    sulfur_count::UInt8
end

getMZ(f::LibraryFragment) = f.frag_mz
getIntensity(f::LibraryFragment) = f.intensity
isyIon(f::LibraryFragment) = f.is_y_ion
getIonIndex(f::LibraryFragment) = f.ion_index
getIonPosition(f::LibraryFragment) = f.ion_position
getFragCharge(f::LibraryFragment) = f.frag_charge
getRank(f::LibraryFragment) = f.rank
sulfurCount(f::LibraryFragment) = f.sulfur
#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), 
zero(UInt8), 
false, 
false,
zero(UInt8), 
zero(UInt8), 
zero(Float32), 
zero(UInt8),
zero(UInt32), 
zero(UInt8),
zero(UInt8)
 )

 struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    mz::T
    total_intensity::T
    base_peak_intensity::T
    isDecoy::Bool
    charge::UInt8
    pep_id::UInt32
    prot_ids::Vector{UInt32}
    accession_numbers::String
    sequence::String
    missed_cleavages::UInt8
    variable_mods::UInt8
    length::UInt8
    sulfur_count::UInt8
end

isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT
getCharge(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.charge
getMz(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.mz
getTotalIntensity(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.total_intensity
getPepID(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.pep_id
getBasePeakInt(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.base_peak_intensity
sulfurCount(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.sulfur_count