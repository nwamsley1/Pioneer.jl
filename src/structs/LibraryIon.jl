abstract type LibraryIon{T<:AbstractFloat} <: Ion{T} end

getPrecCharge(f::LibraryIon)::UInt8 = f.prec_charge

struct LibraryPrecursorIon{T<:AbstractFloat} <: LibraryIon{T}
    irt::T
    mz::T

    is_decoy::Bool

    proteome_identifiers::String
    accession_numbers::String
    sequence::String
    structural_mods::Union{Missing, String}
    isotopic_mods::Union{Missing, String}

    prec_charge::UInt8
    missed_cleavages::UInt8
    length::UInt8
    sulfur_count::UInt8
end

getIRT(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.irt
isDecoy(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.is_decoy
getAccessionNumbers(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.accession_numbers
getSequence(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.sequence
getMissedCleavages(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.missed_cleavages
getVariableMods(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.variable_mods
getLength(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.length
sulfurCount(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.sulfur_count
getCharge(p::LibraryPrecursorIon{T}) where {T<:AbstractFloat} = p.prec_charge

abstract type LibraryFragmentIon{T<:AbstractFloat} <: LibraryIon{T} end

getPrecID(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_id
getPrecMZ(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_mz
getIntensity(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.intensity

function reset!(lf::Vector{L}, last_non_empty::Int64) where {L<:LibraryFragmentIon{Float32}}
    for i in range(1, last_non_empty)
        lf[i] = L()
    end
end

struct SimpleFrag{T<:AbstractFloat} <: LibraryFragmentIon{T}
    mz::T
    prec_id::UInt32
    prec_mz::T
    prec_irt::T
    prec_charge::UInt8
    score::UInt8
end
ArrowTypes.arrowname(::Type{SimpleFrag}) = :SimpleFrag
ArrowTypes.JuliaType(::Val{:SimpleFrag}) = SimpleFrag


getScore(pbi::SimpleFrag)::UInt8 = pbi.score
getIRT(pbi::SimpleFrag{T}) where {T<:AbstractFloat} = pbi.prec_irt


struct PioneerFrag
    mz::Float32
    intensity::Float16
    ion_type::UInt16 
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_axcz::Bool
    has_neutral_diff::Bool
    frag_index::UInt8 #posiiton of fragment
    charge::UInt8
    isotope::UInt8
    internal::Bool
    immonium::Bool
    internal_ind::Tuple{UInt8, UInt8} #If an internal ion, the start and stop. 0,0 if not internal
    sulfur_count::UInt8
end
ArrowTypes.arrowname(::Type{PioneerFrag}) = :PioneerFrag
ArrowTypes.JuliaType(::Val{:PioneerFrag}) = PioneerFrag

getIntensity(pf::PioneerFrag) = pf.intensity
getMZ(pf::PioneerFrag) = pf.mz
getType(pf::PioneerFrag) = pf.ion_type
getIndex(pf::PioneerFrag) = pf.frag_index
getCharge(pf::PioneerFrag) = pf.charge
getSulfurCount(pf::PioneerFrag) = pf.sulfur_count
isY(pf::PioneerFrag) = pf.is_y


struct PioneerSplineFrag{N}
    mz::Float32
    spl_coef::NTuple{N, Float32}
    intensity::Float16
    ion_type::UInt16 
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_axcz::Bool
    has_neutral_diff::Bool
    frag_index::UInt8 #posiiton of fragment
    charge::UInt8
    isotope::UInt8
    internal::Bool
    immonium::Bool
    internal_ind::Tuple{UInt8, UInt8} #If an internal ion, the start and stop. 0,0 if not internal
    sulfur_count::UInt8
end
ArrowTypes.arrowname(::Type{PioneerSplineFrag}) = :PioneerSplineFrag
ArrowTypes.JuliaType(::Val{:PioneerSplineFrag}) = PioneerSplineFrag

getIntensity(pf::PioneerSplineFrag) = pf.intensity
getMZ(pf::PioneerSplineFrag) = pf.mz
getType(pf::PioneerSplineFrag) = pf.ion_type
getIndex(pf::PioneerSplineFrag) = pf.frag_index
getCharge(pf::PioneerSplineFrag) = pf.charge
getSulfurCount(pf::PioneerSplineFrag) = pf.sulfur_count
isY(pf::PioneerSplineFrag) = pf.is_y

#Need this information for each distinct fragment type
#Need this information for each distinct fragment type
struct PioneerFragAnnotation
    base_type::Char
    frag_index::UInt8
    charge::UInt8
    isotope::UInt8
    internal::Bool
    immonium::Bool
    neutral_diff::Bool
    sulfur_diff::Int8
end
ArrowTypes.arrowname(::Type{PioneerFragAnnotation}) = :PioneerFragAnnotation
ArrowTypes.JuliaType(::Val{:PioneerFragAnnotation}) = PioneerFragAnnotation
getBaseType(pfa::PioneerFragAnnotation) = pfa.base_type

abstract type AltimeterFragment{T} <: LibraryFragmentIon{T} end

getPID(f::AltimeterFragment) = f.prec_id
getMz(f::AltimeterFragment) = f.mz
getIonType(f::AltimeterFragment) = f.ion_type
isY(f::AltimeterFragment) = f.is_y
isB(f::AltimeterFragment) = f.is_b
isP(f::AltimeterFragment) = f.is_p
isIso(f::AltimeterFragment) = f.is_isotope
getFragCharge(f::AltimeterFragment) = f.frag_charge
getIonPosition(f::AltimeterFragment) = f.ion_position
getPrecCharge(f::AltimeterFragment) = f.prec_charge
getRank(f::AltimeterFragment) = f.rank
sulfurCount(f::AltimeterFragment) = f.sulfur_count

struct DetailedFrag{T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32

    mz::T
    intensity::Float16
    
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_isotope::Bool

    frag_charge::UInt8
    ion_position::UInt8
    prec_charge::UInt8
    rank::UInt8
    sulfur_count::UInt8
end
# Example usage
function save_detailed_frags(filename::String, data::Vector{<:DetailedFrag})
    jldsave(filename; data)
end

function load_detailed_frags(filename::String)
    return load(filename, "data")
end
ArrowTypes.arrowname(::Type{DetailedFrag{Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:DetailedFrag}) = DetailedFrag

#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
DetailedFrag{T}() where {T<:AbstractFloat} = DetailedFrag(

                zero(UInt32), #prec_id

                zero(T), #mz
                zero(Float16), #intensity

                zero(UInt16),
                false,
                false,
                false,
                false,
                
                zero(UInt8), #frag_charge
                zero(UInt8), #ion_position
                zero(UInt8), #prec_charge
                zero(UInt8), #rank
                zero(UInt8)  #sulfur_count
 )
 
 struct SplineDetailedFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32

    mz::T
    intensity::NTuple{N, T}
    
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_isotope::Bool

    frag_charge::UInt8
    ion_position::UInt8
    prec_charge::UInt8
    rank::UInt8
    sulfur_count::UInt8
end
# Example usage
function save_detailed_frags(filename::String, data::Vector{<:SplineDetailedFrag})
    jldsave(filename; data)
end
function load_detailed_frags(filename::String)
    return load(filename, "data")
end
ArrowTypes.arrowname(::Type{SplineDetailedFrag{4, Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:SplineDetailedFrag}) = DetailedFrag
#This needs to be eddited for the splines
function getIntensity(
    pf::SplineDetailedFrag{N, T},
    knots::NTuple{M, T},
    degree::Int64,
    nce::T
    ) where {N,M,T<:AbstractFloat}
    return splevl(nce, knots, pf.intensity, degree)::T
end

#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
SplineDetailedFrag{N,T}() where {N,T<:AbstractFloat} = SplineDetailedFrag(

                zero(UInt32), #prec_id

                zero(T), #mz
                NTuple{N, T}(undef), #intensity

                zero(UInt16),
                false,
                false,
                false,
                false,
                
                zero(UInt8), #frag_charge
                zero(UInt8), #ion_position
                zero(UInt8), #prec_charge
                zero(UInt8), #rank
                zero(UInt8)  #sulfur_count
 )

abstract type LibraryFragmentLookup end

struct StandardFragmentLookup{T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{DetailedFrag{T}}
    prec_frag_ranges::Vector{UInt64}
end
getFrag(lfp::StandardFragmentLookup{<:AbstractFloat}, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::StandardFragmentLookup{<:AbstractFloat}) = lfp.frags
getPrecFragRange(lfp::StandardFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} = range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))

struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    nce::Base.RefValue{T}
end
getFrag(lfp::SplineFragmentLookup, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::SplineFragmentLookup) = lfp.frags
getPrecFragRange(lfp::SplineFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} = range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))
getNCE(lfp::SplineFragmentLookup) = lfp.nce[]
getKnots(lfp::SplineFragmentLookup) = lfp.knots
"""
    PrecursorBinItem{T<:AbstractFloat}

Item in a precursor bin. Minimal information required to know if the fragment has a correct precursor mass, and if so, 
what the precursor ID is.  

### Fields

- prec_id::UInt32 -- Unique identifier for the precursor
- prec_mz::T -- m/z of the precursor

### Examples

### GetterMethods

- getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
- getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
"""
struct PrecursorBinFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32
    prec_mz::T #Only need to tell if the peptide is in the quad isolation window
    score::UInt8 
    charge::UInt8
end

getScore(pbi::PrecursorBinFragment{T}) where {T<:AbstractFloat} = pbi.score

