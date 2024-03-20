abstract type LibraryIon{T<:AbstractFloat} <: Ion{T} end

getPrecCharge(f::LibraryIon)::UInt8 = f.prec_charge

struct LibraryPrecursorIon{T<:AbstractFloat} <: LibraryIon{T}
    irt::T
    mz::T

    is_decoy::Bool

    accession_numbers::String
    sequence::String

    prec_charge::UInt8
    missed_cleavages::UInt8
    variable_mods::UInt8
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
getScore(pbi::SimpleFrag)::UInt8 = pbi.score

struct DetailedFrag{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32

    mz::T
    intensity::Float16
    
    is_y_ion::Bool
    is_isotope::Bool

    frag_charge::UInt8
    ion_position::UInt8
    prec_charge::UInt8
    rank::UInt8
    sulfur_count::UInt8
end

ArrowTypes.arrowname(::Type{DetailedFrag{Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:DetailedFrag}) = DetailedFrag


getIntensity(f::DetailedFrag) = f.intensity
isyIon(f::DetailedFrag) = f.is_y_ion
getIonIndex(f::DetailedFrag) = f.ion_index
getIonPosition(f::DetailedFrag) = f.ion_position
getFragCharge(f::DetailedFrag) = f.frag_charge
getRank(f::DetailedFrag) = f.rank
sulfurCount(f::DetailedFrag) = f.sulfur
#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
DetailedFrag{T}() where {T<:AbstractFloat} = DetailedFrag(

                zero(UInt32), #prec_id

                zero(T), #mz
                zero(Float16), #intensity

                false, #is_y_ion
                false, #is_isotope
                
                zero(UInt8), #frag_charge
                zero(UInt8), #ion_position
                zero(UInt8), #prec_charge
                zero(UInt8), #rank
                zero(UInt8)  #sulfur_count
 )
#=
struct PrecFragRange
    lb::UInt64
    ub::UInt64
end
ArrowTypes.arrowname(::Type{PrecFragRange}) = :PrecFragRange
ArrowTypes.JuliaType(::Val{:PrecFragRange}) = PrecFragRange


struct LibraryFragmentLookup{T<:AbstractFloat}
    frags::Arrow.Struct{DetailedFrag, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{Float16, Vector{Float16}}, Arrow.BoolVector{Bool}, Arrow.BoolVector{Bool}, Vararg{Arrow.Primitive{UInt8, Vector{UInt8}}, 5}}, (:prec_id, :mz, :intensity, :is_y_ion, :is_isotope, :frag_charge, :ion_position, :prec_charge, :rank, :sulfur_count)}
    prec_frag_ranges:: Arrow.Struct{PrecFragRange, Tuple{Arrow.Primitive{UInt64, Vector{UInt64}}, Arrow.Primitive{UInt64, Vector{UInt64}}}, (:lb, :ub)}
end

getFrag(lfp::LibraryFragmentLookup{<:AbstractFloat}, prec_idx::Integer) = lfp.frags[prec_idx]

function getPrecFragRange(lfp::LibraryFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} 
    prec_frag_range = lfp.prec_frag_ranges[prec_idx]
    return prec_frag_range.lb:prec_frag_range.ub
end
=#


struct LibraryFragmentLookup{T<:AbstractFloat}
    frags::Vector{DetailedFrag{T}}
    prec_frag_ranges::Vector{UnitRange{UInt32}}
end

getFrag(lfp::LibraryFragmentLookup{<:AbstractFloat}, prec_idx::Integer) = lfp.frags[prec_idx]

getPrecFragRange(lfp::LibraryFragmentLookup, prec_idx::Integer)::UnitRange{UInt32} = lfp.prec_frag_ranges[prec_idx]

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

