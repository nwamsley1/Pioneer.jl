"""
Base fragment type containing core information about a peptide fragment
"""
abstract type AbstractKoinaFragment end

getIntensity(pf::AbstractKoinaFragment) = pf.intensity
getMZ(pf::AbstractKoinaFragment) = pf.mz
getType(pf::AbstractKoinaFragment) = pf.ion_type
getIndex(pf::AbstractKoinaFragment) = pf.frag_index
getCharge(pf::AbstractKoinaFragment) = pf.charge
getSulfurCount(pf::AbstractKoinaFragment) = pf.sulfur_count
isY(pf::AbstractKoinaFragment) = pf.is_y


"""
Standard fragment with intensity information
"""
struct PioneerFrag <: AbstractKoinaFragment
    mz::Float32
    intensity::Float16
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_other::Bool
    has_neutral_diff::Bool
    fragment_index::UInt8
    charge::UInt8
    isotope::UInt8
    is_internal::Bool
    is_immonium::Bool
    sequence_bounds::Tuple{UInt8,UInt8}
    sulfur_count::UInt8
end

ArrowTypes.arrowname(::Type{PioneerFrag}) = :PioneerFrag
ArrowTypes.JuliaType(::Val{:PioneerFrag}) = PioneerFrag

"""
Fragment with spline coefficients for intensity modeling
"""
struct PioneerSplineFrag{N} <: AbstractKoinaFragment
    mz::Float32
    coefficients::NTuple{N,Float32}
    intensity::Float16
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_other::Bool
    has_neutral_diff::Bool
    fragment_index::UInt8
    charge::UInt8
    isotope::UInt8
    is_internal::Bool
    is_immonium::Bool
    sequence_bounds::Tuple{UInt8,UInt8}
    sulfur_count::UInt8
end

ArrowTypes.arrowname(::Type{PioneerSplineFrag}) = :PioneerSplineFrag
ArrowTypes.JuliaType(::Val{:PioneerSplineFrag}) = PioneerSplineFrag

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

# src/fragments/fragment_parser.jl

