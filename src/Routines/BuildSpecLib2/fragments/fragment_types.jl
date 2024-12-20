"""
Base fragment type containing core information about a peptide fragment
"""
abstract type AbstractFragment end

"""
Standard fragment with intensity information
"""
struct PioneerFrag <: AbstractFragment
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

"""
Fragment with spline coefficients for intensity modeling
"""
struct PioneerSplineFrag{N} <: AbstractFragment
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

# src/fragments/fragment_parser.jl

