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
