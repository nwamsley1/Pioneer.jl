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

struct GeneralGaussQuadParams{T <: AbstractFloat}
    a::T
    b::T
    overhang::T
end

function (rqm::GeneralGaussQuadParams)(x::T) where {T<:AbstractFloat}
        x = abs(x)
        return 1/(1 + abs(x/(rqm.a + rqm.overhang))^(2*rqm.b))
        #return one(T)
end


struct GeneralGaussModel{T<:AbstractFloat} <: QuadTransmissionModel
    b::T
    overhang::T
end

struct GeneralGaussFunction{T<:AbstractFloat} <: QuadTransmissionFunction 
    min_mz::T
    max_mz::T
    center_mz::T
    params::GeneralGaussQuadParams{T}
end


function getQuadTransmissionBounds(ggm::GeneralGaussModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    return T(centerMz - isolationWidthMz/2), T(centerMz + isolationWidthMz/2)
end

function getPrecMinBound(ggf::GeneralGaussFunction{T}) where {T<:AbstractFloat}
    return ggf.min_mz
end

function getPrecMaxBound(ggf::GeneralGaussFunction{T}) where {T<:AbstractFloat}
    return ggf.max_mz
end


function (ggm::GeneralGaussFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    ggm.params(ionMz-ggm.center_mz)
end

function getQuadTransmissionFunction(ggm::GeneralGaussModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    half_width_mz = isolationWidthMz/T(2)
    GeneralGaussFunction(
        centerMz - half_width_mz,
        centerMz + half_width_mz,
        centerMz,
        GeneralGaussQuadParams(
            half_width_mz,
            ggm.b,
            ggm.overhang
        )
    )
end

