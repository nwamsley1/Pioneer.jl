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

struct SquareQuadModel{T<:AbstractFloat} <: QuadTransmissionModel
    overhang::T
end

struct SquareQuadFunction{T<:AbstractFloat} <: QuadTransmissionFunction 
    min_mz::T
    max_mz::T
    center_mz::T
end


function getQuadTransmissionFunction(qtm::SquareQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    SquareQuadFunction(
        centerMz - isolationWidthMz/T(2.0) - qtm.overhang,
        centerMz + isolationWidthMz/T(2.0) + qtm.overhang,
        centerMz
    )
end

function getQuadTransmissionBounds(sqm::SquareQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    return T(centerMz - isolationWidthMz/2), T(centerMz + isolationWidthMz/2)
end


function getPrecMinBound(sqf::SquareQuadFunction{T}) where {T<:AbstractFloat}
    return sqf.min_mz
end

function getPrecMaxBound(sqf::SquareQuadFunction{T}) where {T<:AbstractFloat}
    sqf.max_mz
end

function (sqf::SquareQuadFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    if ionMz > sqf.max_mz
        return zero(U)
    elseif ionMz < sqf.min_mz
        return zero(U)
    else
        return one(U)
    end
end