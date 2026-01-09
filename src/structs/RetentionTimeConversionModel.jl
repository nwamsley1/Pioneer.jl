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

#Uniform Cubic Smoothing Spline
#Given number of evenly spaced control points and data (time, value)
#solve coefficients for B-spline basis. Then build a speedy implementation. 
struct UniformSpline{N, T<:AbstractFloat} 
    coeffs::SVector{N, T}
    degree::Int64
    first::T
    last::T
    bin_width::T
end

abstract type RtConversionModel end
getModel(rt::RtConversionModel) = rt.model

struct SplineRtConversionModel <: RtConversionModel
    model::UniformSpline
end
(s::SplineRtConversionModel)(x::AbstractFloat) = s.model(x)

struct InterpolationRtConversionModel <: RtConversionModel
    model::Interpolations.Extrapolation
end
(s::InterpolationRtConversionModel)(x::AbstractFloat) = s.model(x)

struct LinearRtConversionModel <: RtConversionModel
    slope::Float32
    intercept::Float32
end

function (m::LinearRtConversionModel)(rt::AbstractFloat)
    return m.intercept + m.slope * rt
end

# Override getModel for LinearRtConversionModel since it doesn't wrap another model
getModel(m::LinearRtConversionModel) = m

struct IdentityModel <: RtConversionModel
    model::Function
    function IdentityModel()
        new(x::Float32 -> x::Float32)
    end
end

(i::IdentityModel)(x::AbstractFloat) = i.model(x)

RtConversionModel() = IdentityModel()