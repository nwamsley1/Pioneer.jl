# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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

# Parametrized on the concrete extrapolation type: an abstract
# `model::Interpolations.Extrapolation` field made `s.model(x)` type-unstable,
# boxing a value on every RT->iRT conversion (millions per file).
struct InterpolationRtConversionModel{M<:Interpolations.Extrapolation} <: RtConversionModel
    model::M
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

struct IdentityModel <: RtConversionModel end

(::IdentityModel)(x::AbstractFloat) = x
getModel(::IdentityModel) = identity

RtConversionModel() = IdentityModel()