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

struct IdentityModel <: RtConversionModel
    model::Function
    function IdentityModel()
        new(x::Float32 -> x::Float32)
    end
end

(i::IdentityModel)(x::AbstractFloat) = i.model(x)

RtConversionModel() = IdentityModel()