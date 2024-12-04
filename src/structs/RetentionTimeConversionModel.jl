
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