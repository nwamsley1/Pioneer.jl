
abstract type RtConversionModel end
getModel(rt::RtConversionModel) = rt.model

struct SplineRtConversionModel <: RtConversionModel
    model::UniformSpline
end

struct IdentityModel <: RtConversionModel
    model::Function
    function IdentityModel()
        new(x::Float32 -> x::Float32)
    end
end

RtConversionModel() = IdentityModel()