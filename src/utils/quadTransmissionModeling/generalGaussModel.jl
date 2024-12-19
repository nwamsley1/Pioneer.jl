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

