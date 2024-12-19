#A Quad model that always gives a precursor transmssion efficienty of 100%

struct NoQuadModel{T<:AbstractFloat} <: QuadTransmissionModel
    overhang::T
end

struct NoQuadFunction{T<:AbstractFloat} <: QuadTransmissionFunction 
    min_mz::T
    max_mz::T
    center_mz::T
end


function getQuadTransmissionFunction(qtm::NoQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    NoQuadFunction(
        centerMz - isolationWidthMz/T(2.0) - qtm.overhang,
        centerMz + isolationWidthMz/T(2.0) + qtm.overhang,
        centerMz
    )
end

function getQuadTransmissionBounds(sqm::NoQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    return T(centerMz - isolationWidthMz/2), T(centerMz + isolationWidthMz/2)
end


function getPrecMinBound(sqf::NoQuadFunction{T}) where {T<:AbstractFloat}
    return sqf.min_mz
end

function getPrecMaxBound(sqf::NoQuadFunction{T}) where {T<:AbstractFloat}
    sqf.max_mz
end

function (sqf::NoQuadFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    return one(U)
end