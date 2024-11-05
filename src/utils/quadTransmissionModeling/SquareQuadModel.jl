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