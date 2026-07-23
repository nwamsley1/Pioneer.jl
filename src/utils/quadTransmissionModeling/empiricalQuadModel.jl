# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Empirical LUT quadrupole transmission model for Sciex ZT scanning DIA.
#
# The ZT physical quad is a ~5 Da window swept across m/z in ~1 Da steps, so the
# *effective* per-scan transmission is far wider than the recorded isolation width
# and is not well described by a square/Gaussian. QuadTuningSearch measures it from
# overlapping-scan fragment envelopes and stores it here as uniform knots
# (peak-normalized, spanning ±half_range about offset 0). Evaluated by linear
# interpolation on the signed offset (ionMz − center_mz); zero outside the knots.
#
# `values` is stored centered: length must be odd, values[(n+1)/2] is offset 0.

struct EmpiricalQuadModel{T<:AbstractFloat} <: QuadTransmissionModel
    values::Vector{T}   # transmission at uniform knots, peak-normalized, centered (odd length)
    step::T             # knot spacing (Da)
    bound_hw::T         # half-width (Da) for isotope-set admission (model support)
end

struct EmpiricalQuadFunction{T<:AbstractFloat} <: QuadTransmissionFunction
    min_mz::T
    max_mz::T
    center_mz::T
    values::Vector{T}
    step::T
    half_range::T       # (length(values)-1)/2 * step
end

function getQuadTransmissionFunction(m::EmpiricalQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    half_range = T((length(m.values) - 1) / 2) * m.step
    EmpiricalQuadFunction(
        centerMz - m.bound_hw,
        centerMz + m.bound_hw,
        centerMz,
        m.values,
        m.step,
        half_range,
    )
end

# Isotope-set admission bounds: the model's own support (where T is non-negligible).
function getQuadTransmissionBounds(m::EmpiricalQuadModel{T}, centerMz::T, isolationWidthMz::T) where {T<:AbstractFloat}
    return (centerMz - m.bound_hw, centerMz + m.bound_hw)
end

getPrecMinBound(f::EmpiricalQuadFunction{T}) where {T<:AbstractFloat} = f.min_mz
getPrecMaxBound(f::EmpiricalQuadFunction{T}) where {T<:AbstractFloat} = f.max_mz

@inline function (f::EmpiricalQuadFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    x = ionMz - f.center_mz
    (x <= -f.half_range || x >= f.half_range) && return zero(U)
    fi = (x + f.half_range) / f.step        # position in [0, n-1]
    lo = floor(Int, fi)
    frac = fi - lo
    i0 = lo + 1
    n = length(f.values)
    i0 >= n && return U(f.values[n])
    return U((one(U) - frac) * f.values[i0] + frac * f.values[i0 + 1])
end
