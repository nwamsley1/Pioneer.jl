# Fragment lookup tables and NCE models.
#
# A LibraryFragmentLookup maps precursor IDs to their fragment ions.
# Two concrete types:
# - StandardFragmentLookup: fixed intensity (DetailedFrag)
# - SplineFragmentLookup: spline-interpolated intensity (SplineDetailedFrag)
#
# NceModel predicts normalized collision energy from precursor m/z and charge.

abstract type LibraryFragmentLookup end

# ============================================================================
# StandardFragmentLookup (fixed intensity)
# ============================================================================

struct StandardFragmentLookup{T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{CompactFrag{T}}
    prec_frag_ranges::Vector{UInt64}
end

getFrag(lfp::StandardFragmentLookup, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::StandardFragmentLookup) = lfp.frags
getPrecFragRange(lfp::StandardFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} =
    range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))

function getSplineData(lfp::StandardFragmentLookup, prec_charge::UInt8, prec_mz::T) where {T<:AbstractFloat}
    return ConstantType()
end
function getSplineData(lfp::StandardFragmentLookup)
    return ConstantType()
end
function getSplineData(lfp::StandardFragmentLookup, nce_model::NceModel{T}, prec_charge::UInt8, prec_mz::T) where {T<:AbstractFloat}
    return ConstantType()
end
function getSplineData(lfp::StandardFragmentLookup, nce_model::NceModel{T}) where {T<:AbstractFloat}
    return ConstantType()
end

# ============================================================================
# NceModel — collision energy prediction
# ============================================================================

# NceModel{T} is defined in fragment_types.jl

"""
    PiecewiseNceModel{T}

Piecewise-linear NCE model with charge dependence.
- When x ≤ breakpoint: f(x,z) = left_slope * x + left_intercept + charge_slope * z
- When x > breakpoint: f(x,z) = right_value + charge_slope * z
"""
struct PiecewiseNceModel{T<:AbstractFloat} <: NceModel{T}
    breakpoint::T
    left_slope::T
    left_intercept::T
    right_value::T
    charge_slope::T
end

# ============================================================================
# SplineFragmentLookup (spline-interpolated intensity)
# ============================================================================

if !@isdefined(SplineFragmentLookup)
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineCompactFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    degree::Int64
end
end

getDegree(lfp::SplineFragmentLookup) = lfp.degree
getKnots(lfp::SplineFragmentLookup) = lfp.knots
getFrag(lfp::SplineFragmentLookup, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::SplineFragmentLookup) = lfp.frags
getPrecFragRange(lfp::SplineFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} =
    range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))

function getSplineData(lfp::SplineFragmentLookup{N,M,T}, nce_model::NceModel{T}, prec_charge::UInt8, prec_mz::T) where {N,M,T<:AbstractFloat}
    return SplineType(getKnots(lfp), nce_model(prec_mz, prec_charge), getDegree(lfp))
end

function getSplineData(lfp::SplineFragmentLookup{N,M,T}, nce_model::NceModel{T}) where {N,M,T<:AbstractFloat}
    return SplineType(getKnots(lfp), nce_model(), getDegree(lfp))
end
