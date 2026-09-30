# Fragment lookup tables and NCE models.
#
# A LibraryFragmentLookup maps precursor IDs to their fragment ions.
# Two concrete types:
# - StandardFragmentLookup: fixed intensity (DetailedFrag)
# - SplineFragmentLookup: spline-interpolated intensity (SplineCompactFrag)
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

struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineCompactFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
end

getKnots(lfp::SplineFragmentLookup) = lfp.knots
getFrag(lfp::SplineFragmentLookup, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::SplineFragmentLookup) = lfp.frags
getPrecFragRange(lfp::SplineFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} =
    range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))

# ============================================================================
# Prepared spline interpolation data
# ============================================================================

struct ConstantSplineIntensityModel{T<:AbstractFloat}
    data::PreparedSplineFractions{T}
end

struct DynamicSplineIntensityModel{M,K}
    nce_model::M
    knots::K
end

struct BinnedSplineIntensityModel{M,T<:AbstractFloat}
    nce_model::M
    data::Vector{PreparedSplineFractions{T}}
    default_data::PreparedSplineFractions{T}
end

prepare_fragment_intensity_model(
    ::StandardFragmentLookup, ::NceModel) = ConstantType()

function prepare_fragment_intensity_model(
        lookup::SplineFragmentLookup, nce_model::NceModel)
    return DynamicSplineIntensityModel(nce_model, getKnots(lookup))
end

# The fifth argument is the scan's collision energy in eV (0 when unknown). Only
# collision-energy-keyed NCE models (CeBinnedNceModel) use it; the four-argument
# form passes 0 for callers without a scan.
@inline getSplineData(lookup::LibraryFragmentLookup, model, prec_charge::UInt8, prec_mz::AbstractFloat) =
    getSplineData(lookup, model, prec_charge, prec_mz, 0f0)

@inline getSplineData(
    ::StandardFragmentLookup, intensity_data::ConstantType,
    ::UInt8, ::AbstractFloat, ::AbstractFloat) = intensity_data

@inline getSplineData(
    ::SplineFragmentLookup, model::ConstantSplineIntensityModel,
    ::UInt8, ::AbstractFloat, ::AbstractFloat) = model.data

@inline function getSplineData(
    ::SplineFragmentLookup, model::DynamicSplineIntensityModel,
        prec_charge::UInt8, prec_mz::AbstractFloat, scan_ev::AbstractFloat)
    nce = model.nce_model(prec_mz, prec_charge, scan_ev)
    return prepare_spline_fractions(nce, model.knots)
end

@inline function getSplineData(
        ::SplineFragmentLookup, model::BinnedSplineIntensityModel,
        prec_charge::UInt8, prec_mz::AbstractFloat, scan_ev::AbstractFloat)
    slot = nce_cache_slot(model.nce_model, prec_mz, prec_charge, scan_ev)
    return slot == 0 ? model.default_data : @inbounds(model.data[slot])
end
