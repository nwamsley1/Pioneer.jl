# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Precomputed linear extrapolation at the boundaries of a spline.
Outside [lo, hi] (configured from training-data quantiles),
evaluation switches from the spline to: val + slope * (t - boundary).
This avoids trusting the spline in sparse tail regions while keeping
the full data for fitting.
"""
struct SplineExtrap{T<:AbstractFloat}
    lo::T        # lower extrapolation boundary
    hi::T        # upper extrapolation boundary
    lo_val::T    # spline(lo)
    lo_slope::T  # spline'(lo)
    hi_val::T    # spline(hi)
    hi_slope::T  # spline'(hi)
end

"""
Intensity-aware mass error model.

Bias is decomposed into m/z- and intensity-dependent smoothing splines (Da):
  bias_Da(mz, I, rt) =
      mz_bias_spline(mz) + intensity_bias_spline(log2(I)) + rt_bias_spline(rt)

Spread is Laplace-distributed with scale parameter **in mDa** modeled as a
regularized cubic B-spline of log2(I), then adjusted by an m/z-dependent
correction factor.

Bias splines are fit on all training data but use configured extrapolation
outside their supported ranges. `conservative_tol_da` is the stopping bound
for 2-arg methods. The active tolerance is determined by `k * σ_mDa / 1e3`.
"""
struct IntensityMassErrorModel{Nmz, NI, Nspread, Nrt, Nmzspread, T<:AbstractFloat} <: AbstractMassErrorModel
    # m/z-dependent bias (Da): binned robust regularized spline
    mz_bias_spline::UniformSpline{Nmz, T}

    # Intensity-dependent bias (Da): convexity-constrained smoothing spline of log2(I)
    intensity_bias_spline::UniformSpline{NI, T}

    # Laplace spread (mDa): regularized spline of log2(I)
    spread_spline::UniformSpline{Nspread, T}

    # RT-dependent additive bias (Da), evaluated on raw RT within this run.
    rt_bias_spline::UniformSpline{Nrt, T}

    # Extrapolation rules outside supported spline ranges
    mz_bias_extrap::SplineExtrap{T}
    int_bias_extrap::SplineExtrap{T}
    spread_extrap::SplineExtrap{T}
    rt_bias_extrap::SplineExtrap{T}

    # Coverage multiplier (k * Laplace_scale_mDa / 1e3 gives half-width in Da)
    k::T

    # Maximum tolerance in Da. Used as conservative stopping bound for 2-arg methods.
    conservative_tol_da::T

    # m/z-dependent spread correction: effective_spread_mDa = b_mda(log2I) * c(mz).
    # Coefficients are retained for logging; evaluation uses mz_spread_spline
    # with constant endpoint extrapolation.
    mz_spread_spline::UniformSpline{Nmzspread, T}
    mz_spread_extrap::SplineExtrap{T}
    mz_spread_α::T
    mz_spread_β::T
    mz_spread_γ::T

    # Default Laplace log-density for missing top-3 fragments.
    # Computed once at fit time: log f at 99th percentile tail using median intensity/mz.
    default_top3_ll::T

    # 99th percentile of |bias-corrected Da error| from training fragments.
    # Diagnostic only — NOT used as a clamp. Tolerance is determined by the
    # Gaussian model (k * σ_mDa / 1e3) with no upper bound.
    max_da_tolerance::T

    # Raw RT range covered by rt_bias_spline.
    rt_min::T
    rt_max::T
end

#==========================================================
Fast log2 approximation (bit manipulation, Float32 only)
==========================================================#

@inline function fast_log2(x::Float32)::Float32
    i = reinterpret(Int32, x)
    return Float32(i) * 1.1920929f-7 - 126.94269504f0
end

#==========================================================
Spline derivative (for precomputing extrapolation slopes)
==========================================================#

# Derivative of a UniformSpline at point t (clamped to [first, last]).
# For cubic: d/dt [c0 + c1*u + c2*u² + c3*u³] = (c1 + 2*c2*u + 3*c3*u²) / bin_width
@inline function spline_derivative(s::UniformSpline, t::U) where {U<:AbstractFloat}
    @fastmath @inbounds begin
        t = max(min(t, s.last), s.first)
        idx = floor(Int32, (t - s.first) / s.bin_width)
        u = (t - (s.first + s.bin_width * idx)) / s.bin_width
        coeff = idx * (s.degree + 1) + 1
        c1 = U(s.coeffs[coeff + 1])
        c2 = U(s.coeffs[coeff + 2])
        c3 = U(s.coeffs[coeff + 3])
        deriv_du = muladd(U(3) * c3, u, U(2) * c2)
        deriv_du = muladd(deriv_du, u, c1)
        return deriv_du / U(s.bin_width)
    end
end

# Build SplineExtrap from a spline and quantile boundaries.
function make_spline_extrap(spline::UniformSpline{N, T}, lo::T, hi::T) where {N, T}
    SplineExtrap{T}(
        lo, hi,
        spline(lo), spline_derivative(spline, lo),
        spline(hi), spline_derivative(spline, hi)
    )
end

#==========================================================
Spline evaluation with linear extrapolation outside stored quantile boundaries
==========================================================#

@inline function _eval_with_extrap(spline::UniformSpline, e::SplineExtrap{T}, t::U) where {T, U<:AbstractFloat}
    @fastmath begin
        if t < U(e.lo)
            return muladd(U(e.lo_slope), t - U(e.lo), U(e.lo_val))
        elseif t > U(e.hi)
            return muladd(U(e.hi_slope), t - U(e.hi), U(e.hi_val))
        else
            return U(spline(t))
        end
    end
end

#==========================================================
Bias computation (all in Da, via spline evaluation)
==========================================================#

@inline function _mz_bias_da(m::IntensityMassErrorModel, mz::T) where {T<:AbstractFloat}
    return _eval_with_extrap(m.mz_bias_spline, m.mz_bias_extrap, mz)
end

@inline function _intensity_bias_da(m::IntensityMassErrorModel, log2I::T) where {T<:AbstractFloat}
    return _eval_with_extrap(m.intensity_bias_spline, m.int_bias_extrap, log2I)
end

@inline function _rt_bias_da(m::IntensityMassErrorModel, rt::T) where {T<:AbstractFloat}
    if !isfinite(rt)
        rt = (T(m.rt_min) + T(m.rt_max)) / T(2)
    end
    return _eval_with_extrap(m.rt_bias_spline, m.rt_bias_extrap, rt)
end

@inline function _total_bias_da(m::IntensityMassErrorModel, mz::T, log2I::T) where {T<:AbstractFloat}
    return _mz_bias_da(m, mz) + _intensity_bias_da(m, log2I)
end

@inline function _total_bias_da(m::IntensityMassErrorModel, mz::T, log2I::T, rt::T) where {T<:AbstractFloat}
    return _total_bias_da(m, mz, log2I) + _rt_bias_da(m, rt)
end

#==========================================================
Spread (Laplace scale in mDa) — spline with linear extrapolation
==========================================================#

@inline function _laplace_scale_mda(m::IntensityMassErrorModel, log2I::T) where {T<:AbstractFloat}
    scale = _eval_with_extrap(m.spread_spline, m.spread_extrap, log2I)
    return max(scale, T(1e-4))  # floor at 0.0001 mDa
end

# Convert Laplace scale (median/ln2) to Gaussian σ (median/0.6745)
# σ = b_laplace * ln(2) / 0.6745 ≈ 1.0269 * b_laplace
const _LAPLACE_TO_GAUSSIAN = Float32(log(2.0) / 0.6744897501960817)

@inline function _gaussian_sigma_mda(m::IntensityMassErrorModel, log2I::T) where {T<:AbstractFloat}
    return _laplace_scale_mda(m, log2I) * T(_LAPLACE_TO_GAUSSIAN)
end

#==========================================================
m/z-dependent spread correction factor (in mDa space)
==========================================================#

@inline function _mz_spread_correction(spline::UniformSpline, e::SplineExtrap{T}, mz::U) where {T, U<:AbstractFloat}
    @fastmath return min(max(_eval_with_extrap(spline, e, mz), U(0.1)), U(10.0))
end

@inline function _mz_spread_correction(m::IntensityMassErrorModel, mz::T) where {T<:AbstractFloat}
    return _mz_spread_correction(m.mz_spread_spline, m.mz_spread_extrap, mz)
end

#==========================================================
2-arg interface (conservative Da bounds — ignores intensity)
==========================================================#

function getCorrectedMz(m::IntensityMassErrorModel, mz::Float32)
    bias_da = _mz_bias_da(m, mz)
    return Float32(mz - bias_da)
end

function getMzBoundsReverse(m::IntensityMassErrorModel, mass::Float32)
    hw = m.conservative_tol_da
    return Float32(mass - hw), Float32(mass + hw)
end

function getMzBounds(m::IntensityMassErrorModel, mass::Float32)
    hw = m.conservative_tol_da
    return Float32(mass - hw), Float32(mass + hw)
end

getRightTol(m::IntensityMassErrorModel) = m.conservative_tol_da
getLeftTol(m::IntensityMassErrorModel) = m.conservative_tol_da
getMassOffset(m::IntensityMassErrorModel) = Float32(m.mz_bias_spline(m.mz_bias_spline.first))
getMassCorrection(m::IntensityMassErrorModel) = getMassOffset(m)

#==========================================================
Intensity/RT-aware interface used by fragment matching.
Spread is in mDa; converted to Da directly.
==========================================================#

function getCorrectedMz(m::IntensityMassErrorModel, mz::Float32, intensity::Float32)
    log2I = fast_log2(intensity)
    bias_da = _total_bias_da(m, mz, log2I)
    return Float32(mz - bias_da)
end

function getCorrectedMz(m::IntensityMassErrorModel, mz::Float32, intensity::Float32, rt::Float32)
    log2I = fast_log2(intensity)
    bias_da = _total_bias_da(m, mz, log2I, rt)
    return Float32(mz - bias_da)
end

function getMzBoundsReverse(m::IntensityMassErrorModel, mass::Float32, log2I::Float32)
    σ_mda = _gaussian_sigma_mda(m, log2I) * _mz_spread_correction(m, mass)
    half_width = m.k * σ_mda / 1f3  # mDa → Da
    return Float32(mass - half_width), Float32(mass + half_width)
end

getMzBoundsReverse(m::IntensityMassErrorModel, mass::Float32, log2I::Float32, ::Float32) =
    getMzBoundsReverse(m, mass, log2I)

"""
    getCorrectedMzAndBounds(m, mz, intensity) -> (corrected, low, high)

Bias-correct an observed peak and compute its Da tolerance window in one call.
Returns the corrected m/z and the range of theoretical values that could match.
Avoids redundant log2 computation vs calling getCorrectedMz + getMzBoundsReverse separately.
"""
@inline function getCorrectedMzAndBounds(m::IntensityMassErrorModel, mz::Float32, intensity::Float32)
    @fastmath begin
        log2I = fast_log2(intensity)
        bias_da = _total_bias_da(m, mz, log2I)
        corrected = mz - bias_da
        σ_mda = _gaussian_sigma_mda(m, log2I) * _mz_spread_correction(m, corrected)
        half_width = m.k * σ_mda * 1f-3  # mDa → Da (multiply faster than divide)
        return corrected, corrected - half_width, corrected + half_width
    end
end

@inline function getCorrectedMzAndBounds(m::IntensityMassErrorModel, mz::Float32, intensity::Float32, rt::Float32)
    @fastmath begin
        log2I = fast_log2(intensity)
        bias_da = _total_bias_da(m, mz, log2I, rt)
        corrected = mz - bias_da
        σ_mda = _gaussian_sigma_mda(m, log2I) * _mz_spread_correction(m, corrected)
        half_width = m.k * σ_mda * 1f-3
        return corrected, corrected - half_width, corrected + half_width
    end
end

#==========================================================
Laplace log-density for fragment mass error scoring
==========================================================#

"""
    laplace_log_density(m, theo_mz, obs_mz, intensity) -> Float32

Compute the Gaussian log-density of a fragment mass error given the model.
Returns log f(x | σ) = -0.5·ln(2πσ²) - x²/(2σ²) where x is the bias-corrected
Da error and σ is the intensity-aware Gaussian scale in Da (= σ_mDa / 1e3).
(Function name kept as laplace_log_density for API compatibility.)
"""
@inline function laplace_log_density(m::IntensityMassErrorModel, theo_mz::Float32, obs_mz::Float32, intensity::Float32)
    log2I = fast_log2(intensity)
    bias_da = _total_bias_da(m, obs_mz, log2I)
    corrected_err = Float64(obs_mz - bias_da) - Float64(theo_mz)
    σ_mda = Float64(_gaussian_sigma_mda(m, log2I) * _mz_spread_correction(m, theo_mz))
    σ_da = max(σ_mda / 1e3, 1e-12)
    return Float32(-0.5 * log(2π * σ_da^2) - corrected_err^2 / (2.0 * σ_da^2))
end

# SimpleMassErrorModel fallback: no intensity model, return 0
laplace_log_density(::SimpleMassErrorModel, ::Float32, ::Float32, ::Float32) = Float32(0)

# Default log-density for missing top-3 fragments
get_default_top3_ll(m::IntensityMassErrorModel) = m.default_top3_ll
get_default_top3_ll(::SimpleMassErrorModel) = Float32(0)

#==========================================================
Spline conversion helper
==========================================================#

function convert_spline_to_f32(s::UniformSpline{N, Float64}) where N
    return UniformSpline{N, Float32}(
        SVector{N, Float32}(Float32.(s.coeffs)),
        s.degree,
        Float32(s.first),
        Float32(s.last),
        Float32(s.bin_width)
    )
end

convert_spline_extrap_to_f32(e::SplineExtrap{Float64}) = SplineExtrap{Float32}(
    Float32(e.lo), Float32(e.hi),
    Float32(e.lo_val), Float32(e.lo_slope),
    Float32(e.hi_val), Float32(e.hi_slope)
)

#==========================================================
ScoutCalibratedMassErrorModel

Wide-scout calibration model: spline-based m/z + intensity bias
correction (in Da) with a CONSTANT mDa tolerance.

Tier 1 (≥1000 PSMs at 1% FDR): full m/z + intensity spline bias.
Tier 2 (<1000 PSMs): set has_intensity_bias=false; intensity spline
  is present but zeroed — only m/z bias is active.

The constant tolerance_da is derived from the MAD of corrected
residuals: tolerance_da = MAD/0.6745 × 2.576 × 1.5 / 1e3.
==========================================================#

struct ScoutCalibratedMassErrorModel{Nmz, NI, T<:AbstractFloat} <: AbstractMassErrorModel
    mz_bias_spline::UniformSpline{Nmz, T}
    mz_bias_extrap::SplineExtrap{T}
    intensity_bias_spline::UniformSpline{NI, T}
    int_bias_extrap::SplineExtrap{T}
    has_intensity_bias::Bool
    tolerance_da::T
end

# ── Bias computation ──

@inline function _scout_mz_bias_da(m::ScoutCalibratedMassErrorModel, mz::T) where {T<:AbstractFloat}
    return _eval_with_extrap(m.mz_bias_spline, m.mz_bias_extrap, mz)
end

@inline function _scout_intensity_bias_da(m::ScoutCalibratedMassErrorModel, log2I::T) where {T<:AbstractFloat}
    m.has_intensity_bias || return zero(T)
    return _eval_with_extrap(m.intensity_bias_spline, m.int_bias_extrap, log2I)
end

@inline function _scout_total_bias_da(m::ScoutCalibratedMassErrorModel, mz::T, log2I::T) where {T<:AbstractFloat}
    return _scout_mz_bias_da(m, mz) + _scout_intensity_bias_da(m, log2I)
end

# ── 2-arg interface (no intensity — conservative, mz-only bias) ──

function getCorrectedMz(m::ScoutCalibratedMassErrorModel, mz::Float32)
    bias_da = _scout_mz_bias_da(m, mz)
    return Float32(mz - bias_da)
end

function getMzBoundsReverse(m::ScoutCalibratedMassErrorModel, mass::Float32)
    return Float32(mass - m.tolerance_da), Float32(mass + m.tolerance_da)
end

function getMzBounds(m::ScoutCalibratedMassErrorModel, mass::Float32)
    return Float32(mass - m.tolerance_da), Float32(mass + m.tolerance_da)
end

getRightTol(m::ScoutCalibratedMassErrorModel) = m.tolerance_da
getLeftTol(m::ScoutCalibratedMassErrorModel) = m.tolerance_da
getMassOffset(m::ScoutCalibratedMassErrorModel) = Float32(m.mz_bias_spline(m.mz_bias_spline.first))
getMassCorrection(m::ScoutCalibratedMassErrorModel) = getMassOffset(m)

# ── 3-arg interface (intensity-aware bias, constant tolerance) ──

function getCorrectedMz(m::ScoutCalibratedMassErrorModel, mz::Float32, intensity::Float32)
    log2I = fast_log2(intensity)
    bias_da = _scout_total_bias_da(m, mz, log2I)
    return Float32(mz - bias_da)
end

getMzBoundsReverse(m::ScoutCalibratedMassErrorModel, mass::Float32, ::Float32) = getMzBoundsReverse(m, mass)
getCorrectedMz(m::ScoutCalibratedMassErrorModel, mz::Float32, intensity::Float32, ::Float32) =
    getCorrectedMz(m, mz, intensity)
getMzBoundsReverse(m::ScoutCalibratedMassErrorModel, mass::Float32, log2I::Float32, ::Float32) =
    getMzBoundsReverse(m, mass, log2I)
laplace_log_density(::ScoutCalibratedMassErrorModel, ::Float32, ::Float32, ::Float32) = Float32(0)
get_default_top3_ll(::ScoutCalibratedMassErrorModel) = Float32(0)

@inline function getCorrectedMzAndBounds(m::ScoutCalibratedMassErrorModel, mz::Float32, intensity::Float32)
    @fastmath begin
        log2I = fast_log2(intensity)
        bias_da = _scout_total_bias_da(m, mz, log2I)
        corrected = mz - bias_da
        return corrected, Float32(corrected - m.tolerance_da), Float32(corrected + m.tolerance_da)
    end
end

@inline function getCorrectedMzAndBounds(m::ScoutCalibratedMassErrorModel, mz::Float32, intensity::Float32, ::Float32)
    return getCorrectedMzAndBounds(m, mz, intensity)
end
