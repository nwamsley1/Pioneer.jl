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

#==========================================================
Intensity-Aware Mass Error Model Fitting Pipeline

Fits an IntensityMassErrorModel from matched fragment data:
  1. m/z-dependent bias: binned robust regularized spline (Da)
  2. intensity-dependent bias: convexity-constrained smoothing spline (Da)
  3. RT-dependent bias: binned robust regularized spline on raw RT (Da)
  4. Laplace spread: regularized spline of log2(I) (mDa)
  5. m/z-dependent spread correction
  6. Empirical coverage multiplier k

Bias is fit in Da space. Spread is fit in mDa space after bias correction,
then adjusted by a fitted m/z-dependent correction factor.
==========================================================#

# Empirical-k coverage parameters. Quantile of |residual|/σ used to set the
# coverage multiplier k; the result is clamped to [EMPK_MIN_K, EMPK_CLAMP_HI].
const EMPK_QUANTILE = 0.95
const EMPK_MIN_K = Float32(quantile(Normal(), (1.0 + EMPK_QUANTILE) / 2.0))
const EMPK_CLAMP_HI = 4.5f0

#==========================================================
Binned Laplace scale estimator
==========================================================#

"""
    binned_laplace_scale(order, vals, errors, bin_sz)

Compute Laplace scale b = median(|err - median(err)|) / ln(2) per quantile bin.
Returns (centers, b_values) vectors.
"""
function binned_laplace_scale(order::Vector{Int}, vals::Vector{Float64},
                               errors::Vector{Float64}, bin_sz::Int)
    ntot = length(order)
    nb = max(ntot ÷ bin_sz, 1)
    centers = Float64[]
    b_vals = Float64[]
    for b in 1:nb
        s = (b-1)*bin_sz + 1
        e = b == nb ? ntot : b*bin_sz
        idx = order[s:e]
        errs = errors[idx]
        mu = median(errs)
        push!(centers, median(vals[idx]))
        push!(b_vals, median(abs.(errs .- mu)) / log(2))
    end
    return centers, b_vals
end

function _binned_bias_medians(x::AbstractVector{<:Real},
                              y::AbstractVector{<:Real},
                              lo::Real,
                              hi::Real;
                              bin_size::Int = TUNING_CALIBRATION_BIN_SIZE)
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))
    lo_f = Float64(lo)
    hi_f = Float64(hi)
    valid = Int[]
    for i in eachindex(x)
        xi = Float64(x[i])
        yi = Float64(y[i])
        if isfinite(xi) && isfinite(yi) && lo_f <= xi <= hi_f
            push!(valid, Int(i))
        end
    end
    length(valid) < 2 && return Float64[], Float64[]

    order = valid[sortperm(Float64.(x[valid]))]
    n = length(order)
    n_bins = max(n ÷ bin_size, 1)
    centers = Float64[]
    vals = Float64[]
    sizehint!(centers, n_bins)
    sizehint!(vals, n_bins)

    for b in 1:n_bins
        s = (b - 1) * bin_size + 1
        e = b == n_bins ? n : b * bin_size
        idx = order[s:e]
        push!(centers, median(Float64.(x[idx])))
        push!(vals, median(Float64.(y[idx])))
    end

    return centers, vals
end

function _theil_sen_slope(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))
    length(x) < 2 && return nothing

    slopes = Float64[]
    for j in 2:length(x)
        for i in 1:(j - 1)
            dx = x[j] - x[i]
            abs(dx) <= eps(Float64) && continue
            push!(slopes, (y[j] - y[i]) / dx)
        end
    end

    isempty(slopes) && return nothing
    return median(slopes)
end

function _edge_point_slope(centers::Vector{Float64},
                           vals::Vector{Float64},
                           side::Symbol,
                           n_edge_points::Int)
    n = length(centers)
    n < 2 && return nothing
    n_take = min(max(n_edge_points, 2), n)
    if side === :lo
        idx = 1:n_take
    elseif side === :hi
        idx = (n - n_take + 1):n
    else
        throw(ArgumentError("side must be :lo or :hi"))
    end
    return _theil_sen_slope(centers[idx], vals[idx])
end

function make_edge_point_spline_extrap(spline::UniformSpline{N, T},
                                       x::AbstractVector{<:Real},
                                       y::AbstractVector{<:Real},
                                       lo::T,
                                       hi::T;
                                       n_edge_points::Int = TUNING_BIAS_EXTRAP_EDGE_POINTS,
                                       bin_size::Int = TUNING_CALIBRATION_BIN_SIZE) where {N, T}
    centers, vals = _binned_bias_medians(x, y, lo, hi; bin_size=bin_size)
    lo_slope = _edge_point_slope(centers, vals, :lo, n_edge_points)
    hi_slope = _edge_point_slope(centers, vals, :hi, n_edge_points)
    extrap_lo = length(centers) >= 2 ? T(first(centers)) : lo
    extrap_hi = length(centers) >= 2 ? T(last(centers)) : hi

    return SplineExtrap{T}(
        extrap_lo, extrap_hi,
        spline(extrap_lo), T(lo_slope === nothing ? 0.0 : lo_slope),
        spline(extrap_hi), T(hi_slope === nothing ? 0.0 : hi_slope)
    )
end

function make_binned_bound_spline_extrap(spline::UniformSpline{N, T},
                                         x::AbstractVector{<:Real},
                                         y::AbstractVector{<:Real},
                                         lo::T,
                                         hi::T;
                                         bin_size::Int = TUNING_CALIBRATION_BIN_SIZE) where {N, T}
    centers, _ = _binned_bias_medians(x, y, lo, hi; bin_size=bin_size)
    extrap_lo = length(centers) >= 2 ? T(first(centers)) : lo
    extrap_hi = length(centers) >= 2 ? T(last(centers)) : hi
    return make_spline_extrap(spline, extrap_lo, extrap_hi)
end

function make_rt_range_edge_spline_extrap(spline::UniformSpline{N, T},
                                          x::AbstractVector{<:Real},
                                          y::AbstractVector{<:Real},
                                          lo::T,
                                          hi::T;
                                          bound_bin_size::Int = TUNING_RT_BIAS_BIN_SIZE) where {N, T}
    bound_centers, _ = _binned_bias_medians(x, y, lo, hi; bin_size=bound_bin_size)
    lo_f = length(bound_centers) >= 2 ? first(bound_centers) : Float64(lo)
    hi_f = length(bound_centers) >= 2 ? last(bound_centers) : Float64(hi)
    span = hi_f - lo_f
    if !(isfinite(span) && span > 0)
        return SplineExtrap{T}(
            lo, hi,
            spline(lo), zero(T),
            spline(hi), zero(T)
        )
    end

    return SplineExtrap{T}(
        T(lo_f), T(hi_f),
        spline(T(lo_f)), zero(T),
        spline(T(hi_f)), zero(T)
    )
end

function _constant_mz_spread_correction_spline(first_mz::Float64,
                                                last_mz::Float64;
                                                value::Float64 = 1.0)
    _first = isfinite(first_mz) ? first_mz : 0.0
    _last = isfinite(last_mz) ? last_mz : _first + 1.0
    _last <= _first && (_last = _first + 1.0)
    value = clamp(value, 0.1, 10.0)
    degree = 3
    n_knots = 3
    bin_width = (_last - _first) / (n_knots - 1)
    coeffs = Float64[]
    for _ in 1:n_knots
        append!(coeffs, (value, 0.0, 0.0, 0.0))
    end
    return UniformSpline{length(coeffs), Float64}(
        SVector{length(coeffs), Float64}(coeffs),
        degree,
        _first,
        _last,
        bin_width
    )
end

function _constant_mz_spread_correction_extrap(spline::UniformSpline{N, Float64}) where N
    return SplineExtrap{Float64}(
        spline.first,
        spline.last,
        spline(spline.first),
        0.0,
        spline(spline.last),
        0.0
    )
end

function _constant_bias_spline(first_x::Float64,
                               last_x::Float64;
                               value::Float64 = 0.0)
    _first = isfinite(first_x) ? first_x : 0.0
    _last = isfinite(last_x) ? last_x : _first + 1.0
    _last <= _first && (_last = _first + 1.0)
    degree = 3
    n_knots = 3
    coeffs = Float64[]
    for _ in 1:n_knots
        append!(coeffs, (value, 0.0, 0.0, 0.0))
    end
    return UniformSpline{length(coeffs), Float64}(
        SVector{length(coeffs), Float64}(coeffs),
        degree,
        _first,
        _last,
        (_last - _first) / (n_knots - 1)
    )
end

function _constant_bias_extrap(spline::UniformSpline{N, Float64}) where N
    return SplineExtrap{Float64}(
        spline.first,
        spline.last,
        spline(spline.first),
        0.0,
        spline(spline.last),
        0.0
    )
end

function _shift_spline(spline::UniformSpline{N, Float64}, offset::Float64) where N
    coeffs = collect(spline.coeffs)
    stride = spline.degree + 1
    for i in 1:stride:length(coeffs)
        coeffs[i] -= offset
    end
    return UniformSpline{N, Float64}(
        SVector{N, Float64}(coeffs),
        spline.degree,
        spline.first,
        spline.last,
        spline.bin_width
    )
end

function _linear_mz_spread_correction_spline(first_mz::Float64,
                                             last_mz::Float64,
                                             intercept::Float64,
                                             slope::Float64)
    _first = isfinite(first_mz) ? first_mz : 0.0
    _last = isfinite(last_mz) ? last_mz : _first + 1.0
    _last <= _first && (_last = _first + 1.0)

    degree = 3
    n_knots = 3
    bin_width = (_last - _first) / (n_knots - 1)
    coeffs = Float64[]
    for k in 0:(n_knots - 1)
        mz = _first + k * bin_width
        append!(coeffs, (intercept + slope * mz, slope * bin_width, 0.0, 0.0))
    end

    return UniformSpline{length(coeffs), Float64}(
        SVector{length(coeffs), Float64}(coeffs),
        degree,
        _first,
        _last,
        bin_width
    )
end

function _fit_robust_linear_mz_spread_correction(
    corr_centers::Vector{Float64},
    corr_vals::Vector{Float64};
    floor_val::Float64 = 0.1,
    ceiling_val::Float64 = 10.0,
)
    n_pts = length(corr_centers)
    n_pts < 2 && return nothing, Float32(1.0), Float32(0.0), Float32(0.0),
                         "constant (too few m/z correction bins)"

    order = sortperm(corr_centers)
    x = corr_centers[order]
    y = clamp.(corr_vals[order], floor_val, ceiling_val)
    _first = minimum(x)
    _last = maximum(x)
    _last <= _first && return nothing, Float32(1.0), Float32(0.0), Float32(0.0),
                         "constant (degenerate m/z correction bins)"

    slope = _theil_sen_slope(x, y)
    slope === nothing && return nothing, Float32(1.0), Float32(0.0), Float32(0.0),
                              "constant (too few m/z correction bins)"
    intercept = median(y .- slope .* x)

    spline = _linear_mz_spread_correction_spline(_first, _last, intercept, slope)
    return spline, Float32(intercept), Float32(slope), Float32(0.0),
           "robust linear (Theil-Sen, bins=$(n_pts))"
end

function _select_regularized_mz_spread_correction(corr_centers::AbstractVector{<:Real},
                                                  corr_vals::AbstractVector{<:Real})
    centers = Float64[]
    vals = Float64[]
    for (center, val) in zip(corr_centers, corr_vals)
        c = Float64(center)
        v = Float64(val)
        if isfinite(c) && isfinite(v) && v > 0.0
            push!(centers, c)
            push!(vals, v)
        end
    end

    if isempty(centers)
        spline = _constant_mz_spread_correction_spline(0.0, 1.0)
        return spline, _constant_mz_spread_correction_extrap(spline),
               Float32(1.0), Float32(0.0), Float32(0.0), "constant (too few fragments)"
    end

    spline, mz_spread_α, mz_spread_β, mz_spread_γ, label =
        _fit_robust_linear_mz_spread_correction(centers, vals)
    if spline === nothing
        value = median(clamp.(vals, 0.1, 10.0))
        spline = _constant_mz_spread_correction_spline(minimum(centers), maximum(centers); value=value)
        mz_spread_α = Float32(value)
        mz_spread_β = Float32(0.0)
        mz_spread_γ = Float32(0.0)
    end

    return spline, _constant_mz_spread_correction_extrap(spline),
           mz_spread_α, mz_spread_β, mz_spread_γ, label
end

function _project_convex_coeffs!(c::Vector{Float64}, convex_up::Bool; floor_val::Union{Nothing, Float64}=nothing)
    for _ in 1:30
        changed = false
        for i in 1:(length(c)-2)
            v = c[i+2] - 2*c[i+1] + c[i]
            if convex_up ? (v < 0) : (v > 0)
                correction = abs(v)
                if convex_up
                    c[i] += correction / 3
                    c[i+1] -= 2 * correction / 3
                    c[i+2] += correction / 3
                else
                    c[i] -= correction / 3
                    c[i+1] += 2 * correction / 3
                    c[i+2] -= correction / 3
                end
                changed = true
            end
        end
        !changed && break
    end

    if floor_val !== nothing
        for i in eachindex(c)
            c[i] = max(c[i], floor_val)
        end
    end
    return c
end

function _project_monotone_coeffs!(c::Vector{Float64}, increasing::Bool)
    for _ in 1:30
        changed = false
        for i in 1:(length(c)-1)
            if increasing ? (c[i+1] < c[i]) : (c[i+1] > c[i])
                avg = (c[i] + c[i+1]) / 2
                c[i] = avg
                c[i+1] = avg
                changed = true
            end
        end
        !changed && break
    end
    return c
end

function _projected_gradient!(
    c::Vector{Float64},
    H::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    project!::F;
    lr::Float64,
    max_iter::Int,
    tol::Float64,
    check_every::Int,
) where {F}
    project!(c)

    # opnorm(H) is the gradient Lipschitz constant of this convex quadratic.
    lr_eff = min(lr, 1.8 / opnorm(H))
    grad = similar(c)
    prev_obj = Inf
    for iter in 1:max_iter
        mul!(grad, H, c)
        grad .-= rhs
        c .-= lr_eff .* grad
        project!(c)
        if iter % check_every == 0
            mul!(grad, H, c)
            obj = 0.5 * dot(c, grad) - dot(rhs, c)
            if abs(prev_obj - obj) / max(abs(obj), 1.0) < tol
                break
            end
            prev_obj = obj
        end
    end
    return c
end

function _projected_irls_coeffs(
    X::AbstractMatrix{Float64},
    y::Vector{Float64},
    n_coeffs::Int,
    project!::F;
    λ::Float64,
    max_iter::Int,
    lr::Float64,
    tol::Float64,
    n_irls::Int,
) where {F}
    w = ones(length(y))
    c = zeros(n_coeffs)
    resid = similar(y)

    for irls_iter in 1:n_irls
        H, rhs = _penalized_spline_system(X, y; λ=λ, weights=w, order=2)
        c = H \ rhs
        _projected_gradient!(c, H, rhs, project!;
                             lr=lr, max_iter=max_iter, tol=tol, check_every=200)

        irls_iter == n_irls && break
        mul!(resid, X, c)
        resid .= y .- resid
        σ_hat = median(abs.(resid)) / 0.6744897501960817
        δ = 1.345 * max(σ_hat, 1e-10)
        for i in eachindex(w, resid)
            r_abs = abs(resid[i])
            w[i] = r_abs <= δ ? 1.0 : δ / r_abs
        end
    end
    return c
end

#==========================================================
Convexity-constrained smoothing spline for bias
==========================================================#

"""
    fit_convex_bias_spline(x, y; convex_up, n_knots, λ, max_iter, lr, tol)

Fit a smoothing spline to (x, y) data with a second-derivative sign constraint.

`convex_up=true`:  f''(x) >= 0 (bowl shape) — coeffs: c[i+2]-2c[i+1]+c[i] >= 0
`convex_up=false`: f''(x) <= 0 (dome shape) — coeffs: c[i+2]-2c[i+1]+c[i] <= 0

Fits to ALL data points via penalized LS + projected gradient descent.
Returns a UniformSpline (Float64), or nothing if fitting fails.
"""
function fit_convex_bias_spline(x::Vector{Float64}, y::Vector{Float64};
                                 convex_up::Bool = true,
                                 n_knots::Int = 10,
                                 λ::Float64 = 10.0,
                                 max_iter::Int = 5000,
                                 lr::Float64 = 0.0001,
                                 tol::Float64 = 1e-10,
                                 n_irls::Int = 3)
    n_pts = length(x)
    n_pts < n_knots && return nothing

    basis = _make_uniform_spline_basis(x, n_knots)
    X = basis.X
    n_coeffs = basis.n_coeffs

    project!(c) = _project_convex_coeffs!(c, convex_up)
    c = _projected_irls_coeffs(X, y, n_coeffs, project!;
                               λ=λ, max_iter=max_iter, lr=lr, tol=tol, n_irls=n_irls)

    return _coeffs_to_spline(c, basis)
end

"""
    fit_best_convex_bias(x, y; n_knots, λ, ...)

Fit both convex-up and convex-down splines, return the one with lower L1 loss.
Returns (spline, label) or (nothing, "failed").
"""
function fit_best_convex_bias(x::Vector{Float64}, y::Vector{Float64};
                               n_knots::Int = 10, λ::Float64 = 10.0,
                               max_iter::Int = 5000, lr::Float64 = 0.0001,
                               tol::Float64 = 1e-10)
    spline_up = fit_convex_bias_spline(x, y; convex_up=true, n_knots=n_knots,
                                        λ=λ, max_iter=max_iter, lr=lr, tol=tol)
    spline_dn = fit_convex_bias_spline(x, y; convex_up=false, n_knots=n_knots,
                                        λ=λ, max_iter=max_iter, lr=lr, tol=tol)

    loss_up = spline_up !== nothing ? sum(abs.(y .- [spline_up(xi) for xi in x])) : Inf
    loss_dn = spline_dn !== nothing ? sum(abs.(y .- [spline_dn(xi) for xi in x])) : Inf

    if loss_up <= loss_dn && spline_up !== nothing
        return spline_up, "convex-up"
    elseif spline_dn !== nothing
        return spline_dn, "convex-down"
    else
        return nothing, "failed"
    end
end

"""
    fit_binned_regularized_bias_spline(x, y; n_knots, λ, bin_size)

Fit an unconstrained cubic smoothing spline to equal-count m/z bins. Each bin
contributes the median m/z and median Da error, which lets local m/z structure
drive the bias curve without forcing one global convex/concave shape.
Returns (spline, label) or (nothing, "failed").
"""
function fit_binned_regularized_bias_spline(
    x::Vector{Float64},
    y::Vector{Float64};
    n_knots::Int = TUNING_MZ_BIAS_KNOTS,
    λ::Float64 = 1.0,
    bin_size::Int = TUNING_MZ_BIAS_BIN_SIZE,
    n_irls::Int = 3,
    min_bins::Int = 5
)
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))

    valid = findall(i -> isfinite(x[i]) && isfinite(y[i]), eachindex(x))
    n_valid = length(valid)
    n_valid == 0 && return nothing, "failed (no valid points)"

    order = valid[sortperm(x[valid])]
    n_bins = max(n_valid ÷ bin_size, 1)
    n_bins < min_bins && return nothing, "failed (too few m/z bins)"

    centers = Vector{Float64}(undef, n_bins)
    vals = Vector{Float64}(undef, n_bins)
    weights = Vector{Float64}(undef, n_bins)
    for b in 1:n_bins
        s = (b - 1) * bin_size + 1
        e = b == n_bins ? n_valid : b * bin_size
        idx = order[s:e]
        centers[b] = median(x[idx])
        vals[b] = median(y[idx])
        weights[b] = Float64(length(idx))
    end

    _first = minimum(centers)
    _last = maximum(centers)
    _last <= _first && return nothing, "failed (degenerate m/z bins)"

    n_fit_knots = min(n_knots, max(3, n_bins - 3))
    basis = _make_uniform_spline_basis(centers, n_fit_knots)
    X = basis.X

    base_weights = weights ./ median(weights)
    fit_weights = copy(base_weights)
    c = zeros(basis.n_coeffs)
    resid = similar(vals)
    for irls_iter in 1:n_irls
        H, rhs = _penalized_spline_system(X, vals; λ=λ, weights=fit_weights, order=2, ridge=1e-8)
        c = H \ rhs

        irls_iter == n_irls && break
        mul!(resid, X, c)
        resid .= vals .- resid
        σ_hat = median(abs.(resid)) / 0.6744897501960817
        δ = 1.345 * max(σ_hat, 1e-10)
        for i in eachindex(fit_weights, base_weights, resid)
            r_abs = abs(resid[i])
            fit_weights[i] = base_weights[i] * (r_abs <= δ ? 1.0 : δ / r_abs)
        end
    end

    spline = _coeffs_to_spline(c, basis)
    return spline, "binned regularized (λ=$(λ), knots=$(n_fit_knots), bins=$(n_bins))"
end

#==========================================================
Monotonicity-constrained smoothing spline for intensity bias
==========================================================#

"""
    fit_monotone_bias_spline(x, y; increasing, n_knots, λ, max_iter, lr, tol)

Fit a smoothing spline to (x, y) data with a first-derivative sign constraint.

`increasing=true`:  f'(x) >= 0 — coeffs: c[i+1] >= c[i]
`increasing=false`: f'(x) <= 0 — coeffs: c[i+1] <= c[i]

Fits to ALL data points via penalized LS + projected gradient descent.
Returns a UniformSpline (Float64), or nothing if fitting fails.
"""
function fit_monotone_bias_spline(x::Vector{Float64}, y::Vector{Float64};
                                   increasing::Bool = true,
                                   n_knots::Int = 10,
                                   λ::Float64 = 1.0,
                                   max_iter::Int = 5000,
                                   lr::Float64 = 0.0001,
                                   tol::Float64 = 1e-10,
                                   n_irls::Int = 3)
    n_pts = length(x)
    n_pts < n_knots && return nothing

    basis = _make_uniform_spline_basis(x, n_knots)
    X = basis.X
    n_coeffs = basis.n_coeffs

    project!(c) = _project_monotone_coeffs!(c, increasing)
    c = _projected_irls_coeffs(X, y, n_coeffs, project!;
                               λ=λ, max_iter=max_iter, lr=lr, tol=tol, n_irls=n_irls)

    return _coeffs_to_spline(c, basis)
end

"""
    fit_best_monotone_bias(x, y; n_knots, λ, ...)

Fit both monotone-increasing and monotone-decreasing splines, return the one
with lower L1 loss. Returns (spline, label) or (nothing, "failed").
"""
function fit_best_monotone_bias(x::Vector{Float64}, y::Vector{Float64};
                                 n_knots::Int = 10, λ::Float64 = 1.0,
                                 max_iter::Int = 5000, lr::Float64 = 0.0001,
                                 tol::Float64 = 1e-10)
    spline_inc = fit_monotone_bias_spline(x, y; increasing=true, n_knots=n_knots,
                                           λ=λ, max_iter=max_iter, lr=lr, tol=tol)
    spline_dec = fit_monotone_bias_spline(x, y; increasing=false, n_knots=n_knots,
                                           λ=λ, max_iter=max_iter, lr=lr, tol=tol)

    loss_inc = spline_inc !== nothing ? sum(abs.(y .- [spline_inc(xi) for xi in x])) : Inf
    loss_dec = spline_dec !== nothing ? sum(abs.(y .- [spline_dec(xi) for xi in x])) : Inf

    if loss_inc <= loss_dec && spline_inc !== nothing
        return spline_inc, "monotone-increasing"
    elseif spline_dec !== nothing
        return spline_dec, "monotone-decreasing"
    else
        return nothing, "failed"
    end
end

# No-op projection: an unconstrained (but still λ-regularized) smoothing spline.
_noop_project!(c::Vector{Float64}) = nothing

"""
    fit_free_bias_spline(x, y; n_knots, λ, ...)

Unconstrained λ-regularized smoothing spline (same penalized-IRLS fit as
`fit_monotone_bias_spline` but with no monotonicity projection), so it can
represent a non-monotone bias such as the rise-then-dip intensity bias seen on
Sciex ZenoTOF ZT-scanning data. Returns a UniformSpline or nothing.
"""
function fit_free_bias_spline(x::Vector{Float64}, y::Vector{Float64};
                              n_knots::Int = 10, λ::Float64 = 1.0,
                              max_iter::Int = 5000, lr::Float64 = 0.0001,
                              tol::Float64 = 1e-10, n_irls::Int = 3)
    n_pts = length(x)
    n_pts < n_knots && return nothing
    basis = _make_uniform_spline_basis(x, n_knots)
    c = _projected_irls_coeffs(basis.X, y, basis.n_coeffs, _noop_project!;
                               λ=λ, max_iter=max_iter, lr=lr, tol=tol, n_irls=n_irls)
    return _coeffs_to_spline(c, basis)
end

"""
    fit_intensity_bias_spline(x, y; ...) -> (spline, label)

Dispatch the intensity-bias fit by acquisition. Sciex ZT scanning DIA
(PIONEER_ZT_METASCAN_K>0) has a non-monotone intensity bias, so use the
unconstrained `fit_free_bias_spline`; all other data keep the robust
monotone-best fit (`fit_best_monotone_bias`).
"""
function fit_intensity_bias_spline(x::Vector{Float64}, y::Vector{Float64};
                                   n_knots::Int = 10, λ::Float64 = 1.0,
                                   lr::Float64 = 0.0001)
    if something(tryparse(Int, get(ENV, "PIONEER_ZT_METASCAN_K", "")), 0) > 0
        sp = fit_free_bias_spline(x, y; n_knots=n_knots, λ=λ, lr=lr)
        return sp, "free (non-monotone, ZT)"
    else
        return fit_best_monotone_bias(x, y; n_knots=n_knots, λ=λ, lr=lr)
    end
end

#==========================================================
Regularized spline for Laplace spread
==========================================================#

"""
    fit_regularized_spread_spline(bin_centers, bin_scales; ...)

Fit an unconstrained regularized cubic B-spline to pre-binned Laplace scale.
Knot count is adaptive: at most n_knots, reduced until each knot interval has
enough binned scale estimates to avoid poorly supported local wiggles.
Takes (bin_centers, bin_scales) from `binned_laplace_scale`.
Returns a UniformSpline (Float64), or nothing.
"""
function fit_regularized_spread_spline(
    bin_centers::Vector{Float64},
    bin_scales::Vector{Float64};
    n_knots::Int = 20,
    λ::Float64 = 0.01,
    floor_val::Float64 = 0.0001
)
    n_pts = length(bin_centers)
    n_pts < 5 && return nothing

    valid = findall(i -> isfinite(bin_centers[i]) && isfinite(bin_scales[i]), eachindex(bin_centers))
    length(valid) < 5 && return nothing

    order = valid[sortperm(bin_centers[valid])]
    t = Float64.(bin_centers[order])
    y = max.(Float64.(bin_scales[order]), floor_val)
    _first = minimum(t); _last = maximum(t)
    _last <= _first && return nothing

    # Adaptive knot count: reduce n_knots until every uniform knot interval
    # contains enough binned Laplace-scale estimates.
    min_pts_per_interval = 3
    while n_knots > 1
        bin_width = (_last - _first) / (n_knots - 1)
        # Count data points in each knot interval
        min_count = typemax(Int)
        for k in 1:(n_knots - 1)
            lo = _first + (k - 1) * bin_width
            hi = lo + bin_width
            # Last interval is closed on both sides
            cnt = k < (n_knots - 1) ? count(x -> lo <= x < hi, t) :
                                       count(x -> lo <= x <= hi, t)
            min_count = min(min_count, cnt)
        end
        min_count >= min_pts_per_interval && break
        n_knots -= 1
    end
    n_knots < 3 && return nothing

    basis = _make_uniform_spline_basis(t, n_knots)
    H, rhs = _penalized_spline_system(basis.X, y; λ=λ, order=2, ridge=1e-8)
    c = H \ rhs
    for i in eachindex(c)
        c[i] = max(c[i], floor_val)
    end

    return _coeffs_to_spline(c, basis)
end

#==========================================================
Orchestrator: fit_intensity_mass_error_model
==========================================================#

"""
    fit_intensity_mass_error_model(fragments, old_model; k, conservative_k)

Fit an IntensityMassErrorModel from matched fragment data.
Bias is fit in Da space with m/z, intensity, and RT components. Spread is fit
in mDa space from binned Laplace scale estimates, with an m/z-dependent spread
correction and empirical coverage multiplier.
"""
function fit_intensity_mass_error_model(
    samples::AbstractVector{MassErrSample},
    old_model::AbstractMassErrorModel;
    k::Float32 = 3.7f0,
    conservative_k::Float32 = 6.0f0
)
    n = length(samples)
    if n < 50
        @user_warn "Too few samples ($n) to fit IntensityMassErrorModel, keeping SimpleMassErrorModel"
        return old_model
    end

    # Extract data in Da
    frag_mzs = Vector{Float32}(undef, n)
    intensities = Vector{Float32}(undef, n)
    rts = Vector{Float32}(undef, n)
    da_errs = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        s = samples[i]
        frag_mzs[i] = s.theoretical_mz
        intensities[i] = s.intensity
        rts[i] = s.rt
        da_errs[i] = Float64(s.observed_mz - s.theoretical_mz)
    end

    valid = (intensities .> 0.0f0) .& isfinite.(rts)
    if sum(valid) < 50
        @user_warn "Too few valid-intensity fragments ($(sum(valid))) for IntensityMassErrorModel"
        return old_model
    end
    frag_mzs = frag_mzs[valid]
    intensities = intensities[valid]
    rts = rts[valid]
    da_errs = da_errs[valid]
    n = length(da_errs)

    log2I = log2.(Float64.(intensities))
    mz_f64 = Float64.(frag_mzs)
    rt_f64 = Float64.(rts)

    # Step 1: m/z bias — binned robust spline, unconstrained in curvature
    mz_spline_f64, mz_label = fit_binned_regularized_bias_spline(
        mz_f64, da_errs; n_knots=TUNING_MZ_BIAS_KNOTS, λ=TUNING_MZ_BIAS_LAMBDA,
        bin_size=TUNING_MZ_BIAS_BIN_SIZE)

    if mz_spline_f64 === nothing
        fallback_spline, fallback_label = fit_best_convex_bias(mz_f64, da_errs;
            n_knots=TUNING_MZ_BIAS_KNOTS, λ=TUNING_MZ_BIAS_LAMBDA, lr=0.0001)
        mz_spline_f64 = fallback_spline
        mz_label = "fallback $(fallback_label)"
    end

    if mz_spline_f64 === nothing
        @user_warn "m/z bias spline fitting failed, keeping SimpleMassErrorModel"
        return old_model
    end

    mz_residuals = da_errs .- [mz_spline_f64(m) for m in mz_f64]

    # Step 2: intensity bias — monotone spline on ALL mz-corrected data
    I_spline_f64, I_label = fit_intensity_bias_spline(log2I, mz_residuals;
        n_knots=10, λ=1.0, lr=0.0001)

    if I_spline_f64 === nothing
        @user_warn "Intensity bias spline fitting failed, keeping SimpleMassErrorModel"
        return old_model
    end

    intensity_residuals = mz_residuals .- [I_spline_f64(li) for li in log2I]

    # Step 3: RT bias — binned robust spline on raw RT.
    rt_min = minimum(rt_f64)
    rt_max = maximum(rt_f64)
    rt_width = rt_max - rt_min
    rt_spline_f64 = _constant_bias_spline(rt_min, rt_max)
    rt_label = "none (degenerate RT range)"

    if rt_width > 0 && n >= TUNING_RT_BIAS_BIN_SIZE * 5
        fitted_rt_spline, fitted_rt_label = fit_binned_regularized_bias_spline(
            rt_f64,
            intensity_residuals;
            n_knots = TUNING_RT_BIAS_KNOTS,
            λ = TUNING_RT_BIAS_LAMBDA,
            bin_size = TUNING_RT_BIAS_BIN_SIZE,
            min_bins = 5,
        )
        if fitted_rt_spline !== nothing
            rt_center = median([fitted_rt_spline(rt) for rt in rt_f64])
            rt_spline_f64 = _shift_spline(fitted_rt_spline, rt_center)
            rt_label = fitted_rt_label
        else
            rt_label = fitted_rt_label
        end
    end

    full_residuals = intensity_residuals .- [rt_spline_f64(rt) for rt in rt_f64]

    # Step 4: Spread — binned MAD → Gaussian σ per intensity bin
    full_residuals_mda = full_residuals .* 1e3  # Da → mDa

    I_order = sortperm(log2I)
    spread_centers, spread_b_vals = binned_laplace_scale(
        I_order, log2I, full_residuals_mda, TUNING_CALIBRATION_BIN_SIZE)

    spread_spline_f64 = fit_regularized_spread_spline(spread_centers, spread_b_vals)

    if spread_spline_f64 === nothing
        @user_warn "Spread spline fitting failed, keeping SimpleMassErrorModel"
        return old_model
    end

    # Convert all splines to Float32
    mz_spline = convert_spline_to_f32(mz_spline_f64)
    I_spline = convert_spline_to_f32(I_spline_f64)
    spread_spline = convert_spline_to_f32(spread_spline_f64)

    # Step 5: m/z-dependent spread correction in mDa space
    # Regularized spline for the correction factor c(mz), with constant endpoint extrapolation.
    mz_spread_α = Float32(1.0)
    mz_spread_β = Float32(0.0)
    mz_spread_γ = Float32(0.0)
    mz_spread_spline_f64 = _constant_mz_spread_correction_spline(minimum(mz_f64), maximum(mz_f64))
    mz_spread_extrap_f64 = _constant_mz_spread_correction_extrap(mz_spread_spline_f64)
    mz_corr_label = "none (too few fragments)"
    if n >= 400
        fitted_b_mda = [max(Float64(spread_spline_f64(Float32(l2i))), 1e-4) for l2i in log2I]
        std_resid = abs.(full_residuals_mda) ./ fitted_b_mda

        mz_corr_bin_sz = TUNING_CALIBRATION_BIN_SIZE
        mz_corr_sort = sortperm(mz_f64)
        n_mz_corr_bins = max(n ÷ mz_corr_bin_sz, 1)
        corr_centers = Vector{Float64}(undef, n_mz_corr_bins)
        corr_vals = Vector{Float64}(undef, n_mz_corr_bins)

        for b in 1:n_mz_corr_bins
            s = (b - 1) * mz_corr_bin_sz + 1
            e = b == n_mz_corr_bins ? n : b * mz_corr_bin_sz
            idx = mz_corr_sort[s:e]
            corr_centers[b] = median(mz_f64[idx])
            corr_vals[b] = median(std_resid[idx]) / log(2)
        end

        mz_spread_spline_f64, mz_spread_extrap_f64,
        mz_spread_α, mz_spread_β, mz_spread_γ, mz_corr_label =
            _select_regularized_mz_spread_correction(corr_centers, corr_vals)
    end

    # Empirical k: fit k so the ±k·σ envelope covers the target fraction of
    # observed bias-corrected residuals. Per-file fit; clamped to
    # [EMPK_MIN_K, EMPK_CLAMP_HI] to avoid pathological under/over-widths.
    laplace_to_gauss = log(2.0) / 0.6744897501960817
    σ_per_frag_mda = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        σ_base = max(Float64(spread_spline_f64(Float32(log2I[i]))), 1e-4) * laplace_to_gauss
        mz_corr = Float64(_mz_spread_correction(mz_spread_spline_f64, mz_spread_extrap_f64, mz_f64[i]))
        σ_per_frag_mda[i] = σ_base * mz_corr
    end
    normalized_residuals = abs.(full_residuals_mda) ./ σ_per_frag_mda
    k_emp = Float32(quantile(normalized_residuals, EMPK_QUANTILE))
    k = clamp(k_emp, EMPK_MIN_K, EMPK_CLAMP_HI)

    # Conservative tolerance in Da = collection tolerance from the SimpleMassErrorModel.
    # SimpleMassErrorModel stores tolerance in ppm; convert to Da at max training m/z.
    max_training_mz = Float32(maximum(mz_f64))
    conservative_tol_da = if old_model isa SimpleMassErrorModel
        tol_ppm = max(getRightTol(old_model), getLeftTol(old_model))
        tol_ppm * (max_training_mz / 1f6)
    else
        # Fallback: use the model's own worst-case mDa tolerance
        max_spread_mda = max(spread_spline(spread_spline.first), Float32(1e-4))
        k * max_spread_mda / 1f3  # mDa → Da
    end

    # Default Gaussian log-density for missing top-3 fragments:
    # evaluated at 99th percentile using median intensity
    median_log2I = Float32(median(log2I))
    σ_mda_default = max(Float64(spread_spline(median_log2I)), 1e-4) * laplace_to_gauss
    σ_da_default = σ_mda_default / 1e3
    z_99 = 2.3263478740408408  # Φ⁻¹(0.99) — 99th percentile of standard normal
    default_top3_ll = Float32(-0.5 * log(2π * σ_da_default^2) - z_99^2 / 2.0)

    # Da tolerance clamp: empirical 99% quantile of bias-corrected |Da| errors.
    # More robust than max (which is driven by interference outliers).
    abs_residuals = abs.(full_residuals)
    max_da_observed = Float32(maximum(abs_residuals))
    max_da_tolerance = Float32(quantile(abs_residuals, 0.99))

    # Precompute bias extrapolation at the first/last binned training points.
    # Splines are fit on all data; outside the binned boundaries, m/z and
    # intensity bias use robust edge slopes and RT bias uses range-binned
    # edge slopes.
    # The spread spline uses constant extrapolation outside its supported range:
    # low intensity is capped at the 5% boundary, while high intensity uses the
    # fitted spline through its right edge.
    mz_qlo, mz_qhi = Float32(minimum(mz_f64)), Float32(maximum(mz_f64))
    I_qlo, I_qhi   = Float32(minimum(log2I)), Float32(maximum(log2I))
    rt_qlo, rt_qhi = Float32(minimum(rt_f64)), Float32(maximum(rt_f64))
    I_lo_spread = Float32(quantile(log2I, 0.05))
    I_hi_spread = spread_spline.last
    mz_extrap = make_edge_point_spline_extrap(
        mz_spline, mz_f64, da_errs, mz_qlo, mz_qhi)
    int_extrap = make_edge_point_spline_extrap(
        I_spline, log2I, mz_residuals, I_qlo, I_qhi)
    rt_spline = convert_spline_to_f32(rt_spline_f64)
    rt_extrap = make_rt_range_edge_spline_extrap(
        rt_spline, rt_f64, intensity_residuals, rt_qlo, rt_qhi;
        bound_bin_size = TUNING_RT_BIAS_BIN_SIZE)
    spread_extrap_raw = make_spline_extrap(spread_spline, I_lo_spread, I_hi_spread)
    spread_extrap = SplineExtrap{Float32}(
        spread_extrap_raw.lo, spread_extrap_raw.hi,
        spread_extrap_raw.lo_val, Float32(0),  # constant extrap on left  (low intensity)
        spread_extrap_raw.hi_val, Float32(0)   # constant extrap on right (high intensity)
    )
    mz_spread_spline = convert_spline_to_f32(mz_spread_spline_f64)
    mz_spread_extrap = convert_spline_extrap_to_f32(mz_spread_extrap_f64)

    model = IntensityMassErrorModel(
        mz_spline,
        I_spline,
        spread_spline,
        rt_spline,
        mz_extrap,
        int_extrap,
        spread_extrap,
        rt_extrap,
        k,
        conservative_tol_da,
        mz_spread_spline,
        mz_spread_extrap,
        mz_spread_α,
        mz_spread_β,
        mz_spread_γ,
        default_top3_ll,
        max_da_tolerance,
        Float32(rt_min),
        Float32(rt_max)
    )

    # Model's unclamped worst-case: k * spread_mDa(leftmost) * max_mz_correction
    max_spread_mda = max(Float64(spread_spline(spread_spline.first)), 1e-4)
    mz_min, mz_max = Float64(minimum(mz_f64)), Float64(maximum(mz_f64))
    _corr(mz) = Float64(_mz_spread_correction(mz_spread_spline_f64, mz_spread_extrap_f64, mz))
    max_mz_corr = max(_corr(mz_min), _corr((mz_min + mz_max) / 2), _corr(mz_max), 0.1)
    model_max_mda = Float64(k) * max_spread_mda * max_mz_corr

    mad_corr_mda = mad(full_residuals_mda, normalize=false)
    @debug_l1 "Fitted IntensityMassErrorModel (bias Da, spread mDa):" *
        "\n  m/z bias: $mz_label spline ($(length(mz_spline.coeffs)) coeffs)" *
        "\n  intensity bias: $I_label spline ($(length(I_spline.coeffs)) coeffs)" *
        "\n  RT bias: $rt_label spline ($(length(rt_spline.coeffs)) coeffs)" *
        "\n  spread: regularized spline ($(length(spread_spline.coeffs)) coeffs)" *
        "\n  m/z spread correction (mDa): $mz_corr_label ($(length(mz_spread_spline.coeffs)) coeffs; α=$(round(mz_spread_α, digits=4)), β=$(round(mz_spread_β, sigdigits=2)), γ=$(round(mz_spread_γ, sigdigits=2)))" *
        "\n  k=$(k), Da clamp (q99)=$(round(max_da_tolerance * 1e3, digits=2)) mDa, max observed=$(round(max_da_observed * 1e3, digits=2)) mDa, model max=$(round(model_max_mda, digits=2)) mDa" *
        "\n  MAD after correction: $(round(mad_corr_mda, digits=2)) mDa" *
        "\n  n_fragments=$(n)"

    return model
end

"""
    fit_scout_calibrated_model(fragments; with_intensity=true)

Fit a ScoutCalibratedMassErrorModel from wide-scout matched fragments.

Tier 1 (with_intensity=true, ≥1000 PSMs): m/z + intensity spline bias, constant mDa tol.
Tier 2 (with_intensity=false, <1000 PSMs): m/z spline bias only, constant mDa tol.

Tolerance = MAD(corrected residuals) / 0.6745 × 2.576 × 1.5 / 1e3 (Gaussian 99%, 1.5× safety).

Returns the model or `nothing` if fitting fails.
"""
function fit_scout_calibrated_model(
    samples::AbstractVector{MassErrSample};
    with_intensity::Bool = true
)
    n = length(samples)
    if n < 50
        @user_warn "Too few samples ($n) to fit ScoutCalibratedMassErrorModel"
        return nothing
    end

    frag_mzs = Vector{Float32}(undef, n)
    intensities = Vector{Float32}(undef, n)
    da_errs = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        s = samples[i]
        frag_mzs[i] = s.theoretical_mz
        intensities[i] = s.intensity
        da_errs[i] = Float64(s.observed_mz - s.theoretical_mz)
    end

    valid = intensities .> 0.0f0
    if sum(valid) < 50
        @user_warn "Too few valid-intensity samples ($(sum(valid))) for ScoutCalibratedMassErrorModel"
        return nothing
    end
    frag_mzs = frag_mzs[valid]
    intensities = intensities[valid]
    da_errs = da_errs[valid]
    n = length(da_errs)

    log2I = log2.(Float64.(intensities))
    mz_f64 = Float64.(frag_mzs)

    # Step 1: m/z bias — robust linear fit (Tukey bisquare)
    mz_label = "robust linear"
    rob_a, rob_b = 0.0, 0.0
    try
        rob_df = DataFrame(da = da_errs, mz = mz_f64)
        rob_model = rlm(@formula(da ~ mz), rob_df, TauEstimator{TukeyLoss}())
        rob_coefs = coef(rob_model)
        rob_a = Float64(rob_coefs[1])
        rob_b = Float64(rob_coefs[2])
        mz_label = "robust linear: $(round(rob_a*1e3, digits=2)) + $(round(rob_b*1e6, digits=2))e-6×mz mDa"
    catch e
        @user_warn "Robust linear m/z bias fit failed ($e), using OLS"
        mean_mz = sum(mz_f64) / n
        mean_err = sum(da_errs) / n
        sxx = sum((mz_f64[i] - mean_mz)^2 for i in 1:n)
        sxy = sum((mz_f64[i] - mean_mz) * (da_errs[i] - mean_err) for i in 1:n)
        rob_b = sxx > 0 ? sxy / sxx : 0.0
        rob_a = mean_err - rob_b * mean_mz
        mz_label = "OLS fallback"
    end

    # Create a 1-segment linear spline from the robust fit coefficients
    mz_lo, mz_hi = Float64(minimum(mz_f64)), Float64(maximum(mz_f64))
    val_lo = rob_a + rob_b * mz_lo
    val_hi = rob_a + rob_b * mz_hi
    mz_spline_f64 = UniformSpline{4, Float64}(
        SVector{4, Float64}(val_lo, (val_hi - val_lo), 0.0, 0.0),
        3, mz_lo, mz_hi, mz_hi - mz_lo)

    mz_residuals = da_errs .- (rob_a .+ rob_b .* mz_f64)

    # Step 2: intensity bias spline (tier 1 only)
    I_spline_f64 = nothing
    I_label = "none (tier 2)"
    full_residuals = mz_residuals

    if with_intensity
        I_spline_f64, I_label = fit_intensity_bias_spline(log2I, mz_residuals;
            n_knots=10, λ=1.0, lr=0.0001)
        if I_spline_f64 !== nothing
            full_residuals = mz_residuals .- [I_spline_f64(li) for li in log2I]
        else
            @user_warn "Intensity bias spline failed, falling back to m/z-only"
            with_intensity = false
        end
    end

    # If no intensity spline, create a zeroed one matching the log2I range
    if I_spline_f64 === nothing
        I_lo, I_hi = Float64(minimum(log2I)), Float64(maximum(log2I))
        I_spline_f64 = UniformSpline{4, Float64}(
            SVector{4, Float64}(0.0, 0.0, 0.0, 0.0), 3, I_lo, I_hi, I_hi - I_lo)
    end

    # Step 3: constant mDa tolerance from corrected residuals
    full_residuals_mda = full_residuals .* 1e3
    abs_resid = abs.(full_residuals_mda .- median(full_residuals_mda))
    mad_mda = median(abs_resid)
    q90_mda = quantile(abs_resid, 0.90)
    sigma_mda = mad_mda / 0.6744897501960817
    tol_mda = max(q90_mda * 3.0, 20.0)
    tolerance_da = Float32(tol_mda / 1e3)

    # Convert splines to Float32
    mz_spline = convert_spline_to_f32(mz_spline_f64)
    I_spline = convert_spline_to_f32(I_spline_f64)

    # Precompute spline-derivative bias extrapolation at binned data endpoints.
    mz_qlo, mz_qhi = Float32(minimum(mz_f64)), Float32(maximum(mz_f64))
    I_qlo, I_qhi = Float32(minimum(log2I)), Float32(maximum(log2I))
    mz_extrap = make_binned_bound_spline_extrap(
        mz_spline, mz_f64, da_errs, mz_qlo, mz_qhi)
    int_extrap = make_binned_bound_spline_extrap(
        I_spline, log2I, mz_residuals, I_qlo, I_qhi)

    model = ScoutCalibratedMassErrorModel(
        mz_spline, mz_extrap,
        I_spline, int_extrap,
        with_intensity,
        tolerance_da
    )

    @debug_l1 "Fitted ScoutCalibratedMassErrorModel:" *
        "\n  m/z bias: $mz_label ($(length(mz_spline.coeffs)) coeffs)" *
        "\n  intensity bias: $(with_intensity ? I_label : "disabled") ($(length(I_spline.coeffs)) coeffs)" *
        "\n  MAD corrected: $(round(mad_mda, digits=2)) mDa, q90: $(round(q90_mda, digits=2)) mDa" *
        "\n  tolerance: ±$(round(tol_mda, digits=2)) mDa (q90×3, floor 20 mDa)" *
        "\n  n_fragments=$n"

    return model
end
