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
  3. Laplace spread: monotone-decreasing convex spline of log2(I) (**ppm**)
  4. Coverage multiplier k (from PEP-based tolerance scanning)

Bias is fit in Da space. Spread is fit in ppm space — after bias correction,
Da residuals are divided by (mz / 1e6) to get ppm, then binned Laplace scale
is computed on those ppm residuals. This yields m/z-proportional tolerance.
==========================================================#

# Empirical-k coverage parameters. Quantile of |residual|/σ used to set the
# coverage multiplier k; the result is clamped to [1.96, EMPK_CLAMP_HI].
const EMPK_QUANTILE = 0.95
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

    return SplineExtrap{T}(
        lo, hi,
        spline(lo), T(lo_slope === nothing ? 0.0 : lo_slope),
        spline(hi), T(hi_slope === nothing ? 0.0 : hi_slope)
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

function _fit_regularized_mz_spread_correction_spline(
    corr_centers::Vector{Float64},
    corr_vals::Vector{Float64};
    n_knots::Int = 8,
    λ_candidates::Tuple = (100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01),
    min_local_ratio::Float64 = 0.85,
    floor_val::Float64 = 0.1,
    ceiling_val::Float64 = 10.0,
)
    n_pts = length(corr_centers)
    n_pts < 5 && return nothing, "constant (too few m/z correction bins)"

    order = sortperm(corr_centers)
    x = corr_centers[order]
    y = clamp.(corr_vals[order], floor_val, ceiling_val)
    _first = minimum(x)
    _last = maximum(x)
    _last <= _first && return nothing, "constant (degenerate m/z correction bins)"

    n_knots = min(n_knots, max(3, n_pts - 3))
    knots = collect(LinRange(_first, _last, n_knots))
    bin_width = (_last - _first) / (n_knots - 1)
    X = _build_numeric_design_matrix(x, knots, bin_width)
    n_coeffs = n_knots + 3
    D2 = build_difference_matrix(n_coeffs, 2)
    P = D2' * D2
    ridge = 1e-6 * Matrix{Float64}(I, n_coeffs, n_coeffs)
    Xty = X' * y
    XtX = X' * X

    best_spline = nothing
    best_label = "regularized spline"
    best_score = Inf
    best_λ = last(λ_candidates)

    for λ in λ_candidates
        H = XtX + λ * P + ridge
        c = H \ Xty
        spline = _coeffs_to_spline(c, knots, bin_width, 3, n_knots, _first, _last)
        pred = clamp.([Float64(spline(xi)) for xi in x], floor_val, ceiling_val)
        rss = sum((pred .- y).^2)
        roughness = sum((D2 * c).^2)
        score = rss + 0.01 * roughness
        if all(pred .>= min_local_ratio .* y)
            return spline, "regularized spline (λ=$(λ), bins=$(n_pts))"
        end
        if score < best_score
            best_score = score
            best_spline = spline
            best_λ = λ
        end
    end

    return best_spline, "regularized spline (λ=$(best_λ), undercoverage guard relaxed, bins=$(n_pts))"
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

    spline, label = _fit_regularized_mz_spread_correction_spline(centers, vals)
    if spline === nothing
        value = median(clamp.(vals, 0.1, 10.0))
        spline = _constant_mz_spread_correction_spline(minimum(centers), maximum(centers); value=value)
    end

    return spline, _constant_mz_spread_correction_extrap(spline),
           Float32(1.0), Float32(0.0), Float32(0.0), label
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

    _first = minimum(x); _last = maximum(x)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = _build_numeric_design_matrix(x, knots, bin_width)
    n_coeffs = n_knots + 3

    D2 = build_difference_matrix(n_coeffs, 2)
    P = D2' * D2

    function project!(c::Vector{Float64})
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
    end

    w = ones(n_pts)
    c = zeros(n_coeffs)

    for irls_iter in 1:n_irls
        WX = Diagonal(w) * X
        H = WX' * X + λ * P
        Xty = WX' * y
        c = H \ Xty
        project!(c)

        # Lipschitz-safe step cap: opnorm(H) is the gradient Lipschitz
        # constant of this convex quadratic, so any lr ≤ 2/L is stable.
        # Without this cap, fixed lr diverges once n_pts is large enough
        # that L grows past 1/lr (LU factorization at the next IRLS
        # iteration then throws on the resulting Inf/NaN coefficients).
        lr_eff = min(lr, 1.8 / opnorm(H))

        prev_obj = Inf
        for iter in 1:max_iter
            grad = H * c - Xty
            c .-= lr_eff .* grad
            project!(c)
            if iter % 200 == 0
                obj = 0.5 * (c' * H * c) - (Xty' * c)
                if abs(prev_obj - obj) / max(abs(obj), 1.0) < tol
                    break
                end
                prev_obj = obj
            end
        end

        irls_iter == n_irls && break
        resid = y .- X * c
        σ_hat = median(abs.(resid)) / 0.6744897501960817
        δ = 1.345 * max(σ_hat, 1e-10)
        w = [abs(r) <= δ ? 1.0 : δ / abs(r) for r in resid]
    end

    degree = 3
    return _coeffs_to_spline(c, knots, bin_width, degree, n_knots, _first, _last)
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
    knots = collect(LinRange(_first, _last, n_fit_knots))
    bin_width = (_last - _first) / (n_fit_knots - 1)
    X = _build_numeric_design_matrix(centers, knots, bin_width)
    n_coeffs = n_fit_knots + 3
    D2 = build_difference_matrix(n_coeffs, 2)
    P = D2' * D2
    ridge = 1e-8 * Matrix{Float64}(I, n_coeffs, n_coeffs)

    base_weights = weights ./ median(weights)
    fit_weights = copy(base_weights)
    c = zeros(n_coeffs)
    for irls_iter in 1:n_irls
        WX = Diagonal(fit_weights) * X
        H = WX' * X + λ * P + ridge
        Xty = X' * (fit_weights .* vals)
        c = H \ Xty

        irls_iter == n_irls && break
        resid = vals .- X * c
        σ_hat = median(abs.(resid)) / 0.6744897501960817
        δ = 1.345 * max(σ_hat, 1e-10)
        huber_weights = [abs(r) <= δ ? 1.0 : δ / abs(r) for r in resid]
        fit_weights = base_weights .* huber_weights
    end

    spline = _coeffs_to_spline(c, knots, bin_width, 3, n_fit_knots, _first, _last)
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

    _first = minimum(x); _last = maximum(x)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = _build_numeric_design_matrix(x, knots, bin_width)
    n_coeffs = n_knots + 3

    D2 = build_difference_matrix(n_coeffs, 2)
    P = D2' * D2

    function project!(c::Vector{Float64})
        for _ in 1:30
            changed = false
            for i in 1:(length(c)-1)
                if increasing ? (c[i+1] < c[i]) : (c[i+1] > c[i])
                    avg = (c[i] + c[i+1]) / 2
                    c[i] = avg; c[i+1] = avg; changed = true
                end
            end
            !changed && break
        end
    end

    w = ones(n_pts)
    c = zeros(n_coeffs)

    for irls_iter in 1:n_irls
        WX = Diagonal(w) * X
        H = WX' * X + λ * P
        Xty = WX' * y
        c = H \ Xty
        project!(c)

        # See fit_convex_bias_spline for rationale.
        lr_eff = min(lr, 1.8 / opnorm(H))

        prev_obj = Inf
        for iter in 1:max_iter
            grad = H * c - Xty
            c .-= lr_eff .* grad
            project!(c)
            if iter % 200 == 0
                obj = 0.5 * (c' * H * c) - (Xty' * c)
                if abs(prev_obj - obj) / max(abs(obj), 1.0) < tol
                    break
                end
                prev_obj = obj
            end
        end

        irls_iter == n_irls && break
        resid = y .- X * c
        σ_hat = median(abs.(resid)) / 0.6744897501960817
        δ = 1.345 * max(σ_hat, 1e-10)
        w = [abs(r) <= δ ? 1.0 : δ / abs(r) for r in resid]
    end

    degree = 3
    return _coeffs_to_spline(c, knots, bin_width, degree, n_knots, _first, _last)
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

#==========================================================
Monotone-decreasing convex spline for Laplace spread
==========================================================#

"""
    fit_monotone_convex_spread_spline(bin_centers, bin_scales; ...)

Fit a monotone-decreasing, convex cubic B-spline to pre-binned Laplace scale.
Knot count is adaptive: at most n_knots, reduced to maintain a 10:1
data-to-coefficient ratio.
Takes (bin_centers, bin_scales) from `binned_laplace_scale`.
Returns a UniformSpline (Float64), or nothing.
"""
function fit_monotone_convex_spread_spline(
    bin_centers::Vector{Float64},
    bin_scales::Vector{Float64};
    n_knots::Int = 20,
    λ::Float64 = 0.01,
    max_iter::Int = 5000,
    lr::Float64 = 0.001,
    tol::Float64 = 1e-8,
    floor_val::Float64 = 0.0001
)
    n_pts = length(bin_centers)
    n_pts < 5 && return nothing

    t = Float64.(bin_centers)
    y = Float64.(bin_scales)
    _first = minimum(t); _last = maximum(t)

    # Adaptive knot count: reduce n_knots until every uniform knot interval
    # contains at least `min_pts_per_interval` data points (Laplace bin centers).
    min_pts_per_interval = 5
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
    n_knots < 1 && return nothing

    spline_bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = _build_numeric_design_matrix(t, knots, spline_bin_width)
    n_coeffs = n_knots + 3

    D2 = build_difference_matrix(n_coeffs, 2)
    P = D2' * D2
    H = X' * X + λ * P
    Xty = X' * y
    c = H \ Xty

    function project!(c::Vector{Float64})
        for _ in 1:30
            changed = false
            # Convexity: f''(x) >= 0
            for i in 1:(length(c)-2)
                v = c[i+2] - 2*c[i+1] + c[i]
                if v < 0
                    c[i] += (-v) / 3; c[i+1] -= 2 * (-v) / 3; c[i+2] += (-v) / 3
                    changed = true
                end
            end
            !changed && break
        end
        for i in eachindex(c); c[i] = max(c[i], floor_val); end
    end

    project!(c)
    # See fit_convex_bias_spline for rationale.
    lr_eff = min(lr, 1.8 / opnorm(H))
    prev_obj = Inf
    for iter in 1:max_iter
        grad = H * c - Xty
        c .-= lr_eff .* grad
        project!(c)
        if iter % 100 == 0
            obj = 0.5 * (c' * H * c) - (Xty' * c)
            if abs(prev_obj - obj) / max(abs(obj), 1.0) < tol; break; end
            prev_obj = obj
        end
    end

    degree = 3
    return _coeffs_to_spline(c, knots, spline_bin_width, degree, n_knots, _first, _last)
end

#==========================================================
Orchestrator: fit_intensity_mass_error_model
==========================================================#

"""
    fit_intensity_mass_error_model(fragments, old_model; k, conservative_k)

Fit an IntensityMassErrorModel from matched fragment data.
Bias in Da space; spread in ppm space. Bias uses convexity-constrained
smoothing splines fit to all per-fragment data; spread uses monotone-decreasing
convex spline fit to binned Laplace scale values of ppm residuals.
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
    da_errs = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        s = samples[i]
        frag_mzs[i] = s.theoretical_mz
        intensities[i] = s.intensity
        da_errs[i] = Float64(s.observed_mz - s.theoretical_mz)
    end

    valid = intensities .> 0.0f0
    if sum(valid) < 50
        @user_warn "Too few valid-intensity fragments ($(sum(valid))) for IntensityMassErrorModel"
        return old_model
    end
    frag_mzs = frag_mzs[valid]
    intensities = intensities[valid]
    da_errs = da_errs[valid]
    n = length(da_errs)

    log2I = log2.(Float64.(intensities))
    mz_f64 = Float64.(frag_mzs)

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
    I_spline_f64, I_label = fit_best_monotone_bias(log2I, mz_residuals;
        n_knots=10, λ=1.0, lr=0.0001)

    if I_spline_f64 === nothing
        @user_warn "Intensity bias spline fitting failed, keeping SimpleMassErrorModel"
        return old_model
    end

    full_residuals = mz_residuals .- [I_spline_f64(li) for li in log2I]

    # Step 3: Spread — binned MAD → Gaussian σ per intensity bin
    full_residuals_mda = full_residuals .* 1e3  # Da → mDa

    I_order = sortperm(log2I)
    spread_centers, spread_b_vals = binned_laplace_scale(
        I_order, log2I, full_residuals_mda, TUNING_CALIBRATION_BIN_SIZE)

    spread_spline_f64 = fit_monotone_convex_spread_spline(spread_centers, spread_b_vals)

    if spread_spline_f64 === nothing
        @user_warn "Spread spline fitting failed, keeping SimpleMassErrorModel"
        return old_model
    end

    # Convert all splines to Float32
    mz_spline = convert_spline_to_f32(mz_spline_f64)
    I_spline = convert_spline_to_f32(I_spline_f64)
    spread_spline = convert_spline_to_f32(spread_spline_f64)

    # Step 4: m/z-dependent spread correction in mDa space
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

    # Empirical k: fit k so the ±k·σ envelope covers ~95% of observed
    # bias-corrected residuals. Replaces the previous static k=1.96 (which
    # only achieves 95% under a Gaussian assumption — Olsen Exploris fragments
    # showed ~86% empirical coverage at k=1.96, evidence the residual tails
    # are heavier than Gaussian). Per-file fit; clamped to [1.96, 4.5] to
    # avoid pathological blow-up on noisy files.
    # 2026-05-15: 98% target was tested and dramatically regressed (−50k prec,
    # −3.7k PG on 23-file Olsen) because BitVec LUT calibration trained on
    # the wider window learns a stricter score cutoff. 95% is the sweet spot.
    laplace_to_gauss = log(2.0) / 0.6744897501960817
    σ_per_frag_mda = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        σ_base = max(Float64(spread_spline_f64(Float32(log2I[i]))), 1e-4) * laplace_to_gauss
        mz_corr = Float64(_mz_spread_correction(mz_spread_spline_f64, mz_spread_extrap_f64, mz_f64[i]))
        σ_per_frag_mda[i] = σ_base * mz_corr
    end
    normalized_residuals = abs.(full_residuals_mda) ./ σ_per_frag_mda
    k_emp = Float32(quantile(normalized_residuals, EMPK_QUANTILE))
    k = clamp(k_emp, 1.96f0, EMPK_CLAMP_HI)

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

    # Precompute extrapolation boundaries at central quantiles of training data.
    # Splines are fit on all data but evaluation switches to extrapolation
    # outside these boundaries.
    # The SPREAD spline uses CONSTANT extrapolation on both ends. The
    # low-intensity tail can be sparsely sampled and the monotone-convex
    # spline's linear extrapolation tended to blow up the spread there,
    # which derails BitVecCalibration (admits more low-intensity noise →
    # fewer bits clear the target/decoy excess threshold → narrower
    # candidate set → fewer PSMs). Constant extrapolation beyond [q05, q95]
    # caps the spread at its value at the boundary.
    bias_q_lo, bias_q_hi = TUNING_BIAS_EXTRAP_QUANTILES
    mz_qlo, mz_qhi = Float32.(quantile(mz_f64, [bias_q_lo, bias_q_hi]))
    I_qlo, I_qhi   = Float32.(quantile(log2I, [bias_q_lo, bias_q_hi]))
    I_q05, I_q95   = Float32.(quantile(log2I, [0.05, 0.95]))
    mz_extrap = make_edge_point_spline_extrap(
        mz_spline, mz_f64, da_errs, mz_qlo, mz_qhi)
    int_extrap = make_edge_point_spline_extrap(
        I_spline, log2I, mz_residuals, I_qlo, I_qhi)
    spread_extrap_raw = make_spline_extrap(spread_spline, I_q05, I_q95)
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
        mz_extrap,
        int_extrap,
        spread_extrap,
        k,
        conservative_tol_da,
        mz_spread_spline,
        mz_spread_extrap,
        mz_spread_α,
        mz_spread_β,
        mz_spread_γ,
        default_top3_ll,
        max_da_tolerance
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
        "\n  spread: monotone-convex spline ($(length(spread_spline.coeffs)) coeffs)" *
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
        I_spline_f64, I_label = fit_best_monotone_bias(log2I, mz_residuals;
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

    # Precompute linear extrapolation at central bias quantiles. Scout keeps
    # endpoint-derivative slopes so the Wide Scout bias diagnostics stay
    # comparable to the pre-edge-slope model.
    bias_q_lo, bias_q_hi = TUNING_BIAS_EXTRAP_QUANTILES
    mz_qlo, mz_qhi = Float32.(quantile(mz_f64, [bias_q_lo, bias_q_hi]))
    I_qlo, I_qhi = Float32.(quantile(log2I, [bias_q_lo, bias_q_hi]))
    mz_extrap = make_spline_extrap(mz_spline, mz_qlo, mz_qhi)
    int_extrap = make_spline_extrap(I_spline, I_qlo, I_qhi)

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
