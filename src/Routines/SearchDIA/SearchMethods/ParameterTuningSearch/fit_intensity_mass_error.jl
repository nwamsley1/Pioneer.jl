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
  1. m/z-dependent bias: convexity-constrained smoothing spline (Da)
  2. intensity-dependent bias: convexity-constrained smoothing spline (Da)
  3. Laplace spread: monotone-decreasing convex spline of log2(I) (**ppm**)
  4. Coverage multiplier k (from PEP-based tolerance scanning)

Bias is fit in Da space. Spread is fit in ppm space — after bias correction,
Da residuals are divided by (mz / 1e6) to get ppm, then binned Laplace scale
is computed on those ppm residuals. This yields m/z-proportional tolerance.
==========================================================#

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

        prev_obj = Inf
        for iter in 1:max_iter
            grad = H * c - Xty
            c .-= lr .* grad
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

        prev_obj = Inf
        for iter in 1:max_iter
            grad = H * c - Xty
            c .-= lr .* grad
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
    prev_obj = Inf
    for iter in 1:max_iter
        grad = H * c - Xty
        c .-= lr .* grad
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

    # Step 1: m/z bias — convex spline on ALL data
    mz_spline_f64, mz_label = fit_best_convex_bias(mz_f64, da_errs;
        n_knots=10, λ=1.0, lr=0.0001)

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
    spread_centers, spread_b_vals = binned_laplace_scale(I_order, log2I, full_residuals_mda, 100)
    # binned_laplace_scale returns Laplace b = MAD / ln(2). Convert to Gaussian σ = MAD / 0.6745
    spread_σ_vals = spread_b_vals .* (log(2.0) / 0.6744897501960817)

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
    # AICc selects among constant (1), linear (α + β·mz), quadratic (α + β·mz + γ·mz²)
    mz_spread_α = Float32(1.0)
    mz_spread_β = Float32(0.0)
    mz_spread_γ = Float32(0.0)
    mz_corr_label = "none (too few fragments)"
    if n >= 400
        fitted_b_mda = [max(Float64(spread_spline_f64(Float32(l2i))), 1e-4) for l2i in log2I]
        std_resid = abs.(full_residuals_mda) ./ fitted_b_mda

        mz_corr_bin_sz = 200
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

        nb = n_mz_corr_bins
        function _aicc(rss::Float64, n_obs::Int, k_params::Int)
            n_obs <= k_params + 1 && return Inf
            aic = n_obs * log(rss / n_obs) + 2.0 * k_params
            return k_params == 0 ? aic : aic + 2.0 * k_params * (k_params + 1) / (n_obs - k_params - 1)
        end

        # Null: c(mz) = 1.0 (constant, 0 free parameters)
        rss_null = sum((corr_vals .- 1.0).^2)
        aicc_null = _aicc(rss_null, nb, 0)

        # Linear: c(mz) = α + β·mz (2 parameters)
        X_lin = hcat(ones(nb), corr_centers)
        coefs_lin = X_lin \ corr_vals
        rss_lin = sum((corr_vals .- X_lin * coefs_lin).^2)
        aicc_lin = _aicc(rss_lin, nb, 2)

        # Quadratic: c(mz) = α + β·mz + γ·mz² (3 parameters)
        X_quad = hcat(ones(nb), corr_centers, corr_centers .^ 2)
        coefs_quad = nb > 3 ? X_quad \ corr_vals : [1.0, 0.0, 0.0]
        rss_quad = nb > 3 ? sum((corr_vals .- X_quad * coefs_quad).^2) : Inf
        aicc_quad = nb > 3 ? _aicc(rss_quad, nb, 3) : Inf

        # Select best model by AICc (require ΔAICc > 4 over null)
        best_aicc = min(aicc_null, aicc_lin, aicc_quad)
        if aicc_quad == best_aicc && aicc_null - aicc_quad > 4.0
            mz_spread_α = Float32(coefs_quad[1])
            mz_spread_β = Float32(coefs_quad[2])
            mz_spread_γ = Float32(coefs_quad[3])
            mz_corr_label = "quadratic (ΔAICc=$(round(aicc_null - aicc_quad, digits=1)))"
        elseif aicc_lin == best_aicc && aicc_null - aicc_lin > 4.0
            mz_spread_α = Float32(coefs_lin[1])
            mz_spread_β = Float32(coefs_lin[2])
            mz_corr_label = "linear (ΔAICc=$(round(aicc_null - aicc_lin, digits=1)))"
        else
            mz_corr_label = "none (best ΔAICc=$(round(aicc_null - best_aicc, digits=1)))"
        end
    end

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
    laplace_to_gauss = log(2.0) / 0.6744897501960817
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

    # Precompute linear extrapolation at [2%, 98%] quantiles of training data.
    # Splines are fit on all data but evaluation switches to linear outside
    # these boundaries to avoid trusting the spline in sparse tail regions.
    mz_q02, mz_q98 = Float32.(quantile(mz_f64, [0.02, 0.98]))
    I_q02, I_q98 = Float32.(quantile(log2I, [0.02, 0.98]))
    mz_extrap = make_spline_extrap(mz_spline, mz_q02, mz_q98)
    int_extrap = make_spline_extrap(I_spline, I_q02, I_q98)
    spread_extrap_raw = make_spline_extrap(spread_spline, I_q02, I_q98)
    spread_extrap = SplineExtrap{Float32}(
        spread_extrap_raw.lo, spread_extrap_raw.hi,
        spread_extrap_raw.lo_val, spread_extrap_raw.lo_slope,
        spread_extrap_raw.hi_val, Float32(0)  # constant extrapolation on right (high intensity)
    )

    model = IntensityMassErrorModel(
        mz_spline,
        I_spline,
        spread_spline,
        mz_extrap,
        int_extrap,
        spread_extrap,
        k,
        conservative_tol_da,
        mz_spread_α,
        mz_spread_β,
        mz_spread_γ,
        default_top3_ll,
        max_da_tolerance
    )

    # Model's unclamped worst-case: k * spread_mDa(leftmost) * max_mz_correction
    max_spread_mda = max(Float64(spread_spline(spread_spline.first)), 1e-4)
    mz_min, mz_max = Float64(minimum(mz_f64)), Float64(maximum(mz_f64))
    _corr(mz) = Float64(mz_spread_α) + Float64(mz_spread_β) * mz + Float64(mz_spread_γ) * mz^2
    max_mz_corr = max(_corr(mz_min), _corr((mz_min + mz_max) / 2), _corr(mz_max), 0.1)
    model_max_mda = Float64(k) * max_spread_mda * max_mz_corr

    mad_corr_mda = mad(full_residuals_mda, normalize=false)
    @debug_l1 "Fitted IntensityMassErrorModel (bias Da, spread mDa):" *
        "\n  m/z bias: $mz_label spline ($(length(mz_spline.coeffs)) coeffs)" *
        "\n  intensity bias: $I_label spline ($(length(I_spline.coeffs)) coeffs)" *
        "\n  spread: monotone-convex spline ($(length(spread_spline.coeffs)) coeffs)" *
        "\n  m/z spread correction (mDa): $mz_corr_label (α=$(round(mz_spread_α, digits=4)), β=$(round(mz_spread_β, sigdigits=2)), γ=$(round(mz_spread_γ, sigdigits=2)))" *
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

    # Precompute linear extrapolation at [2%, 98%] quantiles
    mz_q02, mz_q98 = Float32.(quantile(mz_f64, [0.02, 0.98]))
    I_q02, I_q98 = Float32.(quantile(log2I, [0.02, 0.98]))
    mz_extrap = make_spline_extrap(mz_spline, mz_q02, mz_q98)
    int_extrap = make_spline_extrap(I_spline, I_q02, I_q98)

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
