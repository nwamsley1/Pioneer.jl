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
Shared utilities for RT (Retention Time) to iRT (indexed Retention Time) alignment.

This module provides functions for fitting RT conversion models using various methods:
- Linear regression (robust M-estimator with Tukey bisquare loss)
- P-spline fitting with optional RANSAC for robustness
- Adaptive knot selection based on data size
- Outlier removal via MAD threshold

These utilities are shared across ParameterTuningSearch, FirstPassSearch, and potentially
other search methods to ensure consistent RT alignment methodology.
"""

"""
    fit_linear_irt_model(psms::DataFrame)

Fits a linear model mapping RT to iRT using robust M-estimation.

Uses RobustModels.jl TauEstimator with Tukey bisquare loss for robust parameter
estimation via iteratively reweighted least squares. Falls back to OLS if robust
regression fails.

# Arguments
- `psms`: DataFrame with columns `:rt` and `:irt_predicted`

# Returns
Tuple of:
- LinearRtConversionModel: Fitted linear model
- Float32: Residual standard deviation
- Int: Number of parameters (always 2 for linear)

# Examples
```julia
psms = DataFrame(rt = [...], irt_predicted = [...])
model, residual_std, n_params = fit_linear_irt_model(psms)
predicted_irt = model(observed_rt)
```
"""
function fit_linear_irt_model(
    psms::DataFrame
)::Tuple{LinearRtConversionModel, Float32, Int}

    rt = psms[!, :rt]
    irt = psms[!, :irt_predicted]
    n = length(rt)

    # Handle degenerate case - all RTs identical or insufficient variance
    if std(rt) < 1e-10
        return (
            LinearRtConversionModel(0.0f0, Float32(mean(irt))),
            Float32(std(irt)),
            2
        )
    end

    try
        # Create DataFrame for RobustModels (requires formula interface)
        fit_data = DataFrame(
            irt = Float64.(irt),  # Response variable
            rt = Float64.(rt)      # Predictor variable
        )

        # Fit robust M-estimator with Tukey bisquare loss
        # MEstimator provides robust parameter estimation via iteratively reweighted
        # least squares, with better stability for small sample sizes than MM-estimator
        robust_fit = rlm(@formula(irt ~ rt), fit_data, TauEstimator{TukeyLoss}())

        # Extract coefficients using coef() function
        coeffs = coef(robust_fit)
        intercept = coeffs[1]  # First coeff is intercept
        slope = coeffs[2]      # Second coeff is slope for 'rt'

        # Calculate residuals
        predicted = intercept .+ slope .* rt
        residuals = irt .- predicted
        residual_std = Float32(std(residuals))

        model = LinearRtConversionModel(Float32(slope), Float32(intercept))

        return (model, residual_std, 2)

    catch e
        # Fallback to OLS if robust regression fails
        @warn "Robust regression failed: $e. Falling back to OLS."

        # Simple least squares fallback
        sum_rt = sum(rt)
        sum_rt2 = sum(rt .^ 2)
        sum_irt = sum(irt)
        sum_rt_irt = sum(rt .* irt)

        det = n * sum_rt2 - sum_rt^2
        if abs(det) < 1e-10
            return (
                LinearRtConversionModel(0.0f0, Float32(mean(irt))),
                Float32(std(irt)),
                2
            )
        end

        intercept = (sum_rt2 * sum_irt - sum_rt * sum_rt_irt) / det
        slope = (n * sum_rt_irt - sum_rt * sum_irt) / det

        predicted = intercept .+ slope .* rt
        residuals = irt .- predicted
        residual_std = Float32(std(residuals))

        model = LinearRtConversionModel(Float32(slope), Float32(intercept))
        return (model, residual_std, 2)
    end
end

"""
    make_spline_monotonic(original_spline, rt_data, irt_data)

Enforces monotonic increasing property on a fitted RT→iRT spline using bidirectional
cumulative max filter from the median, then creates a linear interpolation.

# Algorithm
1. Sample original spline at uniform grid (~N PSMs points)
2. Find median RT in original data
3. Apply cumulative max filter moving right (enforce increasing)
4. Apply cumulative max filter moving left (enforce increasing in reverse)
5. Create linear interpolation of filtered points with constant extrapolation

# Arguments
- `original_spline::UniformSpline`: Fitted spline that may violate monotonicity
- `rt_data::Vector{Float32}`: Original RT values from PSM data
- `irt_data::Vector{Float32}`: Original iRT values from PSM data

# Returns
- `Interpolations.Extrapolation`: Linear interpolation guaranteed to be monotonic

# Examples
```julia
monotonic_interp = make_spline_monotonic(original_spline, rt_psms, irt_psms)
```
"""
function make_spline_monotonic(
    original_spline::UniformSpline,
    rt_data::Vector{Float32},
    irt_data::Vector{Float32}
)::Interpolations.Extrapolation

    # 1. Sample from original spline at uniform grid
    rt_min, rt_max = extrema(rt_data)
    n_sample_points = length(rt_data) - 1
    rt_grid = collect(LinRange(Float32(rt_min), Float32(rt_max), n_sample_points))
    irt_grid = [original_spline(r) for r in rt_grid]

    # 2. Find median split point
    median_rt = median(rt_data)
    median_idx = argmin(abs.(rt_grid .- median_rt))

    # 3. Filter right side (median → end): enforce monotonic increasing
    for i in (median_idx+1):n_sample_points
        if irt_grid[i] < irt_grid[i-1]
            irt_grid[i] = irt_grid[i-1]
        end
    end

    # 4. Filter left side (median → start): enforce monotonic increasing in reverse
    for i in (median_idx-1):-1:1
        if irt_grid[i] > irt_grid[i+1]
            irt_grid[i] = irt_grid[i+1]
        end
    end

    # 5. Create linear interpolation of filtered monotonic data
    # Linear interpolation exactly preserves monotonicity and requires no knot selection
    final_interp = linear_interpolation(
        rt_grid,          # x values (sorted RT)
        irt_grid,         # y values (filtered monotonic iRT)
        extrapolation_bc = Interpolations.Flat()  # Constant extrapolation outside range
    )

    # 6. Validate monotonicity (optional diagnostic)
    violations = sum(diff(irt_grid) .< 0)
    if violations > 0
        @user_info "Monotonic filter had $violations violations (edge corrections applied)"
    end

    return final_interp
end

"""
    fit_irt_model(psms::DataFrame; kwargs...)

Fits retention time alignment model between library and empirical retention times.

Uses RANSAC-based spline fitting for limited data (<ransac_threshold PSMs) and
standard spline fitting for abundant data. Linear model used only as error fallback.

# Arguments
- `psms`: DataFrame with columns `:rt` and `:irt_predicted`

# Keyword Arguments
- `lambda_penalty::Float32 = 0.1f0`: P-spline smoothing penalty
- `ransac_threshold::Int = 1000`: PSM count below which RANSAC is used
- `min_psms::Int = 10`: Minimum PSMs required for any RT model
- `spline_degree::Int = 3`: Degree of B-spline basis
- `max_knots::Int = 7`: Maximum number of knots allowed
- `outlier_threshold::Float32 = 5.0f0`: MAD multiplier for outlier removal

# Process
1. Determines fitting strategy based on PSM count and configuration
2. Performs initial spline fit (with RANSAC if limited data)
3. Removes outliers based on MAD threshold
4. Refits spline model on cleaned data
5. Applies monotonic enforcement with linear interpolation

# Returns
Tuple containing:
- Final RT conversion model (InterpolationRtConversionModel, LinearRtConversionModel, or IdentityModel)
- Valid RT values (Float32[])
- Valid iRT values (Float32[])
- iRT median absolute deviation (Float32)

# Examples
```julia
psms = DataFrame(rt = [...], irt_predicted = [...])

# Use defaults
model, rt, irt, mad = fit_irt_model(psms)

# Custom parameters
model, rt, irt, mad = fit_irt_model(psms,
    lambda_penalty = 0.05f0,
    ransac_threshold = 500,
    min_psms = 20
)
```
"""
function fit_irt_model(
    psms::DataFrame;
    lambda_penalty::Float32 = Float32(0.1),
    ransac_threshold::Int = 1000,
    min_psms::Int = 30,
    spline_degree::Int = 3,
    max_knots::Int = 7,
    outlier_threshold::Float32 = Float32(5.0)
)::Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}

    n_psms = nrow(psms)

    # Calculate adaptive knots: 1 per 100 PSMs, minimum 3
    @user_info "max_knots set to: $max_knots"
    @user_info "min_psms set to: $min_psms"
    lambda_penalty = Float32(0.1)
    @user_info "lambda_penalty set to $lambda_penalty \n"
    min_knots = max(min(10, n_psms ÷ 20), 3) # At least 3 knots, up to 1 per 20 PSMs, capped at 10
    n_knots = min(max(min_knots, Int(floor(n_psms / 100))), max_knots)
    @user_info "n_knots set to: $n_knots \n n_psms: $n_psms \n"
    # Early exit for insufficient data
    if n_psms < min_psms
        @user_info "Too few PSMs ($n_psms) for RT alignment (need ≥$min_psms), using identity model"
        return (IdentityModel(), Float32[], Float32[], 0.0f0)
    end

    @user_info "Using $n_knots knots for RT spline ($n_psms PSMs)"

    # For small datasets (< 100 PSMs), use simple linear model instead of splines
    # Linear models are more robust and appropriate for sparse data
    if n_psms < 30
        @user_info "Small dataset ($n_psms PSMs < 100): Using linear model instead of spline"
        linear_model, linear_std, _ = fit_linear_irt_model(psms)

        # Calculate MAD for consistency with spline path
        predicted = [linear_model(rt) for rt in psms[!, :rt]]
        residuals = psms[!, :irt_predicted] .- predicted
        linear_mad = median(abs.(residuals .- median(residuals)))

        return (
            linear_model,
            psms[!, :rt],
            psms[!, :irt_predicted],
            Float32(linear_mad)
        )
    end

    # Single unified spline fitting path with conditional RANSAC
    try
        # Determine if we should use RANSAC based on PSM count
        use_ransac = n_psms < ransac_threshold

        @user_info if use_ransac
            "Limited data ($n_psms PSMs): Using RANSAC + penalty (λ=$lambda_penalty)"
        else
            "Abundant data ($n_psms PSMs): Using standard fitting (λ=$lambda_penalty)"
        end

        # Fit spline (with or without RANSAC)
        rt_to_irt_map = if use_ransac
            # Use RANSAC for robustness with limited data
            fit_spline_with_ransac(
                psms[!, :irt_predicted],
                psms[!, :rt],
                spline_degree,
                n_knots,
                lambda_penalty,
                2,  # 2nd order penalty
                ransac_iterations=50,
                ransac_sample_size=30,
                ransac_threshold_factor=Float32(2.0)
            )
        else
            # Standard P-spline fitting for abundant data
            UniformSplinePenalized(
                psms[!, :irt_predicted],
                psms[!, :rt],
                spline_degree,
                n_knots,
                lambda_penalty,
                2  # 2nd order penalty
            )
        end

        # Calculate residuals
        predicted_irt = [rt_to_irt_map(rt) for rt in psms[!, :rt]]
        residuals = psms[!, :irt_predicted] .- predicted_irt

        # Remove outliers
        irt_mad = mad(residuals, normalize=true)::Float32

        outlier_limit = outlier_threshold * irt_mad
        inlier_mask = abs.(residuals) .<= outlier_limit
        valid_psms = psms[inlier_mask, :]

        rt_to_irt_map = if use_ransac
            # Use RANSAC for robustness with limited data
            fit_spline_with_ransac(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                spline_degree,
                n_knots,
                lambda_penalty,
                2,  # 2nd order penalty
                ransac_iterations=50,
                ransac_sample_size=30,
                ransac_threshold_factor=Float32(2.0)
            )
        else
            # Standard P-spline fitting for abundant data
            UniformSplinePenalized(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                spline_degree,
                n_knots,
                lambda_penalty,
                2  # 2nd order penalty
            )
        end

        # Apply monotonic enforcement to prevent backwards slopes at edges
        @user_info "Applying monotonic enforcement with bidirectional cumulative max filter"
        final_map_monotonic = make_spline_monotonic(
            rt_to_irt_map,
            valid_psms[!, :rt],
            valid_psms[!, :irt_predicted]
        )
        final_model = InterpolationRtConversionModel(final_map_monotonic)
        return (final_model, psms[!, :rt], psms[!, :irt_predicted], irt_mad)

    catch e
        # Unexpected failure during spline fitting - fall back to linear model
        @user_warn "RT spline fitting failed unexpectedly ($e), falling back to linear model \n"
        @user_warn "Full stack trace:"
        for (exc, bt) in Base.catch_stack()
            showerror(stderr, exc, bt)
            println(stderr)
        end
        linear_model, linear_std, _ = fit_linear_irt_model(psms)

        # Calculate MAD for linear model (σ-estimate)
        predicted = [linear_model(rt) for rt in psms[!, :rt]]
        residuals = psms[!, :irt_predicted] .- predicted
        linear_mad = mad(residuals, normalize=true)::Float32

        return (
            linear_model,
            psms[!, :rt],
            psms[!, :irt_predicted],
            linear_mad
        )
    end
end
