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

# Returns
Tuple containing:
- Final RT conversion model (SplineRtConversionModel, LinearRtConversionModel, or IdentityModel)
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
    min_psms::Int = 10,
    spline_degree::Int = 3,
    max_knots::Int = 7,
    outlier_threshold::Float32 = Float32(5.0)
)::Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}

    n_psms = nrow(psms)

    # Calculate adaptive knots: 1 per 100 PSMs, minimum 3
    n_knots = min(max(3, Int(floor(n_psms / 100))), max_knots)

    # Early exit for insufficient data
    if n_psms < min_psms
        @debug_l2 "Too few PSMs ($n_psms) for RT alignment (need ≥$min_psms), using identity model"
        return (IdentityModel(), Float32[], Float32[], 0.0f0)
    end

    @debug_l2 "Using $n_knots knots for RT spline ($n_psms PSMs)"

    # Single unified spline fitting path with conditional RANSAC
    try
        # Determine if we should use RANSAC based on PSM count
        use_ransac = n_psms < ransac_threshold

        @debug_l2 if use_ransac
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
        valid_mask = abs.(residuals) .< (irt_mad * outlier_threshold)
        valid_psms = psms[valid_mask, :]

        # Check if enough PSMs remain
        if nrow(valid_psms) < min_psms
            @debug_l2 "Too few PSMs after outlier removal ($(nrow(valid_psms))), using identity model"
            return (IdentityModel(), Float32[], Float32[], 0.0f0)
        end

        # Refit with filtered PSMs (same method as initial fit)
        n_knots_final = min(n_knots, max(3, Int(floor(nrow(valid_psms) / 100))))

        final_map = if use_ransac
            fit_spline_with_ransac(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                spline_degree,
                n_knots_final,
                lambda_penalty,
                2,
                ransac_iterations=50,
                ransac_sample_size=30,
                ransac_threshold_factor=Float32(2.0)
            )
        else
            UniformSplinePenalized(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                spline_degree,
                n_knots_final,
                lambda_penalty,
                2
            )
        end

        final_model = SplineRtConversionModel(final_map)

        return (final_model, valid_psms[!, :rt], valid_psms[!, :irt_predicted], irt_mad)

    catch e
        # Only fallback: use linear model on error
        @user_warn "RT spline fitting failed ($e), falling back to linear model"

        linear_model, linear_std, _ = fit_linear_irt_model(psms)

        # Calculate MAD for linear model
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
end
