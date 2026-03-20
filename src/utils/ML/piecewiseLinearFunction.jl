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

function PiecewiseNceModel(x::T) where {T<:AbstractFloat}
    PiecewiseNceModel(zero(T), one(T), zero(T), x, zero(T))
end

"""
   fit_nce_model(pwlm::PiecewiseNceModel{T}, x::AbstractVector, y::AbstractVector, 
                charge::AbstractVector, breakpoint::Real) where {T<:AbstractFloat}

Fit a piecewise model with charge dependence to normalized collision energy data, ensuring continuity at the breakpoint.

# Arguments
- `pwlm::PiecewiseNceModel{T}`: Template for the model type to fit
- `x::AbstractVector`: Mass-to-charge ratio values
- `y::AbstractVector`: Observed collision energy values
- `charge::AbstractVector`: Charge states
- `breakpoint::Real`: Point where model transitions from linear to constant behavior

# Returns
- `PiecewiseNceModel{T}`: Fitted model with parameters:
   - left_slope: Slope of linear region (x ≤ breakpoint)
   - left_intercept: Intercept of linear region
   - right_value: Constant value for x > breakpoint (= left_slope * breakpoint + left_intercept)
   - charge_slope: Linear charge dependence coefficient

# Model
For mass-to-charge ratio x and charge state z:
- When x ≤ breakpoint: f(x,z) = left_slope * x + left_intercept + charge_slope * z
- When x > breakpoint: f(x,z) = right_value + charge_slope * z

The model is fit by minimizing the sum of squared residuals while maintaining continuity 
at the breakpoint through the constraint: right_value = left_slope * breakpoint + left_intercept.

# Notes
- Initial parameter estimates are obtained using linear regression on the left region
- Optimization is performed using the LBFGS algorithm
- Charge dependence is linear and consistent across both regions
- Continuity at breakpoint is enforced by construction

See also: [`PiecewiseNceModel`](@ref)
"""
function fit_nce_model(
    pwlm::PiecewiseNceModel{T},
    x::AbstractVector,
    y::AbstractVector,
    charge::AbstractVector,
    breakpoint::Real
) where {T<:AbstractFloat}
    # Define the objective function to minimize sum of squared residuals
    function objective(params)
        # params[1] = left_slope
        # params[2] = left_intercept
        # params[3] = charge_slope
        
        # Calculate the right_value from continuity constraint
        right_value = params[1] * breakpoint + params[2]
        
        # Calculate residuals for left side
        left_mask = x .<= breakpoint
        y_pred_left = params[1] .* x[left_mask] .+
                     params[3] .* charge[left_mask] .+
                     params[2]
        residuals_left = y_pred_left .- y[left_mask]
        
        # Calculate residuals for right side
        right_mask = x .> breakpoint
        y_pred_right = fill(right_value, sum(right_mask)) .+
                     params[3] .* charge[right_mask]
        residuals_right = y_pred_right .- y[right_mask]
        
        # Total sum of squared residuals
        return sum(residuals_left.^2) + sum(residuals_right.^2)
    end
    
    # Initial guess using simple linear regression
    left_mask = x .<= breakpoint
    X_left = [ones(sum(left_mask)) x[left_mask] charge[left_mask]]
    initial_coef = X_left \ y[left_mask]
    
    initial_guess = [
        initial_coef[2],  # left_slope
        initial_coef[1],  # left_intercept
        initial_coef[3]   # charge_slope
    ]
    
    # Optimize
    result = optimize(objective, initial_guess, LBFGS())
    optimal_params = Optim.minimizer(result)
    
    # Calculate right_value using continuity constraint
    right_value = optimal_params[1] * breakpoint + optimal_params[2]
    
    return PiecewiseNceModel(
        T(float(breakpoint)),
        T(optimal_params[1]),  # left_slope
        T(optimal_params[2]),  # left_intercept
        T(right_value),        # right_value from continuity
        T(optimal_params[3])   # charge_slope
    )
end
 
function (f::PiecewiseNceModel)()
    return f.right_value
end

function (f::PiecewiseNceModel)(x::AbstractFloat, charge::Integer)
base = if x <= f.breakpoint
    f.left_slope * x + f.left_intercept
else
    f.right_value
end
return base + f.charge_slope * charge
end
 
 # Add method for vectors
 function (f::PiecewiseNceModel)(x::AbstractVector, charge::AbstractVector)
    return map((x,c) -> f(x,c), x, charge)
 end