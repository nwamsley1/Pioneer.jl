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

#=
function (f::NceModelContainer)(c::UInt8, x::AbstractFloat)
    #Suppose the charge-state bounds are 2-4
    #Then the indices for the 2, 3, and 4 charge state models
    #are 1, 2, and 3 respectively. To evaluate at a charge outside 
    #the bound, default to the nearest. So 1 -> 2. 6 -> 4 in this example. 
    if c < first(f.charge_bounds)
        first(f.nce_models)(x)
    elseif c > last(f.charge_bounds)
        last(f.nce_models)(x)
    else
        f.nce_models[c-first(f.charge_bounds)+1](x)
    end
end
=#

 
function PiecewiseNceModel(x::T) where {T<:AbstractFloat}
    PiecewiseNceModel(zero(T), one(T), zero(T), x, zero(T))
end

 """
   PiecewiseNceModel{T<:AbstractFloat} <: NceModel{T}

A piecewise model for normalized collision energy prediction that includes charge dependence.

# Type Parameters
- `T`: Floating point precision type

# Fields
- `breakpoint::T`: x-value where model transitions from linear to constant
- `left_slope::T`: Slope of linear component for x ≤ breakpoint
- `left_intercept::T`: Intercept of linear component for x ≤ breakpoint
- `right_value::T`: Constant value for x > breakpoint
- `charge_slope::T`: Linear charge dependence coefficient

# Model
For a given mass-to-charge ratio x and charge state z:
- When x ≤ breakpoint: f(x,z) = left_slope * x + left_intercept + charge_slope * z
- When x > breakpoint: f(x,z) = right_value + charge_slope * z

# Methods
   (f::PiecewiseNceModel)(x::AbstractFloat, charge::Integer)
   (f::PiecewiseNceModel)(x::AbstractVector, charge::AbstractVector)
   fit_nce_model(pwlm::PiecewiseNceModel, x, y, charge, breakpoint)

See also: [`fit_nce_model`](@ref)
"""
#=
 function fit_nce_model(
    pwlm::PiecewiseNceModel{T},
    x::AbstractVector,
    y::AbstractVector,
    charge::AbstractVector,
    breakpoint::Real) where {T<:AbstractFloat}
    # Define the objective function to minimize sum of squared residuals
    function objective(params)
        # params[1] = left_slope
        # params[2] = left_intercept
        # params[3] = right_value
        # params[4] = charge_slope
        
        # Calculate residuals for left side
        left_mask = x .<= breakpoint
        y_pred_left = params[1] .* x[left_mask] .+
                     params[4] .* charge[left_mask] .+
                     params[2]
        residuals_left = y_pred_left .- y[left_mask]
        
        # Calculate residuals for right side
        right_mask = x .> breakpoint
        y_pred_right = fill(params[3], sum(right_mask)) .+
                     params[4] .* charge[right_mask]
        residuals_right = y_pred_right .- y[right_mask]
        
        # Total sum of squared residuals
        return sum(residuals_left.^2) + sum(residuals_right.^2)
    end
    
    # Initial guess using simple linear regression
    left_mask = x .<= breakpoint
    X_left = [ones(sum(left_mask)) x[left_mask] charge[left_mask]]
    initial_coef = X_left \ y[left_mask]
    
    # Initial guess for right value
    right_mask = x .> breakpoint
    initial_right = mean(y[right_mask])
    
    initial_guess = [
        initial_coef[2],  # left_slope
        initial_coef[1],  # left_intercept
        initial_right,    # right_value
        initial_coef[3]   # charge_slope
    ]
    
    # Optimize
    result = optimize(objective, initial_guess, LBFGS())
    optimal_params = Optim.minimizer(result)
    
    return PiecewiseNceModel(
        T(float(breakpoint)),
        T(optimal_params[1]), # left_slope
        T(optimal_params[2]), # left_intercept
        T(optimal_params[3]), # right_value (now fitted directly)
        T(optimal_params[4])  # charge_slope
    )
end
=#
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



 #=


function PiecewiseConstantRhsNceModel(x::T) where {T<:AbstractFloat}
    return PiecewiseConstantRhsNceModel(
        NCE_MODEL_BREAKPOINT #Global constant
        , zero(T), x, x
    )
end

"""
    fit_nce_model(pwlm::PiecewiseConstantRhsNceModel, x::AbstractVector, y::AbstractVector, breakpoint::Real)

Fit a piecewise linear-constant model to data with fixed breakpoint.

# Arguments
- `pwlm`: Model template
- `x`: Input values
- `y`: Target values
- `breakpoint`: Fixed x-value for transition point

# Returns
New `PiecewiseConstantRhsNceModel` with optimized parameters minimizing squared error.

# Implementation
Uses LBFGS optimization to find left slope and intercept that minimize squared 
residuals. Right-hand value is determined by continuity at breakpoint.
"""
function fit_nce_model(pwlm::PiecewiseConstantRhsNceModel{T}, x::AbstractVector, y::AbstractVector, breakpoint::Real) where {T<:AbstractFloat}
    # Define the objective function to minimize sum of squared residuals
    function objective(params)
        # params[1] = left_slope
        # params[2] = left_intercept
        
        # Calculate residuals for left side
        left_mask = x .<= breakpoint
        y_pred_left = params[1] .* x[left_mask] .+ params[2]
        residuals_left = y_pred_left .- y[left_mask]
        
        # Right value is determined by continuity at breakpoint
        right_value = params[1] * breakpoint + params[2]
        
        # Calculate residuals for right side
        right_mask = x .> breakpoint
        y_pred_right = fill(right_value, sum(right_mask))
        residuals_right = y_pred_right .- y[right_mask]
        
        # Total sum of squared residuals
        return sum(residuals_left.^2) + sum(residuals_right.^2)
    end
    
    # Initial guess using simple linear regression on left side
    left_mask = x .<= breakpoint
    X_left = [ones(sum(left_mask)) x[left_mask]]
    initial_coef = X_left \ y[left_mask]
    initial_guess = [initial_coef[2], initial_coef[1]]  # [slope, intercept]
    
    # Optimize
    result = optimize(objective, initial_guess, LBFGS())
    optimal_params = Optim.minimizer(result)
    
    # Calculate right value
    right_value = optimal_params[1] * breakpoint + optimal_params[2]
    
    return PiecewiseConstantRhsNceModel(
        T(float(breakpoint)),
        T(optimal_params[1]),  # left_slope
        T(optimal_params[2]),  # left_intercept
        T(right_value)
    )
end

# Make the struct callable
function (f::PiecewiseConstantRhsNceModel)(x::AbstractFloat)
    if x <= f.breakpoint
        return f.left_slope * x + f.left_intercept
    else
        return f.right_value
    end
end

function (f::PiecewiseConstantRhsNceModel)()
    return f.left_intercept
end


# Add method for vectors
function (f::PiecewiseConstantRhsNceModel)(x::AbstractVector)
    return map(f, x)
end

function fit_nce_model(
    pwlm::PiecewiseConstantRhsNceModel, 
    x::AbstractVector, 
    y::AbstractVector,
    charge::AbstractVector,
    breakpoint::Real
)
    # Define the objective function to minimize sum of squared residuals
    function objective(params)
        # params[1] = left_slope
        # params[2] = left_intercept
        # params[3] = charge_slope
        
        # Calculate residuals for left side
        left_mask = x .<= breakpoint
        y_pred_left = params[1] .* x[left_mask] .+ 
                     params[3] .* charge[left_mask] .+ 
                     params[2]
        residuals_left = y_pred_left .- y[left_mask]
        
        # Right value is determined by continuity at breakpoint (excluding charge effect)
        right_value = params[1] * breakpoint + params[2]
        
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
    initial_guess = [initial_coef[2], initial_coef[1], initial_coef[3]] # [slope, intercept, charge_slope]
    
    # Optimize
    result = optimize(objective, initial_guess, LBFGS())
    optimal_params = Optim.minimizer(result)
    
    # Calculate right value (base value without charge effect)
    right_value = optimal_params[1] * breakpoint + optimal_params[2]
    
    # Need to modify the PiecewiseConstantRhsNceModel struct to include charge_slope
    return PiecewiseConstantRhsNceModel(
        float(breakpoint),
        optimal_params[1],  # left_slope
        optimal_params[2],  # left_intercept
        right_value,
        optimal_params[3]   # charge_slope
    )
end


 =#