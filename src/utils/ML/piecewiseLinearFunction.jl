struct PiecewiseLinear
    breakpoint::Float64
    left_slope::Float64
    left_intercept::Float64
    right_slope::Float64
    right_intercept::Float64
end

function fit_piecewise_linear(x::AbstractVector, y::AbstractVector, breakpoint::Real)
    # Split data into left and right of breakpoint
    left_mask = x .<= breakpoint
    right_mask = x .> breakpoint
    
    # Fit left segment
    X_left = [ones(sum(left_mask)) x[left_mask]]
    left_coef = X_left \ y[left_mask]
    left_intercept, left_slope = left_coef
    
    # Fit right segment
    #X_right = [ones(sum(right_mask)) x[right_mask]]
    #=
    X_right = [ones(sum(right_mask)) x[right_mask]]
    right_coef = X_right \ y[right_mask]
    right_intercept, right_slope = right_coef
    =#
    right_intercept = mean(y[right_mask])
    right_slope = zero(eltype(y))
    # Create and return callable struct
    return PiecewiseLinear(
        float(breakpoint),
        left_slope,
        left_intercept,
        right_slope,
        right_intercept
    )
end

# Make the struct callable
function (f::PiecewiseLinear)(x::AbstractFloat)
    if x <= f.breakpoint
        return f.left_slope * x + f.left_intercept
    else
        return f.right_slope * x + f.right_intercept
    end
end

# Add method for vectors
function (f::PiecewiseLinear)(x::AbstractVector)
    return map(f, x)
end