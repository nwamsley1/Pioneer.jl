function residual2(
    t
)
end

m(t, p) = p[1] * exp.(p[2] * t[1] + t[2])
p0 = [0.5, 0.5]
tdata = collect(zip(randn(100), randn(100)))
tdata = hcat(randn(100), randn(100))
ydata = [0.3*exp(0.4*t[1] + t[2]) for t in eachrow(tdata)]
fit = curve_fit(m, tdata, ydata, p0)



function fit_linear_regression(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
    # Our loss function: sum of squared residuals
    function loss(params)
        m, b = params  # slope and intercept
        # Compute predictions
        predictions = m .* x .+ b
        # Return sum of squared residuals
        return sum((predictions .- y).^2)
    end
    
    # Initial guess for parameters [slope, intercept]
    initial_guess = [0.0f0, 0.0f0]
    
    # Optimize using LBFGS (good for smooth functions)
    result = optimize(loss, initial_guess, LBFGS())
    
    # Extract optimized parameters
    m, b = Optim.minimizer(result)
    
    return (slope=m, intercept=b)
end

x = randn(100)
y = 1 .+ (x .* 10) .+ (0.01).*randn(100)

fit_linear_regression(x, y)




function fit_ratio_function(x0::Vector{Float32}, x1::Vector{Float32}, y::Vector{Float32};
                          initial_guess::Vector{Float32}=[1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.0f0])
    
    # Helper function to compute f(x) for a single point
    function f(x, params)
        h, w, l, r, c = params
        xl = c - w/2
        xr = c + w/2
        σl = l/sqrt(2*log(2))
        σr = r/sqrt(2*log(2))
        
        if x < xl
            return h * exp(-(x - xl)^2 / (2*σl^2))
        elseif x > xr
            return h * exp(-(x - xr)^2 / (2*σr^2))
        else
            return h
        end
    end
    
    # Loss function computing log of ratios to avoid numerical issues
    function loss(params)
        h, w, l, r, c = params
        
        # Ensure parameters are valid
        if h <= 0 || h > 1 || w <= 0 || l <= 0 || r <= 0
            return Inf32
        end
        
        total_error = 0.0f0
        
        for i in 1:length(x0)
            # Compute log ratio of function values
            log_ratio = log2(f(x1[i], params)) - log2(f(x0[i], params))
            
            # Compare with target log ratio
            error = (log_ratio - y[i])^2
            total_error += error
        end
        
        return total_error
    end
    
    # Set up optimization
    lower = Float32[0.0, 0.0, 0.0, 0.0, -Inf]  # h, w, l, r, c
    upper = Float32[1.0, Inf, Inf, Inf, Inf]
    
    # Use Optimization with bounds
    result = optimize(loss, lower, upper, initial_guess, Fminbox(LBFGS()))
    
    # Extract optimized parameters
    h, w, l, r, c = Optim.minimizer(result)
    
    return (
        h = h,
        w = w,
        l = l,
        r = r,
        c = c,
        converged = Optim.converged(result),
        minimum = Optim.minimum(result),
        iterations = Optim.iterations(result)
    )
end

# Function to compute predictions using optimized parameters
function compute_predictions(x0::Vector{Float32}, x1::Vector{Float32}, params)
    function f(x, params)
        h, w, l, r, c = params.h, params.w, params.l, params.r, params.c
        xl = c - w/2
        xr = c + w/2
        σl = l/sqrt(2*log(2))
        σr = r/sqrt(2*log(2))
        
        if x < xl
            return h * exp(-(x - xl)^2 / (2*σl^2))
        elseif x > xr
            return h * exp(-(x - xr)^2 / (2*σr^2))
        else
            return h
        end
    end
    
    predictions = similar(x0)
    for i in 1:length(x0)
        predictions[i] = log2(f(x1[i], params)) - log2(f(x0[i], params))
    end
    return predictions
end
#=

#x0 = Float32[-1.0, -0.5, 0.0, 0.5, 1.0]
#x1 = Float32[0.0, 0.5, 1.0, 1.5, 2.0]
#y = Float32[0.5, 0.7, 0.9, 1.1, 1.3]
result = fit_ratio_function(x0, x1, exp2.(yt))
predictions = compute_predictions(x0, x1, result)
residuals = predictions .- y

=#