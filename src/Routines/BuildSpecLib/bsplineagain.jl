#=
"""
Direct port of FITPACK's B-spline evaluation routines from __fitpack.cc
"""

"""
    find_interval(t::Vector{Float64}, k::Int, x::Float64, prev_l::Int, extrapolate::Bool)

Find an interval such that t[interval] <= xval < t[interval+1].
Direct port of FITPACK's _find_interval function.
"""
function find_interval(t::Vector{Float64}, k::Int, x::Float64, prev_l::Int, extrapolate::Bool)
    n = length(t) - k - 1
    tb = t[k+1]
    te = t[n+1]
    
    # Handle nan
    if isnan(x)
        return -1
    end
    
    # Handle out of bounds
    if ((x < tb) || (x > te)) && !extrapolate
        return -1
    end
    
    # Initial interval guess
    l = (k < prev_l < n) ? prev_l : k
    
    # Search for interval where t[l] <= x < t[l+1]
    while (x < t[l+1]) && (l != k)
        l -= 1
    end
    
    l += 1
    while (x >= t[l+1]) && (l != n)
        l += 1
    end
    
    return l
end

"""
    deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int)

Core B-spline evaluation routine.
Direct port of FITPACK's _deBoor_D function.

On completion returns the k+1 non-zero values of beta^(m)_i,k(x) 
for i=ℓ, ℓ-1, ℓ-2, ℓ-k where t[ℓ] <= x < t[ℓ+1].
"""
function deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int)
    # Allocate work arrays
    h = zeros(k + 1)
    hh = zeros(k + 1)
    
    # Initial value
    h[1] = 1.0
    
    # Perform k-m "standard" deBoor iterations
    for j in 1:k-m
        # Save current values
        copyto!(hh, 1, h, 1, j)
        h[1] = 0.0
        
        for n in 1:j
            ind = ℓ + n
            xb = t[ind+1]
            xa = t[ind-j+1]
            
            if xb == xa
                h[n+1] = 0.0
                continue
            end
            
            w = hh[n] / (xb - xa)
            h[n] += w * (xb - x)
            h[n+1] = w * (x - xa)
        end
    end
    
    # Now do m "derivative" recursions
    for j in (k-m+1):k
        # Save current values
        copyto!(hh, 1, h, 1, j)
        h[1] = 0.0
        
        for n in 1:j
            ind = ℓ + n
            xb = t[ind+1]
            xa = t[ind-j+1]
            
            if xb == xa
                h[m+1] = 0.0
                continue
            end
            
            w = j * hh[n] / (xb - xa)
            h[n] -= w
            h[n+1] = w
        end
    end
    
    return h[1:k+1]
end

"""
    evaluate_spline(t::Vector{Float64}, c::Matrix{Float64}, k::Int, x::Vector{Float64}, 
                   nu::Int, extrapolate::Bool)

Evaluate a B-spline or its derivatives.
Direct port of FITPACK's spline evaluation logic.
"""
function evaluate_spline(t::Vector{Float64}, c::Matrix{Float64}, k::Int, x::Vector{Float64}, 
                        nu::Int, extrapolate::Bool)
    out = zeros(length(x), size(c, 2))
    interval = k
    
    for i in 1:length(x)
        xval = x[i]
        
        # Find correct interval
        interval = find_interval(t, k, xval, interval, extrapolate)
        
        if interval < 0
            out[i, :] .= NaN
            continue
        end
        
        # Get k+1 non-zero basis functions
        basis = deboor_d(t, xval, k, interval, nu)
        
        # Form linear combinations
        for j in 1:size(c, 2)
            out[i, j] = 0.0
            for a in 1:k+1
                out[i, j] += c[interval - k + a, j] * basis[a]
            end
        end
    end
    
    return out
end

"""
Scipy-compatible BSpline class implementation
"""
struct BSpline
    t::Vector{Float64}  # knots
    c::Matrix{Float64}  # coefficients 
    k::Int             # degree
    extrapolate::Union{Bool,String}  # extrapolation mode
    
    function BSpline(t::Vector{Float64}, c::Vector{Float64}, k::Int; 
                    extrapolate::Union{Bool,String}=true)
        # Input validation
        if k < 0
            throw(ArgumentError("Spline order cannot be negative."))
        end
        
        if length(t) < 2k + 2
            throw(ArgumentError("Need at least $(2k + 2) knots for degree $k"))
        end
        
        if length(c) < length(t) - k - 1
            throw(ArgumentError("Too few coefficients"))
        end
        
        if any(diff(t) .< 0)
            throw(ArgumentError("Knots must be in non-decreasing order"))
        end
        
        # Reshape coefficients to matrix form if needed
        c_mat = reshape(c, :, 1)
        
        new(t, c_mat, k, extrapolate)
    end
end

# Evaluation function
function (spl::BSpline)(x::Union{Real,AbstractVector}, nu::Int=0)
    x_arr = x isa AbstractVector ? collect(Float64, x) : [Float64(x)]
    
    # Handle periodic extrapolation 
    if spl.extrapolate == "periodic"
        n = length(spl.t) - spl.k - 1
        period = spl.t[n] - spl.t[spl.k+1]
        x_arr = @. spl.t[spl.k+1] + mod(x_arr - spl.t[spl.k+1], period)
        extrapolate = false
    else
        extrapolate = spl.extrapolate
    end
    
    result = evaluate_spline(spl.t, spl.c, spl.k, x_arr, nu, extrapolate)
    
    return x isa AbstractVector ? vec(result) : result[1]
end
=#
"""
    evaluate_spline(t::Vector{Float64}, c::Matrix{Float64}, k::Int, x::Vector{Float64}, 
                   nu::Int, extrapolate::Bool)

Evaluate a B-spline or its derivatives.
Direct port of FITPACK's spline evaluation logic.
"""
function evaluate_spline(t::Vector{Float64}, c::Matrix{Float64}, k::Int, x::Vector{Float64}, 
                        nu::Int, extrapolate::Bool)
    out = zeros(length(x), size(c, 2))
    interval = k
    
    for i in 1:length(x)
        xval = x[i]
        
        # Find interval where t[interval] <= x < t[interval+1]
        interval = find_interval(t, k, xval, interval, extrapolate)
        
        if interval < 0
            out[i, :] .= NaN
            continue
        end
        
        # Get k+1 non-zero basis functions 
        basis = deboor_d(t, xval, k, interval, nu)
        
        # For interval ℓ, the basis functions are ordered:
        # β_{ℓ-k,k}(x), β_{ℓ-k+1,k}(x), ..., β_{ℓ,k}(x)
        for j in 1:size(c, 2)
            val = 0.0
            for a in 0:k
                # Start from ℓ-k and go up to ℓ
                idx = interval - k + a
                println("Using c[", idx, "] with basis[", a+1, "]")
                val += c[idx, j] * basis[a+1]
            end
            out[i,j] = val
        end
    end
    
    return out
end
# Rest of the code remains the same
"""
Direct port of FITPACK's _find_interval function.
"""

function find_interval(t::Vector{Float64}, k::Int, x::Float64, prev_l::Int, extrapolate::Bool)
    n = length(t) - k - 1
    tb = t[k+1]
    te = t[n+1]
    
    if isnan(x)
        return -1
    end
    
    if ((x < tb) || (x > te)) && !extrapolate
        return -1
    end
    
    # Start at previous interval if valid
    l = (k < prev_l < n) ? prev_l : k
    
    # First search left
    while (x < t[l+1]) && (l != k)
        l -= 1
    end
    
    # Then increment and search right
    l += 1
    while (x >= t[l+1]) && (l != n)
        l += 1
    end
    
    return l
end
"""
    deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int)

Core B-spline evaluation routine, direct port of FITPACK's _deBoor_D function.
"""

function deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int)
    # Following FITPACK exactly:
    # Returns the k+1 non-zero values of beta^(m)_i,k(x) for i=ℓ,ℓ-1,...,ℓ-k
    
    # Allocate a single work array like FITPACK 
    work = zeros(2k + 2)
    h = view(work, 1:k+1)        # first k+1 elements 
    hh = view(work, k+2:2k+2)    # last k+1 elements
    
    # Initial value
    h[1] = 1.0
    
    # Perform k-m "standard" deBoor iterations
    for j in 1:k-m
        # Save current values
        hh[1:j] .= h[1:j]
        h[1] = 0.0
        
        for n in 1:j
            ind = ℓ + n
            xb = t[ind+1]
            xa = t[ind-j+1]
            
            if xb == xa
                h[n+1] = 0.0
                continue
            end
            
            w = hh[n] / (xb - xa)
            h[n] += w * (xb - x)
            h[n+1] = w * (x - xa)
        end
    end
    
    # Now do m "derivative" recursions
    for j in (k-m+1):k
        hh[1:j] .= h[1:j]
        h[1] = 0.0
        
        for n in 1:j
            ind = ℓ + n
            xb = t[ind+1]
            xa = t[ind-j+1]
            
            if xb == xa
                h[m+1] = 0.0
                continue
            end
            
            w = j * hh[n] / (xb - xa)
            h[n] -= w
            h[n+1] = w
        end
    end
    
    return h[1:k+1]
end


struct BSpline
    t::Vector{Float64}  # knots
    c::Matrix{Float64}  # coefficients 
    k::Int             # degree
    extrapolate::Union{Bool,String}  # extrapolation mode
    
    function BSpline(t::Vector{Float64}, c::Vector{Float64}, k::Int; 
                    extrapolate::Union{Bool,String}=true)
        # Input validation
        if k < 0
            throw(ArgumentError("Spline order cannot be negative."))
        end
        
        if length(t) < 2k + 2
            throw(ArgumentError("Need at least $(2k + 2) knots for degree $k"))
        end
        
        if length(c) < length(t) - k - 1
            throw(ArgumentError("Too few coefficients"))
        end
        
        if any(diff(t) .< 0)
            throw(ArgumentError("Knots must be in non-decreasing order"))
        end
        
        # Reshape coefficients to matrix form if needed
        c_mat = reshape(c, :, 1)
        
        new(t, c_mat, k, extrapolate)
    end
end

# Evaluation function
function (spl::BSpline)(x::Union{Real,AbstractVector}, nu::Int=0)
    x_arr = x isa AbstractVector ? collect(Float64, x) : [Float64(x)]
    
    # Handle periodic extrapolation 
    if spl.extrapolate == "periodic"
        n = length(spl.t) - spl.k - 1
        period = spl.t[n] - spl.t[spl.k+1]
        x_arr = @. spl.t[spl.k+1] + mod(x_arr - spl.t[spl.k+1], period)
        extrapolate = false
    else
        extrapolate = spl.extrapolate
    end
    
    result = evaluate_spline(spl.t, spl.c, spl.k, x_arr, nu, extrapolate)
    
    return x isa AbstractVector ? vec(result) : result[1]
end

knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coeffs = [1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
degree = 3

# Create BSpline
spl = BSpline(knots, coeffs, degree)

# Evaluate at x = 37.0
result = spl(37.0)
