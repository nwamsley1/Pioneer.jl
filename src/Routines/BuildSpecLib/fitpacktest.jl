"""
    fpbspl(t::Vector{Float64}, n::Int, k::Int, x::Float64, l::Int)

Evaluates the (k+1) non-zero B-splines of degree k at x using the efficient stable recurrence 
relation of de Boor and Cox. 
This is an exact port of the FITPACK fpbspl subroutine.
"""
function fpbspl(t::Vector{Float64}, n::Int, k::Int, x::Float64, l::Int)
    h = zeros(Float64, 6)  # Assumes max degree of 5
    hh = zeros(Float64, 5)
    
    h[1] = 1.0
    
    for j in 1:k
        for i in 1:j
            hh[i] = h[i]
        end
        h[1] = 0.0
        
        for i in 1:j
            li = l + i
            lj = li - j
            f = hh[i] / (t[li] - t[lj])
            h[i] += f * (t[li] - x)
            h[i+1] = f * (x - t[lj])
        end
    end
    
    return h
end

"""
    splev(t::Vector{Float64}, c::Vector{Float64}, k::Int, x::Vector{Float64}, ext::Int=0)

Evaluate a B-spline at specified points.
This is an exact port of the FITPACK splev subroutine.

Parameters:
- t: knot vector
- c: B-spline coefficients
- k: degree of the spline
- x: points where to evaluate the spline
- ext: extrapolation mode (0:extrapolate, 1:return zeros, 2:raise error, 3:return boundary)

Returns:
- y: evaluated values
- ier: error code (0:success, 1:error)
"""
function splev(t::Vector{Float64}, c::Vector{Float64}, k::Int, x::Vector{Float64}, ext::Int=0)
    n = length(t)
    m = length(x)
    y = zeros(Float64, m)
    
    # Input validation
    if m < 1
        return y, 10  # Error
    end
    
    # Get boundaries of approximation interval
    k1 = k + 1
    nk1 = n - k1
    tb = t[k1]
    te = t[nk1 + 1]
    
    l = k1
    l1 = l + 1
    
    # Main loop for the different points
    for i in 1:m
        # Fetch new x-value arg
        arg = x[i]
        
        # Check if arg is in the support
        if arg < tb || arg > te
            if ext == 0
                # Do nothing, evaluate anyway by extrapolation
            elseif ext == 1
                y[i] = 0.0
                continue
            elseif ext == 2
                return y, 1  # Error
            elseif ext == 3
                # Set arg to nearest boundary
                arg = arg < tb ? tb : te
            end
        end
        
        # Find knot interval t[l] <= arg < t[l+1]
        while arg ≥ t[l+1] && l != nk1
            l = l1
            l1 = l + 1
        end
        
        # Evaluate the non-zero B-splines at arg
        h = fpbspl(t, n, k, arg, l)
        
        # Compute value at x = arg
        sp = 0.0
        ll = l - k1
        for j in 1:k1
            ll += 1
            sp += c[ll] * h[j]
        end
        y[i] = sp
    end
    
    return y, 0  # Success
end

"""
    evaluate_bspline(knots::Vector{Float64}, coeffs::Vector{Float64}, degree::Int, 
                    x::Union{Float64, Vector{Float64}}, ext::Int=0)

High-level function to evaluate B-spline with SciPy-compatible interface.

Parameters:
- knots: knot vector from SciPy
- coeffs: coefficients from SciPy
- degree: spline degree from SciPy
- x: evaluation point(s)
- ext: extrapolation mode (same as SciPy)

Returns:
- Evaluated spline value(s)
"""
function evaluate_bspline(knots::Vector{Float64}, coeffs::Vector{Float64}, degree::Int, 
                         x::Union{Float64, Vector{Float64}}, ext::Int=0)
    x_vec = x isa Vector ? x : [x]
    y, ier = splev(knots, coeffs, degree, x_vec, ext)
    if ier == 1
        throw(ArgumentError("x value out of bounds and ext=2"))
    elseif ier == 10
        throw(ArgumentError("Invalid input data"))
    end
    return x isa Vector ? y : y[1]
end