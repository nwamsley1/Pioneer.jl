"""
    find_interval(t::Vector{Float64}, k::Int, x::Float64, prev_l::Int, extrapolate::Bool) -> Int

Port of FITPACK's _find_interval function. Finds an interval such that t[interval] <= x < t[interval+1].

Parameters:
- t: knot vector
- k: spline degree
- x: point to find interval for
- prev_l: interval where previous value was located (use k if unknown)
- extrapolate: whether to extrapolate beyond boundaries

Returns:
- interval index or -1 if x is out of bounds or NaN
"""
function find_interval(t::Vector{Float64}, k::Int, x::Float64, prev_l::Int, extrapolate::Bool)
    n = length(t) - k - 1
    tb = t[k+1]
    te = t[n+1]
    
    # Handle nan
    if isnan(x)
        return -1
    end
    
    # Check bounds
    if ((x < tb) || (x > te)) && !extrapolate
        return -1
    end
    
    # Start at either prev_l or k, whichever is valid
    l = (k < prev_l < n) ? prev_l : k
    
    # Search for interval
    while (x < t[l+1]) && (l != k)
        l -= 1
    end
    
    l += 1
    while (x >= t[l+1]) && (l != n)
        l += 1
    end
    
    return l-1  # Convert to 0-based index for compatibility
end

"""
    deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int) -> Vector{Float64}

Port of FITPACK's _deBoor_D function. Evaluates the (k+1) non-zero B-splines of degree k at x.

Parameters:
- t: knot vector
- x: point to evaluate at
- k: spline degree
- ℓ: interval where t[ℓ] <= x < t[ℓ+1]
- m: derivative order (0 for function value)

Returns:
- Vector of k+1 non-zero basis function values
"""
function deboor_d(t::Vector{Float64}, x::Float64, k::Int, ℓ::Int, m::Int)
    # Allocate single work array like FITPACK
    work = zeros(2k + 2)
    
    # Use views into work array for h and hh
    h = view(work, 1:k+1)
    hh = view(work, k+2:2k+2)
    
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
    
    # Do m "derivative" recursions
    for j in (k-m+1):k
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

function evaluate_spline(t::Vector{Float64}, c::Matrix{Float64}, k::Int, 
                        x::Vector{Float64}, nu::Int, extrapolate::Bool)
    out = zeros(length(x), size(c, 2))
    interval = k
    
    for ip in 1:length(x)
        xval = x[ip]
        # Find correct interval
        interval = find_interval(t, k, xval, interval, extrapolate)
        
        if interval < 0
            out[ip, :] .= NaN
            continue
        end
        
        # Get (k+1) non-zero B-splines
        basis = deboor_d(t, xval, k, interval, nu)
        
        # Form linear combinations exactly as in scipy/_bspl.pyx
        for jp in 1:size(c, 2)
            out[ip, jp] = 0.0
            for a in 0:k
                # Note: original is c[interval + a - k, jp]
                idx = interval + a - k + 1  # +1 for Julia indexing
                out[ip, jp] += c[idx, jp] * basis[a+1]
            end
        end
    end
    
    return out
end

struct BSpline
    t::Vector{Float64}
    c::Matrix{Float64}
    k::Int
    extrapolate::Union{Bool,String}
    
    function BSpline(t::Vector{Float64}, c::Vector{Float64}, k::Int; 
                    extrapolate::Union{Bool,String}=true)
        # ... validation code ...
        c_mat = reshape(c, :, 1)
        new(t, c_mat, k, extrapolate)
    end
end

function (spl::BSpline)(x::Union{Real,AbstractVector}, nu::Int=0)
    x_arr = x isa AbstractVector ? collect(Float64, x) : [Float64(x)]
    
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

# Example usage:
function example(x::Float64)
    # Define a cubic spline (k=3)
    k = 3
    t = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
    
    # Find interval containing x=37.0
    #x = 37.0
    interval = find_interval(t, k, x, k, true)
    println("Found interval: ", interval)
    
    # Compute basis functions
    basis = deboor_d(t, x, k, interval, 0)
    println("Basis functions: ", basis)
    
    return interval, basis
end

function test_eval(x::Float64, knots::Vector{Float64}, coeffs::Vector{Float64}, k::Int64)
    interval = find_interval(knots, k, x, k, true)
    # Compute basis functions
    basis = deboor_d(knots, x, k, interval, 0)
    return sum(coeffs.*basis)
end

function test_bspline()
    # Same test data
    knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
    coeffs = [1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
    k = 3
    x = 37.0

    # Create spline
    spl = BSpline(knots, coeffs, k)
    
    # Test each step
    println("\nTesting each step:")
    println("1. Find interval:")
    interval = find_interval(knots, k, x, k, true)
    println("   interval = ", interval)
    
    println("\n2. Basis functions:")
    basis = deboor_d(knots, x, k, interval, 0)
    println("   basis = ", basis)
    
    println("\n3. Coefficients used:")
    for a in 0:k
        idx = interval + a - k + 1
        println("   c[", idx, "] = ", coeffs[idx], " * basis[", a+1, "] = ", basis[a+1])
    end
    
    println("\n4. Final evaluation:")
    result = spl(x)
    println("   spl(", x, ") = ", result)
    
    return result
end

# Run the test
result = test_bspline()
