"""
    A minimal implementation of B-spline evaluation matching SciPy's behavior.
"""
struct BSpline
    t::Vector{Float64}  # knots
    c::Vector{Float64}  # coefficients 
    k::Int             # degree
end

"""
    evaluate_basis(t::Vector{Float64}, k::Int, x::Float64, idx::Int)

Evaluate a B-spline basis function value.
Based on the de Boor-Cox recursive formula.
"""
function evaluate_basis(t::Vector{Float64}, k::Int, x::Float64, idx::Int)
    if k == 0
        return Float64(t[idx] ≤ x < t[idx + 1])
    end
    
    # Handle boundary cases for the last basis function
    if idx == length(t) - k - 1 && x == t[end]
        return 1.0
    end
    
    # Standard de Boor-Cox formula
    value = 0.0
    
    # First term
    den1 = t[idx + k] - t[idx]
    if den1 != 0
        value += (x - t[idx]) / den1 * evaluate_basis(t, k - 1, x, idx)
    end
    
    # Second term
    den2 = t[idx + k + 1] - t[idx + 1]
    if den2 != 0
        value += (t[idx + k + 1] - x) / den2 * evaluate_basis(t, k - 1, x, idx + 1)
    end
    
    return value
end

"""
    evaluate(spline::BSpline, x::Float64)

Evaluate a B-spline at point x.
"""
function evaluate(spline::BSpline, x::Float64)
    # Get valid domain
    domain_start = spline.t[spline.k + 1]
    domain_end = spline.t[end - spline.k]
    
    # Clamp x to domain (this matches SciPy's behavior)
    x = clamp(x, domain_start, domain_end)
    
    # Calculate the value
    result = 0.0
    n = length(spline.t) - spline.k - 1
    
    for i in 1:n
        basis_val = evaluate_basis(spline.t, spline.k, x, i)
        result += spline.c[i] * basis_val
    end
    
    return result
end

# Convenience function to evaluate at multiple points
function evaluate(spline::BSpline, x::Vector{Float64})
    return [evaluate(spline, xi) for xi in x]
end

# Helper function to convert from SciPy format
function from_scipy_format(knots::Vector{Float64}, coeffs::Vector{Float64}, degree::Int)
    return BSpline(knots, coeffs, degree)
end

# Your SciPy data
knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coeffs = [1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
degree = 3

# Create BSpline
spline = from_scipy_format(knots, coeffs, degree)

# Evaluate at x = 37.0
result = plot(tbins, [evaluate(spline, t) for t in tbins])
