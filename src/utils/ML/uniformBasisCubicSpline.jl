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

function UniformSpline(
                        u::Vector{T}, 
                        t::Vector{T}, 
                        degree::I, #Degree of the piecewise polynomials
                        n_knots::I, #Number of control points
                        ) where {I<:Integer, T<:AbstractFloat}
    if degree != 3
        error("Non-cubic splines not yet implemented. Use a degree of 3")
    end
    if n_knots < 3
        error("need at least 3 knots")
    end
    if length(u) != length(t)
        error("length(u) is not equal to length(t)")
    end

    #Uniform B Spline basis for the given degree
    #only implemented for d=3 but coule expand in the future 
    function getSplineBasis(degree::I)
        return NTuple{4, Polynomial}([
            Polynomial([0, 0, 0, 1])/6, #b0
            Polynomial([1, 3, 3, -3])/6, #b1
            Polynomial([4, 0, -6, 3])/6, #b2
            Polynomial([1, -3, 3, -1])/6, #b3
        ])
    end

    function buildDesignMat(t::Vector{T}, #location of data points
                            knots::Vector{T},
                            bin_width::T,
                            spline_basis::NTuple{4, Polynomial}
                    ) where {T<:AbstractFloat}

        function fillDesignMatRow!(X::Matrix{T}, 
                                    row::Int,
                                    knot_idx::Int,
                                    u::T,
                                    spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
            i = length(spline_basis)
            #println(" t - knot_val: ", t - knot_val)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i](u)
                i -= 1
            end
        end

        X = zeros(T, (length(t), length(knots) + 3))
        for (i, t) in enumerate(t)
            #knot_idx = min(Int64((t - first(knots))÷bin_width)+1, length(knots))
            knot_idx = min(
                            floor(Int32, (t - first(knots))/bin_width)+one(Int32), 
                            length(knots)
                            )
            fillDesignMatRow!(
                X,
                i,
                knot_idx,
                (t-knots[knot_idx])/bin_width,
                spline_basis
            )
        end

        return X
    end

    function buildPieceWise(
                    knots::Vector{T},
                    bin_width::T,
                    spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}

        function fillDesignMatRow!(X::Matrix{Polynomial}, 
                                    row::Int,
                                    knot_idx::Int,
                                    spline_basis::NTuple{4, Polynomial})
            i = length(spline_basis)
            #println(" t - knot_val: ", t - knot_val)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i]
                i -= 1
            end
        end

        X = zeros(Polynomial, (length(knots), length(knots) + length(spline_basis) - 1))
        for (i, t) in enumerate(knots)
            t = t + bin_width/2
            knot_idx = min(
                floor(Int32, (t - first(knots))/bin_width)+one(Int32), 
                length(knots)
                )
            fillDesignMatRow!(
                X,
                i,
                knot_idx,
                spline_basis
            )
        end

        return X
    end

    # insure the input is sorted so the spline fitting is stable under
    # different numbers of threads
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first)/(n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, collect(knots), bin_width, spline_basis)
    c = X\sorted_u
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly*c
    n_coeffs = n_knots*(degree + 1)
    coeffs = SVector{n_coeffs}(vcat([polynomial.coeffs for polynomial in piecewise_polynomials]...))


    UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end

function (s::UniformSpline)(t::U) where {U<:AbstractFloat}
    t = max(min(t, s.last), s.first)
    idx = floor(Int32,
                    (t - s.first)/s.bin_width
                )
    u = (t - (s.first + s.bin_width*(idx)))/s.bin_width
    x = zero(U)
    coeff = idx*(s.degree + 1) + 1
    c = one(U)
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    return x
end

"""
    build_difference_matrix(n_coeffs::Int, order::Int=2)

Build k-th order difference matrix for penalized B-splines (P-splines).

The difference matrix D is used to penalize roughness in spline coefficients:
- Order 1: Penalizes differences (cᵢ₊₁ - cᵢ)
- Order 2: Penalizes curvature (cᵢ₊₂ - 2cᵢ₊₁ + cᵢ) [RECOMMENDED]
- Order 3: Penalizes jerk (cᵢ₊₃ - 3cᵢ₊₂ + 3cᵢ₊₁ - cᵢ)

# Arguments
- `n_coeffs::Int`: Number of spline coefficients
- `order::Int`: Order of differences (1, 2, or 3; default 2)

# Returns
- `D::Matrix`: Difference matrix of size (n_coeffs - order) × n_coeffs

# Example
```julia
# For 8 coefficients with 2nd order penalty
D = build_difference_matrix(8, 2)  # Returns 6×8 matrix
P = D' * D  # 8×8 penalty matrix
```

# References
Eilers & Marx (1996). "Flexible smoothing with B-splines and penalties."
Statistical Science, 11(2), 89-121.
"""
function build_difference_matrix(n_coeffs::Int, order::Int=2)
    if order < 1 || order > 3
        error("Difference order must be 1, 2, or 3")
    end

    if order == 1
        # First order difference: D[i,:] = [0, ..., -1, 1, ..., 0]
        n_rows = n_coeffs - 1
        D = zeros(Float64, n_rows, n_coeffs)
        for i in 1:n_rows
            D[i, i] = -1.0
            D[i, i+1] = 1.0
        end
        return D

    elseif order == 2
        # Second order: recursively compute D² = D₁(n-1) × D₁(n)
        D1_n = build_difference_matrix(n_coeffs, 1)  # (n-1) × n
        D1_n_minus_1 = build_difference_matrix(n_coeffs - 1, 1)  # (n-2) × (n-1)
        return D1_n_minus_1 * D1_n  # (n-2) × n

    else  # order == 3
        # Third order: D³ = D₁(n-2) × D²(n)
        D2 = build_difference_matrix(n_coeffs, 2)  # (n-2) × n
        D1_n_minus_2 = build_difference_matrix(n_coeffs - 2, 1)  # (n-3) × (n-2)
        return D1_n_minus_2 * D2  # (n-3) × n
    end
end

"""
    UniformSplinePenalized(u, t, degree, n_knots, λ=0.1, order=2)

Fit uniform B-spline with difference penalty regularization (P-splines).

This uses the P-spline approach (Eilers & Marx, 1996) which adds a penalty term
to the least squares objective:

    minimize: ||Xc - u||² + λ||D^k c||²

The penalty encourages smooth coefficient sequences, which may improve monotonicity
by reducing oscillations caused by noise or outliers.

# Arguments
- `u::Vector{T}`: Target values (iRT in RT alignment context)
- `t::Vector{T}`: Input values (RT in RT alignment context)
- `degree::I`: Polynomial degree (must be 3 for cubic splines)
- `n_knots::I`: Number of knots (typically 5)
- `λ::T`: Penalty parameter (default 0.1)
  - λ = 0: No penalty (standard least squares)
  - λ small (0.01-0.1): Mild smoothing
  - λ medium (0.1-1.0): Moderate smoothing [RECOMMENDED]
  - λ large (1.0-10.0): Strong smoothing
- `order::Int`: Order of difference penalty (default 2)
  - Order 2 is recommended (penalizes curvature changes)

# Returns
- `UniformSpline{N, T}`: Fitted spline with difference penalty

# Notes
- Does NOT guarantee monotonicity, but improves it significantly
- For 5 knots: 8 coefficients, penalty matrix is 6×8 (2nd order)

# Example
```julia
# Basic usage with defaults
spline = UniformSplinePenalized(irt, rt, 3, 5)

# Custom penalty strength
spline = UniformSplinePenalized(irt, rt, 3, 5, 0.5, 2)
```

# References
Eilers, P. H. C., & Marx, B. D. (1996). Statistical Science, 11(2), 89-121.
"""
function UniformSplinePenalized(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
    λ::T = T(0.1),
    order::Int = 2
) where {I<:Integer, T<:AbstractFloat}

    # Input validation
    if degree != 3
        error("Non-cubic splines not yet implemented. Use degree = 3")
    end
    if n_knots < 3
        error("Need at least 3 knots")
    end
    if length(u) != length(t)
        error("length(u) must equal length(t)")
    end
    if λ < 0
        error("Penalty parameter λ must be non-negative")
    end

    # Sort data for numerical stability
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    # Build B-spline basis (reuse existing internal functions)
    function getSplineBasis(degree::I)
        return NTuple{4, Polynomial}([
            Polynomial([0, 0, 0, 1])/6,
            Polynomial([1, 3, 3, -3])/6,
            Polynomial([4, 0, -6, 3])/6,
            Polynomial([1, -3, 3, -1])/6,
        ])
    end

    function buildDesignMat(t::Vector{T}, knots::Vector{T}, bin_width::T, spline_basis)
        function fillDesignMatRow!(X, row, knot_idx, u, spline_basis)
            i = length(spline_basis)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i](u)
                i -= 1
            end
        end

        X = zeros(T, (length(t), length(knots) + 3))
        for (i, t_val) in enumerate(t)
            knot_idx = min(
                floor(Int32, (t_val - first(knots))/bin_width) + one(Int32),
                length(knots)
            )
            fillDesignMatRow!(
                X, i, knot_idx,
                (t_val - knots[knot_idx])/bin_width,
                spline_basis
            )
        end
        return X
    end

    function buildPieceWise(knots, bin_width, spline_basis)
        function fillDesignMatRow!(X, row, knot_idx, spline_basis)
            i = length(spline_basis)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i]
                i -= 1
            end
        end

        X = zeros(Polynomial, (length(knots), length(knots) + length(spline_basis) - 1))
        for (i, t_val) in enumerate(knots)
            t_adj = t_val + bin_width/2
            knot_idx = min(
                floor(Int32, (t_adj - first(knots))/bin_width) + one(Int32),
                length(knots)
            )
            fillDesignMatRow!(X, i, knot_idx, spline_basis)
        end
        return X
    end

    # Setup spline parameters
    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, knots, bin_width, spline_basis)

    # Compute coefficients with penalty
    n_coeffs = n_knots + degree  # = 8 for n_knots=5, degree=3

    if λ == 0
        # No penalty: standard least squares
        c = X \ sorted_u
    else
        # Penalized least squares
        D = build_difference_matrix(n_coeffs, order)
        P = D' * D  # Penalty matrix (n_coeffs × n_coeffs)

        # Solve: (X'X + λP)c = X'u
        c = (X'X + T(λ) * P) \ (X'sorted_u)
    end

    # Build piecewise polynomials
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly * c
    n_total_coeffs = n_knots * (degree + 1)
    coeffs = SVector{n_total_coeffs}(vcat([poly.coeffs for poly in piecewise_polynomials]...))

    return UniformSpline{n_total_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end
