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

#############################################################################
# Shared helpers for uniform cubic B-spline construction
#############################################################################

"""
Cubic B-spline basis polynomials (b0..b3) for uniform knots.
"""
const _CUBIC_SPLINE_BASIS = NTuple{4, Polynomial}([
    Polynomial([0, 0, 0, 1])/6,   # b0
    Polynomial([1, 3, 3, -3])/6,  # b1
    Polynomial([4, 0, -6, 3])/6,  # b2
    Polynomial([1, -3, 3, -1])/6, # b3
])

"""
    _n_spline_coeffs(n_knots) -> Int

Number of cubic B-spline basis functions over `n_knots` uniform knots, and hence
the column count of the design matrix.

`n_knots` knots bound `n_knots - 1` spans, and a cubic basis over those spans has
`(n_knots - 1) + 3` functions. Callers that need to know how much data is required
to identify a fit should use this rather than open-coding the arithmetic.
"""
_n_spline_coeffs(n_knots::Integer) = Int(n_knots) + 2

"""
    _find_knot_index(t_val, knots_first, bin_width, n_knots)

Find the knot span index for a data point. Clamped to [1, n_knots].
"""
@inline function _find_knot_index(t_val::T, knots_first::T, bin_width::T, n_knots::Int) where {T}
    return min(floor(Int32, (t_val - knots_first) / bin_width) + one(Int32), n_knots)
end

"""
    _build_numeric_design_matrix(t, knots, bin_width)

Build the numeric design matrix X where X[i,j] = basis_j(t[i]).
Rows are data points, columns are B-spline coefficients.
"""
function _build_numeric_design_matrix(t::Vector{T}, knots::Vector{T},
                                       bin_width::T) where {T<:AbstractFloat}
    n_data = length(t)
    n_knots = length(knots)
    # n_knots knots bound n_knots-1 spans, and a cubic B-spline basis over those
    # spans has (n_knots-1)+3 = n_knots+2 functions -- not n_knots+3.
    #
    # The extra column was structurally unreachable. `_find_knot_index` clamps to
    # n_knots, so the single point at t = max(t) landed one span past the last,
    # at u = 0, and wrote basis[1](0) = 0 into column n_knots+3. That column was
    # therefore always zero, every fit was rank deficient by one, and when the
    # system happened to be square `\` dispatched to `lu` and raised
    # SingularException on the zero pivot. See dev_docs/quant_spline_guard.
    #
    # Clamping one span lower puts that point at u = 1 of the last real span,
    # which writes the same three non-zero values into the same columns: the
    # matrix is the old one with the empty column removed, agreeing elsewhere to
    # within 1e-16. Fitted values are unchanged; the system is now full rank.
    n_spans = max(n_knots - 1, 1)
    n_cols = _n_spline_coeffs(n_knots)
    X = zeros(T, n_data, n_cols)
    knots_first = first(knots)
    basis = _CUBIC_SPLINE_BASIS

    for row in 1:n_data
        knot_idx = min(_find_knot_index(t[row], knots_first, bin_width, n_knots), n_spans)
        u = (t[row] - knots[knot_idx]) / bin_width
        i = 4
        for col in knot_idx:(knot_idx + 3)
            X[row, col] = basis[i](u)
            i -= 1
        end
    end
    return X
end

"""
    _build_piecewise_matrix(knots, bin_width)

Build the symbolic piecewise polynomial matrix. Each row corresponds to a knot
and produces a Polynomial in the local coordinate. Used to convert B-spline
coefficients into piecewise polynomial coefficients.
"""
function _build_piecewise_matrix(knots::Vector{T}, bin_width::T) where {T<:AbstractFloat}
    n_knots = length(knots)
    basis = _CUBIC_SPLINE_BASIS
    n_cols = n_knots + length(basis) - 1
    X = zeros(Polynomial, n_knots, n_cols)
    knots_first = first(knots)

    for row in 1:n_knots
        t_adj = knots[row] + bin_width / 2
        knot_idx = _find_knot_index(t_adj, knots_first, bin_width, n_knots)
        i = 4
        for col in knot_idx:(knot_idx + 3)
            X[row, col] = basis[i]
            i -= 1
        end
    end
    return X
end

"""
    _coeffs_to_spline(c, knots, bin_width, degree, n_knots, _first, _last)

Convert B-spline coefficients to a callable UniformSpline by building
piecewise polynomials and extracting their coefficients.
"""
function _coeffs_to_spline(c::Vector{T}, knots::Vector{T}, bin_width::T,
                            degree::Int, n_knots::Int,
                            _first::T, _last::T) where {T<:AbstractFloat}
    XPoly = _build_piecewise_matrix(knots, bin_width)
    piecewise_polynomials = XPoly * c
    n_total_coeffs = n_knots * (degree + 1)

    # Extract coefficients — Polynomial type drops trailing zeros, so pad them back
    expected_length = degree + 1
    coeff_vecs = [
        let pc = poly.coeffs
            if length(pc) < expected_length
                vcat(pc, zeros(T, expected_length - length(pc)))
            else
                pc[1:expected_length]
            end
        end
        for poly in piecewise_polynomials
    ]
    coeffs = SVector{n_total_coeffs}(vcat(coeff_vecs...))

    return UniformSpline{n_total_coeffs, T}(coeffs, degree, _first, _last, bin_width)
end

"""
    _validate_spline_inputs(u, t, degree, n_knots)

Common input validation for all spline constructors.
"""
function _validate_spline_inputs(u, t, degree, n_knots)
    degree != 3 && error("Non-cubic splines not yet implemented. Use degree = 3")
    n_knots < 3 && error("Need at least 3 knots")
    length(u) != length(t) && error("length(u) must equal length(t)")
end

"""
    _sort_inputs(u, t)

Return (sorted_u, sorted_t), sorting by t if not already sorted.
"""
function _sort_inputs(u::Vector{T}, t::Vector{T}) where {T}
    if issorted(t)
        return u, t
    else
        perm = sortperm(t)
        return u[perm], t[perm]
    end
end

"""
    _setup_knots(sorted_t, n_knots)

Compute knot positions, bin width, and domain bounds from sorted data.
Returns (knots, bin_width, _first, _last).
"""
function _setup_knots(sorted_t::Vector{T}, n_knots::Int) where {T}
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    return knots, bin_width, _first, _last
end

function _make_uniform_spline_basis(t::Vector{T}, n_knots::Integer) where {T<:AbstractFloat}
    n_knots_int = Int(n_knots)
    knots, bin_width, _first, _last = _setup_knots(t, n_knots_int)
    X = _build_numeric_design_matrix(t, knots, bin_width)
    return (
        X = X,
        knots = knots,
        bin_width = bin_width,
        first = _first,
        last = _last,
        n_knots = n_knots_int,
        # Read off the matrix rather than recomputed, so the coefficient count
        # cannot drift from the number of columns callers actually solve against
        # (fit_intensity_mass_error.jl sizes its coefficient vector from this).
        n_coeffs = size(X, 2),
    )
end

function _penalized_spline_system(
    X::AbstractMatrix{T},
    y::AbstractVector{T};
    λ::Real = zero(T),
    weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
    order::Int = 2,
    ridge::Real = 0,
) where {T<:AbstractFloat}
    size(X, 1) == length(y) || throw(DimensionMismatch("X row count must equal length(y)"))

    if weights === nothing
        H = X' * X
        rhs = X' * y
    else
        length(weights) == size(X, 1) || throw(DimensionMismatch("length(weights) must equal number of rows in X"))
        w = T.(weights)
        H = (X .* reshape(w, :, 1))' * X
        rhs = X' * (w .* y)
    end

    if λ != 0
        D = build_difference_matrix(size(X, 2), order)
        H = H + T(λ) * (D' * D)
    end

    if ridge != 0
        H = H + T(ridge) * I
    end

    return H, rhs
end

function _coeffs_to_spline(c::Vector{T}, basis) where {T<:AbstractFloat}
    # The design matrix now has n_knots+2 columns, but `_build_piecewise_matrix`
    # still emits n_knots segments over n_knots+3 coefficients -- the final segment
    # covers only the single point t = last, where the dropped coefficient enters
    # through b0(0) = 0 and so contributes nothing. Pad it back as an explicit zero
    # rather than reshaping the piecewise conversion, which is shared with
    # `build_temp_spline_from_coeffs` and the evaluation path's segment indexing.
    n_expected = basis.n_knots + 3
    padded = length(c) == n_expected ? c : vcat(c, zeros(T, n_expected - length(c)))
    return _coeffs_to_spline(padded, basis.knots, basis.bin_width, 3, basis.n_knots, basis.first, basis.last)
end

#############################################################################
# UniformSpline — unpenalized least squares fit
#############################################################################

function UniformSpline(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
) where {I<:Integer, T<:AbstractFloat}

    _validate_spline_inputs(u, t, degree, n_knots)
    sorted_u, sorted_t = _sort_inputs(u, t)
    basis = _make_uniform_spline_basis(sorted_t, n_knots)
    return _coeffs_to_spline(convert(Vector{T}, basis.X \ sorted_u), basis)
end

#############################################################################
# Callable evaluation (Horner's method, unrolled for degree=3)
#############################################################################

@inline function (s::UniformSpline)(t::U) where {U<:AbstractFloat}
    @fastmath @inbounds begin
        t = max(min(t, s.last), s.first)
        idx = floor(Int32, (t - s.first) / s.bin_width)
        u = (t - (s.first + s.bin_width * idx)) / s.bin_width
        coeff = idx * (s.degree + 1) + 1
        # Horner's method: c0 + u*(c1 + u*(c2 + u*c3))
        x = muladd(u, s.coeffs[coeff + 3], s.coeffs[coeff + 2])
        x = muladd(u, x, s.coeffs[coeff + 1])
        x = muladd(u, x, s.coeffs[coeff])
        return x
    end
end

#############################################################################
# build_temp_spline_from_coeffs — used by SplineQuadModel
#############################################################################

"""
    build_temp_spline_from_coeffs(c, knots, bin_width, spline_basis, degree, _first, _last, n_knots)

Build a UniformSpline from pre-computed B-spline coefficients.
Used by SplineQuadModel for residual calculation during quad fitting.
"""
function build_temp_spline_from_coeffs(
    c::Vector{T},
    knots::Vector{T},
    bin_width::T,
    spline_basis::NTuple{4, Polynomial},  # kept for API compatibility
    degree::Int,
    _first::T,
    _last::T,
    n_knots::Int
) where {T<:AbstractFloat}
    return _coeffs_to_spline(c, knots, bin_width, degree, n_knots, _first, _last)
end

#############################################################################
# build_difference_matrix — P-spline penalty
#############################################################################

"""
    build_difference_matrix(n_coeffs::Int, order::Int=2)

Build k-th order difference matrix for penalized B-splines (P-splines).

- Order 1: Penalizes differences (cᵢ₊₁ - cᵢ)
- Order 2: Penalizes curvature (cᵢ₊₂ - 2cᵢ₊₁ + cᵢ) [RECOMMENDED]
- Order 3: Penalizes jerk (cᵢ₊₃ - 3cᵢ₊₂ + 3cᵢ₊₁ - cᵢ)

Reference: Eilers & Marx (1996), Statistical Science 11(2), 89-121.
"""
function build_difference_matrix(n_coeffs::Int, order::Int=2)
    if order < 1 || order > 3
        error("Difference order must be 1, 2, or 3")
    end

    if order == 1
        n_rows = n_coeffs - 1
        D = zeros(Float64, n_rows, n_coeffs)
        for i in 1:n_rows
            D[i, i] = -1.0
            D[i, i+1] = 1.0
        end
        return D
    elseif order == 2
        D1_n = build_difference_matrix(n_coeffs, 1)
        D1_nm1 = build_difference_matrix(n_coeffs - 1, 1)
        return D1_nm1 * D1_n
    else
        D2 = build_difference_matrix(n_coeffs, 2)
        D1_nm2 = build_difference_matrix(n_coeffs - 2, 1)
        return D1_nm2 * D2
    end
end

#############################################################################
# UniformSplinePenalized — P-spline with difference penalty
#############################################################################

"""
    UniformSplinePenalized(u, t, degree, n_knots, λ=0.1, order=2)

Fit uniform B-spline with difference penalty regularization (P-splines).
Minimizes: ||Xc - u||² + λ||D^k c||²

# Arguments
- `u::Vector{T}`: Target values (iRT in RT alignment context)
- `t::Vector{T}`: Input values (RT)
- `degree::I`: Polynomial degree (must be 3)
- `n_knots::I`: Number of knots (typically 5)
- `λ::T`: Penalty parameter (0=no penalty, 0.1-1.0 recommended)
- `order::Int`: Difference penalty order (default 2, penalizes curvature)

Reference: Eilers & Marx (1996), Statistical Science 11(2), 89-121.
"""
function UniformSplinePenalized(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
    λ::T = T(0.1),
    order::Int = 2
) where {I<:Integer, T<:AbstractFloat}

    _validate_spline_inputs(u, t, degree, n_knots)
    λ < 0 && error("Penalty parameter λ must be non-negative")

    sorted_u, sorted_t = _sort_inputs(u, t)
    basis = _make_uniform_spline_basis(sorted_t, n_knots)
    c = if λ == 0
        basis.X \ sorted_u
    else
        H, rhs = _penalized_spline_system(basis.X, sorted_u; λ=λ, order=order)
        H \ rhs
    end
    return _coeffs_to_spline(convert(Vector{T}, c), basis)
end

#############################################################################
# calculate_adaptive_iterations — RANSAC iteration count
#############################################################################

"""
    calculate_adaptive_iterations(outlier_ratio, sample_size, confidence=0.99)

Calculate RANSAC iterations needed: k = log(1-p) / log(1-(1-ε)^s).
Result clamped to [20, 200].
"""
function calculate_adaptive_iterations(
    outlier_ratio::T,
    sample_size::Int,
    confidence::T = T(0.99)
) where T<:AbstractFloat
    if outlier_ratio <= 0
        return 20
    end
    prob_all_inliers = (1 - outlier_ratio)^sample_size
    if prob_all_inliers < 1e-10
        return 200
    end
    k = log(1 - confidence) / log(1 - prob_all_inliers)
    return max(20, min(200, ceil(Int, k)))
end

#############################################################################
# fit_spline_with_ransac — robust fitting with RANSAC + P-spline
#############################################################################

"""
    fit_spline_with_ransac(u, t, degree, n_knots, λ=0.1, order=2; kwargs...)

Fit uniform B-spline with RANSAC robust fitting and P-spline regularization.
Handles datasets with up to 50% outlier contamination.

Algorithm:
1. Initial P-spline fit → compute MAD-based inlier threshold
2. Randomly sample subsets, fit P-splines, count inliers
3. Select model with most inliers
4. Refit on all inliers from best model

Reference: Fischler & Bolles (1981), CACM 24(6), 381-395.
"""
function fit_spline_with_ransac(
    u::Vector{T},
    t::Vector{T},
    degree::Int,
    n_knots::Int,
    λ::T = T(0.1),
    order::Int = 2;
    ransac_iterations::Int = 50,
    ransac_sample_size::Int = 30,
    ransac_threshold_factor::T = T(2.0),
    ransac_seed::Union{Int,Nothing} = nothing
) where {T<:AbstractFloat}

    _validate_spline_inputs(u, t, degree, n_knots)

    spline_initial = UniformSplinePenalized(u, t, degree, n_knots, λ, order)
    sorted_u, sorted_t = _sort_inputs(u, t)

    if ransac_seed !== nothing
        Random.seed!(ransac_seed)
    end

    n_data = length(sorted_u)
    n_coeffs = n_knots + degree
    sample_size = max(n_coeffs + 2, min(ransac_sample_size, div(n_data, 2)))

    # Adaptive threshold from initial fit
    predicted_initial = [spline_initial(ti) for ti in sorted_t]
    residuals_initial = abs.(sorted_u .- predicted_initial)
    inlier_threshold = median(residuals_initial) * ransac_threshold_factor

    # Adaptive iteration count
    outlier_ratio = sum(residuals_initial .> inlier_threshold) / n_data
    n_iters = min(ransac_iterations,
                  calculate_adaptive_iterations(T(outlier_ratio), sample_size))

    # RANSAC main loop
    best_inlier_count = 0
    best_spline = spline_initial
    best_inlier_mask = trues(n_data)

    for _ in 1:n_iters
        sample_idx = randperm(n_data)[1:sample_size]

        spline_sample = try
            UniformSplinePenalized(sorted_u[sample_idx], sorted_t[sample_idx],
                                   degree, n_knots, λ, order)
        catch
            continue
        end

        residuals_sample = abs.(sorted_u .- [spline_sample(ti) for ti in sorted_t])
        inlier_mask = residuals_sample .<= inlier_threshold
        inlier_count = sum(inlier_mask)

        if inlier_count > best_inlier_count
            best_inlier_count = inlier_count
            best_spline = spline_sample
            best_inlier_mask = inlier_mask
        end

        inlier_count / n_data > 0.95 && break
    end

    # Refit on all inliers
    if best_inlier_count >= n_coeffs + 2
        return try
            UniformSplinePenalized(sorted_u[best_inlier_mask], sorted_t[best_inlier_mask],
                                   degree, n_knots, λ, order)
        catch
            best_spline
        end
    else
        return best_spline
    end
end
