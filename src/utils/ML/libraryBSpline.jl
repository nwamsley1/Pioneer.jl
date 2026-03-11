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

"""
    splevl(x, knots, c, k)

Evaluate cubic B-spline at point x using the iterative de Boor algorithm.
Knots `knots`, coefficients `c`, degree `k` (must be 3).

Replaces the recursive Cox-de Boor approach (B() + sum) with a single-pass
triangular iteration: O(k²) = 6 multiply-adds vs 32 recursive calls.
"""
@inline function splevl(x::T, knots::NTuple{N,T}, c::NTuple{M,T}, k::Int) where {M,N,T<:AbstractFloat}
    # Find knot span: j where knots[j] ≤ x < knots[j+1]
    j = 0
    @inbounds for idx in 1:(N-1)
        if knots[idx] ≤ x < knots[idx+1]
            j = idx
            break
        end
    end

    # x outside all knot spans → all basis functions are zero
    j == 0 && return zero(T)

    # Safe coefficient access: zero outside [1, M]
    @inline _getc(i) = (1 ≤ i ≤ M) ? @inbounds(c[i]) : zero(T)
    # Safe knot access: clamp to [1, N]
    # Out-of-bounds accesses only occur when multiplied by zero coefficients
    @inline _getk(i) = @inbounds knots[clamp(i, 1, N)]

    # Initialize d[1..4] with active coefficients (zero-padded outside range)
    d1 = _getc(j - 3)
    d2 = _getc(j - 2)
    d3 = _getc(j - 1)
    d4 = _getc(j)

    # de Boor triangular iterations, fully unrolled for k=3
    # r = 1: three updates
    denom = _getk(j + 3) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    denom = _getk(j + 2) - _getk(j - 1)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 1)) / denom)
    d3 = (one(T) - α) * d2 + α * d3

    denom = _getk(j + 1) - _getk(j - 2)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 2)) / denom)
    d2 = (one(T) - α) * d1 + α * d2

    # r = 2: two updates
    denom = _getk(j + 2) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    denom = _getk(j + 1) - _getk(j - 1)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j - 1)) / denom)
    d3 = (one(T) - α) * d2 + α * d3

    # r = 3: one update
    denom = _getk(j + 1) - _getk(j)
    α = ifelse(denom == zero(T), zero(T), (x - _getk(j)) / denom)
    d4 = (one(T) - α) * d3 + α * d4

    return d4
end

"""
    getSplineQuadrature(dtype::Type)

Get Quadrature nodes and weights as static vector to integrate B-Splines
"""
function getSplineQuadrature(dtype::Type)
    x, w = gausslegendre(20)
    return SVector{20}(dtype.(x)), SVector{20}(dtype.(w))
end

"""
    splevl(knotw, c, d, gqw, gqx)

Use Fast Gaussian Quadrature to compute the definite integral of a B-spline
over its domain
"""
function splint(knots::NTuple{N, T},
                c::NTuple{M, T},
                d::Int,
                gqx::SVector{20, T},
                gqw::SVector{20, T}) where {M,N,T<:AbstractFloat}
    i_eval = zero(T)
    @inbounds for i in 1:20
        i_eval += splevl(gqx[i], knots, c, d)*gqw[i]
    end
    return i_eval
end

"""
    getSplineQuadrature(dtype::Type, x0::AbstractFloat, x1::AbstractFloat)

Get quadrature weights and nodes as static arrays. Defaults to 20. 
"""
function getSplineQuadrature(dtype::Type, x0::AbstractFloat, x1::AbstractFloat)
    x, w = gausslegendre(20)
    ws = (x1 - x0)/2
    for i in range(1, length(x))
        x[i] = x0 + (x[i] + 1)*(x1- x0)/2
        w[i] = ws*w[i]
    end
    return SVector{20}(dtype.(x)), SVector{20}(dtype.(w))
end

function getSplineAreas(knots::NTuple{N, T}, 
                        coefficients::AbstractVector{NTuple{M, T}},
                        d::Int,
                        gqx::SVector{20, T},
                        gqw::SVector{20, T}) where {M,N,T<:AbstractFloat}
    intensities = Vector{T}(undef, length(coefficients))
    for i in range(1, length(coefficients))
        intensities[i] = (
            splint(knots, coefficients[i], d, gqx, gqw)
        )
    end
    return intensities
end
#=
knots = (6.0f0, 13.0f0, 20.0f0, 27.0f0, 34.0f0, 41.0f0, 48.0f0, 55.0f0)
c = (1.6181915f-6, 7.382022f-6, 7.887343f-5, 0.00023642876f0)
t = 37.0f0
degree = 3
@btime splevl(t, knots, c, degree)
gqx, gqw = getSplineQuadrature(Float32, 20.0, 40.0);
@btime splint(knots, c, 3, gqx, gqw)
=#