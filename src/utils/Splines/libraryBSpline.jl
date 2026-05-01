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

