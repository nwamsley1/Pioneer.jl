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
    prepare_spline_fractions(x, knots)

Find the cubic spline span and calculate its six de Boor interpolation
fractions once for reuse across fragment coefficient tuples. The active knot
spans must have nonzero widths.
"""
@inline function prepare_spline_fractions(
        x::T, knots::NTuple{N,T}) where {N,T<:AbstractFloat}
    j = 0
    @inbounds for idx in 1:(N - 1)
        if knots[idx] ≤ x < knots[idx + 1]
            j = idx
            break
        end
    end
    j == 0 && return PreparedSplineFractions(UInt8(0), ntuple(_ -> zero(T), Val(6)))

    @inline getk(i) = @inbounds knots[clamp(i, 1, N)]
    @inline alpha(num, denom) = num / denom

    return PreparedSplineFractions(
        UInt8(j),
        (
            alpha(x - getk(j),     getk(j + 3) - getk(j)),
            alpha(x - getk(j - 1), getk(j + 2) - getk(j - 1)),
            alpha(x - getk(j - 2), getk(j + 1) - getk(j - 2)),
            alpha(x - getk(j),     getk(j + 2) - getk(j)),
            alpha(x - getk(j - 1), getk(j + 1) - getk(j - 1)),
            alpha(x - getk(j),     getk(j + 1) - getk(j)),
        ),
    )
end


"""
    splevl_prepared(coefficients, prepared)

Evaluate four spline coefficients using a previously prepared knot span and
interpolation fractions. The cubic de Boor iteration is fully unrolled.
"""
@inline function splevl_prepared(
        c::NTuple{4,T}, prepared::PreparedSplineFractions{T}) where {T<:AbstractFloat}
    j = Int(prepared.span)
    j == 0 && return zero(T)
    a = prepared.alpha

    @inline getc(i) = (1 ≤ i ≤ 4) ? @inbounds(c[i]) : zero(T)
    d1 = getc(j - 3)
    d2 = getc(j - 2)
    d3 = getc(j - 1)
    d4 = getc(j)

    d4 = (one(T) - a[1]) * d3 + a[1] * d4
    d3 = (one(T) - a[2]) * d2 + a[2] * d3
    d2 = (one(T) - a[3]) * d1 + a[3] * d2
    d4 = (one(T) - a[4]) * d3 + a[4] * d4
    d3 = (one(T) - a[5]) * d2 + a[5] * d3
    return (one(T) - a[6]) * d3 + a[6] * d4
end
