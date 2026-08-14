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
    FragBoundModel(low_mass, high_mass[, lo_limit, hi_limit])

The fragment m/z window a precursor's fragments must fall in, as a pair of
polynomials evaluated at the precursor m/z.

`lo_limit` and `hi_limit` clamp the result. They default to `(-Inf, Inf)`, which
is a no-op, so a two-argument model behaves exactly as it always has.

The clamp matters as soon as a bound has a slope. A ceiling of `2.0·x + 10`
evaluated at a 1250 m/z precursor gives 2510 m/z, and no Orbitrap records there
— without a limit, a wide precursor range silently admits fragments above
anything the instrument can see. Fitting a line on one file and evaluating it
outside that file's range is a mistake this codebase has made before: an
identical ceiling fit in `DownloadSpecLib/fit.jl` once reported a 4626 m/z
ceiling, fixed in `f35564b52` by holding the line flat beyond its domain.

The auto path is safe from this only by coincidence — it derives the precursor
range from the same file it fits, so the polynomial is never evaluated outside
its own domain. A manually configured slope has no such coupling, so the clamp
lives here rather than at the call site.
"""
struct FragBoundModel
    low_mass::ImmutablePolynomial{Float32}
    high_mass::ImmutablePolynomial{Float32}
    lo_limit::Float32
    hi_limit::Float32

    # Deliberately untyped in the element: the auto path fits in Float64 and
    # relied on the default constructor converting on the way in. `new` still
    # converts to the field types, so this keeps that.
    function FragBoundModel(
        low_mass::ImmutablePolynomial,
        high_mass::ImmutablePolynomial,
        lo_limit::Real = -Inf32,
        hi_limit::Real = Inf32,
    )
        return new(low_mass, high_mass, Float32(lo_limit), Float32(hi_limit))
    end
end

function (fbm::FragBoundModel)(x::AbstractFloat)
    lo = max(fbm.low_mass(x), fbm.lo_limit)
    hi = min(fbm.high_mass(x), fbm.hi_limit)
    return lo, hi
end
