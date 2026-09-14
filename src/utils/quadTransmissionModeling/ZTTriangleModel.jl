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
    ZTTriangleModel{T}

Quadrupole transmission for a SWEPT quad (scanning DIA / "ZT"). Replaces the Razo model on such
files: Razo is fit from the isotope ratio of a precursor seen in two abutting windows, which has
no analogue here — on a swept quad a precursor is seen in ~2k+1 bins of one meta-scan, and its
deconvolved weight in each bin is already a direct sample of transmission.

Transmission is a symmetric triangle falling linearly to zero at `±half_width_mz`:

    T(ionMz) = max(0, 1 - |ionMz - centerMz| / half_width_mz)

`half_width_mz` ("h") is measured per file: per meta-scan least squares of
`w = a - b*|Δ|` gives `h = a/b`, which is scale-free because `a` and `b` share the precursor's
unknown abundance. Measured on Sciex ZT 5 Da/5 min: h = 6.17 Da, and `h / bin_step = 6.04` —
which is where `metascan_k = 6` comes from physically, and why k = 6 also holds on the 10 Da
method (its Q1 step and window both double, leaving the ratio unchanged).

The apex is NOT at the monoisotopic m/z. Transmission peaks at the intensity-weighted centre of
the transmitted isotope envelope, which sits above M0 and moves further as peptides get heavier:
measured apex offset runs -0.02 Da at m/z 400 to +0.31 Da at m/z 900. Callers pass a
precursor-specific `center_mz` already carrying that offset.
"""
struct ZTTriangleModel{T<:AbstractFloat} <: QuadTransmissionModel
    half_width_mz::T
end

struct ZTTriangleFunction{T<:AbstractFloat} <: QuadTransmissionFunction
    center_mz::T
    half_width_mz::T
end

function getQuadTransmissionFunction(qtm::ZTTriangleModel{T}, centerMz::T,
                                     isolationWidthMz::T) where {T<:AbstractFloat}
    ZTTriangleFunction(centerMz, qtm.half_width_mz)
end

function getQuadTransmissionBounds(qtm::ZTTriangleModel{T}, centerMz::T,
                                   isolationWidthMz::T) where {T<:AbstractFloat}
    return T(centerMz - isolationWidthMz/2), T(centerMz + isolationWidthMz/2)
end

# Support is exactly the triangle's base: beyond +/-h transmission is zero, so there is no
# reason to consider a precursor there. This is what lets `metascan_k` be DERIVED rather than
# configured -- the expansion needs to reach exactly as far as the support does.
getPrecMinBound(f::ZTTriangleFunction{T}) where {T<:AbstractFloat} = f.center_mz - f.half_width_mz
getPrecMaxBound(f::ZTTriangleFunction{T}) where {T<:AbstractFloat} = f.center_mz + f.half_width_mz

@inline function (f::ZTTriangleFunction{T})(ionMz::U) where {T,U<:AbstractFloat}
    d = abs(U(ionMz) - U(f.center_mz))
    h = U(f.half_width_mz)
    return d >= h ? zero(U) : one(U) - d / h
end
