# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
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

"Normalised Gaussian over -h:h with h = ceil(extent * sigma); sigma <= 0 gives [1.0]."
function gauss_kernel(sigma::Real, extent::Real = 3.0)
    sigma <= 0 && return [1.0]
    h = ceil(Int, extent * sigma)
    k = [exp(-(x^2) / (2sigma^2)) for x in -h:h]
    k ./ sum(k)
end
"""
IM kernel of a level: Gaussian scaled to sum to the stride when `sum_scale`.

Why scale: a slice stands for `stride` scans. With weights summing to 1 its intensity would be a per-scan average;
summing to `stride` makes it approximately the ion's total over those scans, so intensities stay on the scale of
raw counts whatever the stride.
"""
im_kernel(lp::LevelParams) = gauss_kernel(lp.im_sigma, lp.kernel_extent) .* (lp.sum_scale ? lp.stride : 1)
mz_kernel(lp::LevelParams) = gauss_kernel(lp.mz_sigma, lp.kernel_extent)
