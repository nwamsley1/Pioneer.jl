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

# src/utils/math.jl
"""
Calculate parts per million difference between two masses.
"""
function get_ppm(mass1::T, mass2::T)::T where T <: AbstractFloat
    return abs(mass1 - mass2) / ((mass1 + mass2) / 2.0) * 1e6
end

"""
Get PPM tolerance window for a given mass.
"""
function get_ppm_window(mass::T, tolerance_ppm::T)::Tuple{T,T} where T <: AbstractFloat
    window = mass * tolerance_ppm / 1e6
    return (mass - window, mass + window)
end
