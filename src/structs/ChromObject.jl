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

abstract type ChromObject end
struct MS2ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    scan_idx::UInt32
    precursor_idx::UInt32
end
struct MS1ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    m0::Bool
    n_iso::UInt8
    scan_idx::UInt32
    precursor_idx::UInt32
end

function growChromObjects!(chromatograms::Vector{ChromObject}, block_size::Int64)
    chromatograms = append!(chromatograms, Vector{ChromObject}(undef, block_size))
end

