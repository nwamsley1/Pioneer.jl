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

abstract type IsotopeTraceType end

struct CombineTraces <: IsotopeTraceType
    min_ft::Float32
end

seperateTraces(itt::CombineTraces) = false

function getPsmGroupbyCols(itt::CombineTraces)
    return [:precursor_idx]
end

struct SeperateTraces <: IsotopeTraceType
end
 
seperateTraces(itt::SeperateTraces) = true

function getPsmGroupbyCols(itt::SeperateTraces)
    return [:precursor_idx,:isotopes_captured]
end




