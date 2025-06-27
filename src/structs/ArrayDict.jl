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

mutable struct ArrayDict{I<:Unsigned,C<:Real}
    keys::Vector{I}
    vals::Vector{C}
    size::Int64
    function ArrayDict(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(zeros(I, size), zeros(C, size), 0)
    end
end

function update!(c::ArrayDict{I,C}, key::I, val::C) where {I,C<:Unsigned}
    c.size += 1
    c.vals[key] = val
    c.keys[c.size] = key
end

function reset!(c::ArrayDict{I,C}) where {I,C<:Unsigned}
    @turbo for i in range(1, c.size)
        c.vals[c.keys[i]] = zero(I)
        c.keys[i] = zero(C)
    end
    c.size = 0
end

function Base.getindex(c::ArrayDict{I,C}, i::Ti) where {I,C<:Unsigned, Ti<:Integer}
    c.vals[i]
end