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

#=
Iteration scheme implementations.
=#

"""
    get_iteration_rounds(scheme::IterationScheme) -> Vector{Int}

Get the number of boosting rounds for each iteration.
"""
function get_iteration_rounds(scheme::FixedIterationScheme)
    return scheme.rounds
end

function get_iteration_rounds(scheme::SinglePassScheme)
    return [scheme.num_rounds]
end

"""
    get_iterations_count(scheme::IterationScheme) -> Int

Get total number of training iterations.
"""
get_iterations_count(scheme::FixedIterationScheme) = length(scheme.rounds)
get_iterations_count(scheme::SinglePassScheme) = 1
