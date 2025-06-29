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

# Core implementation for MassErrEstimationStrategy
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::MassErrEstimationStrategy,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    library_fragment_lookup::L,  # Generic type parameter
    scan_to_prec_idx::UnitRange{Int64},
    precursors_passed_scoring::Vector{UInt32};
    max_rank::Int64 = 5,
    block_size::Int64 = 10000) where {L <: LibraryFragmentLookup}

    for i in scan_to_prec_idx
        prec_idx = precursors_passed_scoring[i]
        for frag_idx in getPrecFragRange(library_fragment_lookup, prec_idx)
            transition_idx = process_fragment!(
                transitions, 
                transition_idx,
                library_fragment_lookup,
                frag_idx,
                max_rank,
                block_size
            )
        end
    end
    return transition_idx, 0
end

# Specialized fragment processing for standard lookup
function process_fragment!(
    transitions::Vector{DetailedFrag{Float32}},
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    frag_idx::Integer,
    max_rank::Int64,
    block_size::Int64)
    frag = getFrag(lookup, frag_idx)
    if getRank(frag) <= max_rank
        transition_idx += 1
        transitions[transition_idx] = convert_to_detailed(frag, getSplineData(lookup))
        ensureTransitionCapacity!(transitions, transition_idx, block_size)
    end
    return transition_idx
end