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

# Define a protocol for transition selection through traits
abstract type TransitionSelectionStrategy end
struct StandardTransitionSelection <: TransitionSelectionStrategy end
struct MassErrEstimationStrategy <: TransitionSelectionStrategy end
struct RTIndexedTransitionSelection <: TransitionSelectionStrategy end
struct QuadEstimationTransitionSelection <: TransitionSelectionStrategy end
"""
    PrecEstimation
Abstract type to specify precursor isotopes estimation strategy.
Used to determine how fragment ion isotope patterns are calculated.
"""
abstract type PrecEstimation end

"""
    PartialPrecCapture <: PrecEstimation

Strategy for calculating fragment isotope patterns when only a portion of the precursor 
isotope pattern is captured in the isolation window. Uses a more detailed calculation
that accounts for the partial capture of precursor isotopologues.
"""
struct PartialPrecCapture <: PrecEstimation end

"""
    FullPrecCapture <: PrecEstimation

Strategy for calculating fragment isotope patterns when the full precursor isotope pattern
is captured in the isolation window. Uses a simplified calculation that assumes complete
capture of relevant precursor isotopologues.
"""
struct FullPrecCapture <: PrecEstimation end


# Helper function to ensure capacity
function ensureTransitionCapacity!(
    transitions::Vector{DetailedFrag{Float32}}, 
    transition_idx::Int64,
    block_size::Int64)
    
    if transition_idx + 1 > length(transitions)
        append!(transitions, [DetailedFrag{Float32}() for _ in 1:block_size])
    end
end

function getPrecursorProperties(
    prec_idx::UInt32,
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8}
)
    return (
        mz = prec_mzs[prec_idx],
        charge = prec_charges[prec_idx],
        sulfur_count = prec_sulfur_counts[prec_idx]
    )
end


# Core function that handles common setup and cleanup
function selectTransitions!(
    transitions::Vector{DetailedFrag{Float32}},
    strategy::TransitionSelectionStrategy,
    prec_estimation_type::PrecEstimation,
    common_args...;
    kwargs...)
    
    # Common setup code here
    transition_idx = 0
    
    # Delegate to specific implementation
    transition_idx, n_precursors = _select_transitions_impl!(transitions, strategy, prec_estimation_type, transition_idx, common_args...; kwargs...)
    
    # Common cleanup/sorting code
    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
          
    return transition_idx, n_precursors
end

# Core function that handles common setup and cleanup
#=
function selectIsotopes!(
    isotopes::Vector{Isotope{Float32}},
    prec_estimation_type::PrecEstimation,
    common_args...;
    kwargs...)
    
    # Common setup code here
    transition_idx = 0
    
    # Delegate to specific implementation
    transition_idx, n_precursors = _select_transitions_impl!(transitions, strategy, prec_estimation_type, transition_idx, common_args...; kwargs...)
    
    # Common cleanup/sorting code
    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
          
    return transition_idx, n_precursors
end
=#