# Define a protocol for transition selection through traits
abstract type TransitionSelectionStrategy end
struct StandardTransitionSelection <: TransitionSelectionStrategy end
struct MassErrEstimationStrategy <: TransitionSelectionStrategy end
struct RTIndexedTransitionSelection <: TransitionSelectionStrategy end
struct QuadEstimationTransitionSelection <: TransitionSelectionStrategy end


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
    transition_idx = _select_transitions_impl!(transitions, strategy, prec_estimation_type, transition_idx, common_args...; kwargs...)
    
    # Common cleanup/sorting code
    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
          
    return transition_idx
end