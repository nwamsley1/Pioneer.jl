# Core implementation for MassErrEstimationStrategy
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::MassErrEstimationStrategy,
    transition_idx::Int64,
    library_fragment_lookup::T,  # Generic type parameter
    scan_to_prec_idx::UnitRange{Int64},
    precursors_passed_scoring::Vector{UInt32};
    max_rank::Int64 = 5,
    block_size::Int64 = 10000) where {T <: LibraryFragmentLookup}

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
    return transition_idx
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
        transitions[transition_idx] = convert_to_detailed(frag, getKnots(lookup), getNCE(lookup))
        ensureTransitionCapacity!(transitions, transition_idx, block_size)
    end
    return transition_idx
end