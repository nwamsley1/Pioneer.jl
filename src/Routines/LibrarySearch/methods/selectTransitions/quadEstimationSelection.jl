"""
    _select_transitions_impl!(transitions, ::QuadEstimationTransitionSelection, ...)

Implementation for quadrupole transmission estimation. For each precursor,
processes multiple isotope indices to estimate transmission characteristics.

Returns (transition_idx, n) where:
- transition_idx: Number of transitions added
- n: Number of precursor-isotope combinations processed
"""
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::QuadEstimationTransitionSelection,
    transition_idx::Int64,
    library_fragment_lookup::LibraryFragmentLookup,
    prec_idxs::AbstractVector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    frag_mz_bounds::Tuple{Float32, Float32};
    isotope_range::UnitRange{Int64} = 0:2,  # Make isotope range configurable
    block_size::Int64 = 10000
)
    n = 0
    
    for prec_idx in prec_idxs
        prec_props = getPrecursorProperties(
            prec_idx, prec_mzs, prec_charges, prec_sulfur_counts
        )
        
        for prec_iso_idx in isotope_range
            n += 1
            
            transition_idx = @inline fillTransitionListForQuadEstimation!(
                transitions,
                getPrecFragRange(library_fragment_lookup, prec_idx),
                getFragments(library_fragment_lookup),
                prec_props.mz,
                prec_props.charge,
                prec_props.sulfur_count,
                transition_idx,
                precursor_transmission,
                isotopes,
                prec_iso_idx,
                iso_splines,
                frag_mz_bounds,
                block_size
            )
        end
    end

    return transition_idx, n
end
