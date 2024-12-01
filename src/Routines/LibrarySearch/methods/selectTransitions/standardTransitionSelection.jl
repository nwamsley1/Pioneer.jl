function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::StandardTransitionSelection,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    scan_to_prec_idx::UnitRange{Int64},
    precursors_passed_scoring::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    abreviate_precursor_calc::Bool,
    iRT::Float32,
    iRT_tol::Float32,
    frag_mz_bounds::Tuple{Float32, Float32};
    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
    block_size::Int64 = 10000
    )

    for i in scan_to_prec_idx
        # Get precursor properties
        prec_idx = precursors_passed_scoring[i]
        prec_charge = prec_charges[prec_idx]
        prec_mz = prec_mzs[prec_idx]
        
        # Check retention time tolerance
        if abs(prec_irts[prec_idx] - iRT) > iRT_tol
            continue
        end
        
        # Handle isotope errors
        mz_low = getPrecMinBound(quad_transmission_func) - NEUTRON*first(isotope_err_bounds)/prec_charge
        mz_high = getPrecMaxBound(quad_transmission_func) + NEUTRON*last(isotope_err_bounds)/prec_charge
        
        if (prec_mz < mz_low) | (prec_mz > mz_high)
            continue
        end
        
        prec_sulfur_count = prec_sulfur_counts[prec_idx]
        nce = getNCE(lookup, prec_charge, prec_mz)
        
        # Fill transition list using spline-specific implementation
        transition_idx = @inline fillTransitionList!(
            transitions,
            getPrecFragRange(lookup, prec_idx),
            getFragments(lookup),
            nce,
            getKnots(lookup),
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            transition_idx,
            quad_transmission_func,
            precursor_transmission,
            isotopes,
            n_frag_isotopes,
            max_frag_rank,
            abreviate_precursor_calc,
            iso_splines,
            frag_mz_bounds,
            block_size
        )
    end
    
    return transition_idx
end
