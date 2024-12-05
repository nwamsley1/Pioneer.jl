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
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
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
        prec_charge = prec_charges[prec_idx]
        prec_mz = prec_mzs[prec_idx]
        prec_sulfur_count = prec_sulfur_counts[prec_idx]
        for prec_iso_idx in isotope_range
            n += 1
            transition_idx = fill_quad_transitions!(
                transitions,
                lookup,
                prec_idx,
                prec_mz,
                prec_charge,
                prec_sulfur_count,
                prec_iso_idx,
                transition_idx,
                precursor_transmission,
                isotopes,
                iso_splines,
                frag_mz_bounds,
                getSplineData(lookup, prec_charge, prec_mz),
                block_size
            )
        end
    end

    return transition_idx
end


function fill_quad_transitions!(
    transitions::Vector{DetailedFrag{Float32}},
    lookup::LibraryFragmentLookup,
    prec_idx::UInt32,
    prec_mz::Float32,
    prec_charge::UInt8,
    prec_sulfur_count::UInt8,
    prec_iso_idx::Int64,
    transition_idx::Int64,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    iso_splines::IsotopeSplineModel,
    frag_mz_bounds::Tuple{Float32, Float32},
    spline_data::G,
    block_size::Int64,
) where {G<:IntensityDataType}
    # Reset and set up precursor transmission
    fill!(precursor_transmission, zero(Float32))
    precursor_transmission[prec_iso_idx + 1] = one(Float32)

    # Calculate isotope factor
    iso_fac = one(Float32) / (
        iso_splines(
            min(Int64(prec_sulfur_count), 5),
            prec_iso_idx,
            prec_mz * prec_charge
        )
    )

    # Process each fragment
    for frag_idx in getPrecFragRange(lookup, prec_idx)
        frag = getFrag(lookup, frag_idx)
        getRank(frag) > 5 && continue

        # Get isotope abundances
        getFragIsotopes!(
            PartialPrecCapture(),
            isotopes,
            precursor_transmission,
            (prec_iso_idx, prec_iso_idx),
            0:2,
            iso_splines,
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            frag,
            spline_data
        )

        # Add transitions for mono and mono+1 isotopes
        transition_idx = add_quad_transitions!(
            transitions,
            transition_idx,
            frag,
            isotopes,
            iso_fac,
            frag_mz_bounds,
            block_size
        )
    end

    return transition_idx
end

function add_quad_transitions!(
    transitions::Vector{DetailedFrag{Float32}},
    transition_idx::Int64,
    frag::F,
    isotopes::Vector{Float32},
    iso_fac::Float32,
    frag_mz_bounds::Tuple{Float32, Float32},
    block_size::Int64
) where {F <: AltimeterFragment}
    NEUTRON = 1.00335f0

    for iso_idx in 0:2
        frag_mz = getMz(frag) + iso_idx * NEUTRON/getFragCharge(frag)
        
        if frag_mz < first(frag_mz_bounds) || frag_mz > last(frag_mz_bounds)
            continue
        end

        transition_idx += 1
        transitions[transition_idx] = DetailedFrag(
            UInt32((getPID(frag) - 1) * 3 + (iso_idx + 1)),
            Float32(frag_mz),
            Float16(isotopes[iso_idx + 1] * iso_fac),
            getIonType(frag),
            isY(frag),
            isB(frag),
            isP(frag),
            iso_idx > 0,
            getFragCharge(frag),
            getIonPosition(frag),
            getPrecCharge(frag),
            getRank(frag),
            getSulfurCount(frag)
        )

        if transition_idx >= length(transitions)
            append!(transitions, [DetailedFrag{Float32}() for _ in 1:block_size])
        end
    end

    return transition_idx
end