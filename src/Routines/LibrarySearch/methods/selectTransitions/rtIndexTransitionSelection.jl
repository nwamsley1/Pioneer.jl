"""
    _select_transitions_impl!(transitions, ::RTIndexedTransitionSelection, ...)

Implement transition selection based on retention time indexing. Handles both filtered 
and unfiltered precursor selection through an optional precursors_passing parameter.

Returns (transition_idx, n, precs_temp_size) where:
- transition_idx: Number of transitions added
- n: Number of precursors processed
- precs_temp_size: Final size of temporary precursor array
"""
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::RTIndexedTransitionSelection,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    precs_temp::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    rt_index::retentionTimeIndex{Float32, Float32},
    rt_start_idx::Int64,
    rt_stop_idx::Int64,
    frag_mz_bounds::Tuple{Float32, Float32};
    precursors_passing::Union{Set{UInt32}, Nothing} = nothing,
    isotope_err_bounds::Tuple{I, I} = (3, 1),
    block_size::Int64 = 10000
) where {I<:Integer}
    n = 0
    min_prec_mz, max_prec_mz = getQuadrupoleBounds(quad_transmission_func)
    precs_temp_size = 0
    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec
        
        # Get precursor range within isolation window
        start, stop = getPrecursorWindowRange(
            precs, min_prec_mz, max_prec_mz, isotope_err_bounds
        )

        for i in start:stop
            prec_idx = first(precs[i])
            
            # Optional precursor filtering
            if !isnothing(precursors_passing) && prec_idx ∉ precursors_passing
                continue
            end

            # Get and validate precursor properties
            
            #prec_props = getPrecursorProperties(
            #    prec_idx, prec_mzs, prec_charges, prec_sulfur_counts
            #)
            prec_charge = prec_charges[prec_idx]
            prec_mz = prec_mzs[prec_idx] 
            prec_sulfur_count = prec_sulfur_counts[prec_idx]
   
            # Check quadrupole bounds
            if !withinQuadrupoleBounds(
                prec_mz, prec_charge, 
                min_prec_mz, max_prec_mz, isotope_err_bounds
            )
                continue
            end

            # Update counters and temp storage
            precs_temp_size += 1
            n += 1
            precs_temp[precs_temp_size] = prec_idx

            nce = getNCE(lookup, prec_charge, prec_mz)
        
            
            transition_idx = @inline fillTransitionList!(
                transitions,
                prec_estimation_type,
                getPrecFragRange(lookup, prec_idx),
                getFragments(lookup),
                getSplineData(lookup, prec_charge, prec_mz),
                prec_mz,
                prec_charge,
                prec_sulfur_count,
                transition_idx,
                quad_transmission_func,
                precursor_transmission,
                isotopes,
                n_frag_isotopes,
                max_frag_rank,
                iso_splines,
                frag_mz_bounds,
                block_size
            )
        end
    end

    return transition_idx, precs_temp_size
end


# Helper functions to improve readability and maintainability
function getQuadrupoleBounds(quad_func::QuadTransmissionFunction)
    return getPrecMinBound(quad_func), getPrecMaxBound(quad_func)
end

function getPrecursorWindowRange(
    precs, min_prec_mz::Float32, max_prec_mz::Float32, 
    isotope_err_bounds::Tuple{I, I}
) where {I<:Integer}
    start = searchsortedfirst(precs, by = x->last(x), 
        min_prec_mz - first(isotope_err_bounds)*NEUTRON/2)
    stop = searchsortedlast(precs, by = x->last(x), 
        max_prec_mz + last(isotope_err_bounds)*NEUTRON/2)
    return start, stop
end

function withinQuadrupoleBounds(
    prec_mz::Float32, prec_charge::UInt8,
    min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{I, I}
)where {I<:Integer}
    mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
    mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge
    return mz_low ≤ prec_mz ≤ mz_high
end

