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
    transition_idx::Int64,
    library_fragment_lookup::LibraryFragmentLookup,
    precs_temp::Vector{UInt32},
    precs_temp_size::Int64,
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
    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
    block_size::Int64 = 10000
)
    n = 0
    min_prec_mz, max_prec_mz = get_quadrupole_bounds(quad_transmission_func)
    
    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec
        
        # Get precursor range within isolation window
        start, stop = get_precursor_window_range(
            precs, min_prec_mz, max_prec_mz, isotope_err_bounds
        )

        for i in start:stop
            prec_idx = first(precs[i])
            
            # Optional precursor filtering
            if !isnothing(precursors_passing) && prec_idx ∉ precursors_passing
                continue
            end

            # Get and validate precursor properties
            prec_props = get_precursor_properties(
                prec_idx, prec_mzs, prec_charges, prec_sulfur_counts
            )
            
            # Check quadrupole bounds
            if !within_quadrupole_bounds(
                prec_props.mz, prec_props.charge, 
                min_prec_mz, max_prec_mz, isotope_err_bounds
            )
                continue
            end

            # Update counters and temp storage
            precs_temp_size += 1
            n += 1
            precs_temp[precs_temp_size] = prec_idx

            # Fill transitions for this precursor
            transition_idx = @inline fillTransitionList!(
                transitions, library_fragment_lookup,
                prec_idx, prec_props.mz, prec_props.charge, prec_props.sulfur_count,
                transition_idx, quad_transmission_func,
                precursor_transmission, isotopes,
                n_frag_isotopes, max_frag_rank, false,
                iso_splines, frag_mz_bounds, block_size
            )
        end
    end

    return transition_idx, n, precs_temp_size
end


# Helper functions to improve readability and maintainability
function get_quadrupole_bounds(quad_func::QuadTransmissionFunction)
    return getPrecMinBound(quad_func), getPrecMaxBound(quad_func)
end

function get_precursor_window_range(
    precs, min_prec_mz::Float32, max_prec_mz::Float32, 
    isotope_err_bounds::Tuple{Int64, Int64}
)
    start = searchsortedfirst(precs, by = x->last(x), 
        min_prec_mz - first(isotope_err_bounds)*NEUTRON/2)
    stop = searchsortedlast(precs, by = x->last(x), 
        max_prec_mz + last(isotope_err_bounds)*NEUTRON/2)
    return start, stop
end

function get_precursor_properties(
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

function within_quadrupole_bounds(
    prec_mz::Float32, prec_charge::UInt8,
    min_prec_mz::Float32, max_prec_mz::Float32,
    isotope_err_bounds::Tuple{Int64, Int64}
)
    mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
    mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge
    return mz_low ≤ prec_mz ≤ mz_high
end
