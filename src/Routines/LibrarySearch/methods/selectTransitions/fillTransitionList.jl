"""
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


"""
    fillTransitionList!(transitions, prec_estimation_type, args...) -> Int64

Fill a transition list with fragment ions and their isotopes based on the specified
precursor estimation strategy.

# Arguments
- `transitions::Vector{DetailedFrag{Float32}}`: Vector to store generated transitions
- `prec_estimation_type::PrecEstimation`: Strategy for isotope pattern calculation
- `precursor_fragment_range::UnitRange{UInt64}`: Range of fragments to process
- `fragment_ions::Vector{AltimeterFragment}`: The library fragment ions
- `nce::Union{Missing, Float32}`: Normalized collision energy (if applicable)
- `knots::Union{Missing, NTuple{M, Float32}}`: Spline knots for intensity prediction
- `prec_mz::Float32`: Precursor m/z
- `prec_charge::UInt8`: Precursor charge
- `prec_sulfur_count::UInt8`: Number of sulfur atoms in precursor
- `transition_idx::Int64`: Current transition index
- `quad_transmission_func::QuadTransmissionFunction`: Quadrupole transmission function
- `precursor_transmission::Vector{Float32}`: Buffer for precursor transmission values
- `isotopes::Vector{Float32}`: Buffer for isotope calculations
- `n_frag_isotopes::Int64`: Number of fragment isotopes to consider
- `max_frag_rank::UInt8`: Maximum fragment rank to include
- `iso_splines::IsotopeSplineModel`: Model for isotope pattern prediction
- `frag_mz_bounds::Tuple{Float32, Float32}`: m/z bounds for fragments
- `block_size::Int64`: Size for array growth when needed

# Returns
- `Int64`: Updated transition index
"""
function fillTransitionList!(transitions::Vector{DetailedFrag{Float32}}, 
                            prec_estimation_type::PrecEstimation,
                            precursor_fragment_range::UnitRange{UInt64},
                            fragment_ions::Vector{F},
                            nce::Union{Missing, Float32},
                            knots::Union{Missing, NTuple{M, Float32}},
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            transition_idx::Int64, 
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32}, 
                            n_frag_isotopes::Int64,
                            max_frag_rank::UInt8,
                            iso_splines::IsotopeSplineModel, 
                            frag_mz_bounds::Tuple{Float32, Float32},
                            block_size::Int64)::Int64 where {M, F <: AltimeterFragment}#where {T,U,V,W<:AbstractFloat,I<:Integer}

    # Calculate precursor isotope transmission and range
    getPrecursorIsotopeTransmission!(precursor_transmission, prec_mz, prec_charge, quad_transmission_func)
    prec_isotope_set = getPrecursorIsotopeSet(prec_mz, prec_charge, quad_transmission_func)
    frag_iso_idx_range = range(0, min(n_frag_isotopes - 1, last(prec_isotope_set)))

    for frag_idx in precursor_fragment_range
        frag = fragment_ions[frag_idx]
        frag.rank > max_frag_rank && continue
        # Calculate isotope pattern based on estimation strategy
        getFragIsotopes!(prec_estimation_type, isotopes, precursor_transmission,
        prec_isotope_set, frag_iso_idx_range, iso_splines,
        prec_mz, prec_charge, prec_sulfur_count, frag, knots, nce)
        # Create transitions for each isotope
        transition_idx = addTransitionIsotopes!(transitions, transition_idx, 
                                                frag, isotopes, frag_iso_idx_range,
                                                frag_mz_bounds, block_size)
    end
    return transition_idx
end


"""
Helper function to add transitions for each isotope of a fragment.
"""
function addTransitionIsotopes!(transitions::Vector{DetailedFrag{Float32}},
                                transition_idx::Int64,
                                frag::AltimeterFragment,
                                isotopes::Vector{Float32},
                                frag_iso_idx_range::UnitRange{Int64},
                                frag_mz_bounds::Tuple{Float32, Float32},
                                block_size::Int64)::Int64
    for iso_idx in frag_iso_idx_range
        frag_mz = Float32(frag.mz + iso_idx * NEUTRON/frag.frag_charge)
        
        # Skip if outside m/z bounds
        (frag_mz < first(frag_mz_bounds) || frag_mz > last(frag_mz_bounds)) && continue

        transition_idx += 1
        transitions[transition_idx] = DetailedFrag(
            frag.prec_id, frag_mz, Float16(isotopes[iso_idx + 1]),
            frag.ion_type, frag.is_y, frag.is_b, frag.is_p, iso_idx > 0,
            frag.frag_charge, frag.ion_position, frag.prec_charge,
            frag.rank, frag.sulfur_count
        )

        ensureTransitionCapacity!(transitions, transition_idx, block_size)
    end
    return transition_idx
end

function getFragIsotopes!(
                            ::PartialPrecCapture,
                            frag_isotopes::Vector{Float32}, 
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::SplineDetailedFrag{N, Float32}, 
                            knots::NTuple{M, Float32},
                            nce::Float32) where {M, N}
    #Reset relative abundances of isotopes to zero 
    fill!(frag_isotopes, zero(eltype(frag_isotopes)))
    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity =  getIntensity(frag, knots, 3, nce)

    getFragAbundance!(
                    frag_isotopes, 
                    precursor_transmition,
                    iso_splines,  
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count, 
                    frag
                    )

    #Estimate abundances of M+n fragment ions relative to the monoisotope
    for i in reverse(range(1, length(frag_isotopes)))
        frag_isotopes[i] = total_fragment_intensity*frag_isotopes[i]
    end
end

function getFragIsotopes!(
                            ::FullPrecCapture,
                            frag_isotopes::Vector{Float32}, 
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel, 
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::SplineDetailedFrag{N, Float32}, 
                            knots::NTuple{M, Float32},
                            nce::Float32) where {M, N}

    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity = getIntensity(frag, knots, 3, nce)
    frag_mz = getMz(frag)
    frag_charge = getCharge(frag)
    frag_nsulfur = Int64(getSulfurCount(frag))
    for i in frag_iso_idx_range
        frag_isotopes[i+1] = iso_splines(
                                min(frag_nsulfur, 5), 
                                iso_idx, 
                                frag_mz*frag_charge
                                )*total_fragment_intensity
    end
    return nothing
end
