# Per-fragment isotope abundance prediction.
#
# Used by run_fused! (fusedMatch.jl) to predict the relative intensity
# of each isotope of a library fragment under the current precursor
# transmission profile. Originally lived in
# selectTransitions/fillTransitionList.jl alongside the deleted
# selectTransitions! dispatcher; survived because the fused path still
# needs it.
#
# Two methods, dispatching on PrecEstimation:
#   - PartialPrecCapture: Goldfarb-style detailed abundance calculation
#     (uses getFragAbundance!) — accounts for partial isotopologue capture
#     when only part of the precursor envelope is in the isolation window.
#   - FullPrecCapture: simplified direct iso_splines evaluation —
#     assumes complete capture.

function getFragIsotopes!(
        ::PartialPrecCapture,
        frag_isotopes::Vector{Float32},
        precursor_transmission::Vector{Float32},
        prec_isotope_set::Tuple{Int64, Int64},
        frag_iso_idx_range::UnitRange{Int64},
        iso_splines::IsotopeSplineModel,
        prec_mz::Float32,
        prec_charge::UInt8,
        prec_sulfur_count::UInt8,
        frag::F,
        spline_data::G) where {F<:AltimeterFragment, G<:IntensityDataType}

    fill!(frag_isotopes, zero(eltype(frag_isotopes)))
    total_fragment_intensity = getIntensity(frag, spline_data)

    getFragAbundance!(
        frag_isotopes,
        precursor_transmission,
        iso_splines,
        prec_mz,
        prec_charge,
        prec_sulfur_count,
        frag)

    for i in reverse(range(1, length(frag_isotopes)))
        frag_isotopes[i] = total_fragment_intensity * frag_isotopes[i]
    end
end

function getFragIsotopes!(
        ::FullPrecCapture,
        frag_isotopes::Vector{Float32},
        precursor_transmission::Vector{Float32},
        prec_isotope_set::Tuple{Int64, Int64},
        frag_iso_idx_range::UnitRange{Int64},
        iso_splines::IsotopeSplineModel,
        prec_mz::Float32,
        prec_charge::UInt8,
        prec_sulfur_count::UInt8,
        frag::F,
        spline_data::G) where {F<:AltimeterFragment, G<:IntensityDataType}

    total_fragment_intensity = getIntensity(frag, spline_data)
    frag_mz = getMz(frag)
    frag_charge = getPrecCharge(frag)
    frag_nsulfur = Int64(getSulfurCount(frag))
    for iso_idx in frag_iso_idx_range
        frag_isotopes[iso_idx + 1] = iso_splines(
            min(frag_nsulfur, 5),
            iso_idx,
            frag_mz * frag_charge) * total_fragment_intensity
    end
    return nothing
end
