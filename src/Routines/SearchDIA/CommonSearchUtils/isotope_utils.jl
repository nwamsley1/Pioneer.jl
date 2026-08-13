# Isotope capture and transmission utilities.
#
# DataFrame-level wrappers around isotopeSplines.jl functions
# (getPrecursorIsotopeSet, getPrecursorFractionTransmitted!).
# Shared by MainSearch and IntegrateChromatogramsSearch.

#############################################################################
# Per-row isotope computation (shared inner loop)
#############################################################################

"""
    _compute_scan_window(scan_id, centerMz, isolationWidthMz)

Extract scan center m/z and half-width for a given scan, handling missing values.
Returns (low_mz, high_mz, scan_mz, window_width).
"""
@inline function _compute_scan_window(scan_id::UInt32,
                                       centerMz::AbstractVector{Union{Missing, Float32}},
                                       isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    scan_mz = coalesce(centerMz[scan_id], zero(Float32))::Float32
    window_width = coalesce(isolationWidthMz[scan_id], zero(Float32))::Float32
    low_mz = Float32(scan_mz - window_width / 2)
    high_mz = Float32(scan_mz + window_width / 2)
    return low_mz, high_mz, scan_mz, window_width
end

"""
    _compute_fraction_transmitted(iso_splines, quad_transmission_model,
                                  mz, charge, sulfur, scan_mz, window_width,
                                  precursor_transmission)

Compute the fraction of precursor signal transmitted through the quadrupole.
Uses the pre-allocated `precursor_transmission` buffer to avoid allocation.
"""
@inline function _compute_fraction_transmitted(
    iso_splines::IsotopeSplineModel{Float32},
    quad_transmission_model::QuadTransmissionModel,
    mz::Float32, charge::UInt8, sulfur::UInt8,
    scan_mz::Float32, window_width::Float32,
    precursor_transmission::Vector{Float32})

    return getPrecursorFractionTransmitted!(
        precursor_transmission,
        iso_splines,
        (1, 5),
        getQuadTransmissionFunction(quad_transmission_model, scan_mz, window_width),
        mz, charge, sulfur)
end

#############################################################################
# Public API
#############################################################################

"""
    get_isotopes_captured!(chroms, quad_transmission_model, search_data,
                           scan_idx, prec_charge, prec_mz, sulfur_count,
                           centerMz, isolationWidthMz;
                           compute_isotope_set=true)

Compute precursor fraction transmitted and optionally which isotopes are
captured in each isolation window. The caller may pass either a PSM table or a
chromatogram-row table.

Adds columns to `chroms`:
- `:precursor_fraction_transmitted` (always)
- `:isotopes_captured` (when `compute_isotope_set=true`)

Chromatogram integration always needs `:precursor_fraction_transmitted` for
transmission-corrected smoothing. It only needs row-level `:isotopes_captured`
when integrating separate isotope traces.

Multi-threaded: partitions rows into chunks and spawns tasks.
"""
function get_isotopes_captured!(chroms::DataFrame,
                                quad_transmission_model::QuadTransmissionModel,
                                search_data::Vector{SimpleLibrarySearch{IsotopeSplineModel{Float32}}},
                                scan_idx::AbstractVector{UInt32},
                                prec_charge::AbstractArray{UInt8},
                                prec_mz::AbstractArray{Float32},
                                sulfur_count::AbstractArray{UInt8},
                                centerMz::AbstractVector{Union{Missing, Float32}},
                                isolationWidthMz::AbstractVector{Union{Missing, Float32}};
                                compute_isotope_set::Bool = true)
    n = size(chroms, 1)
    precursor_fraction_transmitted = Vector{Float32}(undef, n)
    isotopes_captured = compute_isotope_set ? Vector{Tuple{Int8, Int8}}(undef, n) : nothing

    # Cache column vector (avoid DataFrame row-indexing in hot loop)
    prec_idx_col = chroms[!, :precursor_idx]

    parallel_foreach!(n) do chunk
        thread_id = (first(chunk) % Threads.nthreads()) + 1
        iso_splines = getIsoSplines(search_data[thread_id])
        precursor_transmission = zeros(Float32, 5)

        for i in chunk
            prec_id = prec_idx_col[i]
            mz = prec_mz[prec_id]
            charge = prec_charge[prec_id]
            sulfur = sulfur_count[prec_id]
            scan_id = scan_idx[i]

            low_mz, high_mz, scan_mz, window_width = _compute_scan_window(
                scan_id, centerMz, isolationWidthMz)

            precursor_fraction_transmitted[i] = _compute_fraction_transmitted(
                iso_splines, quad_transmission_model,
                mz, charge, sulfur, scan_mz, window_width,
                precursor_transmission)

            if compute_isotope_set
                isotopes_captured[i] = getPrecursorIsotopeSet(mz, charge, low_mz, high_mz)
            end
        end
    end

    chroms[!, :precursor_fraction_transmitted] = precursor_fraction_transmitted
    if compute_isotope_set
        chroms[!, :isotopes_captured] = isotopes_captured
    end
    return nothing
end
