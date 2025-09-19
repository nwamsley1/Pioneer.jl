# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    FragmentMassError

Structure to hold individual fragment mass error data for analysis.

Contains all the information needed to investigate mass error distributions
relative to fragment intensity, retention time, and other factors.
"""
struct FragmentMassError{T<:AbstractFloat}
    # Mass error information
    ppm_error::T                    # PPM mass error: (observed - theoretical) / (theoretical / 1e6)
    theoretical_mz::T               # Theoretical m/z of the fragment
    observed_mz::T                  # Observed m/z of the fragment

    # Intensity information
    fragment_intensity::T           # Intensity of this specific fragment
    library_intensity::T            # Library/predicted intensity of this fragment
    max_fragment_intensity::T       # Intensity of the most intense fragment for this precursor in this scan

    # Context information
    precursor_q_value::T            # Q-value of the precursor PSM this fragment belongs to
    retention_time::T               # Retention time of the scan
    scan_idx::UInt32                # Scan index
    ms_file_idx::UInt32             # MS file index
    precursor_idx::UInt32           # Precursor index

    # Fragment metadata
    fragment_charge::UInt8          # Fragment charge
    ion_type::UInt8                 # Ion type (b, y, p)
    fragment_number::UInt8          # Fragment position/number
    is_isotope::Bool               # Whether this is an isotope peak
end

"""
    FragmentMassError{Float32}()

Default constructor for FragmentMassError.
"""
FragmentMassError{Float32}() = FragmentMassError{Float32}(
    zero(Float32),  # ppm_error
    zero(Float32),  # theoretical_mz
    zero(Float32),  # observed_mz
    zero(Float32),  # fragment_intensity
    zero(Float32),  # library_intensity
    zero(Float32),  # max_fragment_intensity
    zero(Float32),  # precursor_q_value
    zero(Float32),  # retention_time
    zero(UInt32),   # scan_idx
    zero(UInt32),   # ms_file_idx
    zero(UInt32),   # precursor_idx
    zero(UInt8),    # fragment_charge
    zero(UInt8),    # ion_type
    zero(UInt8),    # fragment_number
    false           # is_isotope
)

"""
    collect_fragment_mass_errors!(fragment_errors::Vector{FragmentMassError{Float32}},
                                  fragment_matches::Vector{FragmentMatch{Float32}},
                                  n_matches::Int,
                                  psm_q_value::Float32,
                                  retention_time::Float32,
                                  scan_idx::UInt32,
                                  ms_file_idx::UInt32,
                                  precursor_idx::UInt32,
                                  current_idx::Int) -> Int

Collects fragment mass errors from fragment matches into the collection vector.

# Arguments
- `fragment_errors`: Vector to store fragment mass error data
- `fragment_matches`: Vector of fragment matches from peak matching
- `n_matches`: Number of valid matches in fragment_matches
- `psm_q_value`: Q-value of the PSM these fragments belong to
- `retention_time`: Retention time of the scan
- `scan_idx`: Scan index
- `ms_file_idx`: MS file index
- `precursor_idx`: Precursor index
- `current_idx`: Current index in fragment_errors vector

# Returns
- Updated index after adding fragment errors

# Notes
- Only collects data from PSMs that pass 1% FDR threshold (q_value <= 0.01)
- Calculates PPM mass error for each fragment
- Finds the maximum fragment intensity for normalization
- Grows the fragment_errors vector as needed
"""
function collect_fragment_mass_errors!(
    fragment_errors::Vector{FragmentMassError{Float32}},
    fragment_matches::Vector{FragmentMatch{Float32}},
    n_matches::Int,
    psm_q_value::Float32,
    retention_time::Float32,
    scan_idx::UInt32,
    ms_file_idx::UInt32,
    precursor_idx::UInt32,
    current_idx::Int
)::Int
    # Only collect from PSMs passing 1% FDR
    if psm_q_value > 0.01f0 || n_matches == 0
        return current_idx
    end

    # Find maximum fragment intensity for this precursor/scan
    max_intensity = zero(Float32)
    for i in 1:n_matches
        fragment = fragment_matches[i]
        if getIntensity(fragment) > max_intensity
            max_intensity = getIntensity(fragment)
        end
    end

    # Collect each fragment mass error
    idx = current_idx
    for i in 1:n_matches
        fragment = fragment_matches[i]

        # Calculate PPM mass error
        theoretical_mz = getMZ(fragment)
        observed_mz = getMatchMZ(fragment)
        ppm_error = (observed_mz - theoretical_mz) / (theoretical_mz / 1e6f0)

        # Grow vector if needed
        if idx + 1 > length(fragment_errors)
            resize!(fragment_errors, length(fragment_errors) + 10000)
        end

        idx += 1
        fragment_errors[idx] = FragmentMassError{Float32}(
            ppm_error,
            theoretical_mz,
            observed_mz,
            getIntensity(fragment),
            getPredictedIntensity(fragment),
            max_intensity,
            psm_q_value,
            retention_time,
            scan_idx,
            ms_file_idx,
            precursor_idx,
            getCharge(fragment),
            getIonType(fragment),
            getFragInd(fragment),
            isIsotope(fragment)
        )
    end

    return idx
end

"""
    write_fragment_mass_errors(fragment_errors::Vector{FragmentMassError{Float32}},
                               n_errors::Int,
                               output_dir::String,
                               ms_file_name::String)

Writes collected fragment mass errors to an Arrow file.

# Arguments
- `fragment_errors`: Vector containing fragment mass error data
- `n_errors`: Number of valid entries in fragment_errors
- `output_dir`: Output directory path
- `ms_file_name`: Name of the MS file for naming the output

# Output File Format
- Name: `fragment_mass_errors_{ms_file_name}_{timestamp}.arrow`
- Contains columns for all FragmentMassError fields
- Suitable for analysis with DataFrames.jl and Arrow.jl
"""
function write_fragment_mass_errors(
    fragment_errors::Vector{FragmentMassError{Float32}},
    n_errors::Int,
    output_dir::String,
    ms_file_name::String
)
    if n_errors == 0
        @user_warn "No fragment mass errors collected for file $ms_file_name"
        return
    end

    # Create output filename with timestamp
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_filename = "fragment_mass_errors_$(ms_file_name)_$(timestamp).arrow"
    output_path = joinpath(output_dir, output_filename)

    # Convert to DataFrame format
    df = DataFrame(
        ppm_error = [fragment_errors[i].ppm_error for i in 1:n_errors],
        theoretical_mz = [fragment_errors[i].theoretical_mz for i in 1:n_errors],
        observed_mz = [fragment_errors[i].observed_mz for i in 1:n_errors],
        fragment_intensity = [fragment_errors[i].fragment_intensity for i in 1:n_errors],
        library_intensity = [fragment_errors[i].library_intensity for i in 1:n_errors],
        max_fragment_intensity = [fragment_errors[i].max_fragment_intensity for i in 1:n_errors],
        precursor_q_value = [fragment_errors[i].precursor_q_value for i in 1:n_errors],
        retention_time = [fragment_errors[i].retention_time for i in 1:n_errors],
        scan_idx = [fragment_errors[i].scan_idx for i in 1:n_errors],
        ms_file_idx = [fragment_errors[i].ms_file_idx for i in 1:n_errors],
        precursor_idx = [fragment_errors[i].precursor_idx for i in 1:n_errors],
        fragment_charge = [fragment_errors[i].fragment_charge for i in 1:n_errors],
        ion_type = [fragment_errors[i].ion_type for i in 1:n_errors],
        fragment_number = [fragment_errors[i].fragment_number for i in 1:n_errors],
        is_isotope = [fragment_errors[i].is_isotope for i in 1:n_errors]
    )

    # Write to Arrow file
    Arrow.write(output_path, df)

    @user_info "Fragment mass errors written to: $output_path"
    @user_info "Collected $(n_errors) fragment mass errors from $(length(unique(df.precursor_idx))) precursors"

    return output_path
end

"""
    collect_fragment_errors_from_psms!(results::FirstPassSearchResults,
                                       spectra::MassSpecData,
                                       search_context::SearchContext,
                                       params::FirstPassSearchParameters,
                                       ms_file_idx::Int64)

Collects fragment mass errors from PSMs that pass 1% FDR threshold.

This function re-performs peak matching for all PSMs with q-value <= 0.01
to collect detailed fragment mass error information for analysis.

# Arguments
- `results`: FirstPassSearchResults to store fragment errors
- `spectra`: Mass spectrometry data
- `search_context`: Current search context
- `params`: Search parameters
- `ms_file_idx`: Index of the MS file being processed

# Process
1. Filters PSMs to those passing 1% FDR (q_value <= 0.01)
2. For each passing PSM, re-performs fragment matching
3. Collects mass error data for each matched fragment
4. Stores data in results.fragment_errors with proper indexing

# Notes
- Only collects from target PSMs (excludes decoys)
- Re-performs peak matching to get detailed fragment match data
- Calculates retention time and intensity information
"""
function collect_fragment_errors_from_psms!(
    results::FirstPassSearchResults,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    ms_file_idx::Int64
)
    psms = results.psms[]

    # Filter to PSMs passing 1% FDR and targets only
    passing_psms = psms[(psms.q_value .<= 0.01f0) .& psms.target, :]

    if nrow(passing_psms) == 0
        @user_info "No PSMs passed 1% FDR for fragment error collection in file $ms_file_idx"
        return
    end

    @user_info "Collecting fragment mass errors from $(nrow(passing_psms)) PSMs passing 1% FDR"

    # Get necessary data structures
    spec_lib = getSpecLib(search_context)
    search_data = getSearchData(search_context)
    fragment_index = getFragmentIndex(spec_lib)
    precursors = getPrecursors(spec_lib)
    ion_list = getFragmentLookupTable(spec_lib)
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_model = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = getIrtErrors(search_context)[ms_file_idx]

    current_error_idx = results.fragment_error_count[]

    # Group PSMs by scan for efficient processing
    psm_groups = groupby(passing_psms, :scan_idx)

    for group in psm_groups
        scan_idx = first(group.scan_idx)

        # Skip invalid scans
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        getMsOrder(spectra, scan_idx) ∉ getSpecOrder(params) && continue

        # Get precursor IDs for this scan
        prec_ids = UInt32.(group.precursor_idx)

        # Select transitions for these precursors
        isotopes = zeros(Float32, 5)
        precursor_transmission = zeros(Float32, 5)

        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data[1]),
            StandardTransitionSelection(),
            getPrecEstimation(params),
            ion_list,
            1:length(prec_ids), prec_ids,  # Mock scan_to_prec_idx range
            getMz(precursors),
            getCharge(precursors),
            getSulfurCount(precursors),
            getIrt(precursors),
            getIsoSplines(search_data[1]),
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            precursor_transmission, isotopes, getNFragIsotopes(params),
            getMaxFragRank(params),
            Float32(getModel(rt_model)(getRetentionTime(spectra, scan_idx))),
            Float32(irt_tol),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        if ion_idx < 2
            continue
        end

        # Match peaks for this scan
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data[1]),
            getIonMisses(search_data[1]),
            getIonTemplates(search_data[1]),
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )

        if nmatches == 0
            continue
        end

        # For each PSM in this scan, collect fragment errors
        for psm_row in eachrow(group)
            precursor_idx = UInt32(psm_row.precursor_idx)
            q_value = Float32(psm_row.q_value)
            retention_time = Float32(psm_row.rt)

            # Filter fragment matches to this precursor
            precursor_matches = Vector{FragmentMatch{Float32}}()
            for i in 1:nmatches
                fragment = getIonMatches(search_data[1])[i]
                if getPrecID(fragment) == precursor_idx
                    push!(precursor_matches, fragment)
                end
            end

            if !isempty(precursor_matches)
                current_error_idx = collect_fragment_mass_errors!(
                    results.fragment_errors,
                    precursor_matches,
                    length(precursor_matches),
                    q_value,
                    retention_time,
                    UInt32(scan_idx),
                    UInt32(ms_file_idx),
                    precursor_idx,
                    current_error_idx
                )
            end
        end
    end

    # Update the count
    results.fragment_error_count[] = current_error_idx

    @user_info "Collected $(current_error_idx - results.fragment_error_count[]) fragment mass errors from file $ms_file_idx"
end