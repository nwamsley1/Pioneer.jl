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
    FirstPassSearch

Initial search method to identify PSMs and establish retention time calibration.

This search:
1. Performs initial PSM identification with learned scoring
2. Calculates retention time indices and FWHM statistics
3. Maps between library and empirical retention times
4. Generates RT calibration curves for subsequent searches

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (1, 0),
    :first_search_params => Dict(
        "n_train_rounds_probit" => 10,
        "max_iter_probit" => 100,
        "max_q_value_probit_rescore" => 0.01,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "max_precursors_passing" => 5000,
        "abreviate_precursor_calc" => false
    ),
    :summarize_first_search_params => Dict(
        "min_inference_points" => 1000,
        "max_q_val_for_irt" => 0.01,
        "max_precursors" => 10000,
        "max_irt_bin_size" => 0.1,
        "max_prob_to_impute" => 0.99
    ),
    :irt_mapping_params => Dict(
        "min_prob" => 0.9
    )
)

# Execute search
results = execute_search(FirstPassSearch(), search_context, params)
```
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

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
Results container for first pass search.
Holds FWHM statistics, PSM file paths, and model updates.
"""
struct FirstPassSearchResults <: SearchResults
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32,mad_fwhm::Float32}}
    psms::Base.Ref{DataFrame}
    ms1_mass_err_model::Base.Ref{<:MassErrorModel}
    ms1_ppm_errs::Vector{Float32}
    ms1_mass_plots::Vector{Plots.Plot}
    qc_plots_folder_path::String
    # Fragment mass error collection
    fragment_errors::Vector{FragmentMassError{Float32}}
    fragment_error_count::Base.Ref{Int}
end

"""
Parameters for first pass search.
Configures PSM identification, scoring, and RT calibration.
"""
struct FirstPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    frag_tol_ppm::Float32
    ms1_tol_ppm::Float32
    frag_err_quantile::Float32
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    match_between_runs::Bool
    relative_improvement_threshold::Float32
    
    # Scoring parameters
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
    max_PEP::Float32
    
    # RT parameters
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    prec_estimation::P

    function FirstPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        first_params = params.first_search
        frag_params = first_params.fragment_settings
        score_params = first_params.scoring_settings
        rt_params = params.rt_alignment
        irt_mapping_params = first_params.irt_mapping
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        # Determine precursor estimation strategy
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            0.0f0,  # No transmission threshold for first pass
            0.0f0,  # No fragment tolerance for first pass
            Float32(params.parameter_tuning.iteration_settings.ms1_tol_ppm),  # MS1 tolerance from config
            Float32(params.parameter_tuning.search_settings.frag_err_quantile),
            # Handle min_score as either single value or array (use first value if array)
            begin
                min_score_raw = frag_params.min_score
                if min_score_raw isa Vector
                    UInt8(first(min_score_raw))
                else
                    UInt8(min_score_raw)
                end
            end,
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1), # max_best_rank
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            global_params.match_between_runs,
            Float32(frag_params.relative_improvement_threshold),
            
            Int64(score_params.n_train_rounds),
            Int64(score_params.max_iterations),
            Float32(score_params.max_q_value_probit_rescore),
            Float32(score_params.max_PEP),
            
            Int64(1000), # Default min_inference_points
            Float32(rt_params.min_probability),
            Float32(rt_params.min_probability),
            Float32(0.1), # Default max_irt_bin_size
            Float32(irt_mapping_params.max_prob_to_impute_irt),  # Default max_prob_to_impute
            Float32(irt_mapping_params.fwhm_nstd),   # Default fwhm_nstd
            Float32(irt_mapping_params.irt_nstd),   # Default irt_nstd
            prec_estimation
        )
    end
end

#==========================================================
Fragment Mass Error Collection Functions
==========================================================#

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
"""
function write_fragment_mass_errors(
    fragment_errors::Vector{FragmentMassError{Float32}},
    n_errors::Int,
    output_dir::String,
    ms_file_name::String
)
    if n_errors == 0
        @warn "No fragment mass errors collected for file $ms_file_name"
        return
    end

    # Create output filename with timestamp
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    output_filename = "fragment_mass_errors_$(ms_file_name)_$(timestamp).arrow"
    output_path = joinpath(output_dir, output_filename)

    # Convert to DataFrame format
    df = DataFrames.DataFrame(
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

    @info "Fragment mass errors written to: $output_path"
    @info "Collected $(n_errors) fragment mass errors from $(length(unique(df.precursor_idx))) precursors"

    return output_path
end

"""
    collect_fragment_errors_from_psms!(results::FirstPassSearchResults,
                                       spectra::MassSpecData,
                                       search_context::SearchContext,
                                       params::FirstPassSearchParameters,
                                       ms_file_idx::Int64)

Collects fragment mass errors from PSMs that pass 1% FDR threshold.
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

    if DataFrames.nrow(passing_psms) == 0
        @info "No PSMs passed 1% FDR for fragment error collection in file $ms_file_idx"
        return
    end

    @info "Collecting fragment mass errors from $(DataFrames.nrow(passing_psms)) PSMs passing 1% FDR"

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
    psm_groups = DataFrames.groupby(passing_psms, :scan_idx)

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
        for psm_row in DataFrames.eachrow(group)
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
    previous_count = results.fragment_error_count[]
    results.fragment_error_count[] = current_error_idx

    @info "Collected $(current_error_idx - previous_count) fragment mass errors from file $ms_file_idx"
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = FirstPassSearchParameters(params)
getMs1MassErrorModel(ptsr::FirstPassSearchResults) = ptsr.ms1_mass_err_model[]
getMs1TolPpm(params::FirstPassSearchParameters) = params.ms1_tol_ppm

function init_search_results(
    ::FirstPassSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
    !isdir(temp_folder) && mkdir(temp_folder)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    ms1_mass_error_plots = joinpath(qc_dir, "ms1_mass_error_plots")
    !isdir(ms1_mass_error_plots ) && mkdir(ms1_mass_error_plots )
    return FirstPassSearchResults(
        Dictionary{Int64, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}(),
        Base.Ref{DataFrame}(),
        Base.Ref{MassErrorModel}(),
        Vector{Float32}(),
        Plots.Plot[],
        qc_dir,
        # Initialize fragment error collection with pre-allocated space
        Vector{FragmentMassError{Float32}}(undef, 100000),
        Ref(0)
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the first pass search.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:FirstPassSearchParameters}

    """
    Perform library search with current parameters.
    """
    function perform_library_search(
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        
        return library_search(spectra, search_context, params, ms_file_idx)
    end

    """
    Process PSMs from library search.
    """
    function process_psms!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)

        """
        Select best PSMs based on criteria.
        """
        function select_best_psms!(
            psms::DataFrame,
            precursor_mzs::AbstractVector,
            params::FirstPassSearchParameters,
            search_context::SearchContext
        )
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
            get_best_psms!(
                psms,
                precursor_mzs,
                max_PEP=params.max_PEP,
                fdr_scale_factor=fdr_scale_factor
            )
        end

        rt_model = getRtIrtModel(search_context, ms_file_idx)
        # Add columns
        add_psm_columns!(psms, spectra, search_context, rt_model, ms_file_idx)
        
        # Score PSMs
        score_psms!(psms, params, search_context)
        # Get best PSMs
        select_best_psms!(
            psms,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            params,
            search_context
        )
        return psms
    end

    """
    Add necessary columns to PSM DataFrame.
    """
    function add_psm_columns!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        rt_model::RtConversionModel,
        ms_file_idx::Int64)
        add_main_search_columns!(
            psms,
            getModel(rt_model),
            getStructuralMods(getPrecursors(getSpecLib(search_context))),
            getMissedCleavages(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getIrt(getPrecursors(getSpecLib(search_context))),
            getCharge(getPrecursors(getSpecLib(search_context))),
            getRetentionTimes(spectra),
            getTICs(spectra),
            getMzArrays(spectra)
        )
        # Calculate RT values
        psms[!, :irt_observed] = rt_model.(psms[!, :rt])
        psms[!, :irt_error] = Float16.(abs.(psms[!, :irt_observed] .- psms[!, :irt_predicted]))
        psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)
        psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    end

    """
    Score PSMs using probit model.
    """
    function score_psms!(
        psms::DataFrame,
        params::FirstPassSearchParameters,
        search_context::SearchContext)
        column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe, :percent_theoretical_ignored,
            :charge2, :poisson, :irt_error, 
            :missed_cleavage, 
            :Mox,
            #:charge, Only works with charge 2 if at least 3 charge states presence. otherwise singular error
            #:b_count, might be good for non-tryptic enzymes
            :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
        ]

        # Avoid singular error if no peaks were ignored
        if maximum(psms.percent_theoretical_ignored) == 0
            deleteat!(column_names, findfirst(==(:percent_theoretical_ignored), column_names))
        end


        # Select scoring columns
        select!(psms, vcat(column_names, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target]))
        # Score PSMs
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        try
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        catch
            column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe,
            :charge2, :poisson, :irt_error, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
            ]
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        end
        # Process scores
       
        select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target])
        get_probs!(psms, psms[!,:score])
        sort!(psms, :rt)
    end

    try
        # Get models and update fragment lookup table
        psms = perform_library_search(spectra, search_context, params, ms_file_idx)
        results.psms[] = process_psms!(psms, spectra, search_context, params, ms_file_idx)

        # Collect fragment mass errors from PSMs passing 1% FDR
        collect_fragment_errors_from_psms!(results, spectra, search_context, params, ms_file_idx)

        temp_psms = results.psms[] 
        temp_psms = temp_psms[temp_psms[!,:q_value].<=0.001,:]
        most_intense = sortperm(temp_psms[!,:log2_summed_intensity], rev = true)
        ms1_errs = vcat(
            mass_error_search(
                spectra,
                temp_psms[most_intense[1:(min(3000, length(most_intense)))],:scan_idx],
                temp_psms[most_intense[1:(min(3000, length(most_intense)))],:precursor_idx],
                UInt32(ms_file_idx),
                getSpecLib(search_context),
                getSearchData(search_context),
                MassErrorModel(
                0.0f0,
                (getMs1TolPpm(params), getMs1TolPpm(params))  # Use MS1 tolerance from JSON config
                ),
                params,
                MS1CHROM()
            )...
        )
        if length(ms1_errs) > 1
            mad_dev= mad(ms1_errs)
            med_errs = median(ms1_errs)
            low_bound, high_bound = med_errs - mad_dev*7, med_errs + mad_dev*7
            filter!(x->(low_bound<x)&(high_bound>x), ms1_errs)
            ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
            results.ms1_mass_err_model[] = ms1_mass_err_model
            append!(results.ms1_ppm_errs, ms1_ppm_errs)
            select!(results.psms[] , Not(:log2_summed_intensity))
        else
            #ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
            #Default to MS2 pattern 
            results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
            append!(results.ms1_ppm_errs, Float32[])
            select!(results.psms[] , Not(:log2_summed_intensity))
        end
    catch e
        # Get file name for debugging
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        
        reason = "FirstPassSearch failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        @user_warn "First pass search failed for MS data file: $file_name. Error type: $(typeof(e)). Creating empty results to continue pipeline."
        
        # Create an empty but properly structured DataFrame to avoid downstream errors
        empty_psms = DataFrame(
            ms_file_idx = UInt32[],
            scan_idx = UInt32[], 
            precursor_idx = UInt32[],
            rt = Float32[],
            irt_predicted = Float32[],
            q_value = Float32[],
            score = Float32[], 
            prob = Float32[],
            scan_count = UInt32[],
            fwhm = Float32[]  # Add this to prevent missing column error
        )
        results.psms[] = empty_psms
        
        # Set default mass error model
        results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
        #rethrow(e)
    end

    return results
end

"""
Initial file processing complete, no additional processing needed.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:FirstPassSearchParameters}
    psms = results.psms[]
    fwhms = skipmissing(psms[!, :fwhm])
    fwhm_points = count(!ismissing, fwhms)
    if fwhm_points >= 1#params.min_inference_points
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = median(fwhms),
            mad_fwhm = mad(fwhms, normalize=true)))
    else
        @user_warn "Insuficient fwhm_points to estimate for $ms_file_idx"
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = 0.2f0,
            mad_fwhm = 0.2f0))
    end
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    temp_path = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms", parsed_fname * ".arrow")
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
            :irt_predicted, :q_value, :score, :prob, :scan_count])
    )
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, temp_path)

    #####
    #MS1 mass tolerance 
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    # Generate mass error plot
    push!(results.ms1_mass_plots, generate_ms1_mass_error_plot(results, parsed_fname))
    # Update models in search context
    setMs1MassErrorModel!(search_context, ms_file_idx, getMs1MassErrorModel(results))

    # Write fragment mass errors to Arrow file
    if results.fragment_error_count[] > 0
        output_dir = getDataOutDir(search_context)
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        write_fragment_mass_errors(
            results.fragment_errors,
            results.fragment_error_count[],
            output_dir,
            parsed_fname
        )
    end
end

"""
No cleanup needed between files.
"""
function reset_results!(results::FirstPassSearchResults)
    empty!(results.psms[])
    resize!(results.ms1_ppm_errs, 0)
    results.fragment_error_count[] = 0
    return nothing
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}
    
    """
    Process precursors and calculate iRT errors.
    """
    function get_best_precursors_accross_runs!(
        search_context::SearchContext,
        results::FirstPassSearchResults,
        params::FirstPassSearchParameters
    )

        # Filter out failed files
        valid_indices = get_valid_file_indices(search_context)
        all_psms_paths = getFirstPassPsms(getMSData(search_context))
        valid_psms_paths = [all_psms_paths[i] for i in valid_indices]
        
        # Create RT-IRT map for valid files only
        all_rt_irt = getRtIrtModel(search_context)
        valid_rt_irt = Dict{Int64, RtConversionModel}(i => all_rt_irt[i] for i in valid_indices if haskey(all_rt_irt, i))
        
        if isempty(valid_psms_paths)
            @user_warn "No valid files for cross-run precursor analysis"
            return Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}}()
        end
        
        # Get best precursors from valid files only
        return get_best_precursors_accross_runs(
            valid_psms_paths,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            valid_rt_irt,
            max_q_val=params.max_q_val_for_irt
        )
    end
    # Map retention times
    map_retention_times!(search_context, results, params)
    # Process precursors
    precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)

    if false==true#params.match_between_runs==true
        #######
        #Each target has a corresponding decoy and vice versa
        #Add the complement targets/decoys to the precursor dict 
        #if the `sibling_peptide_scores` parameter is set to true
        #In the target/decoy scoring (see SearchMethods/ScoringSearch)
        #the maximum score for each target/decoy pair is shared accross runs
        #in an iterative training scheme. 
        precursors = getPrecursors(getSpecLib(search_context))
        i = 1
        for (pid, val) in pairs(precursor_dict)
            i += 1
            setPredIrt!(search_context, pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
            partner_pid = getPartnerPrecursorIdx(precursors)[pid]
            if ismissing(partner_pid)
                continue
            end

            # If the partner needs to be added, then give it the irt of the currently identified precursor
            # Otherwise if the partner was ID'ed, it should keep its original predicted iRT
            if !haskey(precursor_dict, partner_pid)
                insert!(precursor_dict, partner_pid, val)
                setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
            else
                setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[partner_pid])
            end
            
        end
    else
        for (pid, val) in pairs(precursor_dict)
            setPredIrt!(search_context, pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
        end
    end

    setPrecursorDict!(search_context, precursor_dict)
    # Calculate RT indices
    create_rt_indices!(search_context, results, precursor_dict, params)
    
    # Merge mass error plots
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    output_path = joinpath(ms1_mass_error_folder, "ms1_mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @user_warn "Could not clear existing file: $e"
    end

    if !isempty(results.ms1_mass_plots)
        save_multipage_pdf(results.ms1_mass_plots, output_path)
        empty!(results.ms1_mass_plots)
    end
end

