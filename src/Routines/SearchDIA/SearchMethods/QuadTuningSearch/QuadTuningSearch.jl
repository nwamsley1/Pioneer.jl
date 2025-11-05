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
    QuadTuningSearch

Search method for optimizing quadrupole transmission models.

This search:
1. Uses fitted MassErrorModel from ParameterTuningSearch for accurate fragment matching
2. Collects PSMs with extended precursor isotope patterns
3. Performs deconvolution to estimate relative isotope abundances
4. Fits transmission model based on isotope ratio deviations
5. Stores optimized models in SearchContext for other methods

Note: Mass tolerances are determined by the fitted MassErrorModel from ParameterTuningSearch,
not by parameter settings. This ensures data-driven, instrument-specific tolerances.

# Configuration
Quadrupole tuning parameters are configured in the parameter_tuning.quad_tuning section:

```json
{
    "parameter_tuning": {
        "quad_tuning": {
            "min_psms_per_thompson": 250,  // PSMs required per Thomson isolation width
            "min_fragments": 3,             // Minimum fragments for PSM acceptance
            "initial_percent": 2.5,         // Initial sampling percentage
            "min_initial_scans": 5000       // Minimum scans for initial sample
        }
    }
}
```

The search automatically:
- Calculates dynamic PSM requirements based on isolation window width
- Uses progressive sampling starting at max(initial_percent, min_initial_scans/total_scans)
- Prioritizes scans by m/z distribution for optimal spectral coverage
"""
struct QuadTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for quadrupole tuning search.
"""
struct QuadTuningSearchResults <: SearchResults
    tuning_results::Vector{Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }}}
    quad_model::Base.Ref{QuadTransmissionModel}
    quad_plot_dir::String
    quad_model_plots::Vector{Plots.Plot}
    quad_data_plots::Vector{Plots.Plot}
end

"""
Parameters for quadrupole tuning search.
Configures deconvolution and quad model fitting.
"""
struct QuadTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    #Fit Quad Model from Data?
    fit_from_data::Bool

    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_log2_matched_ratio::Float32
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    max_frag_rank::UInt8
    n_frag_isotopes::Int64
    irt_tol::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    
    # Deconvolution parameters
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Quad tuning specific parameters
    min_quad_tuning_fragments::Int64
    min_quad_tuning_psms_per_thompson::Int64
    initial_percent::Float32
    prec_estimation::P

    function QuadTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        tuning_params = params.parameter_tuning
        frag_params = tuning_params.fragment_settings
        search_params = tuning_params.search_settings
        deconv_params = params.optimization.deconvolution
        
        # Always use partial capture for quad tuning
        prec_estimation = PartialPrecCapture()

        # Extract Quad tuning parameters with backwards compatibility
        quad_tuning_params = get(tuning_params, :quad_tuning, nothing)

        min_quad_tuning_fragments = if quad_tuning_params !== nothing && haskey(quad_tuning_params, :min_fragments)
            Int64(quad_tuning_params.min_fragments)
        else
            Int64(get(search_params, :min_quad_tuning_fragments, 3))  # Backwards compatibility
        end

        min_quad_tuning_psms_per_thompson = if quad_tuning_params !== nothing && haskey(quad_tuning_params, :min_psms_per_thompson)
            Int64(quad_tuning_params.min_psms_per_thompson)
        else
            Int64(get(search_params, :min_quad_tuning_psms_per_thompson, 250))  # Backwards compatibility
        end

        initial_percent = if quad_tuning_params !== nothing && haskey(quad_tuning_params, :initial_percent)
            Float32(quad_tuning_params.initial_percent)
        else
            Float32(2.5)  # Default value
        end


        new{typeof(prec_estimation)}(
            # Search parameters
            params.acquisition.quad_transmission.fit_from_data,
            (UInt8(0), UInt8(0)),  # Fixed isotope bounds for quad tuning
            # Handle min_score as either single value or array (use maximum value if array)
            begin
                min_score_raw = frag_params.min_score
                if min_score_raw isa Vector
                    UInt8(maximum(min_score_raw))  # Use most lenient (highest) score for Quad tuning
                else
                    UInt8(min_score_raw)
                end
            end,
            typemin(Float32),  # Minimum possible value for min_log2_matched_ratio
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1),  # Fixed max_best_rank
            UInt8(frag_params.max_rank),
            Int64(1),  # Fixed n_frag_isotopes for quad tuning
            typemax(Float32),  # Maximum possible value for irt_tol
            Set{Int64}([2]),
            Float32(frag_params.relative_improvement_threshold),
            
            # Deconvolution parameters
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.bisection_iters),
            Int64(deconv_params.outer_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),
            
            # Quad tuning specific parameters
            min_quad_tuning_fragments,
            min_quad_tuning_psms_per_thompson,
            initial_percent,
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

# Getters
getQuadModel(q::QuadTuningSearchResults) = q.quad_model[]

# Setters
function setQuadModel(q::QuadTuningSearchResults, model::Q) where {Q<:QuadTransmissionModel}
    q.quad_model[] = model
end

get_parameters(::QuadTuningSearch, params::Any) = QuadTuningSearchParameters(params)

function init_search_results(::QuadTuningSearchParameters, search_context::SearchContext)
    # Initialize empty tuning results vector in each search data structure
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)

    qpp = joinpath(qc_dir, "quad_transmission_model")
    !isdir(qpp) && mkdir(qpp)
    quad_data = joinpath(qpp, "quad_data")
    !isdir(quad_data) && mkdir(quad_data)
    quad_models = joinpath(qpp, "quad_models")
    !isdir(quad_models) && mkdir(quad_models)
    temp_data = Vector{Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}}()#(undef, length(getSearchData(search_context)))
    for i in range(1, length(getSearchData(search_context)))
        push!(temp_data, Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}())
    end
    return QuadTuningSearchResults(
        temp_data,
        Ref{QuadTransmissionModel}(),
        qpp,
        Plots.Plot[],
        Plots.Plot[]
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Main file processing method for quad tuning search.
"""
function process_file!(
    results::QuadTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData) where {P<:QuadTuningSearchParameters}

    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(1.0f0))

    # Get file name for debugging
    file_name = try
        getFileIdToName(getMSData(search_context), ms_file_idx)
    catch
        "file_$ms_file_idx"
    end

    # Check if quad model fitting is enabled BEFORE any expensive operations
    if params.fit_from_data == false
        return nothing
    end

    # Log fitted mass error model from ParameterTuningSearch
    fitted_model = getMassErrorModel(search_context, ms_file_idx)

    try
        # Check if file has any scans
        if length(spectra) == 0
            @user_warn "\nSkipping quad tuning for $file_name - file contains no scans"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end
        
        # Check for inconsistent array types (data quality issue)
        try
            # Test array access - this will fail if there are type mismatches
            test_scan_idx = findfirst(i -> getMsOrder(spectra, i) == 2, 1:length(spectra))
            if test_scan_idx !== nothing
                _ = getMzArray(spectra, test_scan_idx)
                _ = getIntensityArray(spectra, test_scan_idx)
            end
        catch type_error
            if isa(type_error, MethodError) || contains(string(type_error), "SubArray")
                @user_warn "\nData type inconsistency detected in $file_name. Array types don't match expected schema. This may be due to dummy/test data with inconsistent typing. Skipping quad tuning."
                setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
                return results
            else
                rethrow(type_error)
            end
        end

        # Adjust arrays for isotope variants
        adjust_precursor_arrays!(search_context)
        
        # Check window widths
        window_widths = check_window_widths(spectra)
        if length(window_widths) != 1
            @user_warn "\nMultiple window sizes detected: $(join(collect(window_widths), ';'))"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end
        window_width = first(window_widths)
        # Build scan priority index (metadata only, no peak data)
        scan_index = build_quad_scan_priority_index(spectra)

        # Calculate minimum PSMs required based on window width
        required_psms = params.min_quad_tuning_psms_per_thompson * parse(Float64, window_width)
        required_psms_int = round(Int, required_psms)

        total_psms, _, _ = progressive_quad_psm_collection!(
            scan_index,
            spectra,
            search_context,
            params,
            ms_file_idx;
            min_psms_required=required_psms_int,
            verbose=false
        )

        # Convert to DataFrame format expected by downstream processing
        if !isempty(total_psms)
            # Apply quad-specific processing if we got PSMs
            total_psms = process_quad_pipeline(total_psms, spectra, search_context, results, params, ms_file_idx, parse(Float64, window_width))
        else
            # Return empty DataFrame with correct schema
            total_psms = DataFrame(
                scan_idx = Int64[],
                precursor_idx = UInt32[],
                center_mz = Union{Float32, Missing}[],
                δ = Union{Float32, Missing}[],
                yt = Union{Float32, Missing}[],
                x0 = Union{Float32, Missing}[],
                x1 = Union{Float32, Missing}[],
                prec_charge = Union{UInt8, Missing}[],
                half_width_mz = Float32[]
            )
        end

        if nrow(total_psms) < required_psms
            @user_warn "Too few PSMs found for quad modeling. required_psms $required_psms and total_psms $(nrow(total_psms)) \n"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end

        # Plot charge states
        push!(results.quad_data_plots, plot_charge_distributions(total_psms, results, getFileIdToName(getMSData(search_context), ms_file_idx)))
        
        # Fit quad model
        window_width = parse(Float64, first(window_widths))
        fitted_model = RazoQuadModel(fit_quad_model(total_psms, window_width))
        setQuadModel(results, fitted_model)
        # Plot quad model
        push!(results.quad_model_plots, plot_quad_model(fitted_model, window_width, results, getFileIdToName(getMSData(search_context), ms_file_idx)))


    catch e
        @user_warn "\nQuad transmission function fit failed for MS data file: $file_name. Error: $e. Using fallback model: GeneralGaussModel(5.0, 0.0)"
        setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
    end

    return results
end

function process_search_results!(
    results::QuadTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:QuadTuningSearchParameters}
    
    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "QuadTuningSearch results processing")
        return nothing  # Return early 
    end
    
    if params.fit_from_data==true
        setQuadTransmissionModel!(search_context, ms_file_idx, getQuadModel(results))
    else
        setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))
    end
end

function summarize_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext
) where {P<:QuadTuningSearchParameters}
    
    models_path = joinpath(results.quad_plot_dir, "quad_models", "quad_model_plots.pdf")
    data_path = joinpath(results.quad_plot_dir, "quad_data", "quad_data_plots.pdf")
    try
        if isfile(models_path)
            safeRm(models_path, nothing)
        end
        if isfile(data_path)
            safeRm(data_path, nothing)
        end
    catch e
        @user_warn "\nCould not clear existing file: $e"
    end

    if !isempty(results.quad_model_plots)
        save_multipage_pdf(results.quad_model_plots, models_path)
        empty!(results.quad_model_plots)
    end

    if !isempty(results.quad_data_plots)
        save_multipage_pdf(results.quad_data_plots, data_path)
        empty!(results.quad_data_plots)
    end

    reset_precursor_arrays!(search_context)
    return nothing
end

function reset_results!(results::QuadTuningSearchResults)
    for r in results.tuning_results
        empty!(r)
    end
    return nothing
end
