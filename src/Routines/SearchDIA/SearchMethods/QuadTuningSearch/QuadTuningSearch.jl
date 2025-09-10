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

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 0),  # Fixed for quad tuning
    :presearch_params => Dict(
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "quad_tuning_sample_rate" => 0.1,
        "min_quad_tuning_fragments" => 3,
        "min_quad_tuning_psms_per_thompson" => 250,
        "abreviate_precursor_calc" => false
    ),
    :deconvolution_params => Dict(
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 100,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    )
)

# Execute search
results = execute_search(QuadTuningSearch(), search_context, params)
```
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
    prec_estimation::P

    function QuadTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        tuning_params = params.parameter_tuning
        frag_params = tuning_params.fragment_settings
        search_params = tuning_params.search_settings
        deconv_params = params.optimization.deconvolution
        
        # Always use partial capture for quad tuning
        prec_estimation = PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            # Search parameters
            params.acquisition.quad_transmission.fit_from_data,
            (UInt8(0), UInt8(0)),  # Fixed isotope bounds for quad tuning
            # Handle min_score as either single value or array (use first value if array)
            begin
                min_score_raw = frag_params.min_score
                if min_score_raw isa Vector
                    UInt8(first(min_score_raw))
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
            Int64(get(search_params, :min_quad_tuning_fragments, 3)),  # Default if not specified
            Int64(get(search_params, :min_quad_tuning_psms_per_thompson, 250)),   # Default if not specified
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
    if params.fit_from_data == false
        @debug_l1 "QuadTuningSearch: Skipping quad model fitting (fit_from_data = false)"
        return nothing
    end

    # Get file name for debugging
    file_name = try
        getFileIdToName(getMSData(search_context), ms_file_idx)
    catch
        "file_$ms_file_idx"
    end
    
    @info "Quad tuning processing file: $file_name"
          
    # Log fitted mass error model from ParameterTuningSearch
    fitted_model = getMassErrorModel(search_context, ms_file_idx)
    
    try
        # Check if file has any scans
        if length(spectra) == 0
            @user_warn "Skipping quad tuning for $file_name - file contains no scans"
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
                @user_warn "Data type inconsistency detected in $file_name. Array types don't match expected schema. This may be due to dummy/test data with inconsistent typing. Skipping quad tuning."
                setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
                return results
            else
                rethrow(type_error)
            end
        end
        
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Adjust arrays for isotope variants
        adjust_precursor_arrays!(search_context)
        
        # Check window widths
        window_widths = check_window_widths(spectra)
        if length(window_widths) != 1
            @user_warn "Multiple window sizes detected: $(join(collect(window_widths), ';'))"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end
        window_width = first(window_widths)
        # Collect and process PSMs
        total_psms = collect_psms(spectra, search_context, results, params, ms_file_idx, parse(Float64, window_width))
        
        required_psms = params.min_quad_tuning_psms_per_thompson * parse(Float64, window_width)
        if nrow(total_psms) < required_psms
            @user_warn "Too few PSMs found for quad modeling. required_psms $required_psms and total_psms $(nrow(total_psms))"
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
        @user_warn "Quad transmission function fit failed for MS data file: $file_name. Error type: $(typeof(e)). Using fallback model: GeneralGaussModel(5.0, 0.0)"
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
        @user_warn "Could not clear existing file: $e"
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
