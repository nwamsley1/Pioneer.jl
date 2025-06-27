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
    NceTuningSearch

Search method for optimizing normalized collision energy (NCE) parameters.

This search:
1. Performs grid search over NCE values
2. Collects PSMs for each NCE value
3. Fits piecewise NCE models based on precursor m/z
4. Stores optimized models in SearchContext for use by other methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 2),
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "sample_rate" => 0.1,
        "abreviate_precursor_calc" => false
    )
)

# Execute search
results = execute_search(NceTuningSearch(), search_context, params)
```
"""
struct NceTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for NCE tuning search.
Holds NCE models and associated PSM data for each file.
"""
struct NceTuningSearchResults <: SearchResults
    nce_models::Dict{Int64, NceModel}
    nce_psms::DataFrame
    nce_plot_dir::String
end

"""
Parameters for NCE tuning search.
Configures NCE grid search and general search behavior.
"""
struct NceTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    
    # NCE specific parameters
    nce_grid::LinRange{Float32, Int64}
    nce_breakpoint::Float32
    max_q_val::Float32
    min_samples::Int64
    prec_estimation::P

    function NceTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        tuning_params = params.parameter_tuning
        frag_params = tuning_params.fragment_settings
        search_params = tuning_params.search_settings
        global_params = params.global_settings

        # Always use partial capture for NCE tuning
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        # Create NCE grid
        nce_grid = LinRange{Float32}(21.0f0, 40.0f0, 30)

        new{typeof(prec_estimation)}(
            # Core parameters
            (UInt8(0), UInt8(0)),  # Fixed isotope bounds for NCE tuning
            UInt8(frag_params.min_score),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1),  # Fixed max_best_rank
            Int64(1),  # Fixed n_frag_isotopes for NCE tuning
            UInt8(frag_params.max_rank),
            Float32(search_params.sample_rate),
            Set{Int64}([2]),
            Float32(frag_params.relative_improvement_threshold),
            
            # NCE specific parameters
            nce_grid,
            NCE_MODEL_BREAKPOINT,  # Assuming this is defined as a constant
            Float32(0.01),  # Fixed max_q_val for NCE tuning
            Int64(search_params.min_samples),
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::NceTuningSearch, params::Any) = NceTuningSearchParameters(params)

function init_search_results(
    ::NceTuningSearchParameters,
    search_context::SearchContext
)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)
    nce_model_plots_path = joinpath(qc_dir, "collision_energy_alignment")
    !isdir(nce_model_plots_path) && mkdir(nce_model_plots_path)
    return NceTuningSearchResults(
        Dict{Int64, NceModel}(),
        DataFrame(),
        nce_model_plots_path
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Main file processing method for NCE tuning search.
Performs grid search and fits NCE model.
"""
function process_file!(
    results::NceTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:NceTuningSearchParameters}

    try
        if typeof(getSpecLib(search_context))==FragmentIndexLibrary
            #No NCE tuning for basic FramgentIndexLibrary
            return nothing
        end
        processed_psms = DataFrame()
        for i in range(1, 10)
            # Perform grid search
            psms = library_search(spectra, search_context, params, ms_file_idx)   
            # Process and filter PSMs
            append!(processed_psms, process_psms!(psms, spectra, search_context, params))
            if size(processed_psms, 1)>params.min_samples
                break
            end
            if i == 10
                n = size(processed_psms, 1)
                @warn "Could not get collect enough psms for nce alignment. In 10 iterations collected $n samples"
            end
        end

        # Fit and store NCE model
        nce_model = fit_nce_model(
            PiecewiseNceModel(0.0f0),
            processed_psms[!, :prec_mz],
            processed_psms[!, :nce],
            processed_psms[!, :charge],
            params.nce_breakpoint
        )
        
        fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        # Create the main plot
        # Create the main plot with adjusted right margin
        p = plot(
            title = "NCE calibration for $fname",
            right_margin = 50Plots.px  # Add extra margin on the right
        )

        # Calculate bin range
        pbins = LinRange(minimum(processed_psms[!,:prec_mz]), maximum(processed_psms[!,:prec_mz]), 100)

        # Extend x-axis range to accommodate annotations
        x_range = maximum(pbins) - minimum(pbins)
        #plot_xlims = (minimum(pbins), maximum(pbins) + x_range * 0.15)  # Add 15% to x-axis

        # Plot each charge state with annotations
        for charge in sort(unique(processed_psms[!,:charge]))
            # Calculate the curve
            curve_values = nce_model.(pbins, charge)
            
            # Plot the line
            plot!(p, pbins, curve_values, 
                label = "+"*string(charge))
            
            # Add annotation at the rightmost point
            last_x = pbins[end]
            last_y = curve_values[end]
            
            # Add text annotation
            annotate!(p, [(last_x + x_range*0.02,  # Slight offset from end
                        last_y,
                        text("$(round(last_y, digits=1))", 
                                :left, 
                                8))])
        end
        # Save the plot
        savefig(p, joinpath(results.nce_plot_dir, "nce_model_"*"$fname"*".pdf"))  
        results.nce_models[ms_file_idx] = nce_model
        append!(results.nce_psms, processed_psms)

    catch e
        @warn "NCE tuning failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

"""
Store results in search context.
"""
function process_search_results!(
    results::NceTuningSearchResults,
    ::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:NceTuningSearchParameters}
    try
        setNceModel!(search_context, ms_file_idx, results.nce_models[ms_file_idx])
    catch
        nothing
    end
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::NceTuningSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:NceTuningSearchParameters}
    try
        output_path = joinpath(results.nce_plot_dir, "nce_alignment_plots.pdf")
        try
            if isfile(output_path)
                rm(output_path)
            end
        catch e
            @warn "Could not clear existing file: $e"
        end
        if !isempty(results.nce_plot_dir)
            merge_pdfs([x for x in readdir(results.nce_plot_dir, join=true) if endswith(x, ".pdf")],
                    output_path, 
                    cleanup=true)
        end
    catch
        nothing
    end
    # Could add NCE model statistics or plots here
    return nothing
end

"""
Reset results containers.
"""
function reset_results!(results::NceTuningSearchResults)
    empty!(results.nce_models)
    empty!(results.nce_psms)
end

#==========================================================
Helper Methods
==========================================================#
