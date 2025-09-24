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
    HuberTuningSearch

Search method to optimize the Huber loss parameter δ for spectral deconvolution.

This search:
1. Collects PSMs from previous search results
2. Tests multiple δ values to find optimal deconvolution performance
3. Analyzes weight distributions to determine best δ
4. Stores optimized δ in SearchContext for other methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :huber_delta_params => Dict(
        "delta_grid" => [10.0f0, 50.0f0, 100.0f0, 500.0f0, 1000.0f0],
        "min_pct_diff" => 10.0f0,
        "q_value_threshold" => 0.01f0
    ),
    :deconvolution_params => Dict(
        "lambda" => 0.1f0,
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 100,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    )
)

# Execute search
results = execute_search(HuberTuningSearch(), search_context, params)
```
"""
struct HuberTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for Huber tuning search.
"""
struct HuberTuningSearchResults <: SearchResults
    huber_delta::Base.Ref{Float32}
    tuning_psms::Vector{DataFrame}
end

"""
Parameters for Huber tuning search.
Configures deconvolution and δ selection.
"""
struct HuberTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    lambda::Float32
    reg_type::RegularizationType
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Huber tuning specific parameters
    delta_grid::Vector{Float32}
    min_pct_diff::Float32
    q_value_threshold::Float32
    prec_estimation::P

    # Override Parameters
    huber_override_bool::Bool
    huber_override_delta::Float32

    function HuberTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        deconv_params = params.optimization.deconvolution
        global_params = params.global_settings
        frag_params = params.quant_search.fragment_settings
        # Calculate Huber delta grid
        delta0 = Float32(deconv_params.ms2.huber_delta)
        delta_exp = Float32(deconv_params.huber_exp)
        delta_iters = Int64(deconv_params.huber_iters)
        huber_δs = Float32[delta0 * (delta_exp^i) for i in range(-4,delta_iters+6)]
        isotope_bounds = global_params.isotope_settings.err_bounds_quant_search

        # Always use partial capture for Huber tuning
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        reg_type = deconv_params.ms2.reg_type
        if reg_type == "none"
            reg_type = NoNorm()
        elseif reg_type == "l1"
            reg_type = L1Norm()
        elseif reg_type == "l2"
            reg_type = L2Norm()
        else
            reg_type = NoNorm()
            @user_warn "Warning. Reg type `$reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Int64(frag_params.n_isotopes),  # Fixed n_frag_isotopes
            UInt8(frag_params.max_rank),  # Using max possible rank
            Set{Int64}([2]),
            
            Float32(deconv_params.ms2.lambda),
            reg_type, 
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.bisection_iters),
            Int64(deconv_params.outer_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),
            
            huber_δs,
            10.0f0,  # Default min_pct_diff
            0.001f0, # Default q_value_threshold
            prec_estimation,

            global_params.huber_override.override_huber_delta_fit,
            Float32(global_params.huber_override.huber_delta)

        )
    end
end

#==========================================================
Interface Implementation  
==========================================================#

get_parameters(::HuberTuningSearch, params::Any) = HuberTuningSearchParameters(params)

function init_search_results(::HuberTuningSearchParameters, search_context::SearchContext)
    return HuberTuningSearchResults(
        Ref{Float32}(),
        Vector{DataFrame}()
    )
end


"""
Process a single file for Huber tuning.
"""
function process_file!(
    results::HuberTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:HuberTuningSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "HuberTuningSearch")
        return results  # Return early with unchanged results
    end

    if params.huber_override_bool==true
        return 
    end
    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Get PSMs to tune on
        best_psms = get_best_psms(search_context, params.q_value_threshold)
        file_psms = filter(row -> row.ms_file_idx == ms_file_idx, best_psms)
        isempty(file_psms) && return results
        
        # Create scan/precursor set and scan indices
        prec_set = Set(zip(file_psms[!, :precursor_idx], file_psms[!, :scan_idx]))
        scan_idxs = Set(file_psms[!, :scan_idx])
        
        # Create scan mapping
        scan_to_prec = get_scan_precursor_mapping(file_psms)
        
        # Perform Huber tuning search
        tuning_results = perform_huber_search(
            spectra,
            scan_to_prec,
            #scan_idxs,
            #prec_set,
            search_context,
            params,
            ms_file_idx
        )
        
        push!(results.tuning_psms, tuning_results)
        
    catch e
        # Handle failures gracefully using helper function
        handle_search_error!(search_context, ms_file_idx, "HuberTuningSearch", e, createFallbackResults!, results)
    end
    
    return results
end

"""
Create fallback results for a failed file in HuberTuningSearch.
Uses default Huber delta value for robust quantification.
"""
function createFallbackResults!(results::HuberTuningSearchResults, ms_file_idx::Int64)
    # Use default conservative Huber delta value
    results.huber_delta[] = 1.0f0  # Conservative default
    
    # Create empty tuning PSMs DataFrame
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        delta = Float32[]
    )
    push!(results.tuning_psms, empty_psms)
end

function process_search_results!(
    ::HuberTuningSearchResults,
    ::P, 
    ::SearchContext,
    ::Int64,
    ::MassSpecData
) where {P<:HuberTuningSearchParameters}
    return nothing
end

"""
Summarize results across all files and set optimal δ.
"""
function summarize_results!(
    results::HuberTuningSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:HuberTuningSearchParameters}
    
    if params.huber_override_bool==true
        default_delta = params.huber_override_delta
        results.huber_delta[] = default_delta
        setHuberDelta!(search_context, default_delta)
        return 
    end
    try
        # Combine all tuning PSMs
        all_psms = vcat(results.tuning_psms...)

        # Process results to get optimal δ
        optimal_delta = estimate_optimal_delta(
            all_psms,
            params.delta_grid,
            params.min_pct_diff
        )
        search_context.deconvolution_stop_tolerance[] = Float32(quantile(all_psms[!,:weight], 0.01)/100000)
        # Store results
        if params.huber_override_bool #If overriding the fit
            optimal_delta = params.huber_override_delta
        end
        results.huber_delta[] = optimal_delta
        setHuberDelta!(search_context, optimal_delta)
        
    catch e
        @user_warn "Failed to determine optimal Huber delta, using default" exception=e
        default_delta = params.huber_override_delta
        results.huber_delta[] = default_delta
        setHuberDelta!(search_context, default_delta)
    end
end

function reset_results!(::HuberTuningSearchResults)
    return nothing
end
