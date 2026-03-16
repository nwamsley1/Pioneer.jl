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
    max_psms_for_huber::Int64
    prec_estimation::P

    # Override Parameters
    huber_override_bool::Bool
    huber_override_delta::Float32

    function HuberTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        frag_params = params.search.fragment_settings
        # Calculate Huber delta grid (hardcoded defaults)
        delta0 = Float32(300)
        delta_exp = Float32(1.5)
        delta_iters = Int64(15)
        huber_δs = Float32[delta0 * (delta_exp^i) for i in range(-4,delta_iters+6)]
        # Hardcoded isotope error bounds (always (1,0))
        isotope_bounds = (1, 0)

        # Always use partial capture for Huber tuning
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()

        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Int64(frag_params.n_isotopes),  # Fixed n_frag_isotopes
            UInt8(frag_params.max_rank),  # Using max possible rank
            Set{Int64}([2]),

            Float32(0.0),     # lambda (no regularization)
            NoNorm(),         # reg_type
            Int64(50),        # max_iter_newton
            Int64(100),       # max_iter_bisection
            Int64(1000),      # max_iter_outer
            Float32(10),      # accuracy_newton
            Float32(10),      # accuracy_bisection
            Float32(0.01),    # max_diff
            
            huber_δs,
            10.0f0,  # Default min_pct_diff
            0.001f0, # Default q_value_threshold
            5000,    # max_psms_for_huber
            prec_estimation,

            false,              # huber_override_bool
            Float32(1055)       # huber_override_delta

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

    # Bypass: no first-pass PSMs available for Huber tuning
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
    
    # Bypass: hardcode large delta (effectively no Huber clipping)
    bypass_delta = Float32(1e6)
    results.huber_delta[] = bypass_delta
    setHuberDelta!(search_context, bypass_delta)
    search_context.deconvolution_stop_tolerance[] = Float32(1e-6)
end

function reset_results!(::HuberTuningSearchResults)
    return nothing
end
