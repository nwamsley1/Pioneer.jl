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
    sample_rate::Float32
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    lambda::Float32
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

    function HuberTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        deconv_params = params.optimization.deconvolution
        global_params = params.global_settings
        frag_params = params.quant_search.fragment_settings
        # Calculate Huber delta grid
        delta0 = Float32(deconv_params.huber_delta)
        delta_exp = Float32(deconv_params.huber_exp)
        delta_iters = Int64(deconv_params.huber_iters)
        huber_δs = Float32[delta0 * (delta_exp^i) for i in 1:delta_iters]
        
        isotope_bounds = global_params.isotope_settings.err_bounds_quant_search

        # Always use partial capture for Huber tuning
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Int64(frag_params.n_isotopes),  # Fixed n_frag_isotopes
            UInt8(frag_params.max_rank),  # Using max possible rank
            1.0f0,  # Full sampling
            Set{Int64}([2]),
            
            Float32(deconv_params.lambda),
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.newton_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),
            
            huber_δs,
            50.0f0,  # Default min_pct_diff
            0.001f0, # Default q_value_threshold
            prec_estimation
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
        @warn "Huber tuning failed for file $ms_file_idx" exception=e
    end
    
    return results
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
    
    try
        # Combine all tuning PSMs
        all_psms = vcat(results.tuning_psms...)
        
        # Process results to get optimal δ
        optimal_delta = estimate_optimal_delta(
            all_psms,
            params.delta_grid,
            params.min_pct_diff
        )
        search_context.deconvolution_stop_tolerance[] = Float32(quantile(all_psms[!,:weight], 0.01)/1000)
        #println("test ", search_context.deconvolution_stop_tolerance[] )
        #println("_")
        # Store results
        results.huber_delta[] = optimal_delta
        setHuberDelta!(search_context, optimal_delta)
        
    catch e
        throw(e)
        @warn "Failed to determine optimal Huber delta, using default" exception=e
        default_delta = 100000.0f0
        results.huber_delta[] = default_delta
        setHuberDelta!(search_context, default_delta)
    end
end

function reset_results!(::HuberTuningSearchResults)
    return nothing
end
