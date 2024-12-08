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
end

"""
Parameters for NCE tuning search.
Configures NCE grid search and general search behavior.
"""
struct NceTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
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
    nce_grid::LinRange{Float32, Int64}
    nce_breakpoint::Float32
    max_q_val::Float32
    prec_estimation::P

    function NceTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            Float32(pp["frag_tol_ppm"]),
            UInt8(pp["min_index_search_score"]),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            Float32(pp["min_log2_matched_ratio"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            Int64(pp["n_frag_isotopes"]),
            UInt8(pp["max_frag_rank"]),
            Float32(pp["sample_rate"]),
            Set(2),
            LinRange(21.0f0, 40.0f0, 15),  
            NCE_MODEL_BREAKPOINT,
            0.01f0,
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
    ::SearchContext
)
    return NceTuningSearchResults(
        Dict{Int64, NceModel}(),
        DataFrame()
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
    spectra::Arrow.Table
) where {P<:NceTuningSearchParameters}

    try
        # Perform grid search
        psms = library_search(spectra, search_context, params, ms_file_idx)
        
        # Process and filter PSMs
        processed_psms = process_psms!(psms, spectra, search_context, params)
        
        # Fit and store NCE model
        nce_model = fit_nce_model(
            PiecewiseNceModel(0.0f0),
            processed_psms[!, :prec_mz],
            processed_psms[!, :nce],
            processed_psms[!, :charge],
            params.nce_breakpoint
        )
        
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
    ::Arrow.Table
) where {P<:NceTuningSearchParameters}
    
    setNceModel!(search_context, ms_file_idx, results.nce_models[ms_file_idx])
end

"""
Summarize results across all files.
"""
function summarize_results!(
    ::NceTuningSearchResults,
    ::P,
    ::SearchContext
) where {P<:NceTuningSearchParameters}
    
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

"""
Process and filter PSMs from search results.
"""
function process_psms!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::NceTuningSearchParameters
)
    # Add columns
    addPreSearchColumns!(
        psms,
        spectra,
        getPrecursors(getSpecLib(search_context))[:is_decoy],
        getPrecursors(getSpecLib(search_context))[:irt],
        getPrecursors(getSpecLib(search_context))[:prec_charge],
        spectra[:retentionTime],
        spectra[:TIC]
    )

    # Score PSMs
    scorePresearch!(psms)
    getQvalues!(psms[!, :prob], psms[!, :target], psms[!, :q_value])

    # Get best PSMs per precursor/scan
    spsms = combine(groupby(psms, [:precursor_idx, :scan_idx])) do group
        max_idx = argmax(group[!, :scribe])
        return group[max_idx:max_idx, :]
    end

    # Apply filters
    filter!(row -> row.target && row.q_value <= params.max_q_val, spsms)
    passing_precs = Set(spsms[!, :precursor_idx])
    filter!(row -> row.precursor_idx ∈ passing_precs, psms)

    # Add precursor info
    psms[!, :prec_mz] = [
        getPrecursors(getSpecLib(search_context))[:mz][pid]
        for pid in psms[!, :precursor_idx]
    ]
    
    # Select best PSMs
    psms[!, :best_psms] .= false
    for group in groupby(psms, :precursor_idx)
        best_idx = argmax(group[!, :scribe])
        group[best_idx, :best_psms] = true
    end
    filter!(row -> row.best_psms, psms)

    return psms
end