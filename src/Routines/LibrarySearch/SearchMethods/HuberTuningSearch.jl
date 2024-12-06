"""
    HuberTuningSearch

Search method for optimizing Huber loss parameters for spectral deconvolution.

This search:
1. Tests multiple Huber delta values on high-confidence PSMs
2. Analyzes weight distributions across delta values
3. Selects optimal delta based on cumulative probability
4. Updates deconvolution parameters in search context

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :deconvolution_params => Dict(
        "lambda" => 0.1,
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 50,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    ),
    :quant_search_params => Dict(
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "min_y_count" => 1
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
    huber_param::Base.Ref{Float32}  # Optimal Huber delta
    weight_distributions::Dict{Int64, DataFrame}  # Weight distributions per file
end

"""
Parameters for Huber tuning search.
"""
struct HuberTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Basic search parameters
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
    
    # Deconvolution parameters
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Huber tuning parameters
    huber_deltas::Vector{Float32}
    min_y_count::Int64
    min_percent_diff::Float32
    cum_prob_threshold::Float32
    bin_rt_size::Float32
    prec_estimation::P

    function HuberTuningSearchParameters(params::Any)
        dp = params[:deconvolution_params]
        qp = params[:quant_search_params]
        prec_estimation = true ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(0)),  # No isotope error for Huber tuning
            UInt8(3),  # Default min index score
            Int64(3),  # Default min fragment count
            0.1f0,    # Default spectral contrast
            -3.0f0,   # Default log2 ratio
            (3, 5),   # Default topN of M
            UInt8(3), # Default max rank
            Int64(qp["n_frag_isotopes"]),
            UInt8(qp["max_frag_rank"]),
            1.0f0,    # Full sampling
            Set(2),   # MS2 only
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            collect(Float32, range(1.0f0, 100.0f0, length=20)),  # Default delta grid
            Int64(qp["min_y_count"]),
            10.0f0,   # Default minimum percent difference
            0.5f0,    # Default cumulative probability threshold
            0.1f0,    # Default RT bin size
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::HuberTuningSearch, params::Any) = HuberTuningSearchParameters(params)

function init_search_results(
    ::HuberTuningSearchParameters,
    ::SearchContext,
    ::Int64
)
    return HuberTuningSearchResults(
        Ref{Float32}(0.0f0),
        Dict{Int64, DataFrame}()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the Huber tuning search.
"""
function process_file!(
    results::HuberTuningSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:HuberTuningSearchParameters}

    try
        # Set up RT index
        rt_index = build_rt_index(
            search_context,
            ms_file_idx,
            params.bin_rt_size
        )

        # Get high confidence PSMs for tuning
        psm_set = get_tuning_psms(
            search_context,
            ms_file_idx
        )

        # Perform Huber grid search
        psms = huber_grid_search(
            spectra,
            rt_index,
            psm_set,
            search_context,
            params,
            ms_file_idx
        )

        # Process results and store weight distributions
        results.weight_distributions[ms_file_idx] = process_huber_results(
            psms,
            params.huber_deltas,
            params.min_percent_diff
        )

    catch e
        @warn "Huber tuning failed for file" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

"""
Store results in search context.
"""
function process_search_results!(
    ::HuberTuningSearchResults,
    ::P,
    ::SearchContext,
    ::Int64
) where {P<:HuberTuningSearchParameters}
    return nothing
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::HuberTuningSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:HuberTuningSearchParameters}
    
    @info "Determining optimal Huber parameter..."
    
    # Combine weight distributions
    combined_dist = combine_weight_distributions(
        results.weight_distributions,
        params.huber_deltas
    )
    
    # Get optimal parameter
    optimal_delta = get_median_huber_delta(
        combined_dist[!, :cum_prob],
        combined_dist[!, :huber_delta],
        cum_prob_threshold=params.cum_prob_threshold
    )
    
    results.huber_param[] = optimal_delta
    @info "Selected Huber delta: $optimal_delta"
end

"""
No cleanup needed between files.
"""
function reset_results!(::HuberTuningSearchResults)
    return nothing
end

#==========================================================
Search Implementation Methods
==========================================================#

"""
Build RT index for a file.
"""
function build_rt_index(
    search_context::SearchContext,
    ms_file_idx::Int64,
    bin_size::Float32
)
    rt_paths = getRtIndexPaths(search_context)
    fname = getParsedFileName(search_context, ms_file_idx)
    rt_df = DataFrame(Arrow.Table(rt_paths[fname]))
    
    return buildRtIndex(rt_df, bin_rt_size=bin_size)
end

"""
Get high confidence PSMs for tuning.
"""
function get_tuning_psms(
    search_context::SearchContext,
    ms_file_idx::Int64
)
    prec_dict = getPrecursorDict(search_context)
    is_decoy = getPrecursors(getSpecLib(search_context))[:is_decoy]
    
    best_psms = getPsmsForHuberEstimation(
        prec_dict,
        is_decoy,
        q_value_threshold=0.01f0
    )
    
    return Set(zip(best_psms[!, :precursor_idx], best_psms[!, :scan_idx])), 
           Set(best_psms[!, :scan_idx])
end

"""
Perform Huber grid search over delta values.
"""
function huber_grid_search(
    spectra::Arrow.Table,
    rt_index::retentionTimeIndex{Float32, Float32},
    psm_set::Tuple{Set{Tuple{UInt32,UInt32}}, Set{UInt32}},
    search_context::SearchContext,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            return huber_tuning_search(
                spectra,
                last(thread_task),
                psm_set[1],  # prec_set
                psm_set[2],  # scan_idxs
                search_context,
                search_data,
                rt_index,
                params,
                ms_file_idx
            )
        end
    end
    
    return vcat(fetch.(tasks)...)
end

"""
Core Huber tuning search implementation.
"""
function huber_tuning_search(
    spectra::Arrow.Table,
    thread_task::Vector{Int64},
    prec_set::Set{Tuple{UInt32,UInt32}},
    scan_idxs::Set{UInt32},
    search_context::SearchContext,
    search_data::SearchDataStructures,
    rt_index::retentionTimeIndex{Float32, Float32},
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    residuals = getResiduals(search_data)
    
    tuning_results = Dict(
        :precursor_idx => UInt32[],
        :scan_idx => UInt32[],
        :weight => Float32[],
        :huber_delta => Float32[]
    )

    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)

    # Process each scan
    for scan_idx in thread_task
        process_scan_for_huber!(
            scan_idx,
            scan_idxs,
            prec_set,
            spectra,
            search_context,
            search_data,
            params,
            ms_file_idx,
            Hs,
            weights,
            residuals,
            isotopes,
            precursor_transmission,
            tuning_results
        )
    end
    
    return DataFrame(tuning_results)
end

"""
Process a single scan for Huber tuning.
"""
function process_scan_for_huber!(
    scan_idx::Int64,
    scan_idxs::Set{UInt32},
    prec_set::Set{Tuple{UInt32,UInt32}},
    spectra::Arrow.Table,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64,
    Hs::SparseArray,
    weights::Vector{Float32},
    residuals::Vector{Float32},
    isotopes::Vector{Float32},
    precursor_transmission::Vector{Float32},
    tuning_results::Dict{Symbol,Vector}
)
    scan_idx ∉ scan_idxs && return
    msn = spectra[:msOrder][scan_idx]
    msn ∉ params.spec_order && return
    
    # Select transitions
    ion_idx = select_huber_transitions!(
        search_data,
        search_context,
        scan_idx,
        spectra,
        isotopes,
        precursor_transmission
    )
    
    nmatches, nmisses = match_peaks!(
        search_data,
        ion_idx,
        spectra,
        scan_idx,
        ms_file_idx,
        search_context
    )
    
    nmatches ≤ 2 && return
    
    perform_huber_grid!(
        Hs,
        weights,
        residuals,
        search_data,
        nmatches,
        nmisses,
        params,
        scan_idx,
        prec_set,
        tuning_results
    )
    
    reset!(getIdToCol(search_data))
    reset!(Hs)
end

#==========================================================
Helper Methods
==========================================================#

"""
Process Huber results to get weight distributions.
"""
function process_huber_results(
    psms::DataFrame,
    huber_deltas::Vector{Float32},
    min_percent_diff::Float32
)
    gpsms = groupby(psms, [:precursor_idx, :scan_idx])
    
    # Process weight curves
    weight_curves = combine(gpsms) do sdf
        process_huber_loss_curve(sdf[!, :weight], sdf[!, :huber_delta])
    end
    
    # Apply filters
    filter!(x -> x.n == length(huber_deltas), weight_curves)
    filter!(x -> x.wdiff > (min_percent_diff/100), weight_curves)
    filter!(x -> !ismissing(x.huber50), weight_curves)
    
    weight_curves[!, :huber50] = ceil.(Int, weight_curves[!, :huber50])
    
    # Get distribution
    huber_hist = combine(groupby(weight_curves, :huber50), nrow)
    sort!(huber_hist, :huber50)
    huber_hist[!, :prob] = huber_hist[!, :nrow] ./ sum(huber_hist[!, :nrow])
    huber_hist[!, :cum_prob] = Float32.(cumsum(huber_hist[!, :prob]))
    
    return huber_hist
end

"""
Get median Huber delta from cumulative distribution.
"""
function get_median_huber_delta(
    cum_prob::Vector{Float32},
    delta::Vector{Int64};
    cum_prob_threshold::Float32=0.5f0
)
    N = length(cum_prob)
    N == 1 && return first(delta)
    
    for i in 1:N-1
        if cum_prob[i+1] >= cum_prob_threshold
            x1, x2 = cum_prob[i], cum_prob[i+1]
            y1, y2 = delta[i], delta[i+1]
            slope = (y2 - y1)/(x2 - x1)
            midpoint = (x2 + x1)/2
            return y1 + slope*(midpoint - x1)
        end
    end
    
    @warn "Could not estimate Huber delta"
    return first(delta)
end

"""
Process weight curve for a single precursor.
"""
function process_huber_loss_curve(
    weights::AbstractVector{Float32},
    huber_deltas::AbstractVector{Float32}
)
    min_w, max_w = extrema(weights)
    w50 = min_w + (max_w - min_w)/2
    huber50 = missing
    
    if length(weights) > 1
        for i in 1:length(weights)-1
            if (w50 >= weights[i]) && (w50 <= weights[i + 1])
                huber50 = huber_deltas[i] + (huber_deltas[i + 1] - huber_deltas[i])/2
            end
        end
    end
    
    return (
        min = min_w,
        max = max_w,
        n = length(weights),
        huber50 = huber50,
        w50 = w50,
        wdiff = (max_w - min_w)/min_w
    )
end