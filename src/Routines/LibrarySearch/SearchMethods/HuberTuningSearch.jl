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

    function HuberTuningSearchParameters(params::Any)
       #hp = params[:huber_delta_params]
        dp = params[:deconvolution_params]
        prec_estimation = PartialPrecCapture()  # Always use partial capture for Huber tuning
        
        delta0 = params[:deconvolution_params]["huber_delta0"];
        delta_exp = params[:deconvolution_params]["huber_delta_exp"];
        delta_iters = params[:deconvolution_params]["huber_delta_iters"];
        huber_δs = Float32[delta0*(delta_exp^i) for i in range(1, delta_iters)];
        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(2)),
            2,
            UInt8(params[:quant_search_params]["max_frag_rank"]),
            1.0f0,
            Set(2),
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            huber_δs,#Float32.(hp["delta_grid"]),
            Float32(50),#Float32(hp["min_pct_diff"]),
            0.001f0,#Float32(hp["q_value_threshold"]),
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

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file for Huber tuning.
"""
function process_file!(
    results::HuberTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table
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
    ::Arrow.Table
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

#==========================================================
Helper Methods
==========================================================#

"""
Get best PSMs from precursor dictionary for tuning.
"""
function get_best_psms(search_context::SearchContext, q_value_threshold::Float32)
    prec_dict = getPrecursorDict(search_context)
    is_decoy = getIsDecoy(getPrecursors(getSpecLib(search_context)))#[:is_decoy]
    
    N = length(prec_dict)
    df = DataFrame(
        precursor_idx = UInt32[],
        ms_file_idx = UInt32[],
        scan_idx = UInt32[],
        best_prob = Float32[]
    )
    
    for (pid, value) in pairs(prec_dict)
        push!(df, (
            pid,
            value[:best_ms_file_idx],
            value[:best_scan_idx],
            value[:best_prob]
        ))
    end
    
    df[!, :target] = .!is_decoy[df.precursor_idx]
    
    # Calculate q-values
    sort!(df, :best_prob, rev=true)
    df[!, :q_value] = zeros(Float32, nrow(df))
    target_count = 0
    
    for i in 1:nrow(df)
        target_count += df[i, :target]
        df[i, :q_value] = (i - target_count) / i
    end
    
    return filter(row -> row.q_value <= q_value_threshold, df)
end

"""
Create mapping between scan indices and precursor IDs.
"""
function get_scan_precursor_mapping(psms::DataFrame)
    scan_to_prec = Dict{UInt32, Vector{UInt32}}()
    for row in eachrow(psms)
        if !haskey(scan_to_prec, row.scan_idx)
            scan_to_prec[row.scan_idx] = UInt32[]
        end
        push!(scan_to_prec[row.scan_idx], row.precursor_idx)
    end
    return scan_to_prec
end

"""
Perform Huber tuning search on a single file.
"""
function perform_huber_search(
    spectra::Arrow.Table,
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    #scan_idxs::Set{UInt32},
    #prec_set::Set{Tuple{UInt32, UInt32}},
    search_context::SearchContext,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            tuning_results = process_scans_for_huber!(
                last(thread_task),
                #scan_idxs,
                scan_to_prec,
                #prec_set,
                spectra,
                search_context,
                search_data,
                params,
                ms_file_idx
            )
            
            return DataFrame(tuning_results)
        end
    end
    
    return vcat(fetch.(tasks)...)
end
#=
"""
Process scans for Huber tuning.
"""
function process_scans_for_huber!(
    scan_range::Vector{Int64},
    scan_idxs::Set{UInt32},
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    spectra::Arrow.Table,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64
)
    tuning_results = Dict(
        :precursor_idx => UInt32[],
        :scan_idx => UInt32[],
        :weight => Float32[],
        :huber_δ => Float32[]
    )
    
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    rt_index = buildRtIndex(DataFrame(Arrow.Table(search_context.rt_index_paths[][ms_file_idx])), bin_rt_size = 0.1)
    for scan_idx in scan_range
        process_scan_for_huber!(
            scan_idx,
            scan_idxs,
            scan_to_prec,
            rt_index,
            spectra,
            search_context,
            search_data,
            params,
            ms_file_idx,
            Hs,
            weights,
            precursor_weights,
            residuals,
            tuning_results
        )
    end
    
    return tuning_results
end
=#

"""
Process scans for Huber tuning with RT bin caching.
"""
function process_scans_for_huber!(
    scan_range::Vector{Int64},
    #scan_idxs::Set{UInt32},
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    #prec_set::Set{Tuple{UInt32, UInt32}},
    spectra::Arrow.Table,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64
)
    tuning_results = Dict(
        :precursor_idx => UInt32[],
        :scan_idx => UInt32[],
        :weight => Float32[],
        :huber_δ => Float32[]
    )
    
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    
    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    ion_idx = 0
    
    # Get RT index
    rt_index = buildRtIndex(
        DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
        bin_rt_size = 0.1)
    
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    for scan_idx in scan_range
        scan_idx ∉ keys(scan_to_prec) && continue
        
        msn = spectra[:msOrder][scan_idx]
        msn ∉ params.spec_order && continue
        
        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(spectra[:retentionTime][scan_idx])
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
        
        # Check for m/z change
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]
        
        # Update transitions if RT window or m/z changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new
            
            ion_idx, _ = selectTransitions!(
                getIonTemplates(search_data),
                RTIndexedTransitionSelection(),
                PartialPrecCapture(),
                getFragmentLookupTable(getSpecLib(search_context)),
                getPrecIds(search_data),
                getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
                getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
                getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
                getIsoSplines(search_data),
                getQuadTransmissionFunction(
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    spectra[:centerMz][scan_idx],
                    spectra[:isolationWidthMz][scan_idx]
                ),
                getPrecursorTransmission(search_data),
                getIsotopes(search_data),
                params.n_frag_isotopes,
                params.max_frag_rank,
                rt_index,
                irt_start,
                irt_stop,
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                block_size = 10000
            )
        end
        
        # Match peaks and process if enough matches found
        nmatches, nmisses = match_peaks_for_huber!(
            search_data,
            ion_idx,
            spectra,
            scan_idx,
            search_context,
            ms_file_idx
        )
        
        nmatches ≤ 2 && continue
        
        # Process delta values for this scan
        process_delta_values!(
            params.delta_grid,
            Hs,
            weights,
            precursor_weights,
            residuals,
            search_data,
            nmatches,
            nmisses,
            params,
            scan_idx,
            scan_to_prec,
            tuning_results
        )
    end
    
    return tuning_results
end

"""
Process single scan for Huber tuning.
"""
function process_scan_for_huber!(
    scan_idx::Int,
    #scan_idxs::Set{UInt32},
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    rt_index::retentionTimeIndex,
    spectra::Arrow.Table,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::HuberTuningSearchParameters,
    ms_file_idx::Int64,
    Hs::SparseArray,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32},
    residuals::Vector{Float32},
    tuning_results::Dict
)
    scan_idx ∉ scan_idxs && return
    
    msn = spectra[:msOrder][scan_idx]
    msn ∉ params.spec_order && return
    
    # Select transitions
    ion_idx = select_transitions_for_huber!(
        search_data,
        search_context,
        scan_idx,
        rt_index,
        ms_file_idx,
        scan_to_prec,
        spectra,
        params
    )
    
    # Match peaks
    nmatches, nmisses = match_peaks_for_huber!(
        search_data,
        ion_idx,
        spectra,
        scan_idx,
        search_context,
        ms_file_idx
    )
    
    nmatches ≤ 2 && return
    
    # Process each δ value
    process_delta_values!(
        params.delta_grid,
        Hs,
        weights,
        precursor_weights,
        residuals,
        search_data,
        nmatches,
        nmisses,
        params,
        scan_idx,
        scan_to_prec,
        tuning_results
    )
end

"""
Estimate optimal delta from tuning results.
"""
function estimate_optimal_delta(
    psms::DataFrame,
    delta_grid::Vector{Float32},
    min_pct_diff::Float32
)
    # Group by precursor/scan
    gpsms = groupby(psms, [:precursor_idx, :scan_idx])
    
    # Process each curve
    curves = combine(gpsms) do sdf
        process_huber_curve(sdf[!, :weight], sdf[!, :huber_δ])
    end
    
    # Filter curves
    filter!(row -> row.n == length(delta_grid), curves)
    filter!(row -> row.wdiff > (min_pct_diff/100), curves)
    filter!(row -> !ismissing(row.huber50), curves)
    
    # Get median delta
    curves[!, :huber50] = ceil.(Int, curves[!, :huber50])
    huber_hist = combine(groupby(curves, :huber50), nrow)
    sort!(huber_hist, :huber50)
    
    huber_hist[!, :prob] = huber_hist[!, :nrow] ./ sum(huber_hist[!, :nrow])
    huber_hist[!, :cum_prob] = Float32.(cumsum(huber_hist[!, :prob]))
    
    return get_median_huber_delta(
        huber_hist[!, :cum_prob],
        huber_hist[!, :huber50]
    )
end



"""
Process a single Huber curve to extract statistics.
"""
function process_huber_curve(
    weights::AbstractVector{Float32},
    huber_δs::AbstractVector{Float32}
)
    min_w, max_w = minimum(weights), maximum(weights)
    huber50 = missing
    w50 = min_w + (max_w - min_w)/2
    
    if length(weights) > 1
        for i in 1:(length(weights)-1)
            if (w50 >= weights[i]) && (w50 <= weights[i + 1])
                huber50 = huber_δs[i] + (huber_δs[i + 1] - huber_δs[i])/2
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

"""
Get median Huber delta value from cumulative probability distribution.
"""
function get_median_huber_delta(
    cum_prob::Vector{Float32},
    δ::Vector{Int64};
    cum_prob_threshold::Float32 = 0.5f0
)
    N = length(cum_prob)
    N == 1 && return first(δ)
    
    for i in 1:(N-1)
        if cum_prob[i+1] >= cum_prob_threshold
            x1, x2 = cum_prob[i], cum_prob[i+1]
            y1, y2 = δ[i], δ[i+1]
            slope = (y2 - y1)/(x2 - x1)
            midpoint = (x2 + x1)/2
            return y1 + slope*(midpoint - x1)
        end
    end
    
    @warn "Could not estimate huber delta"
    return first(δ)
end

"""
Process a range of delta values for a single scan.
"""
function process_delta_values!(
    delta_grid::Vector{Float32},
    Hs::SparseArray,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32},
    residuals::Vector{Float32},
    search_data::SearchDataStructures,
    nmatches::Int,
    nmisses::Int,
    params::HuberTuningSearchParameters,
    scan_idx::Int64,
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    tuning_results::Dict
)
    # Build design matrix
    buildDesignMatrix!(
        Hs,
        getIonMatches(search_data),
        getIonMisses(search_data),
        nmatches,
        nmisses,
        getIdToCol(search_data)
    )
    
    # Process each delta value
    for δ in delta_grid
        # Resize arrays if needed
        if getIdToCol(search_data).size > length(weights)
            new_entries = getIdToCol(search_data).size - length(weights) + 1000
            resize!(weights, length(weights) + new_entries)
            resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
            append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
        end
        
        # Initialize weights
        for i in 1:getIdToCol(search_data).size
            weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
                precursor_weights[getIdToCol(search_data).keys[i]]
        end
        
        # Solve deconvolution problem
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            δ,
            params.lambda,
            params.max_iter_newton,
            params.max_iter_bisection,
            params.max_iter_outer,
            params.accuracy_newton,
            params.accuracy_bisection,
            10.0,
            params.max_diff
        )
        
        # Record results
        for i in 1:getIdToCol(search_data).size
            id = getIdToCol(search_data).keys[i]
            colid = getIdToCol(search_data)[id]
            
            # Update precursor weights
            precursor_weights[id] = weights[colid]
            
            # Record if this is a target PSM
            if id ∈ scan_to_prec[scan_idx]
                push!(tuning_results[:precursor_idx], id)
                push!(tuning_results[:weight], weights[colid])
                push!(tuning_results[:huber_δ], δ)
                push!(tuning_results[:scan_idx], UInt32(scan_idx))
            end
        end
    end
    
    reset!(getIdToCol(search_data))
    reset!(Hs)
end

"""
Select transitions for Huber tuning.
"""
function select_transitions_for_huber!(
    search_data::SearchDataStructures,
    search_context::SearchContext,
    scan_idx::Int,
    rt_index::retentionTimeIndex,
    ms_file_idx::Int,
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    spectra::Arrow.Table,
    params::HuberTuningSearchParameters
)

    return selectTransitions!(
        getIonTemplates(search_data),
        RTIndexedTransitionSelection(),
        PartialPrecCapture(),
        getFragmentLookupTable(getSpecLib(search_context)),
        getPrecIds(search_data),
        getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
        getIsoSplines(search_data),
        getQuadTransmissionFunction(
            getQuadTransmissionModel(search_context, ms_file_idx),
            spectra[:centerMz][scan_idx],
            spectra[:isolationWidthMz][scan_idx]
        ),
        getPrecursorTransmission(search_data),
        getIsotopes(search_data),
        params.n_frag_isotopes,
        params.max_frag_rank,
        rt_index,#getRtIndex(search_context),
        getIrtStart(search_context),
        getIrtStop(search_context),
        (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
        block_size = 10000
    )
end

"""
Match peaks for Huber tuning.
"""
function match_peaks_for_huber!(
    search_data::SearchDataStructures,
    ion_idx::Int,
    spectra::Arrow.Table,
    scan_idx::Int,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    return matchPeaks!(
        getIonMatches(search_data),
        getIonMisses(search_data),
        getIonTemplates(search_data),
        ion_idx,
        spectra[:mz_array][scan_idx],
        spectra[:intensity_array][scan_idx],
        getMassErrorModel(search_context, ms_file_idx),
        spectra[:highMz][scan_idx],
        UInt32(scan_idx),
        UInt32(ms_file_idx)
    )
end
