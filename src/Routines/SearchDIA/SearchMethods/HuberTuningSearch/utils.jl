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

#==========================================================
PSMs processing
==========================================================#
"""
    get_best_psms(search_context::SearchContext, q_value_threshold::Float32) -> DataFrame

Retrieve high-confidence PSMs for Huber parameter tuning.

# Arguments
- `search_context`: Contains spectral library and precursor data
- `q_value_threshold`: Q-value cutoff for PSM filtering

# Process
1. Extracts PSMs from precursor dictionary
2. Adds target/decoy status
3. Calculates q-values
4. Filters by q-value threshold

# Returns
DataFrame of filtered PSMs containing:
- precursor_idx, ms_file_idx, scan_idx
- best_prob, target status, q_value
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
    get_scan_precursor_mapping(psms::DataFrame) -> Dict{UInt32, Vector{UInt32}}

Create mapping between scan indices and their corresponding precursor IDs.

# Arguments
- `psms`: DataFrame containing PSMs with scan_idx and precursor_idx columns

# Returns
Dictionary where:
- Keys: Scan indices
- Values: Vector of precursor IDs for that scan
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

#==========================================================
Core Search Functions 
==========================================================#
"""
    perform_huber_search(spectra::MassSpecData, scan_to_prec::Dict{UInt32, Vector{UInt32}},
                        search_context::SearchContext, params::HuberTuningSearchParameters,
                        ms_file_idx::Int64) -> DataFrame

Execute Huber tuning search across MS data.

# Arguments
- `spectra`: MS/MS spectral data
- `scan_to_prec`: Mapping of scans to precursors
- `search_context`: Search context
- `params`: Huber tuning parameters
- `ms_file_idx`: MS file index

# Process
1. Partitions scans across threads
2. Processes each scan batch in parallel
3. Combines results into single DataFrame

# Returns
DataFrame containing tuning results across all scans
"""
function perform_huber_search(
    spectra::MassSpecData,
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

"""
    process_scans_for_huber!(scan_range::Vector{Int64}, scan_to_prec::Dict{UInt32, Vector{UInt32}},
                            spectra::MassSpecData, search_context::SearchContext,
                            search_data::SearchDataStructures, params::HuberTuningSearchParameters,
                            ms_file_idx::Int64) -> Dict

Process a batch of scans for Huber tuning with RT bin caching.

# Arguments
- `scan_range`: Range of scans to process
- `scan_to_prec`: Scan to precursor mapping
- `spectra`: MS/MS data
- `search_context`: Search context
- `search_data`: Thread-specific search data
- `params`: Search parameters
- `ms_file_idx`: MS file index

# Returns
Dictionary of tuning results with keys:
- precursor_idx, scan_idx, weight, huber_δ
"""
function process_scans_for_huber!(
    scan_range::Vector{Int64},
    #scan_idxs::Set{UInt32},
    scan_to_prec::Dict{UInt32, Vector{UInt32}},
    #prec_set::Set{Tuple{UInt32, UInt32}},
    spectra::MassSpecData,
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

    # Log delta grid info
    delta_grid = params.delta_grid
    n_deltas = length(delta_grid)
    @user_info "\nTesting $n_deltas delta values: [$(delta_grid[1]), ..., $(delta_grid[end])]"

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
        
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && continue
        
        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
        
        # Check for m/z change
        prec_mz_string_new = string(getCenterMz(spectra, scan_idx))
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]
        
        # Update transitions if RT window or m/z changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || (prec_mz_string_new != prec_mz_string)
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new

            ion_idx, _ = select_transitions_for_huber!(
                search_data, search_context, scan_idx,
                ms_file_idx, spectra, irt_start, irt_stop, rt_index, params
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

#==========================================================
Transition and Peak Matching Functions
==========================================================#
"""
    select_transitions_for_huber!(search_data::SearchDataStructures, search_context::SearchContext,
                                scan_idx::Int, rt_index::retentionTimeIndex, ms_file_idx::Int,
                                scan_to_prec::Dict{UInt32, Vector{UInt32}}, spectra::MassSpecData,
                                params::HuberTuningSearchParameters)

Select ion transitions for Huber tuning analysis.

# Arguments
- Standard search parameters and data structures
- RT index for efficient transition selection
- Current scan and file indices

# Returns
Index of selected transitions and success status
"""
function select_transitions_for_huber!(
    search_data::SearchDataStructures,
    search_context::SearchContext,
    scan_idx::Int,
    ms_file_idx::Int,
    spectra::MassSpecData,
    irt_start::Int64, irt_stop::Int64, rt_index::Any,
    params::HuberTuningSearchParameters
)
    nce_model = getNceModel(search_context, ms_file_idx)

    return selectTransitions!(
        getIonTemplates(search_data),
        RTIndexedTransitionSelection(),
        params.prec_estimation,
        getFragmentLookupTable(getSpecLib(search_context)),
        nce_model,
        getPrecIds(search_data),
        getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
        getIsoSplines(search_data),
        getQuadTransmissionFunction(
            getQuadTransmissionModel(search_context, ms_file_idx),
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx)
        ),
        getPrecursorTransmission(search_data),
        getIsotopes(search_data),
        params.n_frag_isotopes,
        params.max_frag_rank,
        rt_index,#getRtIndex(search_context),
        irt_start,
        irt_stop,
        (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
        block_size = 10000
    )
end

"""
    match_peaks_for_huber!(search_data::SearchDataStructures, ion_idx::Int,
                          spectra::MassSpecData, scan_idx::Int, search_context::SearchContext,
                          ms_file_idx::Int64)

Match theoretical peaks to observed peaks for Huber tuning.

# Returns
Tuple of (number of matches, number of misses)
"""
function match_peaks_for_huber!(
    search_data::SearchDataStructures,
    ion_idx::Int,
    spectra::MassSpecData,
    scan_idx::Int,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    return matchPeaks!(
        getIonMatches(search_data),
        getIonMisses(search_data),
        getIonTemplates(search_data),
        ion_idx,
        getMzArray(spectra, scan_idx),
        getIntensityArray(spectra, scan_idx),
        getMassErrorModel(search_context, ms_file_idx),
        getHighMz(spectra, scan_idx),
        UInt32(scan_idx),
        UInt32(ms_file_idx)
    )
end


"""
    process_delta_values!(delta_grid::Vector{Float32}, Hs::SparseArray, weights::Vector{Float32},
                         precursor_weights::Vector{Float32}, residuals::Vector{Float32},
                         search_data::SearchDataStructures, nmatches::Int, nmisses::Int,
                         params::HuberTuningSearchParameters, scan_idx::Int64,
                         scan_to_prec::Dict{UInt32, Vector{UInt32}}, tuning_results::Dict)
                         -> (Float64, Float64, Dict{Symbol, Int})

Process multiple delta values for a single scan.

# Process
1. Builds design matrix
2. For each delta value:
   - Initializes weights
   - Solves deconvolution problem
   - Records results for target PSMs

# Returns
Nothing (mutates tuning_results in place)
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
    # Build design matrix (once per scan)
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
            params.max_diff,
            NoNorm()
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

    return nothing
end

#==========================================================
Huber delta value optim
==========================================================#
"""
    estimate_optimal_delta(psms::DataFrame, delta_grid::Vector{Float32},
                         min_pct_diff::Float32) -> Float32

Estimate optimal Huber delta parameter from tuning results.

# Process
1. Groups PSMs by precursor/scan
2. Processes each weight curve
3. Filters curves based on criteria
4. Calculates median delta value
"""
function estimate_optimal_delta(
    psms::DataFrame,
    delta_grid::Vector{Float32},
    min_pct_diff::Float32
)
    # Group by precursor/scan
    gpsms = groupby(psms, [:precursor_idx, :scan_idx])
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
    process_huber_curve(weights::AbstractVector{Float32},
                       huber_δs::AbstractVector{Float32}) -> NamedTuple

Process a single Huber curve to extract statistics.

# Returns
NamedTuple containing:
- min/max weights
- number of points
- huber50 (interpolated δ at 50% weight)
- w50 (50% weight value)
- wdiff (relative weight difference)
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
            elseif (w50 <= weights[i]) && (w50 >= weights[i + 1])
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
    get_median_huber_delta(cum_prob::Vector{Float32}, δ::Vector{Int64};
                          cum_prob_threshold::Float32=0.5f0) -> Float32

Calculate median Huber delta from cumulative probability distribution.

Uses linear interpolation to find delta value at specified probability threshold.
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
    
    @user_warn "Could not estimate huber delta"
    return first(δ)
end

