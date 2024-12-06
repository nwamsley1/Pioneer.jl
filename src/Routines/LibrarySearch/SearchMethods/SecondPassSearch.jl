"""
    SecondPassSearch

Second pass search method implementing quantification using optimized parameters.

This search:
1. Uses tuned parameters from previous searches (Huber, NCE, Quad)
2. Performs detailed PSM identification with quantification
3. Applies isotope pattern analysis
4. Stores results for protein-level quantification

# Example Implementation
```julia
params = Dict(
    :deconvolution_params => Dict(
        "huber_delta" => 50.0,  # Will be overridden by tuned value
        "lambda" => 0.1,
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 50,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    ),
    :quant_search_params => Dict(
        "min_y_count" => 1,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10
    ),
    :isotope_params => Dict(
        "trace_type" => "full",  # or "simplified"
        "bin_rt_size" => 0.1
    )
)

results = execute_search(SecondPassSearch(), search_context, params)
```
"""
struct SecondPassSearch <: SearchMethod end

"""
Results container for second pass quantitative search.

# Fields
- `psms_paths::Dict{String, String}`: Maps file identifiers to paths of stored PSM results
- `summary_stats::Dict{String, Any}`: Summary statistics per file including coverage, quantification success rates
- `file_metrics::Dict{String, DataFrame}`: Quality metrics for each processed file
"""
struct SecondPassSearchResults <: SearchResults
    psms_paths::Dict{String, String}
    summary_stats::Dict{String, Any}
    file_metrics::Dict{String, DataFrame}
end

"""
Parameters for second pass search with quantification.

# Fields
## Search Parameters
- `isotope_err_bounds::Tuple{UInt8, UInt8}`: Bounds for isotope error consideration
- `min_index_search_score::UInt8`: Minimum score for index search matches
- `min_frag_count::Int64`: Minimum number of matching fragments
- `min_spectral_contrast::Float32`: Minimum spectral contrast score
- `min_log2_matched_ratio::Float32`: Minimum log2 ratio of matched intensity
- `min_topn_of_m::Tuple{Int64, Int64}`: Minimum top N of M matches
- `max_best_rank::UInt8`: Maximum rank for best matches
- `n_frag_isotopes::Int64`: Number of fragment isotopes to consider
- `max_frag_rank::UInt8`: Maximum fragment rank to consider
- `spec_order::Set{Int64}`: Set of MS levels to analyze

## Deconvolution Parameters
- `lambda::Float32`: Lambda parameter for Huber loss
- `max_iter_newton::Int64`: Maximum Newton iterations
- `max_iter_bisection::Int64`: Maximum bisection iterations
- `max_iter_outer::Int64`: Maximum outer loop iterations
- `accuracy_newton::Float32`: Newton method accuracy threshold
- `accuracy_bisection::Float32`: Bisection method accuracy threshold
- `max_diff::Float32`: Maximum difference threshold

## Quantification Parameters
- `min_y_count::Int64`: Minimum number of y ions
- `bin_rt_size::Float32`: RT bin size for indexing
- `isotope_trace_type::IsotopeTraceType`: Type of isotope trace analysis
- `prec_estimation::P`: Precursor estimation method
"""
struct SecondPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Quantification parameters
    min_y_count::Int64
    bin_rt_size::Float32
    isotope_trace_type::IsotopeTraceType
    prec_estimation::P

    function SecondPassSearchParameters(params::Any)
        dp = params[:deconvolution_params]
        qp = params[:quant_search_params]
        ip = params[:isotope_params]
        
        # Determine precursor estimation method
        prec_estimation = qp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        # Determine isotope trace type
        trace_type = ip["trace_type"] == "full" ? FullIsotopeTrace() : SimplifiedIsotopeTrace()
        
        new{typeof(prec_estimation)}(
            (UInt8(3), UInt8(0)),  # isotope_err_bounds for quantification
            UInt8(3),              # min_index_search_score
            Int64(qp["min_frag_count"]),
            Float32(qp["min_spectral_contrast"]),
            Float32(qp["min_log2_matched_ratio"]),
            (Int64(first(qp["min_topn_of_m"])), Int64(last(qp["min_topn_of_m"]))),
            UInt8(qp["max_best_rank"]),
            Int64(qp["n_frag_isotopes"]),
            UInt8(qp["max_frag_rank"]),
            Set(2),               # MS2 only
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            Int64(qp["min_y_count"]),
            Float32(ip["bin_rt_size"]),
            trace_type,
            prec_estimation
        )
    end
end

# Interface implementation
get_parameters(::SecondPassSearch, params::Any) = SecondPassSearchParameters(params)

"""
Initialize search results and create necessary output directories.

# Arguments
- `search_parameters::SecondPassSearchParameters`: Search parameters
- `search_context::SearchContext`: Search context
- `ms_file_idx::Int64`: MS file index

# Returns
- `SecondPassSearchResults`: Initialized results container

Creates output directories for:
- Second pass quantification results
- Quality control plots
- Summary statistics
"""
function init_search_results(
    ::SecondPassSearchParameters,
    search_context::SearchContext,
    ::Int64
)
    # Create output directories
    base_dir = getDataOutDir(search_context)
    quant_dir = joinpath(base_dir, "second_quant")
    qc_dir = joinpath(base_dir, "second_qc")
    
    # Ensure directories exist
    for dir in [quant_dir, qc_dir]
        !isdir(dir) && mkdir(dir)
    end
    
    return SecondPassSearchResults(
        Dict{String, String}(),
        Dict{String, Any}(),
        Dict{String, DataFrame}()
    )
end


#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the second pass quantitative search.

# Arguments
- `results::SecondPassSearchResults`: Results container
- `params::SecondPassSearchParameters`: Search parameters
- `search_context::SearchContext`: Search context
- `ms_file_idx::Int64`: MS file index
- `spectra::Arrow.Table`: MS spectra data

# Returns
- Modified `results` with updated PSMs and metrics

# Processing Steps
1. Build RT index for precursor matching
2. Update NCE model from context
3. Perform quantitative search
4. Process and filter PSMs
5. Add isotope pattern analysis
6. Save results to disk
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:SecondPassSearchParameters}

    try
        # Get file identifier
        fname = getParsedFileName(search_context, ms_file_idx)
        
        # Build RT index from stored indices
        rt_df = DataFrame(Arrow.Table(getRtIndexPaths(search_context)[fname]))
        rt_index = buildRtIndex(rt_df, bin_rt_size=params.bin_rt_size)
        
        # Update NCE model for fragment predictions
        frag_lookup = getFragmentLookupTable(getSpecLib(search_context))
        updateNceModel(frag_lookup, getNceModelModel(search_context, ms_file_idx))

        # Store file metrics
        results.file_metrics[fname] = DataFrame(
            metric = String[],
            value = Float64[]
        )

        # Perform quantitative search
        psms = perform_quant_search(
            spectra,
            rt_index,
            search_context,
            params,
            ms_file_idx
        )
        
        # Early exit if no PSMs found
        if isempty(psms)
            @warn "No PSMs found for file" ms_file_idx
            return results
        end

        # Process PSMs and add metadata
        psms = add_psm_metadata!(
            psms,
            spectra,
            search_context,
            params,
            ms_file_idx
        )

        # Process isotope patterns
        psms = process_isotope_patterns!(
            psms,
            spectra,
            search_context,
            params,
            ms_file_idx
        )

        # Calculate summary scores
        psms = calculate_summary_scores!(
            psms,
            params.isotope_trace_type
        )

        # Save results
        save_path = save_quant_results!(
            psms,
            search_context,
            ms_file_idx
        )
        
        # Update results container
        results.psms_paths[fname] = save_path
        update_summary_stats!(results, psms, fname)

    catch e
        @warn "Second pass search failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

"""
No additional per-file processing needed.
"""
function process_search_results!(
    ::SecondPassSearchResults,
    ::P,
    ::SearchContext,
    ::Int64
) where {P<:SecondPassSearchParameters}
    return nothing
end

"""
Calculate and store summary statistics across all files.

Computes:
- Total PSMs identified
- Quantification success rate
- Distribution of PSM scores
- Isotope pattern coverage
"""
function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}
    
    @info "Calculating summary statistics..."
    
    total_stats = Dict{String, Any}()
    
    # Combine PSMs from all files
    all_psms = DataFrame[]
    for path in values(results.psms_paths)
        push!(all_psms, DataFrame(Arrow.Table(path)))
    end
    
    if !isempty(all_psms)
        combined_psms = vcat(all_psms...)
        
        # Calculate global statistics
        total_stats["total_psms"] = nrow(combined_psms)
        total_stats["quant_success_rate"] = mean(combined_psms[!, :best_scan])
        total_stats["median_score"] = median(combined_psms[!, :score])
        
        # Store in results
        results.summary_stats["global"] = total_stats
    end
    
    @info "Summary complete" total_psms=get(total_stats, "total_psms", 0)
end

#==========================================================
Search Implementation Methods
==========================================================#

"""
Perform quantitative search across all scans.

# Arguments
- `spectra::Arrow.Table`: MS spectra data
- `rt_index::retentionTimeIndex`: RT index for precursor matching
- `search_context::SearchContext`: Search context
- `params::SecondPassSearchParameters`: Search parameters
- `ms_file_idx::Int64`: MS file index

# Returns
- DataFrame containing PSMs with quantification

Implements parallel processing across MS scans using thread tasks.
"""
function perform_quant_search(
    spectra::Arrow.Table,
    rt_index::retentionTimeIndex{Float32, Float32},
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    # Create thread tasks for parallel processing
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    
    # Launch parallel tasks
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            return second_pass_search(
                spectra,
                last(thread_task),
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx
            )
        end
    end
    
    # Combine results
    return vcat(fetch.(tasks)...)
end

"""
Core second pass search implementation per thread.

# Arguments
- `spectra::Arrow.Table`: MS spectra data
- `thread_task::Vector{Int64}`: Scan indices for this thread
- `rt_index::retentionTimeIndex`: RT index for precursor matching
- `search_context::SearchContext`: Search context
- `search_data::SearchDataStructures`: Thread-specific search data
- `params::SecondPassSearchParameters`: Search parameters
- `ms_file_idx::Int64`: MS file index

# Returns
- DataFrame containing PSMs found in assigned scans

Handles:
1. RT window calculation
2. Ion matching and scoring
3. Deconvolution using tuned Huber parameter
4. PSM collection and initial filtering
"""
function second_pass_search(
    spectra::Arrow.Table,
    thread_task::Vector{Int64},
    rt_index::retentionTimeIndex{Float32, Float32},
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    residuals = getResiduals(search_data)
    isotopes = getIsotopes(search_data)
    precursor_transmission = getPrecursorTransmission(search_data)
    
    # Initialize tracking variables
    last_val = 0
    cycle_idx = 0
    irt_start = irt_stop = 1
    prec_mz_string = ""
    
    # Get models and parameters
    rt_model = getRtIrtModel(search_context, ms_file_idx)
    huber_delta = getHuberDelta(search_context)
    mass_err_model = getMassErrorModel(search_context, ms_file_idx)
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    
    # Process each scan
    for scan_idx in thread_task
        # Basic validation
        (scan_idx == 0 || scan_idx > length(spectra[:mz_array])) && continue
        
        # Check MS level
        msn = spectra[:msOrder][scan_idx]
        msn ∉ params.spec_order && continue
        
        # Update cycle index for MS1 scans
        cycle_idx += (msn < 2)
        
        # Calculate RT windows
        irt = rt_model(spectra[:retentionTime][scan_idx])
        irt_err = getIrtErrors(search_context)[getParsedFileName(search_context, ms_file_idx)]
        
        new_irt_start = max(
            searchsortedfirst(
                rt_index.rt_bins,
                irt - irt_err,
                lt=(r,x)->r.lb<x
            ) - 1,
            1
        )
        
        new_irt_stop = min(
            searchsortedlast(
                rt_index.rt_bins,
                irt + irt_err,
                lt=(x,r)->r.ub>x
            ) + 1,
            length(rt_index.rt_bins)
        )
        
        # Get updated precursor m/z string
        new_prec_mz = string(spectra[:centerMz][scan_idx])
        new_prec_mz = new_prec_mz[1:min(length(new_prec_mz), 6)]
        
        # Check if we need to update transitions
        if new_irt_start != irt_start || 
           new_irt_stop != irt_stop || 
           new_prec_mz != prec_mz_string
            
            irt_start = new_irt_start
            irt_stop = new_irt_stop
            prec_mz_string = new_prec_mz
            
            # Select transitions
            ion_idx = selectTransitions!(
                getIonTemplates(search_data),
                RTIndexedTransitionSelection(),
                params.prec_estimation,
                getFragmentLookupTable(getSpecLib(search_context)),
                Vector{UInt32}(undef, 50000),  # Temporary precursor buffer
                getPrecursors(getSpecLib(search_context))[:mz],
                getPrecursors(getSpecLib(search_context))[:prec_charge],
                getPrecursors(getSpecLib(search_context))[:sulfur_count],
                getIsoSplines(search_data),
                getQuadTransmissionFunction(
                    quad_model,
                    spectra[:centerMz][scan_idx],
                    spectra[:isolationWidthMz][scan_idx]
                ),
                precursor_transmission,
                isotopes,
                params.n_frag_isotopes,
                params.max_frag_rank,
                rt_index,
                irt_start,
                irt_stop,
                (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                isotope_err_bounds=params.isotope_err_bounds,
                block_size=10000
            )
        end
        
        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            spectra[:mz_array][scan_idx],
            spectra[:intensity_array][scan_idx],
            mass_err_model,
            spectra[:highMz][scan_idx],
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )
        
        nmatches ≤ 2 && continue
        
        # Sort matches
        sort!(
            @view(getIonMatches(search_data)[1:nmatches]),
            by=x->(x.peak_ind, x.prec_id),
            alg=QuickSort
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
        
        # Resize arrays if needed
        if getIdToCol(search_data).size > length(weights)
            new_entries = getIdToCol(search_data).size - length(weights) + 1000
            resize!(weights, length(weights) + new_entries)
            resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
            append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
        end
        
        # Initialize weights from previous values
        for i in 1:getIdToCol(search_data).size
            weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
                search_data.precursor_weights[getIdToCol(search_data).keys[i]]
        end
        
        # Perform deconvolution
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            huber_delta,
            params.lambda,
            params.max_iter_newton,
            params.max_iter_bisection,
            params.max_iter_outer,
            params.accuracy_newton,
            params.accuracy_bisection,
            10.0,
            params.max_diff
        )
        
        # Update precursor weights
        for i in 1:getIdToCol(search_data).size
            search_data.precursor_weights[getIdToCol(search_data).keys[i]] = 
                weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]]
        end
        
        # Calculate spectral metrics
        getDistanceMetrics(weights, residuals, Hs, getSpectralScores(search_data))
        
        # Score matches
        ScoreFragmentMatches!(
            getUnscoredPsms(search_data),
            getIdToCol(search_data),
            getIonMatches(search_data),
            nmatches,
            mass_err_model,
            last(params.min_topn_of_m)
        )
        
        # Record PSMs
        last_val = Score!(
            getScoredPsms(search_data),
            getUnscoredPsms(search_data),
            getSpectralScores(search_data),
            weights,
            getIdToCol(search_data),
            cycle_idx,
            nmatches/(nmatches + nmisses),
            last_val,
            Hs.n,
            Float32(sum(spectra[:intensity_array][scan_idx])),
            scan_idx;
            min_spectral_contrast=params.min_spectral_contrast,
            min_log2_matched_ratio=params.min_log2_matched_ratio,
            min_y_count=params.min_y_count,
            min_frag_count=params.min_frag_count,
            max_best_rank=params.max_best_rank,
            min_topn=first(params.min_topn_of_m),
            block_size=500000
        )
        
        # Reset arrays
        for i in 1:Hs.n
            getUnscoredPsms(search_data)[i] = eltype(getUnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end
    
    return DataFrame(@view(getScoredPsms(search_data)[1:last_val]))
end


#==========================================================
PSM Processing Methods
==========================================================#

"""
Add metadata to PSMs from search results.

Adds:
- Retention time
- Charge states
- Decoy status
- CV fold assignment
"""
function add_psm_metadata!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    addSecondSearchColumns!(
        psms,
        spectra[:retentionTime],
        getPrecursors(getSpecLib(search_context))[:prec_charge],
        getPrecursors(getSpecLib(search_context))[:is_decoy],
        getPrecursors(getSpecLib(search_context))[:cv_fold]
    )
    
    psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)
    psms[!, :best_scan] = zeros(Bool, nrow(psms))
    
    return psms
end

"""
Process isotope patterns for quantification.

Analyzes:
- Isotope capture patterns
- Pattern completeness
- Peak overlap
"""
function process_isotope_patterns!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    getIsotopesCaptured!(
        psms,
        params.isotope_trace_type,
        getQuadTransmissionModel(search_context, ms_file_idx),
        psms[!, :scan_idx],
        getPrecursors(getSpecLib(search_context))[:prec_charge],
        getPrecursors(getSpecLib(search_context))[:mz],
        spectra[:centerMz],
        spectra[:isolationWidthMz]
    )
    
    filter!(row -> first(row.isotopes_captured) < 2, psms)
    
    return psms
end

"""
Calculate summary scores for PSM groups.

Computes:
- Goodness of fit
- Matched ratios
- Distance metrics
- Spectral contrast
"""
function calculate_summary_scores!(
    psms::DataFrame,
    isotope_trace_type::IsotopeTraceType
)
    initSummaryColumns!(psms)
    
    for sub_psms in groupby(psms, getPsmGroupbyCols(isotope_trace_type))
        getSummaryScores!(
            sub_psms,
            sub_psms[!, :weight],
            sub_psms[!, :gof],
            sub_psms[!, :matched_ratio],
            sub_psms[!, :fitted_manhattan_distance],
            sub_psms[!, :fitted_spectral_contrast],
            sub_psms[!, :y_count]
        )
    end
    
    filter!(row -> row.best_scan, psms)
    
    return psms
end

#==========================================================
Result Storage Methods
==========================================================#

"""
Save quantification results to file.

Creates Arrow file containing:
- PSM identifications
- Quantification values
- Quality metrics
"""
function save_quant_results!(
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    # Initialize probability columns
    for col in [:prob, :max_prob, :mean_prob, :min_prob]
        psms[!, col] = zeros(Float32, nrow(psms))
    end
    
    # Save to file
    fname = getParsedFileName(search_context, ms_file_idx)
    out_path = joinpath(getDataOutDir(search_context), "second_quant", fname * ".arrow")
    
    Arrow.write(out_path, psms)
    
    return out_path
end

"""
Update summary statistics with results from a file.

Tracks:
- PSM counts
- Quantification success rates
- Quality metrics
"""
function update_summary_stats!(
    results::SecondPassSearchResults,
    psms::DataFrame,
    fname::String
)
    file_stats = Dict{String, Any}()
    
    # Basic counts
    file_stats["total_psms"] = nrow(psms)
    file_stats["quantified_psms"] = count(psms[!, :best_scan])
    
    # Quality metrics
    if !isempty(psms)
        file_stats["median_score"] = median(psms[!, :score])
        file_stats["quant_success_rate"] = mean(psms[!, :best_scan])
        file_stats["median_fragments"] = median(psms[!, :n_matches])
    end
    
    results.summary_stats[fname] = file_stats
end

#==========================================================
Utility Methods
==========================================================#

"""
Get grouping columns for PSM summaries based on isotope trace type.
"""
function getPsmGroupbyCols(isotope_trace_type::IsotopeTraceType)
    if isotope_trace_type isa FullIsotopeTrace
        return [:precursor_idx, :scan_idx]
    else
        return [:precursor_idx]
    end
end

"""
Initialize summary columns in PSM DataFrame.
"""
function initSummaryColumns!(psms::DataFrame)
    for col in [:gof, :matched_ratio, :fitted_manhattan_distance, 
                :fitted_spectral_contrast, :y_count]
        psms[!, col] = zeros(Float32, nrow(psms))
    end
end

"""
Calculate summary scores for a group of PSMs.
"""
function getSummaryScores!(
    psms::AbstractDataFrame,
    weights::AbstractVector{Float32},
    gof::AbstractVector{Float32},
    matched_ratio::AbstractVector{Float32},
    fitted_manhattan::AbstractVector{Float32},
    fitted_contrast::AbstractVector{Float32},
    y_count::AbstractVector{Float32}
)
    # Find best scan based on weights and scores
    best_idx = argmax(weights .* psms[!, :score])
    psms[best_idx, :best_scan] = true
    
    # Calculate group-level metrics
    mean_weight = mean(weights)
    for i in eachindex(weights)
        gof[i] = calculate_gof(weights[i], mean_weight)
        matched_ratio[i] = calculate_matched_ratio(psms[i, :])
        fitted_manhattan[i] = calculate_fitted_manhattan(psms[i, :])
        fitted_contrast[i] = calculate_fitted_contrast(psms[i, :])
        y_count[i] = count_y_ions(psms[i, :])
    end
end

# Helper calculation functions
calculate_gof(weight::Float32, mean_weight::Float32) = 1.0f0 - abs(weight - mean_weight) / max(weight, mean_weight)
calculate_matched_ratio(psm::DataFrameRow) = psm.n_matches / (psm.n_matches + psm.n_misses)
calculate_fitted_manhattan(psm::DataFrameRow) = sum(abs.(psm.theoretical_intensities .- psm.observed_intensities))
calculate_fitted_contrast(psm::DataFrameRow) = cor(psm.theoretical_intensities, psm.observed_intensities)
count_y_ions(psm::DataFrameRow) = count(x -> startswith(x, "y"), psm.matched_ion_labels)
