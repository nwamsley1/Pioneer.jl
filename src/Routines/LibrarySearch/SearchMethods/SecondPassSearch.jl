"""
    SecondPassSearch

Second pass search method using optimized parameters from initial searches.

This search:
1. Uses optimized parameters from first pass search
2. Performs PSM identification with previously calculated Huber delta
3. Tracks retention time windows for efficient searching
4. Records chromatogram information for later analysis

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :second_pass_params => Dict(
        "min_y_count" => 1,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10
    )
)

# Execute search
results = execute_search(SecondPassSearch(), search_context, params)
```
"""
struct SecondPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for second pass search.
"""
struct SecondPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}          # PSMs for each file
end

"""
Parameters for second pass search.
"""
struct SecondPassSearchParameters{P<:PrecEstimation,I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
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

    # PSM filtering
    min_y_count::Int64
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::Int64
    
    # Precursor estimation strategy
    isotope_tracetype::I
    prec_estimation::P

    function SecondPassSearchParameters(params::Any) 
        sp = params[:quant_search_params]
        dp = params[:deconvolution_params]
        _ISOTOPE_TRACE_TYPE_ = nothing
        if params[:quant_search_params]["combine_isotope_traces"]
            _ISOTOPE_TRACE_TYPE_ = CombineTraces(Float32(params_[:quant_search_params]["min_fraction_transmitted"]))
            @warn "Combine Traces"
        else
            _ISOTOPE_TRACE_TYPE_ = SeperateTraces()
            @warn "Seperate Traces"
        end

        new{typeof(PartialPrecCapture()),typeof(_ISOTOPE_TRACE_TYPE_)}(
            (UInt8(3), UInt8(0)),  # Fixed isotope bounds
            Int64(sp["n_frag_isotopes"]),
            UInt8(sp["max_frag_rank"]),
            1.0f0,
            Set(2),
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            Int64(sp["min_y_count"]),
            Int64(sp["min_frag_count"]),
            Float32(sp["min_spectral_contrast"]),
            Float32(sp["min_log2_matched_ratio"]),
            (Int64(first(sp["min_topn_of_m"])), Int64(last(sp["min_topn_of_m"]))),
            Int64(sp["max_best_rank"]),
            _ISOTOPE_TRACE_TYPE_,
            PartialPrecCapture()
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
#==========================================================
Interface Implementation
==========================================================#

get_parameters(::SecondPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    second_pass_psms = joinpath(getDataOutDir(search_context), "second_pass_psms")
    !isdir(second_pass_psms) && mkdir(second_pass_psms)
    return SecondPassSearchResults(
        DataFrame()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file for second pass search.
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:SecondPassSearchParameters}

    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Get RT index
        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
            bin_rt_size = 0.1)

        # Perform second pass search
        psms = perform_second_pass_search(
            spectra,
            rt_index,
            search_context,
            params,
            ms_file_idx
        )

        results.psms[] = psms

    catch e
        @warn "Second pass search failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:SecondPassSearchParameters}

    try
        psms = results.psms[]
        
        addSecondSearchColumns!(psms, 
            spectra[:retentionTime],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge], 
            getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
            getPrecursors(getSpecLib(search_context))
            );
        # Add required columns
        psms[!, :charge2] = Vector{UInt8}(psms[!, :charge] .== 2)
        # Get isotope captures 
        getIsotopesCaptured!(
            psms,
            getIsotopeTraceType(params),#.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            spectra[:centerMz],
            spectra[:isolationWidthMz]
        )

        # Filter isotopes
        filter!(row -> first(row.isotopes_captured) < 2, psms)

        psms[!,:best_scan] = zeros(Bool, size(psms, 1));
        initSummaryColumns!(psms);
        for (key, gpsms) in pairs(groupby(psms, getPsmGroupbyCols(getIsotopeTraceType(params))))
            getSummaryScores!(
                gpsms, 
                gpsms[!,:weight],
                gpsms[!,:gof],
                gpsms[!,:matched_ratio],
                gpsms[!,:fitted_manhattan_distance],
                gpsms[!,:fitted_spectral_contrast],
                gpsms[!,:y_count]
            );
        end
        filter!(x->x.best_scan, psms);
        # Add post-integration features
        addPostIntegrationFeatures!(
            psms,
            getSequence(getPrecursors(getSpecLib(search_context))),#[:sequence],
            getStructuralMods(getPrecursors(getSpecLib(search_context))),#[:structural_mods],
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            getIrt(getPrecursors(getSpecLib(search_context))),#[:irt],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
            getMissedCleavages(getPrecursors(getSpecLib(search_context))),#[:missed_cleavages],
            spectra[:TIC],
            spectra[:mz_array],
            ms_file_idx,
            getRtIrtModel(search_context, ms_file_idx),
            getPrecursorDict(search_context)
        )

        # Save updated results
        temp_path = joinpath(
            getDataOutDir(search_context),
            "second_pass_psms",
            getParsedFileName(search_context, ms_file_idx) * ".arrow"
        )
        psms[!,:prob], psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
        Arrow.write(temp_path, psms)
        setSecondPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
    catch e
        @warn "Failed to process search results" ms_file_idx exception=e
        rethrow(e)
    end

    return nothing
end

"""
Perform second pass search for a file.
"""
function perform_second_pass_search(
    spectra::Arrow.Table,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            
            return process_scans!(
                last(thread_task),
                spectra,
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx
            )
        end
    end
    
    return vcat(fetch.(tasks)...)
end

"""
Process scans with RT bin caching.
"""
function process_scans!(
    scan_range::Vector{Int64},
    spectra::Arrow.Table,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    # Get working arrays
    Hs = getHs(search_data)
    weights = getTempWeights(search_data)
    precursor_weights = getPrecursorWeights(search_data)
    residuals = getResiduals(search_data)
    last_val = 0

    # RT bin tracking state
    irt_start, irt_stop = 1, 1
    prec_mz_string = ""
    ion_idx = 0
    cycle_idx = 0

    irt_tol = getIrtErrors(search_context)[ms_file_idx]

    for scan_idx in scan_range
        ((scan_idx < 1) || scan_idx > length(spectra[:mz_array])) && continue
        msn = spectra[:msOrder][scan_idx]
        if msn < 2
            cycle_idx += 1
        end
        msn ∉ params.spec_order && continue

        # Calculate RT window
        irt = getRtIrtModel(search_context, ms_file_idx)(spectra[:retentionTime][scan_idx])
        irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, lt=(r,x)->r.lb<x) - 1, 1)
        irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        # Check for m/z change
        prec_mz_string_new = string(spectra[:centerMz][scan_idx])
        prec_mz_string_new = prec_mz_string_new[1:min(length(prec_mz_string_new), 6)]

        # Update transitions if window changed
        if (irt_start_new != irt_start) || (irt_stop_new != irt_stop) || 
           (prec_mz_string_new != prec_mz_string)
            
            irt_start = irt_start_new
            irt_stop = irt_stop_new
            prec_mz_string = prec_mz_string_new

            ion_idx, _ = selectTransitions!(
                getIonTemplates(search_data),
                RTIndexedTransitionSelection(),
                params.prec_estimation,
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

        # Match peaks
        nmatches, nmisses = matchPeaks!(
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

        nmatches ≤ 2 && continue

        # Process scan
        buildDesignMatrix!(
            Hs,
            getIonMatches(search_data),
            getIonMisses(search_data),
            nmatches,
            nmisses,
            getIdToCol(search_data)
        )

        # Handle weights
        if getIdToCol(search_data).size > length(weights)
            resize_arrays!(search_data, weights)
        end

        initialize_weights!(search_data, weights, precursor_weights)
        
        # Solve deconvolution problem
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            getHuberDelta(search_context),
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
        update_precursor_weights!(search_data, weights, precursor_weights)

        # Score PSMs
        getDistanceMetrics(weights, residuals, Hs, getComplexSpectralScores(search_data))
        
        ScoreFragmentMatches!(
            getComplexUnscoredPsms(search_data),
            getIdToCol(search_data),
            getIonMatches(search_data),
            nmatches,
            getMassErrorModel(search_context, ms_file_idx),
            last(params.min_topn_of_m)
        )

        last_val = Score!(
            getComplexScoredPsms(search_data),
            getComplexUnscoredPsms(search_data),
            getComplexSpectralScores(search_data),
            weights,
            getIdToCol(search_data),
            cycle_idx,
            nmatches/(nmatches + nmisses),
            last_val,
            Hs.n,
            Float32(sum(spectra[:intensity_array][scan_idx])),
            scan_idx;
            min_spectral_contrast = params.min_spectral_contrast,
            min_log2_matched_ratio = params.min_log2_matched_ratio,
            min_y_count = params.min_y_count,
            min_frag_count = params.min_frag_count,
            max_best_rank = params.max_best_rank,
            min_topn = first(params.min_topn_of_m),
            block_size = 500000
        )

        # Reset arrays
        reset_arrays!(search_data, Hs)
    end
    return DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val]))
end

#==========================================================
Helper Methods
==========================================================#

"""
Resize arrays when needed.
"""
function resize_arrays!(search_data::SearchDataStructures, weights::Vector{Float32})
    new_entries = getIdToCol(search_data).size - length(weights) + 1000
    resize!(weights, length(weights) + new_entries)
    resize!(getComplexSpectralScores(search_data), length(getComplexSpectralScores(search_data)) + new_entries)
    append!(getComplexUnscoredPsms(search_data), [eltype(getComplexUnscoredPsms(search_data))() for _ in 1:new_entries])
end

"""
Initialize weights for deconvolution.
"""
function initialize_weights!(
    search_data::SearchDataStructures,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32}
)
    for i in 1:getIdToCol(search_data).size
        weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
            precursor_weights[getIdToCol(search_data).keys[i]]
    end
end

"""
Update precursor weights after deconvolution.
"""
function update_precursor_weights!(
    search_data::SearchDataStructures,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32}
)
    for i in 1:getIdToCol(search_data).size
        id = getIdToCol(search_data).keys[i]
        colid = getIdToCol(search_data)[id]
        precursor_weights[id] = weights[colid]
    end
end

"""
Reset arrays between scans.
"""
function reset_arrays!(search_data::SearchDataStructures, Hs::SparseArray)
    for i in 1:Hs.n
        getComplexUnscoredPsms(search_data)[i] = eltype(getComplexUnscoredPsms(search_data))()
    end
    reset!(getIdToCol(search_data))
    reset!(Hs)
end

"""
Reset results containers.
"""
function reset_results!(results::SecondPassSearchResults)
    results.psms[] = DataFrame()
end

function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    #best_psms = samplePSMsForXgboost(quant_psms_folder, params_[:xgboost_params]["max_n_samples"]);
   
    return nothing
end