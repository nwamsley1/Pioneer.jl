"""
    ScoringSearch

Search method for post-processing second pass results to get final protein scores.
This includes XGBoost model training, trace scoring, and protein group analysis.
"""
struct ScoringSearch <: SearchMethod end

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for scoring search.
"""
struct ScoringSearchResults <: SearchResults
    # Paths to results
    best_traces::Dict{Int64, Float32}
    precursor_pep_spline::Base.Ref{<:UniformSpline}  # Spline for PEP calculation
    precursor_qval_interp::Base.Ref{Any} # Interpolation for q-values
    pg_qval_interp::Base.Ref{Any}       # Protein group q-value interpolation
    merged_quant_path::String # Path to merged quantification results
end

"""
Parameters for scoring search.
"""
struct ScoringSearchParameters <: SearchParameters
    # XGBoost parameters
    max_n_samples::Int64
    min_best_trace_prob::Float32
    precursor_prob_spline_points_per_bin::Int64
    precursor_q_value_interpolation_points_per_bin::Int64
    pg_q_value_interpolation_points_per_bin::Int64

    function ScoringSearchParameters(params::Any)
        xp = params[:xgboost_params]
        
        new(
            Int64(xp["max_n_samples"]),
            Float32(xp["min_best_trace_prob"]),
            Int64(xp["precursor_prob_spline_points_per_bin"]),
            Int64(xp["precursor_q_value_interpolation_points_per_bin"]),
            Int64(xp["pg_q_value_interpolation_points_per_bin"])
        )
    end
end

#==========================================================
Interface Implementation  
==========================================================#

get_parameters(::ScoringSearch, params::Any) = ScoringSearchParameters(params)

function init_search_results(::ScoringSearchParameters, search_context::SearchContext)
    return ScoringSearchResults(
        Dict{Int64, Float32}(),  # best_traces
        Ref{UniformSpline}(),  # precursor_pep_spline
        Ref(undef),  # precursor_qval_interp
        Ref(undef),  # pg_qval_interp
        joinpath(getDataOutDir(search_context), "merged_quant.arrow")
    )
end

function process_file!(
    results::ScoringSearchResults,
    params::ScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
)
    # No per-file processing needed
    return results
end

function process_search_results!(
    results::ScoringSearchResults,
    params::ScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
)
    # No per-file results processing needed
    return nothing
end

function reset_results!(results::ScoringSearchResults)
    return nothing
end

"""
Process all results to get final protein scores.
"""
function summarize_results!(
    results::ScoringSearchResults,
    params::ScoringSearchParameters,
    search_context::SearchContext
)
    temp_folder = getDataOutDir(search_context)
    
    # Create necessary folders
    second_pass_folder = joinpath(temp_folder, "second_pass_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    
    for folder in [passing_psms_folder, passing_proteins_folder]
        !isdir(folder) && mkdir(folder)
    end

    try
        # Train XGBoost models
        @info "Training XGBoost models..."
        best_psms = samplePSMsForXgboost(second_pass_folder, params.max_n_samples)
        models = scoreTraces!(
            best_psms,
            getSecondPassPsms(getMSData(search_context)),
            getPrecursors(getSpecLib(search_context))
        )
        best_psms = nothing
        GC.gc()

        # Get best traces
        @info "Finding best traces..."
        best_traces = getBestTraces(
            getSecondPassPsms(getMSData(search_context)),
            params.min_best_trace_prob
        )

        # Sort and filter quantification tables
        @info "Processing quantification results..."
        sortAndFilterQuantTables(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path,
            best_traces
        )

        # Merge scores
        @info "Merging PSM scores..."
        mergeSortedPSMScores(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path
        )

        # Calculate error probabilities
        @info "Calculating error probabilities..."
        results.precursor_pep_spline[] = getPEPSpline(
            results.merged_quant_path,
            :prob,
            min_pep_points_per_bin = params.precursor_prob_spline_points_per_bin,
            n_spline_bins = 5
        )

        results.precursor_qval_interp[] = getQValueSpline(
            results.merged_quant_path,
            :prob,
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin
        )

        # Get passing PSMs
        @info "Filtering passing PSMs..."
        getPSMsPassingQVal(
            getPassingPsms(getMSData(search_context)),
            passing_psms_folder,
            getSecondPassPsms(getMSData(search_context)),
            results.precursor_pep_spline[],
            results.precursor_qval_interp[],
            0.01f0
        )

        # Score protein groups
        @info "Scoring protein groups..."
        sorted_pg_score_path = getProteinGroups(
            getPassingPsms(getMSData(search_context)),
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            temp_folder,
            getPrecursors(getSpecLib(search_context))
        )

        search_context.pg_score_to_qval[] = getQValueSpline(
            sorted_pg_score_path,
            :max_pg_score,
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
        )

        # Store best traces
        #results.best_traces = best_traces
        best_traces = nothing
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end