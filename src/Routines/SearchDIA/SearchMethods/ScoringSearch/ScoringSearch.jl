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
    max_psms_in_memory::Int64
    min_best_trace_prob::Float32
    precursor_prob_spline_points_per_bin::Int64
    precursor_q_value_interpolation_points_per_bin::Int64
    pg_prob_spline_points_per_bin::Int64  # Added based on original struct
    pg_q_value_interpolation_points_per_bin::Int64  # Added based on original struct
    match_between_runs::Bool
    min_peptides::Int64

    function ScoringSearchParameters(params::PioneerParameters)
        # Extract machine learning parameters from optimization section
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings 
        protein_inference_params = params.protein_inference
        
        new(
            Int64(ml_params.max_psms_in_memory),
            Float32(ml_params.min_trace_prob),
            Int64(ml_params.spline_points),
            Int64(ml_params.interpolation_points),
            Int64(ml_params.spline_points),        # Using same value for protein groups
            Int64(ml_params.interpolation_points), # Using same value for protein groups
            Bool(global_params.match_between_runs),
            Int64(protein_inference_params.min_peptides)
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
    spectra::MassSpecData
)
    # No per-file processing needed
    return results
end

function process_search_results!(
    results::ScoringSearchResults,
    params::ScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
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
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    
    # Set up output folders for different stages
    second_pass_folder = joinpath(temp_folder, "second_pass_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    
    for folder in [passing_psms_folder, passing_proteins_folder]
        !isdir(folder) && mkdir(folder)
    end

    try
        # Step 1: Train XGBoost Models
        @info "Training XGBoost models..."
        # Sample PSMs for training to avoid memory consumption issues
        psms_count = get_psms_count(getSecondPassPsms(getMSData(search_context)))

        if psms_count > params.max_psms_in_memory #Use out-of-memory algorithm
            #Sample psms for xgboost training. Only the sampled psms are stored in-memory
            best_psms = sample_psms_for_xgboost(second_pass_folder, psms_count, params.max_psms_in_memory)#params.max_n_samples)
            models = score_precursor_isotope_traces_out_of_memory!(
                best_psms,
                getSecondPassPsms(getMSData(search_context)),
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs
            )
        else #In memory algorithm
            best_psms = load_psms_for_xgboost(second_pass_folder)#params.max_n_samples)
            models = score_precursor_isotope_traces_in_memory!(
                best_psms,
                getSecondPassPsms(getMSData(search_context)),
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs
            )
        end
        best_psms = nothing
        GC.gc()

        # Step 2: Find Best Isotope Traces
        @info "Finding best traces..."
        best_traces = get_best_traces(
            getSecondPassPsms(getMSData(search_context)),
            params.min_best_trace_prob
        )

        # Step 3: Process Quantification Results
        @info "Processing quantification results..."
        # Filter to best traces and sort tables
        sort_and_filter_quant_tables(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path,
            best_traces
        )

        # Step 4: Merge PSM Scores
        @info "Merging PSM scores..."
        merge_sorted_psms_scores(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path
        )

        # Step 5: Calculate Error Probabilities
        @info "Calculating error probabilities..."
        # Create PEP spline
        results.precursor_pep_spline[] = get_pep_spline(
            results.merged_quant_path,
            :prob,
            min_pep_points_per_bin = params.precursor_prob_spline_points_per_bin,
            n_spline_bins = 5
        )
        # Create q-value interpolation
        results.precursor_qval_interp[] = get_qvalue_spline(
            results.merged_quant_path,
            :prob,
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin
        )

        # Step 6: Filter PSMs
        @info "Filtering passing PSMs..."
        # Apply q-value threshold and store passing PSMs
        get_psms_passing_qval(
            getPassingPsms(getMSData(search_context)),
            passing_psms_folder,
            getSecondPassPsms(getMSData(search_context)),
            results.precursor_pep_spline[],
            results.precursor_qval_interp[],
            0.01f0
        )
        # Step 7: Score Protein Groups
        @info "Scoring protein groups..."
        # Create protein groups and calculate scores
        sorted_pg_score_path, protein_inference_dict = get_protein_groups(
            getPassingPsms(getMSData(search_context)),
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            temp_folder,
            getPrecursors(getSpecLib(search_context)),
            min_peptides = params.min_peptides
        )

        add_protein_inferrence_col(
            getPassingPsms(getMSData(search_context)),
            protein_inference_dict,
            getSequence(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context)))
        )
        # Create protein group q-value interpolation
        search_context.pg_score_to_qval[] = get_qvalue_spline(
            sorted_pg_score_path,
            :pg_score,
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
        )


        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end