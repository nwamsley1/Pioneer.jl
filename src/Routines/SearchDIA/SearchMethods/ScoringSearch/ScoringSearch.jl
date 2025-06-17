"""
    ScoringSearch

Search method for post-processing second pass results to get final protein scores.
This includes XGBoost model training, trace scoring, and protein group analysis.
"""
struct ScoringSearch <: SearchMethod end

# Export new types for use in other modules
export ProteinKey, PeptideKey, ProteinFeatures, ProteinGroup, 
       ProteinGroupBuilder, InferenceResult, FileMapping,
       add_peptide!, finalize, to_protein_key, to_namedtuple, to_peptide_key

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
    precursor_global_qval_interp::Base.Ref{Any} # Interpolation for global q-values
    precursor_qval_interp::Base.Ref{Any} # Interpolation for run-specific q-values
    pg_qval_interp::Base.Ref{Any}       # Protein group q-value interpolation
    merged_quant_path::String # Path to merged quantification results
end

"""
Parameters for scoring search.
"""
struct ScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    # XGBoost parameters
    max_psms_in_memory::Int64
    min_best_trace_prob::Float32
    precursor_prob_spline_points_per_bin::Int64
    precursor_q_value_interpolation_points_per_bin::Int64
    pg_prob_spline_points_per_bin::Int64  # Added based on original struct
    pg_q_value_interpolation_points_per_bin::Int64  # Added based on original struct
    match_between_runs::Bool
    min_peptides::Int64
    max_q_value_xgboost_rescore::Float32
    max_q_value_xgboost_mbr_rescore::Float32
    q_value_threshold::Float32
    isotope_tracetype::I

    function ScoringSearchParameters(params::PioneerParameters)
        # Extract machine learning parameters from optimization section
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings 
        protein_inference_params = params.protein_inference

        # Determine isotope trace type
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) && global_params.isotope_settings.combine_traces
            CombineTraces(0.0f0)  # Default min_fraction_transmitted
        else
            SeperateTraces()
        end
        
        new{typeof(isotope_trace_type)}(
            Int64(ml_params.max_psms_in_memory),
            Float32(ml_params.min_trace_prob),
            Int64(ml_params.spline_points),
            Int64(ml_params.interpolation_points),
            Int64(ml_params.spline_points),        # Using same value for protein groups
            Int64(ml_params.interpolation_points), # Using same value for protein groups
            Bool(global_params.match_between_runs),
            Int64(protein_inference_params.min_peptides),
            Float32(ml_params.max_q_value_xgboost_rescore),
            Float32(ml_params.max_q_value_xgboost_mbr_rescore),
            Float32(global_params.scoring.q_value_threshold),
            isotope_trace_type
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
        Ref(undef),  # precursor_global_qval_interp
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
                params.match_between_runs,
                params.max_q_value_xgboost_rescore,
                params.max_q_value_xgboost_mbr_rescore
            )
        else #In memory algorithm
            best_psms = load_psms_for_xgboost(second_pass_folder)#params.max_n_samples)
            models = score_precursor_isotope_traces_in_memory!(
                best_psms,
                getSecondPassPsms(getMSData(search_context)),
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs,
                params.max_q_value_xgboost_rescore,
                params.max_q_value_xgboost_mbr_rescore
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
            params.isotope_tracetype,
            :global_prob,
            best_traces
        )

        # Step 4: Merge PSM Scores by max_prob
        @info "Merging PSM scores for global q-value estimation..."
        merge_sorted_psms_scores(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path,
            :global_prob
        )
        # Step 5: Create q-value interpolation
        @info "Calculating global precursor q-values..."
        # Create PEP spline
        #results.precursor_pep_spline[] = get_pep_spline(
        #    results.merged_quant_path,
        #    :global_prob,
        #    min_pep_points_per_bin = params.precursor_prob_spline_points_per_bin,
        #    n_spline_bins = 5
        #)
        results.precursor_global_qval_interp[] = get_qvalue_spline(
            results.merged_quant_path,
            :global_prob,
            true; # use unique precursors
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        )
        # Step 6: Merge PSM Scores by prob
        @info "Merging PSM scores for experiment-wide q-value estimation..."
        #Individual PSMs files need to be re-rosted by prec_prob 
        sort_quant_tables(
            getSecondPassPsms(getMSData(search_context)),
            results.merged_quant_path,
            :prec_prob
        )
        merge_sorted_psms_scores(
             getSecondPassPsms(getMSData(search_context)),
             results.merged_quant_path,
             :prec_prob
        )

        # Step 7: Create q-value interpolation
        @info "Calculating experiment-wide precursor q-values..."
        # Create PEP spline
        #results.precursor_pep_spline[] = get_pep_spline(
        #    results.merged_quant_path,
        #    :prob,
        #    min_pep_points_per_bin = params.precursor_prob_spline_points_per_bin,
        #    n_spline_bins = 5
        #)
        # Create q-value interpolation
        results.precursor_qval_interp[] = get_qvalue_spline(
            results.merged_quant_path,
            :prec_prob,
            false; # use all precursors
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        )

        # Step 6: Filter PSMs by global q-value
        @info "Filtering passing PSMs by global q-value and experiment-wide q-value threshold..."
        # Apply q-value threshold and store passing PSMs
        passing_psms_paths = getSecondPassPsms(getMSData(search_context))
        get_psms_passing_qval(
            getPrecursors(getSpecLib(search_context)),
            getPassingPsms(getMSData(search_context)),
            passing_psms_folder,
            passing_psms_paths,
            #results.precursor_pep_spline[],
            results.precursor_global_qval_interp[],
            results.precursor_qval_interp[],
            :global_prob,
            :prec_prob,
            :global_qval,
            :run_specific_qval,
            params.q_value_threshold,
        )

        # Step 10: Count protein peptides
        # This is useful for counting some of the protein-group features
        # for target-decoy discriminatin (e.g. the fraction of observed vs. possible peptides
        #in the library for each protein group)
        @info "Counting protein peptides..."
        protein_to_possible_peptides = count_protein_peptides(
            getPrecursors(getSpecLib(search_context))
        )

        # Step 11: Perform protein inference and initial scoring
        # Also links each psms table path to it's corresponding pg group table path. 
        @info "Performing protein inference and initial scoring..."
        psm_to_pg_path = perform_protein_inference(
            getPassingPsms(getMSData(search_context)),
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            getPrecursors(getSpecLib(search_context)),
            protein_to_possible_peptides,
            min_peptides = params.min_peptides
        )

        # Step 12: Perform protein probit regression
        @info "Performing protein probit regression..."
        qc_folder = joinpath(dirname(temp_folder), "qc_plots")
        !isdir(qc_folder) && mkdir(qc_folder)
        perform_protein_probit_regression(
            getPassingProteins(getMSData(search_context)),
            params.max_psms_in_memory,
            qc_folder
        )

        # Step 18: Merge protein groups by experiment-wide pg_score
        # Merge protein groups by global pg_score (which contains probit scores after regression)
        sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
        @info "Merging protein group scores for experiment-wide q-value estimation..."
        merge_sorted_protein_groups(
            passing_proteins_folder,
            sorted_pg_scores_path,
            :pg_score,
            N = 1000000
        )

        # Step 19: Create experiment-wide q-value interpolation
        @info "Calculating experiment-wide q-values for protein groups..."
        # Create protein group run-specific q-value interpolation
        search_context.pg_score_to_qval[] = get_qvalue_spline(
            sorted_pg_scores_path,
            :pg_score,
            false; # use all protein groups
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
        )

        # Step 13: Calculate global protein scores
        @info "Calculating global protein scores..."
        acc_to_max_pg_score = calculate_global_protein_scores(
            getPassingProteins(getMSData(search_context))
        )

        # Merge protein groups by run-specific prob
        if isfile(sorted_pg_scores_path)

            #rm(sorted_pg_scores_path)
            if Sys.iswindows()
                writeArrow(
                    sorted_pg_scores_path,
                    DataFrame()
                )
            else
                rm(sorted_pg_scores_path)
            end
        end

        # Step 17: Sort protein groups by experiment-wide pg_score
        @info "Sorting protein group tables by experiment-wide pg_score..."
        sort_protein_tables(
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            :global_pg_score
        )

        # Step 15: Merge protein group scores by global_pg_score
        @info "Merging protein group scores for global q-value estimation..."
        # Merge protein groups by global pg_score (which contains probit scores after regression)
        sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
        merge_sorted_protein_groups(
            passing_proteins_folder,
            sorted_pg_scores_path,
            :global_pg_score,
            N = 1000000
        )

        # Step 16: Create global q-value interpolation
        @info "Calculating global q-values for protein groups..."
        # Create protein group global q-value interpolation
        search_context.global_pg_score_to_qval[] = get_qvalue_spline(
            sorted_pg_scores_path,
            :global_pg_score,
            true; # use unique protein groups
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
        )

        # Filter proteins by global q-value
        get_proteins_passing_qval(
            passing_proteins_folder,
            search_context.global_pg_score_to_qval[],
            search_context.pg_score_to_qval[],
            :global_pg_score,
            :pg_score,
            :global_pg_qval,
            :pg_qval,
            params.q_value_threshold,
        )

        # Step 14: Update PSMs with probit-scored pg_score and global scores
        @info "Updating PSMs with probit-scored protein group scores..."
        update_psms_with_probit_scores(
            getPrecursors(getSpecLib(search_context)),
            psm_to_pg_path,
            acc_to_max_pg_score,
            search_context.pg_score_to_qval[],
            search_context.global_pg_score_to_qval[]
        )

        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end