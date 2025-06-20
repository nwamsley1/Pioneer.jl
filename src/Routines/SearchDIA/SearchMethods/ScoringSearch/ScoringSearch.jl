"""
    ScoringSearch

Search method for post-processing second pass results to get final protein scores.
This includes XGBoost model training, trace scoring, and protein group analysis.
"""
struct ScoringSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

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
        score_precursor_isotope_traces(
            second_pass_folder,
            getSecondPassPsms(getMSData(search_context)),
            getPrecursors(getSpecLib(search_context)),
            params.match_between_runs,
            params.max_q_value_xgboost_rescore,
            params.max_q_value_xgboost_mbr_rescore,
            params.max_psms_in_memory
        )

        # Step 2: Find Best Isotope Traces
        @info "Finding best traces..."
        best_traces = get_best_traces(
            getSecondPassPsms(getMSData(search_context)),
            params.min_best_trace_prob
        )

        # Create references for second pass PSMs
        second_pass_paths = getSecondPassPsms(getMSData(search_context))
        second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

        # Step 3: Process Quantification Results
        @info "Processing quantification results..."
        
        # Define columns we need to keep for quantification
        necessary_cols = get_quant_necessary_columns()
        
        # Build explicit pipeline - each operation is clear and testable
        quant_processing_pipeline = TransformPipeline() |>
            add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
            rename_column(:prob, :trace_prob) |>
            select_columns(vcat(necessary_cols, :best_trace)) |>
            filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
            remove_columns(:best_trace) |>
            sort_by([:global_prob, :target], rev=[true, true])
        
        # Apply pipeline to all second pass files
        @info "Applying quantification processing pipeline to $(length(second_pass_refs)) files"
        for ref in second_pass_refs
            if exists(ref)
                apply_pipeline!(ref, quant_processing_pipeline)
            end
        end
        
        # Use the processed references for downstream operations
        filtered_refs = second_pass_refs

        # Step 4: Merge PSM Scores by max_prob
        @info "Merging PSM scores for global q-value estimation..."
        merge_psm_files(
            filtered_refs,
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
        # Re-sort filtered files by prec_prob using references
        for ref in filtered_refs
            sort_file_by_keys!(ref, :prec_prob; reverse=true)
        end
        
        # Merge using new abstraction
        merge_psm_files(
            filtered_refs,
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

        # Step 8: Filter PSMs by global q-value
        @info "Filtering passing PSMs by global q-value and experiment-wide q-value threshold..."
        # Apply q-value threshold using new abstractions
        passing_refs = filter_psms_by_qvalue(
            filtered_refs,
            passing_psms_folder,
            getPrecursors(getSpecLib(search_context)),
            results.precursor_global_qval_interp[],
            results.precursor_qval_interp[],
            params.q_value_threshold
        )
        
        # Update the search context with new passing PSM paths
        passing_psm_paths = [file_path(ref) for ref in passing_refs]
        # Store the paths in the search context
        for (idx, path) in enumerate(passing_psm_paths)
            setPassingPsms!(getMSData(search_context), idx, path)
        end

        # Step 10: Count protein peptides
        # This is useful for counting some of the protein-group features
        # for target-decoy discriminatin (e.g. the fraction of observed vs. possible peptides
        #in the library for each protein group)
        @info "Counting protein peptides..."
        protein_to_possible_peptides = count_protein_peptides(
            getPrecursors(getSpecLib(search_context))
        )

        # Step 11: Perform protein inference and initial scoring
        @info "Performing protein inference and initial scoring..."
        # For now, still use the path-based function but with our reference paths
        psm_to_pg_path = perform_protein_inference(
            passing_psm_paths,  # Use the paths from our references
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            getPrecursors(getSpecLib(search_context)),
            protein_to_possible_peptides,
            min_peptides = params.min_peptides
        )

        # Create paired file references for internal use only
        @info "Creating file references for internal processing..."
        paired_files = PairedSearchFiles[]
        for (ms_file_idx, ref) in enumerate(passing_refs)
            psm_path = file_path(ref)
            if haskey(psm_to_pg_path, psm_path)
                pg_path = psm_to_pg_path[psm_path]
                # Store the protein group path in the search context
                setPassingProteins!(getMSData(search_context), ms_file_idx, pg_path)
                # Create paired reference for internal use
                push!(paired_files, PairedSearchFiles(psm_path, pg_path, ms_file_idx))
            end
        end
        
        if isempty(paired_files)
            error("No protein groups created during protein inference")
        end
        @info "Created $(length(paired_files)) paired references for internal processing"
        
        # Create internal scoring refs (not stored in SearchContext)
        scoring_refs = ScoringSearchResultRefs(paired_files)

        # Step 12: Perform protein probit regression
        @info "Performing protein probit regression..."
        qc_folder = joinpath(dirname(temp_folder), "qc_plots")
        !isdir(qc_folder) && mkdir(qc_folder)
        # Use references from scoring_refs
        pg_refs = get_protein_refs(scoring_refs)
        perform_protein_probit_regression(
            pg_refs,  # Now uses references
            params.max_psms_in_memory,
            qc_folder
        )

        # Step 18: Merge protein groups by experiment-wide pg_score
        # Merge protein groups by global pg_score (which contains probit scores after regression)
        sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
        @info "Merging protein group scores for experiment-wide q-value estimation..."
        
        # Use the protein group references we already have from scoring_refs
        # pg_refs already defined above
        
        # Use reference-based merge
        merge_protein_groups_by_score(pg_refs, sorted_pg_scores_path, batch_size=1000000)
        @info "Merged $(length(pg_refs)) protein group files using reference-based approach"

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
        # Use the protein group references directly
        acc_to_max_pg_score = calculate_and_add_global_scores!(pg_refs)

        # Merge protein groups by run-specific prob
        #if isfile(sorted_pg_scores_path)
        

            # Clear the file by writing an empty DataFrame
        #    temp_ref = ProteinGroupFileReference(sorted_pg_scores_path)
        #    write_arrow_file(temp_ref, DataFrame())
        #end

        # Step 17: Sort protein groups by experiment-wide pg_score
        @info "Sorting protein group tables by experiment-wide pg_score..."
        # Use references directly
        for ref in pg_refs
            sort_file_by_keys!(ref, :global_pg_score; reverse=true)
        end

        # Step 15: Merge protein group scores by global_pg_score
        @info "Merging protein group scores for global q-value estimation..."
        # Merge protein groups by global pg_score (which contains probit scores after regression)
        sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
        
        # Reuse the protein group references we already have
        pg_refs_global = pg_refs  # Already have these from earlier
        
        # Use reference-based merge for global scores (sorted by global_pg_score)
        stream_sorted_merge(pg_refs_global, sorted_pg_scores_path, :global_pg_score;
                           batch_size=1000000, reverse=true)
        @info "Merged $(length(pg_refs_global)) protein group files by global score"

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
        # Use the protein group references we already have
        get_proteins_passing_qval_refs(
            pg_refs,
            search_context.global_pg_score_to_qval[],
            search_context.pg_score_to_qval[],
            :global_pg_score,
            :pg_score,
            :global_pg_qval,
            :pg_qval,
            params.q_value_threshold
        )

        @info "DEBUG"
        test_pgs = Arrow.Table(readdir(passing_proteins_folder, join =true))
        ew_count = (test_pgs[:target].==true).&(test_pgs[:pg_qval].<=params.q_value_threshold) |> sum #|> @info "DEBUG: Number of passing protein groups: "
        global_count = (test_pgs[:target].==true).&(test_pgs[:global_pg_qval].<=params.q_value_threshold) |> sum #|> @info "DEBUG: Number of passing protein groups (global): "
        ew_global_count = (test_pgs[:target].==true).&(test_pgs[:global_pg_qval].<=params.q_value_threshold).&(test_pgs[:pg_qval].<=params.q_value_threshold) |> sum #|> @info "DEBUG: Number of passing protein groups (global and pg): " 
        @info "DEBUG: Number of passing protein groups: $ew_count"
        @info "DEBUG: Number of passing protein groups (global): $global_count"
        @info "DEBUG: Number of passing protein groups (global and pg): $ew_global_count"
        # Step 14: Update PSMs with probit-scored pg_score and global scores
        @info "Updating PSMs with probit-scored protein group scores..."
        # Get the paired references from SearchContext
        update_psms_with_probit_scores_refs(
            scoring_refs.paired_files,
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