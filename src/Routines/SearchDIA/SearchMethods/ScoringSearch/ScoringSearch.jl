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
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin
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
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin
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

        # Step 10: Score Protein Groups
        @info "Scoring protein groups..."
        # Create protein groups and calculate scores
        protein_inference_dict = get_protein_groups(
            getPassingPsms(getMSData(search_context)),
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            temp_folder,
            getPrecursors(getSpecLib(search_context)),
            min_peptides = params.min_peptides,
            max_psms_in_memory = params.max_psms_in_memory
        )

        add_protein_inference_col(
            getPassingPsms(getMSData(search_context)),
            protein_inference_dict,
            getSequence(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getEntrapmentGroupId(getPrecursors(getSpecLib(search_context)))
        )

        # Step 11: Merge PSM Scores by max_prob
        @info "Merging protein group scores for global q-value estimation..."
        # Merge protein groups by global prob
        sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
        merge_sorted_protein_groups(
            passing_proteins_folder,
            sorted_pg_scores_path,
            :global_pg_score,
            N = 1000000
        )

        # Step 12: Create q-value interpolation
        @info "Calculating global q-values for protein groups..."
        # Create protein group global q-value interpolation
        search_context.global_pg_score_to_qval[] = get_qvalue_spline(
            sorted_pg_scores_path,
            :global_pg_score,
            true; # use unique protein groups
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
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

        # Step 12: Create q-value interpolation
        @info "Sorting protein group tables by experiment-wide pg_score..."
        sort_protein_tables(
            getPassingProteins(getMSData(search_context)),
            passing_proteins_folder,
            :pg_score
        )
        # Step 13: Merge PSM Scores by prob
        @info "Merging protein group scores for experiment-wide q-value estimation..."
        merge_sorted_protein_groups(
            passing_proteins_folder,
            sorted_pg_scores_path,
            :pg_score,
            N = 1000000
        )

        # Step 14: Create q-value interpolation
        @info "Calculating experiment-wide q-values for protein groups..."
        # Create protein group run-specific q-value interpolation
        search_context.pg_score_to_qval[] = get_qvalue_spline(
            sorted_pg_scores_path,
            :pg_score,
            false; # use all protein groups
            min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
        )

        # Load and summarize protein groups before filtering
        @info "Loading protein groups for summary statistics..."
        protein_files = [path for path in readdir(passing_proteins_folder, join=true) if endswith(path, ".arrow")]
        all_proteins = DataFrame()
        for file_path in protein_files
            append!(all_proteins, DataFrame(Arrow.Table(file_path)))
        end
        
        # Display summary statistics
        #@info "Protein group statistics before filtering:"
        #@info "Total protein groups: $(nrow(all_proteins))"
        #@info "Target protein groups: $(sum(all_proteins.target))"
        #@info "Decoy protein groups: $(sum(.!all_proteins.target))"
        
        # Column-wise statistics
        for col in names(all_proteins)
            col_data = all_proteins[!, col]
            if eltype(col_data) <: Number && !all(ismissing.(col_data))
                non_missing = skipmissing(col_data)
                #@info "Column '$col': min=$(minimum(non_missing)), max=$(maximum(non_missing)), mean=$(round(mean(non_missing), digits=4))"
            end
        end
        
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



        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end