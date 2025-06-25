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

#==========================================================
Q-value Spline Wrapper Functions
==========================================================#

"""
Create global precursor q-value spline (unique precursors only).
"""
function get_precursor_global_qval_spline(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    return get_qvalue_spline(
        merged_path, :global_prob, true;
        min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

"""
Create experiment-wide precursor q-value spline (all precursors).
"""
function get_precursor_qval_spline(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    return get_qvalue_spline(
        merged_path, :prec_prob, false;
        min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

"""
Create global protein q-value spline (unique protein groups only).
"""
function get_protein_global_qval_spline(merged_path::String, params::ScoringSearchParameters)
    return get_qvalue_spline(
        merged_path, :global_pg_score, true;
        min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
    )
end

"""
Create experiment-wide protein q-value spline (all protein groups).
"""
function get_protein_qval_spline(merged_path::String, params::ScoringSearchParameters)
    return get_qvalue_spline(
        merged_path, :pg_score, false;
        min_pep_points_per_bin = params.pg_q_value_interpolation_points_per_bin
    )
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
    
    # Set up output folders
    second_pass_folder = joinpath(temp_folder, "second_pass_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    
    for folder in [passing_psms_folder, passing_proteins_folder]
        !isdir(folder) && mkdir(folder)
    end

    try
        # Step 1: Train XGBoost Models
        @info "Step 1: Training XGBoost models..."
        step1_time = @elapsed begin
            score_precursor_isotope_traces(
                second_pass_folder,
                getSecondPassPsms(getMSData(search_context)),
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs,
                params.max_q_value_xgboost_rescore,
                params.max_q_value_xgboost_mbr_rescore,
                params.max_psms_in_memory
            )
        end
        @info "Step 1 completed in $(round(step1_time, digits=2)) seconds"

        # Step 2: Find Best Isotope Traces
        @info "Step 2: Finding best isotope traces..."
        step2_time = @elapsed begin
            best_traces = get_best_traces(
                getSecondPassPsms(getMSData(search_context)),
                params.min_best_trace_prob
            )
        end
        @info "Step 2 completed in $(round(step2_time, digits=2)) seconds"

        # Create references for second pass PSMs
        second_pass_paths = getSecondPassPsms(getMSData(search_context))
        second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

        # Step 3: Process Quantification Results
        @info "Step 3: Processing quantification results..."
        step3_time = @elapsed begin
            necessary_cols = get_quant_necessary_columns()
            
            quant_processing_pipeline = TransformPipeline() |>
                add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
                rename_column(:prob, :trace_prob) |>
                select_columns(vcat(necessary_cols, :best_trace)) |>
                filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
                remove_columns(:best_trace) |>
                sort_by([:global_prob, :target], rev=[true, true])
            
            apply_pipeline!(second_pass_refs, quant_processing_pipeline)
            
            filtered_refs = second_pass_refs
        end
        @info "Step 3 completed in $(round(step3_time, digits=2)) seconds"

        # Step 4: Merge PSMs by global_prob for global q-values
        @info "Step 4: Merging PSM scores by global_prob..."
        step4_time = @elapsed begin
            stream_sorted_merge(filtered_refs, results.merged_quant_path, :global_prob, :target;
                               batch_size=10_000_000, reverse=[true,true])
        end
        @info "Step 4 completed in $(round(step4_time, digits=2)) seconds"

        # Step 5: Calculate global precursor q-values
        @info "Step 5: Calculating global precursor q-values..."
        step5_time = @elapsed begin
            results.precursor_global_qval_interp[] = get_precursor_global_qval_spline(results.merged_quant_path, params, search_context)
        end
        @info "Step 5 completed in $(round(step5_time, digits=2)) seconds"

        # Step 6: Merge PSMs by prec_prob for experiment-wide q-values
        @info "Step 6: Re-sorting and merging PSMs by prec_prob..."
        step6_time = @elapsed begin
            sort_file_by_keys!(filtered_refs, :prec_prob, :target; reverse=[true,true])
            stream_sorted_merge(filtered_refs, results.merged_quant_path, :prec_prob, :target;
                               batch_size=10_000_000, reverse=[true,true])
        end
        @info "Step 6 completed in $(round(step6_time, digits=2)) seconds"

        # Step 7: Calculate experiment-wide precursor q-values
        @info "Step 7: Calculating experiment-wide precursor q-values..."
        step7_time = @elapsed begin
            results.precursor_qval_interp[] = get_precursor_qval_spline(results.merged_quant_path, params, search_context)
        end
        @info "Step 7 completed in $(round(step7_time, digits=2)) seconds"

        # Step 8: Filter PSMs by q-value thresholds
        @info "Step 8: Filtering PSMs by q-value thresholds..."
        step8_time = @elapsed begin
            qvalue_filter_pipeline = TransformPipeline() |>
                add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
                add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
                filter_by_multiple_thresholds([
                    (:global_qval, params.q_value_threshold),
                    (:qval, params.q_value_threshold)
                ])
            
            passing_refs = apply_pipeline_batch(
                filtered_refs,
                qvalue_filter_pipeline,
                passing_psms_folder
            )
        end
        @info "Step 8 completed in $(round(step8_time, digits=2)) seconds"
        
        # Update search context with passing PSM paths
        for (idx, ref) in enumerate(passing_refs)
            setPassingPsms!(getMSData(search_context), idx, file_path(ref))
        end

        # Step 9: Count protein peptides
        @info "Step 9: Counting protein peptides for feature calculation..."
        step9_time = @elapsed begin
            protein_to_possible_peptides = count_protein_peptides(
                getPrecursors(getSpecLib(search_context))
            )
        end
        @info "Step 9 completed in $(round(step9_time, digits=2)) seconds"

        # Step 10: Perform protein inference and initial scoring
        @info "Step 10: Performing protein inference and initial scoring..."
        step10_time = @elapsed begin
            pg_refs, psm_to_pg_mapping = perform_protein_inference_pipeline(
                passing_refs,
                passing_proteins_folder,
                getPrecursors(getSpecLib(search_context)),
                protein_to_possible_peptides,
                min_peptides = params.min_peptides
            )
            
            paired_files = [PairedSearchFiles(psm_path, pg_path, idx) 
                           for (idx, (psm_path, pg_path)) in enumerate(psm_to_pg_mapping)]
            
            isempty(paired_files) && error("No protein groups created during protein inference")
        end
        @info "Step 10 completed in $(round(step10_time, digits=2)) seconds"

        # Step 11: Perform protein probit regression
        @info "Step 11: Performing protein probit regression..."
        step11_time = @elapsed begin
            qc_folder = joinpath(dirname(temp_folder), "qc_plots")
            !isdir(qc_folder) && mkdir(qc_folder)
            
            perform_protein_probit_regression(
                pg_refs,
                params.max_psms_in_memory,
                qc_folder
            )
        end
        @info "Step 11 completed in $(round(step11_time, digits=2)) seconds"

        # Step 12: Calculate global protein scores
        @info "Step 12: Calculating global protein scores..."
        step12_time = @elapsed begin
            acc_to_max_pg_score = calculate_and_add_global_scores!(pg_refs)
        end
        @info "Step 12 completed in $(round(step12_time, digits=2)) seconds"

        # Step 13: Sort protein groups by global_pg_score
        @info "Step 13: Sorting protein groups by global_pg_score..."
        step13_time = @elapsed begin
            sort_file_by_keys!(pg_refs, :global_pg_score, :target; reverse=[true, true])
        end
        @info "Step 13 completed in $(round(step13_time, digits=2)) seconds"
        
        # Step 14: Merge protein groups for global q-values
        @info "Step 14: Merging protein groups for global q-value calculation..."
        step14_time = @elapsed begin
            sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
            stream_sorted_merge(pg_refs, sorted_pg_scores_path, :global_pg_score, :target;
                               batch_size=1000000, reverse=[true,true])
        end
        @info "Step 14 completed in $(round(step14_time, digits=2)) seconds"

        # Step 15: Calculate global protein q-values
        @info "Step 15: Calculating global protein q-values..."
        step15_time = @elapsed begin
            search_context.global_pg_score_to_qval[] = get_protein_global_qval_spline(sorted_pg_scores_path, params)
        end
        @info "Step 15 completed in $(round(step15_time, digits=2)) seconds"

        # Step 16: Sort protein groups by pg_score
        @info "Step 16: Sorting protein groups by pg_score..."
        step16_time = @elapsed begin
            sort_file_by_keys!(pg_refs, :pg_score, :target; reverse=[true, true])
        end
        @info "Step 16 completed in $(round(step16_time, digits=2)) seconds"
        
        # Step 17: Merge protein groups for experiment-wide q-values
        @info "Step 17: Merging protein groups for experiment-wide q-value calculation..."
        step17_time = @elapsed begin
            stream_sorted_merge(pg_refs, sorted_pg_scores_path, :pg_score, :target;
                               batch_size=1000000, reverse=[true,true])
        end
        @info "Step 17 completed in $(round(step17_time, digits=2)) seconds"

        # Step 18: Calculate experiment-wide protein q-values
        @info "Step 18: Calculating experiment-wide protein q-values..."
        step18_time = @elapsed begin
            search_context.pg_score_to_qval[] = get_protein_qval_spline(sorted_pg_scores_path, params)
        end
        @info "Step 18 completed in $(round(step18_time, digits=2)) seconds"

        # Step 19: Add q-values to protein groups
        @info "Step 19: Adding q-values and passing flags to protein groups..."
        step19_time = @elapsed begin
            protein_qval_pipeline = TransformPipeline() |>
                add_interpolated_column(:global_pg_qval, :global_pg_score, search_context.global_pg_score_to_qval[]) |>
                add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
                add_column(:passes_qval, df -> 
                    (df.global_pg_qval .<= params.q_value_threshold) .& 
                    (df.pg_qval .<= params.q_value_threshold))

            apply_pipeline!(pg_refs, protein_qval_pipeline)
        end
        @info "Step 19 completed in $(round(step19_time, digits=2)) seconds"

        # Step 20: Update PSMs with final protein scores
        @info "Step 20: Updating PSMs with final protein scores..."
        step20_time = @elapsed begin
            update_psms_with_probit_scores_refs(
                paired_files,
                acc_to_max_pg_score,
                search_context.pg_score_to_qval[],
                search_context.global_pg_score_to_qval[]
            )
        end
        @info "Step 20 completed in $(round(step20_time, digits=2)) seconds"

        # Summary of all step times
        total_time = step1_time + step2_time + step3_time + step4_time + step5_time + 
                    step6_time + step7_time + step8_time + step9_time + step10_time + 
                    step11_time + step12_time + step13_time + step14_time + step15_time + 
                    step16_time + step17_time + step18_time + step19_time + step20_time
        @info "ScoringSearch completed - Total time: $(round(total_time, digits=2)) seconds"
        @info "Breakdown: XGBoost($(round(step1_time, digits=1))s) + Best_Traces($(round(step2_time, digits=1))s) + Quant_Processing($(round(step3_time, digits=1))s) + Merging($(round(step4_time + step6_time + step14_time + step17_time, digits=1))s) + Q-values($(round(step5_time + step7_time + step15_time + step18_time, digits=1))s) + Protein_Inference($(round(step10_time, digits=1))s) + Probit_Regression($(round(step11_time, digits=1))s) + Other($(round(step3_time + step8_time + step9_time + step12_time + step13_time + step16_time + step19_time + step20_time, digits=1))s)"

        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end