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

"""
    ScoringSearch

Search method for post-processing second pass results to get final protein scores.
This includes LightGBM model training, trace scoring, and protein group analysis.
"""
struct ScoringSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for scoring search.
"""
struct ScoringSearchResults <: SearchResults
    # Paths to results
    best_traces::Dict{Int64, Float32}
    precursor_global_qval_interp::Base.Ref{Any} # Interpolation for global q-values
    precursor_qval_interp::Base.Ref{Any} # Interpolation for run-specific q-values
    precursor_pep_interp::Base.Ref{Any}  # Interpolation for experiment-wide PEPs
    pg_qval_interp::Base.Ref{Any}       # Protein group q-value interpolation
    merged_quant_path::String # Path to merged quantification results
end

"""
Parameters for scoring search.
"""
struct ScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    # LightGBM parameters
    max_psms_in_memory::Int64
    min_best_trace_prob::Float32
    precursor_prob_spline_points_per_bin::Int64
    precursor_q_value_interpolation_points_per_bin::Int64
    pg_prob_spline_points_per_bin::Int64  # Added based on original struct
    pg_q_value_interpolation_points_per_bin::Int64  # Added based on original struct
    match_between_runs::Bool
    min_peptides::Int64
    max_q_value_lightgbm_rescore::Float32
    max_q_value_mbr_itr::Float32
    min_PEP_neg_threshold_itr::Float32
    max_MBR_false_transfer_rate::Float32
    q_value_threshold::Float32
    isotope_tracetype::I

    # Quantile binning parameter
    n_quantile_bins::Int64

    # Model comparison parameters
    enable_model_comparison::Bool
    validation_split_ratio::Float64
    qvalue_threshold_comparison::Float64
    min_psms_for_comparison::Int64
    max_psms_for_comparison::Int64

    # MS1 scoring parameter
    ms1_scoring::Bool

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
            Float32(global_params.scoring.q_value_threshold),
            Float32(ml_params.max_q_value_mbr_itr),
            Float32(ml_params.min_PEP_neg_threshold_itr),
            Float32(global_params.scoring.q_value_threshold),
            Float32(global_params.scoring.q_value_threshold),
            isotope_trace_type,

            # Quantile binning parameter
            Int64(ml_params.n_quantile_bins),

            # Model comparison parameters with defaults
            Bool(get(ml_params, :enable_model_comparison, true)),
            Float64(get(ml_params, :validation_split_ratio, 0.2)),
            Float64(get(ml_params, :qvalue_threshold, 0.01)),
            Int64(get(ml_params, :min_psms_for_comparison, 1000)),
            Int64(get(ml_params, :max_psms_for_comparison, 100000)),

            # MS1 scoring parameter
            Bool(global_params.ms1_scoring)
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
        Ref(undef),  # precursor_global_qval_interp
        Ref(undef),  # precursor_qval_interp
        Ref(undef),  # precursor_pep_interp
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
    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "ScoringSearch")
        return results  # Return early with unchanged results
    end
    
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
    # Use MBR-boosted scores when MBR is enabled
    score_col = params.match_between_runs ? :MBR_boosted_global_prob : :global_prob
    return get_qvalue_spline(
        merged_path, score_col, true;
        min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

"""
Create experiment-wide precursor q-value spline (all precursors).
"""
function get_precursor_qval_spline(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    # Use MBR-boosted scores when MBR is enabled
    score_col = params.match_between_runs ? :MBR_boosted_prec_prob : :prec_prob
    return get_qvalue_spline(
        merged_path, score_col, false;
        min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

"""
Create experiment-wide precursor PEP interpolation (all precursors).
"""
function get_precursor_pep_interpolation(merged_path::String, params::ScoringSearchParameters, search_context::SearchContext)
    # Use MBR-boosted scores when MBR is enabled
    score_col = params.match_between_runs ? :MBR_boosted_prec_prob : :prec_prob
    return get_pep_interpolation(
        merged_path, score_col;
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context),
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
Create experiment-wide protein group PEP interpolation (all protein groups).
"""
function get_protein_pep_interpolation(merged_path::String, params::ScoringSearchParameters)
    return get_pep_interpolation(
        merged_path, :pg_score;
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
        # Step 1: Train LightGBM Models
        ##@debug_l1 "Step 1: Training LightGBM models..."
        # Filter to only include valid (non-failed) files
        valid_file_data = get_valid_file_paths(search_context, getSecondPassPsms)
        valid_file_indices = [idx for (idx, _) in valid_file_data]
        valid_second_pass_psms = [path for (_, path) in valid_file_data]
        
        # Check if any valid files remain
        if isempty(valid_second_pass_psms)
            @user_warn "No valid files for ScoringSearch - all files failed in previous search methods"
            return nothing
        end
        
        step1_time = @elapsed begin
            score_precursor_isotope_traces(
                second_pass_folder,
                valid_second_pass_psms,
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs,
                params.max_q_value_lightgbm_rescore,
                params.max_q_value_mbr_itr,
                params.min_PEP_neg_threshold_itr,
                params.max_psms_in_memory,
                params.n_quantile_bins,
                params.q_value_threshold,
                params.ms1_scoring
            )
        end
        #@debug_l1 "Step 1 completed in $(round(step1_time, digits=2)) seconds"

        # Create references for second pass PSMs (only valid files)
        second_pass_paths = valid_second_pass_psms
        second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

        # Step 2: Compute precursor probabilities and FTR filtering
        # Step 2: Compute precursor probabilities and FTR filtering
        #@debug_l1 "Step 2: Computing precursor probabilities..."
        step2_time = @elapsed begin
            # Merge all second pass PSMs for experiment-wide calculations
            merged_scores_path = joinpath(temp_folder, "merged_trace_scores.arrow")
            sort_file_by_keys!(second_pass_refs, :trace_prob, :target; reverse=[true, true])
            stream_sorted_merge(second_pass_refs, merged_scores_path, :trace_prob, :target;
                               reverse=[true, true])

            merged_df = DataFrame(Arrow.Table(merged_scores_path))
            sqrt_n_runs = floor(Int64, sqrt(length(getFilePaths(getMSData(search_context)))))

            if params.match_between_runs
                # Apply MBR filter to MBR_boosted_trace_prob column (modifies in place)
                apply_mbr_filter!(merged_df, params)

                # Aggregate MBR-boosted scores to precursor level
                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :MBR_boosted_trace_prob => (p -> begin
                               prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(prob)
                           end) => :MBR_boosted_prec_prob)

                transform!(groupby(merged_df, :precursor_idx),
                           :MBR_boosted_prec_prob => (p -> logodds(p, sqrt_n_runs)) => :MBR_boosted_global_prob)

                # Also aggregate non-MBR scores for protein inference (steps 11+)
                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :trace_prob => (p -> begin
                               prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(prob)
                           end) => :prec_prob)

                transform!(groupby(merged_df, :precursor_idx),
                           :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
            else
                # No MBR: only aggregate base probabilities
                transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
                           :trace_prob => (p -> begin
                               prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                               prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                               Float32(prob)
                           end) => :prec_prob)

                transform!(groupby(merged_df, :precursor_idx),
                           :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
            end
            # Write updated data back to individual files
            for (file_idx, ref) in zip(valid_file_indices, second_pass_refs)
                sub_df = merged_df[merged_df.ms_file_idx .== file_idx, :]
                write_arrow_file(ref, sub_df)
            end
        end
        #@debug_l1 "Step 2 completed in $(round(step2_time, digits=2)) seconds"

        # Step 3: Find Best Isotope Traces
        #@debug_l1 "Step 3: Finding best isotope traces..."
        step3_time = @elapsed begin
            best_traces = get_best_traces(
                valid_second_pass_psms,
                params.min_best_trace_prob
            )
        end
        #@debug_l1 "Step 3 completed in $(round(step3_time, digits=2)) seconds"

        # Step 4: Process Quantification Results
        #@debug_l1 "Step 4: Processing quantification results..."
        step4_time = @elapsed begin
            necessary_cols = get_quant_necessary_columns(params.match_between_runs)

            if params.match_between_runs
                quant_processing_pipeline = TransformPipeline() |>
                    add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
                    select_columns(vcat(necessary_cols, :best_trace)) |>
                    filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
                    remove_columns(:best_trace) |>
                    sort_by([:MBR_boosted_global_prob, :target], rev=[true, true])
            else
                quant_processing_pipeline = TransformPipeline() |>
                    add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
                    select_columns(vcat(necessary_cols, :best_trace)) |>
                    filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
                    remove_columns(:best_trace) |>
                    sort_by([:global_prob, :target], rev=[true, true])
            end

            apply_pipeline!(second_pass_refs, quant_processing_pipeline)
            filtered_refs = second_pass_refs
        end
        #@debug_l1 "Step 4 completed in $(round(step4_time, digits=2)) seconds"

        # Step 5: Merge PSMs by global_prob for global q-values
        #@debug_l1 "Step 5: Merging PSM scores by global_prob..."
        step5_time = @elapsed begin
            if params.match_between_runs
                stream_sorted_merge(filtered_refs, results.merged_quant_path, :MBR_boosted_global_prob, :target;
                                   batch_size=10_000_000, reverse=[true,true])
            else
                stream_sorted_merge(filtered_refs, results.merged_quant_path, :global_prob, :target;
                                   batch_size=10_000_000, reverse=[true,true])
            end
        end
        #@debug_l1 "Step 5 completed in $(round(step5_time, digits=2)) seconds"

        # Step 6: Calculate global precursor q-values
        #@debug_l1 "Step 6: Calculating global precursor q-values..."
        step6_time = @elapsed begin
            results.precursor_global_qval_interp[] = get_precursor_global_qval_spline(results.merged_quant_path, params, search_context)
        end
        #@debug_l1 "Step 6 completed in $(round(step6_time, digits=2)) seconds"

        # Step 7: Merge PSMs by prec_prob for experiment-wide q-values
        #@debug_l1 "Step 7: Re-sorting and merging PSMs by prec_prob..."
        step7_time = @elapsed begin
            if params.match_between_runs
                sort_file_by_keys!(filtered_refs, :MBR_boosted_prec_prob, :target; reverse=[true,true])
                stream_sorted_merge(filtered_refs, results.merged_quant_path, :MBR_boosted_prec_prob, :target;
                                   batch_size=10_000_000, reverse=[true,true])
            else
                sort_file_by_keys!(filtered_refs, :prec_prob, :target; reverse=[true,true])
                stream_sorted_merge(filtered_refs, results.merged_quant_path, :prec_prob, :target;
                                   batch_size=10_000_000, reverse=[true,true])
            end
        end
        #@debug_l1 "Step 7 completed in $(round(step7_time, digits=2)) seconds"

        # Step 8: Calculate experiment-wide precursor q-values and PEPs
        #@debug_l1 "Step 8: Calculating experiment-wide precursor q-values and PEPs..."
        step8_time = @elapsed begin
            qval_interp = get_precursor_qval_spline(results.merged_quant_path, params, search_context)
            results.precursor_qval_interp[] = qval_interp
            results.precursor_pep_interp[]  = get_precursor_pep_interpolation(results.merged_quant_path, params, search_context)
        end
        #@debug_l1 "Step 8 completed in $(round(step8_time, digits=2)) seconds"

        # Step 9: Filter PSMs by q-value thresholds
        #@debug_l1 "Step 9: Filtering PSMs by q-value thresholds..."
        step9_time = @elapsed begin
            if params.match_between_runs
                qvalue_filter_pipeline = TransformPipeline() |>
                    add_interpolated_column(:MBR_boosted_global_qval, :MBR_boosted_global_prob, results.precursor_global_qval_interp[]) |>
                    add_interpolated_column(:MBR_boosted_qval, :MBR_boosted_prec_prob, results.precursor_qval_interp[]) |>
                    add_interpolated_column(:pep, :MBR_boosted_prec_prob, results.precursor_pep_interp[]) |>
                    filter_by_multiple_thresholds([
                        (:MBR_boosted_global_qval, params.q_value_threshold),
                        (:MBR_boosted_qval, params.q_value_threshold)
                    ])
            else
                qvalue_filter_pipeline = TransformPipeline() |>
                    add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
                    add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
                    add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
                    filter_by_multiple_thresholds([
                        (:global_qval, params.q_value_threshold),
                        (:qval, params.q_value_threshold)
                    ])
            end

            passing_refs = apply_pipeline_batch(
                filtered_refs,
                qvalue_filter_pipeline,
                passing_psms_folder
            )
        end

        #@debug_l1 "Step 9 completed in $(round(step9_time, digits=2)) seconds"

        #@debug_l1 "Step 10: Re-calculate q-values using non-MBR scores..."
        step10_time = @elapsed begin
            if params.match_between_runs
                # Part A: Recalculate MBR_boosted_qval using non-MBR prec_prob
                # Sort by non-MBR prec_prob
                sort_file_by_keys!(passing_refs, :prec_prob, :target; reverse=[true,true])

                # Merge by non-MBR prec_prob
                stream_sorted_merge(passing_refs, results.merged_quant_path, :prec_prob, :target;
                                    batch_size=10_000_000, reverse=[true,true])

                # Calculate q-value spline using non-MBR precursor scores
                qval_interp_nonMBR = get_qvalue_spline(
                    results.merged_quant_path, :prec_prob, false;
                    min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
                    fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
                )

                # Recalculate MBR_boosted_qval column with non-MBR scores
                # (overwrites the MBR-boosted qval with unbiased qval for protein inference)
                recalculate_qvalue_pipeline = TransformPipeline() |>
                    add_interpolated_column(:MBR_boosted_qval, :prec_prob, qval_interp_nonMBR)

                passing_refs = apply_pipeline_batch(
                    passing_refs,
                    recalculate_qvalue_pipeline,
                    passing_psms_folder
                )

                # Part B: Recalculate MBR_boosted_global_qval using non-MBR global_prob
                # Sort by non-MBR global_prob
                sort_file_by_keys!(passing_refs, :global_prob, :target; reverse=[true,true])

                # Merge by non-MBR global_prob
                stream_sorted_merge(passing_refs, results.merged_quant_path, :global_prob, :target;
                                    batch_size=10_000_000, reverse=[true,true])

                # Calculate global q-value spline using non-MBR global scores
                global_qval_interp_nonMBR = get_qvalue_spline(
                    results.merged_quant_path, :global_prob, true;
                    min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
                    fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
                )

                # Recalculate MBR_boosted_global_qval column with non-MBR scores
                recalculate_global_qvalue_pipeline = TransformPipeline() |>
                    add_interpolated_column(:MBR_boosted_global_qval, :global_prob, global_qval_interp_nonMBR)

                passing_refs = apply_pipeline_batch(
                    passing_refs,
                    recalculate_global_qvalue_pipeline,
                    passing_psms_folder
                )
            else
                # No Step 10 needed when MBR is disabled
            end
        end
        #@debug_l1 "Step 10 completed in $(round(step10_time, digits=2)) seconds"

        # Update search context with passing PSM paths
        for (file_idx, ref) in zip(valid_file_indices, passing_refs)
            setPassingPsms!(getMSData(search_context), file_idx, file_path(ref))
        end

        # Step 11: Count protein peptides
        step11_time = @elapsed begin
            protein_to_possible_peptides = count_protein_peptides(
                getPrecursors(getSpecLib(search_context))
            )
        end

        # Step 12: Perform protein inference and initial scoring
        step12_time = @elapsed begin
            pg_refs, psm_to_pg_mapping = perform_protein_inference_pipeline(
                passing_refs,
                passing_proteins_folder,
                getPrecursors(getSpecLib(search_context)),
                protein_to_possible_peptides,
                min_peptides = params.min_peptides
            )
            
            paired_files = [PairedSearchFiles(psm_path, pg_path, file_idx)
                           for (file_idx, (psm_path, pg_path)) in zip(valid_file_indices, psm_to_pg_mapping)]
            
            isempty(paired_files) && error("No protein groups created during protein inference")
        end
        # Step 13: Build protein CV fold mapping from PSMs
        step13_time = @elapsed begin
            # Get PSM paths from passing_refs (these are the high-quality PSMs)
            psm_paths = [file_path(ref) for ref in passing_refs]
            protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, getPrecursors(getSpecLib(search_context)))
        end

        # Step 14: Perform protein probit regression
        step14_time = @elapsed begin

            qc_folder = joinpath(dirname(temp_folder), "qc_plots")
            !isdir(qc_folder) && mkdir(qc_folder)
            
            perform_protein_probit_regression(
                pg_refs,
                params.max_psms_in_memory,
                qc_folder,
                getPrecursors(getSpecLib(search_context));
                protein_to_cv_fold = protein_to_cv_fold,
                ms1_scoring = params.ms1_scoring
            )
        end

        # Step 15: Calculate global protein scores
        step15_time = @elapsed begin
            acc_to_max_pg_score = calculate_and_add_global_scores!(pg_refs)
        end

        # Step 16: Sort protein groups by global_pg_score
        step16_time = @elapsed begin
            sort_file_by_keys!(pg_refs, :global_pg_score, :target; reverse=[true, true])
        end

        # Step 17: Merge protein groups for global q-values
        step17_time = @elapsed begin
            sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
            stream_sorted_merge(pg_refs, sorted_pg_scores_path, :global_pg_score, :target;
                               batch_size=1000000, reverse=[true,true])
        end

        # Step 18: Calculate global protein q-values
        step18_time = @elapsed begin
            search_context.global_pg_score_to_qval[] = get_protein_global_qval_spline(sorted_pg_scores_path, params)
        end

        # Step 19: Sort protein groups by pg_score
        step19_time = @elapsed begin
            sort_file_by_keys!(pg_refs, :pg_score, :target; reverse=[true, true])
        end
        
        # Step 20: Merge protein groups for experiment-wide q-values
        step20_time = @elapsed begin
            stream_sorted_merge(pg_refs, sorted_pg_scores_path, :pg_score, :target;
                               batch_size=1000000, reverse=[true,true])
        end

        # Step 21: Calculate experiment-wide protein q-values and PEPs
        step21_time = @elapsed begin
            search_context.pg_score_to_qval[] = get_protein_qval_spline(sorted_pg_scores_path, params)
            search_context.pg_score_to_pep[]  = get_protein_pep_interpolation(sorted_pg_scores_path, params)
        end

        # Step 22: Add q-values and passing flags to protein groups
        step22_time = @elapsed begin
            protein_qval_pipeline = TransformPipeline() |>
                add_interpolated_column(:global_pg_qval, :global_pg_score, search_context.global_pg_score_to_qval[]) |>
                add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
                add_interpolated_column(:pg_pep, :pg_score, search_context.pg_score_to_pep[]) |> 
                filter_by_multiple_thresholds([
                    (:global_pg_qval, params.q_value_threshold),
                    (:pg_qval, params.q_value_threshold)
                ]) 
            apply_pipeline!(pg_refs, protein_qval_pipeline)
        end

        sort_file_by_keys!(pg_refs, :pg_score, :target; reverse=[true, true])
        stream_sorted_merge(pg_refs, sorted_pg_scores_path, :pg_score, :target;
                    batch_size=1000000, reverse=[true,true])

        search_context.pg_score_to_qval[] = get_protein_qval_spline(sorted_pg_scores_path, params)

        recalculate_experiment_wide_qvalue = TransformPipeline() |>
            add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[])

        pg_refs = apply_pipeline_batch(
                pg_refs,
                recalculate_experiment_wide_qvalue,
                passing_proteins_folder
            )

        # Step 23: Update PSMs with final protein scores
        step23_time = @elapsed begin
            update_psms_with_probit_scores_refs(
                paired_files,
                acc_to_max_pg_score,
                search_context.pg_score_to_qval[],
                search_context.global_pg_score_to_qval[]
            )
        end
        # Summary of all step times
        total_time = step1_time + step2_time + step3_time + step4_time + step5_time +
                    step6_time + step7_time + step8_time + step9_time + step10_time +
                    step11_time + step12_time + step13_time + step14_time + step15_time +
                    step16_time + step17_time + step18_time + step19_time + step20_time +
                    step21_time + step22_time + step23_time
    
        @user_info "ScoringSearch completed - Total time: $(round(total_time, digits=2)) seconds"

        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end