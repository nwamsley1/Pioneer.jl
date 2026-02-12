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
    precursor_global_qval_dict::Base.Ref{Dict{UInt32, Float32}} # Dictionary mapping precursor_idx to global q-values
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

    # OOM scoring parameters
    force_oom::Bool
    max_training_psms::Int64
    max_mbr_training_candidates::Int64

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
            Bool(global_params.ms1_scoring),

            # OOM scoring parameters
            Bool(get(ml_params, :force_oom, false)),
            Int64(get(ml_params, :max_training_psms, ml_params.max_psms_in_memory)),
            Int64(get(ml_params, :max_mbr_training_candidates, 1_000_000))
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
        Ref(Dict{UInt32, Float32}()),  # precursor_global_qval_dict
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

        # Check if any valid files remain
        if isempty(valid_file_data)
            @user_warn "No valid files for ScoringSearch - all files failed in previous search methods"
            return nothing
        end

        # Get all fold-split file paths for LightGBM training and downstream processing
        # SecondPassSearch now writes separate files per CV fold: *_fold0.arrow, *_fold1.arrow
        valid_fold_paths = get_valid_fold_file_paths(search_context)

        if isempty(valid_fold_paths)
            @user_warn "No valid fold-split PSM files found for ScoringSearch"
            return nothing
        end

        step1_time = @elapsed begin
            score_precursor_isotope_traces(
                second_pass_folder,
                valid_fold_paths,
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs,
                params.max_q_value_lightgbm_rescore,
                params.max_q_value_mbr_itr,
                params.min_PEP_neg_threshold_itr,
                params.max_psms_in_memory,
                params.n_quantile_bins,
                params.q_value_threshold,
                params.ms1_scoring,
                params.force_oom,
                params.max_training_psms
            )
        end
        #@debug_l1 "Step 1 completed in $(round(step1_time, digits=2)) seconds"

        # Step 1b: Merge fold files back into single files per MS run
        # After ML scoring, we merge fold0 and fold1 files back together
        # This simplifies downstream processing which expects one file per MS run
        merged_psm_paths = String[]
        for (idx, base_path) in valid_file_data
            fold0_path = "$(base_path)_fold0.arrow"
            fold1_path = "$(base_path)_fold1.arrow"
            merged_path = "$(base_path).arrow"

            # Collect data from both folds
            fold_dfs = DataFrame[]
            if isfile(fold0_path)
                push!(fold_dfs, DataFrame(Arrow.Table(fold0_path)))
            end
            if isfile(fold1_path)
                push!(fold_dfs, DataFrame(Arrow.Table(fold1_path)))
            end

            if !isempty(fold_dfs)
                # Merge and write combined file
                combined_df = vcat(fold_dfs...)
                writeArrow(merged_path, combined_df)
                push!(merged_psm_paths, merged_path)

                # Update search context with merged path
                setSecondPassPsms!(getMSData(search_context), idx, merged_path)

                # Delete fold files to save disk space (merged file now contains all data)
                isfile(fold0_path) && rm(fold0_path)
                isfile(fold1_path) && rm(fold1_path)
            end
        end

        # Create references for second pass PSMs (now using merged files)
        second_pass_paths = merged_psm_paths
        second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

        # Step 2: Apply MBR filtering and calculate precursor probabilities (per-file OOM)
        step2_time = @elapsed begin
            apply_mbr_filter_and_aggregate_per_file!(second_pass_refs, valid_file_indices, params)
        end
        #@debug_l1 "Step 2 completed in $(round(step2_time, digits=2)) seconds"

        # Step 3: Find Best Isotope Traces
        #@debug_l1 "Step 3: Finding best isotope traces..."
        step3_time = @elapsed begin
            best_traces = get_best_traces(
                second_pass_paths,
                params.min_best_trace_prob
            )
        end
        #@debug_l1 "Step 3 completed in $(round(step3_time, digits=2)) seconds"

        # Step 4: Process Quantification Results
        #@debug_l1 "Step 4: Processing quantification results..."
        step4_time = @elapsed begin
            necessary_cols = get_quant_necessary_columns(params.match_between_runs)

            # Note: Removed sorting by global_prob - it will be calculated in Step 5
            quant_processing_pipeline = TransformPipeline() |>
                add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
                select_columns(vcat(necessary_cols, :best_trace)) |>
                filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
                remove_columns(:best_trace)

            apply_pipeline!(second_pass_refs, quant_processing_pipeline)
            filtered_refs = second_pass_refs
        end
        #@debug_l1 "Step 4 completed in $(round(step4_time, digits=2)) seconds"

        # Steps 5-10 (combined): Build dictionaries + sidecar splines + single pipeline pass
        # Replaces 6 separate sort-merge-load-split cycles with:
        #   - Streaming dict accumulation for global_prob (reads ~12 bytes/row)
        #   - In-memory q-value computation from dicts (no I/O)
        #   - Lightweight 2-column sidecar files for spline computation
        #   - Single per-file pipeline combining all column additions + filtering
        step5_10_time = @elapsed begin
            sqrt_n_runs = floor(Int64, sqrt(length(getFilePaths(getMSData(search_context)))))
            fdr_scale = getLibraryFdrScaleFactor(search_context)
            has_mbr = params.match_between_runs

            # Pre-allocation size from spectral library
            n_precursors = length(getPrecursors(getSpecLib(search_context)))

            # A1: Stream per-file to build global_prob dictionaries (~12 bytes/row read)
            global_prob_dict, mbr_global_prob_dict, target_dict =
                build_precursor_global_prob_dicts(filtered_refs, sqrt_n_runs, has_mbr, n_precursors)

            # A2: Compute global q-value dict from global_prob dict (NO file I/O)
            score_dict_for_qval = has_mbr ? mbr_global_prob_dict : global_prob_dict
            global_qval_dict = build_global_qval_dict_from_scores(score_dict_for_qval, target_dict, fdr_scale)
            results.precursor_global_qval_dict[] = global_qval_dict

            # A3-A5: Sidecar lifecycle → q-value spline + PEP interpolation
            score_col = has_mbr ? :MBR_boosted_prec_prob : :prec_prob
            spline_result = build_qvalue_spline_from_refs(filtered_refs, score_col, results.merged_quant_path;
                compute_pep=true, min_pep_points_per_bin=params.precursor_q_value_interpolation_points_per_bin,
                fdr_scale_factor=fdr_scale, temp_prefix="qval_sidecar")
            qval_spline = spline_result.qval_spline
            results.precursor_qval_interp[] = qval_spline
            results.precursor_pep_interp[] = spline_result.pep_interp

            # Phase B — Single per-file pipeline combining Steps 5+10
            global_qval_col = has_mbr ? :MBR_boosted_global_qval : :global_qval
            qval_col = has_mbr ? :MBR_boosted_qval : :qval

            combined_pipeline = TransformPipeline() |>
                add_dict_column(:global_prob, :precursor_idx, global_prob_dict)

            if has_mbr
                combined_pipeline = combined_pipeline |>
                    add_dict_column(:MBR_boosted_global_prob, :precursor_idx, mbr_global_prob_dict)
            end

            combined_pipeline = combined_pipeline |>
                add_dict_column(global_qval_col, :precursor_idx, global_qval_dict) |>
                add_interpolated_column(qval_col, score_col, qval_spline) |>
                add_interpolated_column(:pep, score_col, results.precursor_pep_interp[]) |>
                filter_by_multiple_thresholds([
                    (global_qval_col, params.q_value_threshold),
                    (qval_col, params.q_value_threshold)
                ])

            passing_refs = apply_pipeline_batch(filtered_refs, combined_pipeline, passing_psms_folder)
        end

        # Step 11: Re-calculate q-values using filtered data (sidecar-based)
        step11_time = @elapsed begin
            # Determine score column from filtered data
            sample_tbl = Arrow.Table(file_path(passing_refs[1]))
            has_mbr_cols = hasproperty(sample_tbl, :MBR_boosted_prec_prob)
            sample_tbl = nothing

            recalc_score_col = has_mbr_cols ? :MBR_boosted_prec_prob : :prec_prob
            recalc_qval_col = has_mbr_cols ? :MBR_boosted_qval : :qval

            # Sidecar lifecycle for new spline (on filtered data)
            spline_result = build_qvalue_spline_from_refs(passing_refs, recalc_score_col, results.merged_quant_path;
                min_pep_points_per_bin=params.precursor_q_value_interpolation_points_per_bin,
                fdr_scale_factor=getLibraryFdrScaleFactor(search_context), temp_prefix="recalc_sidecar")
            if spline_result === nothing
                @user_warn "No non-empty files for q-value recalculation — skipping Step 11"
            else
                recalc_pipeline = TransformPipeline() |>
                    add_interpolated_column(recalc_qval_col, recalc_score_col, spline_result.qval_spline)
                passing_refs = apply_pipeline_batch(passing_refs, recalc_pipeline, passing_psms_folder)
            end
        end

        # Update search context with passing PSM paths
        for (file_idx, ref) in zip(valid_file_indices, passing_refs)
            setPassingPsms!(getMSData(search_context), file_idx, file_path(ref))
        end

        # Step 12: Count protein peptides
        step12_time = @elapsed begin
            protein_to_possible_peptides = count_protein_peptides(
                getPrecursors(getSpecLib(search_context))
            )
        end

        # Step 13: Perform protein inference and initial scoring
        step13_time = @elapsed begin
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
        # Step 14: Build protein CV fold mapping from PSMs
        step14_time = @elapsed begin
            # Get PSM paths from passing_refs (these are the high-quality PSMs)
            psm_paths = [file_path(ref) for ref in passing_refs]
            protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, getPrecursors(getSpecLib(search_context)))
        end

        # Step 15: Perform protein probit regression
        step15_time = @elapsed begin

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

        # Steps 16-23 (combined): Build protein dicts + sidecar splines + single pipeline pass
        # Replaces 4 separate sort-merge-load-split cycles with:
        #   - Streaming dict accumulation for global_pg_score (reads ~20 bytes/row)
        #   - In-memory q-value computation from dicts (no I/O)
        #   - Lightweight 2-column sidecar files for spline computation
        #   - Single per-file pipeline combining all column additions + filtering
        step16_23_time = @elapsed begin
            sqrt_n_runs = floor(Int, sqrt(length(pg_refs)))

            # Pre-allocation size from spectral library
            n_proteins = length(getProteins(getSpecLib(search_context)))

            # D1: Stream per-file PGs to build global_pg_score dict (~20 bytes/row read)
            global_pg_score_dict, pg_name_to_global_pg_score =
                build_protein_global_score_dicts(pg_refs, sqrt_n_runs, n_proteins)
            search_context.pg_name_to_global_pg_score[] = pg_name_to_global_pg_score

            # D2: Compute global PG q-value dict from score dict (NO file I/O)
            global_pg_qval_dict = build_protein_global_qval_dict(global_pg_score_dict)
            search_context.global_pg_score_to_qval_dict[] = global_pg_qval_dict

            # D3-D5: Sidecar lifecycle → PG spline + PEP interpolation
            sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
            spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
                batch_size=1_000_000, compute_pep=true,
                min_pep_points_per_bin=params.pg_q_value_interpolation_points_per_bin, temp_prefix="pg_sidecar")
            search_context.pg_score_to_qval[] = spline_result.qval_spline
            search_context.pg_score_to_pep[] = spline_result.pep_interp

            # Phase B — Single per-file pipeline combining Steps 16+23
            protein_combined_pipeline = TransformPipeline() |>
                add_dict_column_composite_key(:global_pg_score, [:protein_name, :target, :entrap_id], global_pg_score_dict) |>
                add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], global_pg_qval_dict) |>
                add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
                add_interpolated_column(:pg_pep, :pg_score, search_context.pg_score_to_pep[]) |>
                filter_by_multiple_thresholds([
                    (:global_pg_qval, params.q_value_threshold),
                    (:pg_qval, params.q_value_threshold)
                ])

            apply_pipeline!(pg_refs, protein_combined_pipeline)

            # Post-filtering recalculation of pg_qval on filtered data
            spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
                batch_size=1_000_000, min_pep_points_per_bin=params.pg_q_value_interpolation_points_per_bin,
                temp_prefix="pg_recalc")
            search_context.pg_score_to_qval[] = spline_result.qval_spline

            recalc_pg_pipeline = TransformPipeline() |>
                add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[])

            pg_refs = apply_pipeline_batch(pg_refs, recalc_pg_pipeline, passing_proteins_folder)
        end

        # Step 24: Update PSMs with final protein scores
        step24_time = @elapsed begin
            update_psms_with_probit_scores_refs(
                paired_files,
                search_context.pg_name_to_global_pg_score[],
                search_context.pg_score_to_qval[],
                search_context.global_pg_score_to_qval_dict[]
            )
        end
        # Summary of all step times
        total_time = step1_time + step2_time + step3_time + step4_time +
                    step5_10_time + step11_time +
                    step12_time + step13_time + step14_time + step15_time +
                    step16_23_time + step24_time

        @user_info "ScoringSearch completed - Total time: $(round(total_time, digits=2)) seconds"

        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize scoring results" exception=e
        rethrow(e)
    end

    return nothing
end