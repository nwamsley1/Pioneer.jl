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
    PrecursorScoringSearch

Search method for precursor rescoring and precursor-level FDR control after the
second pass search has completed. This step performs ML-based precursor
rescoring, MBR-aware precursor aggregation, and final precursor filtering for
downstream chromatogram integration.
"""
struct PrecursorScoringSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for precursor scoring.
"""
struct PrecursorScoringSearchResults <: SearchResults
    best_traces::Dict{Int64, Float32}
    precursor_global_qval_dict::Base.Ref{Dict{UInt32, Float32}}
    precursor_qval_interp::Base.Ref{Any}
    precursor_pep_interp::Base.Ref{Any}
    merged_precursor_scores_path::String
end

"""
Parameters for precursor scoring.
"""
struct PrecursorScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    # Core memory and trace-selection parameters
    max_in_memory_table_mb::Float64
    min_best_trace_prob::Float32
    precursor_prob_spline_points_per_bin::Int64
    q_value_interpolation_points_per_bin::Int64

    # MBR and iterative rescoring parameters
    match_between_runs::Bool
    max_q_value_lightgbm_rescore::Float32
    max_q_value_mbr_itr::Float32
    min_PEP_neg_threshold_itr::Float32
    max_MBR_false_transfer_rate::Float32

    # Final precursor filtering parameters
    q_value_threshold::Float32
    isotope_tracetype::I
    n_quantile_bins::Int64

    # Model comparison parameters
    enable_model_comparison::Bool
    validation_split_ratio::Float64
    qvalue_threshold_comparison::Float64
    min_psms_for_comparison::Int64
    max_psms_for_comparison::Int64

    # OOM scoring parameters
    ms1_scoring::Bool
    force_oom::Bool
    max_mbr_training_candidates::Int64

    function PrecursorScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings

        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) && global_params.isotope_settings.combine_traces
            CombineTraces(0.0f0)
        else
            SeperateTraces()
        end

        new{typeof(isotope_trace_type)}(
            Float64(ml_params.max_in_memory_table_mb),
            Float32(ml_params.min_trace_prob),
            Int64(ml_params.spline_points),
            Int64(ml_params.q_value_interpolation_points_per_bin),
            Bool(global_params.match_between_runs),
            Float32(global_params.scoring.q_value_threshold),
            Float32(ml_params.max_q_value_mbr_itr),
            Float32(ml_params.min_PEP_neg_threshold_itr),
            Float32(global_params.scoring.q_value_threshold),
            Float32(global_params.scoring.q_value_threshold),
            isotope_trace_type,
            Int64(ml_params.n_quantile_bins),
            Bool(get(ml_params, :enable_model_comparison, true)),
            Float64(get(ml_params, :validation_split_ratio, 0.2)),
            Float64(get(ml_params, :qvalue_threshold, 0.01)),
            Int64(get(ml_params, :min_psms_for_comparison, 1000)),
            Int64(get(ml_params, :max_psms_for_comparison, 100000)),
            Bool(global_params.ms1_scoring),
            Bool(get(ml_params, :force_oom, false)),
            Int64(get(ml_params, :max_mbr_training_candidates, 1_000_000))
        )
    end
end

#==========================================================
Interface Implementation  
==========================================================#

get_parameters(::PrecursorScoringSearch, params::Any) = PrecursorScoringSearchParameters(params)

function init_search_results(::PrecursorScoringSearchParameters, search_context::SearchContext)
    return PrecursorScoringSearchResults(
        Dict{Int64, Float32}(),
        Ref(Dict{UInt32, Float32}()),
        Ref(undef),
        Ref(undef),
        joinpath(getDataOutDir(search_context), "merged_precursor_scores.arrow")
    )
end

function process_file!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "PrecursorScoringSearch")
        return results  # Return early with unchanged results
    end
    
    # No per-file processing needed
    return results
end

function process_search_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # No per-file results processing needed
    return nothing
end

function reset_results!(results::PrecursorScoringSearchResults)
    return nothing
end

#==========================================================
Precursor Score Calibration Helpers
==========================================================#

function get_precursor_global_qval_spline(merged_path::String, params::PrecursorScoringSearchParameters, search_context::SearchContext)
    score_col = params.match_between_runs ? :MBR_boosted_global_prob : :global_prob
    return get_qvalue_spline(
        merged_path, score_col, true;
        min_pep_points_per_bin = params.q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

function get_precursor_qval_spline(merged_path::String, params::PrecursorScoringSearchParameters, search_context::SearchContext)
    score_col = params.match_between_runs ? :MBR_boosted_prec_prob : :prec_prob
    return get_qvalue_spline(
        merged_path, score_col, false;
        min_pep_points_per_bin = params.q_value_interpolation_points_per_bin,
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    )
end

function get_precursor_pep_interpolation(merged_path::String, params::PrecursorScoringSearchParameters, search_context::SearchContext)
    score_col = params.match_between_runs ? :MBR_boosted_prec_prob : :prec_prob
    return get_pep_interpolation(
        merged_path, score_col;
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context),
    )
end

#==========================================================
Memory Estimation Helper
==========================================================#

"""
    estimate_max_rows(memory_mb::Float64, sample_file::String)

Estimate the maximum number of rows that fit in `memory_mb` given the schema of `sample_file`.
"""
function estimate_max_rows(memory_mb::Float64, sample_file::String)
    tbl = Arrow.Table(sample_file)
    bytes_per_row = 0
    for name in Tables.columnnames(tbl)
        col = Tables.getcolumn(tbl, name)
        T = eltype(col)
        T_inner = Base.nonmissingtype(T)
        if isbitstype(T_inner)
            bytes_per_row += sizeof(T_inner)
        else
            bytes_per_row += 64  # conservative estimate for strings / non-bits types
        end
        if T !== T_inner
            bytes_per_row += 1  # missing indicator
        end
    end
    bytes_per_row = max(bytes_per_row, 1)
    return max(floor(Int64, memory_mb * 1024 * 1024 / bytes_per_row), 1000)
end

#==========================================================
Main Precursor Scoring Pipeline
==========================================================#

"""
Process all results to produce final precursor scores and precursor-level q-values.
"""
function summarize_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    
    # Set up output folders
    second_pass_folder = joinpath(temp_folder, "second_pass_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    
    for folder in [passing_psms_folder]
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
            @user_warn "No valid files for precursor scoring - all files failed in previous search methods"
            return nothing
        end

        # Get all fold-split file paths for LightGBM training and downstream processing
        # SecondPassSearch now writes separate files per CV fold: *_fold0.arrow, *_fold1.arrow
        valid_fold_paths = get_valid_fold_file_paths(search_context)

        if isempty(valid_fold_paths)
            @user_warn "No valid fold-split PSM files found for precursor scoring"
            return nothing
        end

        step1_time = @elapsed begin
            max_psms = estimate_max_rows(params.max_in_memory_table_mb, first(valid_fold_paths))
            @user_info "Memory budget $(params.max_in_memory_table_mb) MB → max_psms = $max_psms"
            score_precursor_isotope_traces(
                second_pass_folder,
                valid_fold_paths,
                getPrecursors(getSpecLib(search_context)),
                params.match_between_runs,
                params.max_q_value_lightgbm_rescore,
                params.max_q_value_mbr_itr,
                params.min_PEP_neg_threshold_itr,
                max_psms,
                params.n_quantile_bins,
                params.q_value_threshold,
                params.ms1_scoring,
                params.force_oom
            )
        end
        #@debug_l1 "Step 1 completed in $(round(step1_time, digits=2)) seconds"

        # Step 1b: Merge fold files back into single files per MS run
        # After ML scoring, we merge fold0 and fold1 files back together
        # This simplifies downstream processing which expects one file per MS run
        merged_psm_paths = String[]
        fold_paths_to_delete = String[]
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

                # Collect fold file paths for batch deletion after loop
                isfile(fold0_path) && push!(fold_paths_to_delete, fold0_path)
                isfile(fold1_path) && push!(fold_paths_to_delete, fold1_path)
            end
        end

        # Release all mmap handles with a single GC, then batch-delete (Windows EACCES fix)
        GC.gc(false)
        for fpath in fold_paths_to_delete
            safeRm(fpath, nothing)
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
            spline_result = build_qvalue_spline_from_refs(filtered_refs, score_col, results.merged_precursor_scores_path;
                compute_pep=true, min_pep_points_per_bin=params.q_value_interpolation_points_per_bin,
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
            spline_result = build_qvalue_spline_from_refs(passing_refs, recalc_score_col, results.merged_precursor_scores_path;
                min_pep_points_per_bin=params.q_value_interpolation_points_per_bin,
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

        # Summary of all step times
        total_time = step1_time + step2_time + step3_time + step4_time +
                    step5_10_time + step11_time

        @user_info "Precursor scoring completed - Total time: $(round(total_time, digits=2)) seconds"

        best_traces = nothing # Free memory
    catch e
        @error "Failed to summarize precursor scoring results" exception=e
        rethrow(e)
    end

    return nothing
end
