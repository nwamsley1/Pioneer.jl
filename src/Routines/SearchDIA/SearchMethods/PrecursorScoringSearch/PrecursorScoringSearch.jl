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

Per-precursor rescoring + FDR control. Trains a second LightGBM on globally-
filtered MainSearch PSMs (20 features), computes experiment-wide + global q-values/PEPs, filters
to a passing-precursor set, and builds RT indices for the chromatogram
integrator. Protein inference + protein-group scoring are now downstream
(`ProteinInferenceSearch`, `ProteinScoringSearch`).
"""
struct PrecursorScoringSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for scoring search.
"""
struct PrecursorScoringSearchResults <: SearchResults
    precursor_global_qval_dict::Base.Ref{Dict{UInt32, Float32}}  # precursor_idx → global q-value
    precursor_qval_interp::Base.Ref{Any}  # Interpolation for experiment-wide q-values
    precursor_pep_interp::Base.Ref{Any}   # Interpolation for experiment-wide PEPs
    merged_quant_path::String             # Path to merged quantification results
end

"""
Parameters for scoring search.
"""
struct PrecursorScoringSearchParameters <: SearchParameters
    # PSMs per bin for empirical q-value/PEP interpolation in get_qvalue_spline.
    # Smaller = finer-grained but noisier per-bin FDR estimates; larger =
    # smoother but coarser.
    pep_bin_size::Int64

    q_value_threshold::Float32

    # When false, skip MBR feature computation, the MBR-boosted second pass,
    # the FTR controller, and the qval bypass. Driven by global.match_between_runs.
    match_between_runs::Bool

    function PrecursorScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings

        # `match_between_runs` now selects one-round vs two-round scoring (its old MBR meaning is gone;
        # the MBR machinery has been deleted). Default true = do the second (cross-run) round.
        mbr = hasproperty(global_params, :match_between_runs) ?
                Bool(global_params.match_between_runs) : true

        new(
            Int64(ml_params.pep_bin_size),
            _resolve_q_value_threshold(global_params),
            mbr,
        )
    end
end

#==========================================================
Interface Implementation  
==========================================================#

get_parameters(::PrecursorScoringSearch, params::Any) = PrecursorScoringSearchParameters(params)

function init_search_results(::PrecursorScoringSearchParameters, search_context::SearchContext)
    return PrecursorScoringSearchResults(
        Ref(Dict{UInt32, Float32}()),  # precursor_global_qval_dict
        Ref(undef),  # precursor_qval_interp
        Ref(undef),  # precursor_pep_interp
        joinpath(getDataOutDir(search_context), "merged_quant.arrow")
    )
end

function process_file!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
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

"""
Process all results to get final protein scores.
"""
function summarize_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")

    # Set up output folders. PSM intermediates all live in `main_search_psms/`
    # since 2026-05-20 — MainSearch writes fold-split files there,
    # PrecursorScoring reads/merges/MBR-folds in place, and only the
    # post-FDR `passing_psms/` is a separate (and strictly smaller) output.
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    !isdir(passing_psms_folder) && mkdir(passing_psms_folder)

    # No outer try/catch — true errors should surface and abort the run.
    # Empty/sparse-input edge cases are handled inline below by early-return
    # checks rather than by swallowing exceptions.

    # Step 1: Train LightGBM Models
    ##@debug_l1 "Step 1: Training LightGBM models..."
    # Iterate every populated SecondPass PSM path
    valid_file_data = get_all_indexed_paths(getSecondPassPsms, search_context)
    valid_file_indices = [idx for (idx, _) in valid_file_data]

    if isempty(valid_file_data)
        @user_warn "No SecondPass PSM files registered for PrecursorScoringSearch"
        return nothing
    end

    # Get all fold-split file paths for LightGBM training and downstream processing
    # MainSearch now writes separate files per CV fold: *_fold0.arrow, *_fold1.arrow
    valid_fold_paths = get_all_fold_file_paths(search_context)

    if isempty(valid_fold_paths)
        @user_warn "No fold-split PSM files found for PrecursorScoringSearch"
        return nothing
    end

    step1_time = @elapsed begin
        score_precursor_isotope_traces(
            valid_fold_paths,
            getPrecursors(getSpecLib(search_context));
            match_between_runs = params.match_between_runs,
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

        # Collect data from both folds via PSMFileReference so the
        # auto-discovered mbr_outputs sidecars (written by
        # merge_mbr_sidecars_into_main!) are joined in. MainSearch
        # initialized :trace_prob to zero in the fold files; the sidecar
        # provides the real Pass-1 OOF score in :trace_prob_prepass, so we
        # set :trace_prob = :trace_prob_prepass here (matching the legacy
        # merge_mbr_sidecars_into_main! semantics).
        fold_dfs = DataFrame[]
        fold_refs = PSMFileReference[]
        function _load_fold(path)
            ref = PSMFileReference(path)
            push!(fold_refs, ref)
            df = load_with_sidecars(ref)
            if hasproperty(df, :trace_prob_prepass)
                df[!, :trace_prob] = df[!, :trace_prob_prepass]
            end
            return df
        end
        isfile(fold0_path) && push!(fold_dfs, _load_fold(fold0_path))
        isfile(fold1_path) && push!(fold_dfs, _load_fold(fold1_path))

        if !isempty(fold_dfs)
            # Merge and write combined file
            combined_df = vcat(fold_dfs...)
            writeArrow(merged_path, combined_df)
            push!(merged_psm_paths, merged_path)

            # Update search context with merged path
            setSecondPassPsms!(getMSData(search_context), idx, merged_path)

            # Collect fold file paths and their sidecars for batch deletion
            for ref in fold_refs
                fp = file_path(ref)
                isfile(fp) && push!(fold_paths_to_delete, fp)
                for s in ref.sidecars
                    isfile(s.path) && push!(fold_paths_to_delete, s.path)
                end
            end
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

    # Step 2: Aggregate trace-level to precursor-level probabilities (per-file)
    step2_time = @elapsed begin
        aggregate_per_file!(second_pass_refs)
    end

    # MainSearch already selects one row per precursor per file,
    # so no best-trace selection is needed. Pass refs through directly.
    filtered_refs = second_pass_refs

    # Steps 5-10 (combined): Build dictionaries + sidecar splines + single pipeline pass
    # Replaces 6 separate sort-merge-load-split cycles with:
    #   - Experiment-wide precursor scoring
    #   - In-memory q-value computation from dicts (no I/O)
    #   - Lightweight 2-column sidecar files for spline computation
    #   - Single per-file pipeline combining all column additions + filtering
    step5_10_time = @elapsed begin
        n_files_total = length(getFilePaths(getMSData(search_context)))
        fdr_scale = getLibraryFdrScaleFactor(search_context)

        # Pre-allocation size from spectral library
        n_precursors = length(getPrecursors(getSpecLib(search_context)))

        # A1: Build experiment-wide precursor scores.
        if n_files_total > 1
            @user_info "Training global precursor scoring model..."
        end
        global_prob_dict, target_dict =
            build_global_precursor_score_dicts(
                filtered_refs,
                n_precursors,
                n_files_total,
            )

        # A2: Compute global q-value AND global PEP dicts from global_prob dict (NO file I/O)
        global_qval_dict = build_global_qval_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        global_pep_dict  = build_global_pep_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        results.precursor_global_qval_dict[] = global_qval_dict

        # A3-A5: Sidecar lifecycle → q-value spline + PEP interpolation
        spline_result = build_qvalue_spline_from_refs(filtered_refs, :prec_prob, results.merged_quant_path;
            compute_pep=true, min_pep_points_per_bin=params.pep_bin_size,
            fdr_scale_factor=fdr_scale, temp_prefix="qval_sidecar")
        qval_spline = spline_result.qval_spline
        results.precursor_qval_interp[] = qval_spline
        results.precursor_pep_interp[] = spline_result.pep_interp

        # Phase B — Single per-file pipeline combining Steps 5+10.
        # MBR Phase 5b: rows with :mbr_recovered=true bypass the per-file
        # :qval filter (their qval is overridden to 0). The cross-run
        # :global_qval threshold still applies to all rows.
        mbr_qval_bypass = "mbr_recovery_qval_bypass" => function(df)
            if hasproperty(df, :mbr_recovered) && hasproperty(df, :qval)
                qv = df[!, :qval]
                rec = df[!, :mbr_recovered]
                @inbounds for i in 1:nrow(df)
                    if rec[i]
                        qv[i] = 0.0f0
                    end
                end
            end
            return df
        end

        # Cross-run filter on (:global_qval ≤ threshold) AND (:qval ≤ threshold).
        # The :pep and :off variants are retained in git history if needed.
        qval_conditions = [(:global_qval, params.q_value_threshold),
                           (:qval,        params.q_value_threshold)]
        combined_pipeline = TransformPipeline() |>
            add_dict_column(:global_prob, :precursor_idx, global_prob_dict) |>
            add_dict_column(:global_qval, :precursor_idx, global_qval_dict) |>
            add_dict_column(:global_pep,  :precursor_idx, global_pep_dict) |>
            add_interpolated_column(:qval, :prec_prob, qval_spline) |>
            mbr_qval_bypass |>
            add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
            filter_by_multiple_thresholds(qval_conditions)

        passing_refs = apply_pipeline_batch(filtered_refs, combined_pipeline, passing_psms_folder)
    end

    # After Step 5-10, the merged main_search_psms files are no longer read by
    # any downstream code (IntegrateChroms/MaxLFQ read passing_psms only). The
    # mid-pipeline cleanup that frees ~120 MB/file is temporarily disabled to
    # keep these files available for diagnostic inspection. Re-enable once the
    # main_search_psms column set has been pruned to a minimal/diagnostic schema.

    # NOTE: no post-filter q-value recalculation. The experiment-wide :qval (and :pep) written in
    # Phase B come from the spline built over the FULL target/decoy distribution (all filtered_refs),
    # which is the correct estimate. Rebuilding the spline on the AND-filter survivors was
    # entrapment-validated to be optimistic (decoy-depleted survivors) and is intentionally removed.

    # Update search context with passing PSM paths
    for (file_idx, ref) in zip(valid_file_indices, passing_refs)
        setPassingPsms!(getMSData(search_context), file_idx, file_path(ref))
    end

    # Build RT indices for IntegrateChromatogramsSearch (all library precursors per file)
    # The per-file count is logged inside build_rt_indices! at @debug_l1.
    build_rt_indices!(search_context, valid_file_indices, passing_refs)

    # Protein inference + protein-group scoring are handled downstream by
    # ProteinInferenceSearch and ProteinScoringSearch. Per-stage timings are
    # in the Performance Report at the end of the run.

    return nothing
end
