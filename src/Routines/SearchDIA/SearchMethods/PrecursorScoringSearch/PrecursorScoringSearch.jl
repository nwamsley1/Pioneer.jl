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
filtered MainSearch PSMs (20 features), recalculates q-values + PEPs, filters
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
    precursor_qval_interp::Base.Ref{Any}  # Interpolation for run-specific q-values
    precursor_pep_interp::Base.Ref{Any}   # Interpolation for experiment-wide PEPs
    merged_quant_path::String             # Path to merged quantification results
end

"""
    FORCE_OOM

Developer toggle: force the out-of-memory scoring path even when PSMs would
fit in memory. Used to test/profile the OOM path without a deliberately
oversized dataset. Production runs leave this `false`.
"""
const FORCE_OOM = false

"""
Parameters for scoring search.
"""
struct PrecursorScoringSearchParameters <: SearchParameters
    # In-memory PSM memory budget; OOM path triggers when bytes-per-row * Q
    # exceeds this.
    max_psm_memory_mb::Float64

    # PSMs per bin for empirical q-value/PEP interpolation in get_qvalue_spline.
    # Smaller = finer-grained but noisier per-bin FDR estimates; larger =
    # smoother but coarser.
    pep_bin_size::Int64

    q_value_threshold::Float32

    # When false, skip MBR feature computation, the MBR-boosted second pass,
    # and the FTR controller. Driven by global.match_between_runs.
    match_between_runs::Bool

    function PrecursorScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings

        mbr = hasproperty(global_params, :match_between_runs) ?
                Bool(global_params.match_between_runs) : true

        new(
            Float64(ml_params.max_psm_memory_mb),
            Int64(ml_params.pep_bin_size),
            _resolve_q_value_threshold(global_params),
            mbr,
        )
    end
end

function _precursor_scoring_qvalue_filter_pipeline(
    global_prob_dict::Dict{UInt32, Float32},
    global_qval_dict::Dict{UInt32, Float32},
    global_pep_dict::Dict{UInt32, Float32},
    qval_spline,
    pep_interp,
    q_value_threshold::Float32,
)
    return TransformPipeline() |>
        add_dict_column(:global_prob, :precursor_idx, global_prob_dict) |>
        add_dict_column(:global_qval, :precursor_idx, global_qval_dict) |>
        add_dict_column(:global_pep,  :precursor_idx, global_pep_dict) |>
        add_interpolated_column(:qval, :prec_prob, qval_spline) |>
        _initial_precursor_qvalue_filter(q_value_threshold) |>
        add_interpolated_column(:pep, :prec_prob, pep_interp)
end

function _score_floor_for_qvalue(qval_spline, q_value_threshold::Float32)
    low = eps(Float32)
    high = 1.0f0 - eps(Float32)
    Float32(qval_spline(low)) <= q_value_threshold && return low
    Float32(qval_spline(high)) > q_value_threshold && return high

    for _ in 1:32
        mid = (low + high) / 2.0f0
        if Float32(qval_spline(mid)) <= q_value_threshold
            high = mid
        else
            low = mid
        end
    end
    return high
end

function remap_mbr_recovered_prec_probs(
    qval_spline,
    q_value_threshold::Float32,
)
    score_ceiling = prevfloat(_score_floor_for_qvalue(qval_spline, q_value_threshold))
    score_floor = eps(Float32)
    score_width = max(score_ceiling - score_floor, 0.0f0)
    desc = "remap_mbr_recovered_prec_probs"

    op = function(df)
        hasproperty(df, :mbr_target_decoy_prob) || return df
        prec_probs = df[!, :prec_prob]::AbstractVector{Float32}
        recovered = df[!, :mbr_recovered]
        target_decoy_probs = df[!, :mbr_target_decoy_prob]

        @inbounds for i in eachindex(prec_probs)
            Bool(recovered[i]) || continue
            score = clamp(Float32(target_decoy_probs[i]), 0.0f0, 1.0f0)
            prec_probs[i] = score_floor + score_width * score
        end
        return df
    end
    return desc => op
end

function remap_mbr_recovered_prec_probs!(
    refs::Vector{PSMFileReference},
    merged_path::String;
    q_value_threshold::Float32,
    min_pep_points_per_bin::Int,
    fdr_scale_factor::Float32,
)
    any(ref -> has_column_anywhere(ref, :mbr_target_decoy_prob), refs) || return refs

    spline_result = build_qvalue_spline_from_refs(refs, :prec_prob, merged_path;
        compute_pep=false,
        min_pep_points_per_bin=min_pep_points_per_bin,
        fdr_scale_factor=fdr_scale_factor,
        temp_prefix="mbr_score_sidecar")
    remap_pipeline = TransformPipeline() |>
        remap_mbr_recovered_prec_probs(spline_result.qval_spline, q_value_threshold)
    apply_pipeline!(refs, remap_pipeline; parallel=false)
    return refs
end

function _initial_precursor_qvalue_filter(q_value_threshold::Float32)
    desc = "initial_precursor_qvalue_filter_with_mbr_recovery"
    op = function(df)
        qvals = df[!, :qval]
        global_qvals = df[!, :global_qval]
        recovered = df[!, :mbr_recovered]

        keep = BitVector(undef, nrow(df))
        @inbounds for i in 1:nrow(df)
            global_pass = !ismissing(global_qvals[i]) && global_qvals[i] <= q_value_threshold
            row_pass = !ismissing(qvals[i]) && qvals[i] <= q_value_threshold
            keep[i] = global_pass && (row_pass || recovered[i])
        end

        df_filtered = df[keep, :]
        empty!(df)
        append!(df, df_filtered)
        return df
    end
    return desc => op
end

function _precursor_scoring_recalculated_qvalue_filter_pipeline(
    qval_spline,
    pep_interp,
    q_value_threshold::Float32,
)
    return TransformPipeline() |>
        add_interpolated_column(:qval, :prec_prob, qval_spline) |>
        filter_by_multiple_thresholds([(:qval, q_value_threshold)]) |>
        add_interpolated_column(:pep, :prec_prob, pep_interp)
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

function _combine_precursor_scoring_fold_dfs(fold_dfs::Vector{DataFrame})
    isempty(fold_dfs) && return DataFrame()
    return vcat(fold_dfs...; cols = :union)
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
    main_search_psms_folder = joinpath(temp_folder, "main_search_psms")
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

    mbr_rescue_fold_paths = params.match_between_runs ?
                            get_existing_mbr_rescue_fold_paths(valid_fold_paths) :
                            String[]
    !isempty(mbr_rescue_fold_paths) && @debug_l1 "MBR rescue pool: found " *
        "$(length(mbr_rescue_fold_paths)) fold files before cross-run scoring"

    max_psms = estimate_max_rows(params.max_psm_memory_mb, first(valid_fold_paths))
    @debug_l1 "Memory budget $(params.max_psm_memory_mb) MB → max_psms = $max_psms"
    score_precursor_isotope_traces(
        main_search_psms_folder,
        valid_fold_paths,
        getPrecursors(getSpecLib(search_context)),
        getFragmentLookupTable(getSpecLib(search_context)),
        max_psms,
        params.q_value_threshold,
        FORCE_OOM;
        match_between_runs = params.match_between_runs,
        mbr_rescue_file_paths = mbr_rescue_fold_paths,
    )

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
        if isfile(fold0_path)
            push!(fold_dfs, _load_fold(fold0_path))
        end
        if isfile(fold1_path)
            push!(fold_dfs, _load_fold(fold1_path))
        end

        for fold in UInt8[0, 1]
            rescue_path = mbr_rescue_recovered_path(mbr_rescue_fold_path("$(base_path)_fold$(fold).arrow"))
            isfile(rescue_path) || continue
            rescue_df = _load_fold(rescue_path)
            nrow(rescue_df) > 0 && push!(fold_dfs, rescue_df)
        end

        if !isempty(fold_dfs)
            # Merge and write combined file
            combined_df = _combine_precursor_scoring_fold_dfs(fold_dfs)
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
    aggregate_per_file!(second_pass_refs)

    # MainSearch already selects one row per precursor per file,
    # so no best-trace selection is needed. Pass refs through directly.
    filtered_refs = second_pass_refs

    # Steps 5-10 (combined): Build dictionaries + sidecar splines + single pipeline pass
    # Replaces 6 separate sort-merge-load-split cycles with:
    #   - Streaming dict accumulation for global_prob (reads ~12 bytes/row)
    #   - In-memory q-value computation from dicts (no I/O)
    #   - Lightweight 2-column sidecar files for spline computation
    #   - Single per-file pipeline combining all column additions + filtering
    n_files_total = length(getFilePaths(getMSData(search_context)))
    sqrt_n_runs = max(1, floor(Int64, sqrt(n_files_total)))
    fdr_scale = getLibraryFdrScaleFactor(search_context)

    # Pre-allocation size from spectral library
    n_precursors = length(getPrecursors(getSpecLib(search_context)))

    remap_mbr_recovered_prec_probs!(
        filtered_refs,
        results.merged_quant_path;
        q_value_threshold=params.q_value_threshold,
        min_pep_points_per_bin=params.pep_bin_size,
        fdr_scale_factor=fdr_scale,
    )

    # A1: Stream per-file to build global_prob dictionaries (~12 bytes/row read)
    global_prob_dict, target_dict =
        build_precursor_global_prob_dicts(filtered_refs, sqrt_n_runs, n_precursors)

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

    # Initial filter: all rows must pass :global_qval. Non-MBR rows also
    # pass the initial row-level :qval; MBR-recovered rows wait for the
    # recalculated row-level q-value below.
    # The :pep and :off variants are retained in git history if needed.
    combined_pipeline = _precursor_scoring_qvalue_filter_pipeline(
        global_prob_dict,
        global_qval_dict,
        global_pep_dict,
        qval_spline,
        results.precursor_pep_interp[],
        params.q_value_threshold,
    )

    passing_refs = apply_pipeline_batch(filtered_refs, combined_pipeline, passing_psms_folder)

    # After Step 5-10, the merged main_search_psms files are no longer read by
    # any downstream code (IntegrateChroms/MaxLFQ read passing_psms only). The
    # mid-pipeline cleanup that frees ~120 MB/file is temporarily disabled to
    # keep these files available for diagnostic inspection. Re-enable once the
    # main_search_psms column set has been pruned to a minimal/diagnostic schema.

    # Step 11: Re-calculate q-values using filtered data (sidecar-based), now
    # including MBR-recovered rows that passed the global q-value filter. This
    # second row-level q-value is applied to every row, including MBR transfers.
    # Sidecar lifecycle for new spline (on filtered data)
    spline_result = build_qvalue_spline_from_refs(passing_refs, :prec_prob, results.merged_quant_path;
        compute_pep=true,
        min_pep_points_per_bin=params.pep_bin_size,
        fdr_scale_factor=getLibraryFdrScaleFactor(search_context), temp_prefix="recalc_sidecar")
    if spline_result === nothing
        @user_warn "No non-empty files for q-value recalculation — skipping Step 11"
    else
        recalc_pipeline = _precursor_scoring_recalculated_qvalue_filter_pipeline(
            spline_result.qval_spline,
            spline_result.pep_interp,
            params.q_value_threshold,
        )
        passing_refs = apply_pipeline_batch(passing_refs, recalc_pipeline, passing_psms_folder)
    end

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
