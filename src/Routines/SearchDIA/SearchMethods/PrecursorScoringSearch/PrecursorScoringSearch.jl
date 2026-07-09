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

    # When false, skip post-qvalue MBR feature computation and the FTR controller.
    # Driven by global.match_between_runs.
    match_between_runs::Bool
    main_search_params::Any

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
            MainSearchParameters(params),
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
        prec_probs = df[!, :prec_prob]
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

function gate_mbr_recovered_by_global_qvalue(
    global_qval_dict::Dict{UInt32, Float32},
    q_value_threshold::Float32,
)
    desc = "gate_mbr_recovered_by_global_qvalue"
    op = function(df)
        hasproperty(df, :mbr_recovered) || return df
        precursors = df[!, :precursor_idx]
        recovered = df[!, :mbr_recovered]
        target_decoy_probs = hasproperty(df, :mbr_target_decoy_prob) ?
                             df[!, :mbr_target_decoy_prob] :
                             nothing

        @inbounds for i in eachindex(recovered)
            Bool(recovered[i]) || continue
            global_qval = get(global_qval_dict, UInt32(precursors[i]), Inf32)
            global_qval <= q_value_threshold && continue
            recovered[i] = false
            target_decoy_probs !== nothing && (target_decoy_probs[i] = NaN32)
        end
        return df
    end
    return desc => op
end

function _initial_precursor_qvalue_filter(q_value_threshold::Float32)
    desc = "initial_precursor_qvalue_filter_with_mbr_recovery"
    op = function(df)
        qvals = df[!, :qval]
        global_qvals = df[!, :global_qval]
        recovered = hasproperty(df, :mbr_recovered) ? df[!, :mbr_recovered] : falses(nrow(df))

        keep = BitVector(undef, nrow(df))
        @inbounds for i in 1:nrow(df)
            global_pass = !ismissing(global_qvals[i]) && Float32(global_qvals[i]) <= q_value_threshold
            row_pass = !ismissing(qvals[i]) && Float32(qvals[i]) <= q_value_threshold
            keep[i] = global_pass && (row_pass || Bool(recovered[i]))
        end
        return df[keep, :]
    end
    return desc => op
end

function _normalize_optional_mbr_columns!(df::DataFrame)
    for col in (:mbr_recovered, :MBR_transfer_candidate)
        hasproperty(df, col) || continue
        values = df[!, col]
        df[!, col] = Bool[!ismissing(x) && Bool(x) for x in values]
    end
    for col in (:mbr_target_decoy_prob, :ftr_qval_true, :ftr_pep_true)
        hasproperty(df, col) || continue
        values = df[!, col]
        df[!, col] = Float32[ismissing(x) ? NaN32 : Float32(x) for x in values]
    end
    return df
end

function remap_mbr_recovered_prec_probs!(
    refs::Vector{PSMFileReference},
    merged_path::String;
    q_value_threshold::Float32,
    min_pep_points_per_bin::Int,
    fdr_scale_factor::Float32,
    global_qval_dict::Union{Nothing, Dict{UInt32, Float32}} = nothing,
)
    any(ref -> has_column_anywhere(ref, :mbr_target_decoy_prob), refs) || return refs

    spline_result = build_qvalue_spline_from_refs(refs, :prec_prob, merged_path;
        compute_pep = false,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
        temp_prefix = "mbr_score_sidecar")
    spline_result === nothing && return refs

    remap_pipeline = TransformPipeline()
    if global_qval_dict !== nothing
        remap_pipeline = remap_pipeline |>
            gate_mbr_recovered_by_global_qvalue(global_qval_dict, q_value_threshold)
    end
    remap_pipeline = remap_pipeline |>
        remap_mbr_recovered_prec_probs(spline_result.qval_spline, q_value_threshold)
    apply_pipeline!(refs, remap_pipeline; parallel = false)
    return refs
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
    main_search_psms_folder = joinpath(temp_folder, "main_search_psms")
    annotated_psms_folder = joinpath(temp_folder, "precursor_scored_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    !isdir(annotated_psms_folder) && mkdir(annotated_psms_folder)
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
        )
    end
    #@debug_l1 "Step 1 completed in $(round(step1_time, digits=2)) seconds"

    # Step 1b: Merge fold files back into single files per MS run
    # After ML scoring, we merge fold0 and fold1 files back together
    # This simplifies downstream processing which expects one file per MS run
    merged_psm_paths = String[]
    fold_paths_to_delete = String[]
    precursor_accessions = getAccessionNumbers(getPrecursors(getSpecLib(search_context)))
    for (idx, base_path) in valid_file_data
        fold0_path = "$(base_path)_fold0.arrow"
        fold1_path = "$(base_path)_fold1.arrow"
        merged_path = "$(base_path).arrow"

        # Collect data from both folds after Pass-1 scoring. MainSearch
        # initialized :trace_prob to zero in the fold files; Pass-1 scoring
        # writes the OOF score to :trace_prob_prepass, so keep :trace_prob
        # aligned with that non-MBR score for the initial q-value filters.
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
            combined_df = vcat(fold_dfs...; cols = :union)
            _normalize_optional_mbr_columns!(combined_df)
            combined_df[!, :accession_numbers] = [
                precursor_accessions[pid] for pid in combined_df[!, :precursor_idx]
            ]
            combined_df[!, :decoy] = combined_df[!, :target] .== false
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
    #   - Streaming dict accumulation for global_prob (reads ~12 bytes/row)
    #   - In-memory q-value computation from dicts (no I/O)
    #   - Lightweight 2-column sidecar files for spline computation
    #   - Single per-file pipeline combining all column additions + filtering
    step5_10_time = @elapsed begin
        n_files_total = length(getFilePaths(getMSData(search_context)))
        sqrt_n_runs = max(1, floor(Int64, sqrt(n_files_total)))
        fdr_scale = getLibraryFdrScaleFactor(search_context)

        # Pre-allocation size from spectral library
        n_precursors = length(getPrecursors(getSpecLib(search_context)))

        # A1: Stream per-file to build global_prob dictionaries (~12 bytes/row read)
        global_prob_dict, target_dict =
            build_precursor_global_prob_dicts(
                filtered_refs,
                sqrt_n_runs,
                n_precursors,
            )

        # A2: Compute global q-value AND global PEP dicts from global_prob dict (NO file I/O)
        global_qval_dict = build_global_qval_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        global_pep_dict  = build_global_pep_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        results.precursor_global_qval_dict[] = global_qval_dict

        remap_mbr_recovered_prec_probs!(
            filtered_refs,
            results.merged_quant_path;
            q_value_threshold = params.q_value_threshold,
            min_pep_points_per_bin = params.pep_bin_size,
            fdr_scale_factor = fdr_scale,
            global_qval_dict = global_qval_dict,
        )
        filtered_refs = [PSMFileReference(file_path(ref)) for ref in filtered_refs]

        # A3-A5: Sidecar lifecycle → q-value spline + PEP interpolation
        spline_result = build_qvalue_spline_from_refs(filtered_refs, :prec_prob, results.merged_quant_path;
            compute_pep=true, min_pep_points_per_bin=params.pep_bin_size,
            fdr_scale_factor=fdr_scale, temp_prefix="qval_sidecar")
        qval_spline = spline_result.qval_spline
        results.precursor_qval_interp[] = qval_spline
        results.precursor_pep_interp[] = spline_result.pep_interp

        # Cross-run filter on (:global_qval ≤ threshold) AND (:qval ≤ threshold).
        # The :pep and :off variants are retained in git history if needed.
        qval_conditions = [(:global_qval, params.q_value_threshold),
                           (:qval,        params.q_value_threshold)]
        scoring_pipeline = TransformPipeline() |>
            add_dict_column(:global_prob, :precursor_idx, global_prob_dict) |>
            add_dict_column(:global_qval, :precursor_idx, global_qval_dict) |>
            add_dict_column(:global_pep,  :precursor_idx, global_pep_dict) |>
            add_interpolated_column(:qval, :prec_prob, qval_spline) |>
            add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[])

        annotated_refs = apply_pipeline_batch(filtered_refs, scoring_pipeline, annotated_psms_folder)

        initial_filter_pipeline = TransformPipeline() |>
            _initial_precursor_qvalue_filter(params.q_value_threshold)

        passing_refs = apply_pipeline_batch(annotated_refs, initial_filter_pipeline, passing_psms_folder)
    end

    if params.match_between_runs && !isempty(passing_refs)
        post_mbr_time = @elapsed begin
            summary = prepare_mbr_after_qvalue_filter!(
                annotated_refs,
                passing_refs,
                getPrecursors(getSpecLib(search_context)),
                getFragmentLookupTable(getSpecLib(search_context)),
                passing_psms_folder;
                q_value_threshold = params.q_value_threshold,
                donor_q_threshold = MBR_DONOR_Q_THRESHOLD,
            )
            @debug_l1 "Post-qvalue MBR candidate staging completed: files=$(summary.n_files), " *
                      "candidates=$(summary.n_candidates), integration_rows=$(summary.n_integration_rows)"
            annotated_refs = [PSMFileReference(file_path(ref)) for ref in annotated_refs]
            passing_refs = summary.integration_refs
        end
        @debug_l1 "Post-qvalue MBR candidate staging elapsed: $(round(post_mbr_time, digits=2)) seconds"
    end

    # After Step 5-10, the merged main_search_psms files are no longer read by
    # any downstream code (IntegrateChroms/MaxLFQ read passing_psms only). The
    # mid-pipeline cleanup that frees ~120 MB/file is temporarily disabled to
    # keep these files available for diagnostic inspection. Re-enable once the
    # main_search_psms column set has been pruned to a minimal/diagnostic schema.

    # Step 11: Re-calculate q-values using filtered data (sidecar-based).
    # For MBR-enabled searches this moves to IntegrateChromatogramSearch, after
    # the MBR FTR controller has used integration-derived features.
    if !params.match_between_runs
        step11_time = @elapsed begin
        # Sidecar lifecycle for new spline (on filtered data)
            spline_result = build_qvalue_spline_from_refs(passing_refs, :prec_prob, results.merged_quant_path;
                compute_pep=true,
                min_pep_points_per_bin=params.pep_bin_size,
                fdr_scale_factor=getLibraryFdrScaleFactor(search_context), temp_prefix="recalc_sidecar")
            if spline_result === nothing
                @user_warn "No non-empty files for q-value recalculation — skipping Step 11"
            else
                results.precursor_qval_interp[] = spline_result.qval_spline
                results.precursor_pep_interp[] = spline_result.pep_interp
                recalc_pipeline = TransformPipeline() |>
                    add_interpolated_column(:qval, :prec_prob, spline_result.qval_spline) |>
                    filter_by_multiple_thresholds([(:qval, params.q_value_threshold)]) |>
                    add_interpolated_column(:pep, :prec_prob, spline_result.pep_interp)
                passing_refs = apply_pipeline_batch(passing_refs, recalc_pipeline, passing_psms_folder)
            end
        end
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
