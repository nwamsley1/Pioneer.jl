# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

function _write_mbr_pass1_sidecars_from_main!(
    file_paths::Vector{String},
)
    started = last_progress = time()
    n_written = 0
    n_rows = 0
    @debug_l1 "MBR Pass-1 sidecar preparation starting: files=$(length(file_paths))"
    for path in file_paths
        main = Arrow.Table(path)
        for column in (
            :precursor_idx,
            :scan_idx,
            :trace_prob_prepass,
            :trace_prob_infold,
        )
            hasproperty(main, column) ||
                error("Post-integration MBR requires column $column in $path")
        end
        writeArrow(
            path * PASS1_SIDECAR_SUFFIX,
            DataFrame(
                precursor_idx = collect(UInt32.(main.precursor_idx)),
                scan_idx = collect(UInt32.(main.scan_idx)),
                trace_prob_prepass =
                    collect(Float32.(main.trace_prob_prepass)),
                trace_prob_infold =
                    collect(Float32.(main.trace_prob_infold)),
            ),
        )
        n_written += 1
        n_rows += length(main.precursor_idx)
        if time() - last_progress >= 60
            @debug_l1 "MBR Pass-1 sidecar preparation: files=$n_written/$(length(file_paths)) rows=$n_rows elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "MBR Pass-1 sidecar preparation complete: files=$n_written rows=$n_rows elapsed=$(round(time() - started, digits=2))s"
    return n_written
end

@inline function _mbr_initial_pass(
    qval,
    global_qval,
    q_value_threshold::Float32,
)
    q = Float32(qval)
    global_q = Float32(global_qval)
    return isfinite(q) && q <= q_value_threshold &&
           isfinite(global_q) && global_q <= q_value_threshold
end

"""
    write_staged_selection(staged_path, source_ref, rows)

Record a staged MBR integration input as `rows` of `source_ref` (with its sidecars) instead of
writing a copy of those rows.
"""
function write_staged_selection(staged_path::String, source_ref::PSMFileReference, rows::Vector{Int})
    metadata = Dict(
        "source" => file_path(source_ref),
        "sidecars" => join((s.path for s in source_ref.sidecars), '\n'),
    )
    # A new file nothing has memory-mapped, so a direct write (with metadata) is safe.
    Arrow.write(staged_path * MBR_SELECTION_SUFFIX, DataFrame(source_row = UInt32.(rows));
                metadata = metadata)
    return nothing
end

"""
    load_staged_psms(path[, cols]) -> DataFrame

The PSM table at `path`, or, while `path` is a staged MBR selection (before integration writes
it), the selected rows of its source table and sidecars. `cols` limits the columns read.
"""
function load_staged_psms(path::String, cols::Union{Nothing, Vector{Symbol}} = nothing)
    selection_path = path * MBR_SELECTION_SUFFIX
    # Read and unmap: these files are replaced or deleted later in the run.
    isfile(selection_path) || return load_arrow_dataframe(path; cols)
    source_path, sidecar_list, rows = with_arrow_table(selection_path) do selection
        metadata = Arrow.getmetadata(selection)
        (String(metadata["source"]), String(metadata["sidecars"]), Int.(selection.source_row))
    end
    source = PSMFileReference(source_path;
        sidecar_paths = isempty(sidecar_list) ? String[] : split(sidecar_list, '\n'))
    return _selected_rows(source, rows, cols)
end

# The selected rows of the source table and its sidecars, gathered column by column from the
# memory-mapped Arrow files: never the whole table (staging keeps 12-31% of rows). Columns follow
# `load_with_sidecars` order (main table, then each sidecar's), restricted to `cols` when given.
# Each file is unmapped once its rows are copied: these files are replaced or deleted later.
function _selected_rows(source::PSMFileReference, rows::Vector{Int}, cols::Union{Nothing, Vector{Symbol}})
    df = DataFrame()
    wanted(c) = cols === nothing || c in cols
    with_arrow_table(file_path(source)) do main
        for c in Tables.columnnames(main)
            wanted(c) && (df[!, c] = _selected_owned(main, c, rows, file_path(source)))
        end
    end
    for s in source.sidecars
        with_arrow_table(s.path) do side
            for c in s.cols
                wanted(c) && (df[!, c] = _selected_owned(side, c, rows, s.path))
            end
        end
    end
    return df
end

# Selected rows of one column, copied so they own their memory (list cells too).
function _selected_owned(tbl, c::Symbol, rows::Vector{Int}, path::String)
    values = _owned_column(Tables.getcolumn(tbl, c)[rows])
    _owns_memory(values) || error("Column $c of $path has type $(typeof(values)), " *
                                  "which may still reference the mapped file")
    return values
end

"""Remove the staged selection at `path` once the real table has been written there."""
clear_staged_selection!(path::String) =
    (isfile(path * MBR_SELECTION_SUFFIX) && safeRm(path * MBR_SELECTION_SUFFIX; force = true); nothing)

# Function barrier: the column types are only known once the table is read.
function _select_mbr_integration_rows(
    precursor_idx, ms_file_idx, qval, global_qval,
    donor_files::Dict{UInt32, Tuple{UInt32, UInt32}},
    q_value_threshold::Float32,
)
    rows = Int[]
    n_candidates = 0
    @inbounds for row in eachindex(precursor_idx)
        baseline = _mbr_initial_pass(qval[row], global_qval[row], q_value_threshold)
        global_q = Float32(global_qval[row])
        global_pass = isfinite(global_q) && global_q <= q_value_threshold
        run_q = Float32(qval[row])
        run_pass = isfinite(run_q) && run_q <= q_value_threshold
        candidate =
            global_pass && !run_pass &&
            _mbr_has_cross_run_donor(
                donor_files,
                UInt32(precursor_idx[row]),
                UInt32(ms_file_idx[row]),
            )
        (baseline || candidate) && push!(rows, row)
        n_candidates += candidate
    end
    return rows, n_candidates
end

# Takes REFS rather than paths: PrecursorScoringSearch attaches :qval/:global_qval/:global_prob/
# :global_pep/:pep as a row-aligned sidecar, so the candidate tables only carry them via their
# PSMFileReference. Rows are chosen from the four columns the rule reads, and the staged input is
# written as a row selection of the candidate table (see MBR_SELECTION_SUFFIX), not as a copy;
# the sidecar is consolidated when integration loads the selection.
function _stage_mbr_integration_inputs!(
    candidate_refs::Vector{PSMFileReference},
    output_folder::String,
    donor_files::Dict{UInt32, Tuple{UInt32, UInt32}},
    q_value_threshold::Float32,
)
    started = last_progress = time()
    @debug_l1 "MBR integration input staging starting: files=$(length(candidate_refs))"
    mkpath(output_folder)
    refs = PSMFileReference[]
    n_rows = 0
    n_candidates = 0
    rows_processed = 0
    for (file_idx, candidate_ref) in enumerate(candidate_refs)
        path = file_path(candidate_ref)
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing MBR Pass-1 sidecar at $pass1_path")
        decide = materialize_columns(candidate_ref, [:precursor_idx, :ms_file_idx, :qval, :global_qval])
        pass1 = load_arrow_dataframe(pass1_path)
        nrow(decide) == nrow(pass1) ||
            error("MBR Pass-1 sidecar row-count mismatch at $pass1_path")

        rows, file_candidates = _select_mbr_integration_rows(
            decide.precursor_idx, decide.ms_file_idx, decide.qval, decide.global_qval,
            donor_files, q_value_threshold,
        )
        n_candidates += file_candidates

        staged_path = joinpath(output_folder, basename(path))
        # The staged path may hold the q-value-filtered table; the selection supersedes it.
        isfile(staged_path) && safeRm(staged_path; force = true)
        write_staged_selection(staged_path, candidate_ref, rows)
        writeArrow(staged_path * PASS1_SIDECAR_SUFFIX, pass1[rows, :])
        push!(refs, PSMFileReference(staged_path))
        n_rows += length(rows)
        rows_processed += nrow(decide)
        if time() - last_progress >= 60
            @debug_l1 "MBR integration input staging: files=$file_idx/$(length(candidate_refs)) rows=$rows_processed retained=$n_rows candidates=$n_candidates elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "MBR integration input staging complete: files=$(length(candidate_refs)) rows=$rows_processed retained=$n_rows candidates=$n_candidates elapsed=$(round(time() - started, digits=2))s"
    return (
        integration_refs = refs,
        n_rows = n_rows,
        n_candidates = n_candidates,
    )
end

"""
    prepare_postintegration_mbr!(candidate_refs, donor_refs, output_folder; ...)

Stage baseline IDs plus globally-supported, run-level failures for chromatogram
integration. Donor availability is evaluated from the frozen pre-MBR scores;
no MBR score can influence the candidate cohort.
"""
function prepare_postintegration_mbr!(
    candidate_refs::Vector{PSMFileReference},
    donor_refs::Vector{PSMFileReference},
    output_folder::String;
    q_value_threshold::Float32,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
)
    candidate_paths = String[
        file_path(ref) for ref in candidate_refs if exists(ref)
    ]
    donor_paths = String[
        file_path(ref) for ref in donor_refs if exists(ref)
    ]
    if isempty(candidate_paths) || isempty(donor_paths)
        return (
            integration_refs = PSMFileReference[],
            n_files = 0,
            n_rows = 0,
            n_candidates = 0,
            donor_score_floor = Inf32,
        )
    end

    _write_mbr_pass1_sidecars_from_main!(
        unique(vcat(candidate_paths, donor_paths)),
    )
    donor_score_floor = _mbr_donor_score_floor(
        donor_paths;
        donor_q_threshold = donor_q_threshold,
    )
    phase_started = time()
    @debug_l1 "MBR donor index starting: files=$(length(donor_paths)) score_floor=$(round(donor_score_floor, digits=4))"
    donor_files = _mbr_preintegration_donor_files(
        donor_paths,
        donor_score_floor,
    )
    @debug_l1 "MBR donor index complete: precursors=$(length(donor_files)) elapsed=$(round(time() - phase_started, digits=2))s"
    staged = _stage_mbr_integration_inputs!(
        PSMFileReference[ref for ref in candidate_refs if exists(ref)],
        output_folder,
        donor_files,
        q_value_threshold,
    )
    phase_started = last_progress = time()
    @debug_l1 "MBR staging sidecar cleanup starting: files=$(length(candidate_paths))"
    staged_paths = Set(
        file_path(ref) for ref in staged.integration_refs
    )
    for (file_idx, candidate_path) in enumerate(candidate_paths)
        if !(candidate_path in staged_paths)
            sidecar_path = candidate_path * PASS1_SIDECAR_SUFFIX
            isfile(sidecar_path) &&
                safeRm(sidecar_path; force = true)
        end
        if time() - last_progress >= 60
            @debug_l1 "MBR staging sidecar cleanup: files=$file_idx/$(length(candidate_paths)) elapsed=$(round(time() - phase_started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "MBR staging sidecar cleanup complete: elapsed=$(round(time() - phase_started, digits=2))s"
    @debug_l1 "Post-integration MBR staging: donor score floor=" *
              "$(round(donor_score_floor, digits=4)), " *
              "rows=$(staged.n_rows), candidates=$(staged.n_candidates)"
    return (
        integration_refs = staged.integration_refs,
        n_files = length(candidate_paths),
        n_rows = staged.n_rows,
        n_candidates = staged.n_candidates,
        donor_score_floor = donor_score_floor,
    )
end

"""
    _write_mbr_recovery_sidecars_from_candidates!(candidates, masks, n_rows, file_paths)

Write the per-file recovery sidecars from a candidates-only frame.

Non-candidate rows carry the defaults `apply_postintegration_mbr_rescoring!` assigns (`false` /
`NaN32`), so each file's columns are filled with those and the candidate values scattered in by
`masks[i]`. Candidates were concatenated in file order, so walking the masks in the same order
consumes them contiguously.
"""
function _write_mbr_recovery_sidecars_from_candidates!(
    candidates::DataFrame,
    masks::Vector{BitVector},
    n_rows::Vector{Int},
    file_paths::Vector{String},
)
    length(masks) == length(file_paths) == length(n_rows) ||
        error("MBR recovery mask/file count mismatch")
    cand_precursor = candidates[!, :precursor_idx]
    cand_scan      = candidates[!, :scan_idx]
    cand_recovered = candidates[!, :mbr_recovered]
    cand_flag      = candidates[!, :MBR_transfer_candidate]
    cand_prob      = candidates[!, :mbr_target_decoy_prob]
    cand_ftr_q     = candidates[!, :ftr_qval_true]
    cand_ftr_pep   = candidates[!, :ftr_pep_true]
    cand_tot_q     = candidates[!, :mbr_total_error_qval_true]
    cand_tot_r     = candidates[!, :mbr_total_error_rate_true]
    cand_cf_prob   = candidates[!, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN]
    cand_cf_idx    = candidates[!, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN]
    cursor = 0
    for (file_idx, path) in enumerate(file_paths)
        # Just the row identifiers, copied, so `path` is not left mapped (it is rewritten later).
        main = load_arrow_dataframe(path; cols = [:precursor_idx, :scan_idx])
        n = n_rows[file_idx]
        length(main.precursor_idx) == n ||
            error("MBR recovery row-count mismatch at $path")
        mask = masks[file_idx]
        recovered = falses(n)
        flag      = falses(n)
        prob      = fill(NaN32, n)
        ftr_q     = fill(NaN32, n)
        ftr_pep   = fill(NaN32, n)
        tot_q     = fill(NaN32, n)
        tot_r     = fill(NaN32, n)
        cf_prob   = fill(NaN32, n)
        cf_idx    = zeros(UInt8, n)
        @inbounds for row in 1:n
            mask[row] || continue
            cursor += 1
            (
                main.precursor_idx[row] == cand_precursor[cursor] &&
                main.scan_idx[row] == cand_scan[cursor]
            ) || error("MBR recovery candidate misalignment at row $row of $path")
            recovered[row] = cand_recovered[cursor]
            flag[row]      = cand_flag[cursor]
            prob[row]      = cand_prob[cursor]
            ftr_q[row]     = cand_ftr_q[cursor]
            ftr_pep[row]   = cand_ftr_pep[cursor]
            tot_q[row]     = cand_tot_q[cursor]
            tot_r[row]     = cand_tot_r[cursor]
            cf_prob[row]   = cand_cf_prob[cursor]
            cf_idx[row]    = cand_cf_idx[cursor]
        end
        writeArrow(path * RECOVERY_SIDECAR_SUFFIX, DataFrame(
            precursor_idx = UInt32.(main.precursor_idx),
            scan_idx = UInt32.(main.scan_idx),
            mbr_recovered = recovered,
            MBR_transfer_candidate = flag,
            mbr_target_decoy_prob = prob,
            ftr_qval_true = ftr_q,
            ftr_pep_true = ftr_pep,
            mbr_total_error_qval_true = tot_q,
            mbr_total_error_rate_true = tot_r,
            mbr_counterfactual_decoy_prob = cf_prob,
            mbr_counterfactual_decoy_index = cf_idx,
        ))
    end
    cursor == nrow(candidates) ||
        error("MBR recovery left $(nrow(candidates) - cursor) candidates unassigned")
    return length(file_paths)
end

function _write_mbr_recovery_sidecars!(
    frame::DataFrame,
    file_paths::Vector{String},
)
    offset = 0
    for path in file_paths
        main = Arrow.Table(path)
        n = length(main.precursor_idx)
        rows = (offset + 1):(offset + n)
        if n > 0
            @inbounds for (local_row, frame_row) in enumerate(rows)
                (
                    main.precursor_idx[local_row] ==
                        frame.precursor_idx[frame_row] &&
                    main.scan_idx[local_row] ==
                        frame.scan_idx[frame_row]
                ) || error("MBR recovery frame misalignment at $path")
            end
        end
        recovery = DataFrame(
            precursor_idx = UInt32.(frame.precursor_idx[rows]),
            scan_idx = UInt32.(frame.scan_idx[rows]),
            mbr_recovered = Bool.(frame.mbr_recovered[rows]),
            MBR_transfer_candidate =
                Bool.(frame.MBR_transfer_candidate[rows]),
            mbr_target_decoy_prob =
                Float32.(frame.mbr_target_decoy_prob[rows]),
            ftr_qval_true = Float32.(frame.ftr_qval_true[rows]),
            ftr_pep_true = Float32.(frame.ftr_pep_true[rows]),
            mbr_total_error_qval_true =
                Float32.(frame.mbr_total_error_qval_true[rows]),
            mbr_total_error_rate_true =
                Float32.(frame.mbr_total_error_rate_true[rows]),
            mbr_counterfactual_decoy_prob = Float32.(
                frame[rows, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN]
            ),
            mbr_counterfactual_decoy_index = UInt8.(
                frame[rows, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN]
            ),
        )
        writeArrow(path * RECOVERY_SIDECAR_SUFFIX, recovery)
        offset += n
    end
    offset == nrow(frame) ||
        error("MBR recovery frame contains unassigned rows")
    return length(file_paths)
end

# Score-remap bounds derived from the frozen pre-MBR q-value spline. Extracted so the fused path in
# _merge_mbr_recoveries! and the standalone _remap_mbr_scores! fallback compute them identically.
function _mbr_remap_bounds(qval_spline, q_value_threshold::Float32)
    score_ceiling = prevfloat(_score_floor_for_qvalue(
        qval_spline,
        q_value_threshold,
    ))
    score_floor = eps(Float32)
    return score_floor, max(score_ceiling - score_floor, 0.0f0)
end

"""
    _assert_recovery_aligned(main_pid, main_scan, rec_pid, rec_scan, n, path)

Check recovery-sidecar alignment using concrete column types in the row loop.
"""
@noinline function _assert_recovery_aligned(main_pid, main_scan, rec_pid, rec_scan,
                                            n::Int, path::AbstractString)
    @inbounds for row in 1:n
        (main_pid[row] == rec_pid[row] && main_scan[row] == rec_scan[row]) ||
            error("MBR recovery sidecar misalignment at row $row of $path")
    end
    return nothing
end

"""
    _mbr_keep_mask(qval, global_qval, recovered, n, q_value_threshold) -> BitVector

Rows to retain after MBR: those that passed the original q-value test, plus those MBR recovered.
Separate method as a function barrier -- the three columns arrive abstractly typed from the DataFrame.
"""
@noinline function _mbr_keep_mask(qval, global_qval, recovered, n::Int,
                                  q_value_threshold::Float32)
    keep = BitVector(undef, n)
    @inbounds for row in 1:n
        keep[row] = _mbr_initial_pass(qval[row], global_qval[row], q_value_threshold) ||
                    Bool(recovered[row])
    end
    return keep
end

"""
    _copy_sidecar_column!(dest, src) -> dest

Fill `dest` from `src`, converting elementwise. Replaces `collect(T.(src))`, which allocated twice --
the broadcast materialises one array and `collect` copies it again. Conversion semantics are unchanged,
including throwing on `missing`.
"""
@noinline function _copy_sidecar_column!(dest::Vector{T}, src) where {T}
    @inbounds for i in eachindex(dest)
        dest[i] = T(src[i])
    end
    return dest
end

function _merge_mbr_recoveries!(
    file_paths::Vector{String},
    q_value_threshold::Float32;
    remap_bounds::Union{Nothing, Tuple{Float32, Float32}} = nothing,
)
    for path in file_paths
        recovery_path = path * RECOVERY_SIDECAR_SUFFIX
        isfile(recovery_path) ||
            error("Missing MBR recovery sidecar at $recovery_path")
        main = load_arrow_dataframe(path)   # unmapped: `path` is rewritten below
        n = nrow(main)
        # Every recovery column is copied into `main`, so the sidecar is unmapped afterwards.
        with_arrow_table(recovery_path) do recovery
            length(recovery.precursor_idx) == n ||
                error("MBR recovery sidecar row-count mismatch at $recovery_path")
            # Resolved once, not once per row -- see _assert_recovery_aligned.
            _assert_recovery_aligned(main[!, :precursor_idx], main[!, :scan_idx],
                                     recovery.precursor_idx, recovery.scan_idx, n, path)

            main[!, :mbr_recovered] =
                _copy_sidecar_column!(Vector{Bool}(undef, n), recovery.mbr_recovered)
            main[!, :MBR_transfer_candidate] =
                _copy_sidecar_column!(Vector{Bool}(undef, n), recovery.MBR_transfer_candidate)
            main[!, :mbr_target_decoy_prob] =
                _copy_sidecar_column!(Vector{Float32}(undef, n), recovery.mbr_target_decoy_prob)
            main[!, :ftr_qval_true] =
                _copy_sidecar_column!(Vector{Float32}(undef, n), recovery.ftr_qval_true)
            main[!, :ftr_pep_true] =
                _copy_sidecar_column!(Vector{Float32}(undef, n), recovery.ftr_pep_true)
            main[!, :mbr_total_error_qval_true] =
                _copy_sidecar_column!(Vector{Float32}(undef, n), recovery.mbr_total_error_qval_true)
            main[!, :mbr_total_error_rate_true] =
                _copy_sidecar_column!(Vector{Float32}(undef, n), recovery.mbr_total_error_rate_true)
            main[!, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN] =
                _copy_sidecar_column!(
                    Vector{Float32}(undef, n),
                    Tables.getcolumn(recovery, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN),
                )
            main[!, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN] =
                _copy_sidecar_column!(
                    Vector{UInt8}(undef, n),
                    Tables.getcolumn(recovery, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN),
                )
        end

        keep = _mbr_keep_mask(main[!, :qval], main[!, :global_qval],
                              main[!, :mbr_recovered], n, q_value_threshold)
        # This table is privately owned, so compact its columns in place.
        deleteat!(main, .!keep)
        if remap_bounds !== nothing
            score_floor, width = remap_bounds
            recovered = main[!, :mbr_recovered]
            transfer = main[!, :mbr_target_decoy_prob]
            prec_prob = main[!, :prec_prob]
            counterfactual = main[!, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN]
            counterfactual_prec_prob = fill(NaN32, nrow(main))
            @inbounds for row in eachindex(recovered)
                if Bool(recovered[row])
                    transfer_score =
                        clamp(Float32(transfer[row]), 0.0f0, 1.0f0)
                    prec_prob[row] = score_floor + width * transfer_score
                end
                counterfactual_score = Float32(counterfactual[row])
                if isfinite(counterfactual_score)
                    counterfactual_prec_prob[row] = score_floor + width *
                        clamp(counterfactual_score, 0.0f0, 1.0f0)
                end
            end
            main[!, MBR_COUNTERFACTUAL_DECOY_PRECURSOR_PROB_COLUMN] =
                counterfactual_prec_prob
        end
        writeArrow(path, main)
    end
    return file_paths
end

function _remap_mbr_scores!(
    refs::Vector{PSMFileReference},
    merged_path::String;
    q_value_threshold::Float32,
    fdr_scale_factor::Float32,
    pre_mbr_qval_spline = nothing,
)
    qval_spline = pre_mbr_qval_spline
    if qval_spline === nothing
        spline_result = build_qvalue_spline_from_refs(
            refs,
            :prec_prob,
            merged_path;
            compute_pep = false,
            fdr_scale_factor = fdr_scale_factor,
            temp_prefix = "mbr_pre_remap",
        )
        spline_result === nothing && return refs
        qval_spline = spline_result.qval_spline
    end
    score_floor, width = _mbr_remap_bounds(qval_spline, q_value_threshold)

    for ref in refs
        path = file_path(ref)
        main = load_arrow_dataframe(path)   # unmapped: `path` is rewritten below
        hasproperty(main, :mbr_recovered) || continue
        has_counterfactual =
            hasproperty(main, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN)
        counterfactual_prec_prob = has_counterfactual ?
            fill(NaN32, nrow(main)) : Float32[]
        @inbounds for row in 1:nrow(main)
            if Bool(main.mbr_recovered[row])
                transfer_score = clamp(
                    Float32(main.mbr_target_decoy_prob[row]),
                    0.0f0,
                    1.0f0,
                )
                main.prec_prob[row] = score_floor + width * transfer_score
            end
            if has_counterfactual
                counterfactual_score = Float32(
                    main[row, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN]
                )
                if isfinite(counterfactual_score)
                    counterfactual_prec_prob[row] = score_floor + width *
                        clamp(counterfactual_score, 0.0f0, 1.0f0)
                end
            end
        end
        if has_counterfactual
            main[!, MBR_COUNTERFACTUAL_DECOY_PRECURSOR_PROB_COLUMN] =
                counterfactual_prec_prob
        end
        writeArrow(path, main)
    end
    return PSMFileReference[
        PSMFileReference(file_path(ref)) for ref in refs
    ]
end

# Sidecar tag for the deferred q-value columns. Deliberately NOT in _cleanup_mbr_sidecars!'s suffix
# list: the sidecar must survive until summarize_results! consumes it.
const MBR_QVAL_SIDECAR_TAG = "mbr_qval"

# register_sidecar! REFUSES a column that already exists in main ("Sidecar column qval already exists
# in main file schema") -- the sidecar mechanism is append-only by design. :qval and :pep are written
# by PrecursorScoringSearch long before this, so the recomputed values are staged under distinct
# names and copied over the originals when summarize_results! consolidates.
const MBR_QVAL_STAGED_COL = :qval_post_mbr
const MBR_PEP_STAGED_COL = :pep_post_mbr

function _recalculate_post_mbr_qvalues!(
    refs::Vector{PSMFileReference},
    merged_path::String;
    q_value_threshold::Float32,
    fdr_scale_factor::Float32,
)
    spline_result = build_qvalue_spline_from_refs(
        refs,
        :prec_prob,
        merged_path;
        compute_pep = true,
        fdr_scale_factor = fdr_scale_factor,
        temp_prefix = "post_integration_mbr",
    )
    spline_result === nothing && return refs, false

    # Stage calibrated scores in aligned sidecars. process_final_psms! applies the
    # row filter while writing final PSMs, avoiding another full-table rewrite.
    qval_spline = spline_result.qval_spline
    pep_interp = spline_result.pep_interp
    for ref in refs
        scores = materialize_columns(ref, Symbol[:prec_prob])[!, :prec_prob]
        n = length(scores)
        qvals = Vector{Float32}(undef, n)
        peps = Vector{Float32}(undef, n)
        @inbounds for row in 1:n
            score = Float32(scores[row])
            qvals[row] = Float32(qval_spline(score))
            peps[row] = Float32(pep_interp(score))
        end
        add_columns_via_sidecar!(
            ref,
            MBR_QVAL_STAGED_COL => qvals,
            MBR_PEP_STAGED_COL => peps;
            tag = MBR_QVAL_SIDECAR_TAG,
        )
    end
    # Preserve the sidecar registrations needed by summarize_results!.
    return refs, true
end

function _cleanup_mbr_sidecars!(file_paths::Vector{String})
    # PIONEER_MBR_KEEP_SIDECARS=1 preserves them so an optimisation can be checked for
    # byte-identity against a reference run. Off by default.
    get(ENV, "PIONEER_MBR_KEEP_SIDECARS", "0") == "1" && return nothing
    for path in file_paths
        for suffix in (
            PASS1_SIDECAR_SUFFIX,
            MBR_SIDECAR_SUFFIX,
            RECOVERY_SIDECAR_SUFFIX,
        )
            sidecar_path = path * suffix
            isfile(sidecar_path) &&
                safeRm(sidecar_path; force = true)
        end
    end
    return nothing
end

# Called while writing final PSMs to avoid a separate table rewrite.
function _drop_internal_mbr_columns!(main::DataFrame)
    internal_columns = Symbol[
        column for column in MBR_INTERNAL_INTEGRATED_COLUMNS
        if hasproperty(main, column)
    ]
    isempty(internal_columns) || select!(main, Not(internal_columns))
    return main
end

# Keep donor evidence and worker closures out of the subsequent fitting phase.
function _prepare_postintegration_mbr_features!(file_paths, precursors;
    run_similarity_atlas, q_value_threshold, donor_q_threshold,
    bitvec_rank_tables_by_file)
    phase_started = time()
    @debug_l1 "Post-integration MBR donor threshold starting: files=$(length(file_paths))"
    donor_score_floor = _mbr_donor_score_floor(
        file_paths;
        donor_q_threshold = donor_q_threshold,
        require_initial_pass = true,
        q_value_threshold = q_value_threshold,
    )
    @debug_l1 "Post-integration MBR donor threshold complete: elapsed=$(round(time() - phase_started, digits=2))s"
    phase_started = time()
    @debug_l1 "Post-integration MBR donor dictionary starting: files=$(length(file_paths))"
    donor_dict = build_mbr_integrated_donor_dict(
        file_paths,
        donor_score_floor;
        q_value_threshold = q_value_threshold,
    )
    @debug_l1 "Post-integration MBR donor dictionary complete: precursors=$(length(donor_dict)) entries=$(sum(length, values(donor_dict); init=0)) elapsed=$(round(time() - phase_started, digits=2))s"
    try
        phase_started = time()
        @debug_l1 "Post-integration MBR LOD thresholds starting: files=$(length(file_paths))"
        lod_thresholds = _mbr_lod_thresholds(
            file_paths,
            donor_score_floor;
            q_value_threshold = q_value_threshold,
        )
        @debug_l1 "Post-integration MBR LOD thresholds complete: elapsed=$(round(time() - phase_started, digits=2))s"
        phase_started = time()
        @debug_l1 "Post-integration MBR receiver clusters starting: files=$(length(file_paths))"
        receiver_run_clusters = build_mbr_receiver_run_clusters(
            file_paths;
            q_value_threshold = q_value_threshold,
        )
        @debug_l1 "Post-integration MBR receiver clusters complete: elapsed=$(round(time() - phase_started, digits=2))s"
        phase_started = time()
        @debug_l1 "Post-integration MBR partner pools starting: files=$(length(file_paths))"
        partner_pools = build_mbr_partner_pools(file_paths, precursors)
        @debug_l1 "Post-integration MBR partner pools complete: elapsed=$(round(time() - phase_started, digits=2))s"
        phase_started = time()
        @debug_l1 "Post-integration MBR counterfactual eligibility starting: files=$(length(file_paths))"
        eligibility = build_mbr_counterfactual_eligibility(
            file_paths;
            q_value_threshold = q_value_threshold,
        )
        @debug_l1 "Post-integration MBR counterfactual eligibility complete: elapsed=$(round(time() - phase_started, digits=2))s"
        feature_started = time()
        feature_progress = Ref((files = 0, rows = 0, candidates = 0,
                                selection_seconds = 0.0, feature_seconds = 0.0,
                                write_seconds = 0.0, logged_at = feature_started))
        feature_progress_lock = ReentrantLock()
        @debug_l1 "Post-integration MBR features starting: files=$(length(file_paths))"
        _run_files = function (chunk)
            for file_position in chunk
                path = file_paths[file_position]
                tbl = Arrow.Table(path)
                file_idx = isempty(tbl.ms_file_idx) ?
                    UInt32(0) :
                    UInt32(first(tbl.ms_file_idx))
                rank_table = bitvec_rank_tables_by_file === nothing ?
                    nothing :
                    get(bitvec_rank_tables_by_file, file_idx, nothing)
                stats = Ref{Any}()
                compute_postintegration_mbr_features!(
                    path,
                    donor_dict,
                    partner_pools,
                    eligibility;
                    run_similarity_atlas = run_similarity_atlas,
                    receiver_run_clusters = receiver_run_clusters,
                    lod_log2_weight_by_file = lod_thresholds.by_file,
                    lod_log2_weight_global = lod_thresholds.global_lod,
                    bitvec_rank_table = rank_table,
                    q_value_threshold = q_value_threshold,
                    stats = stats,
                )
                file_stats = stats[]
                lock(feature_progress_lock) do
                    state = feature_progress[]
                    now = time()
                    done = state.files + 1
                    rows = state.rows + file_stats.rows
                    candidates = state.candidates + file_stats.candidates
                    log_progress = now - state.logged_at >= 60
                    if log_progress
                        @debug_l1 "Post-integration MBR features: files=$done/$(length(file_paths)) rows=$rows candidates=$candidates elapsed=$(round(now - feature_started, digits=2))s"
                    end
                    feature_progress[] = (
                        files = done, rows = rows, candidates = candidates,
                        selection_seconds = state.selection_seconds + file_stats.selection_seconds,
                        feature_seconds = state.feature_seconds + file_stats.feature_seconds,
                        write_seconds = state.write_seconds + file_stats.write_seconds,
                        logged_at = log_progress ? now : state.logged_at,
                    )
                end
            end
        end
        parallel_foreach!(length(file_paths)) do chunk
            _run_files(chunk)
        end

        feature_totals = feature_progress[]
        @debug_l1 "Post-integration MBR features complete: files=$(length(file_paths)) rows=$(feature_totals.rows) candidates=$(feature_totals.candidates) elapsed=$(round(time() - feature_started, digits=2))s"
        @debug_l1 "Post-integration MBR features cumulative worker time: selection=$(round(feature_totals.selection_seconds, digits=2))s features=$(round(feature_totals.feature_seconds, digits=2))s write=$(round(feature_totals.write_seconds, digits=2))s"
    finally
        # The store is memory-mapped; release it once no receiver needs donors.
        close_mbr_donor_store!(donor_dict)
    end
    return nothing
end

"""
    finalize_postintegration_mbr!(integrated_paths, precursors; ...)

Build integrated donors, generate paired real/counterfactual evidence, train
the out-of-fold transfer model, and apply the combined precursor error budget.
"""
function finalize_postintegration_mbr!(
    integrated_paths::Vector{String},
    precursors::LibraryPrecursors;
    run_similarity_atlas::Union{Nothing, RunSimilarityAtlas},
    q_value_threshold::Float32,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
    fdr_scale_factor::Float32,
    ion_mobility::Bool = true,
    merged_path::String,
    pre_mbr_qval_spline = nothing,
    bitvec_rank_tables_by_file::Union{
        Nothing,
        Dict{UInt32, Vector{UInt16}},
    } = nothing,
)
    started = time()
    file_paths = String[
        path for path in integrated_paths
        if isfile(path) && isfile(path * PASS1_SIDECAR_SUFFIX)
    ]
    isempty(file_paths) && return (
        n_files = 0,
        n_candidates = 0,
        n_recovered = 0,
        base_targets = 0,
        base_decoys = 0,
        baseline_error_rate = 0.0f0,
        mbr_targets = 0,
        mbr_decoys = 0,
        mbr_false_transfers = 0,
        internal_ftr_targets = 0,
        internal_ftr_errors = 0,
        internal_ftr_estimate = NaN32,
        total_targets = 0,
        total_errors = 0,
        combined_error_rate = 0.0f0,
    )

    _prepare_postintegration_mbr_features!(file_paths, precursors;
        run_similarity_atlas, q_value_threshold, donor_q_threshold,
        bitvec_rank_tables_by_file)
    GC.gc()
    # Candidates only (~10% of rows). See load_postintegration_mbr_candidates for why this is safe.
    phase_started = time()
    @debug_l1 "Post-integration MBR candidate loading starting: files=$(length(file_paths))"
    loaded, frame, summary = mktempdir(dirname(first(file_paths)); prefix=".mbr_features_") do dir
        store = _MBRFeatureStore(joinpath(dir, "features.bin"))
        try
            loaded = load_postintegration_mbr_candidates(file_paths, q_value_threshold; feature_store=store)
            frame = loaded.candidates
            @debug_l1 "Post-integration MBR candidate loading complete: candidates=$(nrow(frame)) elapsed=$(round(time() - phase_started, digits=2))s"
            phase_started = time()
            @debug_l1 "Post-integration MBR rescoring starting: candidates=$(nrow(frame))"
            summary = apply_postintegration_mbr_rescoring!(frame;
                alpha=q_value_threshold, q_value_threshold,
                baseline_counts=(loaded.base_targets, loaded.base_decoys),
                frame_is_candidates=true, feature_source=store, ion_mobility)
            return loaded, frame, summary
        finally
            close(store)
        end
    end
    @debug_l1 "Post-integration MBR rescoring complete: elapsed=$(round(time() - phase_started, digits=2))s"
    phase_started = time()
    @debug_l1 "Post-integration MBR recovery sidecars starting: files=$(length(file_paths))"
    _write_mbr_recovery_sidecars_from_candidates!(
        frame, loaded.masks, loaded.n_rows, file_paths)
    frame = DataFrame()
    GC.gc(false)
    @debug_l1 "Post-integration MBR recovery sidecars complete: elapsed=$(round(time() - phase_started, digits=2))s"
    # When the caller supplies the frozen pre-MBR spline (it always does; see
    # PrecursorScoringSearch.jl:358 -> IntegrateChromatogramsSearch.jl:628) the score remap needs no
    # global pass, so it rides along with the merge instead of costing its own read-modify-write of
    # every file. The standalone path is kept for the documented spline === nothing contract.
    remap_bounds = pre_mbr_qval_spline === nothing ?
        nothing :
        _mbr_remap_bounds(pre_mbr_qval_spline, q_value_threshold)
    phase_started = time()
    @debug_l1 "Post-integration MBR recovery merge starting: files=$(length(file_paths))"
    _merge_mbr_recoveries!(
        file_paths,
        q_value_threshold;
        remap_bounds = remap_bounds,
    )

    @debug_l1 "Post-integration MBR recovery merge complete: elapsed=$(round(time() - phase_started, digits=2))s"
    refs = PSMFileReference[PSMFileReference(path) for path in file_paths]
    if remap_bounds === nothing
        phase_started = time()
        @debug_l1 "Post-integration MBR score remapping starting: files=$(length(refs))"
        refs = _remap_mbr_scores!(
            refs,
            merged_path;
            q_value_threshold = q_value_threshold,
            fdr_scale_factor = fdr_scale_factor,
            pre_mbr_qval_spline = pre_mbr_qval_spline,
        )
        @debug_l1 "Post-integration MBR score remapping complete: elapsed=$(round(time() - phase_started, digits=2))s"
    end
    phase_started = time()
    @debug_l1 "Post-integration MBR q-value recalculation starting: files=$(length(refs))"
    refs, qval_deferred = _recalculate_post_mbr_qvalues!(
        refs,
        merged_path;
        q_value_threshold = q_value_threshold,
        fdr_scale_factor = fdr_scale_factor,
    )
    @debug_l1 "Post-integration MBR q-value recalculation complete: elapsed=$(round(time() - phase_started, digits=2))s"
    phase_started = time()
    @debug_l1 "Post-integration MBR sidecar cleanup starting: files=$(length(file_paths))"
    _cleanup_mbr_sidecars!(file_paths)
    @debug_l1 "Post-integration MBR sidecar cleanup complete: elapsed=$(round(time() - phase_started, digits=2))s"
    @debug_l1 "Post-integration MBR finalization complete: files=$(length(file_paths)) elapsed=$(round(time() - started, digits=2))s"
    # refs carry the deferred :qval/:pep sidecars; qval_deferred says whether the q-value filter
    # still has to be applied downstream (false when no spline could be built).
    return merge(
        summary,
        (
            n_files = length(file_paths),
            mbr_refs = refs,
            qval_deferred = qval_deferred,
            qval_threshold = q_value_threshold,
        ),
    )
end
