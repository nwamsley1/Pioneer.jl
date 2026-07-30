# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

function _write_mbr_pass1_sidecars_from_main!(
    file_paths::Vector{String},
)
    n_written = 0
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
    end
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

# Takes REFS rather than paths: PrecursorScoringSearch now attaches :qval/:global_qval/:global_prob/
# :global_pep/:pep as a row-aligned sidecar instead of materialising precursor_scored_psms, so the
# candidate tables only carry them via their PSMFileReference. This function already rewrites every
# row it keeps, so consolidating the sidecar here costs nothing extra.
function _stage_mbr_integration_inputs!(
    candidate_refs::Vector{PSMFileReference},
    output_folder::String,
    donor_files::Dict{UInt32, Vector{Tuple{UInt32, Float32}}},
    q_value_threshold::Float32,
)
    mkpath(output_folder)
    refs = PSMFileReference[]
    n_rows = 0
    n_candidates = 0
    for candidate_ref in candidate_refs
        path = file_path(candidate_ref)
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing MBR Pass-1 sidecar at $pass1_path")
        main = load_with_sidecars(candidate_ref)
        pass1 = DataFrame(Tables.columntable(Arrow.Table(pass1_path)))
        nrow(main) == nrow(pass1) ||
            error("MBR Pass-1 sidecar row-count mismatch at $pass1_path")

        keep = falses(nrow(main))
        @inbounds for row in 1:nrow(main)
            baseline = _mbr_initial_pass(
                main.qval[row],
                main.global_qval[row],
                q_value_threshold,
            )
            global_qval = Float32(main.global_qval[row])
            global_pass =
                isfinite(global_qval) &&
                global_qval <= q_value_threshold
            run_qval = Float32(main.qval[row])
            run_pass =
                isfinite(run_qval) &&
                run_qval <= q_value_threshold
            candidate =
                global_pass && !run_pass &&
                _mbr_has_cross_run_donor(
                    donor_files,
                    UInt32(main.precursor_idx[row]),
                    UInt32(main.ms_file_idx[row]),
                )
            keep[row] = baseline || candidate
            n_candidates += candidate
        end

        staged_path = joinpath(output_folder, basename(path))
        writeArrow(staged_path, main[keep, :])
        writeArrow(staged_path * PASS1_SIDECAR_SUFFIX, pass1[keep, :])
        push!(refs, PSMFileReference(staged_path))
        n_rows += count(keep)
    end
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
    donor_files = _mbr_preintegration_donor_files(
        donor_paths,
        donor_score_floor,
    )
    staged = _stage_mbr_integration_inputs!(
        PSMFileReference[ref for ref in candidate_refs if exists(ref)],
        output_folder,
        donor_files,
        q_value_threshold,
    )
    staged_paths = Set(
        file_path(ref) for ref in staged.integration_refs
    )
    for candidate_path in candidate_paths
        candidate_path in staged_paths && continue
        sidecar_path = candidate_path * PASS1_SIDECAR_SUFFIX
        isfile(sidecar_path) &&
            safeRm(sidecar_path, nothing; force = true)
    end
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
        )
        writeArrow(path * RECOVERY_SIDECAR_SUFFIX, recovery)
        offset += n
    end
    offset == nrow(frame) ||
        error("MBR recovery frame contains unassigned rows")
    return length(file_paths)
end

function _merge_mbr_recoveries!(
    file_paths::Vector{String},
    q_value_threshold::Float32,
)
    for path in file_paths
        recovery_path = path * RECOVERY_SIDECAR_SUFFIX
        isfile(recovery_path) ||
            error("Missing MBR recovery sidecar at $recovery_path")
        main = DataFrame(Tables.columntable(Arrow.Table(path)))
        recovery = Arrow.Table(recovery_path)
        n = nrow(main)
        length(recovery.precursor_idx) == n ||
            error("MBR recovery sidecar row-count mismatch at $recovery_path")
        @inbounds for row in 1:n
            (
                main.precursor_idx[row] == recovery.precursor_idx[row] &&
                main.scan_idx[row] == recovery.scan_idx[row]
            ) || error("MBR recovery sidecar misalignment at row $row of $path")
        end

        main[!, :mbr_recovered] = collect(Bool.(recovery.mbr_recovered))
        main[!, :MBR_transfer_candidate] =
            collect(Bool.(recovery.MBR_transfer_candidate))
        main[!, :mbr_target_decoy_prob] =
            collect(Float32.(recovery.mbr_target_decoy_prob))
        main[!, :ftr_qval_true] =
            collect(Float32.(recovery.ftr_qval_true))
        main[!, :ftr_pep_true] =
            collect(Float32.(recovery.ftr_pep_true))
        main[!, :mbr_total_error_qval_true] =
            collect(Float32.(recovery.mbr_total_error_qval_true))
        main[!, :mbr_total_error_rate_true] =
            collect(Float32.(recovery.mbr_total_error_rate_true))

        keep = BitVector(undef, n)
        @inbounds for row in 1:n
            keep[row] =
                _mbr_initial_pass(
                    main.qval[row],
                    main.global_qval[row],
                    q_value_threshold,
                ) ||
                Bool(main.mbr_recovered[row])
        end
        writeArrow(path, main[keep, :])
    end
    return file_paths
end

function _remap_mbr_scores!(
    refs::Vector{PSMFileReference},
    merged_path::String;
    q_value_threshold::Float32,
    min_pep_points_per_bin::Int,
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
            min_pep_points_per_bin = min_pep_points_per_bin,
            fdr_scale_factor = fdr_scale_factor,
            temp_prefix = "mbr_pre_remap",
        )
        spline_result === nothing && return refs
        qval_spline = spline_result.qval_spline
    end
    score_ceiling = prevfloat(_score_floor_for_qvalue(
        qval_spline,
        q_value_threshold,
    ))
    score_floor = eps(Float32)
    width = max(score_ceiling - score_floor, 0.0f0)

    for ref in refs
        path = file_path(ref)
        main = DataFrame(Tables.columntable(Arrow.Table(path)))
        hasproperty(main, :mbr_recovered) || continue
        @inbounds for row in 1:nrow(main)
            Bool(main.mbr_recovered[row]) || continue
            transfer_score = clamp(
                Float32(main.mbr_target_decoy_prob[row]),
                0.0f0,
                1.0f0,
            )
            main.prec_prob[row] = score_floor + width * transfer_score
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
    min_pep_points_per_bin::Int,
    fdr_scale_factor::Float32,
)
    spline_result = build_qvalue_spline_from_refs(
        refs,
        :prec_prob,
        merged_path;
        compute_pep = true,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
        temp_prefix = "post_integration_mbr",
    )
    spline_result === nothing && return refs, false

    # This used to be an apply_pipeline! over every ref, which loads all 147 columns of each file and
    # rewrites them (measured 107.9 MB per file; 6 files x load+write = 1.26 GB) purely to append two
    # Float32 columns and drop rows. :qval and :pep are NEW columns, so they can ride in a row-aligned
    # sidecar instead -- and the row filter is deferred to the process_final_psms! loop in
    # summarize_results!, which already loads and rewrites each table.
    #
    # Order note: the old pipeline was add(:qval) -> filter -> add(:pep), so :pep was computed only on
    # surviving rows. Interpolation is pointwise, so computing it for every row and filtering
    # afterwards gives identical values for the rows that survive.
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
    # Return the SAME refs: rebuilding them (as this previously did) would discard the sidecar
    # registration that summarize_results! depends on.
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
                safeRm(sidecar_path, nothing; force = true)
        end
    end
    return nothing
end

# Drops the internal MBR evidence columns from an already-loaded table. This used to own its own
# full materialise-and-rewrite pass over every file; it is now called from the process_final_psms!
# loop in summarize_results!, which already reads and writes each table.
function _drop_internal_mbr_columns!(main::DataFrame)
    internal_columns = Symbol[
        column for column in MBR_INTERNAL_INTEGRATED_COLUMNS
        if hasproperty(main, column)
    ]
    isempty(internal_columns) || select!(main, Not(internal_columns))
    return main
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
    min_pep_points_per_bin::Int,
    fdr_scale_factor::Float32,
    merged_path::String,
    pre_mbr_qval_spline = nothing,
    bitvec_rank_tables_by_file::Union{
        Nothing,
        Dict{UInt32, Vector{UInt16}},
    } = nothing,
)
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

    # DIAGNOSTIC (PIONEER_MBR_PHASE_DIAG=1). Measured on Olsen 6-file: this function is ~60 GB /
    # ~36 s of MBR's ~73 GB / ~49 s total cost — the bulk — and it runs once inside
    # summarize_results!, outside any per-file instrumentation. Phase probes below locate it.
    _fdiag = get(ENV, "PIONEER_MBR_PHASE_DIAG", "0") == "1"
    _fstate = Ref((time(), Base.gc_bytes()))
    _mark = function (key::Symbol)
        _fdiag || return nothing
        t0, a0 = _fstate[]
        MBR_FINAL_DIAG[Symbol(key, :_bytes)] =
            get(MBR_FINAL_DIAG, Symbol(key, :_bytes), 0) + (Base.gc_bytes() - a0)
        MBR_FINAL_DIAG[Symbol(key, :_ms)] =
            get(MBR_FINAL_DIAG, Symbol(key, :_ms), 0) + round(Int, (time() - t0) * 1000)
        _fstate[] = (time(), Base.gc_bytes())
        return nothing
    end

    donor_score_floor = _mbr_donor_score_floor(
        file_paths;
        donor_q_threshold = donor_q_threshold,
        require_initial_pass = true,
        q_value_threshold = q_value_threshold,
    )
    _mark(:donor_floor)
    donor_dict = build_mbr_integrated_donor_dict(
        file_paths,
        donor_score_floor;
        q_value_threshold = q_value_threshold,
    )
    _mark(:donor_dict)
    lod_thresholds = _mbr_lod_thresholds(
        file_paths,
        donor_score_floor;
        q_value_threshold = q_value_threshold,
    )
    _mark(:lod_thresholds)
    receiver_run_clusters = build_mbr_receiver_run_clusters(
        file_paths;
        q_value_threshold = q_value_threshold,
    )
    _mark(:run_clusters)
    partner_pools = build_mbr_partner_pools(file_paths, precursors)
    _mark(:partner_pools)
    eligibility = build_mbr_counterfactual_eligibility(
        file_paths;
        q_value_threshold = q_value_threshold,
    )
    @debug_l1 "Post-integration MBR donors: precursors=$(length(donor_dict)), " *
              "entries=$(sum(length, values(donor_dict)))"

    _mark(:eligibility)
    # Base.gc_bytes() is process-global, so the row-loop probes inside
    # compute_postintegration_mbr_features! are meaningless when files run concurrently — each
    # thread's window absorbs every other thread's allocations. Run serially when diagnosing.
    _serial_diag = get(ENV, "PIONEER_MBR_ROW_DIAG", "0") == "1"
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
            )
        end
    end
    if _serial_diag
        _run_files(1:length(file_paths))
    else
        parallel_foreach!(length(file_paths)) do chunk
            _run_files(chunk)
        end
    end

    _mark(:compute_features)
    frame = load_postintegration_mbr_frame(file_paths)
    _mark(:load_frame)
    summary = apply_postintegration_mbr_rescoring!(
        frame;
        alpha = q_value_threshold,
        q_value_threshold = q_value_threshold,
    )
    _mark(:rescoring)
    _write_mbr_recovery_sidecars!(frame, file_paths)
    frame = DataFrame()
    GC.gc(false)
    _mark(:write_sidecars)
    _merge_mbr_recoveries!(file_paths, q_value_threshold)

    _mark(:merge_recoveries)
    refs = PSMFileReference[PSMFileReference(path) for path in file_paths]
    refs = _remap_mbr_scores!(
        refs,
        merged_path;
        q_value_threshold = q_value_threshold,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
        pre_mbr_qval_spline = pre_mbr_qval_spline,
    )
    _mark(:remap_scores)
    refs, qval_deferred = _recalculate_post_mbr_qvalues!(
        refs,
        merged_path;
        q_value_threshold = q_value_threshold,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
    )
    _mark(:recalc_qvalues)
    # The internal-column drop used to be its own full materialise-and-rewrite pass over every file
    # (_drop_internal_mbr_columns!). It is pure per-file work with no cross-file dependency, so it is
    # now folded into the process_final_psms! loop in summarize_results!, which already reads and
    # writes each table. One fewer full pass; peak memory is unchanged (still one file at a time).
    _cleanup_mbr_sidecars!(file_paths)
    _mark(:cleanup)
    _fdiag && _mbr_final_diag_report()
    # refs carry the deferred :qval/:pep sidecars; qval_deferred says whether the q-value filter
    # still has to be applied downstream (false when no spline could be built, matching the old
    # behaviour of leaving rows unfiltered in that case).
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


const MBR_FINAL_DIAG = Dict{Symbol, Int}()

function _mbr_final_diag_report()
    d = MBR_FINAL_DIAG
    keys_ordered = [:donor_floor, :donor_dict, :lod_thresholds, :run_clusters, :partner_pools,
                    :eligibility, :compute_features, :load_frame, :rescoring, :write_sidecars,
                    :merge_recoveries, :remap_scores, :recalc_qvalues, :cleanup]
    tot = sum(get(d, Symbol(k, :_bytes), 0) for k in keys_ordered)
    totms = sum(get(d, Symbol(k, :_ms), 0) for k in keys_ordered)
    gb(x) = round(x / 2^30, digits = 2)
    lines = ["finalize_postintegration_mbr! phase diagnostic:"]
    for k in keys_ordered
        b = get(d, Symbol(k, :_bytes), 0); m = get(d, Symbol(k, :_ms), 0)
        (b == 0 && m == 0) && continue
        pc = tot > 0 ? round(100 * b / tot, digits = 1) : 0.0
        rate = m > 0 ? round((b / 2^30) / (m / 1000), digits = 2) : 0.0
        push!(lines, "  $(rpad(string(k), 22)) $(lpad(gb(b), 7)) GB  $(lpad(m, 7)) ms  ($(lpad(pc, 5))%)  $(rate) GB/s")
    end
    push!(lines, "  $(rpad("TOTAL", 22)) $(lpad(gb(tot), 7)) GB  $(lpad(totms, 7)) ms")
    @user_info join(lines, "\n")
    return nothing
end
