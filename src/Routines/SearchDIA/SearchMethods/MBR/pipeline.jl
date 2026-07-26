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

function _stage_mbr_integration_inputs!(
    candidate_paths::Vector{String},
    output_folder::String,
    donor_files::Dict{UInt32, Vector{Tuple{UInt32, Float32}}},
    q_value_threshold::Float32,
)
    mkpath(output_folder)
    refs = PSMFileReference[]
    n_rows = 0
    n_candidates = 0
    for path in candidate_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing MBR Pass-1 sidecar at $pass1_path")
        main = DataFrame(Tables.columntable(Arrow.Table(path)))
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
        candidate_paths,
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
    spline_result === nothing && return refs
    pipeline = TransformPipeline() |>
        add_interpolated_column(
            :qval,
            :prec_prob,
            spline_result.qval_spline,
        ) |>
        filter_by_multiple_thresholds([
            (:qval, q_value_threshold),
        ]) |>
        add_interpolated_column(
            :pep,
            :prec_prob,
            spline_result.pep_interp,
        )
    apply_pipeline!(refs, pipeline; parallel = false)
    return PSMFileReference[
        PSMFileReference(file_path(ref)) for ref in refs
    ]
end

function _cleanup_mbr_sidecars!(file_paths::Vector{String})
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

function _drop_internal_mbr_columns!(file_paths::Vector{String})
    for path in file_paths
        main = DataFrame(Tables.columntable(Arrow.Table(path)))
        internal_columns = Symbol[
            column for column in MBR_INTERNAL_INTEGRATED_COLUMNS
            if hasproperty(main, column)
        ]
        isempty(internal_columns) || select!(main, Not(internal_columns))
        writeArrow(path, main)
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

    donor_score_floor = _mbr_donor_score_floor(
        file_paths;
        donor_q_threshold = donor_q_threshold,
        require_initial_pass = true,
        q_value_threshold = q_value_threshold,
    )
    donor_dict = build_mbr_integrated_donor_dict(
        file_paths,
        donor_score_floor;
        q_value_threshold = q_value_threshold,
    )
    partner_pools = build_mbr_partner_pools(file_paths, precursors)
    eligibility = build_mbr_counterfactual_eligibility(
        file_paths;
        q_value_threshold = q_value_threshold,
    )
    @debug_l1 "Post-integration MBR donors: precursors=$(length(donor_dict)), " *
              "entries=$(sum(length, values(donor_dict)))"

    parallel_foreach!(length(file_paths)) do chunk
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
                bitvec_rank_table = rank_table,
            )
        end
    end

    frame = load_postintegration_mbr_frame(file_paths)
    summary = apply_postintegration_mbr_rescoring!(
        frame;
        alpha = q_value_threshold,
        q_value_threshold = q_value_threshold,
    )
    _write_mbr_recovery_sidecars!(frame, file_paths)
    frame = DataFrame()
    GC.gc(false)
    _merge_mbr_recoveries!(file_paths, q_value_threshold)

    refs = PSMFileReference[PSMFileReference(path) for path in file_paths]
    refs = _remap_mbr_scores!(
        refs,
        merged_path;
        q_value_threshold = q_value_threshold,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
        pre_mbr_qval_spline = pre_mbr_qval_spline,
    )
    refs = _recalculate_post_mbr_qvalues!(
        refs,
        merged_path;
        q_value_threshold = q_value_threshold,
        min_pep_points_per_bin = min_pep_points_per_bin,
        fdr_scale_factor = fdr_scale_factor,
    )
    _drop_internal_mbr_columns!(file_paths)
    _cleanup_mbr_sidecars!(file_paths)
    return merge(summary, (n_files = length(file_paths),))
end
