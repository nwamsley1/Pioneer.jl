# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Donor pool construction for MBRRecoverySearch.
#
# Reads the post-PrecursorScoringSearch per-file PSM tables
# (temp_data/main_search_psms/{file}.arrow, ~117 cols including q_value,
# trace_prob, weight, log2_intensity_explained, irt_obs, irt_pred, rt,
# log_by_ratio_m0, frag{1..8}_smoothed_intensity, target, cv_fold).
#
# Emits a dict pid → Vector{_MBRDonorEntry} keeping the top
# `max_donors_per_pid` entries by trace_prob from distinct files. Reuses
# the existing _MBRDonorEntry struct from PrecursorScoringSearch/mbr_streaming.jl
# so the inline MBR feature compute can share semantics.

"""
    list_post_scoring_psm_files(search_context) -> Vector{String}

Lists the per-file PSM tables in `temp_data/main_search_psms/` that
PrecursorScoringSearch has annotated with per-file qval and trace_prob.
Filters out fold-split intermediates ({file}_fold0.arrow, {file}_fold1.arrow)
which are stripped after MainSearch's summarize_results! step.
"""
function list_post_scoring_psm_files(search_context::SearchContext)
    dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    isdir(dir) || return String[]
    files = readdir(dir, join=true)
    keep = filter(files) do f
        endswith(f, ".arrow") &&
        !occursin("_fold0", basename(f)) &&
        !occursin("_fold1", basename(f)) &&
        !endswith(basename(f), ".mbr_sidecar.arrow") &&
        !endswith(basename(f), ".pass1_sidecar.arrow") &&
        !endswith(basename(f), ".recovery_sidecar.arrow") &&
        !endswith(basename(f), RECOVERY_SEED_SIDECAR_SUFFIX)
    end
    sort!(keep)
    return keep
end

"""
    _read_smoothed_frag_sqrt(tbl, row_idx) -> NTuple{8, Float32}

Reads the eight per-row `frag{i}_smoothed_intensity` cells, L1-normalizes them,
and returns sqrt of the normalized fractions. Mirrors
`_mbr_smoothed_spectrum_sqrt_tuple` from mbr_streaming.jl but reads directly
from a row of an Arrow table.
"""
@inline function _read_smoothed_frag_sqrt(tbl, i::Integer)
    f1 = Float32(tbl.frag1_smoothed_intensity[i])
    f2 = Float32(tbl.frag2_smoothed_intensity[i])
    f3 = Float32(tbl.frag3_smoothed_intensity[i])
    f4 = Float32(tbl.frag4_smoothed_intensity[i])
    f5 = Float32(tbl.frag5_smoothed_intensity[i])
    f6 = Float32(tbl.frag6_smoothed_intensity[i])
    f7 = Float32(tbl.frag7_smoothed_intensity[i])
    f8 = Float32(tbl.frag8_smoothed_intensity[i])
    s  = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
    s > 0.0f0 || return MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    inv = 1.0f0 / s
    return (sqrt(f1 * inv), sqrt(f2 * inv), sqrt(f3 * inv), sqrt(f4 * inv),
            sqrt(f5 * inv), sqrt(f6 * inv), sqrt(f7 * inv), sqrt(f8 * inv))
end

"""
    build_recovery_donor_pool(psm_paths, spec_lib; qval_threshold, max_donors)
        -> Dict{UInt32, Vector{_MBRDonorEntry}}

Streams each per-file post-scoring PSM table once, collecting rows whose
q_value clears `qval_threshold`. For each pid, retains the top `max_donors`
entries by trace_prob from distinct files (sorted desc).

Includes both targets and decoys naturally — qval ≤ 0.01 lets through ~1 %
decoy donors which feed FTR's target/decoy balance through the paired
counterfactual machinery.
"""
function build_recovery_donor_pool(
    psm_paths::Vector{String},
    spec_lib::SpectralLibrary;
    qval_threshold::Float32,
    max_donors::Int,
)
    pool = Dict{UInt32, Vector{_MBRDonorEntry}}()
    n_donor_rows_seen = 0
    n_files_with_donors = 0
    precursors = getPrecursors(spec_lib)
    prec_irt_pred = getIrt(precursors)

    for (file_idx, path) in enumerate(psm_paths)
        tbl = Arrow.Table(path)
        n = length(tbl.precursor_idx)
        n == 0 && continue

        # The post-scoring main_search_psms file's :q_value column is not
        # populated (always 0.0) — recompute per-file qvals from trace_prob
        # + target via get_qvalues!. fdr_scale_factor=1.0 matches the
        # library FDR scale we saw earlier (3.4M target / 3.4M decoy).
        probs = Float32.(tbl.trace_prob)
        labels = Vector{Bool}(tbl.target)
        qvals = zeros(Float32, n)
        get_qvalues!(probs, labels, qvals; doSort=true, fdr_scale_factor=1.0f0)

        file_had_donor = false
        @inbounds for i in 1:n
            qvals[i] <= qval_threshold || continue
            pid = UInt32(tbl.precursor_idx[i])
            entry = _MBRDonorEntry(
                probs[i],
                pid,
                Float32(tbl.weight[i]),
                Float32(tbl.log2_intensity_explained[i]),
                Float32(prec_irt_pred[pid] - tbl.irt_obs[i]),  # irt_residual
                Float32(tbl.irt_obs[i]),
                Float32(tbl.log_by_ratio_m0[i]),
                Float32(tbl.rt[i]),
                _read_smoothed_frag_sqrt(tbl, i),
                UInt32(file_idx),
                Bool(!tbl.target[i]),
            )
            entries = get!(() -> _MBRDonorEntry[], pool, pid)
            _insert_sorted_donor_entry!(entries, entry; max_kept=max_donors)
            n_donor_rows_seen += 1
            file_had_donor = true
        end
        file_had_donor && (n_files_with_donors += 1)
    end

    @user_info "MBRRecovery donor pool: $(length(pool)) unique pids " *
               "from $n_donor_rows_seen rows in $n_files_with_donors files " *
               "(qval ≤ $qval_threshold)"
    return pool
end

"""
    main_search_pid_set(path) -> Set{UInt32}

Returns the set of precursor_idx values present in a single post-scoring PSM
table. Used to skip cells where the donor pid was already natively identified
in the receiver file.
"""
function main_search_pid_set(path::String)
    tbl = Arrow.Table(path)
    s = Set{UInt32}()
    sizehint!(s, length(tbl.precursor_idx))
    @inbounds for pid in tbl.precursor_idx
        push!(s, UInt32(pid))
    end
    return s
end
