# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR Phase 2: cross-run donor tracking + per-PSM MBR features.
#
# After a first-pass LightGBM produces `:trace_prob_prepass`, this module
# scans all PSMs and builds a per-`pair_id` top-2 donor dict. Top-2 is
# necessary because for a row in file=k, the donor must come from a different
# file; if the top-1 is from file=k, we fall back to top-2.
#
# Computes four cross-run features (one per row):
#   :MBR_max_pair_prob          — donor.trace_prob_prepass
#   :MBR_log2_weight_ratio      — log2(this.weight / donor.weight)
#   :MBR_log2_explained_ratio   — this.log2_intensity_explained − donor.log2_intensity_explained
#   :MBR_is_missing             — true when no donor exists in another file
# Plus a label-side flag used by the FTR controller in Phase 4:
#   :MBR_is_best_decoy          — donor.is_decoy
#
# v0.6.6's MBR_rv_coefficient is intentionally skipped here — it requires
# per-row chromatogram vectors (weights, irts) that develop's MainSearch
# reduction discards before second_pass_psms is written. Phase 4 may revisit.

using Statistics

struct _MBRDonorEntry
    trace_prob::Float32
    weight::Float32
    log2_intensity_explained::Float32
    irt_residual::Float32  # irt_pred − irt_obs of the donor row
    irt_obs::Float32       # raw observed iRT of the donor row
    log_by_ratio::Float32  # log(b_int+1) - log(y_int+1) of the donor row
    rt_obs::Float32        # literal scan RT (minutes) of the donor row
    ms_file_idx::UInt32
    is_decoy::Bool
end
const _MBR_DONOR_SENTINEL = _MBRDonorEntry(typemin(Float32), 0f0, 0f0, 0f0, 0f0, 0f0,
                                            0f0, UInt32(0), false)

"""
    compute_mbr_features!(psms; score_col=:trace_prob_prepass) -> NamedTuple

Compute 4 MBR features and `:MBR_is_best_decoy` on `psms` in-place. The
donor for each row at (pair_id, file=k) is the top-prepass-score entry
across files ≠ k for the same pair_id (top-2 tracked for fallback).

Sentinel values when no donor exists:
- numeric features: `-1.0f0`
- `:MBR_is_missing` = true
- `:MBR_is_best_decoy` = false

Returns a NamedTuple with diagnostic counts.
"""
function compute_mbr_features!(
    psms::DataFrame;
    score_col::Symbol = :trace_prob_prepass,
)
    t0 = time()

    pair_id_col = psms[!, :pair_id]::Vector{UInt32}
    score_v     = Float32.(psms[!, score_col])
    weight_v    = Float32.(psms[!, :weight])
    log2ie_v    = Float32.(psms[!, :log2_intensity_explained])
    file_v      = UInt32.(psms[!, :ms_file_idx])
    decoy_v     = Bool.(psms[!, :decoy])
    has_irt     = hasproperty(psms, :irt_pred) && hasproperty(psms, :irt_obs)
    irt_res_v   = if has_irt
        Float32.(psms[!, :irt_pred]) .- Float32.(psms[!, :irt_obs])
    else
        zeros(Float32, nrow(psms))
    end
    n = nrow(psms)

    # ── Build top-2 donor entries per pair_id ──
    top2 = Dict{UInt32, NTuple{2, _MBRDonorEntry}}()
    sizehint!(top2, n)
    @inbounds for i in 1:n
        pid = pair_id_col[i]
        e = _MBRDonorEntry(score_v[i], weight_v[i], log2ie_v[i], irt_res_v[i],
                           0f0, 0f0, 0f0,
                           file_v[i], decoy_v[i])
        cur = get(top2, pid, (_MBR_DONOR_SENTINEL, _MBR_DONOR_SENTINEL))
        if e.trace_prob > cur[1].trace_prob
            top2[pid] = (e, cur[1])
        elseif e.trace_prob > cur[2].trace_prob
            top2[pid] = (cur[1], e)
        end
    end

    # ── Compute features per row ──
    out_max_pair_prob = fill(-1f0, n)
    out_log2_w_ratio  = fill(-1f0, n)
    out_log2_e_ratio  = fill(-1f0, n)
    out_best_irt_diff = fill(-1f0, n)
    out_is_missing    = trues(n)
    out_is_best_decoy = falses(n)

    n_donor_present = 0
    @inbounds for i in 1:n
        pid = pair_id_col[i]
        my_file = file_v[i]
        cur = get(top2, pid, nothing)
        cur === nothing && continue
        donor = if cur[1].trace_prob > typemin(Float32) && cur[1].ms_file_idx != my_file
            cur[1]
        elseif cur[2].trace_prob > typemin(Float32) && cur[2].ms_file_idx != my_file
            cur[2]
        else
            nothing
        end
        donor === nothing && continue

        n_donor_present += 1
        out_max_pair_prob[i]   = donor.trace_prob
        if donor.weight > 0 && weight_v[i] > 0
            out_log2_w_ratio[i] = log2(weight_v[i] / donor.weight)
        end
        out_log2_e_ratio[i]    = log2ie_v[i] - donor.log2_intensity_explained
        if has_irt
            out_best_irt_diff[i] = abs(irt_res_v[i] - donor.irt_residual)
        end
        out_is_missing[i]      = false
        out_is_best_decoy[i]   = donor.is_decoy
    end

    psms[!, :MBR_max_pair_prob]      = out_max_pair_prob
    psms[!, :MBR_log2_weight_ratio]  = out_log2_w_ratio
    psms[!, :MBR_log2_explained_ratio] = out_log2_e_ratio
    psms[!, :MBR_best_irt_diff]      = out_best_irt_diff
    psms[!, :MBR_is_missing]         = out_is_missing
    psms[!, :MBR_is_best_decoy]      = out_is_best_decoy

    elapsed = time() - t0

    @user_info "MBR Phase 2 — donor features:"
    @user_info "  unique pair_ids tracked: $(length(top2))"
    pct = round(100 * n_donor_present / max(n, 1), digits=1)
    @user_info "  PSMs with cross-run donor: $n_donor_present / $n ($pct%)"
    if n_donor_present > 0
        valid = .!out_is_missing
        med_prob = median(out_max_pair_prob[valid])
        n_decoy_donor = sum(out_is_best_decoy)
        @user_info "  MBR_max_pair_prob median (donors found): $(round(med_prob, digits=4))"
        @user_info "  donor was a decoy on $n_decoy_donor rows"
    end
    @user_info "  compute_mbr_features! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_pair_ids       = length(top2),
        n_with_donor     = n_donor_present,
        n_missing        = n - n_donor_present,
        elapsed_s        = elapsed,
    )
end

"""
    compute_mbr_features_dual!(psms; score_col=:trace_prob_prepass) -> NamedTuple

MBR Batch F: compute parallel `_true` and `_false` MBR features per row.

For each row at (precursor_idx=p, file=k, counterfactual_partner_pid=q):
- `*_true`  donor = top score for precursor `p` in files ≠ k (same as
  the original Phase-2 design but keyed on precursor_idx, not pair_id).
- `*_false` donor = top score for precursor `q` in files ≠ k.

Both donors track top-2 entries for same-file fallback. Sentinels for
missing donors (numeric = -1.0, `*_is_missing` = true).

Adds 10 columns to `psms`:
- `MBR_max_pair_prob_true`,        `MBR_max_pair_prob_false`
- `MBR_log2_weight_ratio_true`,    `MBR_log2_weight_ratio_false`
- `MBR_log2_explained_ratio_true`, `MBR_log2_explained_ratio_false`
- `MBR_best_irt_diff_true`,        `MBR_best_irt_diff_false`
- `MBR_is_missing_true`,           `MBR_is_missing_false`

Requires `:counterfactual_partner_pid` populated by
`regenerate_counterfactual_partners!`.
"""
function compute_mbr_features_dual!(
    psms::DataFrame;
    score_col::Symbol = :trace_prob_prepass,
)
    t0 = time()

    pid_col     = psms[!, :precursor_idx]::Vector{UInt32}
    partner_col = psms[!, :counterfactual_partner_pid]::Vector{UInt32}
    score_v     = Float32.(psms[!, score_col])
    weight_v    = Float32.(psms[!, :weight])
    log2ie_v    = Float32.(psms[!, :log2_intensity_explained])
    file_v      = UInt32.(psms[!, :ms_file_idx])
    has_irt     = hasproperty(psms, :irt_pred) && hasproperty(psms, :irt_obs)
    irt_obs_v   = has_irt ? Float32.(psms[!, :irt_obs]) : zeros(Float32, nrow(psms))
    irt_res_v   = has_irt ?
        (Float32.(psms[!, :irt_pred]) .- irt_obs_v) :
        zeros(Float32, nrow(psms))
    has_logby   = hasproperty(psms, :log_by_ratio_m0)
    log_by_v    = has_logby ? Float32.(psms[!, :log_by_ratio_m0]) : zeros(Float32, nrow(psms))
    has_rt      = hasproperty(psms, :rt)
    rt_obs_v    = has_rt ? Float32.(psms[!, :rt]) : zeros(Float32, nrow(psms))
    n = nrow(psms)

    # Top-2 per precursor_idx is enough for the same-file-fallback case in
    # `_donor_for`: each (pid, file) PSM is unique after MainSearch's
    # best-per-precursor selection + paircomp, so per pid the top-1 and top-2
    # entries are guaranteed to be from different files. K=2 covers the rare
    # case where the global best donor is from the row's own file.
    #
    # Earlier the code tracked top-K donors (K = ceil(sqrt(N_files))) and
    # exported MBR_top_n_median_score_* / MBR_top_n_irt_diff_* features (median
    # over the top-K). 2026-05-16 A/B (PIONEER_MBR_DONOR_MODE=topk vs top1) on
    # Olsen 23-file + MTAC 6-file (with Astral entrap library) showed the
    # top-K aggregation has no material effect: |ΔIDs|≤0.7%, ΔPGs marginally
    # +0.06–0.25% in top-1's favor, entrapment EFDR MAE byte-identical.
    # Median features removed; top-K aggregation collapsed to K=2.
    K = 2
    all_entries = Dict{UInt32, Vector{_MBRDonorEntry}}()
    sizehint!(all_entries, max(64, n ÷ 8))
    @inbounds for i in 1:n
        pid = pid_col[i]
        e = _MBRDonorEntry(score_v[i], weight_v[i], log2ie_v[i], irt_res_v[i],
                           irt_obs_v[i], log_by_v[i], rt_obs_v[i],
                           file_v[i], false)
        push!(get!(() -> _MBRDonorEntry[], all_entries, pid), e)
    end
    # Sort each pid's entries by trace_prob desc, keep top-K (= 2).
    @inbounds for (_pid, entries) in all_entries
        sort!(entries; by = e -> -e.trace_prob)
        if length(entries) > K
            resize!(entries, K)
        end
    end

    # Helper: top-prepass donor for pid with file != my_file.
    @inline function _donor_for(pid::UInt32, my_file::UInt32)
        entries = get(all_entries, pid, nothing)
        entries === nothing && return nothing
        @inbounds for e in entries
            if e.ms_file_idx != my_file
                return e
            end
        end
        return nothing
    end

    out_max_t = fill(-1f0, n); out_max_f = fill(-1f0, n)
    out_lw_t  = fill(-1f0, n); out_lw_f  = fill(-1f0, n)
    out_le_t  = fill(-1f0, n); out_le_f  = fill(-1f0, n)
    out_ir_t  = fill(-1f0, n); out_ir_f  = fill(-1f0, n)
    out_miss_t = trues(n);     out_miss_f = trues(n)
    out_log_by_t    = fill(-1f0, n); out_log_by_f    = fill(-1f0, n)
    # rtv1 (2026-05-13): literal donor-vs-recipient scan RT diff. iRT is
    # spline-calibrated and absorbs gradient drift between files; raw RT
    # carries independent species-discrimination signal.
    out_rt_t        = fill(-1f0, n); out_rt_f        = fill(-1f0, n)

    # Per-thread present-counter buffers (avoid atomic contention; sum at end).
    # Size by maxthreadid(), not nthreads(): Julia 1.9+ may run @threads tasks
    # on the interactive thread (threadid > nthreads), and @inbounds writes to
    # an under-sized buffer corrupt the heap. Bus errors / abort traps / TypeErrors
    # observed when launched without `julia -t N` (nthreads()=1, maxthreadid()=2).
    nt = Threads.maxthreadid()
    n_true_per_thread  = zeros(Int, nt)
    n_false_per_thread = zeros(Int, nt)

    # Threaded per-PSM walk. Each row writes to disjoint output indices and
    # only reads (without modifying) the shared `all_entries`, `med_*`, and
    # input arrays — no synchronization needed.
    Threads.@threads :static for i in 1:n
        tid = Threads.threadid()
        @inbounds begin
            my_file = file_v[i]
            my_pid  = pid_col[i]
            my_partner = partner_col[i]
            # _true donor
            donor_t = _donor_for(my_pid, my_file)
            if donor_t !== nothing
                n_true_per_thread[tid] += 1
                out_max_t[i] = donor_t.trace_prob
                if donor_t.weight > 0 && weight_v[i] > 0
                    out_lw_t[i] = log2(weight_v[i] / donor_t.weight)
                end
                out_le_t[i] = log2ie_v[i] - donor_t.log2_intensity_explained
                if has_irt
                    out_ir_t[i] = abs(irt_res_v[i] - donor_t.irt_residual)
                end
                if has_logby
                    out_log_by_t[i] = log_by_v[i] - donor_t.log_by_ratio
                end
                if has_rt
                    out_rt_t[i] = abs(rt_obs_v[i] - donor_t.rt_obs)
                end
                out_miss_t[i] = false
            end
            # _false donor
            donor_f = _donor_for(my_partner, my_file)
            if donor_f !== nothing
                n_false_per_thread[tid] += 1
                out_max_f[i] = donor_f.trace_prob
                if donor_f.weight > 0 && weight_v[i] > 0
                    out_lw_f[i] = log2(weight_v[i] / donor_f.weight)
                end
                out_le_f[i] = log2ie_v[i] - donor_f.log2_intensity_explained
                if has_irt
                    out_ir_f[i] = abs(irt_res_v[i] - donor_f.irt_residual)
                end
                if has_logby
                    out_log_by_f[i] = log_by_v[i] - donor_f.log_by_ratio
                end
                if has_rt
                    out_rt_f[i] = abs(rt_obs_v[i] - donor_f.rt_obs)
                end
                out_miss_f[i] = false
            end
        end
    end
    n_true_present  = sum(n_true_per_thread)
    n_false_present = sum(n_false_per_thread)

    psms[!, :MBR_max_pair_prob_true]         = out_max_t
    psms[!, :MBR_max_pair_prob_false]        = out_max_f
    psms[!, :MBR_log2_weight_ratio_true]     = out_lw_t
    psms[!, :MBR_log2_weight_ratio_false]    = out_lw_f
    psms[!, :MBR_log2_explained_ratio_true]  = out_le_t
    psms[!, :MBR_log2_explained_ratio_false] = out_le_f
    psms[!, :MBR_best_irt_diff_true]         = out_ir_t
    psms[!, :MBR_best_irt_diff_false]        = out_ir_f
    psms[!, :MBR_is_missing_true]            = out_miss_t
    psms[!, :MBR_is_missing_false]           = out_miss_f
    psms[!, :MBR_log_by_diff_true]           = out_log_by_t
    psms[!, :MBR_log_by_diff_false]          = out_log_by_f
    # rtv1 (2026-05-13)
    psms[!, :MBR_best_rt_diff_true]          = out_rt_t
    psms[!, :MBR_best_rt_diff_false]         = out_rt_f

    elapsed = time() - t0
    @user_info "MBR Batch F — dual donor features:"
    @user_info "  unique precursor_idx tracked: $(length(all_entries))"
    @user_info "  PSMs with TRUE  donor: $n_true_present / $n  ($(round(100*n_true_present/max(n,1), digits=1))%)"
    @user_info "  PSMs with FALSE donor: $n_false_present / $n  ($(round(100*n_false_present/max(n,1), digits=1))%)"
    if n_true_present > 0
        med_t = median(out_max_t[.!out_miss_t])
        @user_info "  MBR_max_pair_prob_true  median: $(round(med_t, digits=4))"
    end
    if n_false_present > 0
        med_f = median(out_max_f[.!out_miss_f])
        @user_info "  MBR_max_pair_prob_false median: $(round(med_f, digits=4))"
    end
    @user_info "  compute_mbr_features_dual! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_true_present  = n_true_present,
        n_false_present = n_false_present,
        elapsed_s       = elapsed,
    )
end
