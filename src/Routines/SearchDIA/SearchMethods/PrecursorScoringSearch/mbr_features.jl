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
    irt_obs::Float32       # raw observed iRT of the donor row (for top-N median)
    log_by_ratio::Float32  # log(b_int+1) - log(y_int+1) of the donor row
    rt_obs::Float32        # literal scan RT (minutes) of the donor row — rtv1.
                           # iRT is spline-calibrated; raw RT carries independent
                           # gradient-bounded signal for species-mismatch transfers.
    # Per-rank M0 fragment intensities of the donor row (for cross-run fragment
    # pattern similarity — the develop-schema analog of v0.6.6's chromatogram
    # RV-coefficient).
    frag1_int::Float32
    frag2_int::Float32
    frag3_int::Float32
    frag4_int::Float32
    frag5_int::Float32
    frag6_int::Float32
    ms_file_idx::UInt32
    is_decoy::Bool
end
const _MBR_DONOR_SENTINEL = _MBRDonorEntry(typemin(Float32), 0f0, 0f0, 0f0, 0f0, 0f0,
                                            0f0,
                                            0f0, 0f0, 0f0, 0f0, 0f0, 0f0,
                                            UInt32(0), false)

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
                           0f0, 0f0, 0f0, 0f0, 0f0, 0f0,
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
    # Per-rank M0 fragment intensities (for fragment-pattern cosine/scribe).
    has_frag = all(c -> hasproperty(psms, c),
                   (:frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int))
    frag_v = has_frag ?
        (Float32.(psms.frag1_int), Float32.(psms.frag2_int), Float32.(psms.frag3_int),
         Float32.(psms.frag4_int), Float32.(psms.frag5_int), Float32.(psms.frag6_int)) :
        (zeros(Float32, nrow(psms)), zeros(Float32, nrow(psms)), zeros(Float32, nrow(psms)),
         zeros(Float32, nrow(psms)), zeros(Float32, nrow(psms)), zeros(Float32, nrow(psms)))
    n = nrow(psms)

    # Top-K per precursor_idx for the consensus aggregates (median score,
    # median irt_obs). K = ceil(sqrt(N_files)).
    # N=1 → K=1; N=2 → K=2; N=5 → K=3; N=10 → K=4; N=20 → K=5; N=100 → K=10.
    # PIONEER_MBR_DONOR_MODE = "topk" (default) | "top1" gates the experiment:
    #   "topk": full K = ceil(sqrt(N_files)) donors, median aggregates over top-K
    #   "top1": K = 2 (only enough for same-file fallback in _donor_for); the
    #           "median" aggregates degenerate to the single best donor's value,
    #           so MBR_top_n_median_score_*/MBR_top_n_irt_diff_* become
    #           duplicates of MBR_max_pair_prob_* / a single-best irt_obs diff.
    #           Tests whether top-K aggregation is worth the complexity.
    donor_mode = lowercase(get(ENV, "PIONEER_MBR_DONOR_MODE", "topk"))
    n_files = Int(maximum(file_v))
    K = donor_mode == "top1" ? 2 :
        max(1, ceil(Int, sqrt(Float32(n_files))))
    # Collect all entries per pid, then trim to top-K by trace_prob descending.
    all_entries = Dict{UInt32, Vector{_MBRDonorEntry}}()
    sizehint!(all_entries, max(64, n ÷ 8))
    @inbounds for i in 1:n
        pid = pid_col[i]
        e = _MBRDonorEntry(score_v[i], weight_v[i], log2ie_v[i], irt_res_v[i],
                           irt_obs_v[i], log_by_v[i], rt_obs_v[i],
                           frag_v[1][i], frag_v[2][i], frag_v[3][i],
                           frag_v[4][i], frag_v[5][i], frag_v[6][i],
                           file_v[i], false)
        push!(get!(() -> _MBRDonorEntry[], all_entries, pid), e)
    end
    # Sort each pid's entries by trace_prob desc, keep top-K.
    @inbounds for (_pid, entries) in all_entries
        sort!(entries; by = e -> -e.trace_prob)
        if length(entries) > K
            resize!(entries, K)
        end
    end

    # Per-pid aggregates. In topk mode = median over top-K entries; in top1
    # mode = the single best donor's score / irt_obs (no aggregation).
    med_score_per_pid    = Dict{UInt32, Float32}()
    med_irt_obs_per_pid  = Dict{UInt32, Float32}()
    sizehint!(med_score_per_pid, length(all_entries))
    sizehint!(med_irt_obs_per_pid, length(all_entries))
    if donor_mode == "top1"
        @inbounds for (pid, entries) in all_entries
            isempty(entries) && continue
            med_score_per_pid[pid]   = entries[1].trace_prob
            med_irt_obs_per_pid[pid] = entries[1].irt_obs
        end
    else
        @inbounds for (pid, entries) in all_entries
            isempty(entries) && continue
            sc = Float32[e.trace_prob for e in entries]
            ir = Float32[e.irt_obs    for e in entries]
            med_score_per_pid[pid]   = Float32(median(sc))
            med_irt_obs_per_pid[pid] = Float32(median(ir))
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
    # New features (2026-05-12 evening additions):
    out_med_score_t = fill(-1f0, n); out_med_score_f = fill(-1f0, n)
    out_top_irt_t   = fill(-1f0, n); out_top_irt_f   = fill(-1f0, n)
    out_log_by_t    = fill(-1f0, n); out_log_by_f    = fill(-1f0, n)
    # Fragment-pattern cosine + scribe (donor vs recipient observed frag-int 6-vector).
    # Re-incarnation of v0.6.6's MBR_rv_coefficient using fragment-intensity
    # fingerprints (chromatogram vectors are no longer available post-MainSearch).
    out_frag_cos_t  = fill(-1f0, n); out_frag_cos_f  = fill(-1f0, n)
    out_frag_scr_t  = fill(-1f0, n); out_frag_scr_f  = fill(-1f0, n)
    # rtv1 (2026-05-13): literal donor-vs-recipient scan RT diff. iRT is
    # spline-calibrated and absorbs gradient drift between files; raw RT
    # carries independent species-discrimination signal.
    out_rt_t        = fill(-1f0, n); out_rt_f        = fill(-1f0, n)
    # rtv2 (2026-05-13): Kendall τ on M0 fragment intensity 6-vector. Rank-based
    # so magnitude-invariant; targets species-mismatch noise patterns directly.
    out_frag_kt_t   = fill(-2f0, n); out_frag_kt_f   = fill(-2f0, n)

    # Per-row scratch: normalized recipient frag vector. Compute on the fly.
    @inline function _norm_l1(a, b, c, d, e, f)
        s = a + b + c + d + e + f
        s > 0f0 ? (a/s, b/s, c/s, d/s, e/s, f/s) : (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    end
    @inline function _cosine6(p::NTuple{6,Float32}, q::NTuple{6,Float32})
        dotpq = p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4] + p[5]*q[5] + p[6]*q[6]
        np2 = p[1]^2 + p[2]^2 + p[3]^2 + p[4]^2 + p[5]^2 + p[6]^2
        nq2 = q[1]^2 + q[2]^2 + q[3]^2 + q[4]^2 + q[5]^2 + q[6]^2
        d = sqrt(np2 * nq2)
        d > 0f0 ? Float32(dotpq / d) : 0f0
    end
    # rtv2 (2026-05-13): Kendall τ over 15 ordered pairs of the 6-vector. Tied
    # pairs (zero diff on either side) excluded from the denominator. Range
    # [-1, +1]; sentinel -2f0 for "no donor or no informative pair." Rank-based
    # so insensitive to magnitude noise — targets the species-mismatch failure
    # mode where matched "fragments" in a HelaOnly file at YEAST m/z targets
    # have random ordering vs the structured donor pattern.
    @inline function _kendall6(p::NTuple{6,Float32}, q::NTuple{6,Float32})
        nc = 0; nd = 0
        @inbounds for i in 1:5
            for j in (i+1):6
                di = p[i] - p[j]; dj = q[i] - q[j]
                si = di > 0f0 ? 1 : (di < 0f0 ? -1 : 0)
                sj = dj > 0f0 ? 1 : (dj < 0f0 ? -1 : 0)
                if si != 0 && sj != 0
                    si == sj ? (nc += 1) : (nd += 1)
                end
            end
        end
        nt = nc + nd
        nt > 0 ? Float32((nc - nd) / nt) : -2f0
    end
    @inline function _scribe6(p::NTuple{6,Float32}, q::NTuple{6,Float32})
        absd = abs(p[1]-q[1]) + abs(p[2]-q[2]) + abs(p[3]-q[3]) +
               abs(p[4]-q[4]) + abs(p[5]-q[5]) + abs(p[6]-q[6])
        Float32(-log2(absd / 2f0 + 1f-6))
    end

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
            # Normalized recipient fragment-intensity 6-vector (lazy; zero-pattern
            # when no frag information available).
            v_self_norm = _norm_l1(frag_v[1][i], frag_v[2][i], frag_v[3][i],
                                    frag_v[4][i], frag_v[5][i], frag_v[6][i])
            v_self_has_signal = (v_self_norm[1] + v_self_norm[2] + v_self_norm[3] +
                                  v_self_norm[4] + v_self_norm[5] + v_self_norm[6]) > 0f0
            # Per-pid aggregates (file-independent; safe to read concurrently).
            med_s_true = get(med_score_per_pid,   my_pid,     -1f0)
            med_i_true = get(med_irt_obs_per_pid, my_pid,     irt_obs_v[i])
            med_s_fls  = get(med_score_per_pid,   my_partner, -1f0)
            med_i_fls  = get(med_irt_obs_per_pid, my_partner, irt_obs_v[i])
            out_med_score_t[i] = med_s_true
            out_med_score_f[i] = med_s_fls
            if has_irt
                out_top_irt_t[i] = abs(irt_obs_v[i] - med_i_true)
                out_top_irt_f[i] = abs(irt_obs_v[i] - med_i_fls)
            end
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
                if has_frag && v_self_has_signal
                    v_donor_norm = _norm_l1(donor_t.frag1_int, donor_t.frag2_int, donor_t.frag3_int,
                                             donor_t.frag4_int, donor_t.frag5_int, donor_t.frag6_int)
                    v_donor_has_signal = (v_donor_norm[1] + v_donor_norm[2] + v_donor_norm[3] +
                                           v_donor_norm[4] + v_donor_norm[5] + v_donor_norm[6]) > 0f0
                    if v_donor_has_signal
                        out_frag_cos_t[i] = _cosine6(v_self_norm, v_donor_norm)
                        out_frag_scr_t[i] = _scribe6(v_self_norm, v_donor_norm)
                        out_frag_kt_t[i]  = _kendall6(v_self_norm, v_donor_norm)
                    end
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
                if has_frag && v_self_has_signal
                    v_donor_norm = _norm_l1(donor_f.frag1_int, donor_f.frag2_int, donor_f.frag3_int,
                                             donor_f.frag4_int, donor_f.frag5_int, donor_f.frag6_int)
                    v_donor_has_signal = (v_donor_norm[1] + v_donor_norm[2] + v_donor_norm[3] +
                                           v_donor_norm[4] + v_donor_norm[5] + v_donor_norm[6]) > 0f0
                    if v_donor_has_signal
                        out_frag_cos_f[i] = _cosine6(v_self_norm, v_donor_norm)
                        out_frag_scr_f[i] = _scribe6(v_self_norm, v_donor_norm)
                        out_frag_kt_f[i]  = _kendall6(v_self_norm, v_donor_norm)
                    end
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
    # New features (2026-05-12 evening additions)
    psms[!, :MBR_top_n_median_score_true]    = out_med_score_t
    psms[!, :MBR_top_n_median_score_false]   = out_med_score_f
    psms[!, :MBR_top_n_irt_diff_true]        = out_top_irt_t
    psms[!, :MBR_top_n_irt_diff_false]       = out_top_irt_f
    psms[!, :MBR_log_by_diff_true]           = out_log_by_t
    psms[!, :MBR_log_by_diff_false]          = out_log_by_f
    psms[!, :MBR_frag_pattern_cosine_true]   = out_frag_cos_t
    psms[!, :MBR_frag_pattern_cosine_false]  = out_frag_cos_f
    psms[!, :MBR_frag_pattern_scribe_true]   = out_frag_scr_t
    psms[!, :MBR_frag_pattern_scribe_false]  = out_frag_scr_f
    # rtv1 + rtv2 (2026-05-13)
    psms[!, :MBR_best_rt_diff_true]          = out_rt_t
    psms[!, :MBR_best_rt_diff_false]         = out_rt_f
    psms[!, :MBR_frag_rank_corr_true]        = out_frag_kt_t
    psms[!, :MBR_frag_rank_corr_false]       = out_frag_kt_f

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
