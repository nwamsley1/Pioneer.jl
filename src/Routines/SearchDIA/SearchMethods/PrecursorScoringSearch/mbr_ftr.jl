# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR Phase 4/5b/6 — False-Transfer-Rate controller as a RECOVERY mechanism
# scoped to the transfer-candidate cohort, with a dedicated FTR LightGBM model.
#
# Phase 6 change vs Phase 5b:
#   * The MBR-boosted score :trace_prob_mbr alone cannot separate good from
#     bad transfers tightly enough at α=0.01 (recovered ≈ 0).
#   * We train a second LightGBM on the candidate cohort with label
#     :is_bad_transfer (i.e. positives = good transfer). Features are the
#     MBR evidence + boosted/prepass scores; we deliberately exclude
#     :target, :decoy, :MBR_is_best_decoy to avoid label leakage.
#   * 2-fold CV on the existing :cv_fold so each candidate gets an
#     out-of-fold FTR score.
#   * `get_ftr_threshold` then runs on the OOF FTR score over candidates.
#
# Architecture unchanged otherwise:
#   * The qval pipeline downstream uses :trace_prob (= :trace_prob_prepass,
#     non-MBR pass-1). Non-candidate rows can never be score-inflated by MBR.
#   * Recovered candidates get :mbr_recovered = true; PrecursorScoringSearch's
#     qval bypass step exempts them from the per-file qval cut.

# Donor q-value threshold for candidate eligibility. Defines the prepass-score
# floor a donor must clear (= score-at-q≤MBR_DONOR_Q_THRESHOLD among targets).
# Tighter (0.001) → only very confident donors → fewer candidates, safer FTR.
# Looser (1.0) → any cross-run PSM is a valid donor → more candidates, broader recovery.
const MBR_DONOR_Q_THRESHOLD = Float32(0.01)

# ─────────────────────────────────────────────────────────────────────────────
# Legacy Phase 6 (single-frame) FTR controller removed 2026-05-18.
# The Batch F paired controller (apply_mbr_filter_paired! below) supplanted it;
# the legacy FTR_FEATURES list + _train_ftr_lgbm + apply_mbr_filter! were dead
# code (no callers). Helpers (_transfer_class_counts, _pct, _log_transfer_counts)
# were only used by the legacy path and were removed with it.
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# MBR Batch F — paired-counterfactual FTR controller with q-value threshold
# ─────────────────────────────────────────────────────────────────────────────

# FTR features for the Batch F controller. Same as Phase 8f's FTR_FEATURES but
# - drops :trace_prob_mbr (Pass 2 no longer exists in Batch F)
# - uses the _true variants for the in-place MBR features (the _false copies
#   are swapped in for the bottom half of the doubled training frame).
# rtv3 (2026-05-13) + rtv4 (2026-05-14) — slim FTR feature set.
# The pre-rtv3 list had ~50 features. The ~40 non-MBR features were
# row-level (identical between top and bottom halves of the doubled
# training frame), so they couldn't directly discriminate real-vs-
# counterfactual. Empirically the FTR LGBM was using them for
# interaction modeling. Dropping them gained +1,015 IDs (rtv3 doc)
# and dropped Dennis FTR by 0.31 pp.
#
# rtv4 dropped `trace_prob_prepass` (OOF) too — the in-fold score
# absorbs ~88% of the combined OOF+infold LGBM gain. The trace_prob_*
# columns are still computed and used downstream (qval pipeline,
# MBR top-K aggregation) — they just aren't exposed to the FTR LGBM
# as features.
#
# Dead features dropped (gain ≈ 0 across all rtv3/rtv4 runs):
#   - MBR_is_missing_true (the -1 sentinel on every other MBR feature
#     already encodes missingness)
#   - MBR_frag_top1_match_true (was redundant with MBR_frag_rank_corr_true;
#     both subsequently dropped along with the other M0 frag-6-vector features)
const FTR_FEATURES_F_TRUE = Symbol[
    :trace_prob_infold,                            # rtv3: in-fold Pass-1 score
    :MBR_max_pair_prob_true,
    :MBR_log2_weight_ratio_true,
    :MBR_log2_explained_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_log_by_diff_true,
    # rtv1 (2026-05-13): literal donor-vs-recipient scan RT diff.
    # Standalone A/B was negative — kept in slim FTR because it gets
    # ~5,900 gain at rank #6 in rtv4 alongside trace_prob_infold.
    # Don't drop without first restoring the in-fold score context.
    :MBR_best_rt_diff_true,
    :MBR_smoothed_frag_hellinger_true,
]
# Dropped 2026-05-16:
#   :MBR_top_n_median_score_true / :MBR_top_n_irt_diff_true
#   A/B on Olsen 23-file + MTAC 6-file (top-K vs top-1 donor): no material ID,
#   PG, or EFDR-MAE change. K trimmed to 2 for same-file fallback only.
# Dropped 2026-05-17:
#   :MBR_frag_pattern_cosine_true / :_scribe_true / :MBR_frag_rank_corr_true
#   A/B on YeastMBR 40-file + HelaOnly 20-file: dropping the three M0 fragment
#   6-vector features lowered Dennis FTR 2.59 % → 2.40 % and raised PGs +0.16 %,
#   at a cost of −0.06 % precursor IDs / −2.1 % MBR recoveries (YeastMBR);
#   HelaOnly pure-FTR 0.252 % → 0.242 %, PGs +0.33 %, IDs −0.65 %. Net quality
#   win; donor entries no longer need frag1_int..frag8_int.

# Same features but with the MBR columns swapped to _false. Used for the
# bottom half of the doubled training frame. Each entry must mirror
# FTR_FEATURES_F_TRUE position-for-position so the LGBM sees the same
# semantic feature columns in the same order. Row-level features
# (trace_prob_infold) pass through the default `f` branch unchanged.
const FTR_FEATURES_F_FALSE = Symbol[
    f === :MBR_max_pair_prob_true        ? :MBR_max_pair_prob_false :
    f === :MBR_log2_weight_ratio_true    ? :MBR_log2_weight_ratio_false :
    f === :MBR_log2_explained_ratio_true ? :MBR_log2_explained_ratio_false :
    f === :MBR_best_irt_diff_true        ? :MBR_best_irt_diff_false :
    f === :MBR_log_by_diff_true          ? :MBR_log_by_diff_false :
    f === :MBR_best_rt_diff_true         ? :MBR_best_rt_diff_false :
    f === :MBR_smoothed_frag_hellinger_true ? :MBR_smoothed_frag_hellinger_false :
    f
    for f in FTR_FEATURES_F_TRUE
]

"""
    apply_mbr_filter_paired!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

MBR Batch F — paired-counterfactual FTR controller.

Algorithm (see BATCH_F_PLAN.md):
1. Candidate gate = (pre_qvals > q_thresh) & (MBR_max_pair_prob_true ≥
   prob_thresh) & (!MBR_is_missing_true) & (!MBR_is_missing_false).
2. Build a row-doubled training frame:
   - top half: candidate rows with FTR_FEATURES_F_TRUE values  → label is_real=true
   - bottom half: SAME candidate rows with FTR_FEATURES_F_FALSE values → label is_real=false
3. Train 2-fold CV LightGBM on the doubled frame with label=is_real.
   OOF scores: higher = more "real" by the model's lights.
4. Compute q-values on the doubled frame via get_qvalues! (target=is_real,
   monotonized). Pick rows whose q ≤ alpha.
5. Recovery: candidate rows whose TOP-HALF (real-MBR) row has q ≤ alpha.

Sets `:mbr_recovered`, `:MBR_transfer_candidate`, `:ftr_qval_true` on
`psms` in-place. Does NOT mutate `:trace_prob`.
"""
function apply_mbr_filter_paired!(
    psms::DataFrame;
    alpha::Float32    = 0.01f0,
    q_thresh::Float32 = 0.01f0,
)
    t0 = time()

    n = nrow(psms)
    if n == 0
        @debug_l1 "MBR Batch F — apply_mbr_filter_paired!: no PSMs to filter"
        psms[!, :mbr_recovered]          = falses(0)
        psms[!, :MBR_transfer_candidate] = falses(0)
        psms[!, :ftr_qval_true]          = Float32[]
        psms[!, :ftr_pep_true]           = Float32[]
        return (n_candidates=0, threshold=Float32(Inf), n_recovered=0, elapsed_s=0.0)
    end

    # ── 1. Pre-MBR q-values from the non-MBR (pass-1) score ──
    pre_score  = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    decoy_col  = Bool.(psms[!, :decoy])
    pre_qvals  = Vector{Float32}(undef, n)
    get_qvalues!(pre_score, target_col, pre_qvals)

    donor_target_pass = (pre_qvals .<= MBR_DONOR_Q_THRESHOLD) .& target_col
    prob_thresh = if any(donor_target_pass)
        minimum(pre_score[donor_target_pass])
    else
        Float32(Inf)
    end

    # ── 2. Candidates: failed q-cut, have BOTH _true and _false donors ──
    mbr_pp_t       = Float32.(psms[!, :MBR_max_pair_prob_true])
    mbr_miss_t     = Bool.(psms[!, :MBR_is_missing_true])
    mbr_miss_f     = Bool.(psms[!, :MBR_is_missing_false])
    candidate_mask = (pre_qvals .> q_thresh) .&
                     (mbr_pp_t  .>= prob_thresh) .&
                     (.!mbr_miss_t) .&
                     (.!mbr_miss_f)
    n_cand = count(candidate_mask)
    cand_idx = findall(candidate_mask)

    @debug_l1 "MBR Batch F — paired FTR recovery:"
    @debug_l1 "  MBR_DONOR_Q_THRESHOLD = $MBR_DONOR_Q_THRESHOLD; prob_thresh = $(round(prob_thresh, digits=4))"
    @debug_l1 "  candidates (pre-MBR q>$q_thresh + both donors): $n_cand / $n ($(round(100*n_cand/n, digits=2))%)"
    @debug_l1 "  α (q-value FTR budget): $alpha"

    if n_cand == 0
        psms[!, :mbr_recovered]          = falses(n)
        psms[!, :MBR_transfer_candidate] = candidate_mask
        psms[!, :ftr_qval_true]          = fill(NaN32, n)
        psms[!, :ftr_pep_true]           = fill(NaN32, n)
        return (n_candidates=0, threshold=Float32(Inf), n_recovered=0, elapsed_s=time()-t0)
    end

    # ── 3. Build the doubled training frame ──
    # Top half: real-MBR rows (FTR_FEATURES_F_TRUE values). Label = true.
    # Bottom half: same candidates with _false MBR values swapped in. Label = false.
    sub = psms[cand_idx, :]
    available_true  = filter(f -> hasproperty(sub, f), FTR_FEATURES_F_TRUE)
    # Position-aligned _false set (only matters that the MBR cols flip; non-MBR
    # cols are identical between TRUE and FALSE lists, so available_true tells
    # us which non-MBR cols actually exist).
    feat_to_idx = Dict{Symbol, Int}()
    for (i, f) in enumerate(FTR_FEATURES_F_TRUE); feat_to_idx[f] = i; end
    available_false = Symbol[]
    for f in available_true
        push!(available_false, FTR_FEATURES_F_FALSE[feat_to_idx[f]])
    end
    X_true  = feature_matrix(sub, available_true)
    X_false = feature_matrix(sub, available_false)
    X       = vcat(X_true, X_false)              # 2 n_cand × n_features
    y       = vcat(trues(n_cand), falses(n_cand))  # label = is_real
    cv_double = vcat(sub[!, :cv_fold], sub[!, :cv_fold])

    # ── 4. 2-fold CV LightGBM on the doubled frame ──
    pos0 = findall(cv_double .== 0)
    pos1 = findall(cv_double .== 1)
    ftr_score_double = fill(NaN32, 2 * n_cand)
    last_cls = nothing
    for (train_pos, test_pos) in ((pos1, pos0), (pos0, pos1))
        (isempty(train_pos) || isempty(test_pos)) && continue
        y_tr = y[train_pos]
        if length(unique(y_tr)) == 1
            ftr_score_double[test_pos] .= y_tr[1] ? 1f0 : 0f0
            continue
        end
        cls = build_lightgbm_classifier(; SHARED_LGBM_HP...)
        LightGBM.fit!(cls, X[train_pos, :], _prepare_labels(y_tr); verbosity = -1)
        raw = LightGBM.predict(cls, X[test_pos, :])
        s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
        ftr_score_double[test_pos] .= Float32.(s)
        last_cls = cls
    end

    # ── 5. Q-value AND PEP on the doubled frame (target=is_real, decoy=is_fake) ──
    qvals_double = Vector{Float32}(undef, 2 * n_cand)
    pep_double   = Vector{Float32}(undef, 2 * n_cand)
    get_qvalues!(ftr_score_double, y, qvals_double)
    get_PEP!(ftr_score_double, y, pep_double)
    # qvals_double[i] = (# fakes ranked ≥ row i) / (# reals ranked ≥ row i),
    # monotonized to be non-increasing as score increases.
    # pep_double[i]   = per-row P(is_false | score), isotonic-regression-derived.

    # ── 6. Recovery: top-half rows (real-MBR) with pep ≤ alpha.
    #    PEP-based recovery at α=0.01 gives ~half the Dennis-FTR species-
    #    mismatch rate vs qval-based at the same α (validated 2026-05-15 on
    #    40-file YeastMBR + 20-file HelaOnly: any-YEAST 9.66% qval → 4.49% pep,
    #    strict-YEAST 8.86% → 3.32%). PEP is more conservative — ~33% fewer
    #    ΔTotal recoveries — but the surviving recoveries are dramatically
    #    cleaner. The qval path is retained in git history if needed.
    qvals_top = qvals_double[1:n_cand]
    pep_top   = pep_double[1:n_cand]
    recovered_in_cand = pep_top .<= alpha
    n_recovered = count(recovered_in_cand)

    mbr_recovered_full = falses(n)
    ftr_qval_full      = fill(NaN32, n)
    ftr_pep_full       = fill(NaN32, n)
    @inbounds for (k, i) in enumerate(cand_idx)
        ftr_qval_full[i] = qvals_top[k]
        ftr_pep_full[i]  = pep_top[k]
        if recovered_in_cand[k]
            mbr_recovered_full[i] = true
        end
    end
    psms[!, :mbr_recovered]          = mbr_recovered_full
    psms[!, :MBR_transfer_candidate] = candidate_mask
    psms[!, :ftr_qval_true]          = ftr_qval_full
    psms[!, :ftr_pep_true]           = ftr_pep_full

    # ── 7. Threshold for log: highest ftr_score on a top-half row with q ≤ α ──
    τ = n_recovered > 0 ?
        minimum(ftr_score_double[1:n_cand][recovered_in_cand]) :
        Float32(Inf)

    # ── 8. Recovered transfer composition ──
    # Use library label + MBR_is_best_decoy-equivalent — donor polarity from
    # the true-donor's target/decoy class. The dual-features file does not
    # currently log donor target/decoy; infer from the precursor itself (T←T
    # for targets, D←D for decoys, since the donor is the SAME precursor).
    # T←D / D←T can only happen via the counterfactual, which is never used
    # for recovery; here counts are simpler: T or D recoveries.
    n_t_rec = count(i -> mbr_recovered_full[i] && target_col[i], 1:n)
    n_d_rec = n_recovered - n_t_rec

    @debug_l1 "  doubled-frame rows: $(2*n_cand)  (top=real, bottom=fake)"
    @debug_l1 "  τ (FTR score, q ≤ α): $(round(τ, digits=4))"
    @debug_l1 "  RECOVERED (top-half q ≤ α): $n_recovered ($(round(100*n_recovered/max(n_cand,1), digits=2))% of candidates)"
    @debug_l1 "  recovered targets: $n_t_rec   recovered decoys: $n_d_rec"

    if last_cls !== nothing
        lgbm_model = LightGBMModel(last_cls, available_true, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            @debug_l1 "  Batch F FTR-model feature importances (gain):"
            for (fname, gain) in sorted_imp
                @debug_l1 "    $(rpad(string(fname), 36)) $(round(gain, digits=2))"
            end
        end
    end

    elapsed = time() - t0
    @debug_l1 "  apply_mbr_filter_paired! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_candidates = n_cand,
        threshold    = τ,
        n_recovered  = n_recovered,
        prob_thresh  = prob_thresh,
        elapsed_s    = elapsed,
    )
end
