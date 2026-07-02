# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR false-transfer controller as a recovery mechanism scoped to the
# transfer-candidate cohort.
#
# The paired frame contains each candidate twice:
#   * top half: candidate features using the same-precursor donor
#   * bottom half: the same receiver with counterfactual donor features
# The FTR label is true for top-half same-precursor transfers and false for
# bottom-half counterfactual transfers. Semi-supervised iterations use only the
# top-half rows passing the 3% counterfactual FTR threshold as positives for the
# next training round and for the next stopping evaluation; unselected top-half
# rows are omitted from the semi-supervised objective.
# The resulting OOF score drives FTR recovery.
#
# Architecture unchanged otherwise:
#   * The qval pipeline downstream uses :trace_prob (= :trace_prob_prepass,
#     non-MBR pass-1). Non-candidate rows can never be score-inflated by MBR.
#   * Recovered candidates get :mbr_recovered = true inside the post-qvalue
#     filtered set; downstream PEP and run-level q-values are recomputed
#     without a recovery bypass.

# Donor q-value threshold for candidate eligibility. Defines the prepass-score
# floor a donor must clear (= score-at-q≤MBR_DONOR_Q_THRESHOLD among targets).
# Tighter (0.001) → only very confident donors → fewer candidates, safer FTR.
# Looser (1.0) → any cross-run PSM is a valid donor → more candidates, broader recovery.
const MBR_DONOR_Q_THRESHOLD = Float32(0.01)
const MBR_SEMISUPERVISED_FTR_THRESHOLD = Float32(0.03)

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
#   - MBR_best_is_missing_true (the -1 sentinel on every other MBR feature
#     already encodes missingness)
#   - MBR_frag_top1_match_true (was redundant with MBR_frag_rank_corr_true;
#     both subsequently dropped along with the other M0 frag-6-vector features)
const FTR_FEATURES_F_TRUE = Symbol[
    :main_search_prob,
    :trace_prob_infold,                            # rtv3: in-fold Pass-1 score
    :MBR_best_pair_prob_true,
    :MBR_worst_pair_prob_true,
    :MBR_log2_weight_lod_ratio,
    :MBR_best_log2_weight_ratio_true,
    :MBR_worst_log2_weight_ratio_true,
    :MBR_best_log2_explained_ratio_true,
    :MBR_worst_log2_explained_ratio_true,
    :MBR_best_abs_n_scans_diff_true,
    :MBR_worst_abs_n_scans_diff_true,
    :MBR_best_log2_n_scans_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_best_observed_irt_diff_true,
    :MBR_worst_observed_irt_diff_true,
    :MBR_single_donor_true,
    :MBR_best_log_by_diff_true,
    :MBR_best_smoothed_frag_hellinger_true,
    :MBR_worst_smoothed_frag_hellinger_true,
    :MBR_best_donor_library_hellinger_true,
    :MBR_worst_donor_library_hellinger_true,
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
    f === :MBR_best_pair_prob_true       ? :MBR_best_pair_prob_false :
    f === :MBR_worst_pair_prob_true      ? :MBR_worst_pair_prob_false :
    f === :MBR_best_log2_weight_ratio_true ? :MBR_best_log2_weight_ratio_false :
    f === :MBR_worst_log2_weight_ratio_true ? :MBR_worst_log2_weight_ratio_false :
    f === :MBR_best_log2_explained_ratio_true ? :MBR_best_log2_explained_ratio_false :
    f === :MBR_worst_log2_explained_ratio_true ? :MBR_worst_log2_explained_ratio_false :
    f === :MBR_best_abs_n_scans_diff_true ? :MBR_best_abs_n_scans_diff_false :
    f === :MBR_worst_abs_n_scans_diff_true ? :MBR_worst_abs_n_scans_diff_false :
    f === :MBR_best_log2_n_scans_ratio_true ? :MBR_best_log2_n_scans_ratio_false :
    f === :MBR_best_irt_diff_true        ? :MBR_best_irt_diff_false :
    f === :MBR_best_observed_irt_diff_true ? :MBR_best_observed_irt_diff_false :
    f === :MBR_worst_observed_irt_diff_true ? :MBR_worst_observed_irt_diff_false :
    f === :MBR_single_donor_true         ? :MBR_single_donor_false :
    f === :MBR_best_log_by_diff_true     ? :MBR_best_log_by_diff_false :
    f === :MBR_best_smoothed_frag_hellinger_true ? :MBR_best_smoothed_frag_hellinger_false :
    f === :MBR_worst_smoothed_frag_hellinger_true ? :MBR_worst_smoothed_frag_hellinger_false :
    f === :MBR_best_donor_library_hellinger_true ? :MBR_best_donor_library_hellinger_false :
    f === :MBR_worst_donor_library_hellinger_true ? :MBR_worst_donor_library_hellinger_false :
    f
    for f in FTR_FEATURES_F_TRUE
]

function _mbr_transfer_training_labels(positive_top::AbstractVector{Bool})
    n = length(positive_top)
    labels = falses(2 * n)
    @inbounds for i in 1:n
        labels[i] = positive_top[i]
    end
    return labels
end

function _mbr_transfer_training_mask(positive_top::AbstractVector{Bool})
    n = length(positive_top)
    mask = trues(2 * n)
    @inbounds for i in 1:n
        mask[i] = positive_top[i]
    end
    return mask
end

function _mbr_transfer_positive_top(
    ftr_qvals_top::AbstractVector{<:Real};
    ftr_threshold::Float32 = MBR_SEMISUPERVISED_FTR_THRESHOLD,
)
    n = length(ftr_qvals_top)
    positive_top = falses(n)
    @inbounds for i in 1:n
        positive_top[i] = Float32(ftr_qvals_top[i]) <= ftr_threshold
    end
    return positive_top
end

function _mbr_transfer_counterfactual_labels(n_cand::Int)
    labels = falses(2 * n_cand)
    @inbounds for i in 1:n_cand
        labels[i] = true
    end
    return labels
end

function _mbr_transfer_iteration_metrics(
    scores_double::AbstractVector{<:Real},
    n_cand::Int;
    ftr_threshold::Float32 = MBR_SEMISUPERVISED_FTR_THRESHOLD,
    eval_labels::AbstractVector{Bool} = _mbr_transfer_counterfactual_labels(n_cand),
    eval_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
)
    n_double = 2 * n_cand
    @assert length(scores_double) == n_double
    @assert length(eval_labels) == n_double

    qvals_double = fill(Inf32, n_double)
    pep_double = fill(Inf32, n_double)
    if eval_mask === nothing
        get_qvalues!(scores_double, eval_labels, qvals_double)
        get_PEP!(scores_double, eval_labels, pep_double)
    else
        @assert length(eval_mask) == n_double
        eval_idx = findall(eval_mask)
        if !isempty(eval_idx)
            scores_eval = Float32.(scores_double[eval_idx])
            labels_eval = eval_labels[eval_idx]
            qvals_eval = Vector{Float32}(undef, length(eval_idx))
            pep_eval = Vector{Float32}(undef, length(eval_idx))
            get_qvalues!(scores_eval, labels_eval, qvals_eval)
            get_PEP!(scores_eval, labels_eval, pep_eval)
            @inbounds for (j, i) in enumerate(eval_idx)
                qvals_double[i] = qvals_eval[j]
                pep_double[i] = pep_eval[j]
            end
        end
    end

    positive_top = _mbr_transfer_positive_top(
        @view(qvals_double[1:n_cand]),
        ftr_threshold = ftr_threshold,
    )

    return (
        qvals_double = qvals_double,
        pep_double = pep_double,
        qvals_top = qvals_double[1:n_cand],
        pep_top = pep_double[1:n_cand],
        positive_top = positive_top,
        n_positive = count(positive_top),
    )
end

"""
    apply_mbr_filter_paired!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

MBR Batch F — paired-counterfactual FTR controller.

Algorithm (see BATCH_F_PLAN.md):
1. Candidate gate = (pre_qvals > q_thresh) & (MBR_best_pair_prob_true ≥
   prob_thresh) & (!MBR_best_is_missing_true) & (!MBR_best_is_missing_false).
2. Build a row-doubled training frame:
   - top half: candidate rows with FTR_FEATURES_F_TRUE values.
   - bottom half: same receivers with FTR_FEATURES_F_FALSE values.
3. Train iterative 2-fold CV LightGBM. Each iteration uses the previous round's
   3% counterfactual-FTR top-half rows as positives, plus counterfactual rows
   as negatives.
4. Compute q-values and PEPs on the doubled frame. Recover candidate rows whose
   top-half row has q-value ≤ alpha.

Sets `:mbr_recovered`, `:MBR_transfer_candidate`, `:ftr_qval_true` on
`psms` in-place. Does NOT mutate `:trace_prob`.
"""
function _mbr_recovery_mask(qvals_top, pep_top, alpha::Float32)
    return Float32.(qvals_top) .<= alpha
end

function apply_mbr_filter_paired!(
    psms::DataFrame;
    alpha::Float32    = 0.01f0,
    q_thresh::Float32 = 0.01f0,
    prob_thresh_override::Union{Nothing, Float32} = nothing,
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
    pre_qvals = if hasproperty(psms, :qval)
        Float32.(psms[!, :qval])
    else
        qvals = Vector{Float32}(undef, n)
        get_qvalues!(pre_score, target_col, qvals)
        qvals
    end
    global_pass = if hasproperty(psms, :global_qval)
        Float32.(psms[!, :global_qval]) .<= q_thresh
    else
        trues(n)
    end

    donor_target_pass = (pre_qvals .<= MBR_DONOR_Q_THRESHOLD) .& target_col
    prob_thresh = if prob_thresh_override !== nothing
        prob_thresh_override
    elseif any(donor_target_pass)
        minimum(pre_score[donor_target_pass])
    else
        Float32(Inf)
    end

    # ── 2. Candidates: failed q-cut, have BOTH _true and _false donors ──
    mbr_pp_t       = Float32.(psms[!, :MBR_best_pair_prob_true])
    mbr_miss_t     = Bool.(psms[!, :MBR_best_is_missing_true])
    mbr_miss_f     = Bool.(psms[!, :MBR_best_is_missing_false])
    candidate_mask = global_pass .&
                     (pre_qvals .> q_thresh) .&
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
    # Top half: true-donor rows, positive under the counterfactual FTR label.
    # Bottom half: same receivers with counterfactual MBR values, always negative.
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
    eval_labels = _mbr_transfer_counterfactual_labels(n_cand)
    eval_mask = trues(2 * n_cand)
    cv_double = vcat(sub[!, :cv_fold], sub[!, :cv_fold])

    # ── 4. Semi-supervised 2-fold CV LightGBM on the doubled frame ──
    pos0 = findall(cv_double .== 0)
    pos1 = findall(cv_double .== 1)
    training_labels = eval_labels
    training_mask = trues(2 * n_cand)
    best_state = nothing
    previous_positive_count = -1

    for iter_idx in 1:SCORING_SEMISUPERVISED_MAX_ITERATIONS
        ftr_score_double = fill(NaN32, 2 * n_cand)
        last_cls_iter = nothing
        n_train_targets = 0
        n_train_decoys = 0

        for (train_pos_all, test_pos) in ((pos1, pos0), (pos0, pos1))
            train_pos = train_pos_all[training_mask[train_pos_all]]
            (isempty(train_pos) || isempty(test_pos)) && continue
            y_tr = training_labels[train_pos]
            n_targets_fold = count(y_tr)
            n_train_targets += n_targets_fold
            n_train_decoys += length(y_tr) - n_targets_fold
            if length(unique(y_tr)) == 1
                ftr_score_double[test_pos] .= y_tr[1] ? 1f0 : 0f0
                continue
            end
            cls = build_lightgbm_classifier(; SHARED_LGBM_HP...)
            LightGBM.fit!(cls, X[train_pos, :], _prepare_labels(y_tr); verbosity = -1)
            raw = LightGBM.predict(cls, X[test_pos, :])
            s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
            ftr_score_double[test_pos] .= Float32.(s)
            last_cls_iter = cls
        end

        metrics = _mbr_transfer_iteration_metrics(
            ftr_score_double,
            n_cand,
            eval_labels = eval_labels,
            eval_mask = eval_mask,
        )
        state = (
            iter = iter_idx,
            scores = ftr_score_double,
            metrics = metrics,
            last_cls = last_cls_iter,
            n_positive = metrics.n_positive,
            n_train_targets = n_train_targets,
            n_train_decoys = n_train_decoys,
        )
        if best_state === nothing || state.n_positive >= best_state.n_positive
            best_state = state
        end

        @debug_l1 "  MBR transfer semi-supervised iter $iter_idx: " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "FTR≤$(round(100 * MBR_SEMISUPERVISED_FTR_THRESHOLD, digits=2))% " *
                   "targets=$(state.n_positive)"

        if iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_positive_count,
            state.n_positive,
        )
            @debug_l1 "  MBR transfer semi-supervised stopping: " *
                       "iter $iter_idx positives=$(state.n_positive) did not improve by " *
                       "$(round(100 * SCORING_SEMISUPERVISED_MIN_TARGET_GAIN, digits=2))% " *
                       "over $previous_positive_count; using iter $(best_state.iter) " *
                       "with positives=$(best_state.n_positive)"
            break
        elseif iter_idx == SCORING_SEMISUPERVISED_MAX_ITERATIONS
            @debug_l1 "  MBR transfer semi-supervised stopping: hit max iterations " *
                       "$SCORING_SEMISUPERVISED_MAX_ITERATIONS; using iter $(best_state.iter) " *
                       "with positives=$(best_state.n_positive)"
            break
        end

        previous_positive_count = state.n_positive
        valid_transfer_labels = _mbr_transfer_training_labels(metrics.positive_top)
        valid_transfer_mask = _mbr_transfer_training_mask(metrics.positive_top)
        training_labels = valid_transfer_labels
        training_mask = valid_transfer_mask
        eval_labels = valid_transfer_labels
        eval_mask = valid_transfer_mask
    end

    # ── 5. Q-value AND PEP on the doubled frame ──
    ftr_score_double = best_state.scores
    qvals_double = best_state.metrics.qvals_double
    pep_double = best_state.metrics.pep_double
    last_cls = best_state.last_cls
    # qvals_double[i] = (# negatives ranked ≥ row i) / (# positives ranked ≥ row i),
    # monotonized to be non-increasing as score increases.
    # pep_double[i]   = per-row P(negative | score), isotonic-regression-derived.

    # ── 6. Recovery: top-half rows (real-MBR) with q-value ≤ alpha.
    # PEP is still retained as a diagnostic/reporting column, but the recovery
    # decision uses the paired FTR q-value.
    qvals_top = qvals_double[1:n_cand]
    pep_top   = pep_double[1:n_cand]
    recovered_in_cand = _mbr_recovery_mask(qvals_top, pep_top, alpha)
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

    @debug_l1 "  doubled-frame rows: $(2*n_cand)  (top=true donor, bottom=counterfactual)"
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
