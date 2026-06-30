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
# The evaluation label is true only for target receivers in the top half. Decoy
# receivers and all counterfactual rows are negatives. Semi-supervised
# iterations use only the 3% FTR ∩ 1% top-half FDR targets as positives for the
# next training round; unselected target receivers are omitted from training.
# The resulting OOF score drives both FTR recovery and downstream recovered-row
# ranking.
#
# Architecture unchanged otherwise:
#   * The qval pipeline downstream uses :trace_prob (= :trace_prob_prepass,
#     non-MBR pass-1). Non-candidate rows can never be score-inflated by MBR.
#   * Recovered candidates get :mbr_recovered = true, bypass only the initial
#     row-level q-value gate, and are filtered by the recalculated row-level
#     q-value after the global q-value pass.

# Donor q-value threshold for candidate eligibility. Defines the prepass-score
# floor a donor must clear (= score-at-q≤MBR_DONOR_Q_THRESHOLD among targets).
# Tighter (0.001) → only very confident donors → fewer candidates, safer FTR.
# Looser (1.0) → any cross-run PSM is a valid donor → more candidates, broader recovery.
const MBR_DONOR_Q_THRESHOLD = Float32(0.01)
const MBR_SEMISUPERVISED_FTR_THRESHOLD = Float32(0.03)
const MBR_SEMISUPERVISED_FDR_THRESHOLD = Float32(0.01)

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
# rtv4 dropped `trace_prob_prepass` (OOF) too — the in-fold score absorbs
# most of the combined OOF+infold LGBM gain. FTR still keeps the per-run
# main-search confidence as a separate feature so rescue rows can provide
# local evidence even though they have no cross-run score.
#
# Dead features dropped (gain ≈ 0 across all rtv3/rtv4 runs):
#   - MBR_is_missing_true (the -1 sentinel on every other MBR feature
#     already encodes missingness)
#   - MBR_frag_top1_match_true (was redundant with MBR_frag_rank_corr_true;
#     both subsequently dropped along with the other M0 frag-6-vector features)
const FTR_FEATURES_F_TRUE = Symbol[
    :main_search_prob,                             # per-run 1 - main_search PEP
    :trace_prob_infold,                            # in-fold cross-run score
    :MBR_max_pair_prob_true,
    :MBR_log2_weight_lod_ratio,
    :MBR_log2_weight_ratio_true,
    :MBR_log2_weight_ratio_worst_true,
    :MBR_log2_explained_ratio_true,
    :MBR_log2_explained_ratio_worst_true,
    :MBR_abs_n_scans_diff_true,
    :MBR_abs_n_scans_diff_worst_true,
    :MBR_log2_n_scans_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_best_irt_diff_worst_true,
    :MBR_log_by_diff_true,
    :MBR_smoothed_frag_hellinger_true,
    :MBR_smoothed_frag_hellinger_worst_true,
    :MBR_donor_library_hellinger_true,
    :MBR_donor_library_hellinger_worst_true,
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
# (main_search_prob, trace_prob_infold) pass through unchanged.
const FTR_FEATURES_F_FALSE = Symbol[
    f === :MBR_max_pair_prob_true        ? :MBR_max_pair_prob_false :
    f === :MBR_log2_weight_ratio_true    ? :MBR_log2_weight_ratio_false :
    f === :MBR_log2_weight_ratio_worst_true ? :MBR_log2_weight_ratio_worst_false :
    f === :MBR_log2_explained_ratio_true ? :MBR_log2_explained_ratio_false :
    f === :MBR_log2_explained_ratio_worst_true ? :MBR_log2_explained_ratio_worst_false :
    f === :MBR_abs_n_scans_diff_true     ? :MBR_abs_n_scans_diff_false :
    f === :MBR_abs_n_scans_diff_worst_true ? :MBR_abs_n_scans_diff_worst_false :
    f === :MBR_log2_n_scans_ratio_true   ? :MBR_log2_n_scans_ratio_false :
    f === :MBR_best_irt_diff_true        ? :MBR_best_irt_diff_false :
    f === :MBR_best_irt_diff_worst_true  ? :MBR_best_irt_diff_worst_false :
    f === :MBR_log_by_diff_true          ? :MBR_log_by_diff_false :
    f === :MBR_smoothed_frag_hellinger_true ? :MBR_smoothed_frag_hellinger_false :
    f === :MBR_smoothed_frag_hellinger_worst_true ? :MBR_smoothed_frag_hellinger_worst_false :
    f === :MBR_donor_library_hellinger_true ? :MBR_donor_library_hellinger_false :
    f === :MBR_donor_library_hellinger_worst_true ? :MBR_donor_library_hellinger_worst_false :
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

function _mbr_transfer_training_mask(
    positive_top::AbstractVector{Bool},
    target_top::AbstractVector{Bool},
)
    n = length(target_top)
    mask = trues(2 * n)
    @inbounds for i in 1:n
        mask[i] = !target_top[i] || positive_top[i]
    end
    return mask
end

function _mbr_transfer_iteration_metrics(
    scores_double::AbstractVector{<:Real},
    eval_labels::AbstractVector{Bool},
    target_top::AbstractVector{Bool},
    n_cand::Int;
    ftr_threshold::Float32 = MBR_SEMISUPERVISED_FTR_THRESHOLD,
    fdr_threshold::Float32 = MBR_SEMISUPERVISED_FDR_THRESHOLD,
)
    qvals_double = Vector{Float32}(undef, 2 * n_cand)
    pep_double = Vector{Float32}(undef, 2 * n_cand)
    get_qvalues!(scores_double, eval_labels, qvals_double)
    get_PEP!(scores_double, eval_labels, pep_double)

    scores_top = scores_double[1:n_cand]
    fdr_qvals_top = Vector{Float32}(undef, n_cand)
    get_qvalues!(scores_top, target_top, fdr_qvals_top)

    positive_top = falses(n_cand)
    @inbounds for i in 1:n_cand
        positive_top[i] = target_top[i] &&
                          pep_double[i] <= ftr_threshold &&
                          fdr_qvals_top[i] <= fdr_threshold
    end

    return (
        qvals_double = qvals_double,
        pep_double = pep_double,
        qvals_top = qvals_double[1:n_cand],
        pep_top = pep_double[1:n_cand],
        fdr_qvals_top = fdr_qvals_top,
        positive_top = positive_top,
        n_positive = count(positive_top),
    )
end

function _mbr_candidate_gate(
    psms::DataFrame;
    q_thresh::Float32 = 0.01f0,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
)
    n = nrow(psms)
    pre_score = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    rescue_mask = hasproperty(psms, :mbr_rescue_candidate) ?
                  Bool.(psms[!, :mbr_rescue_candidate]) :
                  falses(n)

    pre_qvals = fill(1.0f0, n)
    normal_idx = findall(.!rescue_mask)
    if !isempty(normal_idx)
        normal_qvals = Vector{Float32}(undef, length(normal_idx))
        get_qvalues!(pre_score[normal_idx], target_col[normal_idx], normal_qvals)
        @inbounds for (j, i) in enumerate(normal_idx)
            pre_qvals[i] = normal_qvals[j]
        end
    end

    donor_target_pass = (.!rescue_mask) .&
                        (pre_qvals .<= donor_q_threshold) .&
                        target_col
    prob_thresh = any(donor_target_pass) ? minimum(pre_score[donor_target_pass]) : Float32(Inf)

    mbr_pp_t       = Float32.(psms[!, :MBR_max_pair_prob_true])
    mbr_miss_t     = Bool.(psms[!, :MBR_is_missing_true])
    mbr_miss_f     = Bool.(psms[!, :MBR_is_missing_false])
    candidate_mask = (pre_qvals .> q_thresh) .&
                     (mbr_pp_t  .>= prob_thresh) .&
                     (.!mbr_miss_t) .&
                     (.!mbr_miss_f)

    return (
        candidate_mask = candidate_mask,
        pre_qvals = pre_qvals,
        prob_thresh = prob_thresh,
        rescue_mask = rescue_mask,
    )
end

"""
    apply_mbr_filter_paired!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

MBR Batch F — paired-counterfactual FTR controller.

Algorithm (see BATCH_F_PLAN.md):
1. Candidate gate = (pre_qvals > q_thresh) & (MBR_max_pair_prob_true ≥
   prob_thresh) & (!MBR_is_missing_true) & (!MBR_is_missing_false).
2. Build a row-doubled training frame:
   - top half: candidate rows with FTR_FEATURES_F_TRUE values.
   - bottom half: same receivers with FTR_FEATURES_F_FALSE values.
3. Train iterative 2-fold CV LightGBM. Each iteration uses the previous round's
   3% FTR ∩ 1% top-half FDR target rows as positives, plus top-half decoys and
   counterfactual rows as negatives.
4. Compute q-values and PEPs on the doubled frame. Recover candidate rows whose
   top-half row has q-value ≤ alpha.

Sets `:mbr_recovered`, `:MBR_transfer_candidate`, `:mbr_target_decoy_prob`,
`:ftr_qval_true`, and `:ftr_pep_true` on `psms` in-place. Does NOT mutate
`:trace_prob`.
"""
function apply_mbr_filter_paired!(
    psms::DataFrame;
    alpha::Float32    = 0.01f0,
    q_thresh::Float32 = 0.01f0,
    pregated::Bool = false,
    pregated_prob_thresh::Float32 = Float32(NaN),
)
    t0 = time()

    n = nrow(psms)
    target_col = Bool.(psms[!, :target])
    if n == 0
        @debug_l1 "MBR Batch F — apply_mbr_filter_paired!: no PSMs to filter"
        psms[!, :mbr_recovered]          = falses(0)
        psms[!, :MBR_transfer_candidate] = falses(0)
        psms[!, :mbr_target_decoy_prob]  = Float32[]
        psms[!, :ftr_qval_true]          = Float32[]
        psms[!, :ftr_pep_true]           = Float32[]
        return (n_candidates=0, threshold=Float32(Inf), n_recovered=0, elapsed_s=0.0)
    end

    # ── 2. Candidates: failed q-cut, have BOTH _true and _false donors ──
    gate = pregated ? nothing : _mbr_candidate_gate(psms; q_thresh = q_thresh)
    candidate_mask = pregated ? trues(n) : gate.candidate_mask
    rescue_mask = hasproperty(psms, :mbr_rescue_candidate) ?
                  Bool.(psms[!, :mbr_rescue_candidate]) :
                  falses(n)
    prob_thresh = pregated ? pregated_prob_thresh : gate.prob_thresh
    n_cand = count(candidate_mask)
    n_rescue = count(rescue_mask)
    n_rescue_cand = count(candidate_mask .& rescue_mask)
    n_normal_cand = n_cand - n_rescue_cand
    cand_idx = findall(candidate_mask)

    @debug_l1 "MBR Batch F — paired FTR recovery:"
    @debug_l1 "  MBR_DONOR_Q_THRESHOLD = $MBR_DONOR_Q_THRESHOLD; prob_thresh = $(round(prob_thresh, digits=4))"
    if pregated
        @debug_l1 "  candidates (pre-gated sparse rows): $n_cand / $n"
    else
        @debug_l1 "  candidates (pre-MBR q>$q_thresh + both donors): $n_cand / $n ($(round(100*n_cand/n, digits=2))%)"
    end
    @debug_l1 "  rescue rows: $n_rescue; candidates normal=$n_normal_cand rescue=$n_rescue_cand"
    @debug_l1 "  α (q-value FTR budget): $alpha"

    if n_cand == 0
        psms[!, :mbr_recovered]          = falses(n)
        psms[!, :MBR_transfer_candidate] = candidate_mask
        psms[!, :mbr_target_decoy_prob]  = fill(NaN32, n)
        psms[!, :ftr_qval_true]          = fill(NaN32, n)
        psms[!, :ftr_pep_true]           = fill(NaN32, n)
        return (
            n_candidates = 0,
            threshold = Float32(Inf),
            n_recovered = 0,
            prob_thresh = prob_thresh,
            elapsed_s = time() - t0,
        )
    end

    # ── 3. Build the doubled training frame ──
    # Top half: true-donor rows. Label = true only for target receivers.
    # Bottom half: same receivers with _false MBR values swapped in. Label = false.
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
    target_top = Bool.(sub[!, :target])
    eval_labels = _mbr_transfer_training_labels(target_top)
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
            eval_labels,
            target_top,
            n_cand,
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
        if best_state === nothing || state.n_positive > best_state.n_positive
            best_state = state
        end

        @debug_l1 "  MBR transfer semi-supervised iter $iter_idx: " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "FTR≤$(round(100 * MBR_SEMISUPERVISED_FTR_THRESHOLD, digits=2))% " *
                   "∩ FDR≤$(round(100 * MBR_SEMISUPERVISED_FDR_THRESHOLD, digits=2))% " *
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
        training_labels = _mbr_transfer_training_labels(metrics.positive_top)
        training_mask = _mbr_transfer_training_mask(metrics.positive_top, target_top)
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
    qvals_top = qvals_double[1:n_cand]
    pep_top   = pep_double[1:n_cand]
    recovered_in_cand = qvals_top .<= alpha
    n_recovered = count(recovered_in_cand)

    mbr_recovered_full = falses(n)
    target_decoy_prob_full = fill(NaN32, n)
    ftr_qval_full      = fill(NaN32, n)
    ftr_pep_full       = fill(NaN32, n)
    ftr_score_top      = ftr_score_double[1:n_cand]
    @inbounds for (k, i) in enumerate(cand_idx)
        ftr_qval_full[i] = qvals_top[k]
        ftr_pep_full[i]  = pep_top[k]
        if recovered_in_cand[k]
            mbr_recovered_full[i] = true
            target_decoy_prob_full[i] = ftr_score_top[k]
        end
    end
    psms[!, :mbr_recovered]          = mbr_recovered_full
    psms[!, :MBR_transfer_candidate] = candidate_mask
    psms[!, :mbr_target_decoy_prob]  = target_decoy_prob_full
    psms[!, :ftr_qval_true]          = ftr_qval_full
    psms[!, :ftr_pep_true]           = ftr_pep_full

    # ── 7. Threshold for log: lowest FTR score on a recovered top-half row ──
    τ = n_recovered > 0 ?
        minimum(ftr_score_top[recovered_in_cand]) :
        Float32(Inf)

    # ── 8. Recovered transfer composition ──
    # Recovery uses the real donor for the same precursor. Counterfactual rows
    # only train/control FTR, so recovered composition is just target vs decoy
    # receiver rows.
    n_t_rec = count(i -> mbr_recovered_full[i] && target_col[i], 1:n)
    n_d_rec = n_recovered - n_t_rec
    n_rescue_recovered = count(i -> mbr_recovered_full[i] && rescue_mask[i], 1:n)
    n_normal_recovered = n_recovered - n_rescue_recovered

    @debug_l1 "  doubled-frame rows: $(2*n_cand)  (top=true donor, bottom=counterfactual)"
    @debug_l1 "  τ (FTR score, q ≤ α): $(round(τ, digits=4))"
    @debug_l1 "  RECOVERED (top-half q ≤ α): $n_recovered ($(round(100*n_recovered/max(n_cand,1), digits=2))% of candidates)"
    @debug_l1 "  recovered normal=$n_normal_recovered rescue=$n_rescue_recovered"
    @debug_l1 "  recovered targets: $n_t_rec   recovered decoys: $n_d_rec"

    if last_cls !== nothing
        lgbm_model = LightGBMModel(last_cls, available_true, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            @debug_l1 "  Batch F transfer-model feature importances (gain):"
            for (fname, gain) in sorted_imp
                @debug_l1 "    $(rpad(string(fname), 36)) $(round(gain, digits=2))"
            end
        end
    end

    elapsed = time() - t0
    @debug_l1 "  apply_mbr_filter_paired! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_candidates = n_cand,
        threshold = τ,
        n_recovered = n_recovered,
        prob_thresh = prob_thresh,
        elapsed_s = elapsed,
    )
end
