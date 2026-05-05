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

# Features for the FTR classifier. Excludes :target, :decoy, :MBR_is_best_decoy
# (these define the label, would leak). :MBR_is_missing is constant on the
# candidate cohort by construction (candidates require .!MBR_is_missing).
const FTR_FEATURES = Symbol[
    :trace_prob_mbr,
    :trace_prob_prepass,
    :MBR_max_pair_prob,
    :MBR_log2_weight_ratio,
    :MBR_log2_explained_ratio,
]

"""
Internal: train a 2-fold CV LightGBM on the candidate cohort to predict
`is_good_transfer = .!is_bad_transfer`. Returns OOF scores aligned with
the input PSM rows (NaN32 outside the candidate set).
"""
function _train_ftr_lgbm(
    psms::DataFrame,
    candidate_mask::AbstractVector{Bool},
    is_bad_transfer::AbstractVector{Bool};
    features::Vector{Symbol} = FTR_FEATURES,
)
    n = nrow(psms)
    full_scores = fill(NaN32, n)
    cand_idx = findall(candidate_mask)
    isempty(cand_idx) && return full_scores, nothing, Symbol[]

    available = filter(f -> hasproperty(psms, f), features)
    if isempty(available)
        @user_warn "FTR LightGBM: no FTR features available; falling back to :trace_prob_mbr"
        @inbounds for i in cand_idx
            full_scores[i] = Float32(psms[i, :trace_prob_mbr])
        end
        return full_scores, nothing, Symbol[]
    end

    sub = psms[cand_idx, :]
    X_sub = feature_matrix(sub, available)
    y_good_sub = .!is_bad_transfer[cand_idx]

    # 2-fold CV on existing cv_fold (positions are within the candidate subset)
    cv_fold_sub = sub[!, :cv_fold]
    pos0 = findall(cv_fold_sub .== 0)
    pos1 = findall(cv_fold_sub .== 1)

    sub_scores = fill(NaN32, length(cand_idx))
    last_cls = nothing
    for (train_pos, test_pos) in ((pos1, pos0), (pos0, pos1))
        (isempty(train_pos) || isempty(test_pos)) && continue
        y_tr = y_good_sub[train_pos]
        if length(unique(y_tr)) == 1
            sub_scores[test_pos] .= y_tr[1] ? 1f0 : 0f0
            continue
        end
        cls = build_lightgbm_classifier(; SHARED_LGBM_HP...)
        LightGBM.fit!(cls, X_sub[train_pos, :], _prepare_labels(y_tr); verbosity = -1)
        raw = LightGBM.predict(cls, X_sub[test_pos, :])
        s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
        sub_scores[test_pos] .= Float32.(s)
        last_cls = cls
    end

    @inbounds for (k, i) in enumerate(cand_idx)
        full_scores[i] = sub_scores[k]
    end
    return full_scores, last_cls, available
end

"""
Internal: count transfer-class composition (T←T, T←D, D←T, D←D) on a mask.
Donor polarity is read from :MBR_is_best_decoy (true = donor is decoy).
Returns a NamedTuple of counts.
"""
function _transfer_class_counts(
    mask::AbstractVector{Bool},
    target::AbstractVector{Bool},
    donor_is_decoy::AbstractVector{Bool},
)
    n_tt = 0; n_td = 0; n_dt = 0; n_dd = 0
    @inbounds for i in eachindex(mask)
        mask[i] || continue
        is_t = target[i]
        donor_d = donor_is_decoy[i]
        if is_t
            donor_d ? (n_td += 1) : (n_tt += 1)
        else
            donor_d ? (n_dd += 1) : (n_dt += 1)
        end
    end
    return (TT = n_tt, TD = n_td, DT = n_dt, DD = n_dd)
end

_pct(num, den) = den == 0 ? 0.0 : round(100 * num / den, digits = 2)

function _log_transfer_counts(label::String, c::NamedTuple)
    total = c.TT + c.TD + c.DT + c.DD
    @user_info "  $label transfers — total=$total  T←T=$(c.TT) ($(_pct(c.TT, total))%)  T←D=$(c.TD) ($(_pct(c.TD, total))%)  D←T=$(c.DT) ($(_pct(c.DT, total))%)  D←D=$(c.DD) ($(_pct(c.DD, total))%)"
end

"""
    apply_mbr_filter!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

Sets `:mbr_recovered::Bool` and `:MBR_transfer_candidate::Bool` on `psms`
in-place. Does NOT mutate `:trace_prob` (which holds the non-MBR score
that drives the downstream qval pipeline).

Phase 6: trains an FTR LightGBM on the candidate cohort (label =
is_bad_transfer; features = FTR_FEATURES) and applies the FTR threshold
to its OOF scores instead of the raw :trace_prob_mbr.

Required columns on `psms`:
- :trace_prob_prepass — non-MBR pass-1 score
- :trace_prob_mbr     — MBR-boosted pass-2 score
- :target, :decoy, :cv_fold
- :MBR_max_pair_prob, :MBR_is_missing, :MBR_is_best_decoy,
  :MBR_log2_weight_ratio, :MBR_log2_explained_ratio
"""
function apply_mbr_filter!(
    psms::DataFrame;
    alpha::Float32   = 0.01f0,
    q_thresh::Float32 = 0.01f0,
)
    t0 = time()

    n = nrow(psms)
    if n == 0
        @user_info "MBR Phase 6 — apply_mbr_filter!: no PSMs to filter"
        psms[!, :mbr_recovered]          = falses(0)
        psms[!, :MBR_transfer_candidate] = falses(0)
        psms[!, :ftr_score]              = Float32[]
        return (n_candidates=0, n_bad=0, threshold=Float32(Inf), n_recovered=0, elapsed_s=0.0)
    end

    # ── 1. Pre-MBR q-values from the non-MBR (pass-1) score ──
    pre_score  = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    decoy_col  = Bool.(psms[!, :decoy])
    pre_qvals  = Vector{Float32}(undef, n)
    get_qvalues!(pre_score, target_col, pre_qvals)

    target_pass = (pre_qvals .<= q_thresh) .& target_col
    prob_thresh = if any(target_pass)
        minimum(pre_score[target_pass])
    else
        Float32(Inf)
    end

    # ── 2. Candidates: failed pre-MBR q-cut AND have a strong cross-run donor ──
    mbr_pp         = Float32.(psms[!, :MBR_max_pair_prob])
    mbr_missing    = Bool.(psms[!, :MBR_is_missing])
    mbr_best_decoy = Bool.(psms[!, :MBR_is_best_decoy])
    candidate_mask = (pre_qvals .> q_thresh) .&
                     (mbr_pp .>= prob_thresh) .&
                     (.!mbr_missing)
    n_candidates = count(candidate_mask)

    # ── 3. is_bad_transfer: donor signals a swap ──
    is_bad_transfer = candidate_mask .&
                      ((target_col .& mbr_best_decoy) .|
                       (decoy_col  .& .!mbr_best_decoy))
    n_bad = count(is_bad_transfer)

    # ── 4. Diagnostic: transfer composition in candidate (training) set ──
    cand_counts = _transfer_class_counts(candidate_mask, target_col, mbr_best_decoy)

    # ── 5. Train FTR LightGBM on candidates (2-fold CV → OOF scores) ──
    ftr_score, ftr_cls, ftr_feats = if n_candidates > 0
        _train_ftr_lgbm(psms, candidate_mask, is_bad_transfer)
    else
        fill(NaN32, n), nothing, Symbol[]
    end
    psms[!, :ftr_score] = ftr_score

    # ── 6. FTR threshold τ at α on the FTR model's OOF score over candidates ──
    τ = if n_candidates > 0
        get_ftr_threshold(ftr_score, is_bad_transfer, alpha; mask = candidate_mask)
    else
        Float32(Inf)
    end

    # ── 7. Recovery: candidate AND ftr_score ≥ τ AND NOT is_bad_transfer ──
    recovered = candidate_mask .& (ftr_score .>= τ) .& (.!is_bad_transfer)
    n_recovered = count(recovered)
    rec_counts = _transfer_class_counts(recovered, target_col, mbr_best_decoy)

    psms[!, :mbr_recovered]          = recovered
    psms[!, :MBR_transfer_candidate] = candidate_mask

    elapsed = time() - t0

    pct_cand = round(100 * n_candidates / max(n, 1), digits = 2)
    pct_bad  = n_candidates == 0 ? 0.0 : round(100 * n_bad / n_candidates, digits = 2)
    pct_rec  = n_candidates == 0 ? 0.0 : round(100 * n_recovered / n_candidates, digits = 2)
    @user_info "MBR Phase 6 — FTR recovery (LightGBM on candidate cohort):"
    @user_info "  α = $alpha;  pre-MBR q_thresh = $q_thresh;  prob_thresh = $(round(prob_thresh, digits = 4))"
    @user_info "  candidates (pre-MBR q>$q_thresh + donor): $n_candidates / $n ($pct_cand% of PSMs)"
    @user_info "  is_bad_transfer in candidates: $n_bad ($pct_bad% of candidates)"
    @user_info "  FTR features used: $(ftr_feats)"
    _log_transfer_counts("training (candidates)", cand_counts)
    @user_info "  τ (FTR-model score, FTR α): $(round(τ, digits = 4))"
    @user_info "  RECOVERED (ftr_score ≥ τ, not bad): $n_recovered ($pct_rec% of candidates)"
    _log_transfer_counts("recovered", rec_counts)

    if ftr_cls !== nothing
        lgbm_model = LightGBMModel(ftr_cls, ftr_feats, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            @user_info "  FTR-model feature importances (gain):"
            for (fname, gain) in sorted_imp
                @user_info "    $(rpad(string(fname), 32)) $(round(gain, digits = 2))"
            end
        end
    end
    @user_info "  apply_mbr_filter! elapsed: $(round(elapsed, digits = 2))s"

    return (
        n_candidates    = n_candidates,
        n_bad           = n_bad,
        threshold       = τ,
        n_recovered     = n_recovered,
        prob_thresh     = prob_thresh,
        cand_counts     = cand_counts,
        recovered_counts = rec_counts,
        elapsed_s       = elapsed,
    )
end
