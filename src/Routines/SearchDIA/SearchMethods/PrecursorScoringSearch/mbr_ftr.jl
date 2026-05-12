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

# Features for the FTR classifier.
# Phase 8f: use the full ADVANCED_FEATURE_SET (already excludes :target /
# :decoy / :MBR_is_best_decoy by construction — those would leak the label).
# Plus the two pass-1/pass-2 scores (:trace_prob_prepass, :trace_prob_mbr)
# which aren't in ADVANCED_FEATURE_SET. The training helper filters to
# whatever columns actually exist on the candidate DataFrame at runtime,
# so listing missing columns is harmless.
const FTR_FEATURES = Symbol[
    :trace_prob_prepass,
    :trace_prob_mbr,
    # All ADVANCED_FEATURE_SET features — label-safe and rich:
    :missed_cleavage, :Mox, :prec_mz, :sequence_length, :charge,
    :irt_pred, :irt_error, :irt_diff,
    :max_y_ions, :y_ions_sum, :longest_y, :y_count, :b_count,
    :total_ions, :total_ions_iso,
    :best_rank, :best_rank_iso, :topn, :topn_iso,
    :gof, :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_residual, :max_unmatched_residual, :max_gof,
    :fitted_spectral_contrast, :spectral_contrast,
    :max_matched_ratio, :err_norm, :poisson, :weight,
    :log2_intensity_explained, :tic, :smoothness,
    :ms1_ms2_rt_diff,
    :gof_ms1, :max_matched_residual_ms1, :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1, :error_ms1, :m0_error_ms1,
    :n_iso_ms1, :big_iso_ms1, :rt_max_intensity_ms1,
    :rt_diff_max_intensity_ms1, :ms1_features_missing,
    :percent_theoretical_ignored, :scribe, :max_scribe, :max_weight,
    :fitted_hellinger, :n_scans,
    # MBR-specific
    :MBR_max_pair_prob, :MBR_log2_weight_ratio,
    :MBR_log2_explained_ratio, :MBR_best_irt_diff, :MBR_is_missing,
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

    # Donor score floor: lowest prepass score among targets passing
    # MBR_DONOR_Q_THRESHOLD. Decouples donor q-cut from candidate q-cut so we
    # can sweep donor stringency independently.
    donor_target_pass = (pre_qvals .<= MBR_DONOR_Q_THRESHOLD) .& target_col
    prob_thresh = if any(donor_target_pass)
        minimum(pre_score[donor_target_pass])
    else
        Float32(Inf)
    end
    @user_info "  MBR_DONOR_Q_THRESHOLD = $MBR_DONOR_Q_THRESHOLD; prob_thresh = $(round(prob_thresh, digits=4))"

    # ── 2. Candidates: failed pre-MBR q-cut AND have a strong cross-run donor ──
    mbr_pp         = Float32.(psms[!, :MBR_max_pair_prob])
    mbr_missing    = Bool.(psms[!, :MBR_is_missing])
    mbr_best_decoy = Bool.(psms[!, :MBR_is_best_decoy])
    candidate_mask = (pre_qvals .> q_thresh) .&
                     (mbr_pp .>= prob_thresh) .&
                     (.!mbr_missing)
    n_candidates = count(candidate_mask)

    # ── 3. is_bad_transfer (Phase 8e — revert to 6/8a definition):
    # only T←D and D←T (cross-class swaps) count as bad. D←D is treated as
    # neutral (same-class noise, FDR-symmetric mirror of T←T).
    is_bad_transfer = candidate_mask .&
                      ((target_col .& mbr_best_decoy) .|     # T←D
                       (decoy_col  .& .!mbr_best_decoy))     # D←T
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

    # ── 7. Recovery (Phase 8c): candidate AND ftr_score ≥ τ.
    # No post-hoc .!is_bad_transfer hard-exclusion — the FTR α budget is the
    # sole gate. Among recovered, the empirical bad fraction is bounded by α
    # by construction of get_ftr_threshold. Lets the FTR LGBM (trained with
    # D←D as bad in Phase 8b) be the sole arbiter of which transfers to keep.
    recovered = candidate_mask .& (ftr_score .>= τ)
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
