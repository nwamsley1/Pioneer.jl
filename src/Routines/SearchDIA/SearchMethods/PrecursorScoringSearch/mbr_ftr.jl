# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR Phase 4 (Phase 5b reframing): False-Transfer-Rate controller as a
# RECOVERY mechanism scoped to the transfer-candidate cohort.
#
# Architecture:
#   * The qval pipeline downstream uses :trace_prob (the NON-MBR pass-1
#     score). Non-candidate rows can never be score-inflated by MBR.
#   * Transfer candidates (rows that fail q≤0.01 in pre-MBR + have a strong
#     cross-run donor) get evaluated by their :trace_prob_mbr (MBR-boosted)
#     against an FTR threshold τ at α=0.01. Recovered candidates get
#     :mbr_recovered = true; PrecursorScoringSearch's qval filter then
#     bypasses the per-file qval cut for these rows.
#   * Recovered candidates carry their :trace_prob (non-MBR, low) into
#     downstream aggregation so their quant weight reflects per-run evidence
#     rather than a phantom MBR boost.

"""
    apply_mbr_filter!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

Sets `:mbr_recovered::Bool` and `:MBR_transfer_candidate::Bool` on `psms`
in-place. Does NOT mutate `:trace_prob` (which holds the non-MBR score
that drives the downstream qval pipeline).

Required columns on `psms`:
- :trace_prob_prepass — non-MBR pass-1 score (also = :trace_prob in Phase 5b)
- :trace_prob_mbr     — MBR-boosted pass-2 score (used by FTR only)
- :target, :decoy
- :MBR_max_pair_prob, :MBR_is_missing, :MBR_is_best_decoy

Returns a NamedTuple with diagnostic counts.
"""
function apply_mbr_filter!(
    psms::DataFrame;
    alpha::Float32   = 0.01f0,
    q_thresh::Float32 = 0.01f0,
)
    t0 = time()

    n = nrow(psms)
    if n == 0
        @user_info "MBR Phase 4 — apply_mbr_filter!: no PSMs to filter"
        psms[!, :mbr_recovered]          = falses(0)
        psms[!, :MBR_transfer_candidate] = falses(0)
        return (n_candidates=0, n_bad=0, threshold=Float32(Inf), n_recovered=0, elapsed_s=0.0)
    end

    # ── 1. Pre-MBR q-values from the non-MBR (pass-1) score ──
    pre_score  = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    decoy_col  = Bool.(psms[!, :decoy])
    pre_qvals  = Vector{Float32}(undef, n)
    get_qvalues!(pre_score, target_col, pre_qvals)

    # Probability threshold = lowest pre-MBR score among targets passing q_thresh.
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

    # ── 3. is_bad_transfer: candidate where the donor signals a swap ──
    is_bad_transfer = candidate_mask .&
                      ((target_col .& mbr_best_decoy) .|
                       (decoy_col  .& .!mbr_best_decoy))
    n_bad = count(is_bad_transfer)

    # ── 4. FTR threshold τ at α on the MBR-BOOSTED score over candidates ──
    boosted = Float32.(psms[!, :trace_prob_mbr])
    τ = if n_candidates > 0
        get_ftr_threshold(boosted, is_bad_transfer, alpha; mask=candidate_mask)
    else
        Float32(Inf)
    end

    # ── 5. Recovery: candidate AND boosted ≥ τ AND NOT is_bad_transfer ──
    recovered = candidate_mask .& (boosted .>= τ) .& (.!is_bad_transfer)
    n_recovered = count(recovered)

    psms[!, :mbr_recovered]          = recovered
    psms[!, :MBR_transfer_candidate] = candidate_mask

    elapsed = time() - t0

    pct_cand = round(100 * n_candidates / max(n, 1), digits=2)
    pct_bad  = n_candidates == 0 ? 0.0 : round(100 * n_bad / n_candidates, digits=2)
    pct_rec  = n_candidates == 0 ? 0.0 : round(100 * n_recovered / n_candidates, digits=2)
    @user_info "MBR Phase 4/5b — FTR recovery:"
    @user_info "  α = $alpha;  pre-MBR q_thresh = $q_thresh;  prob_thresh = $(round(prob_thresh, digits=4))"
    @user_info "  candidates (pre-MBR q>$q_thresh + donor): $n_candidates / $n ($pct_cand% of PSMs)"
    @user_info "  is_bad_transfer in candidates: $n_bad ($pct_bad% of candidates)"
    @user_info "  τ (MBR-boosted score, FTR α): $(round(τ, digits=4))"
    @user_info "  RECOVERED (boosted ≥ τ, not bad): $n_recovered ($pct_rec% of candidates)"
    @user_info "  apply_mbr_filter! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_candidates = n_candidates,
        n_bad        = n_bad,
        threshold    = τ,
        n_recovered  = n_recovered,
        prob_thresh  = prob_thresh,
        elapsed_s    = elapsed,
    )
end
