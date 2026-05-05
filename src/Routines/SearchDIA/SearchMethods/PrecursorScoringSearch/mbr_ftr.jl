# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR Phase 4: False-Transfer-Rate (FTR) controller on the MBR-boosted score.
#
# After the second-pass LightGBM has produced :trace_prob (the MBR-boosted
# score) and compute_mbr_features! populated :MBR_max_pair_prob /
# :MBR_is_best_decoy / :MBR_is_missing, this module:
#
# 1. Defines `:MBR_transfer_candidate` from q-values of the *non-MBR*
#    `:trace_prob_prepass` (NOT the boosted score — see
#    docs/src/advanced/mbr_in_memory_candidate_labeling.md from v0.6.6 for
#    the bug history that made this distinction load-bearing).
# 2. Defines `is_bad_transfer = candidate & ((target & best_decoy) |
#    (decoy & !best_decoy))` — captures T←D and D←T swaps.
# 3. Computes the FTR threshold τ at α (default 0.01) over the candidate
#    set scored by the boosted `:trace_prob`.
# 4. Zeros `:trace_prob` for rows that are candidates AND (score < τ OR
#    is_bad_transfer). This must happen BEFORE aggregate_per_file! so the
#    zeros propagate into prec_prob.

"""
    apply_mbr_filter!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

In-place FTR control on `psms[:, :trace_prob]` (the MBR-boosted score).
See module docstring for details.

Required columns on `psms`:
- :trace_prob, :trace_prob_prepass, :target, :decoy
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
        return (n_candidates=0, n_bad=0, threshold=Float32(Inf), n_filtered=0, elapsed_s=0.0)
    end

    # ── 1. Pre-MBR q-values from the first-pass score ──
    pre_score = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    decoy_col  = Bool.(psms[!, :decoy])
    pre_qvals = Vector{Float32}(undef, n)
    get_qvalues!(pre_score, target_col, pre_qvals)

    # Probability threshold = lowest pre-MBR score among targets passing q_thresh.
    target_pass = (pre_qvals .<= q_thresh) .& target_col
    prob_thresh = if any(target_pass)
        minimum(pre_score[target_pass])
    else
        Float32(Inf)
    end

    # ── 2. Candidates: failed pre-MBR q-cut AND have a strong cross-run donor ──
    mbr_pp        = Float32.(psms[!, :MBR_max_pair_prob])
    mbr_missing   = Bool.(psms[!, :MBR_is_missing])
    mbr_best_decoy = Bool.(psms[!, :MBR_is_best_decoy])
    candidate_mask = (pre_qvals .> q_thresh) .&
                     (mbr_pp .>= prob_thresh) .&
                     (.!mbr_missing)
    n_candidates = count(candidate_mask)

    # ── 3. is_bad_transfer: candidate where the donor signals a swap ──
    #   T←D : target row whose donor was a decoy (false positive transfer)
    #   D←T : decoy row whose donor was a target (decoy got "rescued" by a target)
    is_bad_transfer = candidate_mask .&
                      ((target_col .& mbr_best_decoy) .|
                       (decoy_col  .& .!mbr_best_decoy))
    n_bad = count(is_bad_transfer)

    # ── 4. FTR threshold τ at α on the MBR-boosted score, restricted to candidates ──
    boosted = Float32.(psms[!, :trace_prob])
    τ = if n_candidates > 0
        get_ftr_threshold(boosted, is_bad_transfer, alpha; mask=candidate_mask)
    else
        Float32(Inf)
    end

    # ── 5. Zero out trace_prob for candidates that fail the gate ──
    #   (a) score < τ  OR  (b) is_bad_transfer
    filter_mask = candidate_mask .& ((boosted .< τ) .| is_bad_transfer)
    n_filtered = count(filter_mask)
    @inbounds for i in 1:n
        if filter_mask[i]
            psms[i, :trace_prob] = 0.0f0
        end
    end

    # Persist the candidate flag for downstream / diagnostics.
    psms[!, :MBR_transfer_candidate] = candidate_mask

    elapsed = time() - t0

    pct_cand = round(100 * n_candidates / max(n, 1), digits=2)
    pct_bad = n_candidates == 0 ? 0.0 : round(100 * n_bad / n_candidates, digits=2)
    pct_filt = n_candidates == 0 ? 0.0 : round(100 * n_filtered / n_candidates, digits=2)
    @user_info "MBR Phase 4 — FTR control:"
    @user_info "  α = $alpha;  pre-MBR q_thresh = $q_thresh;  prob_thresh = $(round(prob_thresh, digits=4))"
    @user_info "  candidates: $n_candidates / $n ($(pct_cand)% of all PSMs)"
    @user_info "  bad transfers in candidates: $n_bad ($(pct_bad)% of candidates)"
    @user_info "  τ (boosted score): $(round(τ, digits=4))"
    @user_info "  zeroed trace_prob: $n_filtered ($(pct_filt)% of candidates)"
    @user_info "  apply_mbr_filter! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_candidates = n_candidates,
        n_bad        = n_bad,
        threshold    = τ,
        n_filtered   = n_filtered,
        prob_thresh  = prob_thresh,
        elapsed_s    = elapsed,
    )
end
