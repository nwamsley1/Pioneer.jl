**Simpler Pairing Strategy — Critique And Alternative Plan**

This document proposes a simpler approach to pairing for MBR that achieves the same downstream effect without cloning rows or regenerating pair_id at scoring time. It also critiques the prior “simple” pairing plan (not in repo here) that introduced additional complexity.

**Critique Of The Prior Plan**

- Complexity and memory cost: Cloning decoy (or target) rows to enforce 1:1 pair_id multiplies the PSM table size and introduces quadratic‑like costs if clones are appended iteratively.
- Moving targets: Re‑pairing at scoring time overrides a stable library invariant. The spectral library already encodes target/decoy pairing (e.g., `pair_id`) upstream; redefining pairs late in the pipeline makes debugging harder and risks drift between runs.
- Coupling to downstream logic: Enforcing 1:1 at the row level only to satisfy a specific grouping in `apply_mbr_filter!` creates tight coupling. A better approach is to compute dominance flags directly and keep `apply_mbr_filter!` agnostic to strict 1:1 groups.
- CV‑fold nuances: Pairing within cv_fold and bins is brittle. The scorer already computes robust MBR features per run; the pairing step is redundant if we express “decoy outranks target” with per‑row flags.

**Goal**

- Preserve the practical outcome (flag and control risky transfers) without cloning or regenerating pair_id.
- Reduce CPU/memory use, avoid write‑amplification, and keep the data model stable.

**Key Observation**

- The MBR filter does not inherently require 1:1 pair groups if we provide per‑row dominance flags. Earlier we already compute:
  - `MBR_max_pair_prob` and `MBR_is_best_decoy` in `summarize_precursors!`.
  - A candidate mask in the scorer.
  - `apply_mbr_filter!` can operate on those signals without grouping.

**Proposed Simpler Design**

1) Do not regenerate pair_id at scoring time
- Keep the library’s `pair_id` (if present) as a passive attribute. Don’t modify it or rely on it for strict 1:1.
- If `pair_id` is missing in some PSMs, that’s okay — the steps below don’t depend on it.

2) Fix candidate labeling in the in‑memory scorer
- File: `src/utils/ML/percolatorSortOf.jl`
- Function: `sort_of_percolator_in_memory!`
- Replace the probability‑based pass mask with q‑value based logic (mirroring `update_mbr_probs!`):

```julia
qvals_prev = similar(nonMBR_estimates)
get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)
pass_mask = (qvals_prev .<= max_q_value_xgboost_rescore) .& psms.target
prob_thresh = any(pass_mask) ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)
psms[!, :MBR_transfer_candidate] .= (qvals_prev .> max_q_value_xgboost_rescore) .&
                                    (psms.MBR_max_pair_prob .>= prob_thresh)
```

3) Compute per‑row dominance without cloning
- `summarize_precursors!` already determines, per run, a “best other run” via `MBR_max_pair_prob` and flags `MBR_is_best_decoy`.
- If we want a single, direct dominance flag for filtering, compute it after final probabilities with a single grouped combine over candidates only (no cloning, no pair regeneration):

```julia
agg = combine(groupby(view(psms, psms.MBR_transfer_candidate, :),
                      [:ms_file_idx, :pair_id, :target]),
              :prob => maximum => :max_prob)
tgt = rename!(agg[agg.target .== true, [:ms_file_idx, :pair_id, :max_prob]], :max_prob => :tmax)
dcy = rename!(agg[agg.target .== false, [:ms_file_idx, :pair_id, :max_prob]], :max_prob => :dmax)
pairmax = outerjoin(tgt, dcy, on=[:ms_file_idx, :pair_id])
pairmax[!, :MBR_paired_decoy_higher] = coalesce.(pairmax.dmax, -Inf32) .> coalesce.(pairmax.tmax, -Inf32)
psms = leftjoin(psms, pairmax[:, [:ms_file_idx, :pair_id, :MBR_paired_decoy_higher]], on=[:ms_file_idx, :pair_id])
psms[!, :MBR_paired_decoy_higher] = ifelse.(psms.target, coalesce.(psms.MBR_paired_decoy_higher, false), false)
```

Notes:
- This is optional if `MBR_is_best_decoy` is already reliable for filtering.
- It avoids cloning and only touches the minimal set of rows.

4) Keep `apply_mbr_filter!` simple
- File: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`
- Function: `apply_mbr_filter!`
- Use the candidate mask and a simple bad‑transfer predicate (original form):

```julia
candidate_mask = merged_df.MBR_transfer_candidate
is_bad_transfer = candidate_mask .& (
    (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .|
    merged_df.decoy
)
τ = get_ftr_threshold(merged_df.prob, merged_df.target,
                      is_bad_transfer, params.max_MBR_false_transfer_rate;
                      mask=candidate_mask)
merged_df._filtered_prob = ifelse.(candidate_mask .& (merged_df.prob .< τ), 0f0, merged_df.prob)
```

This design works whether `pair_id` exists or not and does not require strict 1:1 pairing.

**Benefits**

- No cloning or re‑pairing: Lower memory/CPU, simpler control flow, smaller risk of regression.
- Stable semantics: Relies on per‑row flags and q‑value logic already present.
- Decoupled filter: `apply_mbr_filter!` needs only a candidate mask and a bad‑transfer mask — no grouping gymnastics.

**Implementation Notes**

- Instrumentation: Log candidate counts and `prob_thresh` after labeling to confirm realistic set sizes.
- Out‑of‑memory parity: The OOM path already uses correct q‑value‑based labeling (`update_mbr_probs!`). With this plan, both paths are aligned semantically.

