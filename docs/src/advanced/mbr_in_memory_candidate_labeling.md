**In‑Memory MBR Candidate Labeling — Bug Analysis and Fix Proposal**

This note explains why the in‑memory path is labeling far too many rows as match‑between‑runs (MBR) transfer candidates and how to correct it. The focus is on the in‑memory scorer in `percolatorSortOf.jl` and how it differs from the out‑of‑memory path that uses the correct logic.

**Where Candidates Are Labeled (In‑Memory)**

- File: `src/utils/ML/percolatorSortOf.jl`
- Function: `sort_of_percolator_in_memory!`
- Location: near the end of the function, after cross‑validation models produce final probabilities and just before the MBR features and final write‑back.

Current code (abridged, as present in this branch):

```julia
# Determine which precursors failed the q-value cutoff prior to MBR
qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)
pass_mask = (nonMBR_estimates .<= max_q_value_xgboost_rescore)
prob_thresh = any(pass_mask) ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)

# Label as transfer candidates only those failing the q-value cutoff but
# whose best matched pair surpassed the passing probability threshold.
psms[!, :MBR_transfer_candidate] .= .!pass_mask .&
                                    (psms.MBR_max_pair_prob .>= prob_thresh)

# Use the final MBR probabilities for all precursors
psms[!, :prob] = MBR_estimates
```

Note: `nonMBR_estimates` are probabilities (higher=better). `max_q_value_xgboost_rescore` is a q‑value threshold (≈0.01) — a different scale.

**What’s Wrong (Scale Mismatch)**

- `pass_mask = (nonMBR_estimates .<= max_q_value_xgboost_rescore)` compares probabilities to a q‑value threshold. Since true positives have large probabilities (e.g., 0.8, 0.9), the condition `p ≤ 0.01` is almost never true.
  - Consequence: `pass_mask` is mostly false; `.!pass_mask` becomes “almost everyone”.
- `prob_thresh = minimum(nonMBR_estimates[pass_mask])` is taken over a tiny set (often empty). If any are present, they tend to be extremely small (≈0), so the condition `MBR_max_pair_prob ≥ prob_thresh` is trivially satisfied by most rows.
- Net effect: `MBR_transfer_candidate` is set to true for the vast majority of rows, which matches your observation:

```
FTR probability threshold: 0.8703515
Num passing candidate transfers: 1811822 out of 1949436
```

This happens because the candidate set (where `MBR_transfer_candidate=true`) includes almost every row, forcing the FTR threshold τ high to keep the empirical ratio below α.

**Correct Intent (Reference Implementation)**

- The out‑of‑memory path uses the correct logic inside `update_mbr_probs!` (same file):

```julia
function update_mbr_probs!(df::AbstractDataFrame, probs::AbstractVector{Float32}, qval_thresh::Float32)
    prev_qvals = similar(df.prob)
    get_qvalues!(df.prob, df.target, prev_qvals)   # compute q-values
    pass_mask = (prev_qvals .<= qval_thresh) .& df.target
    prob_thresh = any(pass_mask) ? minimum(df.prob[pass_mask]) : typemax(Float32)
    df[!, :MBR_transfer_candidate] = (prev_qvals .> qval_thresh) .&   # use q-values
                                     (df.MBR_max_pair_prob .>= prob_thresh)
    df[!, :prob] = probs
    return df
end
```

- Two key differences from the in‑memory code above:
  1) `pass_mask` is computed with q‑values `prev_qvals`, not with probabilities.
  2) `prob_thresh` is derived from the probabilities of the passing set (as intended), but the passing set is defined by `qval_thresh`, not a probability threshold.

**Impact of the Bug**

- Candidate set is massively inflated in the in‑memory path, causing:
  - Very high τ even at α = 0.01.
  - Large number of “candidate transfers” reported, which is misleading and slows down filtering.
  - Downstream clamping affects almost all rows, reducing the discriminative power of MBR.

**Minimal Fix (In‑Memory Path)**

Replace the probability‑based `pass_mask` with a q‑value‑based mask (as in `update_mbr_probs!`). Proposed snippet drop‑in for `sort_of_percolator_in_memory!`:

```julia
# 1) Compute q-values of non-MBR predictions
qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)

# 2) Build pass/fail by q-value threshold
pass_mask = (qvals_prev .<= max_q_value_xgboost_rescore) .& psms.target
prob_thresh = any(pass_mask) ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)

# 3) Label transfer candidates: failed q-value but paired to a strong donor
psms[!, :MBR_transfer_candidate] .= (qvals_prev .> max_q_value_xgboost_rescore) .&
                                    (psms.MBR_max_pair_prob .>= prob_thresh)
```

This aligns the in‑memory behavior with the out‑of‑memory `update_mbr_probs!` and ensures the candidate set is limited to plausible transfers.

**Verification Steps**

- Add a short diagnostic around labeling (in‑memory):
  - Count of `pass_mask`, `candidate_count = sum(MBR_transfer_candidate)`, and `prob_thresh`.
  - Sanity: candidate_count should be a minority of all rows; `prob_thresh` should be a realistic boundary (e.g., near the minimum probability among passing targets).

**Summary**

- The in‑memory path mistakenly compares probabilities to a q‑value threshold, inflating `MBR_transfer_candidate`.
- The out‑of‑memory path uses the correct q‑value based mask.
- Switching the in‑memory logic to use q‑values (as shown above) restores consistency and reduces the candidate set, leading to reasonable τ at α=0.01.

