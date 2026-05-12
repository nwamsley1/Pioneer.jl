# Transfer-Candidacy Requirements (Batch F, `apply_mbr_filter_paired!`)

A row in `best_psms` is a **transfer candidate** if it passes **all four**
of these conditions simultaneously. From `mbr_ftr.jl`:

```julia
candidate_mask = (pre_qvals .> q_thresh) .&        # (1)
                 (mbr_pp_t  .>= prob_thresh) .&    # (2)
                 (.!mbr_miss_t) .&                  # (3)
                 (.!mbr_miss_f)                     # (4)
```

## (1) The PSM **failed** the non-MBR per-PSM q-cut

`pre_qvals[i] > q_thresh`, where:

- `pre_qvals` are q-values from `get_qvalues!` on `:trace_prob_prepass`
  (the Pass-1 non-MBR LightGBM score), with `:target` as the positive
  label.
- `q_thresh = 0.01` (default arg to `apply_mbr_filter_paired!`).

A row that already passes q ≤ 0.01 on its own evidence is **not** a
transfer candidate — it doesn't need rescue. We only consider rows that
would otherwise be filtered out.

## (2) The TRUE donor is **confident enough**

`MBR_max_pair_prob_true[i] >= prob_thresh`, where `prob_thresh` is
computed once per run as:

```julia
donor_target_pass = (pre_qvals .<= MBR_DONOR_Q_THRESHOLD) .& target_col
prob_thresh = minimum(pre_score[donor_target_pass])
```

i.e., **the lowest `:trace_prob_prepass` among target rows that
themselves passed q ≤ `MBR_DONOR_Q_THRESHOLD`** (default `0.01`,
module-level constant in `mbr_ftr.jl`).

In other words, `prob_thresh` is the prepass-score floor corresponding
to "donor q ≤ 0.01 among targets". A row's TRUE donor (best score for
the same precursor in a file ≠ k) must clear this floor to count as
confident.

The donor q-threshold (`MBR_DONOR_Q_THRESHOLD`) is **decoupled** from
the candidate q-threshold (`q_thresh`) so they can be swept
independently. Today both are 0.01.

## (3) The TRUE donor **exists** (cross-run)

`MBR_is_missing_true[i] == false`, meaning: the same precursor
(`precursor_idx`) has at least one PSM in some file ≠ k (where k is
this row's `ms_file_idx`). Otherwise there is nothing to "transfer"
from.

## (4) The FALSE (counterfactual) donor **exists** (cross-run)

`MBR_is_missing_false[i] == false`, meaning: the counterfactual partner
precursor (`counterfactual_partner_pid[i]`, the opposite-class partner
assigned at pair-regen time) has at least one PSM in some file ≠ k.

This is the **Batch F addition** vs the live (pre-Batch-F) design,
which only required (1) + (2) + (3). We add (4) because the FTR LGBM
trains on the row-doubled frame; a candidate without a counterfactual
would have a missing "fake" training row, breaking the doubling
invariant.

Practical consequence: (4) tightens the candidate set. On 2-file Olsen
it drops candidates from ~28k (live design) to ~3.8k. On 8-file the
effect is much milder (~73k vs ~150k) because `_false`-donor coverage
is ~92% there.

---

## What the candidate set means downstream

- The candidate set defines the **training cohort for the FTR LGBM**:
  doubled rows (top half = real-MBR row with `_true` values; bottom
  half = same row with `_false` MBR values swapped in; non-MBR features
  identical between paired rows).
- The candidate set is also the **only set eligible for recovery**.
  Non-candidate rows can never have `:mbr_recovered = true`.
- Recovery within the candidate set is gated by the FTR q-value
  ≤ α (default 0.01), computed via `get_qvalues!` on the doubled
  frame's OOF LGBM scores.
- Non-candidate rows fall through to the regular `:qval ≤ 0.01` filter
  downstream; their fate is unchanged by the MBR / FTR pipeline.

## Candidates can be targets OR decoys

`candidate_mask` doesn't filter on target/decoy. A decoy whose own
evidence is weak (q > 0.01) but happens to have a high-scoring
same-precursor donor in another file will also be a candidate. The FTR
LGBM trains on both, and recovered decoys contribute to the global
q-value's empirical FDR estimate, which keeps downstream FDR honest.

On the 8-file Olsen entr-library run, of the 56,989 candidates → 54,983
recovered; the global EFDR stayed within 0.0005 of the decoy FDR (see
`BATCH_F_PLAN.md` → "Entrapment FDR validation").
