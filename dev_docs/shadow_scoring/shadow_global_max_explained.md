---
title: "Shadow-Decoy + `global_max`: A Two-Round Scoring Variant"
author: "N. Wamsley"
date: "2026-07-17"
geometry: margin=1in
fontsize: 11pt
---

# Motivation

Pioneer's single-pass LGBM (round 1) scores each PSM using per-run
features. On its own it has no way to say "the same peptide was
also confidently identified in another run." A two-round scheme adds
a **cross-run consistency feature** and re-trains a second LGBM.

Two design choices control how aggressive that feature is, and
whether the LGBM is allowed to lean on it as a standalone signal:

1. **What is the cross-run feature?** — a summary of the same
   precursor's round-1 scores in *other* runs.
2. **How much can the LGBM trust it?** — do we regularise so the
   marginal on the feature alone can't discriminate, forcing joint
   use with the base features?

`shadow + global_max` is one specific answer to both.

# Pipeline

```
round-1 LGBM  (per-PSM features)     -> per-run s1 score (OOF via CV folds)
       |
       v
write cross-run column: global_max   -> mbr_score per (precursor, run)
       |
       v
inject symmetric shadow rows         -> double the row count per fold-file
       |
       v
round-2 LGBM  (all features + mbr_score)
       |
       v
compute FDR-A (real-only) and FDR-B (shadow-inclusive)
       |
       v
remove shadows                       -> downstream sees the same schema
```

Both new steps sit between the round-1 and round-2 LGBM invocations
in `_score_precursor_isotope_traces_no_mbr` (see
`two_round_scoring.jl`). They only fire when the environment variable
`TWO_ROUND=1` is set AND `match_between_runs` is `false` in config.

# The `global_max` feature

Let `s1(p, r)` = round-1 out-of-fold LGBM score for precursor $p$ in
run $r$. For each PSM at $(p, r)$, define:

$$
\text{global\_max}(p, r)
   = \max_{r' \ne r}\ s1(p, r')
$$

with $0$ when $p$ doesn't appear in any other run. This is written as
a new column into every per-file fold Arrow. It's a strong "does this
precursor appear anywhere else in the experiment?" signal: **one**
confident hit in any other run yields a high value.

Two design contrasts against the earlier `knn5_mean`:

| aspect | `knn5_mean` | `global_max` |
|---|---|---|
| scope | top-5 most similar runs (cosine on s1 profile) | every other run |
| aggregation | mean | max |
| leakage risk | low if top-5 are all same-condition | higher — 1 cross-condition run is enough |

# Shadow-decoy injection

After `mbr_score` (= `global_max` here) is written, but **before**
round-2 training, we inject shadow rows symmetrically.

Per fold-file:

- For every real target $T$: create a **shadow_decoy** row
  = nearest-iRT decoy's columns, with `mbr_score` **overwritten**
  by $T$'s value. Labeled decoy.
- For every real decoy $D$: create a **shadow_target** row
  = nearest-iRT target's columns, with `mbr_score` overwritten
  by $D$'s value. Labeled target.

An `is_shadow` boolean column flags the injected rows (false on
originals). Row count exactly doubles when both classes are present
in a file.

Each shadow inherits its **source row's** `precursor_idx` and
`cv_fold`, so the precursor-keyed CV invariant (all instances of a
precursor share one fold — required for the round-1 scores to be OOF)
remains intact.

## Why symmetric

At every `mbr_score` value that came from a target $T$, we now have:

- 1 target-labeled row (the real $T$)
- 1 decoy-labeled row ($T$'s shadow_decoy)

and analogously for every value from a real decoy. **The marginal
on `mbr_score` alone is uninformative** — the target/decoy ratio at
each `mbr_score` value is exactly 1:1.

The round-2 LGBM is therefore forced to use `mbr_score` **jointly**
with the base features: "high `mbr_score` AND target-like base
features → target". A high `mbr_score` on a row with decoy-like base
features (as in the shadow_decoy) doesn't get promoted.

Empirically, the LGBM's total split-gain on `mbr_score` drops from
**~86%** (no shadows) to **~13%** (with shadows) on the FTR-27
benchmark. The base features re-emerge as the dominant signal.

# Round-2 training and FDR

The round-2 trainer is unchanged — it just sees a fold-file whose
row count is roughly $2\times$ what it was, with the correct
target/decoy labels on every row. The `is_shadow` column is present
in the schema but not in the feature list, so the LGBM ignores it.

After round-2 has written the new OOF scores back to the sidecar,
two FDR variants are logged as a diagnostic before shadow rows are
removed.

## FDR-A — real-only

Filter `is_shadow == true` out first. Compute the standard
target/decoy q-value on the real subset:

$$
\text{FDR}_A(s) = \frac{\#\{\text{real decoys with score} \geq s\}}
                         {\#\{\text{real targets with score} \geq s\}}
$$

Report the maximum count of real targets seen at any $s$ where
$\text{FDR}_A(s) \leq q^*$ (default $q^* = 0.01$).

## FDR-B — shadow-inclusive

Do the ratio on **all** rows including shadows on both sides:

$$
\text{FDR}_B(s) = \frac{\#\{\text{all decoys with score} \geq s\}}
                         {\#\{\text{all targets with score} \geq s\}}
$$

Then count how many **real** targets sit above the score where
$\text{FDR}_B \leq q^*$.

FDR-B is the literal answer to *"even though we doubled the decoy
count, do more real targets still pass the 1% bar?"* If the round-2
LGBM is well-calibrated on the enlarged set, shadow_target rows
scoring similarly to real targets can compensate for the
shadow_decoys in the ratio, and FDR-B can even exceed FDR-A.

## Downstream cleanup

After the diagnostics fire, `remove_shadow_decoys!` filters out
`is_shadow == true` rows from every fold-file AND from its sidecar
(row order must stay aligned). The schema returns to what it was
before injection. Chromatogram integration, protein inference, and
the output writers see the same data they always did.

# Empirical result (FTR-27, MBR off)

Same 27 files, same library, three-species mixture (PA / SAu / UP
high, 9 replicates each). Reported: correct IDs (target passes q ≤
0.01 with expected-species membership) and cross-species false
transfers (target passes but with a wrong-species assignment).

| variant | correct | wrong (cross-species) | LGBM gain on `mbr_score` |
|---|---:|---:|---:|
| develop (no round-2) | 136,878 | 881 | — |
| shadow + `knn5_mean` | 153,393 | 1,192 | 12.5% |
| **shadow + `global_max`** | **159,644** | **1,823** | **13.0%** |
| `k=5 mean` (no shadow) | 168,776 | 1,283 | 85.6% |

`shadow + global_max` recovers **+16.6%** correct IDs vs the develop
baseline while keeping the LGBM's reliance on the cross-run feature
below 15%. It out-recovers `shadow + knn5_mean` (+4% more correct)
but at ~53% higher cross-species FTR — global_max is a stronger
transfer signal AND a stronger leakage signal.

On this dataset the unregularised `k=5 mean` still wins the raw
count contest. The hypothesis under active test is that
`shadow + global_max` will pull ahead on datasets where the top-K
nearest neighbors bleed cross-condition — the regularisation
protects against precisely the leakage that `global_max` invites.

# Config knobs (all in `two_round_scoring.jl`)

```julia
const KNN_SCOPE      = :all_runs      # :nearest_k | :all_runs
const KNN_AGG        = :max           # :mean | :median | :quantile | :max
const KNN_K          = 5              # only when KNN_SCOPE === :nearest_k
const KNN_Q          = 0.75f0         # only when KNN_AGG === :quantile
const SHADOW_DECOY_MODE = :symmetric  # :none | :symmetric
```

Column name auto-derives from these to keep round-2 log lines
labeled correctly (`:global_max`, `:knn5_mean`, `:knn10_median`,
etc.). Setting `SHADOW_DECOY_MODE = :none` turns off shadow
injection while leaving the cross-run feature intact.
