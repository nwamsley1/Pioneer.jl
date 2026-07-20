# Two-Round Experiment-Wide Scoring — Cross-Run Consistency Features

A second round of PSM scoring whose new inputs are **cross-run consistency features** derived
from the first round: how confidently the same precursor was seen in *other* runs, both
experiment-wide and within a cluster of similar (same-condition) runs. Round 2 combines these
with the base features under a **shadow-decoy regularizer** that forces the model to use them
*jointly* with the base evidence rather than as a standalone "seen-elsewhere" shortcut — which is
what keeps cross-run breadth from turning into false transfer.

Enabled with `ENV["TWO_ROUND"] == "1"`. The learned global-aggregation model is a separate opt-in
(`ENV["GLOBAL_MODEL"] == "1"`). Code: `two_round_scoring.jl` (features + shadows) and
`scoring_interface.jl` (global model).

> **History.** Earlier prototypes explored a single-run `twin_score` (cosine-nearest run),
> `cluster_consensus`, a KNN/Gram cosine neighbor scheme, spectral-agreement gates
> (`shadow_hel`/`gated_hel`), and a `MULTI_FEATURE` grab-bag. All were superseded and removed
> (2026-07) by the GROUP formulation below, which is the single shipped path. The protein-group
> analogue of this machinery was evaluated and **abandoned** — it does not improve protein-group
> *total* IDs (see the offline protein harness); two-round is a **precursor-level** win only.

## 1. Round 1 → OOF scores

The round-1 Pass-1 OOM trainer scores every PSM out-of-fold and writes each file's
`.pass1_sidecar.arrow` with the OOF score `s1 = trace_prob_prepass`. Assumes one row per
`precursor_idx` per file (the MainSearch invariant), so `s1[file, precursor]` is a single value.

## 2. Per-fold run clustering

Clustering is done **independently per CV fold** (no information shared across folds). For a fold's
runs, build a sparse presence matrix (`run → sorted precursor-ids with s1 ≥ GROUP_PRESENCE_S1 =
0.9`) and cluster the runs:

- **`_rsvd_run_embed`** — randomized truncated SVD (Halko) on the sparse presence matrix,
  `O(nnz·ℓ)`, **no R×R Gram** (that build is `O(R²·density)` and dies past a few thousand runs).
- **`_kmeans_small`** on the embedding (seeded, multi-restart → deterministic), then
  **`_split_oversized!`** so no cluster exceeds the target size.
- **`_cluster_runs(run_present, R, k)`** ties these together; `k ≤ 0` → one cluster (all runs).

Cluster-size schedule **`_group_size(R)`**: `R < 9 → 0` (clustering OFF — cluster features mirror
global), else `min(9, R÷2)`. So small experiments degrade gracefully to pure global features.

## 3. GROUP cross-run features (`GROUP_COLS`)

For each `(precursor, run)`, computed in **one `O(N)` per-fold pass** via a per-precursor top-4
`(score, run)` record (`TopK4`) with **leave-one-out** (the precursor's own run is excluded from
its own feature — `_tk_loo_max`, `_tk_loo_top3`):

| feature | over | meaning |
|---|---|---|
| `global_max` | all other runs | max other-run `s1` |
| `global_top3_logodds` | all other runs | Σ positive logits of the top-3 other-run `s1` (`_pos_logit`) |
| `global_n_present` | all other runs | # other runs with `s1 ≥ 0.9` |
| `global_n_confident` | all other runs | # other runs with `s1 ≥ GROUP_CONF_S1 = 0.99` |
| `cluster_max`, `cluster_top3_logodds`, `cluster_n_present`, `cluster_n_confident` | the run's fold-cluster | same four, scoped to the run's cluster |

`top3_logodds` rewards **reproducibility** (confident in 3 runs ≫ in 1) and is positive-only, so
weak/absent runs never subtract. The **cluster** features are condition-local: a precursor absent
from a run's own condition-cluster gets no cluster boost even if some *other*-condition run scores
it — the mechanism that limits cross-condition false transfer. `delta_irt = |irt_obs − ref_irt|`
(`ref_irt` = `irt_obs` of the precursor's highest-`s1` instance) is also written; it uses the best
instance's **iRT only, never its score**.

## 4. Shadow-decoy regularization

`inject_shadow_decoys!` adds, for every real row, one **symmetric shadow** of the opposite class:

```
real target T → shadow_decoy  = nearest-iRT decoy's columns,  but GROUP features grafted from T
real decoy  D → shadow_target = nearest-iRT target's columns, but GROUP features grafted from D
```

Pairing is by nearest `irt_obs` (a **score-independent** property) within the same fold-file; the
shadow inherits its *source's* `precursor_idx`/`cv_fold` (so precursor-keyed folds stay valid) and
carries an `is_shadow` flag. Grafting each GROUP feature from the real row forces its target/decoy
marginal to **1:1**: at every grafted value there is one target-labeled and one decoy-labeled row,
so round-2 LGBM cannot trust a GROUP feature standalone — it must use it *jointly* with the base
features. `delta_irt` is **not** grafted (each shadow keeps its source row's own iRT
self-consistency). `remove_shadow_decoys!` restores the original schema before downstream steps;
`log_shadow_fdr_diagnostics` reports the shadow-included vs real-only FDR as a diagnostic.

This is a **tree** regularizer (the 1:1 marginal is inert for a linear model). It is why the GROUP
features add breadth without leaking, and is the piece that does **not** transfer to the (linear,
few-row-per-entity) protein-group scorer.

## 5. Round 2

```
features2 = [ ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt ]
```

Trained by the same OOM trainer on the shadow-augmented set, **same CV folds** as round 1, scoring
OOF. The round-2 OOF score overwrites each `.pass1_sidecar.arrow`, becoming the experiment-wide
score fed to the global aggregation.

## 6. Global aggregation (optional learned model)

By default the experiment-wide score per precursor is the fixed **top-√N log-odds** of its per-run
scores. With `GLOBAL_MODEL=1`, `build_precursor_global_model_dict` (scoring_interface.jl) instead
fits a small LightGBM on per-precursor aggregate features (n, top-1/2/3, min, mean, std, frac>0.99,
range, top-√N log-odds) with **2-fold precursor-keyed CV** (train one fold → score the other;
clean OOF because each precursor is in one fold), using streaming indexed accumulators (no
per-precursor `Vector`/sort). This is independent of §§2–5 and can be A/B'd on its own.

## 7. CV protocol & leakage safety

The "model" is two model sets (`{M1_f}`, `{M2_f}`) plus the per-fold cluster maps. Scores are OOF
through **all** of them, guaranteed by two invariants:

- **Invariant A — precursor-consistent folds.** CV fold is keyed to `precursor_idx`: every PSM of a
  precursor (across all runs) shares one fold. *Required*, because a PSM's cross-run features look up
  other instances of the **same precursor**, and consistent folds keep those instances in the same
  held-out fold, so their `s1` is OOF-consistent. Random per-PSM folds leak through these lookups.
- **Invariant B — identical folds across rounds.** A PSM's fold never changes between rounds.

Run clustering is computed **per fold** from that fold's OOF presence, so cluster assignments never
carry information across folds either.

## 8. Production finding (EWZ) — filter ordering is decisive

Two-round's benefit depends on the FDR filter ordering. On EWZ (40 files: 20 GO113 human-only +
20 GO114 human+yeast; entrapment false-transfer = yeast IDs in GO113):

- **Global-first** (filter `global_q ≤ 1%` first, then recompute/​filter `ew_q`) — the ew filter is
  inert (survivors already ~0.5% FDR), so the round-2 ew improvement filters nothing → roughly
  neutral.
- **AND** (`global_q ≤ 1%` **AND** pre-global `ew_q ≤ 1%` simultaneously) — the ew filter now bites,
  and **two-round is the ID-recovery lever that only expresses under the AND**, recovering a large
  band of global-passing-but-ew-failing instances at a small false-transfer cost.

**Takeaway:** the AND filter is the false-transfer lever; two-round is the ID-recovery lever that
only expresses under it. Judge by **total** IDs at controlled false transfer (not unique).

## 9. Reproduction

```bash
# Two-round GROUP features + shadow-decoys (+ learned global model), AND filter:
TWO_ROUND=1 GLOBAL_MODEL=1 EWFULL_AND_THRESH=0.01 \
  julia --project=. --threads 10 --gcthreads 5,1 \
  -e 'using Pioneer; SearchDIA("<config>.json")'
# baselines: drop TWO_ROUND for round-1 only; drop GLOBAL_MODEL for fixed top-√N logodds.
```

Metrics from `<results>/precursors_long.arrow`: passing = `target & qval ≤ 0.01 & global_qval ≤
0.01`; false transfer = `species == "YEAST" & file_name` starts with `GO113`. Round-2 feature gains
print at `@user_info` via `log_pass_importance`.
