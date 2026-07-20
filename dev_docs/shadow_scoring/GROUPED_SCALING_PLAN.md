# Grouped local cross-run scoring — scaling plan

**Status: design / not implemented.** Goal: make the cross-run feature (`global_max` /
`global_top3_logodds`) scale to large experiments (hundreds–thousands of runs) *and* stay
condition-local, by aggregating within fixed size-bounded run groups instead of globally or
per-run-kNN.

## 1. Why global doesn't scale, and why per-run kNN is wasteful

Let **R** = #runs, **N** = total PSMs (rows), group size target **s** (~10), K = top-K used by the
aggregate (3).

| scheme | per-row work | build cost | total time | notes |
|---|---|---|---|---|
| **global** (`KNN_SCOPE=:all_runs`, current) | scan all R−1 other runs | — | **O(N·R)** | each row does R−1 dict lookups; dies at large R. Also borrows across conditions (leak). |
| **per-run kNN** (`_cosine_topk`, MULTI cluster_*) | scan this run's K neighbors | Gram O(R²)+O(N) | **O(N·K + R²)** | fixes R→K, but every run has a DIFFERENT neighbor set → the per-precursor aggregate can't be shared; K lookups per row. |
| **fixed groups** (proposed) | 1 group-record lookup + O(K) LOO | partition + O(N) | **O(N + R²_partition)** | a group's per-precursor aggregate is computed ONCE and reused by all ~s runs in it. Removes the K factor. |

The key inefficiency of per-run kNN: two runs in the same neighborhood recompute overlapping
aggregates. A **hard disjoint partition** shares one aggregate per group → the per-precursor
top-K is built once per group, then every member run reads it with a cheap leave-one-out.

**Crucial observation that makes s≈10 nearly free of breadth cost:** `global_top3_logodds` only
uses the *top 3* corroborating runs. A same-condition group of ~10 almost always already contains
the 3 best corroborators, so the top-3 drawn from a 10-run group ≈ the top-3 drawn from the whole
40-run condition. We get near-global breadth at O(N) cost.

## 2. Cost/memory model of the proposed scheme

**Build (per group, streamed):** one pass over the group's rows; for each `(precursor, s1)` update
that group's per-precursor record = the **top-(K+1)** `(value, run_id)` pairs. O(N) total, O(K)
per row. Store top-(K+1)=top-4 so LOO of any single run still leaves a valid top-3.

**Apply:** for each row `(p, run r ∈ g)`, read group g's record for p; if r is among the stored
top-4, drop it and take the top-3 of the remainder, else the top-3 is unchanged; sum positive
logits (`_pos_logit`). O(N), O(K) per row.

**Memory:** process **group-by-group** → peak extra memory = one group's per-precursor records
≈ O(N·s/R) = O(N/G). Scores are already in RAM from Pass A (`file_pid`/`file_s1`), so no extra I/O;
each fold Arrow is still read once (Pass A) and written once (feature column). I/O identical to today.

## 3. The grouping — must be condition-pure, then size-bounded

Safety invariant (validated across datasets, see main memory): **never pool runs from different
conditions in one reference group** (that's what lets a WT run vouch for a KO-only precursor). Group
SIZE doesn't matter for safety; PURITY does; over-splitting a same-condition set is always safe.

**Recommended (Approach A — subdivide natural clusters):**
1. Run the existing robust run-clustering (presence-cosine or s1-profile → SVD embed →
   modularity/kmeans + silhouette 0.2 floor; = the `cluster_max` infra) to get natural
   condition clusters.
2. For any cluster with size > s: subdivide into ⌈size/s⌉ sub-groups of ≤ s. Since all its runs
   are same-condition, ANY subdivision is safe — cheap options: chunk by the first within-cluster
   SVD/PCA coordinate, or k-means with k=⌈size/s⌉.
3. Clusters ≤ s stay whole; singletons → own group (feature inert / fall back to twin).

Rejected: direct balanced k-means with k=⌈R/s⌉ — if k < #conditions a group can straddle
conditions → leak. Cluster-first, subdivide-second avoids that.

## 4. Aggregation & feature

Reuse the winner: **`group_top3_logodds`** = LOO top-3 positive-logit sum within the run's group
(identical math to `_top3_logodds`, restricted to group members). Optionally also `group_max`.
Grafted onto shadow-decoys exactly like `global_max` (keeps the 1:1-marginal guarantee).

## 5. Relationship to existing code

- `_top3_logodds` / `_pos_logit` (this branch) — the within-group aggregate, unchanged math.
- `cluster_max` prior art: `feat/cluster-max-scoring` (commit b1e9fee2) already did a modularity
  **hard-partition + LOO-max** — closest existing implementation; this plan = (a) bound cluster size
  to ~s, (b) swap max → top3-logodds, (c) the streamed group-by-group build for O(N).
- MULTI `cluster_*` uses `_cosine_topk` = per-run top-K = the expensive kind we are replacing.

## 5b. Refinements (user, 2026-07-20) — folded into the design

- **Cluster PER CV FOLD.** The grouping is part of feature computation, so it must respect the same
  precursor-keyed CV discipline as the OOF `s1` it aggregates: build the run grouping separately for
  each fold, from that fold's presence. Validated (bench_v5): per-fold homogeneity = the global
  grouping (EWZ 1.00, APMS 0.30@min6); cross-fold ARI EWZ 0.68 / APMS 0.91 — folds agree on the
  condition-level structure (EWZ's 0.68 is only the arbitrary within-condition sub-split, both stay
  1.00 pure). Proper and free.
- **Minimum group size, not strict equal size.** Don't force equal groups (that hurts — splits
  conditions). Instead floor small groups: merge any group < min-size (~6) into its nearest by
  embedding centroid. HONEST TRADEOFF (bench_v5): helps re-merge FRAGMENTS of large conditions (EWZ
  min6 keeps homog 1.00) but can LOWER purity where small WHOLE conditions dominate (APMS min6 drops
  0.50→0.30 — merging pure 3-run bait triplets mixes baits, adds no real corroborators since absent
  members contribute 0 to top-3). Use a MODEST floor (~6), not large; it's a mild safeguard against
  corroborator-starved fragments, not always a win. Full report: `grouping_report.pdf`.

## 6. Open decisions (need user input before implementing)

1. **Group size `s`** — fixed 10? tunable const? (Bigger = marginally more breadth, more leak risk,
   more compute; smaller = safer/cheaper. Top-3 saturation argues ~10 is near-optimal.)
2. **Replace or augment** — does `group_top3_logodds` REPLACE `global_max`, or sit alongside it?
   For large experiments the global feature is both slow and leaky, so likely replace above some R.
3. **Auto-switch by run count** — use plain global/top3 when R ≤ threshold (e.g. ≤ ~2s), grouped when
   larger? (Small experiments: one group = all runs = today's behavior, no code path divergence.)
4. **Grouping input** — presence (on/off, best for knockouts) vs s1-profile vs abundance-Pearson
   (best for up/down). Presence-cosine is the validated default; large experiments may need the
   cheaper embedding+chunk.
5. **Aggregation** — top3-logodds only, or also group_max / group_mean for the model to weigh?

## 7. Suggested first experiment

Since KEAP1 (6 runs) and EWZ (40 runs) are small, one group ≈ global there — grouping won't change
results, only prove correctness/equivalence. The scaling payoff needs a many-run set (APMS 60 /
Pyr) where global is slow and conditions are many. Plan: implement behind `GROUP_SCORING=1` with
`GROUP_SIZE=10`, verify on EWZ that grouped ≈ global (equivalence at R<2s), then benchmark
wall-clock + IDs + leak on the largest available experiment.
