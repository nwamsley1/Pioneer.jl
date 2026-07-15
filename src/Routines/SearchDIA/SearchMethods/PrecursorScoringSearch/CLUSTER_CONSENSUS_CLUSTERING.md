# Run Clustering for the Cluster-Consensus Feature

Exploratory findings (offline, on existing dumps) for how to group runs into
replicate/condition clusters **fully unsupervised** — the reference set for the
`cluster_consensus` generalization of `twin_score` (see `TWO_ROUND_SCORING.md`).
No experiment labels are ever available in production; the clustering must recover
structure from the data alone.

## Recommended recipe

1. **Input = binary presence/absence, NOT scores.** A precursor is "present" in a run
   if its main-search LGBM score passes a stringent per-run FDR (q < 0.001 / PEP < 0.01).
   Scores are never used as continuous values — they only binarize presence. (Scores
   conflate composition with abundance, which mis-groups; presence captures composition.)
2. **SVD-embed** the (L2-normalized) presence matrix (top ~20 components). The SVD
   denoises — outliers can't dominate, and the group structure lives in the top components.
3. **KMeans + silhouette-on-embedding** for auto-K: sweep K, pick the silhouette peak.
4. **Singletons need no fallback:** an outlier / low-quality run has no close neighbor and
   lands in its own size-1 cluster; define its "twin" as **itself**, so `cluster_consensus`
   = its own score (feature is neutral for it). No "low-quality" threshold, no removal.
5. `cluster_consensus(p)` = mean `s1` over the run's cluster (self if alone), leave-one-out.

## Validation against known experimental designs (ground truth from filenames)

| dataset | design | true K | auto-K | ARI | homog | notes |
|---|---|---:|---:|---:|---:|---|
| EWZ | 2 conditions ×20 | 2 | **2** | **1.00** | 1.00 | bulk composition differs (yeast) — trivial |
| MTAC 3P | Sample A/B ×3 | 2 | **2** | **1.00** | 1.00 | recovers the A/B design |
| APMS | 20 baits ×3 | 20 | **20** | **0.60** | 0.87 | subtle (shared proteome); silhouette peak lands exactly on 20 |
| OlsenExploris 3P | 6 E/H/Y ratios | 6 | 2 | 0.28 | 0.39 | **continuous gradient**, not discrete — coarse split is correct |

Silhouette-on-SVD-embedding peaks at the true K for EWZ / APMS / MTAC-3P (2, 20, 2).

## Key findings

- **Presence + SVD/Jaccard, not scores.** APMS silhouette rose 0.39→0.63 switching score→presence.
- **Over-cluster is safe, under-cluster is dangerous.** Homogeneity (are clusters pure?) is the
  safety metric, not ARI — over-clustering keeps homogeneity ≈1 (splits a group into pure pieces);
  under-clustering merges distinct groups (mixes baits → penalizes real IDs). Bias K upward.
- **Outliers self-isolate — no removal.** The 20-ID APMS run A24R_02 becomes a size-1 cluster on
  its own under SVD+KMeans; "low-quality run" ≡ "singleton" ≡ "twin = self". Removing runs by an
  ID threshold is ill-defined and unnecessary.
- **Absolute Jaccard is NOT comparable across datasets.** APMS same-bait similarity (0.52) is
  *lower* than EWZ across-condition (0.69) — so no fixed distance threshold works; use relative
  methods (SVD/silhouette/density).
- **Method scaling (clustering step, ms):** KMeans/GMM/HDBSCAN scale (linear / N log N); Spectral
  and AffinityProp blow up (O(N²)–O(N³)) beyond ~1–2k runs. The real cost is the shared upfront
  SVD/Gram (O(N²·n_precursors)). For typical N (<~1k) clustering time is negligible; pick on accuracy.
- **Method comparison** (presence matrix): SVD+KMeans+silhouette is the most *consistent* (right K
  on all: 20/2/2). HDBSCAN(Jaccard) is perfect on clean/small (EWZ, 3P) but under-clusters subtle
  APMS (K=14). HDBSCAN(SVD-embed) fixes APMS (19) but fails small/clean. HDBSCAN auto-flags noise
  (≡ singleton). Avoid Euclidean/Hamming on the sparse binary (joint-absences dominate) — use Jaccard.
- **Continuous gradients** (OlsenExploris ratios) give a coarse split (rank ~2), not the fine labels —
  a distinct, benign behavior mode; ARI-vs-fine-labels understates usefulness there.

Plots (in `RegressionTestsLite/DIANN_vs_Pioneer_PEP/`): `svd_kmeans_silhouette.png`,
`presence_jaccard_clustering.png`, `cluster_recovery_viz.png`, `cluster_quality_viz.png`,
`olsen3p_clustering.png`.

## Status / next

Clustering method chosen (SVD+KMeans+silhouette on stringent-presence). **Not yet wired into
`write_two_round_feature_columns!`.** Next: implement `cluster_consensus` (cluster → LOO mean `s1`,
self-twin for singletons) and run the AND-metric on EWZ vs single-twin (+94k baseline).
