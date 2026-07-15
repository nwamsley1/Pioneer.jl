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

## Real-world cases: the K=1 floor and abundance-blindness

Two single-species cohorts (searched on develop, human lib) extend the picture to the
homogeneous / subtle end:

- **Pyr (20-patient before/after; tested 6-patient subset):** silhouette peak **0.13** — the
  weakest of any dataset; SVD ≈ rank-1. Patient signal is faint (within-patient Jaccard 0.60 vs
  across 0.54, sep +0.066) and there is **no PRE/POST axis** (sep −0.011: treatment doesn't shift
  *which* precursors are ID'd). Homogeneous tissue → **≈1 cluster**. Broad consensus is the right,
  safe reference here.
- **KEAP1KO vs WT (H292, 3+3):** silhouette **0.044**; within/across-condition Jaccard 0.736/0.727
  (sep +0.009 ≈ zero). **Presence clustering is BLIND to KEAP1 KO** because it's an *abundance*
  change, not composition — the same ~90k precursors are ID'd in both (only 16 KEAP1 precursors
  genuinely differ: present in WT, 0 in KO = 0.012% of the list).

**The K=1 floor (silhouette < ~0.2 → pool to one cluster).** Interpretation of silhouette:
`s = (b−a)/max(a,b)` = how much closer a run is to its own cluster than the nearest other, as a
fraction of local scale. >0.5 strong, 0.25–0.5 modest, **<0.2 essentially one blob**. Datasets on
the SVD-embedding: APMS 0.61, EWZ 0.33, Olsen-3P 0.16, **Pyr 0.13, KEAP1 0.044**. A **0.2 floor**
sits cleanly between EWZ (real 2 → keep) and Pyr/KEAP1 (→ collapse to 1). Asymmetry: over-split
(2 when 1) is cheap (still-homogeneous sub-group); under-merge (1 when 2) is the dangerous error,
so a LOW floor (biases toward splitting) is safer — and EWZ's 0.33 clears it comfortably.

**Why the floor is SAFE for the consensus.** The dangerous merge only bites when compositions
differ — a precursor present in one condition, absent in the other, that a merged consensus would
under-vouch. KEAP1 IS such a case (16 precursors), but (a) 0.012% is below ANY clustering's
resolution, and (b) those 16 are already q≈0 by spectral evidence, so a soft consensus feature
won't drop them. When a composition difference is big enough to actually harm (EWZ: thousands of
yeast precursors, sep +0.11), it's big enough to clear the floor and split. **Rule: cluster by
COMPOSITION (presence), not biology — abundance-only condition differences correctly merge.**
(Caveat: this safety is specific to the presence-based false-transfer purpose; for detecting KO
biology or quant normalization you'd need abundance/weights, and KEAP1 KO/WT would then separate.)

## Status / next

Clustering method chosen (SVD+KMeans+silhouette on stringent-presence). **Not yet wired into
`write_two_round_feature_columns!`.** Next: implement `cluster_consensus` (cluster → LOO mean `s1`,
self-twin for singletons) and run the AND-metric on EWZ vs single-twin (+94k baseline).
