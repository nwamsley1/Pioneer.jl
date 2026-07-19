# Grouping benchmark — inputs × methods, effectiveness & scaling

Unsupervised fixed-size-k run grouping for the grouped local cross-run feature (GROUPED_SCALING_PLAN.md).
Bench: `RegressionTestsLite/grouping_bench/bench_grouping.py` + `plot_grouping.py` → `grouping_benchmark.png`.
Effectiveness scored on 6 real datasets with condition labels (labels used ONLY to score; production
is label-free). Speed/memory on synthetic presence matrices up to 6,000 runs.

Datasets (conditions × runs, runs/condition): KEAP1 2×6 (3), Pyr 6×12 (2), Olsen3P 6×23 (~4),
EWZ 2×40 (20), APMS 20×60 (3), YEAST_KO 20×60 (3). k=9 and 12.

Metric = **within-group condition homogeneity** (mean fraction of a group in its majority condition),
reported vs a **random-grouping baseline** because the ceiling depends on k vs condition size.

## Findings

**1. Input barely matters — use presence @ q<0.001.** Mean homogeneity lift: presence@1e-3 ≈
score-profile ≈ presence@1e-2 (kmeans 0.144 / 0.143 / 0.138). presence@1e-3 is the sparsest matrix
(fastest SVD, least memory) and the validated on/off signal → the default. Continuous score profile
adds nothing and is denser.

**2. Method ranking (effectiveness): balanced-kmeans > kmeans > rec-bisect > sort-chunk.**
- **balanced-kmeans** (kmeans centers + capacity-k assignment) is best everywhere and hits **perfect
  1.00 homogeneity on EWZ** (2 conditions × 20 runs, k=9 — cleanly splits the 40 into pure groups).
- **kmeans** (n_clusters=round(R/k)) is close behind (EWZ 0.67) and much cheaper.
- **rec-bisect** (recursive median split on top SVD coord) is cheap but weaker.
- **sort-chunk** (1-D SVD sort → chunk) ≈ random — a single coordinate loses the structure. Reject.

**3. The k-vs-condition-size ceiling is fundamental, not a method failure.** Where conditions are
SMALLER than k (APMS/YEAST 3 runs/condition, k=9), NO method can be pure — a group of 9 must contain
≥3 conditions, so homogeneity caps ~0.33. balanced-kmeans reaches that cap (it packs whole 3-run
baits). Pyr is a true null (no structure → all ≈ random, correctly). EWZ (conditions ≥ k) is where
methods separate and balanced-kmeans wins big.
  - **Why this is still safe:** the top-3 aggregation is presence-gated — a bait-specific precursor
    only draws corroboration from runs where it's actually present (its own bait); the other baits in
    a mixed group contribute 0. So a mixed group leaks only if a precursor is FALSELY confident in a
    wrong-condition member — the same rare event as global, but now limited to an ~8-run pool instead
    of all R. **Grouping is strictly safer than global regardless of purity.**

**4. Speed/memory: the SVD embed dominates; clustering is cheap.**
- SVD embed: ~0.1 s @100 runs → 3.5 s @6,000 runs (P=40k); linear in nnz (total IDs). Peak memory
  ~712 MB @6,000 runs — the memory driver. Real P is larger but cost scales with nnz, not P.
- Clustering on the embedding: sort-chunk ~µs, rec-bisect ~10 ms @6k, kmeans ~1 s @6k. **balanced-
  kmeans OOM/NaN at 6k** — my impl built the dense R×ng distance matrix (O(R²/k)). Fixable by chunked
  greedy capacity assignment (assign points to nearest center under capacity, streamed) — then
  balanced scales too. The blowup is an implementation artifact, not intrinsic.

## Recommendation

- **Input:** presence at q < 0.001, per-run L2-normalized, precursors present in ≥2 runs.
- **Embed:** truncated (randomized) SVD, d≈16. This is the cost/memory driver — cap d and it's
  seconds + ~GB at thousands of runs. Tractable.
- **Group:** **kmeans (n_clusters = round(R/k))** on the embedding as the scalable default, then split
  any cluster > ~1.5k by one bisection to enforce the size bound. If we want balanced-kmeans's edge at
  scale, swap the dense assignment for a **chunked greedy capacity assignment** (never materialize the
  full R×ng distance matrix).
- **k = 9** (divisible by 3, ≤ R; top-3 saturates so 9 gives ample candidates). Fall back to one
  group = all runs when R ≤ ~2k (small experiments = today's global behavior, no divergence).
- **Don't chase purity on few-run-condition data** (APMS/YEAST): the presence-gated aggregation makes
  impure groups safe; homogeneity is a conservative proxy, real safety = measured leak in a search.

## Open follow-ups
- Confirm on a real search: grouped(k=9) top3-logodds vs global top3 on EWZ (expect ≈, R<2k) and on
  a large set (APMS/YEAST 60 runs) for the scaling + leak payoff.

---

# v2/v3 efficiency bake-off (autonomous session, plot `grouping_final.png`)

Scripts `bench_v2.py` (methods) + `bench_v3.py` (embedding paths). Tested on EWZ + APMS + synthetic
scaling to 20,000 runs. **This corrects the earlier "R×R Gram is the efficient path" claim — it is
NOT, for large R.**

## Effectiveness (condition homogeneity, k=9)

| method | EWZ (2c/40r) | APMS (20c/60r) |
|---|---|---|
| **kmeans** (round(R/k), split oversized) | **1.00** | **0.50** |
| **spectral** (precomputed affinity) | **1.00** | **0.50** |
| agglom-ward | 1.00 | 0.483 |
| balanced-chunk (strict size-k) | 0.95 | 0.35 |
| rec-bisect | 0.85 | 0.30 |
| random baseline | 0.64 | 0.21 |

- **kmeans and spectral tie for best** — perfect on EWZ, at the ceiling on APMS.
- **Strict size-balancing HURTS** (balanced-chunk 0.35 vs kmeans 0.50 on APMS): forcing exactly-k
  splits natural condition groups. Use `kmeans(round(R/k))` and only split clusters that exceed k.

## The efficient embedding — SVD-on-sparse, NOT the Gram (v3)

All three embedding paths give **identical clustering** (EWZ 1.00, APMS 0.50) — so choose on speed:

| R | (A) TruncatedSVD on sparse | (B/C) R×R Gram + eig |
|---|---|---|
| 1,000 | 1.3 s | 6.5 s |
| 5,000 | 6.3 s | 85 s |
| 10,000 | **12.6 s** | OOM / too slow |
| 20,000 | **24.9 s** | — |

Randomized `TruncatedSVD` on the frequency-filtered sparse presence matrix is **O(nnz·d)** and scales
linearly (25 s @ 20k runs). Building the dense R×R Gram is **O(R²·density)** and dominates — eigsh vs
eigh makes no difference because the *build*, not the eig, is the bottleneck. The Gram is a fine
mental model and fine for R ≤ ~2k, but SVD-on-sparse is the implementation that scales.

## Clustering-step scaling (on the embedding)

kmeans ~3 s @10k, balanced-chunk ~3 s @10k, rec-bisect ~50 ms @10k (trivial). **spectral and
agglom-ward do NOT scale** (O(R²) affinity / O(R³) linkage → OOM past ~3–5k runs) — so despite tying
kmeans on effectiveness, they're out for large experiments. **kmeans is the only method that is both
top-effectiveness AND scalable.**

## FINAL RECOMMENDED RECIPE

1. Input: **presence @ q<0.001**, sparse CSR (= `file_pid` per run).
2. **Frequency-filter columns**: drop precursors in <2 runs and in >90% of runs.
3. **L2-normalize rows** (depth control).
4. **Randomized TruncatedSVD, d≈16, on the sparse matrix** — O(nnz·d); 12.6 s @10k, 24.9 s @20k.
   Do NOT build the R×R Gram.
5. **kmeans(n_clusters = round(R/k))** on the embedding; split any cluster > k. Do NOT force exact
   size-k (strict balance drops purity).
6. **k = 9**. Fall back to one group = all runs when R ≤ ~2k (= today's global behavior).

Total ≈ 30 s at 20,000 runs (embed + cluster), memory dominated by the sparse matrix (~2 GB nnz at
20k). Fully tractable for thousands of runs.
