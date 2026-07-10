# Sciex ZT Chromatogram Features — Part 3: Efficiency Plan

*Branch `feat/sciex-zt-scanning` · draft for review · 2026-07-09*

> **Update 2026-07-09 — A1 tested and REJECTED (reverted).** Dropping
> `frag_corr_best_shape` + its consensus loop gave: precursors 75,870 → 75,403
> (−0.6 %), **PGs 12,398 → 12,198 (−1.6 %)** — a real loss beyond noise — for a
> collapse saving of only **2.3 s/file** (30.9 → 28.6 s), not the predicted ~9.6 s.
> The cost attribution in §2 was **wrong**: the O(56) consensus loop is cheap
> (length-13 Pearsons); the shape-block cost is dominated by the irreducible
> per-meta-PSM **fragment gather** (`Fbuf` accumulation) + the O(8) frag-vs-weight
> correlations, which we can't remove without losing the *valuable* shape features.
> And `frag_corr_best_shape`'s low main-search gain (58) was misleading — it pulls
> real weight in the experiment-wide 2nd-pass, especially for protein groups.
> **Conclusion: Tier A is a dead end. The real efficiency prize is Tier B/C.**
> Sections below are kept as originally written for the record; read §2–§3 Tier A
> through this correction.

## 0. Where Parts 1–2 left us (measured, not asserted)

The chromatogram-feature redesign is **accuracy-positive and now committed** (`cf30cc4bd`
+ ancestors). Verified numbers:

| Metric | Pre-redesign baseline (`2e9602bee`) | Refactored (`f641960d2`) | Δ |
|---|---|---|---|
| Precursors (3-file, K=6, MBR off) | 73,657 | **75,870** | **+3.0%** |
| Protein groups | 12,228 | **12,398** | **+1.4%** |
| Total runtime | 7.02 min | 7.36 min | **+4.8%** |

The feature split is doing real work in *both* models (1-file gains):

* **Main-search LGBM**: `frag_corr_strength` (elution/develop-original, post-collapse) **#3**,
  `n_scans` **#6** (was the useless 13×-cycle artifact), `frag_corr_strength_shape` **#7**.
* **2nd-pass LGBM**: `frag_corr_strength_shape` **#6** — the within-metascan shape axis is
  the discriminator on survivors, where the fit features (`gof`, `hellinger`) saturate.

So the goal of Part 3 is narrow: **reclaim the +4.8% runtime without giving back the
+3.0% IDs.** This is not a general ZT-perf pass; it targets the overhead the redesign
itself introduced, plus a couple of adjacent, low-risk trims.

## 1. Where the time actually goes

Per-stage 3-file timing shows the entire regression is in **Main Search
(225.3 s → 254.2 s, +28.9 s / +12.8 %)**; every other stage is flat (±0.3 s).

Instrumented split of the ZT post-deconv work (1-file REP1, K=6):

```
collapse_to_metascans      30.9 s
add_chromatogram_features!   0.7 s   ← post-collapse, on 3.2M rows
```

Two conclusions:

1. **Moving `add_chromatogram_features!` after the collapse was a big win already** —
   develop's chromatogram math on 3.2M meta-PSMs costs **0.7 s**, versus running it on
   ~35M raw rows pre-collapse. That part is effectively free now.
2. **All new cost lives in the collapse.** The baseline (`2e9602bee`) also collapses
   (~21 s/file — the ZT pipeline always has), but did **not** compute the shape-feature
   block. The +9.6 s/file delta *is* the shape block.

## 2. Root cause: the expensive computation feeds the cheap feature

Inside `collapse_to_metascans` (metascan_collapse.jl:171–205), per meta-PSM:

* **O(8)** `_frag_pcor(fragment_profile, weight_profile)` → produces
  `frag_corr_strength_shape`, `frag_corr_effective_n_shape`, `n_correlated_fragments_shape`,
  `frag_apex_dispersion_shape`, `n_correlated_fragments_bitvec_rank_shape`.
  **These are the valuable features.** ~8 Pearson calls.
* **O(8×7)=56** `_frag_pcor(Fbuf[b], Fbuf[b2])` pairwise **best-consensus** loop
  (lines 190–200) → produces **only** `frag_corr_best_shape`.
  **This is ~87 % of the Pearson calls in the block, and its output is the lowest-gain
  shape feature** (main-search gain 58 of ~5800 top; 2nd-pass 2587 of ~319k top).

Compounding it, `_frag_pcor` (features.jl:1253) recomputes both means and both centered
norms on every call, so each `Fbuf[b]` is re-standardized ~7× inside the consensus loop.

**The single highest-leverage move: kill or cheapen the consensus loop.** It costs the
most and returns the least.

## 3. Efficiency levers (ordered by reward ÷ risk)

### Tier A — reclaim the redesign overhead (the +9.6 s/file)

**A1. Drop `frag_corr_best_shape` (remove the O(56) consensus loop).** *Highest ROI.*
Removes ~87 % of the block's Pearson work. Expected: reclaim the bulk of the +9.6 s/file
(→ roughly runtime-neutral vs baseline) at essentially zero ID cost (it's the weakest
shape feature in both models). *Experiment:* 3-file A/B, drop the feature + loop, confirm
IDs within noise of 75,870 and Main Search back toward ~226 s.
Verify: `frag_corr_best_shape` truly low-value — check it isn't quietly propping up the
2nd-pass (gain 2587 is non-zero). If IDs dip, fall back to A2.

**A2. If we keep consensus, vectorize it.** Precompute per signal fragment `b`: mean `μ_b`,
centered vector, and norm `ν_b`; then `pcor(b,b2) = ⟨x̃_b, x̃_b2⟩ / (ν_b ν_b2)` over the
**upper triangle only** (28 pairs, not 56), reusing the already-standardized vectors.
This is exactly the pattern develop's `_add_fragment_chromatogram_features!`
(features.jl:1441–1465) already uses. Expected: ~4× fewer FLOPs in the consensus path.
Strictly a fallback to A1 — only worth it if A1 shows `frag_corr_best_shape` earns its keep.

**A3. Standardize once for the frag-vs-weight pass too.** Precompute `w`'s mean/norm once
per meta-PSM instead of inside each of the 8 `_frag_pcor(fp, w)` calls. Small, free, do it
alongside A2.

### Tier B — the base collapse cost (~21 s/file, pre-existing)

Not a regression, but it's the largest single Main-Search sub-cost and worth attacking
since we're in here.

**B1. `meta = psms[center_rows, :]` (line 211)** copies **every column** for 3.2M center
rows out of the 35M-row frame — a full wide-row gather. Profile it; if it's a large slice
of the 21 s, build `meta` by selecting only the columns the downstream models consume
(the union of `PRESCORE_FEATURES`, `ADVANCED_FEATURE_SET`, ZT columns, and the id/score
bookkeeping columns) rather than all ~60+.

**B2. Column materialization (lines 214–222)** builds 13 weight columns + 9 profile + 6
shape columns via `Float32[... for m in 1:nm]` comprehensions — 28 fresh allocations of
length 3.2M. Preallocate and fill, or write directly into the `meta` frame. Minor, but
free once B1 is touched.

**B3. Window gather bound.** `lo = max(blk_lo, r - L); hi = min(blk_hi, r + L)` scans up to
`2L+1=27` rows per center to bin ±k neighbors. Since rows are sorted by `(pid, scan)` and
bins are contiguous, the true neighborhood is exactly the ±k rows around `r` when the
metascan is dense; the `±L` half-width is a safety margin for gaps. Measure whether
tightening to `±k` (with a correctness guard for missing bins) is safe and helps. Low
priority — the gather is likely cheap relative to the Pearson work.

### Tier C — structural (the real prize, higher risk; propose, don't commit yet)

The collapse runs **after** deconvolution + spectral scoring have already produced ~35M
rows, then discards ~90 % of them. Level 1 (`2e9602bee`) already skips *scoring* on
non-center bins, but the rows are still **materialized, sorted, and carried** through the
Main-Search DataFrame (this is the 191 GB cumulative alloc in Main Search).

**C1. Collapse earlier / never materialize discarded non-center feature rows.** The only
thing we need from non-center bins is their per-bin `weight` and `frag*_int` (to form the
transmission profile + shape features). If we can compute the meta-scan profile directly
off the deconvolution output — before building the full 35M-row scored DataFrame — we skip
building and sorting ~32M rows we throw away. This is the largest potential win but touches
the deconv→scoring→features boundary and needs its own design note + careful byte-identity
proof. **Out of scope for the first Part-3 pass; flagged for a follow-up.**

### Trim — stop computing model-unused features (small, tidy)

The post-collapse `add_chromatogram_features!` computes `ms1_corr_weight_m0`,
`ms1_corr_m0_m1`, and the `other_window_*` / `n_frags_detected_union` features that show
**0 gain** on single-condition ZT data (and the `ms1_corr_*` pair isn't even in the model
feature lists). Only 0.7 s total, so low priority — but if we touch that pass, gate the
dead computations off for ZT. Revisit if multi-condition ZT data changes their gain.

## 4. Correctness guardrails (every change)

* **Byte-identity for untouched features.** After each change, re-run the 1-file
  `center_only` comparison (`zt_cmp_dumps.jl`) and confirm all non-target columns are
  bit-identical. A perf change must not move a single existing feature value.
* **Non-ZT path frozen.** `_zt_k == 0` must be untouched; confirm a non-ZT dataset
  (any Astral/Exploris config) is byte-identical before/after.
* **ID/PG regression gate.** 3-file A/B after each accepted change; require precursors ≥
  75,870 − (run noise ≈ ±0.3 %) and PGs ≥ 12,398 − noise. Any real dip reverts the change.
* **Bit-identical shape features where math is unchanged.** A1 removes a feature; A2/A3
  must leave `frag_corr_strength_shape` et al. bit-identical (pure reassociation of the
  same sums — verify, since float reassociation can drift; if it drifts, document the
  tolerance and confirm IDs unaffected).

## 5. Recommended experiment order

1. **A1** — drop `frag_corr_best_shape` + consensus loop. 3-file A/B. *(Expected: ~neutral
   runtime vs baseline, IDs ≈ 75,870.)* Decision point: if IDs hold, this may be the whole
   Part-3 win and A2 is unnecessary.
2. **A3 + (A2 if A1 regressed)** — standardize-once / upper-triangle. Byte-identity check on
   shape features, 3-file A/B.
3. **B1 profile** — instrument `psms[center_rows, :]` and column materialization; if either
   is a meaningful slice of the 21 s, apply column-narrowing (B1) + preallocation (B2).
4. **Trim** — gate dead post-collapse features for ZT (only if touching that pass anyway).
5. **Write C1 design note** — separate deep-dive; do **not** implement inline here.

## 6. What this plan deliberately does *not* do

* No accuracy/feature-design changes — that's the *next* thread (craft better shape
  features + deep correctness check), kept separate so perf and accuracy don't confound.
* No `K` retuning, no deconv-solver changes, no Level-2 scoring changes.
* No speculative abstraction. Each step is a measured before/after with a revert gate.
