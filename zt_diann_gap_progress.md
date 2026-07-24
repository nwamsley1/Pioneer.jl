# ZT scanning-DIA vs DIA-NN gap — progress & state (2026-07-24)

Branch: `feat/zt-batched-collapse`. Worktree: `/Users/n.t.wamsley/Projects/pio-sciex-zt`.
Test data/workspace: `/Users/n.t.wamsley/SciexZT_5min` (A_REP1 = `arrow_a1/...Condition_A_REP1.arrow`,
library `Pioneer_HYE_canon_std.poin`, config `cfg_a1_*.json`).

## Goal
Close the ~14% per-file precursor-ID gap between Pioneer and DIA-NN on Sciex ZenoTOF
ZT *scanning*-DIA (swept ~5 Da quad reconstructed into ~1.02 Da contiguous Q1 bins;
each precursor's fragments diluted across ±k metascan bins).

## Current best result
**wide-emit tol=2 (`PIONEER_ZT_CANDIDACY_TOL=2.0`, with `PIONEER_ZT_METASCAN_K=6`)**
= **30,359 IDs** on A_REP1 at the default bitvec threshold (θ=0.03) — vs baseline
27,234 and DIA-NN direct 30,478 (**99.6% parity**). DIA-NN overlap (peptidoform-aware,
seq+mods+charge): 24,395 shared / 5,964 Pioneer-only / 6,352 DIA-NN-only.
**The recovery is real; the problem is that it is not tractable** (see Tractability).

## Env knobs added (ALL default-off / unchanged behavior)
- `PIONEER_ZT_CANDIDACY_TOL` (Da half-width) — wide-emit candidacy. Widens the
  fragment-index precursor box; `map_any_hit_to_center!` re-anchors any-bin hits to
  the center bin. tol=2.0 best; tol=4.0 saturates (23k, cross-talk collapse).
- `PIONEER_ZT_MERGE_K` (bins) — no-reset candidacy merge (score ±K same-cycle
  neighbors into the same bitvec counter without resetting; bit-identical to peak-merge
  via idempotent `or!`). **Inferior to wide-emit — abandoned** (see below).
- `PIONEER_PROFILE_DECONV=1` — dump a flat deconv self-time profile (heavyweight
  sampler; pathologically slow on the long threaded deconv — prefer coarse timers).
- Step timers (always on, DEBUG1): `map_any_hit_to_center!` / `filter_to_center_bin!`
  and `expand_to_metascans!` now timed ("ZT candidacy step timing").

Reminder: **the whole ZT path is gated behind `PIONEER_ZT_METASCAN_K`** — forget it and
all ZT features (wide-emit, merge, metascan expand) silently no-op to the non-ZT baseline.

## Mechanisms tried — verdicts
| mechanism | A_REP1 IDs | verdict |
|---|---:|---|
| baseline (tight center-bin, k=6) | 27,234 | control |
| lower bitvec θ 0.03→0.01 | 29,861 | works but noisier |
| **wide-emit tol=2** | **30,359** | **WINNER (parity)** |
| wide-emit tol=4 | 23,013 | saturates — cross-talk collapse |
| no-reset merge K=1 (single-cal) | 29,485 | < wide-emit |
| no-reset merge K=1 + merged-LUT-cal | 24,292 | over-prunes (merged decoys inflate noise bins → 15/256 pass) |

Merge track abandoned: merging inflates count_ones for decoys too, so it either ties
wide-emit (single-cal) or over-prunes (merged-cal). All merge knobs remain harmless.

## ★ Key diagnostic — the remaining gap is mostly SCORING, not extraction
Decomposed the 6,352 DIA-NN-only precursors (6,276 in Pioneer's library):
- **22% (1,378): never a candidate** — genuine recall gap (frag-index/bitvec/metascan).
- **78% (4,898): scored but below the FDR bar** — median `trace_prob` 0.285 vs passing
  floor 0.982; **but 98% of Pioneer's best-PSM RT is within 0.1 min of DIA-NN's RT
  (median |ΔRT| = 0.000)** — same peptide, same feature, just under-scored.

**Conclusion: ~78% of the DIA-NN-only gap is scoring-sensitivity (right feature, right
RT, ranked under 1% FDR); only ~22% is recall.** Next lever = ML score / features / FDR
calibration, NOT more aggressive candidacy.

Verified in code: `main_search_psms/1.arrow` `q_value` column is an all-zeros
placeholder (filled later in ScoringSearch) — use `trace_prob`, final q = `precursors_long.qval`.

Library equivalence verified: Pioneer & DIA-NN libs share the same FASTAs (HYE canonical
UniProt 2026-06-26), digestion (Trypsin/P, 1 missed cleavage, len 7–30, charge 1–4),
mods (fixed CAM/C, variable Ox/M). Differences: predicted by different models
(Altimeter vs DIA-NN DL) and Pioneer prec-m/z cap 359–1034. Not a coverage artifact.

## Tractability (wide-emit cost)
Step timing (A_REP1, `search_a1_wprof`): frag_index 126 s, **`map_any_hit_to_center!`
181.7 s(!)**, expand_to_metascans! 21.5 s, deconv ~300 s. Deconv load =
center_candidates(40M) × metascan_bins(13) = 519M entries.
- **`map_any_hit_to_center!` (181 s) is the fattest single step** — serial
  `Dict{Int,Set}` + `(cycle,prec)` dedup over 60M emissions. Prime optimization
  (parallel / bucket arrays), no ID risk.
- MainSearch deconv solver = **PoissonMM** (`solvePoissonMM_fast!`, MainSearch/types.jl),
  NOT Huber (Huber is chromatogram_integration only).

## Open next steps (parked)
1. Optimize `map_any_hit_to_center!` (181 s → parallel/bucket-array).
2. `PIONEER_ZT_METASCAN_K=3` — ±3 expansion (~half deconv load; box ~6 Da closer to
   true 5 Da; recall unchanged).
3. Coarse `@elapsed` split of deconv build (`run_fused!`) vs solve (`post_design_matrix!`).
4. **Higher-leverage: improve scoring sensitivity** to recover the 78% right-RT
   under-scored precursors (better features / model / FDR calibration).
5. Entrapment library to validate the 5,964 Pioneer-only IDs' FDR.

## Commit trail
`138ee1a2b` checkpoint (revert point) → `885764301` wide-emit → `aa3617dc7` bench →
`caaa88186` no-reset merge → `5ba658578` merged-LUT-cal → `39da415fb` step timers +
deconv profiler. Revert experiments: `git reset --hard 138ee1a2b`.

Explainer: `/Users/n.t.wamsley/SciexZT_5min/wide_emit_explained.pdf` (5-step wide-emit).
Analysis scripts (in SciexZT_5min): `overlap_pepform.jl`, `diann_only_score.jl`,
`rt_compare.jl` (DIA-NN RT via Parquet2 in scratchpad/pqenv).
