# Two-Stage FDR Analysis: What We Tested and What It Means

## What the script does

The script (`scripts/two_stage_fdr.jl`) reads the **FirstPass prescore Arrow files** — the raw, unfiltered per-file LightGBM scores written at `FirstPassSearch.jl:197` before any global aggregation. Each file has columns: `precursor_idx`, `lgbm_prob`, `target`.

It then runs two stages:

### Stage 1: Replicate the current global FDR filter

This is the same logic as `aggregate_prescore_globally!()` in the pipeline:

1. **PEP-calibrate** each file's `lgbm_prob` via isotonic regression → `calibrated_prob = 1 - PEP`
2. For each precursor, collect its calibrated probs across all files it appears in
3. **Log-odds combine** the top-√n files (here √3 = 1, so it just takes the max across files)
4. Compute **global q-values** on these combined scores
5. Mark precursors passing at 1% global FDR → **188,409 target precursors**

### Stage 2: The new experiment-wide FDR

1. Go back to the per-file data
2. **Remove rows** where the precursor did NOT pass Stage 1
3. **Concatenate** remaining rows from all 3 files into one big list
4. Compute **q-values** on this concatenated list (using each row's per-file score)
5. Count how many rows per file pass 1% FDR at this experiment-wide threshold

## Results

```
Stage 1 passing target precursors (unique): 188,409

Stage 2 concatenated: 558,089 rows (552,546 targets, 5,543 decoys)
  → 99.0% of rows are targets

Stage 2 per-file target counts passing q ≤ 0.01:
  lgbm_prob:       187,406 / 177,524 / 186,896  (551,826 total)
  calibrated_prob: 187,410 / 178,238 / 184,264  (549,912 total)

Current pipeline:   187,410 / 178,238 / 186,898  (same as calibrated)
```

Stage 2 removes essentially nothing — at most ~720 observations out of 558K.

## Why Stage 2 has almost no effect

The reason is purely **mechanical**, not biological:

Once you filter to precursors that passed 1% global FDR, you've already removed almost all decoys. The concatenated list is **99% targets, 1% decoys**. When you compute q-values on this list, even the worst-scoring row only needs to beat a handful of decoys. The decoy/target ratio is so skewed that essentially every row passes 1% FDR trivially.

In concrete terms: 5,543 decoys spread across 558,089 rows. You'd need >554 decoys to accumulate above 1 target to even reach 1% FDR at any score threshold. The decoys are so diluted that the experiment-wide q-value filter has no discriminative power.

## What this does NOT tell us about the SecondPass LightGBM

**This script does not evaluate the SecondPass model at all.** Here's why:

The current pipeline after Stage 1 does NOT just count FirstPass observations. It does:

1. Global FDR filter (Stage 1 — what we replicated) → passing precursors
2. **SecondPassSearch** — a completely new fragment-index search using calibrated RT/mass-error/quad models, restricted to passing precursors
3. **ScoringSearch** — trains a new LightGBM with **20 features** (poisson, gof, spectral contrast, etc.) on the SecondPass PSMs
4. Final FDR on the rescored SecondPass PSMs

The 20-feature LightGBM you showed (poisson, total_ions, gof, err_norm, ...) operates on **SecondPass PSMs**, which have richer features from the refined search. Our script never touches those — it only uses FirstPass `lgbm_prob` scores.

So the question "does the SecondPass LightGBM do anything?" is not answered here. To answer that, you'd need to compare:
- (A) SecondPass PSMs scored by the 20-feature LightGBM → final FDR
- (B) SecondPass PSMs scored by something simpler (e.g., just spectral contrast) → final FDR

## What this DOES tell us

**The experiment-wide q-value idea (applying a second FDR filter on per-file scores after global filtering) adds no value when using the same FirstPass scores.** The global filter is already so effective that a second pass of q-value computation on the survivors is redundant.

This makes intuitive sense: if you've already identified which precursors are real at 1% FDR globally, their per-file observations are overwhelmingly correct. There aren't enough surviving decoys to create meaningful FDR pressure.

## Possible next steps

If the goal is to determine whether the SecondPass re-search + 20-feature LightGBM adds value over FirstPass alone, a more informative experiment would be:

1. **Skip SecondPass entirely** — use FirstPass global-passing precursors directly for quantification
2. **Compare final protein/precursor counts** against the full pipeline
3. This would show whether the expensive SecondPass search + rescoring step improves identification beyond what FirstPass already achieves

(This is essentially what the `feature/bypass-first-pass` branch is exploring.)
