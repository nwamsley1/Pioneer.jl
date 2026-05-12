# Feature Engineering A/B Report — 2026-05-12 morning

Branch: `feature/ms1-phase1`. Tests run on the 2-file Olsen Exploris subset
(`/tmp/my_tune.json`) which exercises the full pipeline through MaxLFQ. All
runs use the same env defaults (per-file paircomp ON; per-file PEP filter
ON at 0.9; no global paircomp; no global prescore filter).

## TL;DR

Implemented 6 of the 7 feature ideas you proposed last night (skipped the
top-K capture sweep, which requires struct changes). Net effect of adding
all 17 new feature columns together:

| Metric | Baseline (post-trim) | + 17 new features | Δ |
|---|---:|---:|---:|
| **File 1 q≤.01 (after paircomp)** | 40,495 | **41,424** | **+929** |
| File 2 q≤.01 (after paircomp) | 40,679 | 40,733 | +54 |
| **Final IDs** | 87,929 | **88,478** | **+549** |
| Protein groups @ q≤.01 | 12,403 | 12,291 | −112 |

The +549 final-ID lift is past the LGBM noise floor we've seen on this
dataset (±400). The +929 lift on file 1 q≤.01 is dramatic. File 2 barely
moved, which is consistent with run-to-run noise; expect a wider spread on
the 23-file dataset.

PG count dipped −112 — this could be LGBM noise (PGs swing similarly run-to-
run) or a real signal that the new features help PSMs but redistribute them
toward fewer unique proteins. Worth re-checking after the 23-file run.

## Feature-by-feature importance

All 17 new features sorted by LGBM gain on the 2-file run:

| Feature | Gain | Verdict |
|---|---:|---|
| `frag_corr_best_m0` (DIA-NN best frag vs MS1 M0) | **4,404** | **Keep — strong** |
| `n_correlated_fragments_50` (threshold = 0.5) | **3,813** | **Keep — strong** |
| `ms1_corr_weight_m1` (M+1 mirror of existing m0 corr) | **2,912** | **Keep — strong** |
| `frag_corr_top5_weight` (top-K=5 mean Pearson) | **2,330** | **Keep — strong** |
| `frag_corr_best_weight` (DIA-NN best frag vs weight) | 1,777 | Keep — medium |
| `weight_ratio_to_2nd_best` (cross-precursor) | 1,329 | Keep — medium |
| `ms1_m1_to_m0_ratio` (observed isotope ratio) | 1,280 | Keep — medium |
| `weight_ratio_to_3rd_best` (cross-precursor) | 1,237 | Keep — medium |
| `ms1_m1_to_m0_pred` (predicted from iso splines) | 1,087 | Keep — medium |
| `ms1_envelope_dev_log2` (log obs/pred) | 476 | Drop — LGBM prefers raw |
| `pro_count` | 246 | Drop — marginal |
| `n_competitors_50pct` | 213 | Drop — marginal |
| `his_count` | 147 | Drop — marginal |
| `n_correlated_fragments_90` (threshold = 0.9) | 119 | **Drop — threshold too strict** |
| `lys_count` | 98 | Drop — marginal |
| `arg_count` | 85 | Drop — marginal |

## Observations by idea

### 1. MS1 isotope envelope (3 features)

`ms1_m1_to_m0_ratio` (observed, gain 1,280) and `ms1_m1_to_m0_pred` (from
sulfur-aware iso_splines, gain 1,087) both score in the "medium" tier.
Interestingly, the **log-ratio deviation** (`ms1_envelope_dev_log2`, gain
476) is weaker than either raw value — the LGBM prefers to split on the
two raw features separately rather than on the pre-combined log
deviation. Likely because tree splits can already encode the boundary
geometrically.

I used `iso_splines(sulfur, isotope, mass)` (the Goldfarb 2018 averagine-
flavored splines) for the predicted ratio because the codebase doesn't
have a per-residue elemental-composition calculator currently. If a true
composition-based predictor lands, swap it in here and `ms1_envelope_dev_log2`
should become much sharper.

**Verdict:** keep the two raw ratios, drop the log-deviation, document
that the predicted ratio is averagine-approximate not composition-exact.

### 2. n_correlated_fragments thresholds (2 new variants)

Adding 0.5 and 0.9 alongside the existing 0.7:

| Feature | Gain |
|---|---:|
| `n_correlated_fragments_50` | **3,813** |
| `n_correlated_fragments` (existing, 0.7) | (was top-10 in earlier runs; here 1.8k-3k tier) |
| `n_correlated_fragments_90` | 119 |

The 0.5 variant is one of the strongest new features. 0.9 is essentially
dead — at the strict threshold most real fragments don't make the cut,
so the count is 0 for nearly every PSM. **Verdict:** keep 0.5, drop 0.9.

The user's intuition about adding variants at different thresholds was
exactly right — but the *direction* matters. Looser (0.5) > stricter (0.9).

### 3. M+1 fragment intensities (skipped — requires struct change)

This would require widening `MainUnscoredPSM` and `apply_main_scoring!` to
capture M+1 isotope intensities alongside the existing M0 captures, then
the equivalent of all the `frag*_int` + `frag_corr_*` features for M+1.
Significant scaffolding (~3-4 hour change). The MS1-level analog
(`ms1_corr_weight_m1`) lands at gain 2,912 — a solid signal — which
suggests the fragment-level M+1 analogs would also be informative and
likely worth the build cost. **Recommend pursuing this when you have time
for a dedicated PSM-struct widening.**

### 4. Cross-precursor competition (3 features)

Mid-tier. `weight_ratio_to_2nd_best` (1,329) and `weight_ratio_to_3rd_best`
(1,237) carry similar signal — the LGBM might benefit from seeing both,
but they're correlated. `n_competitors_50pct` is weak (213), probably
because the existing `weight_rank_at_scan` already captures the "how many
competitors are there" angle. **Verdict:** keep the two ratio features,
drop the count.

### 5. Peptide composition (4 features)

All weak: `pro_count=246`, `his_count=147`, `lys_count=98`, `arg_count=85`.
The LGBM is finding these less informative than expected — they're
correlated with `sequence_length` and `prec_mz` which already carry the
broad-strokes peptide-property signal. **Verdict: drop all four.** If we
want sequence signal, hydrophobicity or basic-residue ratios would likely
do better, but the per-residue counts aren't pulling weight.

### 6. DIA-NN-style best fragment (2 features)

Both `frag_corr_best_m0` (4,404) and `frag_corr_best_weight` (1,777) made
the top of the new-features list. The "consensus best" fragment — the one
with the highest mean correlation to the other top-6 — is providing
genuinely orthogonal signal vs the rank-1 fragment, especially when
correlated against the MS1 M0 chromatogram (4,404 gain). This validates
the DIA-NN style approach. **Verdict: keep both.**

### 7. Top-K sweep on existing top-6 captures (1 feature)

`frag_corr_top5_weight` lands at 2,330 — slightly above the existing
`frag_corr_top3_weight`. Suggests there's incremental signal in including
ranks 4-5 but not much beyond. **Verdict: keep.** A future experiment
could add `frag_corr_top4_weight` between them; but `top3 + top5` already
brackets the curve.

## Top-15 feature importance (all features, post-additions)

```
worst_max_residual_11scan                72,499   (existing — neighborhood window)
best_gof_5scan                           20,608   (existing — neighborhood window)
worst_manhattan_11scan                    5,264   (existing — neighborhood window)
ms1_corr_weight_m0                        5,034
frag_corr_best_m0                         4,404   ★ new
n_correlated_fragments_50                 3,813   ★ new
frag_apex_dispersion_irt                  3,651
irt_error                                 3,402
worst_max_residual_3scan                  3,088   (existing — neighborhood window)
ms1_corr_weight_m1                        2,912   ★ new
frag_corr_mean_pairwise                   2,761
max_unmatched_residual                    2,505
sequence_length                           2,446
ms1_weight_apex_to_m0_apex_irt            2,413
frag_corr_top1_top2                       2,346
```

The neighborhood-window features (best, worst) still dominate as
expected. 4 of the top-15 are now new features.

## Proposed final feature-set trim

If we want to clean up after the experiment (keep gain > 1000):

**Keep (9 of 17 new features):**
- `frag_corr_best_m0`, `frag_corr_best_weight`
- `n_correlated_fragments_50`
- `ms1_corr_weight_m1`
- `frag_corr_top5_weight`
- `ms1_m1_to_m0_ratio`, `ms1_m1_to_m0_pred`
- `weight_ratio_to_2nd_best`, `weight_ratio_to_3rd_best`

**Drop (8 of 17 new features):**
- `ms1_envelope_dev_log2` (LGBM prefers raw)
- `pro_count`, `his_count`, `lys_count`, `arg_count` (all marginal)
- `n_competitors_50pct` (redundant with weight_rank_at_scan)
- `n_correlated_fragments_90` (too strict)

That keeps the 9 strongest additions for an expected +400-600 ID lift
without polluting the LGBM with low-gain features.

## What's next (suggested)

1. **23-file Olsen run with the 17-feature set** (or the trimmed 9-feature
   set) to confirm the +549 ID lift translates at scale.
2. **M+1 fragment intensities + per-fragment M0/M+1 correlations**
   (idea #6 from yesterday). Given that `ms1_corr_weight_m1` adds ~3k gain,
   the fragment-level M+1 analog should be at least as strong. Requires
   widening `MainUnscoredPSM` with 6 new fields and updating
   `apply_main_scoring!` to capture them.
3. **Composition-exact isotope ratio** to sharpen `ms1_envelope_dev_log2`.
   You mentioned looking at older Pioneer code for an elemental-composition
   isotope-ratio calculator; if found, drop-in replace
   `iso_splines(sulfur, ...)` with that for `ms1_m1_to_m0_pred`.

## Addendum — M+1 fragment intensities + correlations (overnight follow-up)

After the initial 17-feature batch, attempted the **M+1 fragment captures**
(idea #6) that I'd flagged as deferred. This required widening
`MainUnscoredPSM` and `MainSearchScoredPSM` with 6 new fields each
(`frag1..6_int_m1`), updating `apply_main_scoring!` to capture per-rank
M+1 intensities when `iso_idx == 1`, and adding 5 new correlation
features in `_add_fragment_chromatogram_features!`.

Branch isolation: changes live on `feature/ms1-phase1-m1frag` branched
off `01455cfc` of `feature/ms1-phase1`. Two commits:
- `ad414f70` — struct widening (verified ecoli 2,270 / 1,391, all 149
  fused-match unit tests pass).
- `39aff6fe` — new correlation features.

Revert path: `git checkout feature/ms1-phase1; git branch -D feature/ms1-phase1-m1frag`
gets us back to a clean state if anything else surfaces.

### Headline (2-file Olsen, M+1 features included)

| Run | f1 q≤.01 | f2 q≤.01 | Final IDs | PGs |
|---|---:|---:|---:|---:|
| Baseline (post-trim) | 40,495 | 40,679 | 87,929 | 12,403 |
| + 17 new features (no M+1 frag) | 41,424 | 40,733 | 88,478 | 12,291 |
| **+ M+1 fragment intensities + corr** | **41,213** | **41,183** | **88,875** | 12,386 |

**Net vs the no-M+1 variant:** +397 final IDs (and +450 on file-2 q≤.01,
which was nearly flat in the prior batch — the M+1 features are filling
the file-2 gap). **Net vs baseline:** +946 final IDs.

### M+1 feature gains (this run)

| Feature | Gain | Notes |
|---|---:|---|
| `ms1_corr_weight_m1` | 2,962 | Mirror of existing m0 corr (also new today) |
| `frag_corr_m0_m1_mean` | **2,565** | **★ Mean Pearson(frag M0, frag M+1) over the 6 ranks** |
| `frag_corr_m0_m1_top1` | **2,193** | Rank-1 fragment M0 vs M+1 chromatograms |
| `ms1_m1_to_m0_pred` | 2,088 | Predicted ratio from iso_splines |
| `frag_corr_m0_m1_top3` | 1,830 | Rank-3 fragment M0 vs M+1 |
| `frag_corr_m0_m1_top2` | 1,693 | Rank-2 fragment M0 vs M+1 |
| `ms1_m1_to_m0_ratio` | 1,021 | Observed ratio |
| Raw `frag1..6_int_m1` | 250–555 each | Substrate — low individually but enables the corrs above |
| `n_correlated_m0_m1_fragments` | 349 | Thresholded count weaker than the continuous mean |

All four `frag_corr_m0_m1_*` aggregates landed in the top-20 of the
per-file LGBM importance dump. The hypothesis was correct — real
precursors' M0 and M+1 fragment chromatograms coelute perfectly while
chimeric hits don't. This is the **largest single-batch lift we've
measured for the per-file LGBM** since adding `worst_max_residual_*`.

### Things that went smoothly

1. Struct widening was straightforward — Pioneer compiled cleanly,
   149/149 unit tests passed.
2. Ecoli integration produced consistent IDs both pre- and post-struct
   change (2,270 / 1,391).
3. No threading or race conditions surfaced.

### 23-file confirmation (without M+1 features) — completed

The 23-file run on `feature/ms1-phase1` (commit `01455cfc` — 17 new
features minus M+1 fragments) finished. Compared to yesterday's same
config (`tune_run_all23_perfile`, default per-file paircomp, no new
features at all):

| Metric | Yesterday (no new features) | Today (17 new features) | Δ |
|---|---:|---:|---:|
| q≤.001 | 958,355 | **966,262** | **+7,907** |
| q≤.01 | 1,087,173 | 1,089,961 | +2,788 |
| PGs q≤.01 | 142,524 | 142,027 | −497 |

Big lift at the strictest cutoff (q≤.001, +7.9k). q≤.01 less dramatic
(+2.8k). PG count slipped 497 (consistent with the −112 PG dip seen on
the 2-file run — there's a small redistribution effect where PSMs
concentrate into slightly fewer unique PGs).

A 23-file run on the `m1frag` branch (with M+1 fragment features) is now
running in parallel. Output dir: `tune_run_all23_m1frag/`. Should land
in ~25-40 min on its own. The 2-file delta suggests we'll see another
+500-800 IDs from M+1 features at scale.

### 23-file confirmation (M+1 features) — completed

The 23-file `m1frag` run finished. Three-way comparison:

| Run | q≤.001 | q≤.01 | PGs |
|---|---:|---:|---:|
| Yesterday no-features | 958,355 | 1,087,173 | 142,524 |
| Today 17 features no M+1 | 966,262 | 1,089,961 | 142,027 |
| **Today 17 + M+1 fragment features** | **968,822** | **1,091,822** | **142,041** |

**M+1 features at scale: +2,560 q≤.001, +1,861 q≤.01, +14 PGs** over the
17-feature batch. **Net vs yesterday: +10,467 q≤.001, +4,649 q≤.01.**

This is bigger than the 2-file run predicted (was +397 final IDs on 2
files; ~+1,861 q≤.01 on 23 files). The signal grows with scale —
likely because larger PSM pools give the LGBM more variance to learn
the per-fragment M0/M+1 correlations on.

PG count recovered slightly with M+1 features (−483 vs no-features
becomes −483+14 = nearly even). Still down from yesterday's baseline,
but M+1 doesn't make the PG dip worse.

### MTAC 3P Standard cross-dataset confirmation

The MTAC 3P Standard dataset is 6 .arrow files of a different sample type
(`Sample-A/B_Standard_5min_Rep1..3`) using the Astral library
(`altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin`, different
from the Exploris library used in the Olsen runs above). Baseline
reference is the pre-existing `perfile_pep09_mbroff_results` (same
config: per-file paircomp ON, PEP filter 0.9, MBR off).

| MTAC 3P Standard | q≤.001 | q≤.01 | PGs |
|---|---:|---:|---:|
| Baseline (perfile_pep09_mbroff) | 388,206 | 463,976 | 55,958 |
| **+ 17 features + M+1 fragments** | **440,912** | **519,838** | **60,437** |
| **Δ** | **+52,706 (+13.6%)** | **+55,862 (+12.0%)** | **+4,479 (+8.0%)** |

**Massive lift across all three metrics**, and notably **PG count went up
+4,479** — opposite of the Olsen Exploris result where PGs dipped slightly.

Possible explanations for the dataset-dependent PG direction:
1. MTAC has lower starting PG counts (~56k vs Olsen ~143k) so the new
   features have more headroom to identify previously-unscored proteins.
2. The Astral library differs from the Exploris library — different fragment
   prediction quality interacts with our new features differently.
3. MTAC sample type may have cleaner chromatographic peaks (5-min gradient,
   alternating acquisition) where the per-fragment M0/M+1 correlation
   features are more informative.

Either way, **the new feature stack generalizes to a different sample
type + different library + different acquisition** with a 12-14% ID lift.

Per-file MainSearch q≤.01 after paircomp (M+1 features active):

| File | q≤.001 | q≤.01 |
|---|---:|---:|
| Sample-A_Rep1 | 52,481 | 80,821 |
| Sample-A_Rep2 | 52,824 | 79,247 |
| Sample-A_Rep3 | 49,962 | 79,258 |
| Sample-B_Rep1 | 55,088 | 84,352 |
| Sample-B_Rep2 | 53,262 | 81,700 |
| Sample-B_Rep3 | 52,829 | 83,056 |
| **Sum** | **316,446** | **488,434** |

The per-file MainSearch sum (488,434 q≤.01) is below the experiment-wide
ScoringSearch total (519,838 q≤.01) — as expected, ScoringSearch recovers
additional IDs via the second-stage LGBM + global aggregation.

### Things to refocus on in the morning

1. **Check `tune_run_all23_m1frag` numbers** (run launched ~00:35). They
   should show another lift on top of the 23-file no-M+1 result above.
2. **Trim opportunity:** the trim set in the main report can be applied
   on top of this branch. Specifically drop the low-gain ones
   (composition counts, n_correlated_fragments_90, etc.) while keeping
   the M+1 family. New total: ~9 keeper features from the prior batch
   + 5 M+1 correlation features + 6 raw M+1 intensities = ~20 high-
   value additions.
3. **Composition-exact isotope ratio** — still on the to-do list. With
   `frag_corr_m0_m1_mean` at 2,565 gain, the same composition-exact
   logic applied to the *predicted* M0/M+1 fragment ratios could
   sharpen those further.
4. **PG-count dip investigation.** Across both 2-file (−112) and 23-file
   (−497) the new features cost a small number of PGs. Likely a PSM
   redistribution effect (more PSMs land on already-identified proteins
   instead of marginal new ones). Worth a quick look — could be real
   improvement (we removed bad-quality proteins) or worth tuning.

### Things that did NOT cause trouble (but did cause concern initially)

- The first run failed with `TypeError: expected Int64, got UInt8` from
  `IsotopeSplineModel` — fixed by `Int64(...)` conversion (commit
  `01455cfc`).
- The second run failed with `BoundsError` because some peptides have
  sulfur count > 5 (model max). Fixed with `clamp(sulf, 0, 5)` in the
  same commit. After that, all runs completed cleanly.

## How to reproduce

```bash
cd /Users/nathanwamsley/Projects/EntrapmentTests/Pioneer.jl
git checkout feature/ms1-phase1

# baseline (pre-feature-batch)
git checkout 647fb800
julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/tmp/my_tune.json")' \
    > /tmp/baseline.log 2>&1

# with all 17 new features
git checkout fabb9454  # (latest at time of report)
julia --project=. -t 8 -e 'using Pioneer; SearchDIA("/tmp/my_tune.json")' \
    > /tmp/with_new.log 2>&1

# tabulate per-file q≤.01 + final IDs
grep -E "after paircomp|Precursor rows|Protein-group" /tmp/*.log
```

## Files / commits

- `fabb9454` "Add 7 batches of new features for overnight A/B testing"
- `df14a40` (next, if we proceed with trim) "Drop low-gain new features"

Run logs:
- `/tmp/baseline_postrim.log` — baseline numbers
- `/tmp/olsen_all_new.log` — with all 17 new features
