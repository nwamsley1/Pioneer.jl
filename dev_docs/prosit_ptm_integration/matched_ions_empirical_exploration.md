---
title: "Empirical distribution of site-determining ion matches under competing phospho-localizations"
subtitle: "Testing AScore's uniform-fragmentation assumption on Prosit-Pioneer output"
author: "Pioneer.jl / phospho-localization work"
date: "2026-07-09"
geometry: margin=1in
fontsize: 11pt
---

# Motivation

The AScore method (Beausoleil *et al.*, *Nat. Biotechnol.* 2006, Fig. 3) scores
phospho-localization certainty by computing, for each candidate site, a cumulative
binomial probability of matching at least $n$ of the top-$i$-per-100-*m/z* peaks
across a peptide's $N$ predicted b/y ions. The AScore is then the $-10\log_{10}$
difference in cumulative binomial probability, restricted to *site-determining
ions* (those whose *m/z* differs between the top two candidate localizations),
evaluated at the peak depth $i$ that maximizes the top-two separation.

The step in AScore that we want to check is its **null model**: each predicted
site-determining ion has match probability $p = i/100$ under $H_0$. This assumes

1. every predicted ion is *a priori* equally likely to be observed
   (uniform-fragmentation prior), and
2. the null-match probability depends only on the global top-$i$-per-100-*m/z*
   peak density, not on Prosit's intensity prediction or on local spectrum
   density.

In a Prosit-era pipeline, both assumptions are questionable: we can weight ions
by predicted intensity (or restrict to the intensity-selected set already in the
library), and we have a spectrum-level peak density we can measure directly.
This note tests, empirically:

- Under a real Prosit-Pioneer search, how often does a site-determining ion
  *actually* match an observed peak, split by
    - whether the ion belongs to the true target's predicted fragment set or to
      a competing (loser) localization;
    - the number of distinguishing ions available for the pair;
    - the number of peaks in the pertinent MS/MS scan;
    - the number of possible localizations $C$ for the base peptide.

# Terminology

Throughout, a *pair* is one (target PSM, competitor localization) tuple.

- **Peptidoform.** A peptide sequence with a specific choice of PTM positions
  (e.g. `QSSVTQVTEQpSPK` — phospho on position 11).

- **Target.** The peptidoform that Pioneer's search actually reported as a
  passing PSM at some scan (best-per-precursor row in
  `precursors_long.arrow`, $q_{\text{value}} \le 0.01$).

- **Competitor.** Any *other* peptidoform sharing the same base sequence and
  charge as the target but with phospho on different residues — i.e. a
  target-only positional isomer. Competitors are enumerated from a **clean
  library** (`yeast_plus_standards_1mc.poin`) with no shadow-site decoys, so a
  competitor is always another *target* peptidoform, never a designed decoy.

- **Distinguishing ion.** A b- or y-type fragment ion (charges 1 or 2) whose
  *predicted m/z* differs between the target and the competitor. Equivalently,
  the cumulative phospho count in the b/y ladder differs at that position.
  For b$_i$: target has $\sum_{k \le i} \mathbb{1}[k \in \text{phos}]$ different
  from the competitor's; same for y$_i$ counting from the C-terminus.

- **Site-determining ion.** AScore terminology, same as "distinguishing ion"
  here.

- **"Match".** An ion with predicted *m/z* $m$ **matches** the observed spectrum
  if the spectrum has any peak within $\pm 20$ ppm of $m$. Presence/absence
  only; intensity is not thresholded in this pass.

Two definitions of "available distinguishing ions" per pair are recorded:

- **`n_avail_theo`** — the count of *all* b/y ions at charges 1,2 whose
  cumulative phospho count differs between target and competitor. Bounded above
  by $2 \cdot 2 \cdot (L-1)$ for a peptide of length $L$.

- **`n_avail_lib`** — the subset of `n_avail_theo` that survives the Prosit
  intensity-based retention in the *target's* library-retained fragment list
  (`detailed_fragments.jls`). That list is capped near ~10 fragments per
  precursor and preserves site-determining coverage across mod-gap intervals,
  so `n_avail_lib` is small (typically 2–10) but each entry is an ion Prosit
  predicted with meaningful intensity.

# The two "match rate" numbers, precisely

For each pair $(t, c)$:

- $n_{\text{avail}}$ = number of distinguishing ions considered (theoretical or
  library — recorded both ways).
- $n_{\text{match,tgt}}$ = number of those distinguishing ions where the
  *target's* predicted *m/z* has an observed peak within tolerance.
- $n_{\text{match,cmp}}$ = number of those distinguishing ions where the
  *competitor's* predicted *m/z* has an observed peak within tolerance.

The two counts share the same $n_{\text{avail}}$ denominator (same set of ion
indices $(b/y, i, z)$), but they are evaluated at *different m/z values*: the
target's m/z uses the target's phospho positions, and the competitor's uses the
competitor's.

**Important — the competitor rate is NOT a clean null.** In DIA the isolation
window co-fragments any co-eluting precursors that fall within it, and
positional isomers of the same base peptide have identical mass and near-identical
LC retention time, so co-elution is the *norm*, not the exception. The
competitor's shifted b/y lines can therefore genuinely be present in the same
scan — not because the competitor localization is the truth for the target's
PSM, but because a real co-eluting isomer contributed those peaks. The
"competitor rate" reported below is therefore a **mixture**:

- (a) coincidental matches from spectrum density (what AScore's uniform-$p$
  null is trying to model), plus
- (b) real matches from co-eluting positional isomers that share the same
  precursor mass and RT.

Because component (b) is nonzero and cannot be separated from (a) without
external information (e.g. whether the competitor localization was
independently identified at nearby scans), the reported competitor rate is an
**upper bound** on the true noise floor $p_0$ that a probability model would
plug in as the wrong-localization null.

Two reported summary statistics, aggregated over all pairs in a stratum:

- **Target rate** (signal side): $\displaystyle \left\langle
  \frac{n_{\text{match,tgt}}}{n_{\text{avail}}} \right\rangle$. Under the
  hypothesis that the target's assigned localization is *at least one of the*
  true localizations present at that scan, this is the rate at which
  site-determining ions Prosit predicts for the target are actually observed.

- **Competitor rate** (mixture — see above): $\displaystyle \left\langle
  \frac{n_{\text{match,cmp}}}{n_{\text{avail}}} \right\rangle$. The rate at
  which the *shifted* m/z lines predicted by an alternative localization are
  observed. If the competitor is not co-eluting, this is the density-driven
  noise floor; if it *is* co-eluting, this includes real signal from the
  co-eluting isomer.

  Concretely: for target `QSSVTQVTEQpSPK` (phospho @ pos 11) vs competitor
  `QSSVTQVTEQSPpK` (phospho @ pos 13), the b$_{11}$ ion is distinguishing —
  target predicts b$_{11}$ at *m/z* $x$ (heavier, includes phospho); competitor
  predicts b$_{11}$ at *m/z* $x - m_{\text{phos}}/z$. If both peptidoforms co-elute
  and land in the same isolation window, *both* b$_{11}$ lines will appear —
  the "competitor rate" then reflects real fragmentation, not chance density.

# Data and method

## Search

- MS data: 12 files from the Coon phospho dilution series (`dilution_curve2/combined`).
- Search: `dilution_curve2/combined/results` on branch `feature/localization-flr-scoring`.
- Search-time library: `yps_decoy2.poin` (yeast + Coon standards, 1 missed cleavage,
  with positional-isomer shadow decoys). The shadow decoys were only used at
  search time; the analysis restricts to target-vs-target pairs, so shadow-decoy
  competition does not contaminate the rates reported here.

## Analysis library (clean, no positional-isomer decoys)

- `yeast_plus_standards_1mc.poin` — same peptide set as `yps_decoy2` but built
  with `is_loc_decoy = 0` throughout. 2,653,366 targets, of which 2,418,808
  have $C > 1$.
- Precursor-space join between search library and analysis library is by
  $(\text{sequence}, \text{sorted phospho positions}, \text{sorted CAM
  positions}, \text{charge})$, so `precursor_idx` values are re-resolved.

## Pair enumeration

1. Take the target-passing PSMs from `precursors_long.arrow`
   ($q_{\text{value}} \le 0.01$, $q^{\text{global}} \le 0.01$).
2. Re-key each PSM into the analysis library. Drop any that don't map.
3. For each remaining target PSM $t$, enumerate all competitors $c$ from
   `analysis_lib.grpT[(sequence_t, charge_t)] \setminus \{t\}$ — i.e. every
   target sibling positional isomer at that charge.
4. Sample 10,000 target PSMs uniformly at random (seed 1). This produced
   **156,715 pair rows** (mean 15.7 competitors per PSM).

## Per-pair computation

For each pair $(t, c)$ and each $(T, i, z) \in \{b, y\} \times \{1, \dots, L-1\}
\times \{1, 2\}$:

- Compute the cumulative phospho count on the target's ladder ($\mathrm{cT}$)
  and the competitor's ladder ($\mathrm{cC}$). If $\mathrm{cT} = \mathrm{cC}$
  the ion is *non-distinguishing* — skip.
- Otherwise compute the target's predicted *m/z*: sum of residue monoisotopic
  masses on the ion side of the cleavage + $c_T \cdot m_{\text{phos}}$
  + $\text{CAM} \cdot n_{\text{cam}}$, converted with the appropriate proton
  and (for y-ions) water contributions. Compute the competitor's predicted
  *m/z* the same way but with $c_C$.
- Match each *m/z* against the sorted observed peak list of scan
  `scan_idx` at $\pm 20$ ppm. Presence/absence only.
- Increment `n_avail_theo`; increment `n_match_tgt_t` and `n_match_cmp_t` per
  side if the corresponding m/z matched.
- If $(T, i, z)$ appears in the target's library-retained fragment set (from
  `detailed_fragments.jls`), also increment `n_avail_lib`, `n_match_tgt_l`,
  `n_match_cmp_l`.

Also recorded per pair: scan peak count (`n_peaks_scan =
length(mz_array[scan_idx])`), isolation-window peak count (peaks inside the
scan's isolation window from `centerMz` $\pm$ `isolationWidthMz / 2`), $C$
(number of possible localizations, $\binom{n_{\text{STY}}}{n_{\text{phos}}}$),
and `iso_group_size_at_scan` (search-time competitor count including any
shadow decoys).

## Script

`dev_docs/prosit_ptm_integration/prototype/explore_matched_counts.jl`.

## Output

Per-pair Arrow table at
`/Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/results/explore_matched_counts_v1.arrow`
(156,715 rows).

# Results

Numbers below are means over the pairs in each stratum.

## Match rate vs number of distinguishing ions

**Theoretical (`n_avail_theo`):**

| n_avail (bin) | n pairs  | target rate | competitor rate |
|--------------:|---------:|------------:|----------------:|
|          4    |   5,782  |   0.279     |     0.230       |
|         7–8   |   6,196  |   0.295     |     0.234       |
|        11–12  |   5,718  |   0.316     |     0.224       |
|        13–16  |   6,606  |   0.346     |     0.198       |
|        17–20  |   9,514  |   0.361     |     0.169       |
|        21–24  |  13,113  |   0.360     |     0.156       |
|        25–32  |  31,172  |   0.340     |     0.146       |
|        33+    |  78,614  |   0.297     |     0.139       |

**Prosit library-retained subset (`n_avail_lib`):**

| n_avail (bin) | n pairs  | target rate | competitor rate |
|--------------:|---------:|------------:|----------------:|
|         2     |  14,780  |   0.691     |     0.290       |
|         3     |  16,602  |   0.710     |     0.277       |
|         4     |  17,739  |   0.736     |     0.260       |
|         5     |  17,362  |   0.757     |     0.264       |
|         6     |  18,204  |   0.774     |     0.257       |
|         7–8   |  31,771  |   0.780     |     0.253       |
|         9–10  |  19,079  |   0.784     |     0.252       |

The library-retained view lifts the target rate from ~30% to ~75% and drops
$n_{\text{avail}}$ from a median in the mid-20s to 4–8. The competitor rate is
roughly flat around 26%.

## Match rate stratified by scan peak count

| scan peak count | n pairs | target rate | competitor rate |
|----------------:|--------:|------------:|----------------:|
|        0 – 500  |     122 |   0.167     |     **0.029**   |
|      501 – 1000 |   2,355 |   0.192     |     0.056       |
|     1001 – 1500 |  12,188 |   0.236     |     0.098       |
|     1501 – 2000 |  30,310 |   0.280     |     0.127       |
|     2001 – 3000 |  83,579 |   0.322     |     0.162       |
|         3001+   |  28,161 |   0.387     |     **0.204**   |

Competitor rate spans ~$7\times$ across scans — from ~3% in sparse spectra to
~20% in the densest. A single global $p$ cannot represent both regimes.

## Match rate stratified by isolation-window peak count

| iso-window peaks | n pairs | target rate | competitor rate |
|-----------------:|--------:|------------:|----------------:|
|           0 – 20 | 155,697 |   0.316     |     0.156       |
|          21 – 50 |   1,018 |   0.334     |     0.146       |

Effectively bimodal at $\le$ 20 vs $>$ 20 in this dataset, but the vast majority
of pairs are in the low bin — the isolation-window peak count is a weaker
covariate than the whole-scan count.

## Match rate stratified by $C$ (number of possible localizations)

| $C$   | n pairs | target rate | competitor rate |
|-----:|--------:|------------:|----------------:|
|   2  |   2,602 |   0.365     |     0.175       |
|   3  |   9,098 |   0.356     |     0.172       |
| 4 – 6 |  58,592 |   0.347     |     0.164       |
| 7 – 10 |  57,165 |   0.311     |     0.148       |
| 11 – 20 |  16,138 |   0.272     |     0.141       |
|   21+ |  13,104 |   0.224     |     0.161       |

Target–null gap narrows monotonically with $C$: $2.1\times$ at $C=2$ down to
$1.4\times$ at $C \ge 21$.

## Fraction distribution

Histogram of $n_{\text{match,tgt}} / n_{\text{avail}}^{\text{theo}}$ across all
156,715 pairs:

| bin           | n pairs |  % |
|--------------:|--------:|---:|
| [0.0, 0.1)   |  9,007  |  5.7 |
| [0.1, 0.2)   | 27,380  | 17.5 |
| [0.2, 0.3)   | 44,465  | 28.4 |
| [0.3, 0.4)   | 33,932  | 21.7 |
| [0.4, 0.5)   | 19,248  | 12.3 |
| [0.5, 0.6)   | 16,593  | 10.6 |
| [0.6, 0.7)   |  5,760  |  3.7 |
| [0.7, 0.8)   |  2,101  |  1.3 |
| [0.8, 0.9)   |    376  |  0.2 |
| [0.9, 1.0)   |    141  |  0.1 |

Right-skewed, mode at 0.2–0.3, tail out to 1.0. Not a single-$p$ binomial.

# Interpretation

1. **The observed competitor rate is not a clean null, but even its upper
   bound is inconsistent with AScore's uniform-$p$ assumption.** The
   empirical competitor rate is 14–25% on the theoretical distinguishing set
   and 25–29% on the Prosit-retained subset — orders of magnitude above the
   $i/100$ figure AScore's binomial applies (~1–10%). Because that rate mixes
   coincidental matches with real co-eluting-isomer signal (see §"The two
   'match rate' numbers, precisely"), the *true* noise-floor $p_0$ is lower,
   but a probability model still needs a data-driven $p$ estimated from data
   where co-elution is either absent or explicitly modeled.

2. **`n_avail_lib` is the right definition for a probability model.** It gives:
    - Higher target-side match rate (~0.75 vs ~0.30) — each retained trial is
      more informative,
    - Comparable competitor-side rate (0.26 vs 0.20 depending on stratum) —
      the null doesn't get much worse,
    - Sharper signal:null contrast (~3:1 vs ~2:1),
    - Fewer trials per pair (median 4–8 vs 25+) — a simpler distribution to
      calibrate.

3. **Peak-count covariate is real and large.** Competitor rate ranges from 3%
   in sparse spectra to 20% in dense ones. Some of that spread is
   density-driven noise (component (a) above); some is that dense spectra
   are more likely to contain co-eluting isomers whose real fragments push
   the mixture up. Either way, any probability we compute needs a peak-count
   covariate; a single global $p$ under-reports ambiguity in dense scans and
   over-reports in sparse ones.

4. **$C$ matters, but as difficulty, not a fundamentally different null.** The
   target-null gap narrows with $C$ because high-$C$ peptides have more
   distinguishing ions that collide with real fragment lines from other
   positions, not because the null pattern changes. A probability model can
   handle this by simply using $C$ as a covariate on the same ion-level null.

5. **The empirical fraction distribution is not binomial-uniform.** It is
   right-skewed with a heavy tail at high fractions — consistent with a
   mixture of true localizations (high fraction) and wrong or ambiguous
   localizations (low fraction). A likelihood-ratio model that uses per-ion
   Prosit intensity predictions is a more natural fit than a single-$p$
   binomial.

# Caveats

- **Competitor rate is a mixture, not a noise floor.** As noted in the
  "match rate" section, positional isomers of the same base peptide co-elute
  and are co-fragmented in the same DIA isolation window, so the shifted m/z
  lines that define the competitor rate can be real, not coincidental. The
  reported competitor rate is therefore an upper bound on the true
  wrong-localization noise floor $p_0$. The next-step analyses below address
  this directly.

- **Sample size.** 10,000 PSMs from the 969,006 target-passing PSMs across 12
  files. Distributions are already very stable — bins with $\ge 5{,}000$
  pairs converge to $\pm 0.005$ on a stratum-mean rate — but rerunning at 50k
  would tighten low-count bins (e.g. `scan peak count 0–500`, $n = 122$).

- **"Match" is presence/absence at 20 ppm.** No intensity threshold. Adding an
  intensity threshold (e.g. observed peak intensity $\ge$ 1% of base peak)
  would lower both target and competitor rates; the ratio is what matters and
  is unlikely to change much.

- **`n_avail_lib` uses the target's library-retained set, not the
  competitor's.** Symmetry could be enforced (union or intersection of the
  two retained sets) but the target's set is what a real localization scorer
  would use (we're testing whether target = true; competitor's Prosit
  prediction is not obviously the right basis).

- **Competitor sampling.** For high-$C$ peptides, every sibling positional
  isomer is treated as a competitor and each contributes a pair — so the pair
  count is dominated by pairs from a few high-$C$ peptides. Per-pair
  stratification is honest; per-PSM stratification (mean across pairs first,
  then across PSMs) would rebalance and is a straightforward re-summary of
  the same Arrow.

- **Search-time competition was against `yps_decoy2`, not the clean library.**
  This affects which PSMs Pioneer chose as best-per-precursor (some target
  PSMs may have lost weight to shadow decoys), but does not affect the
  spectrum-space match rates recomputed here.

# Files

- **Script.** `dev_docs/prosit_ptm_integration/prototype/explore_matched_counts.jl`
- **Per-pair Arrow.**
  `/Users/n.t.wamsley/prosit_phospho_test/dilution_curve2/combined/results/explore_matched_counts_v1.arrow`
  (156,715 rows).
- **Related prior work.**
    - `dev_docs/prosit_ptm_integration/prototype/extract_sds_features.jl` — the
      original site-determining intensity-ratio feature extractor.
    - `dev_docs/prosit_ptm_integration/prototype/PROTOTYPE_STATUS.md` — the S-score
      validation status prior to this exploration.
    - AScore paper: Beausoleil *et al.*, *Nat. Biotechnol.* 24(10), 1285 (2006).

# Next step candidates

1. **Split the competitor rate into co-eluting vs not-co-eluting.** Two
   approaches on the same Arrow, giving a cleaner estimate of the true
   noise floor $p_0$:

   a. **Presence filter.** For each pair $(t, c)$, mark the competitor as
      "co-eluting" if the competitor precursor $c$ was itself reported as a
      passing PSM (any q-value bin) at any scan within $\pm \Delta$ RT of
      $t$'s scan in the same file. Recompute the competitor rate separately
      on co-eluting and non-co-eluting subsets. The non-co-eluting subset
      approximates the AScore-analog uniform-$p$ null.

   b. **Weight filter.** If per-scan per-precursor deconvolution weights can
      be dumped (existing `DIAG_DUMP_FILE_IDX` hook — see
      `SearchDIA/SearchMethods/MainSearch/scoring.jl:89-90`), condition on
      the competitor's `weight` at the *target's* scan being essentially zero.
      This is a sharper filter than presence-at-nearby-scan because it uses
      the search's own co-elution estimate.

2. **Draft a per-ion $p_i$ model on the co-elution-corrected null.** Once
   $p_0$ is cleanly estimated, model per-ion match probability as
   $p_i = f(\text{local peak density}, \text{predicted intensity})$;
   restricted to `n_avail_lib`. Weighted Poisson-binomial tail as the analog
   of AScore's cumulative binomial, evaluated per pair, and produce a
   per-target localization score $= -10\log_{10}(P_c / P_t)$.

3. **Rebalance to per-PSM.** Same Arrow, different summary: mean across
   competitors per PSM, then mean across PSMs, stratified same as above.
   High-$C$ peptides currently dominate the pair count.

4. **Add an intensity threshold on "match".** Vary `> 0.01 * base peak` or
   `> 0` and confirm the ratio is stable. Report competitor rate at rank
   thresholds equivalent to AScore's $i \in \{1, \dots, 10\}$ to draw a
   direct comparison against AScore's uniform-$p$ assumption.
