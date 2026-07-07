---
title: "Consistent Target-Decoy Localization Scoring"
subtitle: "Score against the target class, not the deconvolution competition"
date: "2026-07-07"
geometry: margin=1in
fontsize: 11pt
---

# 1. The problem this fixes (point 7 recap)

For a target-decoy localization FLR to be valid, reals and decoys must be scored by the **same**
procedure. But a real peptidoform's natural competitor is a **real sibling** isomer (biological
ambiguity, un-labelable), while a decoy's is a wrong non-acceptor placement - different error modes.
And if "the competitor" is the top-weight candidate, the competitor set is **deconvolution-lambda
dependent**. We want a scoring rule that is (a) symmetric across targets and decoys and (b) independent
of which isomers happened to win weight in the scan.

# 2. The proposal (made precise)

For each `(precursor_idx, scan_idx)` with `C(c,n)>1`, using that scan's observed spectrum:

- The **competitor set is all LIBRARY positional isomers of the sequence** (every target and every
  decoy peptidoform sharing the precursor m/z), **not** the deconvolution candidates. Fixed -> no
  lambda dependency. (We can compute distinguishing-ion presence/absence for any isomer because we know
  its predicted fragment m/z; it need not have been a candidate.)
- **Aggregate by MIN** over competitors: a placement must be distinguished from the *hardest* (best-
  supported) competitor. This is what makes the truly-present isomer dominate (Section 3).
- **Two scores:**
  - `S_vs_targets` = min over the **target** isomers (excluding self) of the competitive ratio
    `s(assigned, competitor)`.
  - `S_vs_decoys` = min over the **decoy** isomers (excluding this instance's own derived decoy) of
    `s(...)`. **Target instances only.**
- A **decoy** instance gets `S_vs_targets` **only** (no decoy-vs-decoy). Its competitor set is *all*
  targets - crucially we do **not** drop its parent target (Section 3 shows why).

So the shared, symmetric axis is `S_vs_targets`: a target is scored vs (targets minus itself), a decoy
vs (all targets). Both are "how well am I distinguished from the target class" - the decoy is evaluated
as an **outsider trying to pass as a target**, never as a deconvolution competitor. That is the
consistency we want.

`s(a, a')` is the presence/absence competitive ratio from the S-feature recipe: over the ions whose
cleavage falls between sites `a` and `a'` (the site-determining ions for that pair), `s = support(a) /
(support(a)+support(a'))`, where `support(x)` sums observed intensity at the ion m/z under placement
`x` (weighted by predicted intensity; neutral-loss channels excluded).

# 3. Worked example

**Peptide** `S A A S A A S A A K` (`L=10`). Phospho acceptors **S1, S4, S7** (`C(3,1)=3` -> three target
isomers). Random non-acceptor **shadow sites 3, 6, 9**. Move-one-random decoys: `pS1->dp3`,
`pS4->dp6`, `pS7->dp9` (each target's phospho moved to a shadow). Isomer group = 3 targets + 3 decoys,
all same precursor m/z.

**b-ion phospho shift** (a `b_i` is +80 iff its mod site `<= i`):

| isomer | site | b-ions shifted (+80) |
|---|---|---|
| pS1 | 1 | b1..b9 |
| pS4 | 4 | b4..b9 |
| pS7 | 7 | b7..b9 |
| dp3 | 3 | b3..b9 |
| dp6 | 6 | b6..b9 |
| dp9 | 9 | b9 |

**The scan actually contains pS4.** Observed b-ions: **b1,b2,b3 unmodified; b4..b9 shifted (+80)**
(y-ions symmetric; omitted for brevity). Distinguishing ions for a pair are the `b_i` with
`min(siteX,siteY) <= i < max(siteX,siteY)`; `support` = observed intensity at each placement's m/z.

## 3a. The true isomer pS4

`S_vs_targets` (competitors pS1, pS7; exclude self pS4):

- vs pS1: differ on b1,b2,b3 (pS1 shifts them, pS4 not). Observed b1-b3 = unmod -> support(pS4)=3,
  support(pS1)=0 -> `s=1.0`.
- vs pS7: differ on b4,b5,b6 (pS4 shifts, pS7 not). Observed b4-b6 = shifted -> support(pS4)=3,
  support(pS7)=0 -> `s=1.0`.
- **min = 1.0** (high, correct).

`S_vs_decoys` (competitors dp3, dp9; exclude own decoy dp6): vs dp3 differ on b3 (obs unmod -> pS4),
vs dp9 differ on b4..b8 (obs shifted -> pS4). Both `s=1.0` -> **min = 1.0**.

## 3b. The decoy dp6 (derived from pS4)

`S_vs_targets` (competitors pS1, pS4, pS7; **no exclusion**):

- **vs pS4** (the present isomer, and dp6's parent): differ on b4,b5 (pS4 shifts, dp6 not). dp6 predicts
  b4,b5 unmodified; observed b4,b5 = **shifted** -> support(dp6)=0, support(pS4)=2 -> `s=0.0`.
- vs pS1: differ on b1..b5. Observed b1-b3 unmod (-> dp6), b4,b5 shifted (-> pS1). support(dp6)=3,
  support(pS1)=2 -> `s=0.6`.
- vs pS7: differ on b6 (dp6 shifts, pS7 not). Observed b6 = shifted -> support(dp6)=1, support(pS7)=0
  -> `s=1.0` (a **coincidental** win: pS4 also shifts b6, which happens to match dp6's prediction).
- **min = 0.0** (low, correct) - dominated by the vs-pS4 comparison.

**This is the crux of the design.** If we took `max`, or if we *excluded the parent* pS4, dp6 would
score `0.6-1.0` (high) off coincidental agreement with absent isomers - a false-confident decoy. **MIN
over all targets, keeping the parent, forces the decoy to beat the truly-present isomer**, which it
cannot. That is exactly the wrong-acceptor error mode we worried the shadow decoys couldn't reach -
here the present real isomer *is* the competitor that kills the decoy.

## 3c. A wrong (non-present) target pS7, if it were a candidate

`S_vs_targets` (competitors pS1, pS4): vs pS4 differ on b4,b5,b6; pS7 predicts unmod, observed shifted
-> support(pS7)=0 -> `s=0.0` -> **min = 0.0** (low). A wrongly-placed *real* peptidoform is filtered
just like a decoy - the desired behavior.

## 3d. Summary

| instance | S_vs_targets | S_vs_decoys | verdict |
|---|---|---|---|
| pS4 (true) | 1.0 | 1.0 | confident |
| pS7 (wrong target) | 0.0 | - | filtered |
| dp6 (decoy) | 0.0 | (n/a) | filtered |

On the shared `S_vs_targets` axis: present-correct targets score high, wrong targets and decoys score
low. The target-decoy FLR (target `S_vs_targets` vs decoy `S_vs_targets`) is therefore valid and
symmetric.

# 4. Why this resolves point 7

- **Symmetric reference class.** Everything is scored "how distinguished am I from the *target* class"
  (min over targets minus self). Targets and decoys use the identical rule -> valid target-decoy FLR.
- **Decoys don't compete in the scoring.** The score is presence/absence of distinguishing ions, not
  deconvolution weight. A decoy is judged as an outsider to the target class, not as a weight rival.
- **No lambda dependency.** Competitor set = all library isomers, fixed regardless of which won weight.
- **Reaches the real error mode.** Because MIN keeps the present isomer in the competitor set, the
  decoy is forced to beat the *actual* real placement - the wrong-acceptor-style comparison the shadow
  decoys otherwise miss.

# 5. Feeding the LGBM + open questions

The clean part: `S_vs_targets` exists for both classes -> it is the **calibration axis** (per-ID
localization FLR = target-decoy q-value over `S_vs_targets`) and a shared model feature. `S_vs_decoys`
exists only for targets.

Questions to settle together:

1. **The asymmetric feature.** `S_vs_decoys` is target-only, so it cannot be a symmetric LGBM feature
   (decoys have no value for it). Options: (a) use it only as a *reported* secondary confidence, not a
   model feature; (b) define a decoy analogue (decoy-vs-*other*-decoys) purely so the column exists -
   but you said we don't want decoy-vs-decoy; (c) drop `S_vs_decoys` from the model and keep only
   `S_vs_targets` + P/W/gap features. Which did you intend?
2. **What the model trains on.** Train LGBM to separate target vs decoy using the *symmetric* features
   (`S_vs_targets`, gap features vs targets, P, W); the FLR then comes from the score's target-decoy
   distribution. Is that the intent, or is `S_vs_decoys` meant to enter somewhere specific?
3. **Reported vs calibrated confidence.** `S_vs_targets` mixes "distinguished from real siblings"
   (biological) and "distinguished from wrong placements" (FLR) into one axis. Do we also want to
   surface a *separate* biological-ambiguity number (min vs *real* siblings only), so a user can see
   "confidently a phosphopeptide, but ambiguous between S4 and S7"?
4. **Competitor scope.** All library isomers of the sequence, or only those within some RT/weight
   plausibility so we don't compare against isomers that could never be here? All-isomers is cleanest
   (no lambda dependency); plausibility-gating trades that for less noise.
5. **Multi-mod.** For di-phospho etc., "competitor isomers" are the `C(c,n)` arrangements and the
   distinguishing ions are the symmetric difference of occupied-site sets; the same min-over-competitors
   rule applies but the enumeration grows - confirm single-move vs full enumeration.

*Companion: `localization_scoring_plan.{md,pdf}` (overall plan), `site_determining_ion_score.{md,pdf}`
(S recipe), `FINDINGS.md` (decoy construction + empirical basis).*
