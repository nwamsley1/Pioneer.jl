---
title: "Per-ID Localization FLR: Opposite-Class Target-Decoy Model"
subtitle: "Target-decoy on site-determining evidence; per-ID FLR"
date: "2026-07-07"
geometry: margin=1in
fontsize: 11pt
---

# 1. Goal

Give every `C(c,n)>1` peptidoform a **per-ID localization FLR** (a localization q-value: the estimated
probability its site assignment is wrong), computed by an LGBM model calibrated against positional-
isomer decoys. `C(c,n)<=1` peptidoforms are unambiguous and carry no localization FLR. No single
global target FLR is imposed; the per-ID value is reported and users filter (typical phospho operating
point ~5-10%).

# 2. The two scores (recap)

For a `(precursor_idx, scan_idx)` instance, using that scan's observed spectrum and the SAME
`permute_fragment_mz` the decoy generator uses (full recipe in `site_determining_ion_score.md` and
`consistent_target_decoy_scoring.md`):

- Compare the instance's assigned placement to a **reference class** of positional isomers of its
  sequence (all such **same-composition** library isomers, not just deconvolution candidates and not
  just single-site moves - no lambda dependency).
- **Distinguishing ions are general (competitors may differ at >1 site).** A fragment distinguishes
  isomers A and B iff their **cumulative modification content differs at that cleavage** --
  `(# A-mods at pos <= i) != (# B-mods <= i)`. Then A and B predict **different m/z** for it, offset by
  `(Delta count) * mod_mass`, which can be `+-1x, +-2x, ...` the mod mass (where mods pile up before a
  cleavage). `support(A)` = observed intensity at A's own predicted m/z, `support(B)` at B's;
  `s(A, B) = support(A)/(support(A)+support(B))` on presence/absence weighted by predicted intensity
  (neutral-loss channels excluded). This replaces the "cleavage between two sites a and a'" special
  case (which only held for single-site differences).
- The distinguishing-fragment set is computed **per competitor pair** (it can be non-contiguous and
  multi-magnitude, e.g. `{S1,S7}` vs `{S4,S10}` distinguishes on `[1,4)` and `[7,10)` separately;
  `{S1,S4}` vs `{S7,S10}` gives a `+2x` shift on `[4,7)`). So a distinguishing ion is NOT attributable
  to one adjacent-site "gap" -- match each isomer at its own m/z, don't bin by fixed gaps.
- Aggregate by **MIN over all competitors** (must beat the hardest one; the hardest need not be a
  single-site move -- a 2-site-different isomer whose distinguishing ions are all low-predicted-
  intensity can be the min, so enumerate every same-composition isomer). Cost is O(L) per pair.

The score (rule: **the decoy must model a *mis-localized* target, which faces its truth** - see 2b(i)):

- **Target** `S_vs_targets` = min over the **other** target isomers (drops itself), `n-1` competitors.
  It is scored only against *real* siblings -- never against a decoy -- so its confidence never rides on
  where a random shadow landed.
- **Decoy** `S_vs_targets` = min over **all** target isomers **including its parent**, `n` competitors.
  The parent is (usually) the present/true isomer, so keeping it makes the decoy get killed by the truth
  -- exactly like a mis-localized target -- rather than escape.
- `S_vs_decoys` (optional secondary feature, target instances only) = min over the **decoy** isomers.
  Not required for the model; `S_vs_targets` is the axis the classifier trains and applies on.

Plus supporting scalars per reference: `n_distinguishing_ions`, `distinguishing_capacity` (total
*predicted* distinguishing intensity vs the hardest competitor - the localizability ceiling),
`coverage`, `min_pair_capacity` (predicted distinguishing intensity against the single hardest
competitor; this is how adjacent/near-tied isomers earn their difficulty rather than being flagged).

# 2b. Two competitor-set details (worked examples)

## (i) Keep the parent for decoys -- a decoy must model a *mis-localized* target

The reference for a decoy is not a correctly-localized target but a **mis-localized** one (both are
wrong placements), and a mis-localized target **faces its truth**:

- Correct target `A` (truth A): vs `{B, C}` -- truth is itself, excluded -> wins -> **high**.
- Mis-localized target `A` (truth B): vs `{B, C}` -- truth **B is in the set** -> killed by B -> **low**.
- Decoy `A'` (truth = present isomer = its parent A): to match the mis-localized target it must **face
  its truth** -> **keep the parent**, `A'` vs `{A, B, C}` -> killed by A -> **low**.

So the decoy lands on the mis-localized target, not on a correct one. Dropping the parent would instead
make `A'` mirror a *correct* target (it would escape to high) -- the wrong reference.

Worked example. Peptide `S A S A S A K`, acceptors **S1, S3, S5** -> targets `A=pS1`, `B=pS3`, `C=pS5`;
shadow sites 2/4/6 -> decoys `A'` (mod->pos2), etc. **The scan contains `A` (pS1)**, so observed b-ions
`b1..b6` are all shifted (+80). A `b_i` is shifted iff mod position <= i:

| isomer (mod pos) | b1 | b2 | b3 | b4 | b5 | b6 |
|---|---|---|---|---|---|---|
| A  (1) | + | + | + | + | + | + |
| A' (2) | - | + | + | + | + | + |
| B  (3) | - | - | + | + | + | + |
| C  (5) | - | - | - | - | + | + |

- `A` (target, vs `{B,C}`): min = **1.0** -> correctly localized vs the other acceptors.
- `A'` (decoy, vs `{A,B,C}`, parent kept): `s(A',A)` distinguishes on **`b1`** (`A'` unmod at pos 2,
  `A` shifted at pos 1; observed `b1` shifted -> supports `A`) = 0.0. min = **0.0** -> `A'` rejected,
  landing on top of a mis-localized target (also 0).

`A` and `A'` differ only at the shadow move (site 1 vs 2), so the one ion that rejects `A'` is `b1`,
which lives in the `A`-vs-`A'` comparison -- keeping the parent is what supplies it. And unlike the
target's competitor set, this does NOT add variance to any *target's* score (a target is only ever
scored against real siblings), so a correct target's confidence never depends on where a random shadow
fell.

Two residual points, neither fixed by the competitor set: (1) the decoy faces one *extra* non-truth
competitor a mis-localized target doesn't, but the `min` is dominated by the truth-kill, so this is a
small, at-most-mildly-optimistic bias. (2) The shadow is a **non-acceptor**, so `A'` is killed by an
easy, implausible `b1`, whereas a real wrong-*acceptor* error is killed by ions between two real (often
farther-apart) sites -- so the shadow decoy is *too easy* and can under-estimate the FLR. **The fix is
decoy construction (short-distance-biased shadow sites near real acceptors), not `C(c,n)`-manipulation**
(bumping decoy `C` would strand the `C=2` stratum with no null). Both are measured on the standards
(Section 5).

## (ii) Accumulate evidence per site, not just per competitor pair

When competitors differ at more than one site, a single pooled `s(A,B)` can be carried entirely by one
site's ions while the other is unconfirmed. Example: peptide `S A A S A A S A A R` (acceptors
`S1,S4,S7`), `A={pS1,pS4}` vs `B={pS4,pS7}` (differ at **S1** and **S7**):

- **b-ions** `b1..b6` are +80 heavier for `A` -> they confirm the **S1** end (`A` has `pS1`, `B` does
  not).
- **y-ions** `y4..` are +80 heavier for `B` -> they confirm the **S7** end.

The two ends are probed by *independent* ion series. If we only observed b-ions, a pooled `s(A,B)`
would read high (S1 confirmed) while **S7 is entirely unconfirmed** -- a false-confident call.

Fix: accumulate **per assigned site** `s`. For each modified site, gather the ions whose m/z depends on
`s`'s occupancy (moving the mod off `s` changes them) and compute `s_site(s)` = "mod at `s`" vs "mod
moved off `s`". Peptidoform confidence = **min over its sites** (weakest link); the per-site vector is
the informative output (and matches per-site-probability reporting conventions).

- A **multi-site-different competitor contributes to multiple sites** (its b-ions to `S1`, its y-ions
  to `S7`).
- **Joint-only ions** -- whose m/z changes only under a *combined* move (e.g. a `+2x` fragment where
  two mods pile up before the cleavage) -- confirm the *arrangement* but neither site individually; a
  site confirmed only by joint ions is flagged **not independently localized**.

New features: `s_site` per site -> `min_site_support` (= the peptidoform localization),
`n_sites_independently_confirmed`, `n_sites_joint_only`.

# 3. The model: target-decoy classifier on `S_vs_targets`

Every instance is scored `S_vs_targets` (2b(i)): a **target** against its `n-1` real siblings, a
**decoy** against **all** targets including its parent. The parent (usually the present/true isomer)
kills the decoy, so a decoy lands exactly where a **mis-localized target** does -- which is precisely
what it must model.

- **Positives:** all **target** positional isomers, `S_vs_targets` (vs siblings).
- **Negatives:** all **decoy** positional isomers, `S_vs_targets` (vs all targets, parent kept).
- **Apply:** score each target's `S_vs_targets` and read off the model probability.

A correctly-localized target beats its siblings -> high; a mis-localized target loses to the truly-
present sibling -> low, landing in the decoy region; the decoys calibrate that region. The apply target
and the mis-localized target it might be are both scored against real siblings, and the decoy (killed
by its truth) sits on the mis-localized target -- so the null is the right one. (The decoy's one extra
non-truth competitor is a small, at-most-mildly-optimistic bias; 2b(i).)

**Feature blocks.** Prior (P: `C(c,n)`, `log C`, site counts, distances, terminal flags) and
competition (W: `iso_group_size_at_scan` = number of competitors in the scan, `weight`, `weight_ratio`,
`iso_weight_fraction`, `mbr`) are intrinsic to the (instance, scan). The **S block** is the
`S_vs_targets` evidence, now **per-site** (2b(ii)): `min_site_support`, `s_site` distribution,
`n_sites_independently_confirmed`, `n_sites_joint_only`, plus `n_distinguishing_ions`,
`distinguishing_capacity`, `coverage`, `min_pair_capacity`.

**One tradeoff to note (label noise).** Because some real (yeast) targets are themselves mis-localized
and we don't know which, a small fraction of the positive labels are "wrong-but-labeled-target." This
is the standard target-decoy situation (a few mislocalized positives; the decoys still define the null)
and is tolerable when correct localizations dominate; we monitor it via the standards. (An earlier
draft avoided this by scoring targets against decoys at train time -- the "opposite-class" variant --
but it does not model the mis-localized-target reference cleanly, so we drop it in favor of this
target-decoy classifier.)

# 4. Per-ID localization FLR

The model outputs a score per instance. Calibrate to a per-ID FLR by target-decoy competition on the
score:

- Score all instances (targets and decoys alike) on `S_vs_targets`; the decoys are the null.
- Sort by score desc; at score `s`, `FLR(s) = D(>=s) / R(>=s)` (`C>1` only).
- Each real ID's **localization FLR = min FLR over cutoffs at or below its score** (monotone q-value).
- **MBR instances are included** on both sides (their decoys are transferred too) -> held to the same
  standard, no special-casing. Report with-MBR and MBR-excluded columns.
- `C<=1`: localization FLR = N/A (unambiguous); still subject to the 1% ID-FDR.

# 5. Validation on the standards (ground truth, no construction)

The standards are `C>1` peptidoforms in the library with **known true sites**; they get a reported
localization and a per-ID FLR like everything else. So:

1. **Correctness label** per detected standard: is the reported isomer's site == the true site
   (`standards_truth`)? This yields the standards' real **acceptor-site** error rate directly - exactly
   the error mode we care about.
2. **Calibration curve:** bin detected standards by predicted per-ID FLR; the empirical error rate in
   each bin should equal the bin's predicted FLR. Calibrated -> the shadow-decoy-based FLR is accurate
   on real localization. Predicted < empirical -> decoys are too easy (under-count) -> the short-
   distance-biased decoy fix (Section 8) is warranted. Predicted > empirical -> conservative (safe),
   consistent with prior runs (decoy ~9.6% vs empirical ~4%).
3. **Ablation:** correctly-localized standard IDs at matched FLR for prior-only vs +W vs +S, to show
   the S block earns its cost.

This is the whole test - no deliberately-mislocalized peptides. The standards' own mis-calls (search
picks the wrong real isomer) are the real wrong-acceptor errors, and the per-ID FLR either predicts
their rate or it doesn't.

# 6. What Phase 1 emits (so we can choose empirically)

The extractor (analysis-only, on the existing combined-search results; reuses `site_determining.jl`)
emits, per `(precursor_idx, scan_idx)` `C>1` instance, the per-site **`S_vs_targets`** block (parent
dropped for decoys) plus P and W. It also emits the optional `S_vs_decoys` block cheaply (pairwise `s`
is already computed) so we can test whether it adds signal, but the model is trained and applied on
`S_vs_targets` alone. Everything is validated by whether the per-ID FLR tracks the standards' empirical
error.

# 7. Phased plan

| Phase | Deliverable | Cost | Gate |
|---|---|---|---|
| 1 | per-site `S_vs_targets` extractor -> per-instance table (P, W, S; reals+decoys) | analysis-only | do S features separate reals from decoys in high-`C`/singleton strata? |
| 2 | LGBM (target-decoy on `S_vs_targets`, parent kept for decoys) + per-ID FLR + standards calibration + ablation | analysis-only | does per-ID FLR track standards' empirical error, and beat prior/W-only? |
| 3 | Wire per-ID localization FLR into reporting (+ tests) | code + tests | holds within `C` strata; MBR handled |
| 4 (opt) | Decoy realism (short-distance-biased shadow sites; decoy iRT perturbation) if under-calibrated | generator + rebuild | predicted FLR tracks empirical more tightly |

# 8. Open questions

1. **Multi-mod** (di-phospho, phospho+ox): *subsumed* by the general cumulative-mod-content definition
   in Section 2 - the competitor set is simply all same-composition isomers, and distinguishing ions
   fall out of the per-pair content difference (no single-move restriction). Remaining question is only
   cost control when `C(c,n)` is large (cap the competitor set by predicted plausibility?).
2. **Separate biological-ambiguity number?** `S_vs_targets` blends "distinguished from real siblings"
   with the FLR calibration. Do we also surface a standalone "confidently a phosphopeptide but ambiguous
   between S4/S7" flag (min vs real siblings only)?
3. **Competitor scope:** all library isomers vs RT/weight-plausibility-gated. All-isomers is cleanest
   (no lambda dependency); gating trades that for less noise from isomers that could never be present.
4. **q-value vs PEP** for the per-ID value (or both).

*Companion: `consistent_target_decoy_scoring.{md,pdf}` (scores + worked example),
`site_determining_ion_score.{md,pdf}` (S recipe), `localization_scoring_plan.{md,pdf}` (broader plan),
`FINDINGS.md` (empirical basis).*
