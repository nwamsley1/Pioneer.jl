---
title: "Phospho Localization Scoring & FLR Control"
subtitle: "Plan v2: an LGBM localization score + per-ID FLR, with site-determining-ion features"
date: "2026-07-07"
geometry: margin=1in
fontsize: 11pt
---

# 0. One-paragraph summary

Train an **LGBM model** to discriminate truly- from falsely-localized peptidoforms, using features
that span three tiers: the **difficulty prior** `C(c,n)`, the existing **competition evidence**
(`iso_weight_fraction`), and a new **site-determining-ion (SDI) spectral evidence** term ("S"
features). Calibrate the model's score against the positional-isomer decoys, target-decoy style, to
assign **each `C(c,n)>1` ID its own localization FLR** (a localization q-value). `C(c,n)<=1` IDs are
unambiguous and carry no localization FLR. Per-`C` stratification remains available as a validation
lens / simpler fallback, but the LGBM (with `C(c,n)` as a feature) is the primary approach. This v2
folds in review comments: LGBM-first, per-ID FLR (no single global target), MBR IDs accepted under the
same FLR standard, a detailed S-feature recipe (correctness over efficiency for now), and an accurate
description of what the decoys currently are.

# 1. Goal & success criteria

**Goal.** Assign every `C>1` peptidoform passing 1% ID-FDR a **localization FLR** (per-ID, like a
q-value): the estimated probability that its site assignment is wrong, calibrated from the decoys.
Users filter on it directly (typical phospho operating point ~5-10% FLR - localization is usually less
stringent than the 1% ID-FDR). No single target is baked in; we *report the per-ID value*.

**Success criteria.**

1. The LGBM localization score, target-decoy-calibrated, gives a monotone **per-ID localization FLR**
   that, at any cutoff, matches the realized decoy-FLR.
2. At a fixed FLR the LGBM recovers **more correctly-localized IDs than the raw `frac`** (and than
   per-`C` `frac` thresholds), i.e. the added features earn their keep — shown by an **ablation**
   (prior-only vs `frac` vs `frac`+`C` vs `frac`+`C`+S).
3. The decoy-FLR is **validated against the standards' empirical FLR** (so far decoy runs conservative:
   ~9.6% decoy vs ~4% empirical at `frac>=0.75`).

# 2. Foundation (what this rests on)

- ~91% of phospho IDs are `C>1` (need a localization FLR); ~9% are `C==1` (unambiguous, exempt).
- `frac` (iso_weight_fraction_at_scan) is a good confidence for **co-eluting** peptidoforms but is
  information-free for **singletons** (`iso_group_size==1` -> `frac=1` trivially).
- `C(c,n)` = product over mod types of `C(sites, present)` is a strong monotone difficulty prior
  (singleton FLR by C: C1 3.2%, C2 6.6%, C3 11.2%, C4-6 20.3%, C7+ 38.0%) and rescues the singleton
  regime `frac` can't see.
- Localization already happens implicitly at two pipeline stages (fragment-index candidate selection,
  deconvolution ridge), both keyed on site-determining ions. The S features make that evidence
  explicit and usable by a model.
- Metric discipline: **per-file localization instances** (never max-frac dedup across files -> spurious
  27%); **MBR accepted but held to the same FLR standard** (its decoys are transferred too, so the D/R
  ratio already prices MBR's effect in).

# 3. The model: LGBM over three feature tiers

One row per (peptidoform, file) localization instance with `C>1`. Label for training = target vs
positional-isomer decoy. Features:

**(P) Difficulty prior - library-only, always defined.**
`C(c,n)`, `log C`, per-type candidate-site count `c`, present-count `n`, min/mean distance between
candidate sites, peptide length, whether the assigned site is terminal (b1/y1/y2 excluded from the
library -> unlocalizable-ish).

**(W) Competition evidence - from the existing search, free.**
`iso_weight_fraction_at_scan`, `iso_group_size_at_scan`, `weight_rank_at_scan`, `weight`, `mbr` flag.

**(S) Site-determining-ion evidence - NEW (Section 4).**
`S_C` (min competitive ratio over alternatives), `S_mean`, `n_sdi_obs` (observed site-determining
ions), `sdi_support_frac`, and the **library-predicted** SDI intensity fraction (a prior for "is this
even localizable" - the sparsity idea, computable without the observed spectrum).

LGBM handles the interactions (e.g. trust `frac` less when `C` is high; lean on `S_C` when
`iso_group_size==1`). Per-`C` stratified thresholds stay as a simpler alternative and as a validation
cross-check (does the model's FLR hold within each stratum?).

# 4. The S (site-determining-ion) features - detailed recipe

Correctness-first; optimize later. Computed per (peptidoform, file) at the instance's best-PSM scan
(`scan_idx`). Reuses `site_determining.jl`, `matchVarMods`/`acceptor_positions`, and the SAME
`permute_fragment_mz` the decoy generator uses.

## 4.1 Inputs

- Peptide `seq` (length `L`), and the assigned modification positions from `structural_mods` (per mod
  type; for a decoy, one assigned position is the shadow site).
- Candidate acceptor sites for the mod type: `acceptor_positions(matchVarMods(seq, var_mods))`. Define
  the comparison candidate set `Cand = (real acceptor sites) U (assigned sites of this peptidoform)`.
  For a real peptidoform `Cand` = real acceptors; for a decoy it also includes the shadow site.
- The peptidoform's **library fragments**: for each fragment, `(base_type in {b,y}, frag_index i,
  frag_charge z, mz, predicted_intensity)` (from `detailed_fragments` via
  `precursor_to_fragment_indices`).
- The **observed MS2 spectrum** at `scan_idx`: peak `(mz, intensity)` arrays from the run's Arrow file.
- Fragment mass tolerance (the search's ppm) for matching.
- `mod_mass` for the mod type (e.g. phospho 79.966331).

## 4.2 Which alternatives to compare

We use **single-move alternatives** (matching how decoys are built): for each assigned site `a` of the
mod type and each *free* alternative site `a' in Cand \ {assigned}`, compare "mod on `a`" vs "mod on
`a'`". `S_C` is the worst (min) over all such pairs - the assigned placement must beat every
alternative. (Multi-mod peptidoforms: iterate the assigned mods; the full `C(c,n)`-way joint
enumeration is a later refinement - single-move covers the dominant ambiguity and is what the decoys
model.)

## 4.3 Identify the site-determining ions for a pair (a, a')

An ion is **site-determining** for `(a,a')` iff its cleavage falls between them (it contains exactly
one of `a`, `a'`), so the mod (+`mod_mass`) shifts its m/z under one placement but not the other. With
`lo=min(a,a')`, `hi=max(a,a')`:

- **b ions:** `b_i` contains residues `1..i`; site-determining iff `lo <= i < hi` (i.e. `i in
  [lo, hi-1]`).
- **y ions:** `y_j` contains residues `L-j+1..L`; site-determining iff `lo <= L-j+1 <= hi-... ` -
  equivalently the complementary cleavages, `j in [L-hi+1, L-lo]`.

There are `~2*(hi-lo)` such ions (a b and a y per interior cleavage), across each available fragment
charge `z`. **Exclude neutral-loss channels** (phospho -98 H3PO4 and NL of mod-bearing fragments):
they are not site-determining and can be isobaric with unmodified ions - including them makes decoys
look artificially supported.

## 4.4 m/z under each placement, and matching

For each site-determining ion (given `base_type`, `i`, `z`):

- m/z under the **assigned** placement `a` is the library-predicted `mz` (the peptidoform is built with
  the mod on `a`).
- m/z under the **alternative** placement `a'` is `permute_fragment_mz(mz, base_type, i, z, L, a, a',
  mod_mass)` - the exact function the decoy generator uses to move a mod's contribution between sites.
- Match each m/z to the observed peaks within tolerance; take the matched peak intensity (0 if no
  match). Call them `I_a(ion)` and `I_{a'}(ion)`.

## 4.5 Support and the competitive ratio

```
support(a  | a') = sum over site-determining ions of I_a(ion)
support(a' | a ) = sum over site-determining ions of I_{a'}(ion)
s(a, a') = support(a|a') / ( support(a|a') + support(a'|a) )        # undefined if denom == 0
```

`s -> 1` when the observed spectrum's discriminating ions support the assigned site `a`; `-> 0` when
they favor the alternative `a'`; `~0.5` when ambiguous. **Undefined (`n_sdi_obs == 0` for this pair):
unlocalizable against `a'` from this spectrum** - flag, don't silently score.

## 4.6 Aggregate to per-instance features

```
S_C          = min over alternatives a' of s(a, a')          # must beat every alternative
S_mean       = mean of s(a, a') over alternatives
n_sdi_obs    = # distinct site-determining ions matched in the observed spectrum (any pair)
sdi_support_frac = (sum I_a over site-determining ions) / (sum intensity of ALL matched frags)
sdi_pred_frac    = (sum predicted intensity of site-determining ions) / (sum predicted intensity all)
```

`sdi_pred_frac` is the **library-only "localizability" prior** (Section 8, sparsity): if the ions that
*would* distinguish the sites are predicted to be near-zero intensity, the peptidoform is inherently
hard regardless of the observed spectrum.

## 4.7 Decoys use the identical recipe

A loc-decoy's assigned site *is* its shadow site; its alternatives are the real acceptors. Its library
fragments already carry the permuted (shadow-placement) m/z. Because the observed spectrum reflects the
*real* peptide, the decoy's site-determining ions (supporting the shadow placement) should fail to
match -> low `S_C`. That gap is exactly what calibrates the FLR. Computing S identically for reals and
decoys is what makes the target-decoy FLR valid.

# 5. Per-ID localization FLR (calibration & reporting)

- Score every `C>1` instance (real and decoy) with the LGBM.
- **Target-decoy q-value:** sort by score desc; at a score `s`, `FLR(s) = D(>=s) / R(>=s)` (decoy vs
  real instance counts, `C>1` only). Each real ID's **localization FLR = min FLR over cutoffs at or
  below its score** (monotone, q-value style). Report it as a per-ID column.
- **MBR IDs are included** in both `R` and (their transferred decoys) `D`, so they are held to the same
  standard automatically - no special-casing.
- **Report with-MBR and MBR-excluded** columns (MBR ~ +2 pts, informative).
- `C==1` IDs: localization FLR = N/A (unambiguous); they still carry the ID-level 1% FDR.
- Cross-check the per-ID FLR against the standards' empirical FLR at matched cutoffs.

# 6. What the decoys currently ARE (accurate description) + realism gaps

Per peptide (grouping its peptidoforms):

1. `real_by = real_sites_by_mod(seq, var_mods)` - the real acceptor positions per mod type.
2. `shadow = shadow_decoy_sites(seq, real_by)` - draws, **per mod type, `k` = (number of real acceptor
   sites) positions UNIFORMLY AT RANDOM from all NON-acceptor residues** (`Random.shuffle(avail)[1:k]`,
   deterministic seed). **NOT biased toward the real sites** - they are random non-acceptor positions
   anywhere in the peptide. (The "short-distance-biased" draw is a *proposed* refinement, not what
   runs today.)
3. For each modified peptidoform with genuine ambiguity (`C(c,n)>1` for some mod type), build ONE
   decoy: move one randomly-chosen mod of an ambiguous type from its real site to a random shadow site
   of that type; other mods stay put.
4. Decoy fragments = parent fragments with the moved site's m/z permuted (`permute_fragment_mz`);
   same sequence, **same precursor m/z** (mass preserved), `is_loc_decoy=true`, fresh `pair_id`.
5. **Decoy predicted iRT == parent iRT, exactly** (copied - Koina can't predict a mod on a non-acceptor
   shadow residue).

**Realism gaps (both make decoys potentially TOO EASY -> conservative/over-estimated FLR):**

- **Shadow sites are non-acceptor residues.** A real mislocalization lands on another *acceptor* (S/T/Y)
  site; a phospho on a random non-STY residue is easier for the search to reject. So the decoy tests a
  slightly different (easier) hypothesis than a true site-swap.
- **Decoy iRT == parent iRT exactly** (confirmed unrealistic). Real positional isomers shift RT; ours
  co-elute perfectly with the parent, which (a) is unphysical and (b) makes MBR over-transfer decoys
  (recovered 29% vs 10% for reals). **Perturbing the decoy iRT** is the concrete fix.

These are Phase-4 levers (Section 9): pursue if the FLR calibration is off, but they don't block
Phases 1-3.

# 7. Validation & ablation

1. **Ablation** (headline): correctly-localized IDs at matched FLR for prior-only, `frac`, `frac`+`C`,
   `frac`+`C`+S. Quantifies each tier's lift and whether S is worth the extraction.
2. Per-`C`-stratum FLR calibration (does the single LGBM FLR hold within strata?).
3. Standards ground-truth cross-check at several cutoffs.
4. Per-dilution-point stability (decoy-FLR was flat ~7-8% across 256x).
5. Coverage-vs-FLR curve vs baselines (global `frac`, winner-take-one).

# 8. Deferred / future (remembered, not now)

- **SDI sparsity as an intrinsic limit.** Some precursors are just hard: if the site-determining ions
  have very low predicted relative intensity in the library spectrum, no method can localize them well.
  Captured now as the feature `sdi_pred_frac`; a later focus is a principled "unlocalizable" call
  (report peptide, not site) based on library predicted intensity.
- **Multi-mod joint alternatives** (full `C(c,n)` enumeration) beyond single-move.
- **Decoy realism fixes** (short-distance-biased shadow sites; decoy iRT perturbation).

# 9. Phased plan

| Phase | Deliverable | Cost | Gate |
|---|---|---|---|
| 1 | S-feature extractor + per-instance feature table (existing results; reals + decoys) | analysis-only | do S features separate reals from decoys in high-`C`/singleton strata? |
| 2 | LGBM training + target-decoy per-ID FLR + **ablation** | analysis-only | does `frac`+`C`+S beat `frac`+`C` at matched FLR? |
| 3 | Wire per-ID localization FLR into reporting (+ tests); stratified cross-check | code + tests | per-ID FLR holds within `C` strata + vs empirical? |
| 4 (opt) | Decoy realism (RT perturb / short-distance shadow) if calibration is off | generator + rebuild | decoy-FLR tracks empirical more tightly? |

# 10. Resolved from review + remaining choices

**Resolved:** LGBM primary with `C(c,n)` as a feature (was Option B); per-ID FLR reported (no single
target; typical 5-10%); MBR accepted under the same FLR standard; decoy realism deferred to Phase 4;
SDI sparsity remembered as `sdi_pred_frac` + future work.

**Remaining for you:**

1. Training scope - train the LGBM **within CV folds on this dataset**, or hold out standards entirely
   as an independent test (cleaner validation, fewer training rows)?
2. Score target - train to predict "decoy vs real" (proxy), or, where we have standards ground truth,
   also fit/validate against **true correct-vs-incorrect** labels?
3. Do we want the per-ID FLR as a **q-value** (monotone, min-over-cutoffs) or a raw **PEP-style**
   local error probability, or both?

*Companion docs: `site_determining_ion_score.{md,pdf}` (SDI worked example), `FINDINGS.md` (empirical
basis), `isomer_decoy_shadow_sites.md` (decoy construction).*
