---
title: "Per-Peptidoform Site-Determining-Ion Localization Score"
subtitle: "Design proposal — the singleton regime for multi-peptidoform reporting"
date: "2026-07-06"
geometry: margin=1in
fontsize: 11pt
---

# 1. Why we need this

We want to report **more than one peptidoform per peptide** (multiple positional isomers can be
genuinely co-present), not just one max-weight winner per isomer group. The data already supports it:
each isomer is its own `precursor_idx`, deconvolved per-scan and FDR'd independently, and
`precursors_long` carries each peptidoform's best-PSM `iso_weight_fraction_at_scan`.

But the last experiment (10k-amol, 3 reps, 54,066 target-target vs 40,814 target-decoy peptidoforms
after the `qval` + `global_qval` gate) showed a hard wall:

| population | frac >= 0.75 | frac >= 0.90 | frac >= 0.95 | floor |
|---|---|---|---|---|
| all peptidoforms | 11.5% | 10.7% | 10.8% | **~10.7%, never < 10%** |
| contested (>=2 competing siblings) | 11.3% | 8.2% | 6.6% | ~6% |

**The `iso_weight_fraction` cutoff cannot push per-peptidoform FLR below ~10%.** The cause is
diagnostic: at frac >= 0.75, ~half the confident calls are **singletons** — the peptidoform is alone in
its isomer group at its apex scan (`iso_group_size_at_scan == 1`), so `frac = 1.0` *automatically* and
carries **zero localization information**. Their target-decoy FLR (~11.6%) sets the floor.

A peptidoform is a singleton at its apex exactly when its sibling isomers elute at a **different RT** —
which is the "genuinely co-present isomers" case we want to report. The sibling fraction measures
**co-elution competition**; RT-separated isomers don't compete, so the fraction is blind to them. A
peptidoform decoy that separates in RT with a clean (permuted) spectrum is indistinguishable from a
real second form *by the fraction alone*.

**Conclusion: the fraction is the right confidence for the co-eluting (contested) regime, but the
RT-resolved (singleton) regime needs a different, spectral score. This document proposes that score.**

# 2. Core idea

For a singleton peptidoform, we can't use co-elution competition, but we *can* ask the spectrum
directly: **do the fragment ions that actually distinguish the assigned site from every alternative
site support the assigned placement?** Those ions are the **site-determining ions**. The score is a
per-peptidoform, intrinsic measure — it needs only the observed spectrum at the peptidoform's apex
scan and the library fragments; it does **not** require the competing isomer to be detected.

Crucially, this is exactly the evidence that ID-level FDR does **not** protect. A loc-decoy passes
`qval` mostly on **shared** ions (b/y ions whose cleavage does not span the moved site are identical
between real and decoy). It should fail on the **site-determining** ions, because its mod mass sits on
the wrong fragments. The site-determining score isolates precisely that discriminating evidence.

# 3. Definitions

## 3.1 Candidate sites

For a peptidoform of length $L$ with a variable mod assigned to position $a$, the **candidate
alternatives** are the other acceptor residues (for phospho: the other S/T/Y positions). Localization
means: the mod is on $a$, not on any alternative $a'$.

## 3.2 Site-determining ions

**General definition (competitors may differ at more than one site).** A fragment distinguishes two
same-composition isomers $A$ and $B$ iff their **cumulative modification content differs at that
cleavage**: $(\#\,A\text{-mods at pos}\le i) \ne (\#\,B\text{-mods}\le i)$. Then $A$ and $B$ predict
**different m/z** for it, offset by $(\Delta\text{count})\times \text{mod\_mass}$ — which may be
$\pm 1\times, \pm 2\times, \dots$ the mod mass where mods pile up before a cleavage. Match each isomer at
its **own** predicted m/z; the distinguishing-fragment set is per **competitor pair** and can be
**non-contiguous and multi-magnitude** (e.g. $\{S_1,S_7\}$ vs $\{S_4,S_{10}\}$ distinguishes on
$[1,4)$ and $[7,10)$ separately; $\{S_1,S_4\}$ vs $\{S_7,S_{10}\}$ gives a $+2\times$ shift on $[4,7)$).
So a distinguishing ion is **not** attributable to a single adjacent-site gap, and the competitor set
is **all** same-composition isomers (not just single-site moves).

**Single-site special case (the common one, used in the worked example below).** When $A$ and $B$
differ at exactly one site — the mod is on $a$ vs $a'$ — the general test reduces to: the ion's cleavage
falls **between** $a$ and $a'$ (it contains exactly one of $\{a, a'\}$), so the mod mass ($+79.966$ for
phospho) shifts the m/z for one placement but not the other.

- $b_i$ (prefix, cleavage after residue $i$) contains positions $1..i$. It is site-determining for
  $\{a,a'\}$ iff exactly one of $a, a' \le i$.
- $y_j$ (suffix) contains positions $(L-j+1)..L$. Site-determining iff exactly one of $a, a'$ is in
  that range.

**Worked example.** Peptide `ABCSDESFGK` ($L=10$), phospho candidate sites $S_4$ and $S_7$.

```
position :  1  2  3  4  5  6  7  8  9 10
residue  :  A  B  C  S  D  E  S  F  G  K
                    ^cand a=4      ^cand a'=7
b-ions   :  b1 b2 b3 b4 b5 b6 b7 b8 b9
y-ions   :     y9 y8 y7 y6 y5 y4 y3 y2 y1   (y_j covers L-j+1..L)
```

Ions whose cleavage lies in the gap $[4,7)$ are site-determining:

- $b_4, b_5, b_6$ — contain $S_4$ but not $S_7$ -> carry $+80$ if the site is $S_4$, unmodified if $S_7$.
- $y_4, y_5, y_6$ — contain $S_7$ but not $S_4$ -> carry $+80$ if the site is $S_7$, unmodified if $S_4$.

So the real peptide (site $S_4$) shows $b_{4,5,6}$ at $+80$ and $y_{4,5,6}$ unmodified. A decoy claiming
$S_7$ predicts the opposite and will **not** match the observed peaks on those six ions. The gap of
distance $\delta = 3$ yields $\sim 2\delta = 6$ site-determining ions — matching the retention design
already in `site_determining.jl`.

## 3.3 Phospho-specific caveat: neutral loss

pS/pT spectra are dominated by the **neutral loss of H3PO4 ($-97.977$)** from the precursor and from
mod-bearing fragments. Neutral-loss peaks are *not* site-determining (the loss can occur from any
placement) and, worse, a $b_i^{-98}$ can be isobaric with an unmodified $b_i$. **The scorer must
exclude neutral-loss channels from the site-determining set** (or model them explicitly), or decoys
will look better than they are. This is the single biggest correctness risk (see §9).

# 4. Score formulations

Let $I(\cdot)$ be the observed intensity of a matched peak (0 if unmatched within mass tolerance),
computed at the peptidoform's apex scan.

**(A) SDI count** — number of matched site-determining ions supporting the assigned site.
Simple, but not abundance-normalized and ignores the alternatives' evidence.

**(B) SDI intensity fraction** —
$$ S_B = \frac{\sum_{\text{ion} \in \text{SDI, supports } a} I(\text{ion})}{\sum_{\text{ion} \in \text{SDI}} I(\text{ion at either placement})} $$
Normalizes for abundance; still pools all alternatives.

**(C) Competitive site-determining ratio (recommended).** Mirror the `iso_weight_fraction` logic but
from spectral evidence. For each alternative $a'$:
$$ s(a, a') = \frac{\text{support}(a)}{\text{support}(a) + \text{support}(a')}, \quad
\text{support}(x) = \!\!\sum_{\text{ion} \in \mathrm{SDI}(a,a')}\!\! I(\text{ion at placement } x) $$
The peptidoform score is the **worst case over all alternatives**:
$$ S_C = \min_{a' \ne a} s(a, a') $$
i.e. "the assigned site must beat *every* alternative on its own discriminating ions." $S_C \to 1$ when
the site-determining ions cleanly support $a$; $S_C \to 0.5$ when ambiguous; $S_C < 0.5$ when they
actually favor an alternative.

**Undefined case.** If some alternative $a'$ has **no observed site-determining ions**
($\text{support}(a)+\text{support}(a') = 0$), the site is **unlocalizable against $a'$** from this
spectrum. This is a real, honest outcome (matches the literature: not every site is localizable). The
scorer should return a distinct `n_sdi = 0` flag so we can *report the peptide but not the site*, rather
than silently scoring it 0 or 1.

Recommendation: compute all three but calibrate on **(C)**, reporting **(A)** `n_sdi` alongside as the
"is this even localizable" guard.

# 5. Two-regime localization confidence

The per-peptidoform localization decision becomes:

```
if iso_group_size_at_scan >= 2:      # contested / co-eluting
    confidence = iso_weight_fraction_at_scan      # existing; ~6.6% FLR at 0.95
else:                                # singleton / RT-resolved
    if n_sdi == 0:  report peptide, flag site "unlocalizable"
    else:           confidence = S_C  (site-determining ratio)
```

Both confidences are then thresholded to a **single target FLR**, calibrated jointly against the
target-decoy population (§6). Longer term, (C) can *replace* the fraction everywhere (it is defined for
contested peptidoforms too), but the two-regime split is the minimal first step and lets us keep the
already-validated fraction where it works.

# 6. Calibration & FLR (unchanged framework)

FLR is estimated exactly as now — **target-target vs target-decoy among sequence-targets only,
sequence decoys excluded** (they only set the upstream `qval`/`global_qval`). Every passing
peptidoform, real or loc-decoy, gets the same $S_C$ computed the same way:

$$ \mathrm{FLR}(\tau) = \frac{\#\{\text{target-decoy singletons}: S_C \ge \tau\}}
{\#\{\text{target-target singletons}: S_C \ge \tau\}} $$

If $S_C$ is a good localization discriminator, target-decoys — whose assigned site is wrong, so their
site-determining ions won't match — collapse to low $S_C$, and FLR should fall well below the fraction's
~11.6% singleton floor. **The central empirical question this prototype answers: how far?**

# 7. What we can compute now vs. what needs extraction

| Input | Source | Status |
|---|---|---|
| apex scan index per peptidoform | `new_best_scan` in `precursors_long` | [have] present |
| passing set (qval, global_qval) | `precursors_long` | [have] present |
| target-target / target-decoy label | `is_loc_decoy` join from library | [have] present |
| candidate sites, b/y positions | `site_determining.jl` (`acceptor_positions`) | [have] built |
| library fragment m/z + ion type + position + mod | library `fragments` table | [have] present |
| **observed peaks (m/z, intensity) at apex scan** | the run's Arrow MS2 spectra | [build] needs a read |
| fragment<->peak matching within tol | new targeted step | [build] to build |

So this is **not** a pure `precursors_long` analysis — it is a **targeted extraction**: for each passing
singleton, load its apex scan's peaks from the run's Arrow file (`new_best_scan`), match the library
fragments, classify site-determining, compute $S_C$. It reuses the existing `site_determining.jl`
machinery and runs offline on the current 10k results — **no re-search, no library rebuild**. Estimated
scope: one focused extraction script (~150–250 lines) + reuse.

# 8. Validation plan

1. **Decoy calibration curve** — $\mathrm{FLR}(\tau)$ vs $\tau$ on 10k singletons; compare to the ~11.6%
   fraction floor. Success = a threshold exists with FLR <= 5% at usable coverage.
2. **Ground-truth check on standards** — the 225 Coon standards have known single sites. For the
   singleton subset of detected standards, compare $S_C$-thresholded **empirical** FLR (site vs truth)
   to the **decoy-estimated** FLR. Tests whether the decoy estimate is honest for singletons.
3. **Coverage cost** — IDs retained vs the fraction-only and winner-take-one baselines, at matched FLR.
4. **Unlocalizable rate** — fraction of singletons with `n_sdi = 0` (report-peptide-not-site); this is a
   real property of the data, not a failure.

# 9. Risks & open questions

- **Neutral loss (biggest).** Must exclude/model H3PO4 loss channels or decoys look too good. §3.3.
- **Shadow decoys may be too easy — again.** Decoy sites are non-acceptor "shadow" residues; their
  site-determining ions may be trivially discriminating, so $S_C$ could **under-estimate** FLR just as
  the fraction did in the winner-take-one comparison. The **short-distance-biased** decoy variant is the
  proper calibrant; SDI scoring and better decoys are complementary, not alternatives.
- **Sparse site-determining ions.** Sites near a terminus, or with the discriminating ions below
  detection, give small `n_sdi`. Expect a non-trivial unlocalizable fraction — this is honest, but it
  caps how many singletons we can confidently report.
- **Fragment charge / isotopes / co-isolation.** Match must respect the fragment charge states in the
  library and avoid crediting a co-isolated precursor's peak. Use the same tolerance the search used.
- **$D/R$ normalization** (carried over): at ~1:1 target-decoy:target-target, $D/R$ is the estimate;
  revisit vs $D/(R+D)$ if we change the decoy multiplicity.
- **Does $S_C$ subsume the fraction?** Possibly — (C) is defined for contested peptidoforms too. If it
  matches or beats the fraction there, we unify on one score. To be tested, not assumed.

# 10. Proposed phased plan

- **Phase 0 (cheap, optional first look).** From `precursors_long` alone: quantify the
  contested-vs-singleton split, how many peptides have >=2 RT-separated real forms, and the target-decoy
  rate among singletons vs contested. Sizes the prize before building the extractor. *(~1 script.)*
- **Phase 1 (this proposal).** Build the targeted SDI extractor; compute $S_C$ + `n_sdi` for all passing
  singletons at 10k; produce the decoy calibration curve (§8.1) and the standards ground-truth check
  (§8.2). **Decision gate:** does a usable threshold reach <=5% FLR?
- **Phase 2 (if Phase 1 passes).** Wire the two-regime confidence (§5) into the reporting path; sweep to
  a target FLR; report multi-peptidoform IDs with per-site confidence and an `unlocalizable` flag.
- **Phase 3 (optional).** Replace the fraction with $S_C$ everywhere if it dominates; and/or fold
  $S_C$, `n_sdi`, the fraction, and PSM score into a single trained localization model calibrated on the
  loc-decoys.

---

*Companion docs: `site_determining_fragment_retention.md` (the ion-retention machinery this reuses),
`FINDINGS.md` (the fraction-floor result motivating this), `isomer_decoy_shadow_sites.md`
(the target-decoy construction).*
