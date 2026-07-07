# FINDINGS — Prosit re-integration

Living record of empirical results per phase. See `PLAN.md` for the design.

## Koina Prosit PTM model landscape (probed 2026-07-06)

Empirically tested via Koina `/v2/models/<m>/infer`. **All phospho-capable Prosit intensity models
cap peptide length at 30** (tokenizer limit -> HTTP 400 above 30). Confirmed phospho(UNIMOD:21) +
ox(UNIMOD:35) support and length cap:

| model | phospho | ox | length cap |
|---|---|---|---|
| Prosit_2024_intensity_PTMs_gl | yes | yes | 30 |
| Prosit_2025_intensity_22PTM | yes | yes | 30 |
| Prosit_2025_intensity_40PTM | yes | yes | 30 |
| Prosit_2025_intensity_ptm2 | **NO** | yes | 60 |

So "phospho + length>30" is UNREACHABLE with Prosit (ptm2 lifts length to 60 but drops phospho ->
would need AlphaPeptDeep/MS2PIP). **Prosit_2025_intensity_40PTM is a format-identical DROP-IN** for
2024 (same y1+1 annotation, same mz/intensity outputs; only values differ; Chronologer RT is model-
independent). Integrated into Pioneer registry (src/Pioneer.jl MODEL_CONFIGS + KOINA_URLS) as
prediction_model="prosit_2025_40ptm" — no parser changes (all dispatch is `isa InstrumentAgnosticModel`,
name is just a URL-key label). 22PTM/40PTM = #PTM TYPES supported, NOT length. Model list:
POST /v2/repository/index.

**A/B (same phospho-std digest MC2/7-30/2-4/+Ox, 12-file combined search; data model_compare/):**
2025_40PTM vs 2024_PTMs_gl are EQUIVALENT — total IDs 95,557 vs 95,403 (+0.2%), C>1 decoy-FLR@0.75
(MBR-excl) 9.6% vs 9.5%, empirical@0.75 4.2% vs 3.8% (small-N noise), stratified@<=5% 40,376 vs
41,915. Newer model = safe modernization, no measurable gain here -> keep 2025_40PTM. **The DIGEST is
the real lever, not the model:** phospho-std (MC2/charge2-4/+Ox, yps_std, 34.8M precs) = 95,557 IDs vs
phospho-only MC1 (yps_decoy2, 10.1M) 64,354 (+48%). With Ox: C==1 share 4%->9% (ox adds unambiguous
single-site cases), C>1 91%; decoy-FLR 7.6%->9.6%@0.75 (richer lib genuinely harder); decoy still
conservative vs empirical (~4%). Library C==1 loc-decoys=0 (multi-mod product rule C(cp,np)*C(co,no)
verified). Non-Prosit long-peptide models (AlphaPeptDeep_ms2_generic, UniSpec) do phospho+len40 but
need new annotation parser + instrument_types input — deferred.

## Combined 12-file phospho search vs C>1 decoy library (yps_decoy2) — 2026-07-06

Library yps_decoy2.poin (C>1-restricted decoys). Single COMBINED search over all 12 good dilution
files (experiment-wide global FDR + MBR), localization+L2. Data: ~/prosit_phospho_test/
dilution_curve2/combined/, analyze_combined.jl. **GLOBAL:**
- Total UNIQUE peptidoform IDs @1% FDR: **64,354**. Breakdown: 0 unmodified (lib min_var_mods=1),
  **2,642 (4%) C(c,n)==1** (unambiguous -> FDR only, no loc score), **61,712 (96%) C(c,n)>1**
  (need localization). This 96/4 split is the (b) payoff — correctly separates the no-loc-needed set.
- **C>1 per-file decoy-FLR: 7.6% @frac>=0.75 (MBR-excluded), 9.6% with MBR.** Matches standards
  **empirical ~7.1%** -> decoy method behaves correctly. STABLE ~7.3-7.8% across all 5 dilution points
  (256x); decoy-FLR flat by design (constant yeast background).
- **Empirical at the PROPER stratified operating point (not flat frac>=0.75):** per-C thresholds for
  target decoy-FLR<=10% (C2@0.5,C3@0.5,drop C4-6/7+) -> 118/201 standards, **empirical 1.7%**;
  <=5% (C2@0.56) -> 76/201, **empirical 0.0%**. => **decoy-FLR is CONSERVATIVE (over-estimates)** —
  a safe upper bound. Empirical drops with frac because frac is a genuine confidence (working).

**Two methodology traps found:** (1) NEVER dedup localization by max-frac across files — cherry-picks,
gave a spurious 27%; use PER-FILE instances (-> 7.6%). (2) MBR adds ~2 pts localization error: loc-
decoys MBR-recovered 29% vs reals 10% (decoy iRT==parent gives every decoy an always-present anchor);
MBR genuinely propagates mislocalizations.

**Global vs C-stratified localized instances (C>1, per-file, MBR-excluded; base FLR 70.8%):**

| target | global 1-threshold | C-stratified |
|---|---|---|
| none | 519,465 | 519,465 |
| <=50% | 388,093 | 373,672 |
| <=20% | 204,899 | 183,739 |
| <=10% | 138,349 | 79,157 |
| <=5% | 0 | 35,490 |
| <=1% | 0 | 0 |

Crossover ~5-10%: at loose targets GLOBAL keeps more IDs, but global "<=10%" is a DILUTED average
(high-C localizations inside are ~20-40% wrong, masked by easy low-C). Stratified guarantees EACH
stratum <=target -> uniform per-localization reliability. Stratification's value = honest uniform FLR,
not raw ID count; it's the only option at <5%.

**Library-design notes:** unmodified peptides (min_var_mods=0) don't affect localization/FLR (different
prec m/z, -80/site) — optional, ID-completeness only, sample is phospho-enriched. Standard phospho
search = fixed CAM(C) [have] + variable Phospho(STY) [have] + variable Ox(M) [MISSING] (+ often N-term
acetyl). Adding Ox(M) is more standard and exercises C(c,n) with a 2nd ambiguous mod type; competes for
max_var_mods=2 budget. Also flagged for phospho: max_frag_rank=10 (+gap-cover retention), b_start=2/
y_start=3 excludes b1/y1/y2 (terminal sites unlocalizable), max_frag_charge=3 & include_neutral_diff=
true permitted. Digest currently MC=1/len7-30/charge2-3 (conservative; phospho standard = MC2/7-40/2-4;
NOT charge 1 — phosphate lowers net charge but tryptic peptides stay >=2+).

## Phase 1 (non-PTM parity) — yeast MTAC 3-min, 2026-07-02

**Result: Prosit performs on par with Altimeter.** Phase-1 definition of done met.

### Setup
- **FASTA:** yeast UP000002311 canonical, 6066 proteins (from RIS `FASTA/UNIPROT_2026-06-26`).
- **Two libraries, identical params, only `prediction_model` (+ Prosit NCE) differ:** tryptic
  `[KR][^P]`, 0 missed cleavages, length 7-30, charge 2-3, no variable mods, fixed
  Carbamidomethyl-C, entrapment_r 0, no contaminants, **top-10 fragments/precursor**,
  manual bounds frag 150-2020 / prec 390-1010. Prosit flat **NCE 27**; Altimeter spline.
  - Both: 355,188 precursors, 3,551,880 fragments (10/prec).
  - Build time: Prosit **2.2 min**, Altimeter **6.5 min** (Altimeter's spline-coefficient
    Koina payload is ~2x heavier; Koina is ~75% of Prosit's fragment-prediction time).
- **Search:** the 3 MTAC `Yeast_Alternating-v2_3min` reps
  (RIS `PaperOct2025/MtacYeastAlternating3M/`), **default develop search params**
  (`GetSearchParams`), clean develop code (bitvec WIP stashed), 8 threads, MBR on, 1% FDR.
- **Prosit searchability:** no search-code change needed — develop's `loadSpectralLibrary.jl`
  detects `CompactFrag` fragments → `StandardFragmentLookup` → `FragmentIndexLibrary`, and
  the matcher treats the stored (total-abundance) intensity via
  `getIntensity(::CompactFrag, ::ConstantType)` → `getFragIsotopes!` redistributes it with
  the same `iso_splines` as Altimeter (build-side mono->total is NOT redone at search).

### Results (target-only, qval <= 0.01, precursors matched by sequence+mods+charge)

| Metric | Prosit | Altimeter | Delta |
|---|---|---|---|
| Precursors rep1/2/3 | 47190 / 46891 / 46154 | 47277 / 47037 / 46616 | Altimeter +0.3-1% / rep |
| Unique precursors (union) | 51,237 | 51,391 | -0.3% |
| Precursor overlap | 48,762 shared (95.2% of Prosit, 94.9% of Altimeter) | | |
| Prosit-only / Altimeter-only | 2,475 / 2,629 | | |
| Protein groups (unique) | 4,401 | 4,386 | Prosit +15 (~0%) |
| Quant CV across 3 reps (median / mean) | 16.1% / 20.7% (n=41147) | 16.6% / 21.3% (n=41444) | Prosit marginally better |
| Search runtime | 4.45 min | 4.55 min | ~same |

### Interpretation
- Precursor IDs within 0.3% (unique) / <1% per rep; protein groups essentially identical;
  95% precursor overlap; quant reproducibility marginally *better* for Prosit.
- Achieved at flat NCE 27 (no per-charge tuning) with only 10 fragments/precursor — so the
  parity is not fragile. The strategy-A mono->total design is validated end-to-end on real data.
- Absolute IDs (~47k/rep) are modest partly because of the 10-fragment cap (search stages
  expect more); both libraries are equally capped, so the comparison is fair. A 25-fragment
  rebuild would likely lift absolute numbers for both.

### Open / next
- Phase 1.4 NCE sweep is optional now — 27 already matches Altimeter; a sweep might gain a
  little but is not needed for parity.
- Per-charge NCE adjustment remains deferred (§5.2) — not needed for parity.
- Ready to proceed to Phase 2 (PTM / `Prosit_2024_intensity_PTMs_gl`) when desired.

### Artifacts (local, not in repo)
Working dir `/Users/n.t.wamsley/prosit_yeast_benchmark/`: `{prosit,altimeter}_yeast.poin`,
`results_{prosit,altimeter}/`, `build_*.json`, `search_*_default.json`, `compare.jl`.

## Phase 2/3 — first phospho search on REAL data (2026-07-03)

**Prosit-PTM phospho pipeline validated end-to-end on real Astral DIA phospho data.**

### Data
- Coon lab Astral phospho benchmark = **MassIVE MSV000093613** (DOI 10.25345/C5D795N33), *not*
  PXD049028 (the benchmark brief mis-attributed it — PXD049028 is HAP1/MGME proteome).
- Benchmark files: `20240306_NML_15min_2Th_250ngYeast_{10000,2500,625,156,39}amol_Stds0{1,2,3}.raw`
  = 250 ng yeast phosphopeptides + synthetic standards (JPT SpikeMix 52/54 + Sigma PhosphoMix
  1/2/3 = "JPTstySig123"), 5 amol dilution points x 3 reps, 15-min 2Th Astral DIA.
- Downloaded the 10000 amol triplicate (~21 GB) via MassIVE HTTPS (~20 MB/s).

### Pipeline (all local, working)
1. `.raw` -> `.arrow` via **PioneerConverter** (~/Projects/PioneerConverter, dotnet 8, built
   `dotnet build -c Release`): 7.4 GB raw -> 3.37 GB arrow in ~57 s; 250,866 scans
   (2,483 MS1 + 248,383 MS2). No mzML step.
2. Search with the yeast **phospho-only** library (`prosit_yeast_phospho_only.poin`, 2.05M
   precursors), default develop params, 8 threads.

### Result (Stds01, single 15-min run, 1% FDR)
- **24,267 phosphopeptide precursors** (100% carry >=1 Unimod:21), **21,284 unique
  phosphopeptides**, **25,802 phospho-sites**, 2,358 protein groups. Runtime 4.82 min.
- Multiplicity: 22,732 mono- + 1,535 di-phospho (0 tri — max_var_mods=2 cap).

### Notes / next
- The **225 synthetic standards** are NOT in the yeast library (they're human/synthetic); their
  sequences are needed for ground-truth ID+localization. Not in the public .xlsx supplements
  (results tables only) — likely in the SI PDF (MOESM1) or the commercial JPT/Sigma kit sheets.
- Next: 3-rep MBR search (all 3 reps converting); site-localization scoring (Phase 3); add the
  225 standards for a real FLR benchmark; dilution series for LOD/FLR curve.

### Update: 225 standards found + 3-rep MBR (2026-07-03)
- **Paper:** "Fast and deep phosphoproteome analysis with the Orbitrap Astral mass spectrometer",
  Nature Communications 2024, DOI 10.1038/s41467-024-51274-0 (Coon lab).
- **225 synthetic standards = Supplementary Data 1 (MOESM3), sheet "Reference Sheet- Peptide
  Level"** (225 peptides, sequence + known Site in Spectronaut notation e.g. Y6; multiplicity
  210 mono / 11 di / 3 tri / 1 tetra). Also "Reference Sheet - Site Level" (245 sites) +
  "Dilution Series Data". Saved to `/Users/n.t.wamsley/prosit_phospho_test/coon_225_standards.tsv`.
  This is the ground truth for ID + false-localization-rate benchmarking.
- **3-rep MBR search** (all 3 x 10000amol reps, phospho-only yeast lib): 75,016 phosphopeptide
  precursors, 7,208 protein groups, 5.85 min.
- Next: add 225 standards to a library (with known sites), search -> ID recovery + FLR; then a
  site-localization scorer (Phase 3) and the dilution series for LOD/FLR curve.

### Ground-truth benchmark: ID + FLR with 225 standards (2026-07-03)
Built a **combined library** = yeast phosphoproteome + synthetic standards, 1 MC, variable
phospho [STY] max_var_mods=2 min_var_mods=1 (all positional isomers so localization can fail):
5,306,732 precursors, 208/212 standard sequences present (built 27 min). Standards FASTA =
217 with <=1 internal K/R (172 fully tryptic + 45 with 1 MC; 8 with >=2 internal deferred).
Truth = coon_225_standards.tsv "Reference Sheet" (sequence + known site).

**Search Stds01 (10000 amol, highest load) vs the combined library, default develop params,
1% FDR** (42,116 total phospho precursors, 5.69 min):
- **ID recovery: 198/212 standards = 93.4%** (mono + a few multi).
- **Localization (best-scoring positional isomer per standard vs known site): 152/198 correct
  -> FLR = 23.2%** (mono 22.3%; multi 2/5). Allowing any valid site for positional-isomer
  standards barely changed it (24.2% -> 23.2%).

**Interpretation:** Prosit-PTM phospho *identification* is strong (~93% ground-truth recovery);
*localization* is the gap (23% FLR vs the field's 1-5%) because Pioneer has **no dedicated
site-localization scorer** — "best isomer wins" is a crude proxy and DIA positional isomers
share most fragments. This is the concrete motivation + harness for **Phase 3 (site-localization
scoring / AScore-style site-determining-ion model)**. Analysis: `analyze_flr2.jl`. Library/data
in `/Users/n.t.wamsley/prosit_phospho_test/` (yeast_plus_standards_1mc.poin, search_std/).

### Experiment: PMM deconvolution L1/L2 penalty sweep (2026-07-03, commit 41f26619c)
Added config-driven L1/L2 penalty to the PoissonMM main-search solver. Swept on the phospho
FLR benchmark (Stds01 vs yeast+standards 1MC lib, 1% FDR):

| reg | lambda | total precursors | standards ID | correct | FLR |
|---|---|---|---|---|---|
| none | 0.0 | 42,116 | 198/212 | 152 | 23.2% |
| L1 | 0.1 | 42,492 | 199/212 | 152 | 23.6% |
| L1 | 1.0 | 42,606 | 199/212 | 152 | 23.6% |
| **L2** | **1.0** | **47,426** | 199/212 | **160** | **19.6%** |

- `none` reproduced the baseline exactly -> the L1/L2 change is bit-identical for NoNorm (no regression).
- **L2 ridge (lambda=1) improved localization: FLR 23.2% -> 19.6%** (real ground-truth gain) AND
  +12.6% total precursors. The extra precursors need FDR validation (entrapment/decoy) before
  claiming a depth gain, but the FLR drop is a clean quality signal.
- L1 had ~no effect (NNLS max(.,0) already induces sparsity).
- Interpretation: L2 shrinkage stabilizes the dense-DIA phospho deconvolution (many co-isolated
  positional isomers), yielding steadier weights -> correct isomer wins more often. Promising
  lever for phospho localization; worth a dose-response + FDR check.

### L2 dose-response (refines the above)
Full L2 sweep on the same benchmark:

| reg | lambda | total precursors | correct | FLR |
|---|---|---|---|---|
| none | 0.0 | 42,116 | 152 | 23.2% |
| L2 | 0.1 | 43,742 | 141 | 28.8% |
| L2 | 1.0 | 47,426 | 160 | 19.6% |
| L2 | 2.0 | 48,619 | 154 | 22.6% |
| L2 | 5.0 | 50,186 | 156 | 21.2% |

- **Total precursors rise MONOTONICALLY with lambda** (42k->50k) — robust effect (L2 smoothing gives
  more precursors non-zero weight); needs FDR (entrapment/decoy) validation before claiming depth.
- **FLR is non-monotonic / noisy**: baseline 23.2%, best 19.6% at lambda=1, but lambda=0.1 is WORSE
  (28.8%). Spread ~9% over ~198 standards (~18 peptides). CORRECTION to the single-point claim
  above: L2 has a lambda~1 localization sweet spot but the FLR gain is NOT robust from one replicate.
- Takeaway: L2 reliably increases ID count and modulates localization with a lambda~1 optimum, but
  a real localization improvement needs replicates + FDR/entrapment checks. Dedicated site-
  localization scoring (Phase 3) remains the principled fix; L2 is a cheap knob worth revisiting.

### KEY RESULT: deconvolution weight localizes far better than q-value (2026-07-03)
On the same NoNorm phospho benchmark + ground truth, re-scoring localization by the *isomer's
deconvolution weight* instead of its q-value (no new model, no re-search — just a different column):

| localization rule | correct | FLR |
|---|---|---|
| best q-value (prior proxy) | 152/198 | 23.2% |
| **highest deconvolution weight** | **187/198** | **5.6%** |
| best weight_rank_at_scan | 175/198 | 11.6% |

**FLR 23.2% -> 5.6%** — into the paper's 1-5% range, using info Pioneer already computes. Rationale:
the q-value compares isomers' OVERALL scores (discriminating signal buried under shared fragments);
the deconvolution weight is the isomer COMPETITION (solver splits shared signal; split driven by
site-determining ions weighted by Prosit-PTM predicted intensities) = a direct, intensity-aware
localization readout. Raw weight > weight_rank (rank is vs all precursors at scan, not just isomers).
The two rules disagree on 51/198; weight wins most. `weight`/`weight_ratio_at_scan`/
`weight_rank_at_scan` are already in precursors_long.arrow.

Headroom (Phase 3 localizer): weight-FRACTION (w_best / sum w_isomers) as a per-site confidence +
cutoff (field uses >=0.75) -> trades recall for lower FLR toward ~1%; integrate weight over the peak;
true shared-apex version (all isomers at the same apex scan). Caveats: 1 replicate; `weight` is each
isomer's own best-scan weight, not strictly the shared apex.

### Phase A isomer-competition feature built + reg sweep (2026-07-03, commit b31302ebf)
Implemented `add_isomer_competition_features!`: per-PSM `iso_weight_fraction_at_scan` = weight /
sum(weights of SIBLING isomers at the same scan), siblings = (sequence, charge, mz-bucket,
is_decoy). Gated on `optimization.localization.enabled`; propagates to precursors_long.arrow. Full
per-standard table: `isomer_competition_standards.tsv` (193 contested groups).

**The fraction is a CONFIDENCE, not a selector** (10k-amol Stds01, 198 IDs):

| localization rule | FLR |
|---|---|
| best q-value | 20.7% |
| **max apex weight** (selector) | **6.6%** |
| max apex fraction (selector) | 18.7% (BAD) |

Max-fraction selection is confounded: a weak WRONG isomer that finds a scan with no co-eluting
sibling gets fraction=1.0 trivially (wrong maxfrac picks have median conf 1.0). So SELECT by max
weight, then cut on the winner's fraction as a filterable confidence:

| tau | FLR | coverage |
|---|---|---|
| 0.0 | 6.6% | 100% |
| 0.75 | 2.5% | 81.3% |

2.5% @ 0.75 is the PTM.SiteProbability>=0.75 analogue and beats the paper's Astral 3.95% @ 0.75
(exact, from Nat Commun 2024 Source Data Fig2I; caveat: uncalibrated fraction vs calibrated
probability, different operating points).

**Deconvolution reg sweep (localization on; matched ~80% coverage) — ridge helps, L1 doesn't:**

| run | FLR | coverage |
|---|---|---|
| no reg (tau=0.75) | 2.5% | 81.3% |
| L2 lambda=0.5 (tau=0.60) | 1.8% | 83.3% |
| **L2 lambda=1.0 (tau=0.60)** | **1.3%** | 80.2% |
| L1 lambda=0.5 (tau=0.75) | 3.0% | 83.3% |

Ridge sharpens the sibling split monotonically with lambda (1.3% FLR at 80% cov, below the paper's
3.95%); L1 neutral-to-worse. Confirms the collinearity hypothesis (sibling columns share all but the
site-determining ions). Caveat: L2 rescales the fraction -> compare at matched coverage, not matched
cutoff. **A light ridge default belongs in the localization path** — and becomes mandatory once
isomer decoys add near-collinear columns.

Residual errors are all ADJACENT-site (AALDHSCSPSK S6-vs-true-S8 @ conf 0.962; DLATVYVDVLK
T4-vs-Y6 @ 1.0): the wrong isomer finds a lonely scan -> high confidence despite being wrong. These
are the ~2.5% the cutoff can't catch -> next step: **isomer decoys + a learned localization
discriminant** (see `localization_decoy_library_plan.md`).

### Isomer-decoy library + FLR-vs-abundance dilution series (2026-07-06)

Built the shadow-site isomer-decoy library (yeast+standards, 10.6M precursors = 5.3M reals + 5.29M
loc-decoys, move-one-random onto count-matched random shadow sites), searched the Coon Astral
phospho **dilution series** (MSV000093613: 250 ng yeast + synthetic standards spiked at 10000 -> 39
amol, 3 reps/point; 3 of 15 raw files are CORRUPT in the MassIVE deposit -> those points use 2
reps). Search = 3-rep MBR + L2(1.0) + localization; decoy-based FLR (loc-decoy wins / real wins per
isomer group) vs empirical FLR (standards ground truth). Reported = max-weight isomer per group,
cut on the winner's iso_weight_fraction.

**FLR-vs-abundance (wt-ratio >= 0.5, reported basis):**

| load | std detected (wt>=0) | reported std IDs | empirical FLR | decoy FLR | global IDs |
|---|---|---|---|---|---|
| 10000 amol | 201 | 160 | 3.8% | 8.8% | 12,401 |
| 2500 amol  | 196 | 165 | 3.6% | 9.3% | 13,400 |
| 625 amol*  | 187 | 158 | 2.5% | 8.9% | 12,602 |
| 156 amol*  | 166 | 143 | 3.5% | 8.7% | 12,639 |
| 39 amol*   | 122 | 108 | 2.8% | 9.6% | 12,551 |

(* 2 reps; the 3rd rep file is corrupt on MassIVE.)

**Three findings:**
1. **Global IDs + decoy-FLR are flat by design** — the yeast background is constant (250 ng in every
   run); only the 225 standards dilute. So ~12-13k global phosphopeptides and ~9% decoy-FLR (both
   yeast-dominated) SHOULD be flat, and are. Sanity check passed.
2. **Abundance hits DETECTION, not localization accuracy of what's detected.** Standards detected
   fall 201 -> 122 (10k -> 39 amol; ~61% retained at 256x dilution, ~ the paper's "half at 39
   amol"), but the empirical FLR of the detected standards stays ~2.5-3.8% across the whole range.
   Localization doesn't degrade with abundance; sensitivity does.
3. **Decoy-FLR (~9%) ~= 3x the standards empirical (~3%) — and that gap IS the decoy method's value.**
   The 225 standards are easy (pristine synthetic peptides); the decoy-FLR estimates the hard, real
   yeast phosphoproteome, which is genuinely ~9% and has NO ground truth. Stable and believable
   across abundance -> exactly the scenario decoys exist for (FLR on the real sample, no standards).

Caveats: decoy design is random/shadow (move-one-random), which UNDER-shot the standards' known FLR
in the single-point comparison (decoys too easy) -> a short-distance-biased draw is the next tune.
And this is per-reported (max-weight winner per group); per-precursor rates are higher.
Pipeline/data: ~/prosit_phospho_test/{run_dilution_curve.sh, dilution_curve/SUMMARY.tsv, yps_decoy.poin}.
See memory: [[localization-isomer-competition]], [[data-download-and-convert]].

---

## Canonical terminology (standardize on this — never write bare "decoy" for localization)

Two orthogonal decoy axes, two questions, two flags already in the code:

| Term | Question | FDR estimated | Code flag | Construction |
|---|---|---|---|---|
| **Sequence decoy** | Is this *peptide* real? | ID / detection FDR | `is_decoy` | reversed / shuffled sequence |
| **Positional isomer decoy** | Did we pick the correct *positional isomer*? | localization FDR (FLR) | `is_loc_decoy` | same sequence/mass, mod on a wrong position |

("Positional isomer decoy" is the precise name; "peptidoform decoy" / "loc-decoy" are older synonyms
for the same `is_loc_decoy` flag.) A positional isomer decoy winning the assignment = we placed the
mod on the wrong position = a localization error; **FLR = decoy-isomer-wins / real-isomer-wins**.

**Existence condition (falls straight out of the definition):** a positional isomer decoy can only
exist if there are >=2 positional isomers to choose among — i.e. for some mod type X, candidate
acceptor sites `c_X` > mods present `n_X`, so `C(c_X, n_X) > 1`. Single-site-per-mod peptidoforms
(1 phospho / 1 S-T-Y; phospho + 1 oxidizable M; 2 phospho / exactly 2 sites) have ONE arrangement,
nothing to mislocalize -> **no decoy**. This is not a special case; it is what "positional isomer
decoy" means. Consequence: the current library over-generates (decoys on single-site peptides), which
only measures ID-level leakage (the 3.1% 1-site floor), not localization error — restrict generation
to `c_X > n_X` peptidoforms.

**Construction tension (open):** a *true* positional isomer keeps the mod on a valid acceptor at a
different position. The current **shadow-site** design puts it on a *non-acceptor* residue — a
"position" that is not a real isomer, chemically implausible, easy for the search to reject -> decoys
too easy -> FLR under-shoots the standards. The ideal positional isomer decoy sits on a *plausible
alternative acceptor-like position* (the short-distance-biased draw reaches toward this). Note the
irreducible difficulty: a *real* alternative acceptor site is itself a legitimate target isomer, so
the decoy must be a position the search will entertain yet that we can label wrong — the shadow-site
is the current compromise.

Peptidoform-decoy **construction methods** (name the axis + method, e.g. "peptidoform decoy
(shadow-site move-one-random)"):
- **shadow-site move-one-random** — CURRENT (`yps_decoy.poin`). Per peptide draw count-matched random
  *non-acceptor* shadow sites (deterministic per-seq seed); per modified peptidoform relocate ONE
  random variable mod to a random shadow site (partial displacement, siblings kept). Verified in the
  library: e.g. `(1,Y)(4,Y)->(1,Y)(15,D)`, `(5,S)(6,C:CAM)(21,S)->(3,A)(5,S)(6,C:CAM)`.
- **adjacent move-one** — retired ±1 draft.
- **short-distance-biased** — proposed next (same shadow scheme, draw biased to the real site-distance
  distribution, to stop under-shooting the standards' known FLR).

## How the localization score / FLR is computed (target-decoy competition)

1. **Group** by `sequence | charge | precursor-m/z`. A peptidoform decoy has identical precursor m/z
   (phospho mass is residue-independent), so it competes INSIDE its reals' group.
2. **Select + score:** winner = `argmax(weight)` in the group (real or decoy). Localization SCORE =
   winner's `iso_weight_fraction_at_scan` (winner weight / Σ sibling isomer weights at apex scan).
   Selection is by weight; the fraction is the winner's confidence (bad selector, good confidence).
3. **Threshold + FLR:** sweep cutoff `t` on the winner fraction. Among groups with winner `frac >= t`:
   `R`=real winners, `D`=decoy winners, **decoy-FLR = D/R**; reported IDs = unique real-winner
   peptidoforms. A decoy winning its group = proxy for a real mislocalization.
   OPEN: normalization is `D/R` at ~1:1 decoy:real — revisit vs `D/(R+D)` / `2D` if the ratio changes.
   Rate is per-reported (one winner/group); per-precursor FLR is higher.
   Code: ~/prosit_phospho_test/analyze_dilution_point.jl.

## Multi-peptidoform reporting: fraction alone floors at ~10.7% FLR (2026-07-06)

Goal: report >1 peptidoform per peptide (multiple positional isomers can be genuinely co-present)
instead of one max-weight winner per group. Data already supports it — each isomer is its own
precursor_idx, deconvolved per-scan + FDR'd independently; `precursors_long` carries each
peptidoform's best-PSM `iso_weight_fraction_at_scan` (its apex value). So: drop the argmax collapse,
filter to real targets passing BOTH `qval<=0.01` (run) AND `global_qval<=0.01` (experiment-wide),
dedup to unique precursor_idx (max frac across reps), sweep the fraction cutoff, FLR = loc-decoys /
reals passing. Ran on 10k-amol (3 reps): 54,066 real + 40,814 loc-decoy peptidoforms post-qval.

| frac >= t | all: real / decoy / FLR | contested(gsz>=2): real / decoy / FLR |
|---|---|---|
| 0.75 | 11,762 / 1,348 / **11.5%** | 6,046 / 683 / 11.3% |
| 0.90 | 7,987 / 851 / 10.7% | 2,271 / 186 / 8.2% |
| 0.95 | 6,831 / 739 / 10.8% | 1,115 / 74 / **6.6%** |
| 1.00 | 5,716 / 665 / 11.6% | 0 / 0 / — |

**Fraction alone CANNOT push per-peptidoform FLR below ~10%** (flat from 0.75->1.0). Cause:
**singletons**. At frac>=0.75, ~half the confident calls (5,716 real + 665 decoy) are singletons
(`iso_group_size_at_scan==1`) — alone in their group that scan, so frac=1.0 AUTOMATICALLY and carries
ZERO localization info; their ~11.6% FLR IS the floor. A peptidoform is a singleton at its apex
exactly when its siblings elute at a DIFFERENT RT — i.e. the "genuinely co-present isomers" case we
want to report. **The sibling-fraction measures co-elution competition; RT-separated isomers don't
compete, so a peptidoform decoy that separates in RT with a clean permuted spectrum is
indistinguishable from a real second form.** Multi-reporting therefore has intrinsically higher FLR
than winner-take-one, and the fraction can't close it. Restricting to contested (gsz>=2) reaches
~6.6% at frac>=0.95 but costs most IDs and still isn't 1%.

**Design implication — two regimes, two scores:** (a) contested/co-eluting peptidoforms -> the
iso_weight_fraction is the right localization confidence (~6.6% at 0.95); (b) singleton/RT-resolved
peptidoforms -> fraction is uninformative, confidence must come from **site-determining-ion evidence**
(b/y ions distinguishing the candidate sites) and/or the PSM score, thresholded against loc-decoys.
The singleton FLR is a detection/PSM-score property (did the search reject the decoy's permuted
spectrum?), not a competition property. Next: prototype a per-peptidoform site-determining-ion score
for the singleton regime — that's the actual bottleneck for low-FLR multi-reporting.
Prototype code: inline (this session); numbers from dilution_curve/amol_10000/results.

## Phase 0: what singletons actually are — deconv collapse + acceptor-count floor (2026-07-06)

Investigated WHY singletons (`iso_group_size_at_scan==1`) exist. `iso_group_size_at_scan` counts
sibling PSM rows in the deconv output at the apex scan; `Score!` (ScoredPSMs.jl:27) DROPS any
precursor with `weight < 1e-6`. So a singleton = the only isomer of its peptide with NONZERO
deconvolution weight at that scan. But positional isomers share every non-site-determining b/y ion,
so siblings are co-selected as fragment-index candidates and enter the SAME design matrix; the L2
ridge collapses the group onto one winner. 4,374 / 5,716 real singletons (77%) are >=2-acceptor-site
peptides where a sibling isomer provably exists. **=> deconvolution is already doing localization;
`gsz==1` mostly means "the ridge collapsed this isomer group to one winner," NOT "nothing competed."
frac=1.0 for a singleton is therefore NOT information-free (earlier claim too strong) — it encodes a
site choice driven by site-determining ions + L2.**

**CANDIDATE CHECK (single-file diagnostic, env-gated dump of per-scan fragment-index candidates,
scan_idx aligned 100% self-hit; hook reverted):** Of 2,654 real singletons, at the gsz scan:
**54% had ZERO sibling positional isomer as a candidate** (median sibling-candidates=0); only 46% had
>=1 sibling candidate (39% among those with a real sibling). So the "deconvolution collapsed a
co-candidate group" story holds for only ~46%; the other ~54% were **genuinely uncontested — the only
positional isomer even proposed at that scan.** WHY siblings weren't candidates is NOT retention time:
predicted-iRT |gap| singleton-vs-non-candidate-sibling is small (median 0.62, only 5% >2 iRT), ~same
as candidate siblings. Since RT puts them in-window, siblings fail at the **fragment-index candidate
score** — their mod-bearing (site-determining) ions are at the wrong m/z, so they don't accumulate
enough matched fragments. **=> localization happens at TWO implicit stages, both keyed on
site-determining ions: (1) candidate selection (the fragment index filters wrong-position isomers
before deconv — handles the 54%), (2) deconvolution ridge (collapses the 46% that still competed).**
Corrects the earlier framing that deconv alone does localization. This STRENGTHENS the SDI-score case:
a per-peptidoform site-determining-ion score just formalizes what both stages already use implicitly.
And it re-explains the shadow-decoy weakness: a positional isomer decoy only registers as a false
localization if it PASSES candidate selection, but shadow (non-acceptor) decoys often fail the
fragment-index score (mod ions at implausible m/z) -> they don't compete -> FLR under-shoots. A useful
decoy must be "matchable enough to be proposed, yet labelable wrong."

## C(c,n) — the localization-multiplicity difficulty prior + decoy-iRT=parent (2026-07-06)

FLR is only needed when >=2 positional isomers exist at a precursor m/z, i.e. `C(c,n) > 1` where
c = candidate acceptor sites (#S/T/Y), n = mods present (#phospho). Single-site (C(1,1)=1) AND
fully-saturated (C(2,2)=1) are UNIQUE at their m/z -> **FDR only, no FLR**. Targets with no
localization question shouldn't be in the FLR at all (neither the real nor its shadow decoy).

**Removing non-localizable (C==1) peptidoforms RAISES the FLR** (they were the easy low-FLR stratum
diluting it): all vs localizable-only at frac>=0.9 = 10.7% -> 12.4%; frac>=0.95 = 10.8% -> 13.0%.
So the honest FLR on genuinely-localizable peptidoforms is ~12-13%, not ~11%. The C==1 removal is a
bias correction, not an improvement.

**C(c,n) is a strong, monotonic difficulty prior — and it rescues the singleton regime** (where
frac=1 is information-free). Singleton FLR by C: C1=3.2%, C2=6.6%, C3=11.2%, C4-6=20.3%, C7+=38.0%.
So a lone uncontested call's trustworthiness is set by C even though its fraction says nothing.
**Proposed FLR-model features (cheap, library-only, no extraction): #sites c, #phosphorylated n,
combinations C(c,n) (or log C = chance-baseline info content 1/C).** c,n not fully redundant with C
(same C, different site-determining geometry). Two multiplicities: C(c,n)=theoretical/prior (always
defined incl. singletons) vs iso_group_size_at_scan=actual competition (evidence); C is the prior,
fraction/SDI the evidence -> combine = localization posterior. Usable as model feature OR stratified
per-C thresholds (near-free). SDI evidence matters most in high-C strata (4-6, 7+) where prior alone
leaves 20-38%.

**Decoy iRT == parent iRT, exactly (100%, 5,285,076/5,285,076, max gap 0.0):** the generator copies
the parent's retention time (Koina can't predict shadow-site phospho iRT). => a positional isomer
decoy is ALWAYS RT-eligible wherever its parent elutes (only its permuted fragments can exclude it —
isolates the localization decision to site-determining ions), but is structurally locked to the
co-eluting regime (can never be an RT-resolved singleton). Real positional isomers DO shift RT
(median ~0.6 iRT between real siblings). For the core FLR question ("could a wrong site win at THIS
peak?") evaluating at the parent RT is correct; but it's on the decoy-construction "reconsider" list
alongside shadow-site and C>1-restriction. Numbers: dilution_curve/amol_10000 (3 reps).

Multi-form reality (per-run, real): 63% of detected phosphopeptides have >=2 peptidoforms passing
FDR, but ~96% CO-ELUTE (median closest-pair |Δrt| = 0.0 min); only ~1,230 peptide-instances (~2%)
are RT-resolved. So genuine "co-present isomer" multi-reporting is a small prize; most multi-form is
co-eluting isomers already in the contested regime.

**Singleton FLR is acceptor-count-driven (the real story behind the ~11.6% floor):**

| # S/T/Y | real singletons | decoy | FLR | note |
|---|---|---|---|---|
| 1 | 1,342 | 41 | **3.1%** | localization trivial (one possible site); 3.1% = pure ID-level leakage of shadow decoys on non-acceptor residues |
| 2 | 1,562 | 103 | 6.6% | |
| 3 | 1,190 | 133 | 11.2% | |
| 4+ | 1,622 | 388 | **23.9%** | genuinely hard; deconv collapse unreliable |

**Implications:** (1) SDI-score target is the >=3-site singletons (11-24%), NOT singletons broadly;
1-site are already ~3% with nothing to localize. (2) Localization difficulty scales cleanly with
acceptor count -> **acceptor-count-stratified per-peptidoform thresholds (or per-stratum FLR
reporting) is a cheap win with no new extraction** — try before the full SDI extractor. (3) The 3.1%
1-site floor is the decoy method's own ID-leakage baseline (shadow decoys occasionally too matchable)
-> relevant to the short-distance-biased-decoy retune. Numbers: dilution_curve/amol_10000 (3 reps).
NOTE: "siblings were co-candidates" is inferred (shared fragments + acceptor-count evidence);
confirm definitively by re-running one file with delete_temp=false + localization and inspecting
design-matrix membership if certainty is needed.
