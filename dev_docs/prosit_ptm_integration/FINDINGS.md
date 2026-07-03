# FINDINGS — Prosit re-integration

Living record of empirical results per phase. See `PLAN.md` for the design.

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
