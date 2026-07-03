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
