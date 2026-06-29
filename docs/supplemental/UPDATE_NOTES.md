# Supplemental.tex Update Tracking

Goal: update `supplemental.tex` to reflect current Pioneer (develop @ `1eb5c1886`, 2026-06-27).
Source doc predates a large set of pipeline changes. Status per section tracked below.

## Substantive PRs since late Jan 2026 (algorithm-affecting), grouped by doc section

### STRUCTURAL — pipeline reorganization (highest priority)
- #375 `bypass-first-pass`, #361 `scoring-refactor` — **FirstPass/SecondPass consolidated into MainSearch**.
  → The doc's "First Pass Search" (§subsec:first-pass-search) and its scoring/aggregation/RT subsections are obsolete as written.
- #300 `oom_percolator` — out-of-memory percolator scoring (streaming, per-file).

### MBR / Target-Decoy / FDR  (§subsec:target-decoy-mbr-fdr — largest rewrite)
- #419 `optimize-flanking-window-features`, #424 `remove-unused-flanking-correlated-features`,
  #428 `flanking-window-main-search-speedup` — flanking-window MBR features (added/pruned/moved).
- #431 `file-aware-mbr-decoys`, #407 `remove-mainsearch-pair-competition` — target-decoy pairing changes.
- #433 `mbr-empirical-spectra-hellinger` — NEW MBR feature (empirical-spectra Hellinger distance).
- #434 `mbr-protein-scoring-features` — MBR features in protein scoring.
- #420 `combined-m1-lgbmc-pass1`, #429 `nombr-oom-streaming` — Pass-1 OOM LGBM trainer / no-MBR streaming.
- #329 `disable-MBR`, #330/#331 `three-nonMBR-iterations` — non-MBR iteration path.
- #297 `deterministic_lightGBM`, #309 `global-pep-threshold` — scoring determinism, global PEP threshold.

### Parameter Tuning  (§subsec:parameter-tuning)
- #409/#410 `rt-aware-mass-calibration` — **mass-error estimation rewrite** (RT-aware; regularized bias/spread splines).
- #411 `fast-ms1-features` — MS1 mass-error fit / features.
- #435 `profile-parameter-tuning` — plots only (no algorithm change; ignore for doc).

### Library / Fragments  (§subsec:fragment-index-search, fragment isotope)
- #299 `serialization-for-fragments` — fragment serialization (`fragments.jls`).
- #432 `build-spec-lib` — spectral library build.
- #408 `c13-isotope-spacing`, #415 `sum-fragment-isotope-intensities`, #416 `add-fragment-detected-counts`.

### Chromatogram / Quant / MaxLFQ  (§subsec:chromatogram-smoothing-integration, maxlfq)
- #365 `maxLFQ-disconnected`, #376 `intensity-spline-lipschitz-step`, #389 `quant-change1-q010`,
  #381 `ols-vs-pmm-comparison` — deconvolution/quant solver changes.
- #400 `document-chromatogram-integration-params` — integration params.

### Protein Inference  (§subsec:protein_inference)
- #379 `global-protein-inference`, #377 `decoys-through-protein-scoring`,
  #378 `protein-group-species-split`, #364 `protein-name-parsing`.

### Perf only — NO algorithm change (do NOT alter doc for these)
- #422 runtime-profiling, #423 combined-allocations, #427 clean-for-develop, #430 integrate-chrom-hoist.

## Section-by-section status
(legend: TODO / IN-PROGRESS / DONE / NO-CHANGE-NEEDED)

| Section | Status | Notes |
|---|---|---|
| Altimeter Training Data (§4.1) | TODO | check vs current library build |
| File Conversion (§4.2) | TODO | #312 converter CLI |
| Spectral Library Sequence Gen (§4.3) | TODO | |
| Intensity-Aware Fragment Index Search (§4.4) | TODO | #299 serialization |
| Parameter Tuning (§4.5) | TODO | #409/#410 mass cal rewrite, RT align |
| **First Pass Search (§4.6)** | TODO | **OBSOLETE — now MainSearch (#375/#361)** |
| Matrix Representation (§4.7) | TODO | |
| Linear Regression of Spectra (§4.8) | TODO | #376/#381 solver |
| **Target-Decoy / MBR / FDR (§4.9)** | TODO | **largest rewrite — MBR PRs above** |
| Chromatogram Smoothing/Integration (§4.10) | TODO | integration changes |
| Protein Inference & Quant (§4.11) | TODO | #379 global inference, maxLFQ |
| Fragment Isotope Correction (§4.12) | TODO | |
| Quadrupole Transmission (§4.13) | TODO | |

## Section progress log
- §1 intro: DONE — 'first-pass searching' -> 'the main search' (consolidation #375/#361). NOTE: rename §4.6 title 'First Pass Search' -> 'Main Search' to match.
- §1.3 Spectral Library Sequence Gen: DONE — corrected the entrapment/decoy reversal claim (code only forbids BOTH being reverse, via check_params.jl:98). OPEN: doc lists 2 decoy methods; code has 3 (shuffle/reverse/diann_mutation) — pending user decision.
- §1.4.1 Index Construction: bitmask scoring DONE (top-8, 2^(r-1), OR); ion-mobility line removed. OPEN: per-fragment fields (P,M,Z,S)->(P,S); partition layer; isotope+precursor-ion removal; within-bin precursor-mz ordering.
- §1.4.2 Score Counter: DONE — W now bitmask indexed by precursor ID, accumulated by bitwise OR (Counter{UInt16,UInt8}). OPEN: 'partition-local' ID wording (coupled to partition decision); title 'Score Counter' could be renamed since value is a bitmask not a score.
- §1.4.1: partition layer ADDED (compact (P_k,S_k) record, 5 Da partitions/65535 split, isotope+precursor-ion removal, dropped within-bin prec-mz ordering). DONE.
- §1.4.3 MS2 Scan Query: REWRITTEN (partition selection, SIMD m/z bin search, bitmask OR, candidate selection via learned LUT, precise m/z/iRT deferred to deconvolution; removed sort+threshold T). DONE.
- §1.4.4 Bit-Vector Calibration: NEW subsection written (256-bin target/decoy tally -> CI pass/fail + coarse-group merge -> single pooled table). DONE. Reconciled: dedicated phase between QuadTuning and MainSearch; ParameterTuning/QuadTuning use rule-based bootstrap LUTs.
- TODO §1.5 Parameter Tuning intro still says 'first pass search' + 'score threshold' for fragment index — update when we reach it.
- §1.5 intro: DONE — first-pass->main; score-threshold reframed (min matched top fragments); MS1 mass error added; scan sampling corrected to RT-stratified (30 bins) round-robin by descending TIC (FilteredMassSpecData.jl:175 get_ms2_scan_priority_order); MS2 mass error -> spline model (IntensityMassErrorModel: mz+intensity+rt bias splines).
- TODO §1.5.2 Mass Error Estimation: still OLD (median bias + 2 exponential tols). Rewrite to IntensityMassErrorModel: bias_Da(mz,I,rt)=mz_bias_spline+intensity_bias_spline+rt_bias_spline; spread spline -> tolerance (k*sigma); SimpleMassErrorModel scalar fallback. Needs deep investigation.
- §1.5 'single retention time bin' presearch index: VERIFIED ACCURATE — presearch index built with rt_bin_tol=typemax(Float32) (build_poin_lib.jl:211) = one RT bin, no binning. No change.
- §1.5.2 Mass Error Estimation: REWRITTEN to IntensityMassErrorModel. bias(m,I,t)=g_m(m)+g_I(log2 I)+g_t(t), each a regularized cubic B-spline by IRLS (Huber); m/z & RT on binned medians, intensity on raw w/ monotonicity; RT term only if enough gradient coverage; linear extrap. Tolerance = k*sigma(log2 I)*c(mz) (intensity spread spline x m/z correction), k~1.96 (~95%). Scalar SimpleMassErrorModel fallback. N=3 frags. DONE.
- §1.5 REMAINING: §1.5.1 RT alignment (light verify: MAD x4 factor, FDR 1%), §1.5.3 NCE (light verify). MS1 mass-error fit not yet described in a subsubsection (only listed in params).
- §1.5.1 RT Alignment: DONE. Tolerance factor 4->3 (TUNING_IRT_TOL_SIGMA=3; irt_mad=mad(normalize=true)). 'least squares' -> penalized B-spline + outlier removal. Direction reconciled: intro bullet now 'empirical->library' to match fit (rt_to_irt; process_scans_fused.jl:75). 1% FDR confirmed (TUNING_MAX_Q_VALUE).
- §1.5 REMAINING: §1.5.3 NCE (light verify). MS1 mass-error fit only in param bullet (no subsubsection).
- §1.5.3 NCE Alignment: REWRITTEN. Selection metric Scribe->deconvolution gof (PTSearch.jl:435). Model piecewise-linear(500 breakpoint)->per-charge m/z-binned median NCE (fit_binned_median_nce; >=50 precursors/bin, adaptive 10/5/3/2/1 bins; default NCE fallback). Removed Scribe cite + piecewise eqns. DONE.
- §1.5.4 MS1 Mass Error Estimation: NEW subsubsection added (ms1_diagnostic.jl: nearest MS1 scan by RT, match M+0/M+1/M+2 within +-100ppm, median bias + 5*MAD tol; scalar not spline). §1.5 COMPLETE.

## §1.6 (Main Search) — verifying every claim against DEVELOP (1eb5c1886)
- §1.6 retitled First Pass Search -> Main Search; intro reframed (no first/second pass). DONE.
- §1.6.1 PSM Scoring: DONE. probit/Percolator math DROPPED -> LightGBM (2-fold CV, ~49 PRESCORE_FEATURES, MAINSEARCH_LGBM_HP 50iter/lr0.20/depth8; probit low-data fallback). best-PSM = highest-WEIGHT among top-scoring (select_best_per_precursor!), not highest score. PEP via wPAVA (get_PEP!/_weighted_pava), PEP<=0.9 (MAIN_PEP_FILTER_THR) — MATCHES old. No Scribe init. Verified on develop.
- TODO §1.6.2 Aggregation: old max-score/top-M -> log-odds combine top-sqrt(N) (aggregate_prescore_globally!/_logodds_combine); global qvals used for RT-tol fit; all PEP-passing written to fold-split. VERIFY on develop before rewrite.
- TODO §1.6.3 RT Alignment+Tolerance: refit prob>0.9 (recalibrate_rt!), + peptide-composition iRT correction (irt_refinement.jl), RT-binned FWHM/cycle tolerance. VERIFY on develop before rewrite.
