# MainSearch Per-File LightGBM — Feature Reference

The 46 features in `PRESCORE_FEATURES`
(`src/Routines/SearchDIA/SearchMethods/MainSearch/features.jl:189`), grouped
by where and how they're computed. Each row of the input matrix corresponds
to one (scan, precursor) PSM that survived the deconvolution + scoring
pipeline; the LightGBM target is `:target` (1 = library target, 0 = decoy).

## A. Match counts and ion statistics

Computed inline by `apply_main_scoring!` during the fused per-precursor scan
loop (`src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl:apply_main_scoring!`),
written to `MainUnscoredPSM`, then surfaced by `Score!` into
`MainSearchScoredPSM` (`src/Routines/SearchDIA/PSMs/ScoredPSMs.jl:Score!`).

| Feature | Type | Description |
|---|---|---|
| `total_ions` | UInt8 | `b_count + y_count` — M0 (monoisotopic) fragment match count. |
| `total_ions_iso` | UInt8 | `isotope_count` — # of matched M1+ (isotopologue) fragments. |
| `y_count` | UInt8 | # of matched M0 y-ions (also a feature on its own). |

## B. Error and Poisson

| Feature | Type | Description |
|---|---|---|
| `poisson` | Float16 | `log(λᵏ·e⁻ᵏ) − logfac(observed)` where `λ = expected_matches`, `k = total_ions + total_ions_iso`. λ is floored at 1e-20 to avoid `log(0)`. Closer to 0 = more PSM matches than expected by chance (`ScoredPSMs.jl:psm_getPoisson`). |
| `err_norm` | Float16 | `min(2^min(error, 15) / max(total_ions, 1), 6e4)` — exponentiation of `log2(sum |ppm_err|)` per match. Penalizes large ppm-error sums and rewards more matches (`features.jl:48`). |

`error` in the unscored struct accumulates `abs(ppm_err)` per (fragment,
isotope) match. `Score!` stores `error_log2 = log2(max(error, 1e-20))` into
the scored PSM, then `err_norm` reverses to a normalized scale.

## C. Spectral fit (from deconvolution residuals)

Computed by `getDistanceMetrics`
(`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl`) on the dense
residual vector `r = x − Hs·w` (observed intensities minus reconstructed).
Each is stored in `SpectralScoresMainSearch` and copied into the scored PSM.

| Feature | Type | Description |
|---|---|---|
| `gof` | Float16 | Goodness of fit: `−log2(Σ|r| / Σ x̂ + 1e-10)`. Larger = smaller relative residual. Sum over all peaks. |
| `max_matched_residual` | Float16 | `−log2(max|r_matched| / Σ x̂_matched + 1e-10)`. Max residual *only* over peaks the precursor's theoretical contributes to. |
| `max_unmatched_residual` | Float16 | `−log2(max|r| / Σ x̂ + 1e-10)`. Max residual over *all* peaks (catches contamination at unmatched m/z). |
| `fitted_manhattan_distance` | Float16 | `−log2(Σ|x−x̂| / Σx + 1e-10)`. L1 fit metric on the reconstructed spectrum vs. observed. |
| `fitted_hellinger` | Float16 | `−log2(max(Hellinger², 1e-10))`. Hellinger² = ½·Σ(√x̂ − √x)². Probabilistic distance between observed and reconstructed spectra after L1 normalization. |
| `weight` | Float32 | Deconvolution weight `w[col]` for this precursor at this scan (the PoissonMM solver output). Negligible weights are zeroed pre-scoring. |

## D. Library / peptide properties

Added by `add_features!` (`features.jl:76-181`). Lookups into the spectral
library indexed by `precursor_idx`.

| Feature | Type | Description |
|---|---|---|
| `irt_error` | Float32 | `|irt_obs − irt_pred|` where `irt_obs = rt_to_irt_interp(rt)` and `irt_pred = library_irt[precursor_idx]`. |
| `missed_cleavage` | UInt8 | `library.missed_cleavages[precursor_idx]`. |
| `Mox` | UInt8 | Count of `Unimod:35` (oxidized methionine) in the structural-modifications string. Cached per-precursor on first use. |
| `spectrum_peak_count` | Float16 | `length(masses[scan_idx])` — total # of peaks in the MS2 spectrum (proxy for spectral complexity / TIC density). |
| `sequence_length` | UInt8 | `library.length[precursor_idx]`. |

## E. Per-scan competition

`add_scan_competition_features!` (`features.jl:708-739`). Groups PSMs by
`scan_idx`; within each scan, ranks precursor candidates by deconvolution
weight.

| Feature | Type | Description |
|---|---|---|
| `weight_ratio_at_scan` | Float32 | `weight / max(weight)` over all PSMs at the same scan. 1.0 = this precursor has the highest weight in the scan. Single-PSM scans get 1.0. |
| `weight_rank_at_scan` | UInt16 | Rank (1 = highest weight) of this PSM among PSMs at the same scan. |

## F. Per-precursor neighborhood (3-scan window)

`add_neighborhood_features!` (`features.jl:647-697`). For each PSM, looks at
the same precursor's PSMs at the *nearest earlier* and *nearest later* MS2
scans (when available — that's the "3-scan" window: self + ≤2 neighbors).
Computes "best" of three quality metrics in the window and the iRT distance
from this scan to the scan that contributed the best.

| Feature | Type | Description |
|---|---|---|
| `best_gof_3scan` | Float32 | `max(gof)` over self + ≤2 chronological neighbors. |
| `best_manhattan_3scan` | Float32 | `min(fitted_manhattan_distance)` over the same window (smaller = better). |
| `best_max_residual_3scan` | Float32 | `min(max_matched_residual)` over the same window. |
| `irt_dist_best_gof_3scan` | Float32 | `|irt_obs(this) − irt_obs(scan_with_best_gof)|`. |
| `irt_dist_best_manhattan_3scan` | Float32 | Analogous for the manhattan-winning scan. |
| `irt_dist_best_max_residual_3scan` | Float32 | Analogous for the max-residual-winning scan. |

Rationale: real precursors elute over a few cycles, so the best neighbor
should be close in iRT. Random hits don't form coherent chromatographic
neighborhoods.

## G. Per-precursor apex distance

`add_apex_distance_feature!` (`features.jl:602-624`). For each precursor,
find the scan with the maximum deconvolution weight (its chromatographic
apex), then report each PSM's iRT distance to that apex.

| Feature | Type | Description |
|---|---|---|
| `irt_dist_to_weight_apex` | Float32 | `|irt_obs(this) − irt_obs(scan_with_max_weight_for_this_precursor)|`. Apex PSMs get 0; far-from-apex scattered hits get large values. |

## H. MS1 point-lookup features

`add_ms1_features!` (`features.jl:217-348`). For each PSM, find the nearest
MS1 scan by RT, then search ±`PIONEER_MS1_PPM_TOL` (default 10 ppm) for the
M0 (monoisotopic) and M+1 peaks of the precursor at its library m/z.

| Feature | Type | Description |
|---|---|---|
| `ms1_m0_intensity` | Float32 | Intensity of the best-matched M0 peak in the nearest MS1 scan, or 0 if no match. |
| `ms1_m1_intensity` | Float32 | Same for M+1 (m/z + NEUTRON / charge). |
| `ms1_m0_mass_err_ppm` | Float32 | `|m0_observed_mz − m0_theoretical_mz| / m0_theoretical_mz · 1e6`, or 0 if no match. |

## I. Per-precursor MS1 chromatogram correlations

`_add_ms1_chromatogram_features!` (`features.jl:371-448`). Groups PSMs by
precursor and treats the per-MS2-scan `ms1_m0_intensity`, `ms1_m1_intensity`,
and deconvolution `weight` columns as parallel chromatograms (sampled at
each MS2 scan where this precursor has a PSM). Pearson correlations require
≥ 2 points; smaller groups get 0.

| Feature | Type | Description |
|---|---|---|
| `ms1_corr_weight_m0` | Float32 | `Pearson(weight_chrom, M0_chrom)` — MS2 deconvolution weight should track MS1 M0 intensity through the precursor's elution. |
| `ms1_corr_m0_m1` | Float32 | `Pearson(M0_chrom, M+1_chrom)` — real isotopologues co-elute; chimeric noise doesn't. |
| `ms1_apex_offset_irt` | Float32 | `|irt_obs(this) − irt_obs(argmax(M0_chrom))|` — distance from this PSM's scan to the MS1 M0 apex. |
| `ms1_weight_apex_to_m0_apex_irt` | Float32 | `|irt_obs(argmax(weight_chrom)) − irt_obs(argmax(M0_chrom))|` — disagreement between MS2 apex and MS1 apex. Real: ~0; chimeric: large. |

## J. Per-precursor MS2 fragment chromatograms

`_add_fragment_chromatogram_features!` (`features.jl:469-589`). The top-6
M0 fragments per precursor have their per-(scan, precursor) intensities
captured by `apply_main_scoring!` (fields `frag1_int..frag6_int` on
`MainUnscoredPSM`) and surfaced into `MainSearchScoredPSM`. Each precursor
thus has 6 across-scan chromatograms. Fragments with `max(F[r]) = 0` are
skipped from correlation computations.

### Per-rank raw intensities

| Feature | Type | Description |
|---|---|---|
| `frag1_int` … `frag6_int` | Float32 each | Summed intensity of rank-1 through rank-6 M0 fragments at this scan. Rank is the library's `getRank(frag)`. |

### Correlation features

| Feature | Type | Description |
|---|---|---|
| `frag_corr_top1_top2` | Float32 | `Pearson(frag1_chrom, frag2_chrom)`. |
| `frag_corr_top1_top3` | Float32 | `Pearson(frag1_chrom, frag3_chrom)`. |
| `frag_corr_top1_weight` | Float32 | `Pearson(frag1_chrom, weight_chrom)`. |
| `frag_corr_mean_pairwise` | Float32 | Mean Pearson across all 15 pairs `(frag_i, frag_j), i<j ∈ 1..6` with signal. |
| `frag_corr_min_pairwise` | Float32 | Min Pearson over the same pairs. A single contaminated fragment drags this low even if the rest co-elute. |
| `frag_corr_top3_weight` | Float32 | Mean `Pearson(frag_r, weight)` for `r ∈ 1..3` with signal. |
| `frag_apex_dispersion_irt` | Float32 | Std-dev of `argmax(frag_r_chrom)` iRT values across the 6 fragments with signal. Real precursors: all 6 fragment apexes cluster; chimeric: scattered. |
| `n_correlated_fragments` | UInt8 | Count of `r ∈ 1..6` such that `Pearson(frag_r_chrom, weight_chrom) > 0.7`. |

## How the LightGBM is trained

1. The full feature matrix `X_all = feature_matrix(psms, PRESCORE_FEATURES)`
   has one row per (scan, precursor) PSM and 46 columns.
2. PSMs are split into two CV folds (column `:cv_fold`, deterministic per
   precursor via `getCvFold(precursors, prec_idx)`).
3. Per fold pair `(train_idx, test_idx)`:
   - `_prepare_labels` converts `:target` Bools to `Int{0,1}`.
   - LightGBM trains with `SHARED_LGBM_HP`
     (`MainSearch/scoring.jl:25` — 200 iter × 0.10 lr, max_depth=8,
     num_leaves=63, min_data_in_leaf=300, feature_fraction=0.8,
     bagging_fraction=0.8, bagging_freq=1, lambda_l1=1.0, lambda_l2=1.0).
   - Predict on the held-out fold and write back to `:lgbm_score`.
4. `select_best_per_precursor!` keeps one row per precursor by either
   highest score or, if ≥ 4 PSMs and ≥ p75 score, highest weight; computes
   per-precursor aggregates (`smoothness`, `irt_fwhm`, `n_above_hm`,
   `rt_fwhm`, `num_scans`, `max_weight`, `max_gof`,
   `max_fitted_manhattan_distance`).
5. Optional per-file pair competition + per-file PEP filter (≤ 0.9 default).
6. Best PSMs written to `temp_data/main_search_psms/<file>_fold{0,1}.arrow`
   for the experiment-wide PrecursorScoringSearch LightGBM to read.

The per-file LightGBM's probabilities (column `:lgbm_score` → renamed
`:lgbm_prob` after best-per-precursor) are also log-odds-combined into the
global prescore aggregation that ScoringSearch later consumes.
