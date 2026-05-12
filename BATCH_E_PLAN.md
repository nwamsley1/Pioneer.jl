# Batch E Plan — alphaDIA-inspired features

Branch: `feature/ms1-phase1-m1frag`. Builds on the A/B/C/D peak-shape work
completed earlier today (see `PEAK_SHAPE_FEATURE_AB.md` for those results).

## A/B/C/D recap (already done)

Cumulative 2-file Olsen result with all 4 batches active:

| Stage | Final IDs | f1 q≤.01 | f2 q≤.01 | PGs q≤.01 |
|---|---:|---:|---:|---:|
| Baseline (m1frag, no peak-shape features) | 88,875 | 41,213 | 41,183 | 12,386 |
| + Batch A (edge-distance) | 88,280 | 41,304 | 40,768 | 12,334 |
| + Batch B (peak-shape stats) | 88,612 | 41,140 | 40,913 | 12,311 |
| + Batch C (Gaussian fit) | 88,217 | 41,101 | 40,825 | 12,260 |
| + Batch D (DIA-NN ports) | 88,798 | 40,951 | 40,800 | **12,506** |

**Net A→D:** +120 PGs (first PG-positive batch this dev cycle); final IDs
flat within LGBM noise (~±400 on 2-file).

### Keepers (gain > 1500, ranked)

| Feature | Batch | Gain | Note |
|---|---|---:|---|
| `weight_chrom_skewness` | B | 2,628 | Top of batch B |
| `apex_to_edge_weight_ratio` | B | 2,243 | Clean apex isolation |
| `weight_chrom_kurtosis` | B | 2,281 | Spike vs broad signal |
| `weight_chrom_gaussian_r2` | C | 1,721 | DIA-NN-style log-quadratic fit |
| `weight_chrom_gaussian_apex_irt` | C | 1,530 | |

### Medium keepers (400–1500)

`weight_chrom_gaussian_sigma` (C, 1,445), `weight_apex_relative_pos` (B, 1,176),
`gof_minus_max_gof_precursor` (D, 915), `shape_0` (D, 795), `shape_p1` (D, 597),
`shape_p2` (D, 583), `shape_m2` (D, 578), `shape_m1` (D, 455),
`dist_to_relative_center` (A, 302), `relative_position` (A, 225).

### Dropped (zero LGBM splits — 6 features)

A1–A4, B5, B7, D7 — all redundant with stronger keepers above.

---

## Batch E — full plan

For each feature: implementation outline, cost, and where in the code it
goes. **Test order: tier 1 (E2, E7, E10, E11, E12) → tier 2 (E6, E1, E14, E3) → defer E9, E15.**

### E7 — `top3_ms2_mass_error_mean`

**What:** Mean absolute ppm error of the top-3 fragments by *library-predicted intensity* (= ranks 1, 2, 3 in our scheme).

**Why:** Existing `err_norm` averages over all matched fragments (including noisy low-intensity matches that drag the mean). Top-3 alone should be cleaner — alphaDIA uses this as its primary mass-error feature.

**Implementation:**
1. **Capture per-rank ppm errors in `apply_main_scoring!`:** Currently captures `frag1..6_int`. Add `frag1..6_ppm_err::Float32` (6 new fields). When rank ≤ 6 and M0, accumulate `ppm_err` (signed or absolute? — start with signed to preserve direction info).
2. **Propagate** through `MainSearchScoredPSM` + `Score!` (mirror of how `frag1_int..frag6_int` is currently propagated).
3. **Feature compute:** `top3_ms2_mass_error_mean = mean(|frag1_ppm_err|, |frag2_ppm_err|, |frag3_ppm_err|)` over matched (non-zero) ranks.

**Cost:** moderate — requires widening `MainUnscoredPSM`/`MainSearchScoredPSM` by 6 fields each. ~30 LOC across 3 files. Same shape as last week's M+1 fragment-intensity capture work.

---

### E2 — `pred_obs_max_intensity_spectral_contrast` (you suggested might beat E1)

**What:** For each precursor, take the **max** of each fragment's intensity across all the precursor's scans. Then compute spectral_contrast (cosine similarity) or scribe between those 6 max-intensities and the library's 6 predicted intensities.

**Why:** Direct test of "do the fragments that the library expects to be strongest actually have the strongest peaks?" Real precursors → similarity ≈ 1; chimeric → low. You noted this might be more informative than E1 (which uses areas).

**Implementation:**
1. **Predicted intensities lookup:** the library has `LibraryFragmentLookup` with per-(precursor, rank) predicted intensities. Need to look up the top-6 predicted intensities for each precursor. **Audit needed:** confirm we can access these from `_add_fragment_chromatogram_features!` (we already pass `search_context` upstream so it should be reachable, but specifics depend on the lookup API).
2. **Observed max:** during the per-precursor pass in `_add_fragment_chromatogram_features!`, we already iterate F[r] for r in 1..6 and have access to all scans. Add `obs_max[r] = maximum(F[r])` for r in 1..6.
3. **Normalize both vectors** (predicted and observed_max) to sum to 1 or to unit L2 norm.
4. **Compute spectral_contrast:**
   ```
   spectral_contrast = dot(pred_norm, obs_norm) / (norm(pred_norm) * norm(obs_norm))
   ```
   (cosine; range -1 to 1, 1 = identical relative intensities).
5. **Compute scribe score** (alternative similarity, more robust):
   ```
   scribe = -log2(sum(|pred_norm - obs_norm|) / 2 + ε)
   ```
   higher = better.
6. **Add features:** `pred_obs_max_spectral_contrast`, `pred_obs_max_scribe` (test both, drop loser).

**Cost:** medium. Library lookup is the only unknown; if accessing the lookup table from features.jl is straightforward, this is ~40 LOC.

---

### E1 — `pred_obs_area_spectral_contrast` (your detailed clarification)

**What:** Same as E2 but use **fragment area** (sum or SG-smoothed integral of intensity across scans in an N-scan window) instead of max.

**Why:** Area is a more stable estimator of total fragment abundance than max (max is sensitive to single-scan spikes). Likely correlated with E2 but more robust.

**Implementation:**
1. **Fragment area:** for each rank r ∈ 1..6, sum `F[r]` across all the precursor's scans (= sum of frag_r_int across this precursor's PSMs). That's the "area" — no smoothing needed because we're already on the deconvolved weight-axis.
   - **Optional SG smoothing:** apply a 5-scan Savitzky-Golay quadratic smoother to F[r] before summing, to reduce per-scan noise. Adds 1 line if we precompute SG weights.
2. **Predicted intensities:** same library lookup as E2.
3. **Normalize + spectral_contrast / scribe:** same as E2.
4. **Add features:** `pred_obs_area_spectral_contrast`, `pred_obs_area_scribe`.

**Cost:** if E2 lands first, this is ~10 additional LOC (same library lookup, different aggregation).

**Decision rule:** test E2 (max) and E1 (area) side-by-side first time, then drop whichever loses head-to-head (per the no-redundancy preference we established last night).

---

### E11 + E12 — `top3_b_ion_correlation` / `top3_y_ion_correlation`

**What:** For each precursor:
1. Identify which of the 6 top-rank fragments are b-ions and which are y-ions.
2. Compute the **median chromatogram** across the 6 fragment chromatograms (per scan, take median of F[1..6][scan]).
3. For each fragment, compute `Pearson(F[r], median_chrom)`.
4. Take the mean of those correlations *restricted to top-3 b-ions* (i.e., the top-3 ranks that are b-ions) and *restricted to top-3 y-ions* (independently).

**Why:** Real peptides have b- AND y-ion families that both coelute. Chimeric false positives often have one ion family co-eluting and the other scattered (because the contaminant precursor only contributes one ion class strongly to the matched ranks). This per-ion-type breakdown of E10 catches that asymmetry.

**Implementation:**
1. **Capture per-rank ion-type bits:** in `apply_main_scoring!`, currently the `frag` object has `isB()`, `isY()` predicates. Add 6 bit-fields (or one UInt8 with 6 bits) recording which of ranks 1..6 are b vs y. Promote through `MainSearchScoredPSM` like the other frag captures.
   - Alternative: look up ion-type from the library at feature-compute time. Less invasive — requires the library accessor. **Audit this first** before adding struct fields.
2. **Median chromatogram:** during the per-precursor loop in `_add_fragment_chromatogram_features!`, for each scan in the precursor's PSM list, compute median(F[1][scan], F[2][scan], ..., F[6][scan]) restricted to ranks with signal. Standard median of up-to-6 values.
3. **Per-fragment correlation to median:** `c_to_median[r] = Pearson(F[r], median_chrom)` for r in 1..6 with signal.
4. **Top-3 by ion type:** take ranks {1..6} ∩ {b-ions}, sort by rank, take first 3. Average c_to_median over those. Same for y. (Note: "top-3" means ranks 1, 2, 3 of each ion type, *not* the rank-1/2/3 overall — alphaDIA's convention is by library predicted intensity within each ion class.)
5. **Add features:** `top3_b_ion_corr_median`, `top3_y_ion_corr_median`, `n_b_ions_top6`, `n_y_ions_top6`.

**Cost:** moderate. The ion-type capture is the main scope question. If library accessor route works, ~30 LOC pure in features.jl.

---

### E10 — `mean_fragment_to_median_corr`

**What:** Mean of `Pearson(F[r], median_chrom)` over all ranks r ∈ 1..6 with signal. Same machinery as E11/E12 but no ion-type partition.

**Why:** Different reference than `frag_corr_best_*` (consensus-best) and `frag_corr_mean_pairwise_spearman` (15 pairs averaged). Median is robust to one or two contaminated fragments dragging the average. Quick comparison test for "which reference (best/median/all-pairs) gives the best chromatographic-quality signal?"

**Implementation:** Tiny additional bookkeeping on top of E11/E12 — compute the mean while you have the per-rank correlations.

**Cost:** ~5 LOC on top of E11/E12.

---

### E6 — `log_b_y_intensity_ratio` (and a M0+M1 variant)

**What:** Two variants of the b/y log-intensity ratio:
1. **`log_by_ratio_m0`** = `log(Σ b_intensity_m0 + 1) − log(Σ y_intensity_m0 + 1)`. M0 only.
2. **`log_by_ratio_m01`** = `log(Σ b_intensity_m0 + Σ b_intensity_m1 + 1) − log(Σ y_intensity_m0 + Σ y_intensity_m1 + 1)`. Includes M+1 isotopes.

**Why (your point on narrow-window DIA):** For wide-isolation DIA the M+1 contribution is largely separated from the target; for narrow-window DIA the isotopes co-occupy the isolation window more, so including them changes the ratio. Test which variant the LGBM prefers on each dataset; if neither hurts, default to M01.

**Implementation:**
1. We already have `b_int` and `y_int` in `MainUnscoredPSM` (M0). Just plumb through to `MainSearchScoredPSM` if not already (already there as part of the existing features).
2. For M+1: we'd need `b_int_m1` and `y_int_m1` additional accumulators in `MainUnscoredPSM` (mirror of `b_int`/`y_int` but only when iso_idx==1). 2 new fields each in Unscored and Scored.
3. **Feature compute:** trivial — one `log(Σ + 1) − log(Σ + 1)` line per variant.

**Cost:** small if we only do M0 variant (the data is already there — just need to compute log-ratio). Moderate if we also do M01 (struct widening, 2 fields each in Unscored/Scored).

**Recommendation:** start with just M0 to confirm signal, then add M01 if M0 is informative.

---

### E14 — `delta_frame_peak_center` (clarification)

**What:** Per precursor:
1. Find the apex iRT of each top-6 fragment chromatogram (argmax of F[r] for r in 1..6 with signal). Call them `apex_irt_r`.
2. Take the **median** of those 6 apex iRTs: `median_apex_irt`.
3. Compute the **center** of the precursor's scan range: `center_irt = (irt_first + irt_last) / 2` (midpoint of the iRT span of all scans for this precursor).
4. `delta_frame_peak_center = median_apex_irt - center_irt`. Signed (positive = fragments peaked late in the scan window, negative = early).

**Why:** Real precursors are searched around the predicted iRT, so the fragment apex should land near the center of the scan range we considered. Chimeric/false hits often have the fragment apex drift toward one end of the window because the contaminant elutes at a different iRT. The signed nature lets the LGBM distinguish "drifted early" from "drifted late".

**Implementation:** simple per-precursor pass in `_add_fragment_chromatogram_features!`. We already compute apex_irts in the existing code; just take median + center.

**Cost:** ~10 LOC.

---

### E3 — `matched_ratio` (your reference to old code)

**What:** Resurrect the historical `matched_ratio` feature. Standard definition:
- `matched_ratio = log2((Σ observed_intensity_of_matched_fragments) / (Σ library_predicted_intensity_of_all_fragments))`
- Higher (less negative) = more of the predicted ion budget actually matched in the spectrum.

**Why:** Captures "completeness of match" — fragmentation should produce the full predicted ladder; missing predicted ions = potential chimeric or poor match.

**Implementation:**
1. Confirm presence in git history (`git log -S "matched_ratio" --all`) — your hint says it's there from last year. Likely lives in old `apply_main_scoring!` or `Score!` paths.
2. If found, port the formula. If not found, implement from scratch: needs `library_total_predicted_intensity` per precursor (one number from the library) + observed matched-intensity sum (we have `b_int + y_int` in unscored).
3. Add as a single feature.

**Cost:** small if we resurrect from history; otherwise ~15 LOC plus a library accessor.

---

### E9 — `n_b_y_overlapping` (deferred per your skepticism)

**What:** Count of fragments where a b-ion at sequence-position `i` and a y-ion at sequence-position `j` both matched with `i + j > peptide_length` (= they cover overlapping residues).

**Why my reasoning was shaky:** I claimed this catches contamination, but you're right to push back. In reality, peptides DO produce overlapping b/y from different fragmentation events (single peptide, many fragment paths, all valid). The alphaDIA code masks "overlapping" but the rationale isn't well-explained in their paper either — could be more about *abundance* of overlapping matches (a chimeric pair would have *more* overlapping than a single peptide's natural ladder).

**Decision:** drop E9 from the plan. Not enough mechanistic justification to spend implementation budget.

---

### E15 — `template_frame_correlation` (deferred for complexity)

**What:** Compute a *predicted* (template) chromatogram shape per fragment — typically a Gaussian centered at the predicted iRT with σ from the per-file FWHM model. Correlate observed fragment chromatograms against this template.

**Why deferred:** Requires building/maintaining a per-precursor template. AlphaDIA does this with their search-window math; for us it'd require accessing the per-file iRT FWHM model and building synthetic Gaussians. ~80 LOC + correctness risk. Revisit if Batch E tier 1 results show big lifts that suggest we're still leaving signal on the table.

---

## Test order

Same approach as A/B/C/D: implement → run → record LGBM gains → keep features with gain > ~300 (rough cutoff from prior work). Each batch ≈ 3-min test run on 2-file Olsen.

**Tier 1 (next session):**
1. **E7** — top3_ms2_mass_error_mean. Independent of fragment-area work; can land first.
2. **E10 + E11 + E12** — fragment-to-median correlations (combined into one batch since they share the per-precursor pass).
3. **E2** — observed-max vs predicted spectral_contrast / scribe. Independent of struct changes (just needs library accessor).

**Tier 2 (after tier 1 lands):**
4. **E1** — area variant of E2. Tested side-by-side with E2 (cheaper since the lookup is already done).
5. **E6** — log_by_ratio (M0 first, M01 if M0 wins).
6. **E14** — delta_frame_peak_center. Standalone, ~10 LOC.
7. **E3** — matched_ratio (if findable in git history, otherwise from scratch).

**Deferred:**
- E9 (no clear mechanism)
- E15 (too much scope vs expected lift)

## Open audits before starting

Two things to confirm before tier 1 implementation:

1. **Library predicted-intensity accessor** — does the `LibraryFragmentLookup` / `SpectralLibrary` API expose per-(precursor, rank) predicted intensities that we can pull from `_add_fragment_chromatogram_features!`? If yes, E2/E1/E10/E11/E12 all become much cheaper. If no, we'd need to capture predicted intensities at scan time (struct widening — ~6 more fields).

2. **Per-rank ion-type bits in `MainUnscoredPSM`** — for E11/E12, simpler to capture during `apply_main_scoring!` (1 bit per rank → UInt8) than to look up from library every PSM. Worth a small struct widening if the alternative is per-PSM library access.

Both are 1-hour audits. I'll do them before implementing tier 1.

---

## Decision: combine or isolate?

So far we've been doing additive batches (A → A+B → A+B+C → A+B+C+D). This makes per-feature attribution cleaner via LGBM gain but accumulates risk. For Batch E I'd suggest:

- **Tier 1 features (3 sub-batches above)** stay additive — each builds on the prior.
- **Tier 2** is also additive on top of tier 1.

If at any point a single sub-batch tanks IDs by >5% (vs LGBM noise), revert that sub-batch's commits and continue forward without it. Each sub-batch gets its own git commit so reverts are mechanical.
