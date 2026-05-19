# Fused Matching and Fragment Index Findings

Updated: 2026-05-08 15:21:35 CDT

## Current Findings

### Fused fragment matching

1. Potential real bug: fused per-fragment peak search advances the next fragment's lower peak bound to the previous fragment's chosen `best_peak`.
   - Location: `src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl`, inside `run_fused!`.
   - Why it matters: old `matchPeaks!` did not advance the empirical peak pointer after a matched transition, so a later theoretical fragment could still match an earlier peak if its window overlapped. The fused code can miss that case when fragment A chooses a later peak and fragment B's valid/best peak is before A's chosen peak.
   - Current assessment: likely worth fixing or at least writing a regression test for close fragment m/z values with overlapping windows.

2. Asymmetric mass tolerance issue is currently latent only.
   - `conservative_half_width` uses the right-side bound and applies it symmetrically.
   - Current assessment: not a current production issue if all active mass error models are symmetric; it is still an API footgun if asymmetric `SimpleMassErrorModel` values are introduced.

3. Corrected-m/z sorting assumption is probably safe for current symmetric/simple models, but not for intensity-dependent corrections.
   - The fused matcher binary-searches `scan_corrected_mz`, which is filled in raw scan order.
   - Current assessment: safe only if correction preserves monotonic order. If `IntensityMassErrorModel` is used, intensity-dependent bias can reorder adjacent peaks and invalidate binary search / linear sweep assumptions.

4. Duplicate precursor IDs within one scan are unlikely in the main fragment-index path.
   - The partitioned fragment index uses a per-partition counter and emits each local precursor once per scan per partition.
   - Current assessment: main path should not duplicate unless a precursor can exist in multiple overlapping partitions. PSM-derived tuning paths are less structurally protected and may still be worth deduping or asserting.

## Fragment Index Investigation Findings

### Precursor m/z partition query

1. High-confidence miss risk: empty precursor m/z partitions can break `get_partition_range`.
   - Build path: `src/structs/SpectralLibrary/PartitionedFragmentIndex/build.jl`.
   - Query path: `src/structs/SpectralLibrary/PartitionedFragmentIndex/types.jl`, `get_partition_range`.
   - Why it matters: `build_partitioned_index_from_lib` creates initial precursor m/z bins across the whole min-to-max range, then calls `_split_balanced!` on every bin. Empty bins are preserved as final partitions. Their bounds are recorded as `(Inf, -Inf)`.
   - Failure mode: `get_partition_range` assumes `partition_bounds` are sorted and binary-searchable by min/max. An empty partition between populated partitions makes the max sequence and min sequence non-monotone, so a scan's precursor m/z query can return an incomplete partition range. Because `_build_partition_scan_mapping` only schedules scans for that returned range, real candidate precursors in skipped partitions are never scored.
   - Example shape: `[(400, 405), (Inf, -Inf), (415, 420)]`. A scan spanning 399-406 can fail to include partition 1, and a scan spanning later m/z can fail in the other direction depending on where the binary search lands.
   - Test gap: existing tests cover ordinary monotone partition ranges and construction of an empty `LocalPartition`, but not an empty partition between populated bounds.
   - Suggested fix: do not push empty precursor partitions into `final_partition_pids`, or make `get_partition_range` use a separate compact, monotone list of non-empty partition bounds. The first option is simpler and likely faster.

2. Conditional miss risk: scan precursor m/z expansion assumes charge 2.
   - Location: `src/structs/SpectralLibrary/PartitionedFragmentIndex/search.jl`, `_precompute_scan_properties`.
   - Current code widens scan precursor m/z bounds by `NEUTRON * iso_bound / 2`.
   - Downstream fused filtering uses `NEUTRON * iso_bound / prec_charge` in `quad_window_with_iso_bounds`.
   - Consequence: for charge 2 this matches; for charge >2 the partition query is wider than necessary; for charge 1 the partition query is narrower than the final filter and can miss charge-1 precursors before scoring.
   - Current assessment: likely not active if libraries are built with default `min_charge = 2`, but it is a real bug for charge-1 libraries or workflows.

### RT bin query

3. Boundary miss risk: the RT-bin sweep excludes bins whose low bound equals the scan iRT upper bound.
   - Location: `src/structs/SpectralLibrary/PartitionedFragmentIndex/search.jl`, `_score_partition_hinted!`.
   - Current loop condition is `getLow(rt_bin) < irt_high`.
   - Downstream precise filter is inclusive: `abs(prec_irt - scan_irt) <= irt_tol`.
   - Consequence: candidates exactly on the upper RT boundary can be omitted by the fragment index even though final filtering would accept them.
   - Suggested fix: use `<= irt_high` if inclusive RT-window semantics are intended.

### Fragment m/z bin query

4. Fragment bin overlap search looks structurally correct under monotone peak-window assumptions.
   - Location: `src/structs/SpectralLibrary/PartitionedFragmentIndex/search.jl`, `queryFragmentHinted!`.
   - It finds the first bin with `high >= frag_mz_min` and scores bins until `low > frag_mz_max`, which is the right inclusive overlap criterion.
   - The upper-bound hint is only used to find the first matching bin; the scoring loop can continue past it until bin lows exceed the peak window.

5. Conditional miss risk: intensity-dependent mass correction can invalidate lower-bound advancement.
   - `_score_partition_hinted!` iterates raw scan peaks in m/z order, computes corrected fragment query windows, and advances the next peak's lower fragment-bin guess from the previous first match.
   - This is safe only if corrected/query lower bounds are nondecreasing with raw m/z. It should hold for the current simple symmetric models, but not necessarily for intensity-dependent corrections.

## Deconvolution / Scoring Investigation Findings

### Design matrix assembly

1. The fused CSC row assembly mostly matches the intended old `buildDesignMatrix!` semantics.
   - Matched rows now use raw scan peak indices instead of remapped dense peak rows. This is valid because the solver and distance metrics only require consistent row IDs.
   - Misses are assigned unique synthetic rows after the real peak rows, so unmatched theoretical signal is penalized without colliding with observed peaks.
   - Duplicate same-column/same-peak matches are deduped by summing predicted intensity. The observed `x` value is shared for that peak, so keeping one `x` is fine.

2. Potential attribution bug if the same precursor ID appears more than once in a scan's precursor list.
   - `run_fused!` creates columns sequentially and writes `id_to_col[prec_idx] = this_col`.
   - `Score!` loops over column-positioned `unscored_psms[i]`, but then re-looks up `scores_idx = IDtoCOL[precursor_idx]`.
   - If the same `precursor_idx` is assigned multiple columns, `IDtoCOL` points to the last one. Earlier rows would be scored with the later column's weight and spectral metrics.
   - Current assessment: this should not happen in the main fragment-index path if each scan's precursor list is unique, but there is no assertion. This is worth guarding because it would silently misattribute scores.

3. Column IDs are stored as `UInt16`.
   - `run_fused!` increments `col::UInt16`; `id_to_col` also stores `UInt16`.
   - If a scan ever produces more than 65,535 active columns, column IDs wrap and attribution becomes invalid.
   - Current assessment: probably impossible in normal DIA windows, but a cheap assertion would make this invariant explicit.

### Predicted intensities / weights

4. Clear bug in `FullPrecCapture`: fragment isotope abundance uses precursor charge instead of fragment charge.
   - Location: `src/Routines/SearchDIA/CommonSearchUtils/getFragIsotopes.jl`, `getFragIsotopes!(::FullPrecCapture, ...)`.
   - Current code sets `frag_charge = getPrecCharge(frag)` and evaluates `iso_splines(..., frag_mz * frag_charge)`.
   - It should use `getFragCharge(frag)`. Otherwise fragment isotope distributions are computed at the wrong neutral mass whenever precursor charge differs from fragment charge.
   - Current assessment: latent for current defaults because MainSearch / ParameterTuning / IntegrateChromatograms default to `PartialPrecCapture`, but real if a user selects `FullPrecCapture`.

5. Conditional stale-weight risk for zero-norm active columns.
   - Columns are created when a match is found, even if the predicted intensity for that match is zero or effectively zero.
   - OLS skips zero-norm columns and leaves their warm-start weight unchanged. PMM also cannot update a column with all-zero design entries.
   - If such a column had a positive `precursor_weights` warm start from a previous scan, it could survive `zero_negligible_weights!` and be emitted as a scored PSM.
   - Suggested guard: after assembly or before scoring, force weights to zero for columns whose column norm / sum of predicted intensity is zero.

### Score calculation

6. Current `Score!` only filters on deconvolution weight.
   - Older complex scoring also filtered on y-ion count, total fragment count, top-N evidence, best rank, fitted spectral contrast, and matched ratio.
   - Current `MainSearchScoredPSM` no longer carries several of those old fields, so this may be intentional. But it means any positive-weight column becomes a row, even with weak fragment evidence.
   - Current assessment: not necessarily a correctness bug, but it changes the candidate population and can amplify any weight-attribution or stale-weight issue.

7. Some old spectral features are now hard-zero or absent.
   - `percent_theoretical_ignored` is always zero in the current single-pass `getDistanceMetrics`.
   - `spectral_contrast`, `fitted_spectral_contrast`, `matched_ratio`, `scribe`, `best_rank`, and `topn` are no longer present in the current main scored PSM rows, and downstream feature selection filters missing features out.
   - Current assessment: likely an intentional feature-set reduction for prescoring, but it is a behavioral change from early-February scoring.
