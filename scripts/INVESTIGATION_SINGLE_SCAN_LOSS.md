# Investigation: Single-Scan Precursor Loss

## Problem Statement

Precursors with excellent first-pass scores (prob=1.0, high scribe, q=0) and excellent second-pass scores (prec_prob≈0.9999) are absent from the final output. This was discovered during FRAGCORR chromatogram browsing.

**Scale of the problem:**
- ~14,000-15,000 target precursors lost per file (~23-25% of scoring-passing targets)
- Of 300 first-pass targets with prob=1.0 missing from global output, **239 (80%) had scan_count=1**

**Exemplar precursor:**
- `precursor_idx=6374682` (DQELYFFHELSPGSCFFLPK 3+, protein P26639)
- trace_prob=0.9997-0.9999 in 5/6 files after LightGBM, 0.936 in file 6
- global_qval=8.1e-5, qval=5.6e-5 to 5.5e-3 across files
- **Survives all of ScoringSearch in all 6 files**
- **Present in final output for only 2/6 files** (E5H50Y45_2 and E5H50Y45_4)

## Root Cause: IntegrateChromatogramSearch

### The loss is NOT in ScoringSearch

Diagnostic tracing (run 2026-03-03) confirmed precursor_idx=6374682 survives every ScoringSearch stage:

| Stage | Result |
|-------|--------|
| Step 1 (LightGBM) | trace_prob 0.936–0.9999 in all 6 files |
| Step 2 (aggregation) | prec_prob ≈ trace_prob (single trace per file) |
| Step 4 (best_trace filter) | Survived all 6 files, isotopes_captured=(0,5) |
| Steps 5-10 (dual q-value) | global_qval=8.1e-5, qval 5.6e-5 to 5.5e-3 — all pass 1% threshold |
| Step 13 (protein inference) | Assigned to P26639 in all 6 files |
| Step 24 (protein scores) | pg_qval=0.0, global_pg_qval=0.0 |

**Hypotheses 1-5 from initial speculation are all RULED OUT for this precursor.**

### The loss IS in IntegrateChromatogramSearch

`IntegrateChromatogramsSearch/utils.jl:701-702` (`process_final_psms!`):
```julia
filter!(row -> !isnan(row.peak_area::Float32), psms)
filter!(row -> row.peak_area::Float32 > 0.0, psms)
```

Precursors with `peak_area == 0.0` are silently removed. `peak_area` is initialized to 0 and only set to a positive value if `integrate_precursors` successfully integrates the chromatographic peak.

### Why integrate_precursors returns peak_area=0

`integrate_precursors` (utils.jl:82-144) has three `continue` bail-outs that leave peak_area at 0:

1. **Line 88/91**: Precursor not found in chromatogram group keys → `continue`
   - `extract_chromatograms` scans MS2 spectra and matches fragments to build chromatogram points
   - If no fragments are matched for this precursor across all scans, no group exists

2. **Line 100**: All chromatogram intensities are zero → `continue`
   - Deconvolution can assign zero weight to a precursor even when fragments match

3. **Line 103-104**: Apex scan_idx not found in trimmed chromatogram → `continue`
   - After trimming to first-to-last positive intensity, the original best scan may fall outside the window
   - For single-scan precursors, if the one positive scan gets a zero weight from deconvolution, nothing survives trimming

### Scale of Integration Loss

Per-file target losses (ScoringSearch output → final output):

| File | Targets passing scoring | In final output | Lost | % Lost |
|------|------------------------|-----------------|------|--------|
| E45H50Y5_2 | 62,623 | 48,107 | 14,516 | 23.2% |
| E45H50Y5_3 | 62,463 | 48,398 | 14,065 | 22.5% |
| E45H50Y5_4 | 61,203 | 47,946 | 13,257 | 21.7% |
| E5H50Y45_1 | 61,623 | 46,005 | 15,618 | 25.3% |
| E5H50Y45_2 | 60,815 | 45,954 | 14,861 | 24.4% |
| E5H50Y45_4 | 60,825 | 45,925 | 14,900 | 24.5% |

**~14,000 targets per file (22-25%) pass ScoringSearch but are lost to failed chromatogram integration.**

Note: not all of these are single-scan. Some may be multi-scan precursors where deconvolution still produces zero peak area. But the 80% single-scan rate from the first-pass analysis suggests the majority are.

In the final output, 5,677 of 282,335 rows (2.0%) have `points_integrated=1` — these are the single-scan precursors that survived integration (barely).

## Why Single-Scan Precursors Fail Integration

For precursor_idx=6374682, the 2 files that survived have `points_integrated=1` and `new_best_scan=57885`. The chromatogram extraction found exactly one scan with positive deconvolved intensity.

For the 4 files where it failed, the likely sequence:
1. `build_chromatograms` scans through MS2 spectra matching fragments
2. For a single-scan peptide, only 1 (or very few) scans have matching fragments
3. Deconvolution competes this precursor's signal against other co-eluting precursors
4. In some files, the competition gives this precursor zero weight (other precursors explain the fragment ions better)
5. With zero weight, there's no positive intensity → bail-out at line 100
6. Or the apex scan gets trimmed → bail-out at line 103-104
7. peak_area stays 0 → filtered out by `process_final_psms!`

## Recommendations

### Short-term: Preserve zero-integration precursors
Instead of discarding peak_area=0 rows, keep them with `peak_area=missing` or `peak_area=0`. This preserves the identification (which passed stringent FDR filtering) even when quantification fails. Downstream tools (MaxLFQ) already handle missing values.

### Medium-term: Report integration failure separately
Add a column like `integration_status` (`:integrated`, `:no_chromatogram`, `:no_apex`, `:zero_weight`) so users can see why quantification failed without losing the identification.

### Long-term: Improve single-scan integration
- Consider using the identified scan's raw intensity as a fallback when deconvolution produces zero weight
- For single-scan precursors, the "chromatographic peak" is a single point — the integration should just use the deconvolved weight from that scan directly rather than requiring a multi-point chromatogram

## Diagnostic Tracing

Diagnostic tracing was added to `ScoringSearch.jl` and `scoring_interface.jl`. Set `TRACE_PRECURSOR_IDX` to a non-zero UInt32 in `ScoringSearch.jl` to trace a specific precursor through the pipeline. Set to `UInt32(0)` to disable (zero-cost).

Checkpoints: After Step 1 (LightGBM), Step 2 (aggregation), Step 4 (best_trace), Steps 5-10 (q-value filter), Step 13 (protein inference). Also in `build_precursor_global_prob_dicts` and `build_global_qval_dict_from_scores`.
