# Fix: Single-Point Chromatogram Integration Bug

## Problem

~14,000 target precursors per file (22-25% of scoring-passing targets) are silently dropped from the final output. These precursors have excellent identification scores (prob=1.0, prec_prob≈0.9999, global_qval≈8e-5) and survive every ScoringSearch stage. They are lost in IntegrateChromatogramSearch because the chromatogram integration code produces NaN peak areas for single-point chromatograms.

## Evidence from Diagnostic Tracing

Precursor_idx=6374682 (DQELYFFHELSPGSCFFLPK 3+, protein P26639) was traced through the entire pipeline:

| File | Chrom Points | Deconvolved Intensities | In Final Output? |
|------|-------------|------------------------|------------------|
| E45H50Y5_2 | **1** | 9.17e9 | **NO** |
| E45H50Y5_3 | **1** | 8.08e9 | **NO** |
| E45H50Y5_4 | **1** | 3.53e9 | **NO** |
| E5H50Y45_1 | **1** | 2.02e9 | **NO** |
| E5H50Y45_2 | 2 | 1.03e8, 1.62e9 | YES (peak_area=2.77e7) |
| E5H50Y45_4 | 2 | 1.31e8, 2.60e8 | YES (peak_area=2.34e6) |

**Pattern: Exactly the files with 1 chromatogram point are lost. Files with 2+ points survive.**

## Root Cause: Division-by-Zero Chain

### Step 1: `avg_cycle_time = 0.0`

**File: `IntegrateChromatogramsSearch/utils.jl:104`**

```julia
avg_cycle_time = (chrom.rt[end] - chrom.rt[1]) / length(chrom.rt)
```

For a single-point chromatogram: `(rt - rt) / 1 = 0.0 / 1 = 0.0`

### Step 2: `x / 0.0 → Inf/NaN` in WHSmooth!

**File: `integrate_chrom.jl:144-145`**

```julia
rt_width = last_rt - start_rt  # = 0.0 for single point
x = x / rt_width               # x / 0.0 → all Inf/NaN!
```

The Whittaker-Henderson smoother then returns these Inf/NaN x-values (it correctly short-circuits the actual smoothing for m≤1 at line 147, but the damage to `x` is already done).

### Step 3: NaN Cascade Through Baseline Subtraction

**File: `integrate_chrom.jl:333-337`**

```julia
x_left = x[li]              # = Inf (from step 2)
x_right = x[ri]             # = Inf
dx = x_right - x_left       # = Inf - Inf = NaN
slope = (rmin - lmin) / dx   # = 0.0 / NaN = NaN
baseline = lmin + (xi - x_left) * slope  # = NaN
u[i] = max(0, u[i] - NaN)   # = NaN  (IEEE 754: max(0, NaN) = NaN)
```

### Step 4: NaN Peak Area

**File: `integrate_chrom.jl:276-278, 436`**

```julia
rt_norm = rt[stop] - start_rt  # = 0.0 (single point)
norm_factor = u[apex_scan+n_pad]  # = NaN (from step 3)
trapezoid_area = rt_norm * norm_factor * integrateTrapezoidal(state, avg_cycle_time)
# = 0.0 * NaN * (0.0 * NaN) = NaN
```

### Step 5: NaN Filtered Out

**File: `utils.jl:719`**

```julia
filter!(row -> !isnan(row.peak_area::Float32), psms)  # Removes NaN rows
```

The precursor is silently dropped with no warning or diagnostic.

## Proposed Fix

### Approach: Bypass `integrate_chrom` for Single-Point Chromatograms

For a chromatogram with only 1 data point, there is no chromatographic peak to smooth, fit, or integrate. The Whittaker-Henderson smoother requires ≥2 points for the difference matrix. Baseline subtraction requires ≥2 points for a meaningful baseline. Trapezoidal integration requires ≥2 points for a meaningful area. For a single point, the deconvolved intensity IS the quantitative measurement.

### Implementation

**File: `IntegrateChromatogramsSearch/utils.jl`**, inside `integrate_precursors`, after line 115 (the chromatogram view trimming) and before the apex scan finding:

```julia
chrom = view(chrom, first_pos:last_pos, :)

# ── Single-point bypass ──────────────────────────────────────────────
# For chromatograms with only 1 data point, the Whittaker-Henderson
# smoother, baseline subtraction, and trapezoidal integration are all
# undefined (division by zero in rt_width and avg_cycle_time).
# Use the deconvolved intensity directly as the peak area.
if size(chrom, 1) == 1
    peak_area[i] = chrom[1, :intensity]
    new_best_scan[i] = chrom[1, :scan_idx]
    points_integrated[i] = UInt32(1)
    if !ismissing(precursor_fraction_transmitted_traces)
        precursor_fraction_transmitted_traces[i], isotopes_captured_traces[i] =
            get_isolated_isotopes_strings(
                chrom[!, :precursor_fraction_transmitted],
                chrom[!, :isotopes_captured])
    end
    continue
end

# (existing code continues: apex scan finding, integrate_chrom call, etc.)
```

### What This Does

- **peak_area** = the deconvolved weight from that single scan. Always positive since we already filtered to `first_pos:last_pos` (positive intensities only).
- **new_best_scan** = the scan_idx of that single scan.
- **points_integrated** = 1.
- Skips `integrate_chrom` entirely via `continue`.

### What This Does NOT Change

- Multi-point chromatograms (≥2 points) follow the existing code path unchanged.
- No changes to `integrate_chrom.jl` itself.
- No changes to `process_final_psms!` filtering logic.
- No changes to scoring, protein inference, or any upstream steps.

### Unit Consistency Consideration

For multi-point chromatograms, `peak_area` has units of intensity × time (area under the chromatographic curve). For single-point chromatograms, the deconvolved intensity has units of intensity only. This means single-point peak_areas are not directly comparable in magnitude to multi-point peak_areas.

This is acceptable because:
1. **MaxLFQ normalizes within precursors** — it computes ratios of the same precursor across files, so the units cancel.
2. **A single-scan "peak" has no meaningful chromatographic width** — any time factor would be arbitrary (e.g., one cycle time ≈ 0.007 min would change magnitude but not relative quantification).
3. **The alternative is losing the identification entirely** — a precursor that passed 1% FDR in both global and per-file scoring should not be silently discarded because the integration math is undefined.
4. The `points_integrated=1` column flags these for downstream analysis.

## Expected Impact

- **~14,000 additional targets per file** retained in final output (from ~48k to ~62k).
- These are currently being lost silently — no warnings, no diagnostics.
- The fix is a simple `continue` bypass, not a change to the integration algorithm.
- Zero impact on existing multi-point chromatogram results.

## Cleanup

- Set `TRACE_PRECURSOR_IDX = UInt32(0)` in `ScoringSearch.jl` to disable diagnostic tracing.
- Rename `pi` → `pts_int` in `IntegrateChromatogramsSearch.jl` POST-INTEGRATE trace (potential conflict with Julia's `pi` constant).
- Update `scripts/INVESTIGATION_SINGLE_SCAN_LOSS.md` with root cause and fix reference.

## Verification Plan

1. Run `SearchDIA` on OlsenEclipse dataset with fix applied
2. Check that precursor_idx=6374682 now appears in all 6 files
3. Compare total precursor counts before/after: expect ~14k increase per file
4. Verify MaxLFQ protein quantification still produces reasonable results
5. Spot-check a few single-point precursors to confirm sensible peak_area values
