# Plan: Improving Poisson-Derived Scoring Metrics

## Status

`fitted_hellinger` is a major success — ranked #1 in LightGBM feature importance on OlsenAstral (108K gain). The other four Poisson metrics (`poisson_deviance`, `poisson_dispersion`, `max_dev_resid_matched`, `max_dev_resid_unmatched`) show weak discriminative power despite strong theoretical motivation.

This document analyzes why and proposes fixes.

## Current Results (OlsenAstral 3-file)

### FirstPass LightGBM (20 features)
```
fitted_hellinger                         108813.14   ← #1, huge success
max_matched_residual                     89086.11
irt_error                                30471.12
poisson                                  26517.67
fitted_manhattan_distance                26206.9
err_norm                                 21804.92
gof                                      5737.25
total_ions                               4197.24
y_count                                  4029.93
max_dev_resid_unmatched                  2330.07    ← modest
...
poisson_deviance                         85.01      ← weak
max_dev_resid_matched                    72.07      ← weak
poisson_dispersion                       6.09       ← nearly useless
```

## Diagnosis: Why Are 3 of 4 Metrics Weak?

### Issue 1: Per-Precursor Shadow vs Fitted Comparison

The current implementation computes Poisson metrics using:
- **μ_i** = `fitted_peak` = `w[col] * H.nzval[i]` — this precursor's model prediction
- **x_i** = `shadow_peak` = `fitted_peak - r[H.rowval[i]]` — the deconvolved observed signal attributed to this precursor

The problem: for well-deconvolved spectra, `shadow_peak ≈ fitted_peak` by construction. The deconvolution solver *minimizes* the residuals, so `r ≈ 0` when the fit is good, making `shadow ≈ fitted`. This means:

- **poisson_deviance**: `Σ [x*log(x/μ) - (x-μ)]` ≈ 0 for all PSMs because x ≈ μ
- **poisson_dispersion**: `Σ (x-μ)²/μ` / (n-1) ≈ 0 for all PSMs because x ≈ μ
- **max_dev_resid_matched**: `max sqrt(2[x*log(x/μ) - (x-μ)])` ≈ 0 for same reason

This is fundamentally different from what the LaTeX paper describes. The paper's metrics compare **raw observed** intensities against the **total model prediction** across all precursors — measuring whether Poisson noise explains the residuals. Our per-precursor computation compares one precursor's contribution against its own deconvolved share.

### Why fitted_hellinger works despite the same inputs

Hellinger distance compares **normalized shapes** via the Bhattacharyya coefficient:
```
H² = 1 - Σ sqrt(μ_i * x_i) / sqrt(Σμ · Σx)
```

Even when `shadow ≈ fitted` in absolute terms, small shape differences (relative intensity profile mismatches) are captured because:
1. The sqrt-transform amplifies differences at low intensities
2. The normalization makes it scale-free — it's comparing *proportions*, not magnitudes
3. True PSMs have nearly identical fitted/shadow *shapes*; false PSMs have shape mismatches even when the deconvolution minimizes total residual

### Issue 2: Possible Float16 Quantization

All spectral scores are stored as `Float16` (via `SpectralScoresFirstPass{Float16}`). Float16 has:
- **10-bit mantissa** → ~3.3 decimal digits precision
- **Range**: ~6e-8 to 65,504
- **Precision at 1.0**: ±0.001

If deviance/dispersion values are clustered in a narrow range near zero, Float16 may quantize many distinct Float32 values to the same Float16 representation, destroying discriminative information before LightGBM ever sees it.

**Tested empirically**: Float16 is NOT the bottleneck. At typical deviance/dispersion values (0.001-1.0), Float16 preserves 3 significant digits — plenty for LightGBM to split on. Even at 0.0005, Float16 distinguishes it from 0.001. The issue is that the *values themselves* are uninformative (all near zero), not that Float16 can't represent them.

**Conclusion**: No need for Float32 storage or log-transform for this reason. The fix is to change what we're measuring (Issue 3), not how we store it.

### Issue 3: The True Chimera Detector Needs Scan-Level Computation

The paper's Poisson dispersion is defined as:
```
φ̂_j = (1/(n-1)) Σ_i (x_i - μ_total_i)² / μ_total_i
```

where:
- **x_i** = raw observed intensity at peak i (`H.x[i]`)
- **μ_total_i** = total model prediction = `Σ_j w_j * H_{ij}` for all precursors j

This measures whether the **overall fit residuals** at peaks belonging to precursor j are consistent with Poisson noise. A chimeric spectrum has extra signal from an un-modeled co-eluting precursor, inflating residuals above Poisson expectation.

Currently we DON'T have `μ_total_i` in the per-precursor loop. But we DO have:
- `H.x[i]` — raw observed
- `r[H.rowval[i]]` = `μ_total - x` at that peak's row

Therefore: `μ_total = H.x[i] + r[H.rowval[i]]`

This means we CAN compute the paper's version within the existing loop:

```julia
mu_total = H.x[i] + r[H.rowval[i]]  # total model prediction at this peak
x_raw = H.x[i]                       # raw observed
if mu_total > zero(T)
    chi2_raw_sum += (x_raw - mu_total)^2 / mu_total
end
```

This would give the *actual* scan-level Poisson dispersion attributed to peaks of this precursor.

## Proposed Changes

### Fix 1: Compute scan-level Poisson dispersion (the real chimera detector)

Replace the current per-precursor dispersion with the paper's version using raw observed vs total model prediction:

```julia
# In the per-peak loop:
mu_total_i = H.x[i] + r[H.rowval[i]]   # = Σ_j w_j H_ij (total model prediction)
x_raw_i = H.x[i]                         # raw observed intensity
if mu_total_i > zero(T)
    chi2_raw_sum += (x_raw_i - mu_total_i)^2 / mu_total_i
end

# After loop:
poisson_dispersion = n_frags > 1 ? chi2_raw_sum / (n_frags - 1) : zero(T)
```

**Interpretation**:
- ≈ 1.0 → residuals consistent with Poisson noise → clean PSM
- >> 1.0 → excess variance → chimeric co-isolation or wrong ID
- << 1.0 → overfitting or degenerate solution

**Cost**: Same as current (1 div/peak), just different inputs.

### Fix 2: Similarly fix poisson_deviance to use raw vs total model

```julia
mu_total_i = H.x[i] + r[H.rowval[i]]
x_raw_i = H.x[i]
if x_raw_i > zero(T) && mu_total_i > zero(T)
    deviance_sum += x_raw_i * log(x_raw_i / mu_total_i) - (x_raw_i - mu_total_i)
elseif mu_total_i > zero(T)
    deviance_sum += mu_total_i
end
```

This measures how well the **total model** explains the observed data at this precursor's peaks. True PSMs will have small deviance (good fit at their peaks), while false IDs or chimeras will have larger deviance.

### Fix 3: Similarly fix deviance residuals

```julia
mu_total_i = H.x[i] + r[H.rowval[i]]
x_raw_i = H.x[i]
if x_raw_i > zero(T) && mu_total_i > zero(T)
    dev_contrib = T(2) * (x_raw_i * log(x_raw_i / mu_total_i) - (x_raw_i - mu_total_i))
else
    dev_contrib = T(2) * mu_total_i
end
dev_resid = sqrt(max(dev_contrib, zero(T)))
```

### Fix 4: Address Float16 quantization

Option A: **Log-transform before storing** — similar to how `gof` uses `-log2(ratio)`:
```julia
poisson_deviance = n_frags > 0 ? -log2(max(deviance_sum / n_frags, T(1e-10))) : zero(T)
poisson_dispersion = n_frags > 1 ? log2(max(chi2_raw_sum / (n_frags - 1), T(1e-10))) : zero(T)
```

This spreads values across Float16's representable range and matches the log2 convention used by other metrics (gof, max_matched_residual, etc).

Option B: **Store as Float32 in `FirstPassScoredPSM`** — use `H` (high precision) type parameter instead of `L` (low precision) for Poisson fields. This is a larger refactor and may not be worth it if log-transform works.

**Recommendation**: Try log-transform first (Option A). It's one-line changes and matches existing conventions. If that doesn't help, Float16 wasn't the bottleneck.

### Fix 5: Keep fitted_hellinger as-is

The current computation using shadow_peak vs fitted_peak with Bhattacharyya coefficient is correct and performs excellently. No changes needed.

## Summary of Changes

| Metric | Current (shadow vs fitted) | Proposed (raw vs total model) | Transform |
|--------|---------------------------|-------------------------------|-----------|
| `fitted_hellinger` | Keep as-is | No change | -log2(H²) |
| `poisson_deviance` | x=shadow, μ=fitted | x=H.x[i], μ=H.x[i]+r[i] | -log2(D/n) |
| `poisson_dispersion` | x=shadow, μ=fitted | x=H.x[i], μ=H.x[i]+r[i] | log2(χ²/(n-1)) |
| `max_dev_resid_matched` | x=shadow, μ=fitted | x=H.x[i], μ=H.x[i]+r[i] | raw (already ~N(0,1)) |
| `max_dev_resid_unmatched` | x=shadow, μ=fitted | x=H.x[i], μ=H.x[i]+r[i] | raw |

## Testing Plan

1. **Float16 diagnostic**: Before changing the computation, add temporary logging to check if Float16 quantization is the bottleneck. Print `length(unique(Float16.(values)))` vs `length(unique(values))` for each metric.

2. **A/B test computation change**: Switch to raw-vs-total computation, run on OlsenAstral, compare feature importances.

3. **A/B test log-transform**: Apply log-transform to deviance and dispersion, compare feature importances.

4. **Regression test**: Verify total unique targets at 1% FDR doesn't decrease with changes.

## Files to Modify

| File | Change |
|------|--------|
| `src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl` | Update getDistanceMetrics FirstPass variant |
| (Optional) `src/Routines/SearchDIA/PSMs/ScoredPSMs.jl` | If changing Float16→Float32 for Poisson fields |

## Safe Checkpoint

Revert to commit `696a1b9b` if changes break things. That commit has the working Poisson metrics (weak but functional) and the strong `fitted_hellinger`.
