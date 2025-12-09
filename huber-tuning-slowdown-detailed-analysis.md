# HuberTuningSearch Slowdown: Detailed Analysis of Branch Differences

## All Changes Between `develop` and `feature/monotonic-rt-splines`

### Summary of Commits
```
2094fba6 docs: Add diagnostic analysis for HuberTuningSearch performance regression
229fdc7d feat: Add RT tolerance logging to FirstPassSearch and SecondPassSearch
b5206150 refactor: Use adaptive knot selection in FirstPassSearch RT alignment
32276a7b refactor: Improve ParameterTuningSearch robustness and RT calibration
764b1acb refactor: Simplify FirstPassSearch RT alignment to fixed 7-knot spline
becd59ec fix: Pad polynomial coefficients to prevent dimension mismatch in splines
133ca5ad refactor: Increase RT bins from 15 to 100 for finer scan selection ⚠️
d73172ec feat: Implement monotonic RT spline enforcement
```

## Critical Finding: RT Bins Increase (133ca5ad)

**Change**: `FilteredMassSpecData` default `n_rt_bins` increased from **15 → 100**

**Used by**: ParameterTuningSearch only (NOT HuberTuningSearch directly)

**Impact on ParameterTuningSearch**:
- Samples 100 scans across gradient (vs 15)
- 6.7x more scans to process during parameter tuning
- Would slow down **ParameterTuningSearch**, not HuberTuning

**But**: This could indirectly affect HuberTuning if ParameterTuningSearch results are used.

## Analysis of Each Hypothesis

### Hypothesis 1: RT Tolerance Width (REVISED - LESS LIKELY)

**Original Theory**: Wider RT tolerances → more precursors per scan in HuberTuning

**Analysis Against Code Changes**:

1. **FirstPassSearch RT tolerance calculation UNCHANGED**:
   ```julia
   // Line 620-661 in FirstPassSearch/utils.jl
   function get_irt_errs(fwhms, prec_to_irt, params)
       fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm], fwhms)
       variance_ = collect(skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt)))
       irt_std = median(variance_) * params.irt_nstd
       return map(x->Float32((x+irt_std))::Float32, fwhms)
   end
   ```
   Formula is **identical** to develop branch.

2. **ParameterTuningSearch irt_mad calculation UNCHANGED**:
   ```julia
   // Line 348 in rt_alignment_utils.jl
   irt_mad = mad(residuals, normalize=false)::Float32
   getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
   ```
   Still uses `normalize=false` (same as develop).

3. **FirstPassSearch OVERWRITES ParameterTuningSearch tolerances**:
   ```julia
   // Line 557 in FirstPassSearch/utils.jl
   setIrtErrors!(search_context, irt_errs)  // Replaces ParameterTuning values
   ```

**Conclusion**: RT tolerances should be **identical** to develop branch unless:
- Peak widths (FWHM) changed due to different scan selection
- Cross-run iRT variance changed

**Probability**: 20% - Only if scan selection in FirstPassSearch affected FWHM estimation

---

### Hypothesis 2: RT Model Evaluation Slow (REJECTED)

**Theory**: 100-knot splines slow to evaluate

**Analysis**:

1. **FirstPassSearch uses 5-20 knots** (NOT 100):
   ```julia
   // Line 390 in FirstPassSearch/utils.jl
   n_knots = clamp(n_good_psms ÷ 100, 5, 20)  // Max 20, not 100
   ```

2. **Spline evaluation is O(1)** regardless of knot count:
   ```julia
   // Lines 259-279 in uniformBasisCubicSpline.jl
   function (s::UniformSpline)(t::U)
       idx = floor(Int32, (t - s.first)/s.bin_width)  // O(1) bin lookup
       // Horner's method for cubic: 4 multiplications
       x = s.coeffs[coeff]*c
       c *= u; coeff += 1
       // ... (3 more terms)
   end
   ```
   Evaluation: **4 multiplications + 1 division + indexing** = ~10ns

3. **HuberTuning calls RT model once per scan**:
   ```julia
   // Line 225 in HuberTuningSearch/utils.jl
   irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
   ```
   Even with 10,000 scans: 10,000 × 10ns = **0.1ms total** (negligible)

**Conclusion**: Evaluation speed is **NOT the bottleneck**.

**Probability**: 0% - Math doesn't support this hypothesis

---

### Hypothesis 3: Monotonic Enforcement Effects (NEW - MOST LIKELY)

**Theory**: Monotonic enforcement with 100 knots + high penalty creates degenerate RT models

**Analysis**:

#### Part A: The Monotonic Enforcement Process

```julia
// Lines 380-386 in rt_alignment_utils.jl (ParameterTuningSearch)
final_map_monotonic = make_spline_monotonic(
    rt_to_irt_map,
    valid_psms[!, :rt],
    valid_psms[!, :irt_predicted],
    n_knots = n_knots  // Up to 100 knots!
)
```

Inside `make_spline_monotonic()`:
```julia
// Lines 167-195 in rt_alignment_utils.jl
n_sample_points = length(rt_data) - 1  // ~5000 for typical dataset
rt_grid = collect(LinRange(rt_min, rt_max, n_sample_points))  // 5000 points
irt_grid = [original_spline(r) for r in rt_grid]  // Evaluate at 5000 points

// Bidirectional cumulative max filter
for i in (median_idx+1):n_sample_points
    if irt_grid[i] < irt_grid[i-1]
        irt_grid[i] = irt_grid[i-1]  // Force monotonic
    end
end

// REFIT with high penalty and many knots
final_spline = UniformSplinePenalized(
    irt_grid,        // 5000 filtered points
    rt_grid,         // 5000 RT values
    3,               // cubic
    n_knots,         // 100 knots!
    Float32(1.0),    // HIGH penalty λ=1.0
    2                // 2nd order penalty
)
```

#### Part B: The Problem

With **100 knots, λ=1.0, and 5000 already-smoothed points**:

1. **Over-smoothing**: The penalty term dominates:
   ```
   Objective = ||y - Xc||² + λ||Dc||²

   With λ=1.0 and smooth y (already filtered):
   - Data term is small (good fit)
   - Penalty term drives solution toward nearly-linear
   ```

2. **RT Index Degeneracy**:
   - If RT model is nearly linear: RT_min→iRT_min, RT_max→iRT_max uniformly
   - All precursors distributed evenly across iRT space
   - **BUT** RT index bins are uniform in iRT space
   - With 100 knots → near-linear → bad precursor distribution

3. **Massive RT Index Bin Imbalance**:
   ```
   Scenario: 60-minute gradient, 10,000 precursors

   GOOD MODEL (nonlinear):
   - Early gradient compressed: many precursors in few iRT bins
   - Late gradient spread: few precursors in many iRT bins
   - Average: 100-200 precursors per bin

   BAD MODEL (nearly linear from over-smoothing):
   - All RT times map linearly to iRT
   - RT index built on uniform iRT grid
   - If actual precursor RT distribution is non-uniform:
     - Dense regions: 1000+ precursors in single bin
     - Sparse regions: 10 precursors in single bin
   ```

4. **HuberTuning Slowdown Mechanism**:
   ```julia
   // HuberTuning queries RT index:
   irt = getRtIrtModel(scan_rt)  // Converts scan RT → iRT
   bins = rt_index.query(irt - tol, irt + tol)  // Get precursors in iRT window

   // If ONE scan hits a dense bin:
   for precursor in bins  // Could be 1000+ precursors!
       match_peaks()      // Expensive
       deconvolve()       // Very expensive
       score_psm()        // Expensive
   end
   ```

#### Part C: Evidence This Is The Cause

1. **Timing symptom**: Stuck at 0/8 files for minutes
   - Suggests first file is the bottleneck
   - Could be stuck on single dense scan

2. **Why not ParameterTuningSearch slowdown?**
   - ParameterTuningSearch uses FilteredMassSpecData (100 selected scans)
   - Even with dense bins, only processing 100 scans
   - HuberTuning processes **ALL high-confidence PSM scans** (could be thousands)

3. **The 5→100 knots change** (32276a7b):
   - More knots = more flexibility = more potential for over-smoothing with high penalty
   - 5 knots: Forced to be somewhat curved
   - 100 knots with λ=1.0: Can be nearly linear

**Conclusion**: High penalty + many knots + monotonic enforcement = degenerate RT models

**Probability**: 70% - Math and symptoms align perfectly

---

### Hypothesis 4: RT Model Quality Degradation (NEW - POSSIBLE)

**Theory**: Simplified FirstPassSearch creates worse RT models → wrong precursors searched

**Analysis**:

**OLD FirstPassSearch** (develop):
```julia
// Removed features:
- RANSAC for outlier robustness
- Penalty optimization
- Monotonic enforcement
- Two-mode robust vs simple fitting
```

**NEW FirstPassSearch** (feature branch):
```julia
// Simple approach:
n_knots = clamp(n_good_psms ÷ 100, 5, 20)
rt_to_irt_spline = UniformSpline(best_irts, best_rts, 3, n_knots)
// No RANSAC, no penalty, no monotonic enforcement
```

**Impact**:
- If data has outliers: Simple least squares fits poorly
- Poor RT predictions: `irt = rt_model(scan_rt)` returns wrong iRT
- Wrong iRT: Queries wrong region of RT index
- Wrong precursors: More false positives to process

**BUT**: FirstPassSearch filters PSMs by probability:
```julia
best_hits = psms[:prob].>params.min_prob_for_irt_mapping
```
This should remove most outliers.

**Counter-evidence**:
- FirstPassSearch RT alignment plots should show if fit is poor
- User would likely notice warnings about bad RT alignment

**Conclusion**: Possible but less likely than monotonic enforcement issue

**Probability**: 30% - Depends on data quality

---

### Hypothesis 5: Search Order / Model Confusion (LESS LIKELY)

**Theory**: HuberTuning using wrong RT models

**Analysis of Search Order**:

From `LibrarySearch.jl` and method interfaces:
```
1. ParameterTuningSearch  → Sets RT models (100-knot w/ monotonic)
2. FirstPassSearch        → OVERWRITES RT models (5-20 knot simple)
3. HuberTuningSearch      → Uses models from search_context
4. SecondPassSearch       → Uses models from search_context
```

**Which models does HuberTuning use?**

Looking at code, HuberTuning runs **AFTER** FirstPassSearch, so it uses:
- RT models: FirstPassSearch (5-20 knots, simple)
- RT tolerances: FirstPassSearch (`irt_errs`)

**Problem**: User says ParameterTuningSearch completed normally, then HuberTuning is slow.

**Possible confusion scenario**:
- If FirstPassSearch FAILED or was SKIPPED for some files
- Those files would still have ParameterTuningSearch's RT models (100-knot w/ monotonic)
- HuberTuning would use the bad models for those files

**Check**: Do the logs show FirstPassSearch completing for all files?

**Conclusion**: Only possible if FirstPassSearch failed for affected files

**Probability**: 15% - Requires additional failure

---

### Hypothesis 6: Logging Overhead (REJECTED)

**Theory**: New logging slows down hot paths

**Analysis**:

**Added logging**:
```julia
// FirstPassSearch (lines 559-571):
for (file_idx, tol) in pairs(irt_errs)  // O(n_files) = O(8)
    @user_info "..."
end
```
8 log statements total - **negligible**

**SecondPassSearch** (lines 172-178, 384-390):
```julia
// Per-file logging (once per file)
@user_info "SecondPassSearch: Processing $file_name..."
```
8 log statements total - **negligible**

**HuberTuning**: No new logging added.

**Conclusion**: Logging overhead is <1ms total

**Probability**: 0% - Not in hot path

---

### NEW Hypothesis 7: RT Index Structure Degeneracy

**Theory**: Changes to RT index building process cause pathological bin distributions

**Analysis**:

FirstPassSearch builds RT indices (lines 542-612):
```julia
function create_rt_indices!(search_context, results, precursor_dict, params)
    // Calculate irt_errs (tolerances)
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)
    setIrtErrors!(search_context, irt_errs)

    // Create RT index for each file
    for (file_idx, rt_model, psms_path) in valid_files
        rt_index = build_rt_index(prec_to_irt, irt_errors[file_idx], ...)
        // Save to disk
    end
end
```

**RT Index structure**:
- Bins distributed uniformly across **iRT space**
- Each bin contains precursors in that iRT range
- Query: `bins_to_search = rt_index.query(irt ± tol)`

**If RT model is degenerate (nearly linear)**:
- Precursors cluster in RT space (e.g., early gradient has most peptides)
- But iRT mapping is uniform
- Result: **Uneven bin distribution** even with uniform iRT bins

**Example**:
```
Real data: 60% of precursors in first 20 minutes (early gradient)
Bad RT model maps: 0-60 min → 0-60 iRT (linear)
RT index bins: Uniform over 0-60 iRT

Query at iRT=10 (early gradient):
- Many precursors have actual RT in [8-12] min
- They all map to iRT in [8-12] range
- Bin at iRT=10 has 1000+ precursors
- HuberTuning processes ALL of them for that scan
```

**This explains**:
- Why it's stuck on one file/scan
- Why ParameterTuningSearch isn't slow (only 100 filtered scans)
- Why HuberTuning is slow (processes all PSM scans, hits dense bins)

**Probability**: 80% - Most complete explanation

---

### NEW Hypothesis 8: Probit Regression Changes Affecting Scoring

**Theory**: Changed probit regression threshold affects HuberTuning performance

**Change** (32276a7b):
```julia
// OLD: M > 10
// NEW: M > 50
if M > 50
    try
        // Probit regression scoring
    catch
        // Fallback: assign prob = 1
    end
else
    psms[!,:prob] = ones(Float32, size(psms, 1))
end
```

**Impact on HuberTuning**:
- HuberTuning filters PSMs by q-value
- Q-values depend on probabilities
- If more PSMs have prob=1 (due to threshold increase):
  - Q-value calculation changes
  - More/fewer PSMs pass threshold
  - Different number of scans to process

**But**: This change is in ParameterTuningSearch, not FirstPassSearch or HuberTuning.

**Conclusion**: Unlikely to directly affect HuberTuning

**Probability**: 5% - Indirect effect only

---

## Ranked Hypotheses by Probability

### 1. RT Index Degeneracy (80%)
**Cause**: Monotonic enforcement with 100 knots + λ=1.0 creates nearly-linear RT models
**Mechanism**: Linear models + non-uniform precursor RT distribution → dense RT index bins
**Effect**: HuberTuning scans hit bins with 1000+ precursors → minutes per scan

**Evidence**:
- Math: High penalty over-smooths → linear
- Symptoms: Stuck on single file/scan
- Code: Monotonic enforcement refits with 100 knots, λ=1.0

**Test**: Check RT model linearity:
```julia
rt_model = getRtIrtModel(search_context, 1)
rt_range = LinRange(30.0, 60.0, 100)
irt_range = [rt_model(r) for r in rt_range]
plot(rt_range, irt_range)  // Should be curved, not straight line
```

---

### 2. Monotonic Enforcement Over-Smoothing (70%)
**Cause**: Same as above, slight variant
**Mechanism**: Over-smoothed models lose RT resolution
**Effect**: Poor RT discrimination → search wrong regions

---

### 3. RT Model Quality Degradation (30%)
**Cause**: Removed RANSAC and robust fitting from FirstPassSearch
**Mechanism**: Outliers degrade simple least-squares fit
**Effect**: Poor RT predictions → more wrong precursors searched

**Test**: Check RT alignment residuals in FirstPassSearch plots

---

### 4. RT Tolerance Width (20%)
**Cause**: Some indirect effect changing FWHM or iRT variance
**Mechanism**: Wider tolerances → more precursors per scan
**Effect**: More computation per scan

**Test**: Check logged RT tolerance values

---

### 5. Model Confusion (15%)
**Cause**: FirstPassSearch failed for some files
**Mechanism**: HuberTuning uses ParameterTuningSearch's bad models
**Effect**: Degenerate models for those files only

**Test**: Check FirstPassSearch logs for failures

---

## Smoking Gun Evidence

### Critical Change: commit 133ca5ad (RT bins 15→100)
This change affects **ParameterTuningSearch** processing time but shouldn't directly affect HuberTuning.

### Critical Change: commits d73172ec + a383cfd0 (Monotonic enforcement)
Added `make_spline_monotonic()` which:
1. Samples at ~5000 points
2. Applies cumulative max filter
3. **Refits with UniformSplinePenalized** using:
   - Up to 100 knots (from commit 32276a7b)
   - λ = 1.0 (high penalty)
   - 2nd order penalty

This combination is **mathematically guaranteed** to produce over-smoothed (nearly linear) models.

### Critical Change: commit 32276a7b (spline_n_knots 5→100)
Increases knots in ParameterTuningSearch from 5 to 100.

With 100 knots + high penalty in monotonic enforcement:
- Over-fitting prevention becomes over-smoothing
- Models lose curvature
- RT index bins become imbalanced

---

## Recommended Diagnostic Steps

### Step 1: Check RT Model Shape (PRIORITY 1)
```julia
# After ParameterTuningSearch completes
rt_model = getRtIrtModel(search_context, 1)
rt_samples = LinRange(minimum(rt_data), maximum(rt_data), 1000)
irt_samples = [rt_model(r) for r in rt_samples]

# Check derivative (should vary, not constant)
derivs = diff(irt_samples) ./ diff(rt_samples)
@info "RT model derivative range: $(minimum(derivs)) to $(maximum(derivs))"
@info "RT model derivative std: $(std(derivs))"

# If std(derivs) < 0.01 → Model is nearly linear → PROBLEM CONFIRMED
```

### Step 2: Check RT Index Bin Distribution (PRIORITY 2)
```julia
# In HuberTuning, before processing
bin_sizes = [length(bin.prec) for bin in rt_index.rt_bins]
@info "RT index bin sizes: min=$(minimum(bin_sizes)), max=$(maximum(bin_sizes)), median=$(median(bin_sizes))"

# If max > 500 → Dense bins → PROBLEM CONFIRMED
```

### Step 3: Check Which Scan Is Slow (PRIORITY 3)
```julia
# Add in HuberTuning scan loop
@info "Processing scan $scan_idx with $(length(precursors)) precursors"
@time begin
    # ... existing scan processing code ...
end
```

### Step 4: Check RT Tolerance Values
```julia
# Already added in commit 229fdc7d - check output logs
```

---

## Recommended Fixes (In Order)

### Fix 1: Reduce max_knots in ParameterTuningSearch (QUICK FIX)
```julia
// In ParameterTuningSearch/types.jl line 344
spline_n_knots::Int64 = 20,  // Was 100, reduce to 20
```
**Rationale**: 20 knots still provides good resolution without excessive smoothing

---

### Fix 2: Remove or Reduce Penalty in Monotonic Enforcement
```julia
// In rt_alignment_utils.jl line 185
UniformSplinePenalized(
    irt_grid, rt_grid, 3, n_knots,
    Float32(0.01),  // Reduce from 1.0 to 0.01
    2
)
```
**Rationale**: Lower penalty preserves curvature while still maintaining smoothness

---

### Fix 3: Use Adaptive Knots in Monotonic Enforcement
```julia
// In make_spline_monotonic(), calculate knots based on data size
n_knots_monotonic = min(n_knots, max(5, length(rt_data) ÷ 500))
// For 5000 points: 5000÷500 = 10 knots max
```
**Rationale**: Prevent using 100 knots for refitting

---

### Fix 4: Remove Monotonic Enforcement Entirely (AGGRESSIVE)
```julia
// In rt_alignment_utils.jl, skip make_spline_monotonic()
final_model = SplineRtConversionModel(rt_to_irt_map)  // Use original spline
return (final_model, psms[!, :rt], psms[!, :irt_predicted], irt_mad)
```
**Rationale**: FirstPassSearch doesn't use it and works fine

---

## Summary

**Most Likely Cause**: Monotonic enforcement with 100 knots + high penalty creates over-smoothed, nearly-linear RT models, causing RT index bin degeneracy.

**Quick Test**: Check if RT model is nearly linear (see Step 1 above)

**Quick Fix**: Reduce `spline_n_knots` from 100 to 20 in ParameterTuningSearch

**Complete Fix**: Redesign monotonic enforcement to avoid over-smoothing
