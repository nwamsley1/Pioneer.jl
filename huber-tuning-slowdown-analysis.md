# HuberTuningSearch Slowdown Analysis

## Symptom

HuberTuningSearch is taking several minutes per file when it should complete quickly:

```
[ Info: Executing Huber Tuning...
0.0%┣                                                                                            ┫ 0/8 [00:00<00:00, -0s/it]
```

Progress bar stuck at 0% for extended period (several minutes), indicating performance regression.

## Context

Recent changes on `feature/monotonic-rt-splines` branch:
1. ParameterTuningSearch: Increased default `spline_n_knots` from 5 → 100 (commit 32276a7b)
2. FirstPassSearch: Simplified to use adaptive knots (5-20) without monotonic enforcement (commit b5206150)
3. Added RT tolerance logging to FirstPassSearch and SecondPassSearch (commit 229fdc7d)
4. ParameterTuningSearch: Uses monotonic RT spline enforcement with 100 knots

## Possible Causes

### 1. RT Index Issues - Wider Tolerances (Most Likely)

**Hypothesis**: RT tolerance changes are causing HuberTuningSearch to select too many precursors per scan.

**Mechanism**:
- FirstPassSearch now calculates RT tolerances differently
- If tolerances are wider → RT index returns more precursors per scan
- HuberTuningSearch line 225: `irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))`
- Lines 226-227: RT index lookup with `irt ± irt_tol`
- More precursors → more transitions → more peak matching → slower

**Evidence to look for**:
- Check logged RT tolerance values from FirstPassSearch
- Compare with previous runs
- Add logging: Number of precursors selected per scan in HuberTuningSearch

**Test**:
```julia
# Add after line 227 in HuberTuningSearch/utils.jl
@user_info "Scan $scan_idx: RT bins [$irt_start_new:$irt_stop_new], width = $(irt_stop_new - irt_start_new)"
```

### 2. RT Model Evaluation Performance

**Hypothesis**: RT models with 100 knots are slow to evaluate or numerically unstable.

**Mechanism**:
- ParameterTuningSearch creates RT models with up to 100 knots
- Monotonic enforcement adds complexity
- Each scan calls `getRtIrtModel(search_context, ms_file_idx)(rt_value)`
- If model evaluation is slow, it compounds across all scans

**Evidence to look for**:
- Profile the `getRtIrtModel()` call
- Check if CPU usage is high (computing) or low (hanging)
- Time individual spline evaluations

**Test**:
```julia
# Before line 225
@time irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
```

### 3. Monotonic Enforcement Side Effects

**Hypothesis**: `make_spline_monotonic()` creates models with poor numerical properties.

**Mechanism**:
- Bidirectional cumulative max filter at ~5000 sample points
- Refit with `UniformSplinePenalized(..., λ=1.0, ...)` using 100 knots
- High penalty + many knots → very smooth (almost linear) splines
- Could cause numerical issues or poor RT predictions
- Might result in all precursors mapping to same RT bin (degeneracy)

**Evidence to look for**:
- Plot the RT models to see if they look reasonable
- Check for near-zero derivatives (too flat)
- Check RT index bin distribution

**Test**:
```julia
# After RT model creation in ParameterTuningSearch
rt_range = LinRange(minimum(valid_rt), maximum(valid_rt), 100)
irt_range = [rt_model(r) for r in rt_range]
@user_info "RT model span: $(minimum(irt_range)) to $(maximum(irt_range))"
```

### 4. Logging Overhead

**Hypothesis**: New `@user_info` logging statements are being called in tight loops.

**Mechanism**:
- Added logging to SecondPassSearch (might be similar pattern in HuberTuningSearch)
- If logging is per-scan or per-precursor, thousands of log calls
- I/O overhead accumulates

**Evidence to look for**:
- Check if log file is growing rapidly
- Count number of log statements in HuberTuningSearch
- Profile I/O operations

**Test**:
```julia
# Comment out all @user_info in HuberTuningSearch hot paths
# Re-run and compare performance
```

### 5. FirstPassSearch RT Models vs ParameterTuningSearch RT Models

**Hypothesis**: Confusion about which RT models should be used.

**Mechanism**:
- FirstPassSearch creates simple 5-20 knot splines (no monotonic enforcement)
- ParameterTuningSearch creates complex 100-knot splines (with monotonic enforcement)
- HuberTuningSearch should use ParameterTuningSearch models
- If using wrong models or models are missing, could cause issues

**Evidence to look for**:
- Check which RT models are actually being retrieved
- Verify ParameterTuningSearch completed successfully
- Check if RT models were overwritten

**Test**:
```julia
# In HuberTuningSearch, before first use
rt_model = getRtIrtModel(search_context, ms_file_idx)
@user_info "RT model type: $(typeof(rt_model))"
```

## Diagnostic Steps

### Step 1: Add Basic Logging
```julia
# In HuberTuningSearch/utils.jl, around line 218
@user_info "Processing scan $scan_idx of $(length(scan_range))"
```

### Step 2: Check Progress
Run with diagnostic logging and observe:
- Does it get past scan 1?
- Is it stuck on all files or just first file?
- CPU usage high (computing) or low (hanging)?

### Step 3: Check RT Tolerances
Review FirstPassSearch output for logged RT tolerances:
```
FirstPassSearch RT tolerances (per file):
  File 1 (...): RT tol = ??? min
```

Compare with previous runs. If significantly larger → Cause #1 likely.

### Step 4: Profile RT Model Evaluation
```julia
# Time a few RT model evaluations
using BenchmarkTools
rt_model = getRtIrtModel(search_context, 1)
@btime $rt_model(50.0)  # Should be nanoseconds
```

If slow (>microseconds) → Cause #2 likely.

### Step 5: Check RT Index Bin Selection
```julia
# After line 227
n_bins_selected = irt_stop_new - irt_start_new
@user_info "Scan $scan_idx: Selected $n_bins_selected RT bins"
```

If consistently large (>50 bins) → Cause #1 likely.

## Quick Experiments

### Experiment 1: Revert spline_n_knots
```julia
# In ParameterTuningSearch/types.jl line 344
spline_n_knots::Int64 = 5,  # Temporarily revert from 100
```

Re-run. If fast → Cause #2 or #3 confirmed.

### Experiment 2: Disable Monotonic Enforcement
```julia
# In rt_alignment_utils.jl line 380-386
# Comment out:
# final_map_monotonic = make_spline_monotonic(...)
# Use rt_to_irt_map directly instead
```

Re-run. If fast → Cause #3 confirmed.

### Experiment 3: Check RT Tolerance Impact
```julia
# In HuberTuningSearch/utils.jl line 217
irt_tol_original = getIrtErrors(search_context)[ms_file_idx]
irt_tol = irt_tol_original * 0.5  # Artificially reduce by half
@user_info "Reduced RT tolerance from $irt_tol_original to $irt_tol"
```

Re-run. If fast → Cause #1 confirmed.

## Questions for User

1. **Which file is it stuck on?** File 0/8, or does it eventually progress?
2. **CPU usage during hang?** High (100%) or low (<20%)?
3. **RT tolerance values logged?** What were the values from FirstPassSearch?
4. **Did ParameterTuningSearch complete?** Were there any warnings/errors?
5. **Log file size?** Growing rapidly or static during hang?

## Expected Performance

HuberTuningSearch should be **fast** because:
- Only processes high-confidence PSMs (q-value < 0.01)
- Typically 100-1000 scans per file
- Simple peak matching and scoring
- No iterative optimization in the hot path
- Should take seconds per file, not minutes

If taking minutes per file → significant performance regression requiring investigation.

## Related Files

- `src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl` (lines 215-260)
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl` (line 344)
- `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl` (lines 270-390)
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl` (lines 553-571)

## Next Steps

1. Add diagnostic logging to narrow down the issue
2. Run one of the quick experiments to test hypotheses
3. Profile the hot path if computational slowdown confirmed
4. Check data integrity if hanging/deadlock suspected
