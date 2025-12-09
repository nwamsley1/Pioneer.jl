# FirstPassSearch RT Alignment Simplification Plan

## Goal

Simplify FirstPassSearch RT alignment to use a fixed 7-knot UniformSpline approach with no double-fitting, no RANSAC, and no monotonic enforcement.

## Current State Problems

1. **Complex two-mode system**: `use_robust_fitting` parameter creates two different code paths
2. **Double-fitting**: Initial fit → monotonic enforcement → refit with high penalty
3. **Adaptive knot selection**: Complex logic for determining number of knots
4. **Multiple spline types**: Uses both `UniformSpline` and `UniformSplinePenalized`
5. **High penalty in monotonic refit**: λ=1.0 causes boundary polynomials to have zero coefficients

## Proposed Simplification

**Single, simple approach**: Use `UniformSpline` with fixed 7 knots, no penalty, no RANSAC, no monotonic enforcement.

### Benefits

1. **Simplicity**: One code path, easy to understand and maintain
2. **Speed**: No RANSAC iterations, no penalty optimization, no double-fitting
3. **Reliability**: Fewer moving parts = fewer failure modes
4. **Consistency**: Same behavior across all datasets
5. **Proven**: This is essentially what the legacy `use_robust_fitting=false` mode does, which works

### Tradeoffs

1. **No outlier robustness**: RANSAC was handling outliers in small datasets
2. **No monotonicity guarantee**: Spline may have small non-monotonic regions (rare)
3. **Fixed knots**: Not adaptive to data size
   - 7 knots works well for typical datasets (500-10000 PSMs)
   - May be over-parameterized for very small datasets (<100 PSMs)
   - May be under-parameterized for very large datasets (>100000 PSMs)

---

## Detailed Changes

### File 1: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

#### Change 1: Simplify RT alignment function (Lines 350-450)

**Current Code:**
```julia
for ms_file_idx in valid_files
    # Skip files that have been marked as failed
    if is_file_failed(search_context, ms_file_idx)
        continue
    end
    psms_path = all_psms_paths[ms_file_idx]
    psms = Arrow.Table(psms_path)
    best_hits = psms[:prob].>params.min_prob_for_irt_mapping
    try
        if params.use_robust_fitting
            # Use robust fitting with RANSAC/penalty
            best_psms_df = DataFrame(
                rt = psms[:rt][best_hits],
                irt_predicted = psms[:irt_predicted][best_hits]
            )

            # Fit RT to iRT model using shared robust fitting
            rt_model, valid_rt, valid_irt, irt_mad = Pioneer.fit_irt_model(
                best_psms_df;
                lambda_penalty = Float32(0.1),
                ransac_threshold = 1000,
                min_psms = 10,
                spline_degree = 3,
                max_knots = 10,
                outlier_threshold = Float32(5.0)
            )

            # Store RT to iRT model
            setRtIrtMap!(search_context, rt_model, ms_file_idx)

            # Fit inverse model (iRT to RT) - swap input/output
            irt_to_rt_df = DataFrame(
                rt = valid_irt,           # Use iRT as input
                irt_predicted = valid_rt  # Use RT as output
            )
            irt_model, _, _, _ = Pioneer.fit_irt_model(
                irt_to_rt_df;
                lambda_penalty = Float32(0.1),
                ransac_threshold = 1000,
                min_psms = 10,
                spline_degree = 3,
                max_knots = 10,
                outlier_threshold = Float32(5.0)
            )
            setIrtRtMap!(search_context, irt_model, ms_file_idx)

            # Optionally generate plots
            if params.plot_rt_alignment
                plot_rt_alignment_firstpass(
                    valid_rt,
                    valid_irt,
                    rt_model,
                    ms_file_idx,
                    getDataOutDir(search_context)
                )
            end
        else
            # Use simple UniformSpline (legacy behavior)
            best_rts = psms[:rt][best_hits]
            best_irts = psms[:irt_predicted][best_hits]

            irt_to_rt_spline = UniformSpline(
                best_rts,
                best_irts,
                3,
                5
            )
            rt_to_irt_spline = UniformSpline(
                best_irts,
                best_rts,
                3,
                5
            )

            #Build rt=>irt and irt=> rt mappings for the file and add to the dictionaries
            setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
            setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)
        end
    catch e
        # ... error handling ...
    end
end
```

**New Code:**
```julia
for ms_file_idx in valid_files
    # Skip files that have been marked as failed
    if is_file_failed(search_context, ms_file_idx)
        continue
    end
    psms_path = all_psms_paths[ms_file_idx]
    psms = Arrow.Table(psms_path)
    best_hits = psms[:prob].>params.min_prob_for_irt_mapping

    try
        # Check if we have enough PSMs for RT alignment
        n_good_psms = sum(best_hits)

        if n_good_psms < 100
            # Get file name for warning
            file_name = try
                getFileIdToName(getMSData(search_context), ms_file_idx)
            catch
                "file_$ms_file_idx"
            end

            @user_warn "Too few PSMs ($n_good_psms) for RT alignment in file: $file_name (need ≥100). Using identity RT model."

            # Use identity mapping as fallback
            identity_model = IdentityModel()
            setRtIrtMap!(search_context, identity_model, ms_file_idx)
            setIrtRtMap!(search_context, identity_model, ms_file_idx)
            continue
        end

        # Extract best PSM data
        best_rts = psms[:rt][best_hits]
        best_irts = psms[:irt_predicted][best_hits]

        # Fit simple 7-knot UniformSpline for RT → iRT conversion
        rt_to_irt_spline = UniformSpline(
            best_irts,    # y values (iRT)
            best_rts,     # x values (RT)
            3,            # degree (cubic)
            7             # n_knots (fixed)
        )

        # Fit simple 7-knot UniformSpline for iRT → RT conversion (inverse)
        irt_to_rt_spline = UniformSpline(
            best_rts,     # y values (RT)
            best_irts,    # x values (iRT)
            3,            # degree (cubic)
            7             # n_knots (fixed)
        )

        # Store models in SearchContext
        setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
        setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)

        # Optionally generate diagnostic plots
        if params.plot_rt_alignment
            plot_rt_alignment_firstpass(
                best_rts,
                best_irts,
                SplineRtConversionModel(rt_to_irt_spline),
                ms_file_idx,
                getDataOutDir(search_context)
            )
        end

    catch e
        # Get file name for debugging
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end

        # Safely compute PSM count to avoid excessive output
        n_good_psms = try
            sum(best_hits)
        catch
            0  # Default to 0 if calculation fails
        end

        @user_warn "RT mapping failed for MS data file: $file_name ($n_good_psms PSMs). Error: $e. Using identity RT model."

        # Use identity mapping as fallback - no RT to iRT conversion
        identity_model = IdentityModel()
        setRtIrtMap!(search_context, identity_model, ms_file_idx)
        setIrtRtMap!(search_context, identity_model, ms_file_idx)
    end
end
```

**Lines to Replace**: 355-449 (entire try-catch block inside the for loop)

**Changes Made**:
1. ✅ Removed `if params.use_robust_fitting` conditional - single code path
2. ✅ Removed calls to `fit_irt_model()` - no RANSAC, no penalty, no monotonic enforcement
3. ✅ Use simple `UniformSpline` with fixed 7 knots
4. ✅ Changed minimum PSM threshold from implicit to explicit (100 PSMs)
5. ✅ Simplified error handling
6. ✅ Kept plotting functionality
7. ✅ Fit both directions (RT→iRT and iRT→RT) with same simple approach

---

### File 2: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/types.jl`

#### Change 2: Remove `use_robust_fitting` parameter

**Current Code:**
```julia
struct FirstPassSearchParameters <: SearchParameters
    # ... other fields ...
    use_robust_fitting::Bool
    # ... other fields ...
end
```

**New Code:**
```julia
struct FirstPassSearchParameters <: SearchParameters
    # ... other fields ...
    # REMOVED: use_robust_fitting::Bool
    # ... other fields ...
end
```

**Lines to Modify**: Find and remove the `use_robust_fitting::Bool` field

**Note**: You'll also need to update the constructor that reads from PioneerParameters to not extract this parameter.

---

### File 3: Parameter JSON Configuration

#### Change 3: Remove `use_robust_fitting` from default parameters

**Files to Update**:
- Any example/default parameter JSON files that specify `use_robust_fitting`

**Action**: Remove or comment out lines like:
```json
{
  "first_pass_search": {
    "use_robust_fitting": true,
    ...
  }
}
```

This parameter will no longer be used.

---

## Impact Assessment

### What This Changes

1. **FirstPassSearch**: Simplified RT alignment, always uses 7-knot UniformSpline
2. **ParameterTuningSearch**: **NOT CHANGED** - still uses robust `fit_irt_model()` with adaptive knots
3. **SecondPassSearch and beyond**: Use models from FirstPassSearch, so indirectly affected

### What This Doesn't Change

1. **ParameterTuningSearch**: Still uses full `fit_irt_model()` with RANSAC and monotonic enforcement
2. **Coefficient padding fix**: Still in place as safety net
3. **Unit tests**: Comprehensive tests remain
4. **Other search methods**: No changes

### Why This Is Safe

1. **FirstPassSearch is exploratory**: Its RT models are used for initial PSM collection
2. **ParameterTuningSearch refines**: More sophisticated RT models built there
3. **7 knots is proven**: Works well in practice (used in legacy `use_robust_fitting=false` mode)
4. **100 PSM threshold**: Conservative, ensures enough data for stable spline fitting
5. **Identity fallback**: Files with insufficient data get identity mapping (safe default)

---

## Testing Strategy

### Unit Tests

No changes needed - existing UniformSpline tests cover this usage.

### Integration Test

Run SearchDIA on test dataset:
```julia
SearchDIA("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/KYSE70_TEST/KYSE70_TEST_A.json")
```

**Expected Results**:
1. FirstPassSearch completes successfully
2. RT alignment plots generated
3. No dimension mismatch errors
4. Downstream searches use the RT models successfully

### Validation Checks

1. **Compare RT models**: Check that FirstPassSearch RT models are reasonable
2. **Check PSM counts**: Verify sufficient PSMs pass to downstream searches
3. **Check quantification**: Ensure final protein quantification is similar
4. **Performance**: Measure runtime improvement (should be faster)

---

## Migration Path

### For Users with `use_robust_fitting = true` (current default)

**Impact**: Will switch from complex robust fitting to simple 7-knot spline
**Action**:
- Remove `use_robust_fitting` from parameter files (will be ignored)
- Results should be very similar for well-behaved data
- May see slight differences in edge cases (rare)

### For Users with `use_robust_fitting = false`

**Impact**: Minimal - already using simple spline, just changing from 5 to 7 knots
**Action**: Remove parameter from config files

---

## Rollback Plan

If issues arise:

1. **Revert commit** - single commit with all changes
2. **Or restore `use_robust_fitting` logic** - add back the conditional
3. **All changes are in FirstPassSearch** - isolated to one module

---

## Summary of Files Modified

| File | Change | Lines | Complexity |
|------|--------|-------|------------|
| `FirstPassSearch/utils.jl` | Replace RT alignment logic | 355-449 | Major simplification |
| `FirstPassSearch/types.jl` | Remove `use_robust_fitting` field | Find field | Minor |
| Example config JSONs | Remove `use_robust_fitting` | Various | Documentation |

**Total Impact**:
- ~95 lines removed/simplified
- 1 parameter removed
- 1 code path instead of 2
- Significant complexity reduction

---

## Open Questions

1. **Should we add a configurable min_psms threshold?**
   - Current plan: Hardcoded 100
   - Alternative: Make it a parameter (e.g., `min_psms_for_rt_alignment`)

2. **Should we keep plotting functionality?**
   - Current plan: Yes, keep it
   - Alternative: Remove if not used

3. **Should we add a maximum PSM limit for UniformSpline?**
   - Current plan: No limit
   - Alternative: Cap at 10000 PSMs for speed (subsample if more)
   - Rationale: UniformSpline with 7 knots can handle large datasets efficiently

4. **Should 7 knots be configurable?**
   - Current plan: Hardcoded 7
   - Alternative: Add `n_knots_firstpass` parameter
   - Rationale: 7 is a good default, making it configurable adds complexity back
