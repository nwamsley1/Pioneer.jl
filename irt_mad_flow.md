# iRT MAD Flow Through Pioneer.jl Codebase

This document traces how the `irt_mad` (indexed Retention Time Median Absolute Deviation) value flows through the Pioneer.jl codebase and affects downstream operations.

## 1. Calculation in `fit_irt_model()`

**File**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
**Lines**: 344-349

```julia
# Fit initial RT to iRT spline
rt_to_irt_map = UniformSplinePenalized(...)

# Calculate residuals
predicted_irt = [rt_to_irt_map(rt) for rt in psms[!, :rt]]
residuals = psms[!, :irt_predicted] .- predicted_irt

# Remove outliers
irt_mad = mad(residuals, normalize=false)::Float32  # ⚠️ BUG: normalize=false
```

**What it represents**:
- The median absolute deviation of RT alignment residuals
- Measures how well the spline model fits the RT-to-iRT relationship
- Should estimate the standard deviation of alignment errors

**Current bug**: `normalize=false` returns raw MAD instead of σ-estimate

## 2. Return from `fit_irt_model()`

**File**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
**Line**: 389 (spline path) or 405 (linear fallback)

```julia
return (final_model, psms[!, :rt], psms[!, :irt_predicted], irt_mad)
#       ^^^^^^^^^^^^  ^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^^  ^^^^^^^^^
#       RT model      RT values     iRT values              MAD value
```

**Return signature**: `Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}`

## 3. ParameterTuningSearch Receives and Processes

### 3.1 Wrapper Function Call

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`
**Lines**: 253-269

```julia
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}
    # Call shared fit_irt_model from CommonSearchUtils/rt_alignment_utils.jl
    return Pioneer.fit_irt_model(
        psms;
        lambda_penalty = Float32(getRtAlignmentLambdaPenalty(params)),
        ransac_threshold = getRtAlignmentRansacThreshold(params),
        min_psms = getRtAlignmentMinPsms(params),
        spline_degree = getSplineDegree(params),
        max_knots = getSplineNKnots(params),
        outlier_threshold = Float32(getOutlierThreshold(params))
    )
end
```

### 3.2 Called During Convergence Check

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`
**Line**: 374

```julia
rt_model_data = fit_irt_model(params, psms_for_rt_model)
#               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#               Returns: (rt_model, rt_vec, irt_vec, irt_mad)

set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
```

### 3.3 Storage in SearchContext

**File**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`
**Lines**: 148-167

```julia
function set_rt_to_irt_model!(
    ptsr::ParameterTuningSearchResults,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64,
    model::Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}
) where {P<:ParameterTuningSearchParameters}

    # Store RT model components in results
    ptsr.rt_to_irt_model[] = model[1]  # RT conversion model
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])          # RT values (for plotting)
    append!(ptsr.irt, model[3])         # iRT values (for plotting)

    # CRITICAL: Store RT model in SearchContext for downstream methods
    setRtIrtMap!(search_context, model[1], ms_file_idx)

    # CRITICAL: Calculate and store iRT tolerance
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
    #                                            ^^^^^^^   ^^^^^^^^^^^^^^^^^^^
    #                                            irt_mad   sigma_tolerance (from JSON)
end
```

**Key calculation**:
```julia
iRT_tolerance = irt_mad * sigma_tolerance
```

Where:
- `irt_mad`: The 4th element from `fit_irt_model()` return tuple
- `sigma_tolerance`: JSON config parameter (typically 4-5), means "number of standard deviations"

**Example**:
- If `irt_mad = 0.3` (minutes) and `sigma_tolerance = 4`
- Then `iRT_tolerance = 0.3 * 4 = 1.2` minutes

## 4. Storage in SearchContext

**Access method**: `getIrtErrors(search_context)[ms_file_idx]`

This dictionary maps each MS file index to its iRT tolerance value, which represents the RT search window width for that file.

## 5. Downstream Usage: RT Index Creation

### 5.1 FirstPassSearch Creates RT Indices

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`
**Lines**: 186-246

```julia
function create_rt_indices!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    precursor_dict::PrecToIrtType,
    params::FirstPassSearchParameters
)
    # Calculate iRT errors combining peak width and cross-run variation
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)

    # Store in SearchContext
    setIrtErrors!(search_context, irt_errs)

    # ...RT index creation uses these tolerances...
end
```

### 5.2 Error Calculation in FirstPassSearch

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`
**Lines**: 264-306

```julia
function get_irt_errs(
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}},
    prec_to_irt::Dictionary{UInt32, ...},
    params::FirstPassSearchParameters
)
    # Get upper bound on peak FWHM using median + n*MAD
    fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm], fwhms)

    # Get variance in iRT across runs
    variance_ = collect(skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt)))
    irt_std = median(variance_) * params.irt_nstd

    # Combine peak width variation + cross-run variation
    return map(x->Float32((x+irt_std))::Float32, fwhms)::Dictionary{Int64, Float32}
    #                      ^   ^^^^^^^
    #                      |     |
    #                      |     Cross-run RT variation
    #                      Peak width (FWHM-based)
end
```

**Note**: FirstPassSearch uses a **more sophisticated tolerance** that combines:
1. Chromatographic peak width (FWHM + n*MAD)
2. Cross-run iRT variation

This is different from ParameterTuningSearch's simpler `irt_mad * sigma_tolerance` approach.

## 6. RT Index Usage in Search Methods

**File**: `src/Routines/SearchDIA/CommonSearchUtils/buildRTIndex.jl`

The RT indices use these tolerances to:
1. Determine which precursors to search in each RT window
2. Define the width of chromatographic extraction windows
3. Control the trade-off between sensitivity (wider windows) and specificity (narrower windows)

### Example RT Index Query:
```julia
# Get precursors within RT tolerance of current scan
rt_window_start = current_irt - irt_tolerance
rt_window_end = current_irt + irt_tolerance

# Binary search RT index for precursors in window
precursors_to_search = query_rt_index(rt_index, rt_window_start, rt_window_end)
```

## 7. The `normalize=false` Bug

### Current Behavior (WRONG):

```julia
# Line 349 in rt_alignment_utils.jl
irt_mad = mad(residuals, normalize=false)  # Returns raw MAD

# Later in ParameterTuningSearch
iRT_tolerance = raw_MAD * sigma_tolerance
```

**Example**:
- If residuals have `raw_MAD = 0.3` minutes
- And `sigma_tolerance = 4` (meaning 4 standard deviations)
- Then `iRT_tolerance = 0.3 * 4 = 1.2` minutes

### Correct Behavior (SHOULD BE):

```julia
# Should be:
irt_mad = mad(residuals, normalize=true)  # Returns σ-estimate = raw_MAD * 1.4826

# Later in ParameterTuningSearch
iRT_tolerance = σ_estimate * sigma_tolerance
```

**Example**:
- If residuals have `raw_MAD = 0.3` minutes
- Then `σ_estimate = 0.3 * 1.4826 = 0.445` minutes
- And `sigma_tolerance = 4` (meaning 4 standard deviations)
- Then `iRT_tolerance = 0.445 * 4 = 1.78` minutes

### Impact:

**Current windows are ~1.48x too narrow!**
- Current: `iRT_tolerance = 1.2 minutes` (4 × raw MAD)
- Should be: `iRT_tolerance = 1.78 minutes` (4 × 1.4826 × raw MAD)
- Ratio: 1.78 / 1.2 = **1.48x**

This means:
- RT search windows are narrower than intended
- May be missing valid precursor matches
- Lower sensitivity in downstream searches

## 8. Comparison with Other MAD Usage

**Other places correctly use `normalize=true`**:

### ParameterTuningSearch Mass Error (Line 1531):
```julia
mad_error = mad(ppm_errs, normalize=true)  # ✓ Correct
```

### FirstPassSearch MS1 Error (Line 401):
```julia
mad_dev = mad(ms1_errs; normalize=true)    # ✓ Correct
```

### FirstPassSearch FWHM (Line 484):
```julia
mad_fwhm = mad(fwhms, normalize=true)      # ✓ Correct
```

**Conclusion**: The RT alignment code is inconsistent with the rest of the codebase.

## 9. Recommended Fix

### Change Required:

**File**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
**Line**: 349

```julia
# FROM:
irt_mad = mad(residuals, normalize=false)::Float32

# TO:
irt_mad = mad(residuals, normalize=true)::Float32  # Returns σ-estimate
```

### Impact of Fix:
- RT search windows will be **1.48x wider** (as originally intended)
- Better matches the semantics of `sigma_tolerance` parameter
- Consistent with rest of codebase
- May improve sensitivity in downstream searches
- Should have minimal performance impact

### Testing:
After fix, monitor:
1. iRT tolerance values in logs (should increase by ~1.48x)
2. Number of PSMs identified in downstream searches (may increase)
3. False discovery rates (should remain stable with proper FDR control)
4. RT alignment quality (should remain similar or improve)

## 10. Summary Flow Diagram

```
fit_irt_model()
    ↓
Calculate residuals from RT alignment
    ↓
irt_mad = mad(residuals, normalize=false)  ← BUG HERE
    ↓
Return (rt_model, rt_vec, irt_vec, irt_mad)
    ↓
ParameterTuningSearch receives tuple
    ↓
iRT_tolerance = irt_mad * sigma_tolerance (from JSON)
    ↓
Store in SearchContext: getIrtErrors(search_context)[ms_file_idx]
    ↓
FirstPassSearch may override with more sophisticated calculation
    ↓
Used by RT indices to define search windows
    ↓
Controls precursor selection in all downstream searches
```

## 11. Files That Need Changes

1. **Primary fix**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
   - Line 349: Change `normalize=false` to `normalize=true`

2. **No other code changes needed** - the value is just passed through

3. **Documentation updates**:
   - Update comments to clarify this returns σ-estimate
   - Update function docstring to mention normalization

## 12. Historical Context

The `normalize=false` has been present since the extraction of RT fitting code to shared utilities (commit 429b112d). It's unclear if this was:
1. An intentional choice with different semantics in mind
2. An oversight during refactoring
3. A placeholder pending further tuning

Given that:
- JSON parameter is named `sigma_tolerance` (implying standard deviations)
- Other MAD calculations use `normalize=true`
- Standard statistical practice uses normalized MAD for σ-estimation

The correct choice is clearly `normalize=true`.
