# Simplify RT Spline Fitting Plan

## Executive Summary

Simplify RT model fitting in ParameterTuningSearch by:
1. **Removing AIC comparison** - Always use spline (no linear fallback except on error)
2. **Conditional RANSAC** - Based on PSM count threshold (default 1000)
3. **Configurable parameters** - Move hardcoded values to defaultSearchParams.json
4. **Refactor UniformSplinePenalized** - Split RANSAC into separate function
5. **Simplify fit_irt_model** - Reduce from 4 calls to 1 UniformSplinePenalized call

---

## Current State Analysis

### fit_irt_model() Current Logic (lines 392-593)

**Two Code Paths:**
1. **Path 1: n_psms < 1000** (lines 419-540)
   - Compares linear vs spline using AIC
   - 2 UniformSplinePenalized calls (initial + final)
   - Complex 7-tuple return from spline fitting
   - Falls back to linear if spline fails or AIC is worse

2. **Path 2: n_psms >= 1000** (lines 542-591)
   - Spline only (no comparison)
   - 2 UniformSplinePenalized calls (initial + final)
   - Falls back to IdentityModel on error

**Hardcoded Values:**
- `comparison_threshold = 1000` (line 417)
- `λ_penalty = 0.1f0` (line 418)
- RANSAC parameters use function defaults

**Problems:**
- Too many code paths (2 main paths + multiple fallbacks)
- Redundant AIC comparison rarely selects linear
- 4 total UniformSplinePenalized calls (2 per path)
- Complex tuple unpacking for spline results
- Hardcoded thresholds not configurable

---

## Proposed Changes

### 1. New Configuration Parameters

**File:** `assets/example_config/defaultSearchParams.json`

Add to `"parameter_tuning"` section:

```json
{
  "parameter_tuning": {
    "rt_alignment": {
      "lambda_penalty": 0.1,
      "ransac_threshold_psms": 1000,
      "min_psms_for_spline": 10
    }
  }
}
```

**Parameters:**
- `lambda_penalty`: P-spline smoothing penalty (default 0.1)
- `ransac_threshold_psms`: PSM count below which RANSAC is enabled (default 1000)
- `min_psms_for_spline`: Minimum PSMs needed for any RT model (default 10)

### 2. Refactor UniformSplinePenalized

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

**Current Issues:**
- Function is 250+ lines with complex RANSAC logic embedded
- Hard to read and maintain
- RANSAC and standard fitting mixed together

**Proposed Split:**

```julia
# NEW: Separate RANSAC wrapper function
function fit_spline_with_ransac(
    u::Vector{T},
    t::Vector{T},
    degree::Int,
    n_knots::Int,
    λ::T = T(0.1),
    order::Int = 2;
    ransac_iterations::Int = 50,
    ransac_sample_size::Int = 30,
    ransac_threshold_factor::T = T(2.0),
    ransac_seed::Union{Int,Nothing} = nothing
) where {T<:AbstractFloat}
    # RANSAC logic here
    # Calls UniformSplinePenalized internally
end

# SIMPLIFIED: Keep UniformSplinePenalized focused on core fitting
function UniformSplinePenalized(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
    λ::T = T(0.1),
    order::Int = 2
) where {I<:Integer, T<:AbstractFloat}
    # Only core P-spline fitting
    # Remove all RANSAC logic
    # ~100 lines instead of 250+
end
```

**Benefits:**
- Separation of concerns (RANSAC vs P-spline fitting)
- Each function has single responsibility
- Easier to test independently
- More readable code

### 3. Simplify fit_irt_model()

**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

**Current:** 200+ lines with 2 paths, 4 spline calls, AIC comparison
**Proposed:** ~50 lines with 1 path, 1 spline call, no comparison

```julia
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}

    n_psms = nrow(psms)

    # Extract parameters
    λ_penalty = Float32(getRtAlignmentLambdaPenalty(params))
    ransac_threshold = getRtAlignmentRansacThreshold(params)
    min_psms = getRtAlignmentMinPsms(params)

    # Calculate knots adaptively
    n_knots = min(max(3, Int(floor(n_psms / 100))), getSplineNKnots(params))

    # Early exit for insufficient data
    if n_psms < min_psms
        @debug_l2 "Too few PSMs ($n_psms) for RT alignment (need ≥$min_psms), using identity model"
        return (IdentityModel(), Float32[], Float32[], 0.0f0)
    end

    @debug_l2 "Using $n_knots knots for RT spline ($n_psms PSMs)"

    # Single unified spline fitting path with conditional RANSAC
    try
        # Determine if we should use RANSAC based on PSM count
        use_ransac = n_psms < ransac_threshold

        @debug_l2 if use_ransac
            "Limited data ($n_psms PSMs): Using RANSAC + penalty (λ=$λ_penalty)"
        else
            "Abundant data ($n_psms PSMs): Using standard fitting (no RANSAC, λ=$λ_penalty)"
        end

        # Fit spline (with or without RANSAC)
        rt_to_irt_map = if use_ransac
            # Hardcoded RANSAC parameters (not configurable)
            fit_spline_with_ransac(
                psms[!, :irt_predicted],
                psms[!, :rt],
                getSplineDegree(params),
                n_knots,
                λ_penalty,
                2,  # 2nd order penalty
                ransac_iterations=50,
                ransac_sample_size=30,
                ransac_threshold_factor=Float32(2.0)
            )
        else
            # Standard P-spline fitting
            UniformSplinePenalized(
                psms[!, :irt_predicted],
                psms[!, :rt],
                getSplineDegree(params),
                n_knots,
                λ_penalty,
                2  # 2nd order penalty
            )
        end

        # Calculate residuals
        predicted_irt = [rt_to_irt_map(rt) for rt in psms[!, :rt]]
        residuals = psms[!, :irt_predicted] .- predicted_irt

        # Remove outliers
        irt_mad = mad(residuals, normalize=false)::Float32
        valid_mask = abs.(residuals) .< (irt_mad * getOutlierThreshold(params))
        valid_psms = psms[valid_mask, :]

        # Check if enough PSMs remain
        if nrow(valid_psms) < min_psms
            @debug_l2 "Too few PSMs after outlier removal ($(nrow(valid_psms))), using identity model"
            return (IdentityModel(), Float32[], Float32[], 0.0f0)
        end

        # Refit with filtered PSMs (same method as initial fit)
        n_knots_final = min(n_knots, max(3, Int(floor(nrow(valid_psms) / 100))))

        final_map = if use_ransac
            fit_spline_with_ransac(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                getSplineDegree(params),
                n_knots_final,
                λ_penalty,
                2,
                ransac_iterations=50,
                ransac_sample_size=30,
                ransac_threshold_factor=Float32(2.0)
            )
        else
            UniformSplinePenalized(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                getSplineDegree(params),
                n_knots_final,
                λ_penalty,
                2
            )
        end

        final_model = SplineRtConversionModel(final_map)

        return (final_model, valid_psms[!, :rt], valid_psms[!, :irt_predicted], irt_mad)

    catch e
        # Only fallback: use linear model on error
        @user_warn "RT spline fitting failed ($e), falling back to linear model"

        linear_model, linear_std, _ = fit_linear_irt_model(psms)

        # Calculate MAD for linear model
        predicted = [linear_model(rt) for rt in psms[!, :rt]]
        residuals = psms[!, :irt_predicted] .- predicted
        linear_mad = median(abs.(residuals .- median(residuals)))

        return (
            linear_model,
            psms[!, :rt],
            psms[!, :irt_predicted],
            Float32(linear_mad)
        )
    end
end
```

**Removed:**
- ❌ AIC calculation (calculate_aic function calls)
- ❌ Model comparison logic (linear vs spline)
- ❌ Separate code paths for <1000 and >=1000 PSMs
- ❌ Complex 7-tuple unpacking from spline fitting
- ❌ min_psms_for_spline variable (now parameter)

**Added:**
- ✅ Single unified code path
- ✅ Configurable parameters from JSON
- ✅ Clear RANSAC decision logic
- ✅ Linear model only as error fallback
- ✅ Simplified logic (~50 lines vs 200+)

### 4. Add Parameter Accessors

**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl` or similar

```julia
# New accessor functions for RT alignment parameters
function getRtAlignmentLambdaPenalty(params::ParameterTuningSearchParameters)
    return get(params, "parameter_tuning", "rt_alignment", "lambda_penalty", 0.1)
end

function getRtAlignmentRansacThreshold(params::ParameterTuningSearchParameters)
    return get(params, "parameter_tuning", "rt_alignment", "ransac_threshold_psms", 1000)
end

function getRtAlignmentMinPsms(params::ParameterTuningSearchParameters)
    return get(params, "parameter_tuning", "rt_alignment", "min_psms_for_spline", 10)
end
```

---

## Implementation Steps

### Step 1: Add Configuration Parameters
**File:** `assets/example_config/defaultSearchParams.json`

1. Add `rt_alignment` section under `parameter_tuning`
2. Set defaults: `lambda_penalty=0.1`, `ransac_threshold_psms=1000`, `min_psms_for_spline=10`

### Step 2: Add Parameter Accessors
**File:** Appropriate location in ParameterTuningSearch module

1. Add `getRtAlignmentLambdaPenalty()`
2. Add `getRtAlignmentRansacThreshold()`
3. Add `getRtAlignmentMinPsms()`

### Step 3: Refactor UniformSplinePenalized
**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

1. Extract RANSAC logic into new `fit_spline_with_ransac()` function
2. Simplify `UniformSplinePenalized()` to core P-spline fitting only
3. Update docstrings for both functions
4. Ensure type conversions are correct (existing fix for Float32/Float64)

### Step 4: Simplify fit_irt_model()
**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

1. Remove entire AIC comparison section (lines 416-540)
2. Remove `calculate_aic()` function if not used elsewhere
3. Implement single unified spline fitting path
4. Use `fit_spline_with_ransac()` when `n_psms < threshold`
5. Use `UniformSplinePenalized()` when `n_psms >= threshold`
6. Keep linear model only as error fallback (catch block)
7. Reduce from 4 UniformSplinePenalized calls to 2 (initial + refit)

### Step 5: Update Call Sites
**Files:** Any other files calling `fit_irt_model()`

1. Verify return type still matches expected signature
2. Update any code expecting AIC-related return values

### Step 6: Testing
1. Test with < 1000 PSMs → should use RANSAC
2. Test with >= 1000 PSMs → should use standard fitting
3. Test with error conditions → should fall back to linear
4. Verify RT alignment quality is maintained
5. Measure performance improvement

---

## Expected Outcomes

### Code Quality
- **Reduced complexity**: ~150 lines removed from fit_irt_model
- **Better organization**: RANSAC separated from core fitting
- **Easier maintenance**: Clear separation of concerns
- **More testable**: Each function has single responsibility

### Performance
- **< 1000 PSMs**: No change (~50-100ms with RANSAC)
- **>= 1000 PSMs**: No change (~2ms standard fitting)
- **Overall**: Slightly faster due to removed AIC overhead

### Configurability
- **Lambda penalty**: User can adjust smoothing strength
- **RANSAC threshold**: User can tune when robust fitting is used
- **Min PSMs**: User can set minimum data requirements

### Maintainability
- **Clearer logic**: Single path easier to understand
- **Fewer edge cases**: Removed complex AIC comparison
- **Better documentation**: Each function has focused purpose

---

## Migration Notes

### Breaking Changes
- ❌ None - Return signature of fit_irt_model remains the same

### Behavior Changes
- ✅ Always prefers spline over linear (unless error)
- ✅ Linear model only used as error fallback
- ✅ More consistent RT alignment (no AIC switching)

### Configuration Changes
- ✅ New parameters in defaultSearchParams.json
- ✅ Users can now tune lambda penalty
- ✅ Users can adjust RANSAC threshold

---

## Risks and Mitigation

### Risk 1: Spline may fail where linear previously succeeded
**Mitigation**: Keep linear model as error fallback

### Risk 2: Users may have optimized for current AIC behavior
**Mitigation**: Document change in release notes, provide migration guide

### Risk 3: Refactoring may introduce bugs
**Mitigation**: Comprehensive testing with various PSM counts

---

## Timeline Estimate

| Task | Time | Cumulative |
|------|------|------------|
| Add config parameters | 10 min | 10 min |
| Add parameter accessors | 10 min | 20 min |
| Extract fit_spline_with_ransac | 30 min | 50 min |
| Simplify UniformSplinePenalized | 20 min | 70 min |
| Rewrite fit_irt_model | 40 min | 110 min |
| Update call sites | 10 min | 120 min |
| Testing and validation | 30 min | 150 min |
| **Total** | **~2.5 hours** | |

---

## Commit Strategy

**Commit 1:** Add configuration parameters
```bash
git commit -m "Add RT alignment configuration parameters

- Add parameter_tuning.rt_alignment section to defaultSearchParams.json
- Add lambda_penalty, ransac_threshold_psms, min_psms_for_spline
- Add parameter accessor functions"
```

**Commit 2:** Refactor spline fitting functions
```bash
git commit -m "Refactor UniformSplinePenalized: extract RANSAC into separate function

- Create fit_spline_with_ransac() wrapper function
- Simplify UniformSplinePenalized() to core P-spline fitting
- Reduce complexity and improve code organization
- Update docstrings for both functions"
```

**Commit 3:** Simplify fit_irt_model
```bash
git commit -m "Simplify RT model fitting: remove AIC comparison, single code path

- Remove AIC comparison between linear and spline models
- Implement single unified spline fitting path
- Use RANSAC conditionally based on PSM count threshold
- Linear model only as error fallback
- Reduce from 200+ lines to ~50 lines
- Reduce from 4 spline calls to 2 (initial + refit)

Breaking changes: None (return signature unchanged)
Behavior changes: Always prefer spline, more consistent RT alignment"
```

---

## References

- **Current RANSAC implementation**: src/utils/ML/uniformBasisCubicSpline.jl:490-633
- **Current fit_irt_model**: src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:382-593
- **Parameter structure**: assets/example_config/defaultSearchParams.json
- **RANSAC plan**: ransac_spline_implementation_plan.md
