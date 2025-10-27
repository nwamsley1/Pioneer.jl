# Plan: AIC-Based Model Selection for RT Alignment with Limited Data

## Overview

When fitting retention time (RT) alignment models with limited data (< 1000 PSMs), we will compare two models and use AIC (Akaike Information Criterion) to automatically select the better one:

1. **UniformSpline** - Current cubic (degree 3) spline approach
2. **Simple Linear Model** - `irt_predicted = intercept + slope * rt`

AIC balances goodness-of-fit against model complexity, penalizing models with more parameters. **Lower AIC = better model.**

## Motivation

### Problem
- With few PSMs (< 1000), cubic splines may overfit
- Splines have many parameters (4 * n_knots coefficients)
- Linear models are more robust and stable with limited data
- Need principled way to choose between models

### Solution
- Fit both models when n_psms < 1000
- Calculate AIC for each
- Automatically select the model with lower AIC
- Transparent logging of model selection

### Benefits
- **Robustness**: Linear model prevents overfitting with limited data
- **Adaptability**: Spline model captures nonlinearity when justified
- **Principled**: AIC provides statistical basis for selection
- **Performance**: Linear fitting is O(n), much faster than spline O(n²)

## When to Apply

Model comparison will be triggered when:
- `n_psms < 1000` (configurable threshold)
- RT alignment is being fitted (not falling back to IdentityModel)

For datasets with ≥ 1000 PSMs, continue using UniformSpline only (current behavior).

## Implementation Plan

### 1. Create Linear Model Type

**Location**: `src/structs/RetentionTimeConversionModel.jl`

Add new model type that implements the RtConversionModel interface:

```julia
struct LinearRtConversionModel <: RtConversionModel
    slope::Float32
    intercept::Float32
end

function (m::LinearRtConversionModel)(rt::AbstractFloat)
    return m.intercept + m.slope * rt
end
```

**Properties:**
- Implements same callable interface as other RT models: `model(rt)`
- 2 parameters: slope and intercept
- Simple, interpretable, fast evaluation

### 2. Implement Linear Model Fitting

**Location**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

Add function using closed-form least squares solution:

```julia
"""
    fit_linear_irt_model(psms::DataFrame)

Fit simple linear model: irt_predicted = intercept + slope * rt

Uses closed-form least squares solution via normal equations.

# Arguments
- `psms`: DataFrame with :rt and :irt_predicted columns

# Returns
Tuple containing:
- `model`: LinearRtConversionModel with fitted parameters
- `residual_std`: Standard deviation of residuals
- `n_params`: Number of parameters (always 2)
"""
function fit_linear_irt_model(
    psms::DataFrame
)::Tuple{LinearRtConversionModel, Float32, Int}

    rt = psms[!, :rt]
    irt = psms[!, :irt_predicted]
    n = length(rt)

    # Least squares: irt = intercept + slope * rt
    # Normal equations: [n, Σrt; Σrt, Σrt²] * [intercept; slope] = [Σirt; Σ(rt*irt)]
    sum_rt = sum(rt)
    sum_rt2 = sum(rt .^ 2)
    sum_irt = sum(irt)
    sum_rt_irt = sum(rt .* irt)

    # Solve 2x2 system
    det = n * sum_rt2 - sum_rt^2

    if abs(det) < 1e-10
        # Degenerate case - all RTs identical
        return (
            LinearRtConversionModel(0.0f0, Float32(mean(irt))),
            Float32(std(irt)),
            2
        )
    end

    intercept = (sum_rt2 * sum_irt - sum_rt * sum_rt_irt) / det
    slope = (n * sum_rt_irt - sum_rt * sum_irt) / det

    # Calculate residuals
    predicted = intercept .+ slope .* rt
    residuals = irt .- predicted
    residual_std = Float32(std(residuals))

    model = LinearRtConversionModel(Float32(slope), Float32(intercept))

    return (model, residual_std, 2)
end
```

**Algorithm Details:**
- Normal equations: Closed-form solution for linear regression
- O(n) complexity (just sums and arithmetic)
- Handles degenerate case (all RTs same)
- Returns model, residual std, and parameter count

### 3. Implement AIC Calculation

**Location**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

Add helper function:

```julia
"""
    calculate_aic(n_observations::Int, n_parameters::Int, residual_variance::Float32)

Calculate Akaike Information Criterion (AIC) for model selection.

AIC = n * log(σ²) + 2k

where:
- n = number of observations
- σ² = residual variance (mean squared error)
- k = number of parameters

Lower AIC indicates better model (balances fit quality vs complexity).

# Arguments
- `n_observations`: Number of data points
- `n_parameters`: Number of model parameters
- `residual_variance`: Variance of residuals (MSE)

# Returns
AIC value (Float32). Lower values indicate better models.
"""
function calculate_aic(
    n_observations::Int,
    n_parameters::Int,
    residual_variance::Float32
)::Float32

    if residual_variance <= 0
        return Inf32
    end

    log_likelihood_term = n_observations * log(residual_variance)
    penalty_term = 2 * n_parameters

    return Float32(log_likelihood_term + penalty_term)
end
```

**AIC Formula:**
- `AIC = n * log(σ²) + 2k`
- First term: Log-likelihood (lower residual variance = better fit)
- Second term: Penalty for complexity (more parameters = higher penalty)
- Lower AIC wins

**Why AIC?**
- Balances fit quality against model complexity
- Penalizes overfitting
- Well-established in statistics
- Simple to compute

### 4. Modify `fit_irt_model` for Model Comparison

**Location**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

**Current function signature:**
```julia
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}
```

**Modification strategy:**
Add comparison logic at the beginning, after checking minimum PSM count.

**Pseudocode:**
```
function fit_irt_model(params, psms)
    n_psms = nrow(psms)

    # Existing checks for minimum PSMs
    if n_psms < 10
        return IdentityModel()
    end

    # NEW: Model comparison for limited data
    comparison_threshold = 1000

    if n_psms < comparison_threshold
        # Fit both models and compare

        # 1. Try spline model
        spline_model = try_fit_spline(psms, params)

        # 2. Fit linear model (always succeeds)
        linear_model, linear_std, linear_n_params = fit_linear_irt_model(psms)

        # 3. Calculate AIC for both
        if spline_model succeeded
            spline_aic = calculate_aic_for_spline(spline_model)
        else
            spline_aic = Inf32
        end

        linear_variance = linear_std^2
        linear_aic = calculate_aic(n_psms, linear_n_params, linear_variance)

        # 4. Select better model
        if linear_aic < spline_aic
            log "Selected linear model"
            return prepare_linear_results(linear_model, psms)
        else
            log "Selected spline model"
            return prepare_spline_results(spline_model, psms)
        end
    else
        # n_psms >= 1000: Use spline only (existing behavior)
        return fit_spline_as_before(psms, params)
    end
end
```

**Detailed implementation:**

```julia
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}

    n_psms = nrow(psms)

    # Calculate knots for spline (needed for both paths)
    max_knots_from_psms = max(2, floor(Int, n_psms / 5))
    configured_knots = getSplineNKnots(params)
    n_knots = min(configured_knots, max_knots_from_psms)

    # Check minimum PSM requirement
    if n_psms < 10
        @debug_l2 "Too few PSMs ($n_psms) for RT alignment (need ≥10), using identity model"
        return (IdentityModel(), Float32[], Float32[], 0.0f0)
    end

    # Log knot reduction if needed
    if n_knots < configured_knots
        @debug_l1 "Reducing RT spline knots from $configured_knots to $n_knots due to limited PSMs ($n_psms)"
    end

    # MODEL COMPARISON PATH: n_psms < 1000
    comparison_threshold = 1000

    if n_psms < comparison_threshold

        # 1. TRY SPLINE MODEL
        spline_result = try
            # Fit initial spline
            rt_to_irt_map = UniformSpline(
                psms[!, :irt_predicted],
                psms[!, :rt],
                getSplineDegree(params),
                n_knots
            )

            # Calculate residuals
            predicted_irt = [rt_to_irt_map(rt) for rt in psms[!, :rt]]
            residuals = psms[!, :irt_predicted] .- predicted_irt

            # Remove outliers
            irt_mad = mad(residuals, normalize=false)::Float32
            valid_mask = abs.(residuals) .< (irt_mad * getOutlierThreshold(params))
            valid_psms = psms[valid_mask, :]

            # Check if enough PSMs remain
            if nrow(valid_psms) < 10
                (nothing, nothing, nothing, nothing)
            else
                # Refit with valid PSMs
                n_valid = nrow(valid_psms)
                max_knots_final = max(2, floor(Int, n_valid / 5))
                n_knots_final = min(n_knots, max_knots_final)

                final_spline = UniformSpline(
                    valid_psms[!, :irt_predicted],
                    valid_psms[!, :rt],
                    getSplineDegree(params),
                    n_knots_final
                )

                # Calculate final residuals for AIC
                final_pred = [final_spline(rt) for rt in valid_psms[!, :rt]]
                final_resid = valid_psms[!, :irt_predicted] .- final_pred
                final_variance = var(final_resid)

                (
                    SplineRtConversionModel(final_spline),
                    valid_psms[!, :rt],
                    valid_psms[!, :irt_predicted],
                    irt_mad,
                    nrow(valid_psms),
                    n_knots_final * 4,  # 4 coefficients per knot for degree 3
                    Float32(final_variance)
                )
            end
        catch e
            @debug_l1 "Spline fitting failed: $e"
            (nothing, nothing, nothing, nothing, nothing, nothing, nothing)
        end

        spline_model, spline_rt, spline_irt, spline_mad,
            spline_n_obs, spline_n_params, spline_variance = spline_result

        # 2. FIT LINEAR MODEL (always succeeds)
        linear_model, linear_std, linear_n_params = fit_linear_irt_model(psms)
        linear_variance = linear_std^2

        # 3. CALCULATE AIC FOR BOTH
        spline_aic = if spline_model !== nothing
            calculate_aic(spline_n_obs, spline_n_params, spline_variance)
        else
            Inf32  # Failed to fit
        end

        linear_aic = calculate_aic(n_psms, linear_n_params, linear_variance)

        # 4. SELECT BETTER MODEL
        @debug_l1 "RT model comparison (n_psms=$n_psms): Linear AIC=$linear_aic, Spline AIC=$spline_aic"

        if linear_aic < spline_aic || spline_model === nothing
            # LINEAR MODEL WINS
            @debug_l1 "Selected linear RT model (AIC: $(round(linear_aic, digits=2)) vs $(round(spline_aic, digits=2)))"

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
        else
            # SPLINE MODEL WINS
            @debug_l1 "Selected spline RT model (AIC: $(round(spline_aic, digits=2)) vs $(round(linear_aic, digits=2)))"

            return (spline_model, spline_rt, spline_irt, spline_mad)
        end

    else
        # STANDARD PATH: n_psms >= 1000, use spline only
        try
            # Existing spline fitting logic (unchanged)
            rt_to_irt_map = UniformSpline(
                psms[!, :irt_predicted],
                psms[!, :rt],
                getSplineDegree(params),
                n_knots
            )

            psms[!, :irt_observed] = rt_to_irt_map.(psms.rt::Vector{Float32})
            residuals = psms[!, :irt_observed] .- psms[!, :irt_predicted]
            irt_mad = mad(residuals, normalize=false)::Float32

            valid_psms = psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :]

            n_valid_psms = nrow(valid_psms)
            if n_valid_psms < 10
                @debug_l2 "Too few PSMs after outlier removal ($n_valid_psms), using identity model"
                return (IdentityModel(), Float32[], Float32[], 0.0f0)
            end

            max_knots_final = max(2, floor(Int, n_valid_psms / 5))
            n_knots_final = min(n_knots, max_knots_final)

            final_model = SplineRtConversionModel(UniformSpline(
                valid_psms[!, :irt_predicted],
                valid_psms[!, :rt],
                getSplineDegree(params),
                n_knots_final
            ))

            return (final_model, valid_psms[!, :rt], valid_psms[!, :irt_predicted], irt_mad)

        catch e
            @debug_l1 "RT spline fitting failed, using identity model exception=$e"
            return (IdentityModel(), Float32[], Float32[], 0.0f0)
        end
    end
end
```

**Key implementation details:**
- Two distinct code paths: comparison (< 1000 PSMs) vs standard (≥ 1000 PSMs)
- Spline fitting wrapped in try-catch (may fail with limited data)
- Linear model fitting always succeeds (closed-form solution)
- AIC calculated on final fitted data (after outlier removal for spline)
- Clear logging of model selection with AIC values
- Return format identical for both models (type compatibility via RtConversionModel)

### 5. Update Type Annotations

**Location**: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl`

**Current type:**
```julia
struct IterationState
    rt_to_irt_map::SplineRtConversionModel
    # ... other fields
end
```

**Change to:**
```julia
struct IterationState
    rt_to_irt_map::RtConversionModel  # Parent type
    # ... other fields
end
```

**Rationale:**
- `RtConversionModel` is the abstract parent type
- `SplineRtConversionModel` and `LinearRtConversionModel` both inherit from it
- Both implement the same callable interface: `model(rt)`
- Change allows either model type to be stored
- No other code changes needed (polymorphism handles the rest)

## Testing Strategy

### 1. Unit Tests

**File**: `test/UnitTests/test_linear_rt_model.jl` (new file)

```julia
using Test
using DataFrames
using Pioneer

@testset "Linear RT Model Fitting" begin
    # Test 1: Perfect linear relationship
    rt = Float32[1.0, 2.0, 3.0, 4.0, 5.0]
    irt = Float32[10.0, 20.0, 30.0, 40.0, 50.0]  # irt = 10 * rt
    psms = DataFrame(rt=rt, irt_predicted=irt)

    model, residual_std, n_params = fit_linear_irt_model(psms)

    @test model.slope ≈ 10.0 atol=1e-4
    @test model.intercept ≈ 0.0 atol=1e-4
    @test n_params == 2
    @test residual_std < 1e-4  # Nearly perfect fit

    # Test 2: With intercept
    irt2 = Float32[15.0, 25.0, 35.0, 45.0, 55.0]  # irt = 10 * rt + 5
    psms2 = DataFrame(rt=rt, irt_predicted=irt2)

    model2, _, _ = fit_linear_irt_model(psms2)

    @test model2.slope ≈ 10.0 atol=1e-4
    @test model2.intercept ≈ 5.0 atol=1e-4

    # Test 3: Model evaluation
    @test model(3.0f0) ≈ 30.0 atol=1e-4
    @test model2(3.0f0) ≈ 35.0 atol=1e-4

    # Test 4: Degenerate case (all RTs same)
    rt_same = Float32[5.0, 5.0, 5.0, 5.0]
    irt_same = Float32[10.0, 12.0, 11.0, 13.0]
    psms_same = DataFrame(rt=rt_same, irt_predicted=irt_same)

    model_same, _, _ = fit_linear_irt_model(psms_same)
    @test model_same.slope == 0.0
    @test model_same.intercept ≈ mean(irt_same) atol=1e-4
end

@testset "AIC Calculation" begin
    # Test 1: More parameters increase AIC
    aic_simple = calculate_aic(100, 2, 1.0f0)
    aic_complex = calculate_aic(100, 20, 1.0f0)

    @test aic_complex > aic_simple  # More parameters = higher AIC

    # Test 2: Better fit (lower variance) decreases AIC
    aic_poor_fit = calculate_aic(100, 2, 10.0f0)
    aic_good_fit = calculate_aic(100, 2, 1.0f0)

    @test aic_good_fit < aic_poor_fit  # Lower variance = lower AIC

    # Test 3: Invalid variance returns Inf
    aic_invalid = calculate_aic(100, 2, 0.0f0)
    @test aic_invalid == Inf32

    # Test 4: Penalty term scaling
    # AIC = n*log(σ²) + 2k
    # Increasing k by 1 adds penalty of 2
    aic1 = calculate_aic(100, 5, 1.0f0)
    aic2 = calculate_aic(100, 6, 1.0f0)
    @test aic2 - aic1 ≈ 2.0 atol=1e-4
end

@testset "Model Type Interface" begin
    # Test that LinearRtConversionModel implements RtConversionModel interface
    model = LinearRtConversionModel(2.0f0, 5.0f0)

    @test model isa RtConversionModel
    @test model(10.0f0) ≈ 25.0f0  # 2*10 + 5
    @test model(0.0f0) ≈ 5.0f0    # intercept
end
```

### 2. Integration Tests

**Scenarios to test:**

**Test 1: Small dataset (< 500 PSMs) - Expect linear**
```julia
# Generate nearly linear RT data with small sample
n = 100
rt = sort(rand(Float32, n) .* 100)
irt = 2.0f0 .* rt .+ 10.0f0 .+ randn(Float32, n) .* 0.5f0
psms = DataFrame(rt=rt, irt_predicted=irt)

model, _, _, _ = fit_irt_model(params, psms)
@test model isa LinearRtConversionModel
```

**Test 2: Medium dataset (500-1000 PSMs) - Depends on linearity**
```julia
# Test with both linear and nonlinear data
# Expect linear model for linear data
# Expect spline model for highly nonlinear data
```

**Test 3: Large dataset (> 1000 PSMs) - Always spline**
```julia
n = 2000
rt = sort(rand(Float32, n) .* 100)
irt = 2.0f0 .* rt .+ 10.0f0 .+ randn(Float32, n)
psms = DataFrame(rt=rt, irt_predicted=irt)

model, _, _, _ = fit_irt_model(params, psms)
@test model isa SplineRtConversionModel  # Should always use spline
```

**Test 4: Spline fitting failure - Fallback to linear**
```julia
# Create pathological case where spline fails
# Verify linear model is selected
```

### 3. Visual Validation

**QC Plots to add:**
- Scatter plot: RT vs iRT with both linear and spline fits overlaid
- Residual plots for both models
- Text annotation showing AIC values and selected model
- Save plots to output folder

## Edge Cases

| Case | Behavior | Rationale |
|------|----------|-----------|
| All RTs identical | Linear model returns mean(irt), slope=0 | Closed-form solution handles degenerate case |
| Perfect linear data | Linear model wins (AIC near -Inf) | Zero residual variance → very low AIC |
| Highly nonlinear | Spline wins despite penalty | Better fit overcomes parameter penalty |
| Very few PSMs (< 10) | IdentityModel fallback | Neither model reliable |
| Spline fitting fails | Auto-select linear model | spline_aic = Inf32 |
| n_psms = 1000 exactly | Use standard spline path | Comparison only for < 1000 |

## Performance Considerations

**Linear model:**
- Fitting: O(n) - just sums
- Memory: ~8 bytes (2 Float32s)
- Evaluation: 2 FLOPs per call

**Spline model:**
- Fitting: O(n²) - matrix operations
- Memory: ~16 * n_knots bytes
- Evaluation: ~6-7 FLOPs per call

**Overhead:**
- Fitting linear model adds ~0.1-1ms for typical datasets
- AIC calculation is negligible
- Total overhead < 5% for n_psms < 1000

## Expected Outcomes

**Small datasets (< 200 PSMs):**
- Linear model selected ~80% of the time
- More stable predictions
- Faster fitting

**Medium datasets (200-1000 PSMs):**
- Mixed results depending on nonlinearity
- Linear for approximately linear RT
- Spline when RT shift is clearly nonlinear

**Large datasets (≥ 1000 PSMs):**
- Always use spline (no comparison)
- Current behavior preserved
- No performance regression

**Quality metrics:**
- Reduced RT prediction errors for low-PSM files
- Fewer IdentityModel fallbacks
- Improved downstream quantification

## Files to Modify

| File | Changes | Lines |
|------|---------|-------|
| `src/structs/RetentionTimeConversionModel.jl` | Add LinearRtConversionModel struct | +7 |
| `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl` | Add fit_linear_irt_model() | +50 |
| | Add calculate_aic() | +25 |
| | Modify fit_irt_model() | +100 |
| `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl` | Change rt_to_irt_map type | ~1 |
| `test/UnitTests/test_linear_rt_model.jl` | New test file | +150 |
| **Total** | | **~333 lines** |

## Configuration (Optional Future Enhancement)

Could add to `tuning_params.json`:
```json
{
  "parameter_tuning": {
    "rt_model_selection": {
      "enable_comparison": true,
      "comparison_threshold": 1000,
      "force_linear_below": 100,  // Auto-select linear if < 100 PSMs
      "aic_logging": true
    }
  }
}
```

## Open Questions

### 1. AIC Variant for Small Samples?
**Question:** Should we use AICc (corrected AIC) for small samples?

**AICc formula:** `AICc = AIC + 2k(k+1)/(n-k-1)`

**When to use:** Generally recommended when n/k < 40

**Analysis:**
- Linear model: k=2, so use AICc when n < 80
- Spline model: k=4*n_knots, e.g., k=12 for 3 knots, use AICc when n < 480

**Recommendation:** Start with standard AIC for simplicity. Can add AICc later if needed.

### 2. Outlier Handling
**Question:** Should outliers be removed before or after AIC calculation?

**Current plan:**
- Spline: Remove outliers, refit, calculate AIC on cleaned data
- Linear: No outlier removal, calculate AIC on all data

**Alternatives:**
1. Remove outliers from both models before AIC
2. Calculate AIC on all data for both models
3. Use robust regression (e.g., Huber loss) for linear model

**Recommendation:** Keep current plan - spline is sensitive to outliers so remove them, linear is more robust so keep all data. If linear model is selected, outliers have less impact anyway.

### 3. Logging Verbosity
**Question:** How much should we log about model selection?

**Current plan:**
- @debug_l1: Model comparison (always shown)
- @debug_l1: Selected model with AIC values
- Could add to summary statistics

**Recommendation:** Keep logging at @debug_l1 level. Add model selection info to RT alignment plots.

## Implementation Order

1. ✅ Plan created
2. Add LinearRtConversionModel struct
3. Implement fit_linear_irt_model() with unit tests
4. Implement calculate_aic() with unit tests
5. Modify fit_irt_model() to add comparison logic
6. Update type annotations in types.jl
7. Add integration tests
8. Test with real data (various PSM counts)
9. Add model comparison info to QC plots
10. Documentation and commit

## Risks and Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Linear model insufficient for nonlinear RT | Poor alignment quality | AIC will select spline if nonlinearity significant |
| Type system changes break downstream | Compilation errors | RtConversionModel is parent type, interface unchanged |
| Performance regression | Slower processing | Only adds ~0.1-1ms for n_psms < 1000 |
| AIC comparison complexity | Maintenance burden | Well-documented, only triggered for limited data |
| Incorrect AIC calculation | Wrong model selected | Unit tests validate AIC properties |

## Success Criteria

- [ ] All unit tests pass
- [ ] Integration tests pass with real data
- [ ] No performance regression for n_psms ≥ 1000
- [ ] Linear model selected for simple linear RT shifts
- [ ] Spline model selected for complex nonlinear RT shifts
- [ ] Clear logging shows model selection rationale
- [ ] QC plots show both models when comparison performed
- [ ] No breaking changes to existing code
