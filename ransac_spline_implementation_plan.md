# RANSAC Spline Fitting Implementation Plan

## Executive Summary

**RANSAC (Random Sample Consensus)** is a robust fitting algorithm that can handle datasets with up to 50% outlier contamination. This document provides a complete implementation plan for adding RANSAC as an alternative fitting mode for RT→iRT spline alignment.

**Key Decision:** RANSAC is more powerful but slower than the current two-pass method. Implement only if current approach proves insufficient.

---

## Background: RANSAC Algorithm

### Core Concept

RANSAC finds the best model by:
1. **Randomly sampling** minimal subsets of data
2. **Fitting models** to each sample
3. **Counting inliers** (points that fit the model well)
4. **Selecting** the model with most inliers
5. **Refitting** using all inliers from best model

### Why RANSAC?

**Advantages over standard fitting:**
- ✅ Extremely robust to outliers (up to 50% contamination)
- ✅ Provable probabilistic guarantees
- ✅ Well-established in computer vision and robotics
- ✅ No assumptions about outlier distribution

**Compared to our current two-pass method:**
- More robust (50% vs 10-15% outlier tolerance)
- Slower (50 fits vs 2 fits)
- Non-deterministic (random sampling)
- More parameters to tune

---

## Mathematical Foundation

### Required Iterations

To guarantee finding a good sample with probability p:

```
k = log(1 - p) / log(1 - (1 - ε)^s)
```

Where:
- **p**: Desired confidence (typically 0.99)
- **ε**: Expected outlier ratio (e.g., 0.2 = 20% outliers)
- **s**: Sample size

**Examples:**
- 99% confidence, 10% outliers, sample=30 → k ≈ 16 iterations
- 99% confidence, 20% outliers, sample=30 → k ≈ 35 iterations
- 99% confidence, 30% outliers, sample=30 → k ≈ 60 iterations

### Inlier Threshold

**Adaptive threshold based on MAD:**
```
threshold = MAD(initial_residuals) × scale_factor
```

Where:
- MAD = Median Absolute Deviation
- scale_factor = 1.5-3.0 (typically 2.0)

**Alternative: Chi-squared based:**
```
threshold = σ × sqrt(χ²(p, df))
```
- For 95% confidence, df=1: threshold ≈ 1.96σ

---

## Implementation Design

### Fitting Mode Architecture

Add RANSAC as a third mode alongside current methods:

```julia
@enum FittingMode begin
    Standard      # Single fit, no outlier handling
    RobustRefit   # Current two-pass method (default)
    RANSAC        # Random sample consensus
end
```

### Function Signature Extension

```julia
function UniformSplinePenalized(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
    λ::T = T(0.1),
    order::Int = 2,
    fitting_mode::Symbol = :robust_refit,
    # RANSAC-specific parameters
    ransac_iterations::Int = 50,
    ransac_sample_size::Int = 30,
    ransac_threshold_factor::T = T(2.0),
    ransac_seed::Union{Int,Nothing} = nothing  # For reproducibility
) where {I<:Integer, T<:AbstractFloat}
```

### Parameter Defaults

**Sample size:**
- Minimum: n_coeffs (8 for 5 knots)
- Recommended: 20-50 PSMs
- Too small → unstable fits
- Too large → less robust (approaches full dataset)

**Iterations:**
- Calculated adaptively based on expected outlier ratio
- Minimum: 20
- Typical: 50-100
- Can terminate early if very good model found

**Threshold factor:**
- Conservative: 1.5 (fewer inliers, more robust)
- Moderate: 2.0 (recommended)
- Permissive: 3.0 (more inliers, less robust)

---

## Step-by-Step Implementation

### Step 1: Add Helper Functions

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

```julia
"""
    build_temp_spline_from_coeffs(c, knots, bin_width, spline_basis, degree, _first, _last, n_knots)

Build temporary UniformSpline from coefficients for residual calculation.
Internal helper for RANSAC algorithm.
"""
function build_temp_spline_from_coeffs(
    c::Vector{T},
    knots::Vector{T},
    bin_width::T,
    spline_basis::NTuple{4, Polynomial},
    degree::Int,
    _first::T,
    _last::T,
    n_knots::Int
) where {T<:AbstractFloat}

    # Build piecewise polynomials
    function buildPieceWise(knots, bin_width, spline_basis)
        function fillDesignMatRow!(X, row, knot_idx, spline_basis)
            i = length(spline_basis)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i]
                i -= 1
            end
        end

        X = zeros(Polynomial, (length(knots), length(knots) + length(spline_basis) - 1))
        for (i, t_val) in enumerate(knots)
            t_adj = t_val + bin_width/2
            knot_idx = min(
                floor(Int32, (t_adj - first(knots))/bin_width) + one(Int32),
                length(knots)
            )
            fillDesignMatRow!(X, i, knot_idx, spline_basis)
        end
        return X
    end

    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly * c
    n_total_coeffs = n_knots * (degree + 1)
    coeffs = SVector{n_total_coeffs}(vcat([poly.coeffs for poly in piecewise_polynomials]...))

    return UniformSpline{n_total_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end

"""
    calculate_adaptive_iterations(outlier_ratio::T, sample_size::Int, confidence::T) where T

Calculate number of RANSAC iterations needed for given parameters.

# Arguments
- `outlier_ratio`: Expected fraction of outliers (0.0 to 0.5)
- `sample_size`: Number of points in each random sample
- `confidence`: Desired confidence (typically 0.99)

# Returns
- Number of iterations required
"""
function calculate_adaptive_iterations(
    outlier_ratio::T,
    sample_size::Int,
    confidence::T = T(0.99)
) where T<:AbstractFloat

    if outlier_ratio <= 0
        return 20  # Minimum
    end

    inlier_ratio = 1 - outlier_ratio
    prob_all_inliers = inlier_ratio^sample_size

    if prob_all_inliers < 1e-10
        return 200  # Maximum
    end

    k = log(1 - confidence) / log(1 - prob_all_inliers)
    return max(20, min(200, ceil(Int, k)))
end
```

### Step 2: Implement RANSAC Core Logic

**Add to `UniformSplinePenalized` after the robust_refit block:**

```julia
# RANSAC fitting mode
elseif fitting_mode == :ransac

    # Set random seed if provided (for reproducibility)
    if ransac_seed !== nothing
        Random.seed!(ransac_seed)
    end

    n_data = length(sorted_u)
    n_coeffs = n_knots + degree

    # Ensure sample size is valid
    sample_size = max(n_coeffs + 2, min(ransac_sample_size, div(n_data, 2)))

    # Step 1: Calculate adaptive threshold from initial fit
    c_initial = if λ == 0
        X \ sorted_u
    else
        (X'X + T(λ) * P) \ (X'sorted_u)
    end

    spline_initial = build_temp_spline_from_coeffs(
        c_initial, knots, bin_width, spline_basis,
        degree, _first, _last, n_knots
    )

    predicted_initial = [spline_initial(ti) for ti in sorted_t]
    residuals_initial = abs.(sorted_u .- predicted_initial)
    mad_residual = median(residuals_initial)
    inlier_threshold = mad_residual * ransac_threshold_factor

    # Step 2: Estimate outlier ratio for adaptive iterations
    outlier_ratio = sum(residuals_initial .> inlier_threshold) / n_data
    adaptive_iterations = calculate_adaptive_iterations(
        T(outlier_ratio),
        sample_size
    )
    n_iterations = min(ransac_iterations, adaptive_iterations)

    # Step 3: RANSAC main loop
    best_inlier_count = 0
    best_c = c_initial
    best_inlier_mask = trues(n_data)

    for iter in 1:n_iterations
        # Randomly sample data points
        sample_idx = sample(1:n_data, sample_size, replace=false)
        X_sample = X[sample_idx, :]
        u_sample = sorted_u[sample_idx]

        # Fit to sample
        c_sample = try
            if λ == 0
                X_sample \ u_sample
            else
                (X_sample'X_sample + T(λ) * P) \ (X_sample'u_sample)
            end
        catch e
            # Singular matrix, skip this sample
            continue
        end

        # Build temporary spline
        spline_sample = build_temp_spline_from_coeffs(
            c_sample, knots, bin_width, spline_basis,
            degree, _first, _last, n_knots
        )

        # Count inliers on ALL data
        predicted_sample = [spline_sample(ti) for ti in sorted_t]
        residuals_sample = abs.(sorted_u .- predicted_sample)
        inlier_mask_sample = residuals_sample .<= inlier_threshold
        inlier_count = sum(inlier_mask_sample)

        # Track best model
        if inlier_count > best_inlier_count
            best_inlier_count = inlier_count
            best_c = c_sample
            best_inlier_mask = inlier_mask_sample
        end

        # Early termination: if >95% are inliers, very likely optimal
        if inlier_count / n_data > 0.95
            break
        end
    end

    # Step 4: Refit using all inliers from best model
    if best_inlier_count >= n_coeffs + 2  # Need enough points
        X_inliers = X[best_inlier_mask, :]
        u_inliers = sorted_u[best_inlier_mask]

        c = if λ == 0
            X_inliers \ u_inliers
        else
            (X_inliers'X_inliers + T(λ) * P) \ (X_inliers'u_inliers)
        end
    else
        # Not enough inliers, use best sample fit
        c = best_c
    end

    # Note: Could add logging here
    # @debug "RANSAC: $best_inlier_count/$n_data inliers after $n_iterations iterations"

end  # End RANSAC mode
```

### Step 3: Update Function Documentation

```julia
"""
    UniformSplinePenalized(u, t, degree, n_knots, λ=0.1, order=2, fitting_mode=:robust_refit, ...)

Fit uniform B-spline with difference penalty and configurable robust fitting.

# Fitting Modes
- `:standard` - Single fit, no outlier handling (fastest)
- `:robust_refit` - Two-pass fitting, masks top 10% residuals (default)
- `:ransac` - Random Sample Consensus (most robust, slowest)

# RANSAC Parameters (only used when fitting_mode=:ransac)
- `ransac_iterations::Int`: Maximum RANSAC iterations (default 50)
  - Automatically adapts based on outlier detection
  - More iterations → higher confidence but slower

- `ransac_sample_size::Int`: PSMs per random sample (default 30)
  - Must be ≥ n_coeffs (8 for 5 knots)
  - Larger → less robust to outliers
  - Smaller → more variance between samples

- `ransac_threshold_factor::T`: Inlier threshold multiplier (default 2.0)
  - threshold = MAD × factor
  - Larger → more inliers (less robust)
  - Smaller → fewer inliers (more robust)

- `ransac_seed::Union{Int,Nothing}`: Random seed for reproducibility
  - Set to integer for deterministic results
  - Leave as nothing for random behavior

# Examples
```julia
# Standard RANSAC with defaults
spline = UniformSplinePenalized(irt, rt, 3, 5, 0.1, 2, :ransac)

# RANSAC with custom parameters
spline = UniformSplinePenalized(
    irt, rt, 3, 5, 0.1, 2, :ransac,
    100,  # More iterations
    25,   # Smaller sample size
    1.5,  # Stricter threshold
    42    # Fixed seed
)

# Compare all modes
s1 = UniformSplinePenalized(irt, rt, 3, 5, 0.1, 2, :standard)
s2 = UniformSplinePenalized(irt, rt, 3, 5, 0.1, 2, :robust_refit)
s3 = UniformSplinePenalized(irt, rt, 3, 5, 0.1, 2, :ransac)
```

# Performance Notes
- Standard: ~1ms per fit
- Robust refit: ~2ms per fit (2× slower)
- RANSAC: ~50-100ms per fit (50-100× slower)

RANSAC is recommended only when:
- Outlier contamination >15%
- Robust refit produces non-monotonic results
- Maximum robustness is required
"""
```

### Step 4: Add Random Dependency

**File:** `src/utils/ML/uniformBasisCubicSpline.jl` (top of file)

```julia
using Random  # For RANSAC sampling
using Statistics  # For median (already imported)
```

### Step 5: Update Call Sites (Optional)

**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

**Current calls use default `:robust_refit` mode automatically.**

**To enable RANSAC, change calls to:**
```julia
rt_to_irt_map = UniformSplinePenalized(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots,
    Float32(0.1),
    2,
    :ransac,  # Enable RANSAC
    50,       # iterations
    30,       # sample_size
    Float32(2.0)  # threshold_factor
)
```

---

## Testing Strategy

### Unit Test 1: Verify RANSAC Handles Severe Outliers

```julia
@testset "RANSAC vs Robust Refit - Severe Outliers" begin
    # Generate clean monotonic data
    n = 500
    t = Float32.(sort(rand(n) * 100))
    u_clean = t .+ 0.1f0 * randn(Float32, n)

    # Add 30% severe outliers
    n_outliers = 150
    outlier_idx = sample(1:n, n_outliers, replace=false)
    u_outliers = copy(u_clean)
    u_outliers[outlier_idx] .+= 10.0f0 * randn(Float32, n_outliers)

    # Fit with all three modes
    spline_standard = UniformSplinePenalized(u_outliers, t, 3, 5, 0.1f0, 2, :standard)
    spline_robust = UniformSplinePenalized(u_outliers, t, 3, 5, 0.1f0, 2, :robust_refit)
    spline_ransac = UniformSplinePenalized(u_outliers, t, 3, 5, 0.1f0, 2, :ransac)

    # Ground truth: fit to clean data
    spline_clean = UniformSplinePenalized(u_clean, t, 3, 5, 0.1f0, 2, :standard)

    # Compare RMSE to ground truth
    test_pts = Float32.(0:5:100)
    rmse_standard = sqrt(mean((spline_standard(tp) - spline_clean(tp))^2 for tp in test_pts))
    rmse_robust = sqrt(mean((spline_robust(tp) - spline_clean(tp))^2 for tp in test_pts))
    rmse_ransac = sqrt(mean((spline_ransac(tp) - spline_clean(tp))^2 for tp in test_pts))

    # RANSAC should be best, robust should be middle, standard worst
    @test rmse_ransac < rmse_robust < rmse_standard
    @test rmse_ransac < 0.5  # Should be close to clean fit
end
```

### Unit Test 2: RANSAC Reproducibility

```julia
@testset "RANSAC Reproducibility with Seed" begin
    n = 200
    t = Float32.(0:0.5:100)
    u = t .+ 0.2f0 * randn(Float32, n)

    # Fit twice with same seed
    spline1 = UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :ransac, 50, 30, 2.0f0, 42)
    spline2 = UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :ransac, 50, 30, 2.0f0, 42)

    # Should be identical
    test_pts = Float32.(0:10:100)
    for tp in test_pts
        @test spline1(tp) ≈ spline2(tp) atol=1e-6
    end

    # Fit without seed - should be different
    spline3 = UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :ransac)

    # Very likely to be different (not guaranteed, but >99.9%)
    diffs = [abs(spline1(tp) - spline3(tp)) for tp in test_pts]
    @test maximum(diffs) > 1e-4  # Should see some difference
end
```

### Unit Test 3: RANSAC Performance

```julia
@testset "RANSAC Performance" begin
    n = 1000
    t = Float32.(sort(rand(n) * 100))
    u = t .+ 0.1f0 * randn(Float32, n)

    # Time all three modes
    t_standard = @elapsed UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :standard)
    t_robust = @elapsed UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :robust_refit)
    t_ransac = @elapsed UniformSplinePenalized(u, t, 3, 5, 0.1f0, 2, :ransac, 50)

    println("Timing comparison:")
    println("  Standard:     $(round(t_standard*1000, digits=2))ms")
    println("  Robust refit: $(round(t_robust*1000, digits=2))ms")
    println("  RANSAC:       $(round(t_ransac*1000, digits=2))ms")

    # RANSAC should be slower but not excessively
    @test t_ransac > t_robust
    @test t_ransac < 1.0  # Should complete within 1 second
end
```

### Integration Test

```julia
@testset "ParameterTuning with RANSAC" begin
    # Would test full pipeline with RANSAC enabled
    # Check that RT models are fitted successfully
    # Verify monotonicity improvement vs robust_refit
end
```

---

## Configuration and Tuning

### Recommended Parameter Sets

**Conservative (Maximum Robustness):**
```julia
fitting_mode = :ransac
ransac_iterations = 100
ransac_sample_size = 20
ransac_threshold_factor = 1.5
```
- Best for datasets with >20% outliers
- Slowest but most robust

**Balanced (Default):**
```julia
fitting_mode = :ransac
ransac_iterations = 50
ransac_sample_size = 30
ransac_threshold_factor = 2.0
```
- Good for 10-30% outliers
- Reasonable speed

**Fast (Minimal Overhead):**
```julia
fitting_mode = :ransac
ransac_iterations = 25
ransac_sample_size = 40
ransac_threshold_factor = 2.5
```
- For <15% outliers
- Faster, less robust

### Adaptive Strategy

Could implement automatic mode selection:

```julia
function fit_with_adaptive_mode(u, t, ...)
    # Try robust_refit first
    spline_robust = UniformSplinePenalized(u, t, ..., :robust_refit)

    # Check quality
    predicted = [spline_robust(ti) for ti in t]
    residuals = abs.(u .- predicted)
    outlier_ratio = sum(residuals .> 3*median(residuals)) / length(u)

    # Use RANSAC if high outlier contamination
    if outlier_ratio > 0.15
        @warn "High outlier ratio ($outlier_ratio), using RANSAC"
        return UniformSplinePenalized(u, t, ..., :ransac)
    else
        return spline_robust
    end
end
```

---

## Performance Analysis

### Time Complexity

**Standard fitting:**
- O(n × m²) where n=data points, m=coefficients
- For n=1000, m=8: ~1ms

**Robust refit:**
- 2 × O(n × m²) = O(n × m²)
- For n=1000, m=8: ~2ms

**RANSAC:**
- k × O(s × m²) + O(n × m²)
  - k = iterations
  - s = sample size
  - Final refit with inliers
- For k=50, s=30, n=1000, m=8: ~50ms

### Memory Usage

All modes have similar memory footprint:
- Design matrix X: n × m floats
- Coefficients: m floats
- Temporary arrays: O(n)

RANSAC adds:
- Inlier mask: n booleans
- Best model storage: m floats

**Total:** Still O(n) memory

### Scaling

| Data Size | Standard | Robust | RANSAC |
|-----------|----------|---------|---------|
| 100 PSMs | 0.5ms | 1ms | 20ms |
| 500 PSMs | 1ms | 2ms | 50ms |
| 1000 PSMs | 2ms | 4ms | 100ms |
| 5000 PSMs | 10ms | 20ms | 500ms |

**Conclusion:** RANSAC feasible even for large datasets

---

## Error Handling

### Potential Issues

**1. Singular matrix in sample**
```julia
c_sample = try
    X_sample \ u_sample
catch e
    # Skip this sample, try next iteration
    continue
end
```

**2. No good samples found**
```julia
if best_inlier_count < n_coeffs + 2
    @warn "RANSAC failed to find sufficient inliers, using initial fit"
    c = c_initial
end
```

**3. All samples fail**
```julia
if best_c === nothing
    @error "RANSAC failed completely, falling back to standard fit"
    c = X \ sorted_u
end
```

---

## When to Use RANSAC

### Use RANSAC when:
- ✅ Outlier contamination >15%
- ✅ Robust refit produces non-monotonic results
- ✅ Data has systematic outliers (e.g., contamination from another species)
- ✅ Maximum robustness required for critical applications

### Don't use RANSAC when:
- ❌ Data is mostly clean (<10% outliers)
- ❌ Speed is critical (real-time processing)
- ❌ Reproducibility is essential (use with seed if needed)
- ❌ Current method works well

### Decision Tree

```
Is monotonicity achieved with current robust_refit?
├─ Yes → DONE (don't add complexity)
└─ No → Are there >15% obvious outliers?
    ├─ Yes → Try RANSAC
    └─ No → Tune robust_refit parameters first
            (try outlier_percentile=0.85)
```

---

## Comparison with Alternatives

| Method | Outlier Tolerance | Speed | Deterministic | Complexity |
|--------|------------------|-------|---------------|------------|
| **Standard** | 0% | ⚡⚡⚡ | ✅ | Very Low |
| **Robust Refit** | 10-15% | ⚡⚡⚡ | ✅ | Low |
| **RANSAC** | 50% | ⚡ | ❌ | Medium |
| **Penalty (Optim)** | 20% | ⚡⚡ | ✅ | Medium |
| **Constrained QP** | 100% mono | ⚡ | ✅ | High |

**Recommendation hierarchy:**
1. Start with robust_refit (default)
2. If failures, try tuning outlier_percentile
3. If still failing, try RANSAC
4. If monotonicity critical, use Constrained QP

---

## Implementation Checklist

- [ ] Add `build_temp_spline_from_coeffs` helper
- [ ] Add `calculate_adaptive_iterations` helper
- [ ] Add `using Random` import
- [ ] Implement RANSAC core logic in `UniformSplinePenalized`
- [ ] Update function documentation
- [ ] Add unit tests for RANSAC
- [ ] Test on synthetic data with outliers
- [ ] Compare performance vs robust_refit
- [ ] Document when to use each mode
- [ ] Add to ParameterTuningSearch (optional)

---

## Timeline Estimate

| Task | Time | Cumulative |
|------|------|------------|
| Add helper functions | 15 min | 15 min |
| Implement RANSAC core | 30 min | 45 min |
| Update documentation | 15 min | 60 min |
| Write unit tests | 20 min | 80 min |
| Test and validate | 20 min | 100 min |
| Integration | 10 min | 110 min |
| **Total** | **~2 hours** | |

---

## Commit Strategy

**Commit 1:** Add helper functions and imports
```bash
git commit -m "Add RANSAC helper functions for robust spline fitting"
```

**Commit 2:** Implement RANSAC core logic
```bash
git commit -m "Implement RANSAC fitting mode for UniformSplinePenalized

Adds Random Sample Consensus algorithm as third fitting mode.
Handles up to 50% outlier contamination.

Usage: UniformSplinePenalized(..., :ransac)

Parameters:
- ransac_iterations: Max iterations (default 50)
- ransac_sample_size: PSMs per sample (default 30)
- ransac_threshold_factor: Inlier threshold (default 2.0)
- ransac_seed: For reproducibility (optional)

Performance: ~50ms per fit (vs 2ms for robust_refit)
Use only when robust_refit insufficient."
```

**Commit 3:** Add tests and documentation
```bash
git commit -m "Add tests and documentation for RANSAC mode"
```

---

## Future Enhancements

### Possible Improvements

1. **Progressive RANSAC:**
   - Start with large sample size, reduce if too slow
   - Adaptive sample size based on data density

2. **MSAC (M-estimator SAC):**
   - Instead of counting inliers, sum squared residuals
   - More continuous objective function

3. **MLESAC (Maximum Likelihood SAC):**
   - Explicit noise model
   - Better theoretical foundation

4. **Guided sampling:**
   - Weight samples by residuals from initial fit
   - Avoid sampling known outliers

5. **Parallel RANSAC:**
   - Run multiple RANSAC trials in parallel
   - Take best result
   - Faster wall-clock time on multi-core

---

## Conclusion

RANSAC is a powerful robust fitting algorithm, but adds complexity.

**Recommendation:**
1. Test current robust_refit implementation first
2. Only implement RANSAC if needed
3. Consider simpler alternatives (tuning outlier_percentile)
4. RANSAC is well-documented here if needed later

The current two-pass robust fitting is likely sufficient for RT data, which typically has mild outliers from measurement noise rather than severe contamination.
