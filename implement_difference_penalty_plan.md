# Implementation Plan: 2nd Order Difference Penalty for RT Splines

## Quick Reference

**Goal:** Add 2nd order difference penalty regularization to RT→iRT spline fitting

**Method:** P-spline approach - solve (X'X + λD'D)c = X'u instead of Xc = u

**Time estimate:** 30-60 minutes implementation + 15 minutes testing

**Files to modify:**
- `src/utils/ML/uniformBasisCubicSpline.jl` - Add new functions
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl` - Use new functions

**Dependencies:** None (pure linear algebra)

**Speed:** Same as unconstrained (~1ms per fit)

**Monotonicity:** Not guaranteed, but likely improved through smoothing

---

## Why Start Here?

✅ **Fastest to implement** - Simple linear algebra, no optimization
✅ **No new dependencies** - Uses existing Julia stdlib
✅ **Same performance** - No speed penalty
✅ **May be sufficient** - If RT data is naturally monotonic
✅ **Easy to test** - Compare before/after monotonicity rates

If this doesn't work well enough, you can upgrade to the penalty method (Optim.jl) or constrained QP (JuMP) later.

---

## Step-by-Step Implementation

### Step 1: Add Helper Function for Difference Matrix

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

**Location:** Add after the UniformSpline struct definition and evaluation function, before the main UniformSpline constructor

**Code to add:**

```julia
"""
    build_difference_matrix(n_coeffs::Int, order::Int=2)

Build k-th order difference matrix for penalized B-splines (P-splines).

The difference matrix D is used to penalize roughness in spline coefficients:
- Order 1: Penalizes differences (cᵢ₊₁ - cᵢ)
- Order 2: Penalizes curvature (cᵢ₊₂ - 2cᵢ₊₁ + cᵢ) [RECOMMENDED]
- Order 3: Penalizes jerk (cᵢ₊₃ - 3cᵢ₊₂ + 3cᵢ₊₁ - cᵢ)

# Arguments
- `n_coeffs::Int`: Number of spline coefficients
- `order::Int`: Order of differences (1, 2, or 3; default 2)

# Returns
- `D::Matrix`: Difference matrix of size (n_coeffs - order) × n_coeffs

# Example
```julia
# For 8 coefficients with 2nd order penalty
D = build_difference_matrix(8, 2)  # Returns 6×8 matrix
P = D' * D  # 8×8 penalty matrix
```

# References
Eilers & Marx (1996). "Flexible smoothing with B-splines and penalties."
Statistical Science, 11(2), 89-121.
"""
function build_difference_matrix(n_coeffs::Int, order::Int=2)
    if order < 1 || order > 3
        error("Difference order must be 1, 2, or 3")
    end

    if order == 1
        # First order difference: D[i,:] = [0, ..., -1, 1, ..., 0]
        n_rows = n_coeffs - 1
        D = zeros(Float64, n_rows, n_coeffs)
        for i in 1:n_rows
            D[i, i] = -1.0
            D[i, i+1] = 1.0
        end
        return D

    elseif order == 2
        # Second order: recursively compute D² = D₁(n-1) × D₁(n)
        D1_n = build_difference_matrix(n_coeffs, 1)  # (n-1) × n
        D1_n_minus_1 = build_difference_matrix(n_coeffs - 1, 1)  # (n-2) × (n-1)
        return D1_n_minus_1 * D1_n  # (n-2) × n

    else  # order == 3
        # Third order: D³ = D₁(n-2) × D²(n)
        D2 = build_difference_matrix(n_coeffs, 2)  # (n-2) × n
        D1_n_minus_2 = build_difference_matrix(n_coeffs - 2, 1)  # (n-3) × (n-2)
        return D1_n_minus_2 * D2  # (n-3) × n
    end
end
```

### Step 2: Add Penalized Spline Fitting Function

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

**Location:** Add after the `build_difference_matrix` function

**Code to add:**

```julia
"""
    UniformSplinePenalized(u, t, degree, n_knots, λ=0.1, order=2)

Fit uniform B-spline with difference penalty regularization (P-spline method).

This uses the P-spline approach (Eilers & Marx, 1996) which adds a penalty term
to the least squares objective:

    minimize: ||Xc - u||² + λ||D^k c||²

The penalty encourages smooth coefficient sequences, which may improve monotonicity
by reducing oscillations caused by noise or outliers.

# Arguments
- `u::Vector{T}`: Target values (iRT in RT alignment context)
- `t::Vector{T}`: Input values (RT in RT alignment context)
- `degree::I`: Polynomial degree (must be 3 for cubic splines)
- `n_knots::I`: Number of knots (typically 5)
- `λ::T`: Penalty parameter (default 0.1)
  - λ = 0: No penalty (standard least squares)
  - λ small (0.01-0.1): Mild smoothing
  - λ medium (0.1-1.0): Moderate smoothing [RECOMMENDED]
  - λ large (1.0-10.0): Strong smoothing
- `order::Int`: Order of difference penalty (default 2)
  - Order 2 is recommended (penalizes curvature changes)

# Returns
- `UniformSpline{N, T}`: Fitted spline with difference penalty

# Notes
- Same computational cost as unconstrained fitting
- Does NOT guarantee monotonicity, but may improve it
- For 5 knots: 8 coefficients, penalty matrix is 6×8 (2nd order)

# Example
```julia
# Basic usage with default penalty
spline = UniformSplinePenalized(irt, rt, 3, 5)

# Custom penalty strength
spline = UniformSplinePenalized(irt, rt, 3, 5, 0.5)  # Stronger smoothing

# Test different penalties
for λ in [0.01, 0.1, 1.0]
    s = UniformSplinePenalized(irt, rt, 3, 5, λ)
    println("λ=$λ: RMSE=$(rmse(s, irt, rt))")
end
```

# References
Eilers, P. H. C., & Marx, B. D. (1996). Statistical Science, 11(2), 89-121.
"""
function UniformSplinePenalized(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I,
    λ::T = T(0.1),
    order::Int = 2
) where {I<:Integer, T<:AbstractFloat}

    # Input validation
    if degree != 3
        error("Non-cubic splines not yet implemented. Use degree = 3")
    end
    if n_knots < 3
        error("Need at least 3 knots")
    end
    if length(u) != length(t)
        error("length(u) must equal length(t)")
    end
    if λ < 0
        error("Penalty parameter λ must be non-negative")
    end

    # Sort data for numerical stability
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    # Build B-spline basis (reuse existing internal functions)
    # These functions should already exist in this file from the original UniformSpline
    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, knots, bin_width, spline_basis)

    # Compute coefficients with penalty
    n_coeffs = n_knots + degree  # = 8 for n_knots=5, degree=3

    if λ == 0
        # No penalty: standard least squares
        c = X \ sorted_u
    else
        # Penalized least squares
        D = build_difference_matrix(n_coeffs, order)
        P = D' * D  # Penalty matrix (n_coeffs × n_coeffs)

        # Solve: (X'X + λP)c = X'u
        c = (X'X + T(λ) * P) \ (X'sorted_u)
    end

    # Build piecewise polynomials (reuse existing function)
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
```

### Step 3: Integrate into ParameterTuningSearch

**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

**Find and replace 4 locations where `UniformSpline` is called:**

#### Location 1: Line ~445 (Comparison path, initial fit)

```julia
# FIND (around line 445):
rt_to_irt_map = UniformSpline(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots
)

# REPLACE WITH:
rt_to_irt_map = UniformSplinePenalized(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots,
    Float32(0.1),  # Default penalty parameter
    2  # 2nd order penalty
)
```

#### Location 2: Line ~467 (Comparison path, refit after outliers)

```julia
# FIND (around line 467):
final_spline = UniformSpline(
    valid_psms[!, :irt_predicted],
    valid_psms[!, :rt],
    getSplineDegree(params),
    n_knots_final
)

# REPLACE WITH:
final_spline = UniformSplinePenalized(
    valid_psms[!, :irt_predicted],
    valid_psms[!, :rt],
    getSplineDegree(params),
    n_knots_final,
    Float32(0.1),
    2
)
```

#### Location 3: Line ~542 (Standard path, initial fit)

```julia
# FIND (around line 542):
rt_to_irt_map = UniformSpline(
    psms[!,:irt_predicted],
    psms[!,:rt],
    getSplineDegree(params),
    n_knots
)

# REPLACE WITH:
rt_to_irt_map = UniformSplinePenalized(
    psms[!,:irt_predicted],
    psms[!,:rt],
    getSplineDegree(params),
    n_knots,
    Float32(0.1),
    2
)
```

#### Location 4: Line ~564 (Standard path, final fit)

```julia
# FIND (around line 564):
final_model = SplineRtConversionModel(UniformSpline(
    valid_psms[!,:irt_predicted],
    valid_psms[!,:rt],
    getSplineDegree(params),
    n_knots_final
))

# REPLACE WITH:
final_model = SplineRtConversionModel(UniformSplinePenalized(
    valid_psms[!,:irt_predicted],
    valid_psms[!,:rt],
    getSplineDegree(params),
    n_knots_final,
    Float32(0.1),
    2
))
```

---

## Testing Strategy

### Unit Test 1: Verify Penalty Matrix Construction

**File:** `test/UnitTests/test_difference_penalty.jl` (create new)

```julia
using Test
using Pioneer
using Pioneer: build_difference_matrix

@testset "Difference Matrix Construction" begin
    @testset "1st order difference (n=5)" begin
        D = build_difference_matrix(5, 1)

        # Should be 4×5 matrix
        @test size(D) == (4, 5)

        # Check structure: each row is [-1, 1, 0, ...]
        @test D[1, :] ≈ [-1, 1, 0, 0, 0]
        @test D[2, :] ≈ [0, -1, 1, 0, 0]
        @test D[3, :] ≈ [0, 0, -1, 1, 0]
        @test D[4, :] ≈ [0, 0, 0, -1, 1]
    end

    @testset "2nd order difference (n=5)" begin
        D2 = build_difference_matrix(5, 2)

        # Should be 3×5 matrix
        @test size(D2) == (3, 5)

        # Check structure: each row is [1, -2, 1, 0, ...]
        @test D2[1, :] ≈ [1, -2, 1, 0, 0]
        @test D2[2, :] ≈ [0, 1, -2, 1, 0]
        @test D2[3, :] ≈ [0, 0, 1, -2, 1]
    end

    @testset "2nd order for 8 coeffs (5 knots)" begin
        D2 = build_difference_matrix(8, 2)

        # Should be 6×8 matrix
        @test size(D2) == (6, 8)

        # Verify first and last rows
        @test D2[1, 1:3] ≈ [1, -2, 1]
        @test D2[6, 6:8] ≈ [1, -2, 1]
    end
end
```

### Unit Test 2: Compare Penalized vs Unpenalized

```julia
using Pioneer
using Pioneer: UniformSpline, UniformSplinePenalized, check_monotonicity

@testset "Penalized Spline Fitting" begin
    @testset "λ=0 matches unconstrained" begin
        # Generate simple monotonic data
        t = Float32.(0:0.1:10)
        u = t .+ 0.1f0 * randn(Float32, length(t))

        # Fit both ways
        spline_unconstrained = UniformSpline(u, t, 3, 5)
        spline_penalized_zero = UniformSplinePenalized(u, t, 3, 5, 0.0f0)

        # Should give nearly identical results
        test_points = Float32.(0:0.5:10)
        max_diff = maximum(abs(spline_unconstrained(tp) - spline_penalized_zero(tp))
                          for tp in test_points)

        @test max_diff < 0.01  # Should be nearly identical
    end

    @testset "Penalty reduces oscillations" begin
        # Generate data with non-monotonic noise
        t = Float32.(0:0.1:10)
        u_base = t  # True monotonic signal
        u_noisy = u_base .+ 0.5f0 * randn(Float32, length(t))  # Add noise

        # Fit with different penalties
        spline_none = UniformSpline(u_noisy, t, 3, 5)
        spline_light = UniformSplinePenalized(u_noisy, t, 3, 5, 0.1f0)
        spline_heavy = UniformSplinePenalized(u_noisy, t, 3, 5, 1.0f0)

        # Compute RMS of second derivative (roughness measure)
        function compute_roughness(spline, t_pts)
            # Approximate second derivative
            dt = 0.01f0
            d2 = [(spline(tp+dt) - 2*spline(tp) + spline(tp-dt)) / dt^2
                  for tp in t_pts[2:end-1]]
            return sqrt(mean(abs2, d2))
        end

        test_pts = Float32.(0.1:0.1:9.9)
        roughness_none = compute_roughness(spline_none, test_pts)
        roughness_light = compute_roughness(spline_light, test_pts)
        roughness_heavy = compute_roughness(spline_heavy, test_pts)

        # Heavier penalty should reduce roughness
        @test roughness_light < roughness_none
        @test roughness_heavy < roughness_light
    end

    @testset "Performance check" begin
        # Large dataset
        n = 1000
        t = sort(rand(Float32, n) * 100)
        u = t .+ 0.1f0 * randn(Float32, n)

        # Time both methods
        t_unconstrained = @elapsed UniformSpline(u, t, 3, 5)
        t_penalized = @elapsed UniformSplinePenalized(u, t, 3, 5, 0.1f0)

        # Should be similar speed
        @test t_penalized < 2 * t_unconstrained  # At most 2x slower

        println("Timing: unconstrained=$(round(t_unconstrained*1000, digits=2))ms, " *
                "penalized=$(round(t_penalized*1000, digits=2))ms")
    end
end
```

### Integration Test: Full Pipeline

```julia
@testset "ParameterTuning with Penalized Splines" begin
    # This would run the full pipeline on test data
    # For now, just check that it runs without errors

    params_path = joinpath(@__DIR__, "../test_config/ecoli_test_params.json")

    @test_nowarn SearchDIA(params_path)

    # Could add more checks here:
    # - Load RT models from results
    # - Check monotonicity rate
    # - Compare PSM counts before/after
end
```

---

## Tuning the Penalty Parameter λ

### Quick Guide

**Start with λ = 0.1** (moderate smoothing)

| λ value | Effect | When to use |
|---------|--------|-------------|
| 0.0 | No penalty (unconstrained) | Baseline comparison |
| 0.01-0.05 | Very mild smoothing | Large datasets (>5000 PSMs) |
| 0.1-0.5 | Moderate smoothing | **RECOMMENDED for most cases** |
| 1.0-5.0 | Strong smoothing | Small datasets (<500 PSMs) |
| >10.0 | Very strong (nearly linear) | Only if severe overfitting |

### Empirical Tuning Procedure

```julia
# Test on one representative file
psms = load_representative_file()  # ~1000 PSMs

for λ in [0.0, 0.01, 0.1, 0.5, 1.0, 5.0]
    spline = UniformSplinePenalized(
        psms[!, :irt_predicted],
        psms[!, :rt],
        3, 5, Float32(λ)
    )

    # Check monotonicity
    is_mono, min_deriv = check_monotonicity(spline)

    # Check fit quality
    residuals = [spline(rt) - irt for (rt, irt) in
                 zip(psms[!, :rt], psms[!, :irt_predicted])]
    rmse = sqrt(mean(abs2, residuals))

    println("λ=$λ: monotonic=$is_mono, min_deriv=$(round(min_deriv, digits=4)), RMSE=$(round(rmse, digits=3))")
end
```

**Look for:**
- Smallest λ where monotonicity is achieved
- RMSE not significantly worse than λ=0
- Good balance of smoothness and fit quality

### Visual Tuning

```julia
using Plots

psms = load_representative_file()

p = plot(layout=(2,3), size=(1200, 800))

for (i, λ) in enumerate([0.0, 0.01, 0.1, 0.5, 1.0, 5.0])
    spline = UniformSplinePenalized(psms[!, :irt_predicted], psms[!, :rt], 3, 5, Float32(λ))

    # Plot data and fit
    scatter!(p[i], psms[!, :rt], psms[!, :irt_predicted],
             alpha=0.1, label="Data", title="λ=$λ")

    rt_range = range(minimum(psms[!, :rt]), maximum(psms[!, :rt]), length=100)
    plot!(p[i], rt_range, [spline(rt) for rt in rt_range],
          lw=2, label="Fit", color=:red)
end

savefig(p, "penalty_tuning.png")
```

Look for smoothness without losing important features.

---

## Success Criteria

### Minimum Requirements

- ✅ Code compiles without errors
- ✅ Unit tests pass
- ✅ Pipeline runs to completion on test data
- ✅ No performance regression (within 20% of original time)

### Improvement Metrics

**Compare before/after on full test dataset:**

| Metric | Baseline (unconstrained) | Target (penalized) |
|--------|-------------------------|-------------------|
| Files with non-monotonic splines | X% | <50% of X |
| Min derivative value (worst case) | Y | >0.9Y |
| RT alignment RMSE | Z | <1.1Z |
| Parameter tuning time per file | T | <1.2T |

**Success = Meaningful reduction in violations without significant performance/quality loss**

### Validation Checklist

- [ ] `build_difference_matrix()` produces correct matrices
- [ ] `UniformSplinePenalized(λ=0)` matches `UniformSpline()`
- [ ] Increasing λ reduces roughness/curvature
- [ ] Performance is acceptable (<2x slower)
- [ ] RT alignment plots look smooth and monotonic
- [ ] PSM counts in FirstPassSearch similar to baseline
- [ ] Overall pipeline completes successfully
- [ ] Monotonicity rate improved (check sample of files)

---

## Troubleshooting

### Issue 1: "UniformSplinePenalized not found"

**Symptom:** `ERROR: UndefVarError: UniformSplinePenalized not defined`

**Solution:**
1. Check that functions were added to `uniformBasisCubicSpline.jl`
2. Restart Julia session to reload module
3. Verify no syntax errors in new functions

### Issue 2: Singular matrix error

**Symptom:** `SingularException: matrix is singular`

**Cause:** Penalty matrix may make system ill-conditioned with very large λ

**Solution:**
```julia
# Add small ridge term for numerical stability
c = (X'X + T(λ) * P + T(1e-8) * I) \ (X'sorted_u)
```

### Issue 3: Performance degradation

**Symptom:** Pipeline takes much longer

**Cause:** Unlikely with this method, but check:

**Solution:**
1. Verify λ is reasonable (<10.0)
2. Check that penalty matrix is cached appropriately
3. Profile to find bottleneck

### Issue 4: No improvement in monotonicity

**Symptom:** Still getting non-monotonic splines

**Diagnosis:**
```julia
# Check if penalty is actually being applied
spline_no_penalty = UniformSpline(u, t, 3, 5)
spline_penalty = UniformSplinePenalized(u, t, 3, 5, 0.1f0)

# Compare coefficients - should be different
println("Coefficient changes: ",
        sum(abs2, spline_penalty.model.coeffs - spline_no_penalty.coeffs))
```

**Solutions:**
1. Increase λ (try 1.0 or 5.0)
2. Check data quality (plot raw RT vs iRT)
3. Consider upgrading to penalty method (Optim.jl) or constrained QP

### Issue 5: Over-smoothing

**Symptom:** RT alignment quality decreased, RMSE increased

**Solution:**
1. Decrease λ (try 0.01 or 0.05)
2. Verify that 5 knots is appropriate for your RT range
3. Check if outlier removal is too aggressive

---

## Rollback Procedure

If this approach doesn't work well:

### Quick Rollback

```bash
# Revert to safe checkpoint commit
git reset --hard ec3137dc  # Original safe checkpoint

# Or just undo the difference penalty changes
git checkout HEAD -- src/utils/ML/uniformBasisCubicSpline.jl
git checkout HEAD -- src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl
```

### Keep Changes but Disable

```julia
# In utils.jl, change penalty to 0.0 (effectively disables it)
rt_to_irt_map = UniformSplinePenalized(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots,
    Float32(0.0),  # ← Change to 0.0 to disable penalty
    2
)
```

### Next Steps After Rollback

If difference penalty insufficient:

1. **Try penalty method** (see `penalty_monotonic_spline_plan.md`)
   - Uses Optim.jl for derivative-based penalty
   - ~99% success rate
   - ~2-5x slower

2. **Try constrained QP** (see `monotonic_spline_plan.md`)
   - Uses JuMP + Ipopt for strict monotonicity
   - 100% guarantee
   - ~5-50x slower

---

## Timeline

| Task | Time | Cumulative |
|------|------|------------|
| Add `build_difference_matrix()` | 5 min | 5 min |
| Add `UniformSplinePenalized()` | 15 min | 20 min |
| Update 4 call sites in utils.jl | 10 min | 30 min |
| Write unit tests | 15 min | 45 min |
| Run tests and fix issues | 15 min | 60 min |
| Test on real data | 15 min | 75 min |
| Tune λ parameter | 15 min | 90 min |
| **Total** | **~90 min** | **with tuning** |

**Core implementation: 30 minutes**
**With testing: 60 minutes**
**With tuning: 90 minutes**

---

## Commit Message

```bash
git add src/utils/ML/uniformBasisCubicSpline.jl
git add src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl
git add test/UnitTests/test_difference_penalty.jl

git commit -m "Add 2nd order difference penalty for RT spline smoothing

Implements P-spline approach (Eilers & Marx 1996) to regularize RT→iRT
spline fitting using 2nd order difference penalty on coefficients.

Changes:
- Add build_difference_matrix() for penalty matrix construction
- Add UniformSplinePenalized() with λ parameter (default 0.1)
- Update all UniformSpline calls to use penalized version
- Add unit tests for penalty matrix and penalized fitting

Benefits:
- Same speed as unconstrained fitting (pure linear algebra)
- May reduce non-monotonic violations through smoothing
- No new dependencies required
- Easy to tune via λ parameter

Does NOT guarantee monotonicity, but may be sufficient in practice.
Can upgrade to penalty method (Optim) or constrained QP (JuMP) if needed.

Default λ=0.1 provides moderate smoothing. Set λ=0.0 to disable penalty.
"
```

---

## Next Steps After Implementation

1. **Test on full dataset:**
   - Run ParameterTuningSearch on all files
   - Collect monotonicity statistics
   - Compare RMSE before/after

2. **Evaluate success:**
   - If >95% files monotonic → Success! Keep this approach.
   - If 80-95% monotonic → Consider tuning λ or upgrading to penalty method
   - If <80% monotonic → Upgrade to penalty method (Optim.jl)

3. **Document results:**
   - Add findings to ParameterTuningSearch/CLAUDE.md
   - Note λ value used and success rate
   - Update user documentation if needed

4. **Optional enhancements:**
   - Make λ configurable via JSON parameters
   - Add adaptive λ selection based on data size
   - Implement cross-validation for λ tuning

---

## Summary

This plan implements the simplest possible approach to improving monotonicity:
- 2nd order difference penalty on spline coefficients
- Pure linear algebra (no optimization)
- Same computational cost as unconstrained
- May be sufficient for naturally monotonic RT data

**Start here.** If it works, you're done in 30-60 minutes. If not, you have two upgrade paths clearly documented in the other plan files.
