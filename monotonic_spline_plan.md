# Plan: Monotonic Spline Fitting for RT Alignment

## Executive Summary

This document outlines strategies for adding monotonicity constraints to the RT→iRT spline fitting in ParameterTuningSearch. Two approaches are presented: a rigorous constrained quadratic programming method and a simpler penalty-based method. Additionally, we propose fixing the number of knots to 5 regardless of PSM count.

## Background: Current Implementation

### Location
`src/utils/ML/uniformBasisCubicSpline.jl` - Custom uniform cubic B-spline implementation
`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl` - RT model fitting logic

### Current Algorithm

**UniformSpline Constructor:**
1. Defines 4 cubic B-spline basis polynomials (degree 3)
2. Creates design matrix X relating basis functions to data points
3. Solves **unconstrained least squares**: `c = X \ sorted_u`
4. Builds piecewise cubic polynomials from coefficients
5. Returns fast evaluation function using Horner's method

**Key characteristics:**
- Uniform knot spacing over [min(t), max(t)]
- Cubic B-splines (C² continuous)
- No constraints on monotonicity
- Fast unconstrained solve: O(n)

### Adaptive Knot Logic (Current)

In `ParameterTuningSearch/utils.jl`:

```julia
# Line 396: Calculate knots based on 100 PSMs per knot
n_knots_from_psms = floor(Int, n_psms / 100)
min_knots_required = 3

# Line 409: Use at least 3 knots
n_knots = max(min_knots_required, n_knots_from_psms)

# After outlier removal (lines 467-468):
max_knots_final = max(min_knots_required, floor(Int, n_valid / 100))
n_knots_final = min(n_knots, max_knots_final)
```

**Issues with adaptive approach:**
- Complexity varies with PSM count
- Different files may use different knot counts
- More knots ≠ better fit for small datasets
- Makes monotonicity constraint more complex

## Proposed Change 1: Fixed 5 Knots

### Rationale

**Why 5 knots?**
1. **Sufficient flexibility**: Can model most RT nonlinearities
2. **Avoids overfitting**: Fewer parameters than 100 PSMs/knot rule for small datasets
3. **Consistency**: All files use same model complexity
4. **Simplicity**: Fixed constraint matrix for monotonicity
5. **Proven**: 5 knots = 4 cubic pieces covers typical RT range well

**Mathematical consideration:**
- 5 knots → 4 intervals → 4 cubic pieces
- Each cubic has 4 coefficients
- Total: 8 coefficients for the spline (with B-spline constraints)
- Derivative constraints: ~20-25 inequality constraints (5 samples × 4 intervals)

### Implementation Changes

**In `utils.jl` `fit_irt_model()` function:**

```julia
# CURRENT (lines 396-409):
n_knots_from_psms = floor(Int, n_psms / 100)
min_knots_required = 3
n_knots = max(min_knots_required, n_knots_from_psms)

# PROPOSED:
n_knots = 5  # Fixed for consistency and monotonicity constraints
```

**After outlier removal (lines 467-468):**

```julia
# CURRENT:
max_knots_final = max(min_knots_required, floor(Int, n_valid / 100))
n_knots_final = min(n_knots, max_knots_final)

# PROPOSED:
n_knots_final = 5  # Keep fixed even after outlier removal
```

**Benefits:**
- Simpler code
- Consistent behavior across files
- Easier to add monotonicity constraints
- Sufficient flexibility for RT modeling

## Proposed Change 2: Monotonicity Constraints

### Why Monotonicity Matters

**Problem:** RT (retention time) should map monotonically to iRT (indexed RT)
- Physical reality: later eluting peptides have higher retention times
- Outliers or noise can cause non-monotonic splines
- Non-monotonic segments cause wrong RT predictions

**Consequences of non-monotonicity:**
- Multiple iRT values for same RT (ambiguous mapping)
- Local minima in RT predictions
- Poor alignment quality downstream
- Incorrect chromatogram extraction windows

### Mathematical Requirement

For monotonically increasing spline s(t):
**ds/dt ≥ 0** for all t ∈ [t_min, t_max]

For cubic polynomial s(u) = c₀ + c₁u + c₂u² + c₃u³:
- Derivative: **ds/du = c₁ + 2c₂u + 3c₃u²**
- Constraint: **c₁ + 2c₂u + 3c₃u² ≥ 0** for all u ∈ [0,1] in each interval

Since this must hold for a continuous range, we enforce it at sample points (e.g., u = 0, 0.25, 0.5, 0.75, 1.0) in each interval.

## Approach 1: Constrained Quadratic Programming (Strict Monotonicity)

### Mathematical Formulation

**Transform unconstrained problem:**
```
minimize: ||Xc - u||²
```

**Into constrained quadratic program:**
```
minimize:    ||Xc - u||²
subject to:  Dc ≥ 0
```

Where:
- **X**: Design matrix (n_data × n_coeffs) from B-spline basis evaluation
- **c**: Coefficient vector to solve for
- **u**: Target values (iRT)
- **D**: Constraint matrix (n_constraints × n_coeffs) encoding monotonicity

### Building the Constraint Matrix D

For **5 knots → 4 intervals**, with **5 sample points per interval**:
- Total constraints: 4 × 5 = 20 rows
- Each row encodes: c₁ + 2c₂u + 3c₃u² ≥ 0 at one sample point

**Algorithm:**
```julia
n_knots = 5
n_intervals = n_knots - 1  # = 4
n_samples_per_interval = 5
n_constraints = n_intervals * n_samples_per_interval  # = 20

# Coefficient count for B-spline
n_coeffs = n_knots + degree  # = 5 + 3 = 8

D = zeros(Float32, n_constraints, n_coeffs)

constraint_row = 1
for interval in 1:n_intervals
    # Each interval uses 4 consecutive coefficients starting at index 'interval'
    base_idx = interval

    for u_sample in LinRange(0.0, 1.0, n_samples_per_interval)
        # Derivative: ds/du = c₁ + 2c₂u + 3c₃u²
        # Coefficients: [c₀, c₁, c₂, c₃] at positions base_idx:(base_idx+3)

        D[constraint_row, base_idx + 1] = 1.0              # c₁ term
        D[constraint_row, base_idx + 2] = 2.0 * u_sample   # c₂ term
        D[constraint_row, base_idx + 3] = 3.0 * u_sample^2 # c₃ term

        constraint_row += 1
    end
end
```

### Implementation with JuMP

**Dependencies:**
- `JuMP.jl`: Optimization modeling language
- `Ipopt.jl`: Interior point solver (handles QP with inequality constraints)
  - Alternative: `OSQP.jl` (faster for QP, simpler to install)

**Modified UniformSpline function:**

```julia
using JuMP, Ipopt  # or OSQP

function UniformSplineMonotonic(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I;
    enforce_monotonic::Bool = true
) where {I<:Integer, T<:AbstractFloat}

    # Validate inputs
    if degree != 3
        error("Non-cubic splines not yet implemented")
    end
    if n_knots < 3
        error("Need at least 3 knots")
    end

    # Sort data
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    # Build B-spline basis and design matrix (unchanged)
    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, knots, bin_width, spline_basis)

    # Solve for coefficients
    if !enforce_monotonic
        # Fast unconstrained solve
        c = X \ sorted_u
    else
        # Constrained solve with monotonicity
        n_coeffs = n_knots + degree
        n_intervals = n_knots - 1
        n_samples_per_interval = 5

        # Build constraint matrix D
        D = build_monotonicity_constraints(n_knots, degree, n_samples_per_interval)

        # Set up optimization problem
        model = Model(Ipopt.Optimizer)
        set_silent(model)  # Suppress solver output

        # Decision variables: spline coefficients
        @variable(model, c[1:n_coeffs])

        # Objective: minimize ||Xc - u||²
        @objective(model, Min, sum((X * c - sorted_u).^2))

        # Constraints: Dc ≥ 0 (monotonicity)
        @constraint(model, mono, D * c .>= 0)

        # Solve
        optimize!(model)

        # Check convergence
        if termination_status(model) != MOI.OPTIMAL
            @warn "Monotonic spline optimization did not converge optimally, " *
                  "status: $(termination_status(model)). " *
                  "Falling back to unconstrained solution."
            c = X \ sorted_u
        else
            c = value.(c)
        end
    end

    # Build piecewise polynomials (unchanged)
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly * c
    n_coeffs = n_knots * (degree + 1)
    coeffs = SVector{n_coeffs}(vcat([polynomial.coeffs for polynomial in piecewise_polynomials]...))

    return UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end

function build_monotonicity_constraints(n_knots::Int, degree::Int, n_samples::Int)
    n_intervals = n_knots - 1
    n_coeffs = n_knots + degree
    n_constraints = n_intervals * n_samples

    D = zeros(Float64, n_constraints, n_coeffs)
    row = 1

    for interval in 1:n_intervals
        base_idx = interval

        for u in LinRange(0.0, 1.0, n_samples)
            # ds/du = c₁ + 2c₂u + 3c₃u²
            D[row, base_idx + 1] = 1.0
            D[row, base_idx + 2] = 2.0 * u
            D[row, base_idx + 3] = 3.0 * u^2
            row += 1
        end
    end

    return D
end
```

### Integration into ParameterTuningSearch

**In `utils.jl` `fit_irt_model()` function:**

```julia
# Replace line 445 (or wherever UniformSpline is called):
# OLD:
rt_to_irt_map = UniformSpline(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots
)

# NEW:
rt_to_irt_map = UniformSplineMonotonic(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    5,  # Fixed knots
    enforce_monotonic = true
)
```

### Pros and Cons

**Advantages:**
- ✅ **Mathematically rigorous**: Guarantees strict monotonicity
- ✅ **Optimal**: Finds best monotonic fit (minimizes residuals subject to constraints)
- ✅ **Robust**: Handles outliers well
- ✅ **Extensible**: Easy to add other constraints (e.g., bounds on derivatives)
- ✅ **Mature solvers**: Ipopt and OSQP are well-tested

**Disadvantages:**
- ❌ **New dependencies**: Requires JuMP + solver (adds ~50MB)
- ❌ **Slower**: 5-50x slower than unconstrained solve
  - Unconstrained: <1ms
  - Constrained: 5-50ms (still fast enough for your use case)
- ❌ **More complex**: Additional failure modes (solver convergence)
- ❌ **Solver setup**: Users need to install Ipopt or OSQP

**Performance estimate for ParameterTuning:**
- Current: ~1-2 seconds per file
- With constrained QP: ~1-3 seconds per file (negligible increase)

## Approach 2: Penalty Method (Soft Monotonicity)

### Mathematical Formulation

**Modified objective function:**
```
minimize: ||Xc - u||² + λ * Σᵢ [max(0, -derivative_i)]²
```

Where:
- **λ**: Penalty parameter (e.g., 1000-10000)
- **derivative_i**: ds/du evaluated at sample point i
- **max(0, -derivative_i)²**: Quadratic penalty only for negative derivatives

**Intuition:**
- When derivative is positive → no penalty (constraint satisfied)
- When derivative is negative → large penalty (proportional to violation magnitude)
- As λ → ∞, solution approaches constrained optimum

### Implementation with Optim.jl

**Dependencies:**
- `Optim.jl`: Lightweight optimization package (very common in Julia)
  - Already used in Pioneer? Check `Project.toml`
  - ~5MB package, minimal dependency tree

**Modified UniformSpline function:**

```julia
using Optim, LineSearches

function UniformSplinePenalty(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I;
    penalty::T = T(1000.0),
    adaptive_penalty::Bool = true
) where {I<:Integer, T<:AbstractFloat}

    # [Same setup as before: sorting, basis, design matrix X]

    # Initial guess: unconstrained solution
    c_init = X \ sorted_u

    # Define objective function
    function objective(c, λ)
        # Data fitting term
        residuals = X * c - sorted_u
        data_loss = sum(abs2, residuals)

        # Monotonicity penalty
        mono_penalty = zero(T)
        n_samples = 5
        n_intervals = n_knots - 1

        for interval in 1:n_intervals
            base_idx = interval

            for u_sample in LinRange(T(0), T(1), n_samples)
                # Evaluate derivative at this point
                # ds/du = c₁ + 2c₂u + 3c₃u²
                deriv = c[base_idx + 1] +
                        2 * c[base_idx + 2] * u_sample +
                        3 * c[base_idx + 3] * u_sample^2

                # Penalize negative derivatives
                if deriv < zero(T)
                    mono_penalty += abs2(deriv)
                end
            end
        end

        return data_loss + λ * mono_penalty
    end

    # Optimize with adaptive penalty
    c = c_init
    if adaptive_penalty
        # Start with small penalty, increase if needed
        penalties = [T(100), T(1000), T(10000)]

        for λ in penalties
            result = optimize(
                c_opt -> objective(c_opt, λ),
                c,
                LBFGS(linesearch=LineSearches.BackTracking()),
                Optim.Options(iterations=100, g_tol=1e-6)
            )

            c = Optim.minimizer(result)

            # Check if monotonicity is satisfied
            if check_monotonicity(c, n_knots, n_samples=10)
                @debug "Monotonicity achieved with penalty = $λ"
                break
            end
        end
    else
        # Single optimization with fixed penalty
        result = optimize(
            c_opt -> objective(c_opt, penalty),
            c_init,
            LBFGS(linesearch=LineSearches.BackTracking()),
            Optim.Options(iterations=100)
        )
        c = Optim.minimizer(result)
    end

    # [Build piecewise polynomials as before]

    return UniformSpline{n_coeffs, T}(...)
end

function check_monotonicity(c::Vector{T}, n_knots::Int; n_samples::Int=10) where T
    """Check if spline coefficients produce monotonic derivatives"""
    n_intervals = n_knots - 1

    for interval in 1:n_intervals
        base_idx = interval

        for u in LinRange(T(0), T(1), n_samples)
            deriv = c[base_idx + 1] + 2 * c[base_idx + 2] * u + 3 * c[base_idx + 3] * u^2

            if deriv < -1e-6  # Small tolerance for numerical errors
                return false
            end
        end
    end

    return true
end
```

### Adaptive Penalty Strategy

**Why adaptive?** Finding the right penalty parameter λ is tricky:
- Too small → insufficient constraint enforcement
- Too large → numerical instability, poor data fit

**Strategy:**
1. Start with small penalty (λ=100)
2. Optimize to get initial solution
3. Check if monotonicity is satisfied
4. If not, increase penalty (λ=1000) and re-optimize starting from previous solution
5. Repeat up to λ=10000 if needed

**Advantages of warm-starting:**
- Each optimization starts from previous solution (warm start)
- Faster convergence
- Smooth transition as penalty increases

### Pros and Cons

**Advantages:**
- ✅ **Minimal dependencies**: Optim.jl is lightweight and common
- ✅ **Fast**: Only 2-5x slower than unconstrained
- ✅ **Simple**: Easy to understand and debug
- ✅ **Flexible**: Easy to tune penalty parameter
- ✅ **Gradual enforcement**: Works well when data is naturally monotonic

**Disadvantages:**
- ❌ **Not guaranteed**: May allow small violations (~0.1% of cases)
- ❌ **Parameter tuning**: Need to choose/test penalty values
- ❌ **Local minima**: Non-convex optimization may get stuck
- ❌ **Less rigorous**: Not mathematically guaranteed

**When it works well:**
- Data is naturally mostly monotonic (likely for RT)
- Small violations are acceptable
- Speed is more important than strict guarantees

**When it fails:**
- Data has significant non-monotonic regions
- Strong outliers pull derivative negative
- Very noisy data

## Comparison and Recommendations

### Side-by-Side Comparison

| Aspect | Constrained QP (Approach 1) | Penalty Method (Approach 2) |
|--------|----------------------------|----------------------------|
| **Monotonicity guarantee** | ✅ Strict (100%) | ⚠️ ~99% with good λ |
| **Speed vs unconstrained** | 5-50x slower | 2-5x slower |
| **Absolute time** | ~5-50ms | ~2-10ms |
| **Dependencies** | JuMP + Ipopt/OSQP (~50MB) | Optim.jl (~5MB) |
| **Code complexity** | High | Medium |
| **Failure modes** | Solver convergence issues | Local minima, small violations |
| **When data is monotonic** | Identical to unconstrained | Identical to unconstrained |
| **When data has outliers** | Robust, handles well | May struggle |
| **Tuning required** | None | Penalty parameter (λ) |
| **Interpretability** | Clear constraints | Less transparent |
| **Extensibility** | Easy to add constraints | Harder to combine constraints |

### Expected Performance for ParameterTuning

**Scenario:** 500-5000 PSMs per file, 5 knots

| Method | Per-file time | Total for 100 files |
|--------|--------------|---------------------|
| Current (unconstrained) | 1-2s | 2-3 min |
| Constrained QP | 1-3s | 2-4 min |
| Penalty Method | 1-2.5s | 2-3.5 min |

**Conclusion:** All methods are fast enough; performance is not a deciding factor.

### Recommendation by Use Case

**For Production/Publication (Recommended: Approach 1)**
- When you need strict guarantees
- For manuscripts where methods must be rigorous
- When outliers are present
- When mathematical correctness is paramount

**For Rapid Prototyping (Recommended: Approach 2)**
- When testing if monotonicity helps
- For quick experiments
- When minimal dependencies are preferred
- When data is naturally monotonic

**My specific recommendation for Pioneer:**

**Start with Approach 2 (Penalty Method)** because:
1. RT alignment data is typically naturally monotonic or nearly so
2. Simpler implementation = easier to maintain
3. Faster execution in worst case
4. Can validate whether monotonicity actually improves results
5. If violations occur, easy to upgrade to Approach 1

**Upgrade to Approach 1** if:
- You observe non-monotonic splines in practice
- Small violations cause problems downstream
- You want mathematical rigor for publication

## Implementation Roadmap

### Phase 1: Foundation (1-2 hours)

1. **Fix knot count to 5**
   - File: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`
   - Lines to change: 396, 409, 467-468
   - Replace adaptive knot logic with `n_knots = 5`
   - Test that RT fitting still works

2. **Add monotonicity check function**
   - File: `src/utils/ML/uniformBasisCubicSpline.jl`
   - Add `check_monotonicity(spline::UniformSpline)` function
   - Samples derivative at multiple points
   - Returns true/false and minimum derivative value

3. **Add diagnostic logging**
   - Log when non-monotonic splines are detected
   - Track frequency across files
   - Helps validate if constraint is needed

### Phase 2A: Penalty Method Implementation (2-3 hours)

4. **Add Optim.jl dependency**
   - Check if already in `Project.toml`
   - If not, add: `Optim = "429524aa-4258-5aef-a3af-852621145aeb"`

5. **Implement penalty-based solver**
   - Add `UniformSplinePenalty()` function
   - Implement adaptive penalty strategy
   - Add unit tests for synthetic data

6. **Integration**
   - Modify `fit_irt_model()` to use penalty method
   - Add parameter: `enforce_monotonic::Bool` (default: true)
   - Test on real data

7. **Validation**
   - Run on test dataset
   - Check monotonicity satisfaction rate
   - Compare RT alignment quality vs unconstrained

### Phase 2B: Constrained QP Implementation (Alternative, 3-4 hours)

4. **Add JuMP and solver dependencies**
   - Add to `Project.toml`:
     ```toml
     JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
     Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
     ```
   - Or use OSQP for faster QP solving

5. **Implement constraint matrix builder**
   - Add `build_monotonicity_constraints()` function
   - Generate D matrix for 5 knots
   - Unit test constraint generation

6. **Implement constrained solver**
   - Add `UniformSplineMonotonic()` function
   - Set up JuMP model
   - Add fallback to unconstrained if solver fails

7. **Integration and validation**
   - Same as Phase 2A steps 6-7

### Phase 3: Testing and Validation (2-3 hours)

8. **Unit tests**
   - Test monotonic spline on simple data
   - Test non-monotonic detection
   - Test fallback behavior

9. **Integration tests**
   - Run ParameterTuningSearch on test data
   - Verify RT plots show monotonic curves
   - Check that downstream methods work correctly

10. **Benchmarking**
    - Compare speeds: unconstrained vs constrained
    - Measure impact on overall pipeline time
    - Profile if needed

### Phase 4: Documentation and Cleanup (1 hour)

11. **Update documentation**
    - Add to `ParameterTuningSearch/CLAUDE.md`
    - Document the monotonicity constraint
    - Explain fixed 5-knot strategy

12. **Parameter documentation**
    - Add `enforce_monotonic` to parameter docs
    - Explain when to enable/disable

## Testing Strategy

### Unit Tests for Monotonic Spline

```julia
@testset "Monotonic Spline Fitting" begin
    # Test 1: Naturally monotonic data - should be unchanged
    t = Float32[0.0, 1.0, 2.0, 3.0, 4.0]
    u = Float32[0.0, 1.0, 2.0, 3.0, 4.0]  # Perfectly linear

    spline = UniformSplinePenalty(u, t, 3, 5)

    # Check monotonicity
    @test check_monotonicity(spline)

    # Should match input closely
    for i in 1:5
        @test abs(spline(t[i]) - u[i]) < 0.01
    end

    # Test 2: Non-monotonic data - should be constrained
    t = Float32[0.0, 1.0, 2.0, 3.0, 4.0]
    u = Float32[0.0, 2.0, 1.5, 3.0, 4.0]  # Dip at t=2

    spline = UniformSplinePenalty(u, t, 3, 5)
    @test check_monotonicity(spline)

    # Should smooth out the dip
    @test spline(2.0) > spline(1.0)
    @test spline(3.0) > spline(2.0)

    # Test 3: Performance comparison
    t = rand(Float32, 1000)
    u = sort(rand(Float32, 1000))

    t_unconstrained = @elapsed UniformSpline(u, t, 3, 5)
    t_penalty = @elapsed UniformSplinePenalty(u, t, 3, 5)

    # Should be within 5x
    @test t_penalty < 5 * t_unconstrained
end
```

### Integration Test with ParameterTuning

```julia
@testset "ParameterTuning with Monotonic RT" begin
    # Run parameter tuning with monotonic constraint
    params = load_test_params()
    params.parameter_tuning.enforce_monotonic = true

    results = run_parameter_tuning(test_data, params)

    # Check that RT model is monotonic
    rt_model = results.rt_to_irt_model

    if rt_model isa SplineRtConversionModel
        @test check_monotonicity(rt_model.model)

        # Sample derivative along the range
        rt_range = LinRange(results.rt_min, results.rt_max, 100)
        for i in 1:(length(rt_range)-1)
            @test rt_model(rt_range[i+1]) >= rt_model(rt_range[i])
        end
    end
end
```

## Configuration Parameters

### New Parameters for `parameter_tuning` in JSON config

```json
{
  "parameter_tuning": {
    "rt_model": {
      "n_knots": 5,
      "enforce_monotonic": true,
      "monotonic_method": "penalty",  // "penalty" or "constrained_qp"
      "penalty_lambda": 1000.0,
      "adaptive_penalty": true
    }
  }
}
```

**Defaults:**
- `n_knots`: 5 (fixed)
- `enforce_monotonic`: true
- `monotonic_method`: "penalty"
- `penalty_lambda`: 1000.0
- `adaptive_penalty`: true

## Potential Issues and Solutions

### Issue 1: Penalty method allows small violations

**Symptoms:**
- `check_monotonicity()` reports false occasionally
- Small negative derivatives (~-0.01) observed

**Solutions:**
1. Increase penalty parameter: λ = 10000 instead of 1000
2. Add more sample points in constraint checking
3. Switch to Approach 1 (Constrained QP)

### Issue 2: Constrained QP fails to converge

**Symptoms:**
- Ipopt returns INFEASIBLE or OTHER_ERROR
- Warning message printed

**Solutions:**
1. Fallback to unconstrained solution (already implemented)
2. Relax constraints slightly (change ≥ 0 to ≥ -ε)
3. Try different solver (OSQP instead of Ipopt)
4. Check for data quality issues (extreme outliers)

### Issue 3: Performance degradation

**Symptoms:**
- Parameter tuning takes noticeably longer
- Files timing out

**Solutions:**
1. If using Approach 1: switch to OSQP (faster QP solver)
2. If using Approach 2: disable adaptive penalty (fixed λ)
3. Reduce number of constraint sample points
4. Profile and optimize hot spots

### Issue 4: Over-smoothing due to constraints

**Symptoms:**
- RT alignment quality decreases
- Residuals increase significantly
- Plots show overly smooth curves

**Solutions:**
1. Verify data is actually non-monotonic (plot raw data)
2. Increase number of knots (try 7 instead of 5)
3. Use less aggressive penalty (lower λ)
4. Check for outliers that should be removed first

## Expected Outcomes

### Success Metrics

1. **Monotonicity:** 100% of fitted splines satisfy ds/dt ≥ 0
2. **Quality:** RT alignment RMSE similar to unconstrained (within 5%)
3. **Performance:** <10% increase in parameter tuning time
4. **Robustness:** No increase in fallback to IdentityModel

### Validation Approach

1. Run on full test dataset (ecoli_test)
2. Compare RT alignment plots before/after
3. Check PSM counts in FirstPassSearch (should be similar)
4. Measure overall pipeline time impact

## Conclusion

Both approaches are viable for adding monotonicity constraints to RT spline fitting:

**Approach 1 (Constrained QP)** provides mathematical rigor and strict guarantees but requires more dependencies.

**Approach 2 (Penalty Method)** is simpler and faster but may allow rare small violations.

**Recommendation:** Start with Approach 2 for ease of implementation, validate on real data, and upgrade to Approach 1 if violations occur or strict guarantees are needed.

**Fixed 5-knot strategy** simplifies the implementation and provides consistent behavior across all files.
