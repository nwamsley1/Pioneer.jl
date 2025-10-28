# Implementation Plan: Penalty-Based Monotonic Spline Fitting

## Quick Reference

**Objective:** Add soft monotonicity constraints to RT→iRT spline fitting using penalty method

**Approach:** Minimize data residuals + penalty on negative derivatives

**Key advantages:**
- ✅ Minimal dependencies (Optim.jl - likely already in Project.toml)
- ✅ Fast (~2-5x slower than unconstrained)
- ✅ Simple implementation and debugging
- ✅ Works well when data is naturally monotonic

**Expected time:** 2-3 hours implementation + 1-2 hours testing

---

## Mathematical Foundation

### Penalty Objective Function

**Standard unconstrained spline:**
```
minimize: ||Xc - u||²
```

**Penalty-constrained spline:**
```
minimize: ||Xc - u||² + λ * Σᵢ [max(0, -derivative_i)]²
```

Where:
- **X**: Design matrix from B-spline basis (n_data × n_coeffs)
- **c**: Spline coefficients to solve for
- **u**: Target iRT values
- **λ**: Penalty parameter (e.g., 1000-10000)
- **derivative_i**: ds/du evaluated at sample point i

### Derivative Constraint

For cubic polynomial s(u) = c₀ + c₁u + c₂u² + c₃u³:

**Derivative:** ds/du = c₁ + 2c₂u + 3c₃u²

**Monotonicity requires:** ds/du ≥ 0 for all u ∈ [0,1] in each interval

**Penalty term:**
```julia
penalty = 0
for each interval:
    for u in [0, 0.25, 0.5, 0.75, 1.0]:
        deriv = c₁ + 2c₂*u + 3c₃*u²
        if deriv < 0:
            penalty += deriv²
return λ * penalty
```

### How It Works

1. **When derivative is positive** → penalty = 0 (constraint satisfied)
2. **When derivative is negative** → penalty grows quadratically with violation
3. **As λ increases** → solution approaches strictly monotonic (but may sacrifice fit quality)
4. **Adaptive strategy** → start small, increase λ until monotonic

---

## Implementation Steps

### Step 1: Check/Add Dependencies

**File:** `Project.toml`

Check if Optim.jl is already present:
```bash
grep -i "Optim" Project.toml
```

If not present, add to dependencies:
```toml
[deps]
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"  # For better convergence
```

### Step 2: Add Monotonicity Check Function

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

Add after the UniformSpline struct definition and (s::UniformSpline)(t) evaluation function:

```julia
"""
    check_monotonicity(s::UniformSpline; n_samples::Int=10, tol::Float64=1e-6)

Check if a spline satisfies monotonicity constraint.

# Arguments
- `s::UniformSpline`: Spline to check
- `n_samples`: Number of sample points per interval to check
- `tol`: Tolerance for numerical errors (derivatives > -tol are considered valid)

# Returns
- `is_monotonic::Bool`: True if all derivatives are ≥ -tol
- `min_derivative::Float64`: Minimum derivative value found
"""
function check_monotonicity(
    s::UniformSpline{N, T};
    n_samples::Int=10,
    tol::T=T(1e-6)
) where {N, T<:AbstractFloat}

    n_intervals = (N ÷ (s.degree + 1))
    min_deriv = T(Inf)

    for interval in 0:(n_intervals-1)
        # Get coefficients for this interval
        base_idx = interval * (s.degree + 1) + 1
        c0, c1, c2, c3 = s.coeffs[base_idx:(base_idx+3)]

        # Sample derivative at multiple points in [0, 1]
        for u in LinRange(T(0), T(1), n_samples)
            # Derivative: ds/du = c1 + 2*c2*u + 3*c3*u²
            deriv = c1 + 2*c2*u + 3*c3*u^2
            min_deriv = min(min_deriv, deriv)

            if deriv < -tol
                return false, min_deriv
            end
        end
    end

    return true, min_deriv
end
```

### Step 3: Implement Penalty-Based Spline Fitting

**File:** `src/utils/ML/uniformBasisCubicSpline.jl`

Add new function after the existing `UniformSpline` constructor:

```julia
using Optim, LineSearches

"""
    UniformSplinePenalty(u, t, degree, n_knots; kwargs...)

Fit uniform B-spline with penalty-based monotonicity enforcement.

# Arguments
- `u::Vector{T}`: Target values (iRT)
- `t::Vector{T}`: Input values (RT)
- `degree::I`: Polynomial degree (must be 3)
- `n_knots::I`: Number of knots (typically 5)

# Keyword Arguments
- `penalty::T`: Penalty parameter λ (default: 1000.0)
- `adaptive_penalty::Bool`: Use adaptive penalty strategy (default: true)
- `max_iterations::Int`: Maximum optimization iterations (default: 100)
- `verbose::Bool`: Print convergence info (default: false)

# Returns
- `UniformSpline{N, T}`: Fitted spline with monotonicity constraint
"""
function UniformSplinePenalty(
    u::Vector{T},
    t::Vector{T},
    degree::I,
    n_knots::I;
    penalty::T = T(1000.0),
    adaptive_penalty::Bool = true,
    max_iterations::Int = 100,
    verbose::Bool = false
) where {I<:Integer, T<:AbstractFloat}

    # Validate inputs
    if degree != 3
        error("Non-cubic splines not yet implemented. Use degree = 3")
    end
    if n_knots < 3
        error("Need at least 3 knots")
    end
    if length(u) != length(t)
        error("length(u) must equal length(t)")
    end

    # Sort data for stability
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    # Build B-spline basis and design matrix (reuse existing functions)
    function getSplineBasis(degree::I)
        return NTuple{4, Polynomial}([
            Polynomial([0, 0, 0, 1])/6,
            Polynomial([1, 3, 3, -3])/6,
            Polynomial([4, 0, -6, 3])/6,
            Polynomial([1, -3, 3, -1])/6,
        ])
    end

    function buildDesignMat(t::Vector{T}, knots::Vector{T}, bin_width::T, spline_basis)
        function fillDesignMatRow!(X, row, knot_idx, u, spline_basis)
            i = length(spline_basis)
            for col in range(knot_idx, knot_idx + length(spline_basis) - 1)
                X[row, col] = spline_basis[i](u)
                i -= 1
            end
        end

        X = zeros(T, (length(t), length(knots) + 3))
        for (i, t_val) in enumerate(t)
            knot_idx = min(
                floor(Int32, (t_val - first(knots))/bin_width) + one(Int32),
                length(knots)
            )
            fillDesignMatRow!(
                X, i, knot_idx,
                (t_val - knots[knot_idx])/bin_width,
                spline_basis
            )
        end
        return X
    end

    # Setup spline parameters
    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first) / (n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, knots, bin_width, spline_basis)

    # Initial guess: unconstrained least squares solution
    c_init = X \ sorted_u

    # Define objective function with penalty term
    function objective(c::Vector{T}, λ::T)
        # Data fitting term
        residuals = X * c - sorted_u
        data_loss = sum(abs2, residuals)

        # Monotonicity penalty term
        mono_penalty = zero(T)
        n_samples = 5  # Sample points per interval
        n_intervals = n_knots - 1
        n_coeffs_per_interval = degree + 1

        for interval in 0:(n_intervals-1)
            base_idx = interval * n_coeffs_per_interval + 1
            c1 = c[base_idx + 1]
            c2 = c[base_idx + 2]
            c3 = c[base_idx + 3]

            for u_sample in LinRange(T(0), T(1), n_samples)
                # Evaluate derivative: ds/du = c1 + 2*c2*u + 3*c3*u²
                deriv = c1 + T(2)*c2*u_sample + T(3)*c3*u_sample^2

                # Penalize negative derivatives
                if deriv < zero(T)
                    mono_penalty += abs2(deriv)
                end
            end
        end

        return data_loss + λ * mono_penalty
    end

    # Optimize with adaptive penalty strategy
    c = c_init
    penalties = adaptive_penalty ? [T(100), T(1000), T(10000)] : [penalty]

    for (i, λ) in enumerate(penalties)
        if verbose
            println("Penalty iteration $i, λ = $λ")
        end

        result = optimize(
            c_opt -> objective(c_opt, λ),
            c,
            LBFGS(linesearch=LineSearches.BackTracking()),
            Optim.Options(
                iterations = max_iterations,
                g_tol = T(1e-6),
                show_trace = verbose
            )
        )

        if !Optim.converged(result)
            @warn "Penalty optimization did not converge at λ=$λ, using current solution"
        end

        c = Optim.minimizer(result)

        # Build temporary spline to check monotonicity
        # (We need to convert c to full piecewise polynomial form)
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
        n_coeffs = n_knots * (degree + 1)
        coeffs = SVector{n_coeffs}(vcat([poly.coeffs for poly in piecewise_polynomials]...))

        temp_spline = UniformSpline{n_coeffs, T}(
            coeffs, degree, _first, _last, bin_width
        )

        is_monotonic, min_deriv = check_monotonicity(temp_spline, n_samples=10)

        if verbose
            println("  Monotonic: $is_monotonic, min derivative: $min_deriv")
        end

        if is_monotonic
            if verbose
                println("  ✓ Monotonicity achieved with λ=$λ")
            end
            break
        end

        if i == length(penalties)
            @warn "Failed to achieve strict monotonicity after trying all penalties. " *
                  "Min derivative: $min_deriv. Using best solution found."
        end
    end

    # Build final piecewise polynomials
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
    n_coeffs = n_knots * (degree + 1)
    coeffs = SVector{n_coeffs}(vcat([poly.coeffs for poly in piecewise_polynomials]...))

    return UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
end
```

### Step 4: Integrate into ParameterTuningSearch

**File:** `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl`

Modify the `fit_irt_model` function to use the penalty-based spline:

**Location 1: Line ~445** (comparison path, initial fit):
```julia
# OLD:
rt_to_irt_map = UniformSpline(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots
)

# NEW:
rt_to_irt_map = UniformSplinePenalty(
    psms[!, :irt_predicted],
    psms[!, :rt],
    getSplineDegree(params),
    n_knots,
    adaptive_penalty = true,
    verbose = false
)
```

**Location 2: Line ~467** (comparison path, refit after outlier removal):
```julia
# OLD:
final_spline = UniformSpline(
    valid_psms[!, :irt_predicted],
    valid_psms[!, :rt],
    getSplineDegree(params),
    n_knots_final
)

# NEW:
final_spline = UniformSplinePenalty(
    valid_psms[!, :irt_predicted],
    valid_psms[!, :rt],
    getSplineDegree(params),
    n_knots_final,
    adaptive_penalty = true,
    verbose = false
)
```

**Location 3: Line ~564** (standard path, ≥1000 PSMs):
```julia
# OLD:
final_model = SplineRtConversionModel(UniformSpline(
    valid_psms[!,:irt_predicted],
    valid_psms[!,:rt],
    getSplineDegree(params),
    n_knots_final
))

# NEW:
final_model = SplineRtConversionModel(UniformSplinePenalty(
    valid_psms[!,:irt_predicted],
    valid_psms[!,:rt],
    getSplineDegree(params),
    n_knots_final,
    adaptive_penalty = true,
    verbose = false
))
```

**Location 4: Line ~542** (standard path, initial fit):
```julia
# OLD:
rt_to_irt_map = UniformSpline(
    psms[!,:irt_predicted],
    psms[!,:rt],
    getSplineDegree(params),
    n_knots
)

# NEW:
rt_to_irt_map = UniformSplinePenalty(
    psms[!,:irt_predicted],
    psms[!,:rt],
    getSplineDegree(params),
    n_knots,
    adaptive_penalty = true,
    verbose = false
)
```

---

## Testing Strategy

### Unit Tests

**File:** `test/UnitTests/test_monotonic_spline.jl` (create new)

```julia
using Test
using Pioneer
using Pioneer: UniformSplinePenalty, check_monotonicity

@testset "Monotonic Spline - Penalty Method" begin

    @testset "Naturally monotonic data" begin
        # Linear increasing data
        t = Float32[0.0, 1.0, 2.0, 3.0, 4.0]
        u = Float32[0.0, 1.0, 2.0, 3.0, 4.0]

        spline = UniformSplinePenalty(u, t, 3, 5)
        is_mono, min_deriv = check_monotonicity(spline)

        @test is_mono
        @test min_deriv > -1e-6

        # Should match input closely
        for i in 1:5
            @test abs(spline(t[i]) - u[i]) < 0.1
        end
    end

    @testset "Non-monotonic data with dip" begin
        # Data with non-monotonic region
        t = Float32[0.0, 1.0, 2.0, 3.0, 4.0]
        u = Float32[0.0, 2.0, 1.5, 3.0, 4.0]  # Dip at t=2

        spline = UniformSplinePenalty(u, t, 3, 5)
        is_mono, min_deriv = check_monotonicity(spline)

        @test is_mono
        @test min_deriv > -1e-6

        # Should smooth out the dip
        @test spline(2.0) >= spline(1.0)
        @test spline(3.0) >= spline(2.0)
    end

    @testset "Large dataset performance" begin
        # Test with realistic dataset size
        n = 1000
        t = sort(rand(Float32, n) * 100)  # RT range 0-100
        u = t .+ 0.1f0 * randn(Float32, n)  # iRT ≈ RT + noise

        # Time unconstrained vs penalty
        t_unconstrained = @elapsed UniformSpline(u, t, 3, 5)
        t_penalty = @elapsed UniformSplinePenalty(u, t, 3, 5)

        @test t_penalty < 5 * t_unconstrained

        spline = UniformSplinePenalty(u, t, 3, 5)
        is_mono, _ = check_monotonicity(spline)
        @test is_mono
    end

    @testset "Comparison with unconstrained" begin
        # For naturally monotonic data, results should be similar
        n = 500
        t = sort(rand(Float32, n) * 100)
        u = sort(t .+ 0.05f0 * randn(Float32, n))  # Monotonic with noise

        spline_unconstrained = UniformSpline(u, t, 3, 5)
        spline_penalty = UniformSplinePenalty(u, t, 3, 5)

        # Sample at multiple points
        test_points = LinRange(minimum(t), maximum(t), 50)
        max_diff = maximum(abs(spline_unconstrained(tp) - spline_penalty(tp)) for tp in test_points)

        @test max_diff < 1.0  # Should be similar for monotonic data
    end
end
```

### Integration Test

**File:** `test/IntegrationTests/test_parameter_tuning_monotonic.jl` (create new)

```julia
using Test
using Pioneer
using Pioneer: check_monotonicity, SplineRtConversionModel

@testset "ParameterTuning with Monotonic RT" begin
    # Load test parameters
    params_path = joinpath(@__DIR__, "../test_config/ecoli_test_params.json")

    # Run parameter tuning
    @test_nowarn SearchDIA(params_path)

    # Check that RT models are monotonic
    # (Would need to extract RT models from results - implementation depends on output structure)
end
```

### Manual Validation

Run on test dataset and visually inspect RT plots:

```bash
julia> using Pioneer
julia> SearchDIA("test/test_config/ecoli_test_params.json")
```

Check RT alignment plots in output folder. All splines should be visually monotonic increasing.

---

## Configuration Parameters

### Add to JSON Config (Optional Enhancement)

```json
{
  "parameter_tuning": {
    "rt_model": {
      "n_knots": 5,
      "enforce_monotonic": true,
      "penalty_lambda": 1000.0,
      "adaptive_penalty": true,
      "verbose_optimization": false
    }
  }
}
```

**For initial implementation**, hardcode these values in the `UniformSplinePenalty` calls (as shown in Step 4).

**Future enhancement**: Extract these from params and pass as keyword arguments.

---

## Troubleshooting Guide

### Issue 1: Penalty method not achieving monotonicity

**Symptoms:**
- Warning: "Failed to achieve strict monotonicity"
- `check_monotonicity()` returns false
- Small negative derivatives observed

**Solutions:**
1. **Increase maximum penalty:**
   ```julia
   penalties = [T(100), T(1000), T(10000), T(100000)]  # Add higher value
   ```

2. **Increase sample points for penalty evaluation:**
   ```julia
   n_samples = 10  # Instead of 5 in objective function
   ```

3. **Check data quality:**
   - Plot raw RT vs iRT data
   - Look for extreme outliers
   - May need better outlier removal first

4. **Last resort - switch to constrained QP:**
   - Implement Approach 1 from main plan
   - Guarantees strict monotonicity

### Issue 2: Optimization not converging

**Symptoms:**
- Warning: "Penalty optimization did not converge"
- Results look reasonable but warning is concerning

**Solutions:**
1. **Increase iteration limit:**
   ```julia
   max_iterations = 200  # Instead of 100
   ```

2. **Relax convergence tolerance:**
   ```julia
   Optim.Options(
       iterations = max_iterations,
       g_tol = T(1e-4),  # Less strict than 1e-6
       ...
   )
   ```

3. **Try different solver:**
   ```julia
   # Instead of LBFGS:
   result = optimize(
       c_opt -> objective(c_opt, λ),
       c,
       ConjugateGradient()  # or Newton(), NelderMead()
   )
   ```

### Issue 3: Performance too slow

**Symptoms:**
- Parameter tuning takes significantly longer
- Each file takes >5 seconds

**Solutions:**
1. **Disable adaptive penalty (use fixed penalty):**
   ```julia
   UniformSplinePenalty(..., adaptive_penalty = false, penalty = T(1000))
   ```

2. **Reduce iteration limit:**
   ```julia
   max_iterations = 50  # Instead of 100
   ```

3. **Reduce sample points:**
   ```julia
   n_samples = 3  # Instead of 5 in objective
   ```

4. **Profile to find bottleneck:**
   ```julia
   using Profile
   @profile UniformSplinePenalty(u, t, 3, 5)
   Profile.print()
   ```

### Issue 4: Results worse than unconstrained

**Symptoms:**
- Higher residuals
- Worse RT alignment quality
- More PSMs falling outside predicted windows

**Solutions:**
1. **Verify data is actually non-monotonic:**
   ```julia
   # Check unconstrained spline
   spline_old = UniformSpline(u, t, 3, 5)
   is_mono, min_deriv = check_monotonicity(spline_old)
   println("Old spline monotonic: $is_mono, min deriv: $min_deriv")
   ```

   If already monotonic → constraint not needed, revert to unconstrained

2. **Reduce penalty to allow slight violations:**
   ```julia
   penalties = [T(10), T(100), T(1000)]  # Lower penalties
   ```

3. **Check if 5 knots is sufficient:**
   - May need more flexibility
   - Try n_knots = 7 temporarily

### Issue 5: Optim.jl not found

**Symptoms:**
- ERROR: ArgumentError: Package Optim not found

**Solutions:**
```julia
using Pkg
Pkg.add("Optim")
Pkg.add("LineSearches")
```

---

## Expected Outcomes

### Success Metrics

After implementation, expect:

1. **Monotonicity rate:** ≥99% of splines should satisfy strict monotonicity
2. **Fit quality:** RMSE within 5% of unconstrained splines
3. **Performance:** Parameter tuning <20% slower overall
4. **Robustness:** No increase in IdentityModel fallbacks

### Validation Checklist

- [ ] All tests pass
- [ ] No new warnings or errors in pipeline
- [ ] RT plots show monotonic curves
- [ ] PSM counts in FirstPassSearch similar to before
- [ ] Overall pipeline time increase <20%
- [ ] check_monotonicity() returns true for >99% of fits

---

## Timeline Estimate

| Task | Time | Notes |
|------|------|-------|
| Check dependencies | 5 min | Optim.jl likely already present |
| Add check_monotonicity() | 15 min | Simple function |
| Implement UniformSplinePenalty() | 60 min | Most complex part |
| Integrate into fit_irt_model | 30 min | 4 replacement sites |
| Write unit tests | 30 min | Multiple test cases |
| Run integration tests | 20 min | Full pipeline test |
| Debug/tune if needed | 30 min | Buffer for issues |
| **Total** | **3 hours** | **Plus 1 hour for unexpected issues** |

---

## Next Steps After Implementation

1. **Commit changes:**
   ```bash
   git add src/utils/ML/uniformBasisCubicSpline.jl
   git add src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl
   git add test/UnitTests/test_monotonic_spline.jl
   git commit -m "Implement penalty-based monotonic spline fitting"
   ```

2. **Run full test suite:**
   ```bash
   julia --project -e 'using Pkg; Pkg.test()'
   ```

3. **Validate on real data:**
   - Run on multiple MS files
   - Check RT plots
   - Compare PSM counts before/after

4. **Document in CLAUDE.md:**
   - Update ParameterTuningSearch/CLAUDE.md
   - Note monotonicity guarantee
   - Add troubleshooting section

5. **Optional enhancements:**
   - Add configuration parameters
   - Make verbose option configurable
   - Add diagnostics logging

---

## Fallback Strategy

If penalty method fails to work well:

1. **Revert to checkpoint:**
   ```bash
   git reset --hard ec3137dc  # Safe checkpoint commit
   ```

2. **Consider Approach 1 (Constrained QP):**
   - See `monotonic_spline_plan.md` for full plan
   - Requires JuMP + Ipopt dependencies
   - Guarantees strict monotonicity

3. **Or adjust strategy:**
   - Maybe monotonicity not actually needed?
   - Investigate if non-monotonic splines causing problems
   - Could use softer constraints (higher tolerance)
