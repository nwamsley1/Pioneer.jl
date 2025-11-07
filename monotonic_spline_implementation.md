# Monotonic RT Spline Implementation

## Selected Approach: Bidirectional Isotonic Regression + Uniform Spline Refit

### Problem Statement
RT alignment splines sometimes violate monotonicity at the edges, particularly during "washout" where outlier peptides elute at incorrect times. This causes the fitted curve to slope backwards (see plot with n=6213 showing downturn at x>100).

**Current behavior**: UniformSplinePenalized with λ=1.0 smooths but doesn't guarantee monotonicity.

### Goal
Enforce monotonic increasing property: `f(x₁) ≤ f(x₂)` whenever `x₁ < x₂`

### Algorithm Steps

```
Input: Original spline fitted with RANSAC/penalty
Output: Monotonic spline guaranteed to be non-decreasing

1. Evaluate original spline at dense grid
   - Create uniform grid: rt_grid = LinRange(rt_min, rt_max, N)
   - Evaluate: irt_grid = [original_spline(r) for r in rt_grid]

2. Find median split point
   - median_rt = median(rt_data)  # From original PSM data
   - median_idx = argmin(abs.(rt_grid .- median_rt))

3. Apply isotonic regression on RIGHT side (median → end)
   - rt_right = rt_grid[median_idx:end]
   - irt_right = irt_grid[median_idx:end]
   - irt_right_isotonic = isotonic_regression(rt_right, irt_right, increasing=true)
   - Enforces: irt[i] ≤ irt[i+1] for all i (left to right)

4. Apply isotonic regression on LEFT side (start → median)
   - rt_left = rt_grid[1:median_idx]
   - irt_left = irt_grid[1:median_idx]
   - Reverse both: rt_left_rev = reverse(rt_left), irt_left_rev = reverse(irt_left)
   - irt_left_rev_isotonic = isotonic_regression(rt_left_rev, irt_left_rev, increasing=true)
   - Reverse back: irt_left_isotonic = reverse(irt_left_rev_isotonic)
   - Enforces: irt[i] ≥ irt[i+1] for all i (right to left)

5. Concatenate corrected halves
   - rt_monotonic = [rt_left; rt_right[2:end]]  # Avoid duplicate median point
   - irt_monotonic = [irt_left_isotonic; irt_right_isotonic[2:end]]

6. Sample uniformly from corrected curve
   - Create uniform sampling: rt_sample = LinRange(rt_min, rt_max, M)
   - Interpolate: irt_sample = linear_interpolate(rt_monotonic, irt_monotonic, rt_sample)

7. Refit UniformSpline to sampled points
   - final_spline = UniformSpline(irt_sample, rt_sample, degree=3, n_knots=K)
   - This gives smooth monotonic curve

8. Return final_spline wrapped in SplineRtConversionModel
```

### Key Design Decisions

**Strict Monotonicity**: Yes - isotonic regression enforces strict non-decreasing property
- Right side: `irt[i] ≤ irt[i+1]` always
- Left side: `irt[i] ≥ irt[i+1]` always (decreasing RT → decreasing iRT)

**No Edge Detection Needed**: The median split naturally handles variable edge positions
- Works whether bad points are at beginning, end, or both
- Automatic adaptation to data characteristics

**Two-Stage Process**:
1. Isotonic correction → guarantees monotonicity
2. Spline refit → recovers smoothness

**Grid Densities**:
- Initial evaluation: `N = 1000` points (captures spline shape)
- Final sampling: `M = 500` points (smooth refit)
- Final spline knots: `K` = adaptive (1 per 100/M points), e.g., K=5 for M=500

### Implementation Location

Add to `rt_alignment_utils.jl` as new function:

```julia
function make_spline_monotonic(
    original_spline::UniformSpline,
    rt_data::Vector{Float32},
    irt_data::Vector{Float32};
    n_eval_points::Int = 1000,
    n_sample_points::Int = 500,
    n_knots::Int = 5
)::UniformSpline
    # Implementation here
end
```

**Integration point**: Call after spline fitting in `fit_irt_model`:
```julia
# Line ~280 in rt_alignment_utils.jl, after fitting final_map
final_map_monotonic = make_spline_monotonic(
    final_map,
    valid_psms[!, :rt],
    valid_psms[!, :irt_predicted]
)
final_model = SplineRtConversionModel(final_map_monotonic)
```

### Required Julia Package

**IsotonicRegression.jl**: Provides efficient isotonic regression
- Installation: `using Pkg; Pkg.add("IsotonicRegression")`
- Usage: `using IsotonicRegression: isotonic_regression`
- Algorithm: Pool Adjacent Violators (PAV) - O(n) complexity

### Validation Strategy

Add diagnostic checks:
```julia
# Check monotonicity of final spline
function check_monotonicity(spline, rt_range)
    rt_test = LinRange(minimum(rt_range), maximum(rt_range), 10000)
    irt_test = [spline(r) for r in rt_test]
    violations = sum(diff(irt_test) .< 0)
    if violations > 0
        @warn "Monotonicity violations detected: $violations"
    end
    return violations == 0
end
```

### Performance Considerations

- Isotonic regression: O(N) with PAV algorithm → ~1ms for N=1000
- Spline evaluation: O(N) → ~1ms for N=1000
- Spline fitting: O(M²) → ~10ms for M=500
- **Total overhead**: ~15-20ms per file (acceptable)

### Advantages of This Approach

1. **Automatic edge handling**: No need to define edge regions
2. **Guaranteed monotonicity**: Mathematical guarantee from isotonic regression
3. **Smooth final curve**: UniformSpline refit preserves smoothness
4. **Preserves center**: Median split keeps well-behaved middle region relatively unchanged
5. **Robust**: Works for beginning outliers, end outliers, or both
6. **Simple logic**: Clean algorithm, easy to debug and validate

### Potential Issues & Mitigations

**Issue 1**: Isotonic regression may over-correct (create flat regions)
- **Mitigation**: Use fine grid (N=1000) to preserve detail
- **Alternative**: Add small regularization to isotonic regression

**Issue 2**: Discontinuity at median point
- **Mitigation**: Both sides include median point, ensuring C⁰ continuity
- **Note**: UniformSpline refit will smooth any kinks

**Issue 3**: Final spline might not be perfectly monotonic after refit
- **Mitigation**: Use high sample density (M=500) to preserve isotonic shape
- **Validation**: Check with `check_monotonicity()` function
- **Fallback**: If violations detected, return isotonic-corrected points with linear interpolation instead

### Implementation Steps

1. Implement `make_spline_monotonic` function
2. Add `IsotonicRegression.jl` dependency to Project.toml
3. Integrate into `fit_irt_model` after final spline fitting
4. Add diagnostic logging and validation checks
5. Test on n=6213 dataset and other problematic cases
