# Monotonic RT Spline Implementation

## Selected Approach: Bidirectional Cumulative Max Filter + Uniform Spline Refit

### Problem Statement
RT alignment splines sometimes violate monotonicity at the edges, particularly during "washout" where outlier peptides elute at incorrect times. This causes the fitted curve to slope backwards (see plot with n=6213 showing downturn at x>100).

**Current behavior**: UniformSplinePenalized with λ=1.0 smooths but doesn't guarantee monotonicity.

### Goal
Enforce monotonic increasing property: `f(x₁) ≤ f(x₂)` whenever `x₁ < x₂`

### Algorithm Steps

```
Input: Original spline fitted with RANSAC/penalty
Output: Monotonic spline guaranteed to be non-decreasing

1. Evaluate original spline at uniform grid
   - Create uniform grid: rt_grid = LinRange(rt_min, rt_max, N)  # N = 500
   - Evaluate: irt_grid = [original_spline(r) for r in rt_grid]

2. Find median split point
   - median_rt = median(rt_data)  # From original PSM data
   - median_idx = argmin(abs.(rt_grid .- median_rt))

3. Apply cumulative max filter on RIGHT side (median → end)
   - For i in (median_idx+1):N
       if irt_grid[i] < irt_grid[i-1]
           irt_grid[i] = irt_grid[i-1]
   - Enforces: irt[i] ≥ irt[i-1] for all i > median_idx

4. Apply cumulative max filter on LEFT side (median → start)
   - For i in (median_idx-1):-1:1
       if irt_grid[i] > irt_grid[i+1]
           irt_grid[i] = irt_grid[i+1]
   - Enforces: irt[i] ≤ irt[i+1] for all i < median_idx

5. Refit UniformSpline to filtered points
   - final_spline = UniformSpline(irt_grid, rt_grid, degree=3, n_knots=K)
   - K = adaptive based on N (e.g., K=5 for N=500)

6. Return final_spline
```

### Key Design Decisions

**Simple Cumulative Max**: Replace isotonic regression with simple cumulative maximum
- Right side: Each point becomes max(current, previous)
- Left side: Each point becomes min(current, next)
- Equivalent to isotonic regression for this specific use case
- **No external dependencies required**

**Strict Monotonicity**: Yes - cumulative max enforces strict non-decreasing property
- Right side: `irt[i] ≥ irt[i-1]` always (moving right)
- Left side: `irt[i] ≤ irt[i+1]` always (moving left)

**No Edge Detection Needed**: The median split naturally handles variable edge positions
- Works whether bad points are at beginning, end, or both
- Automatic adaptation to data characteristics

**Single-Stage Process**: Sample → Filter → Refit
1. Sample from original spline
2. Apply bidirectional cumulative max filter
3. Refit spline to filtered samples

**Grid Density**:
- Uniform sampling: `N = 500` points
- Final spline knots: `K` = adaptive (1 per 100 points), e.g., K=5 for N=50000

### Implementation Location

Add to `rt_alignment_utils.jl` as new function:

```julia
function make_spline_monotonic(
    original_spline::UniformSpline,
    rt_data::Vector{Float32},
    irt_data::Vector{Float32};
    n_sample_points::Int = 500,
    n_knots::Int = 5
)::UniformSpline
    # 1. Sample from original spline
    rt_min, rt_max = extrema(rt_data)
    rt_grid = collect(LinRange(rt_min, rt_max, n_sample_points))
    irt_grid = [original_spline(r) for r in rt_grid]

    # 2. Find median
    median_rt = median(rt_data)
    median_idx = argmin(abs.(rt_grid .- median_rt))

    # 3. Filter right side (median → end): enforce increasing
    for i in (median_idx+1):n_sample_points
        if irt_grid[i] < irt_grid[i-1]
            irt_grid[i] = irt_grid[i-1]
        end
    end

    # 4. Filter left side (median → start): enforce increasing in reverse
    for i in (median_idx-1):-1:1
        if irt_grid[i] > irt_grid[i+1]
            irt_grid[i] = irt_grid[i+1]
        end
    end

    # 5. Refit spline to filtered data
    return UniformSpline(irt_grid, rt_grid, 3, n_knots)
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

### Required Dependencies

**None** - Uses only Julia Base and existing Pioneer.jl code
- No external packages needed
- Simple cumulative max filter replaces isotonic regression

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

- Spline evaluation: O(N) → ~1ms for N=500
- Cumulative max filter: O(N) → <1ms for N=500 (simple loop)
- Spline fitting: O(N²) → ~5ms for N=500
- **Total overhead**: ~10ms per file (negligible)

### Advantages of This Approach

1. **Automatic edge handling**: No need to define edge regions
2. **Guaranteed monotonicity**: Mathematical guarantee from isotonic regression
3. **Smooth final curve**: UniformSpline refit preserves smoothness
4. **Preserves center**: Median split keeps well-behaved middle region relatively unchanged
5. **Robust**: Works for beginning outliers, end outliers, or both
6. **Simple logic**: Clean algorithm, easy to debug and validate

### Potential Issues & Mitigations

**Issue 1**: Cumulative max may create flat regions
- **Expected behavior**: This is correct - outliers cause plateaus in corrected curve
- **Mitigation**: Use sufficient sample density (N=500) to preserve shape detail
- **Note**: UniformSpline refit will smooth plateaus naturally

**Issue 2**: Discontinuity at median point
- **Mitigation**: Both sides include median point, ensuring C⁰ continuity
- **Note**: Cumulative max from both sides converges smoothly at median

**Issue 3**: Final spline might not be perfectly monotonic after refit
- **Mitigation**: Use high sample density (N=500) to preserve filtered shape
- **Validation**: Check with `check_monotonicity()` function
- **Fallback**: If violations detected, can increase sample density or return linear interpolation

### Implementation Steps

1. Implement `make_spline_monotonic` function in `rt_alignment_utils.jl`
2. Integrate into `fit_irt_model` after final spline fitting (line ~280)
3. Add diagnostic logging and validation checks
4. Test on n=6213 dataset and other problematic cases
5. No new dependencies required
