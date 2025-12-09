# FirstPassSearch RT Alignment: Current vs Develop Branch

## Overview

The user wants to understand how RT to iRT alignment is done in FirstPassSearch and what has changed from the develop branch.

## Current State (feature/monotonic-rt-splines branch)

### FirstPassSearch RT Alignment (utils.jl:364-432)

Has **two modes** controlled by `use_robust_fitting` parameter:

#### Mode 1: Robust Fitting (`use_robust_fitting = true`)
```julia
rt_model, valid_rt, valid_irt, irt_mad = Pioneer.fit_irt_model(
    best_psms_df;
    lambda_penalty = Float32(0.1),
    ransac_threshold = 1000,
    min_psms = 10,
    spline_degree = 3,
    max_knots = 10,              # ← Changed from 7 to 10
    outlier_threshold = Float32(5.0)
)
```

Calls `fit_irt_model()` which:
1. Fits `UniformSplinePenalized` with RANSAC (if < 1000 PSMs) or without (if ≥ 1000 PSMs)
2. Removes outliers based on MAD
3. Refits the spline
4. **NEW: Applies `make_spline_monotonic()` enforcement** ← **KEY ADDITION**

#### Mode 2: Simple Fitting (`use_robust_fitting = false`)
```julia
rt_to_irt_spline = UniformSpline(best_irts, best_rts, 3, 5)
irt_to_rt_spline = UniformSpline(best_rts, best_irts, 3, 5)
```

Uses plain `UniformSpline` with 5 knots - no RANSAC, no penalty, no monotonic enforcement.

---

## Develop Branch State

### FirstPassSearch RT Alignment

**Identical structure** with same two modes, but key differences:

#### Mode 1: Robust Fitting
```julia
rt_model, valid_rt, valid_irt, irt_mad = Pioneer.fit_irt_model(
    best_psms_df;
    lambda_penalty = Float32(0.1),
    ransac_threshold = 1000,
    min_psms = 10,
    spline_degree = 3,
    max_knots = 7,               # ← Was 7
    outlier_threshold = Float32(5.0)
)
```

Calls `fit_irt_model()` which:
1. Fits `UniformSplinePenalized` with RANSAC (if < 1000 PSMs) or without
2. Removes outliers based on MAD
3. Refits the spline
4. **Does NOT apply monotonic enforcement** ← **KEY DIFFERENCE**

#### Mode 2: Simple Fitting
Same as current branch - plain `UniformSpline` with 5 knots.

---

## Key Differences Summary

| Aspect | Develop Branch | Current Branch (feature/monotonic-rt-splines) |
|--------|---------------|-------------------------------------------|
| **Monotonic Enforcement** | ❌ None | ✅ `make_spline_monotonic()` applied |
| **max_knots** | 7 | 10 |
| **min_psms** | 10 | 30 (in fit_irt_model default) |
| **Adaptive knots calculation** | `max(3, n_psms/100)` | `max(min(10, n_psms/20), 3)` then `max(that, n_psms/100)` |
| **Logging** | `@debug_l2` | `@user_info` (more verbose) |

### The Monotonic Enforcement Addition

**This is the major change.** On the current branch, after fitting the RT spline, the code applies:

```julia
# In fit_irt_model (line 380-386 in current branch)
final_map_monotonic = make_spline_monotonic(
    rt_to_irt_map,
    valid_psms[!, :rt],
    valid_psms[!, :irt_predicted],
    n_knots = n_knots
)
```

`make_spline_monotonic()` algorithm:
1. Samples the fitted spline at N grid points (N = n_psms - 1)
2. Finds median RT
3. Applies cumulative max filter going right from median (enforces increasing)
4. Applies cumulative max filter going left from median (enforces increasing backward)
5. **Refits using `UniformSplinePenalized`** with:
   - `degree = 3`
   - `n_knots = n_knots` (same as original)
   - `lambda = 1.0` (high penalty for smoothness)
   - `penalty_order = 2`

**This is where the bug occurred!** The high penalty (`λ=1.0`) with many knots (up to 59 in your case) + already-smooth monotonic-filtered data → drove boundary polynomial cubic coefficients to exactly zero → Julia's `Polynomial` dropped them → dimension mismatch.

---

## Impact on FirstPassSearch

### When `use_robust_fitting = true` (default in most configs):

**Develop Branch:**
- RT spline may have non-monotonic regions (rare but possible)
- Fewer knots (max 7 vs 10)
- Slightly less smooth at boundaries

**Current Branch:**
- RT spline is guaranteed monotonic
- More knots allowed (max 10 vs 7)
- Smoother due to high penalty in monotonic refit
- **More prone to the coefficient padding bug** (now fixed)

### When `use_robust_fitting = false`:

No differences - both use simple `UniformSpline(data, data, 3, 5)`.

---

## User's Question: "How do we simplify the fitting?"

The **complication** is in the monotonic enforcement step:

```julia
make_spline_monotonic(
    original_spline,     # Already fitted spline
    rt_data, irt_data,   # Original data
    n_knots = n_knots    # Same n_knots
)
```

This creates a **double-fitting** workflow:
1. Fit `UniformSplinePenalized` with RANSAC/standard method
2. Sample it at N grid points
3. Apply monotonic filtering
4. **Refit `UniformSplinePenalized` again** with high penalty (λ=1.0)

### Potential Simplifications:

#### Option 1: Remove Monotonic Enforcement Entirely
- Go back to develop branch behavior
- Simpler, faster, less prone to bugs
- **Tradeoff**: May get non-monotonic regions (though rare with good data)

#### Option 2: Keep Monotonic but Simplify Refit
Instead of refitting with `UniformSplinePenalized(λ=1.0)`, use simpler approach:
```julia
# After filtering, just use UniformSpline (no penalty)
final_spline = UniformSpline(irt_grid, rt_grid, 3, n_knots)
```
- Avoids the high-penalty over-smoothing that triggered the bug
- Faster (no penalized least squares solve)
- **Tradeoff**: May not be as smooth at boundaries

#### Option 3: Keep Monotonic but Reduce Penalty
Keep `UniformSplinePenalized` but use lower penalty:
```julia
final_spline = UniformSplinePenalized(
    irt_grid, rt_grid, 3, n_knots,
    Float32(0.1),  # Instead of 1.0
    2
)
```
- Less aggressive smoothing = less likely to produce exact zeros
- Still benefits from penalty regularization
- **Tradeoff**: May not enforce smoothness as strongly

#### Option 4: Use Different Knot Count for Monotonic Refit
Use fewer knots in the monotonic refit to avoid boundary issues:
```julia
final_spline = UniformSplinePenalized(
    irt_grid, rt_grid, 3,
    min(5, n_knots),  # Cap at 5 knots for refit
    Float32(1.0), 2
)
```
- Fewer knots = fewer boundary polynomials = less prone to zero coefficients
- **Tradeoff**: Less flexible spline

---

## Recommendation

Given that:
1. The coefficient padding bug is now **fixed**
2. Monotonic enforcement adds value (guarantees monotonicity)
3. The bug was rare (only with specific data patterns)

**Recommendation**: Keep the monotonic enforcement but **reduce the penalty** (Option 3):

```julia
# In make_spline_monotonic, line 192-199
final_spline = UniformSplinePenalized(
    irt_grid, rt_grid, 3, n_knots,
    Float32(0.1),   # Reduce from 1.0 to 0.1
    2
)
```

This:
- ✅ Maintains monotonic enforcement
- ✅ Reduces likelihood of exact-zero coefficients
- ✅ Faster convergence (less penalty = easier optimization)
- ✅ Keeps padding fix as safety net
- ✅ More consistent with the penalty used in the initial fit (λ=0.1)

### Alternative: Make Penalty Configurable

Add parameter to `make_spline_monotonic`:
```julia
function make_spline_monotonic(
    original_spline, rt_data, irt_data;
    n_knots::Int = 5,
    lambda::Float32 = Float32(0.1)  # Make configurable
)
```

Then users can tune based on their data characteristics.
