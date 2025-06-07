# SpectralLinearRegression.jl Optimization Analysis

## Overview
This document analyzes the adaptive convergence criteria implementation in `spectralLinearRegression.jl`, specifically focusing on the relative tolerance calculation for the Huber solver.

## Current Implementation Analysis

### Relative Tolerance Formula (Line 540)
```julia
rel_tol = T(10^(-2 - (log10(min(X₁[col], max_x)) - log10(min_x))))
```

Let's break this down mathematically:

1. **Dynamic Range**: `min_x = max_x / 1e4`, so `log10(min_x) = log10(max_x) - 4`

2. **Current coefficient**: `x = X₁[col]` (clamped to max_x)

3. **Formula expansion**:
   - `log10(min(X₁[col], max_x)) - log10(min_x)`
   - = `log10(min(X₁[col], max_x)) - (log10(max_x) - 4)`
   - = `log10(min(X₁[col], max_x)) - log10(max_x) + 4`

4. **Full formula**:
   - `rel_tol = 10^(-2 - (log10(min(X₁[col], max_x)) - log10(max_x) + 4))`
   - = `10^(-2 - 4 - (log10(min(X₁[col], max_x)) - log10(max_x)))`
   - = `10^(-6 - (log10(min(X₁[col], max_x)) - log10(max_x)))`

### Analysis of Current Behavior

When `X₁[col] = max_x`:
- `log10(min(X₁[col], max_x)) - log10(max_x) = 0`
- `rel_tol = 10^(-6) = 1e-6` ✓ (close to desired 1e-7)

When `X₁[col] = min_x = max_x/1e4`:
- `log10(min_x) - log10(max_x) = -4`
- `rel_tol = 10^(-6 - (-4)) = 10^(-2) = 1e-2` ✗ (should be 1e-3)

## Issues with Current Implementation

1. **Range Mismatch**: The formula gives 1e-6 to 1e-2 instead of the desired 1e-7 to 1e-3
2. **Complexity**: The formula is harder to understand and debug
3. **Clamping**: Using `min(X₁[col], max_x)` prevents handling coefficients larger than max_x

## Suggested Alternatives

### Option 1: Direct Linear Interpolation in Log Space
```julia
# Map coefficient magnitude to tolerance
coeff_magnitude = abs(X₁[col])
min_coeff = max_x / T(1e4)

if coeff_magnitude <= min_coeff
    newton_rel_tol = T(1e-3)  # Least precise for small coefficients
elseif coeff_magnitude >= max_x
    newton_rel_tol = T(1e-7)  # Most precise for large coefficients
else
    # Linear interpolation in log space
    t = (log10(coeff_magnitude) - log10(min_coeff)) / 4  # Since log10(max_x/min_coeff) = 4
    t = clamp(t, T(0), T(1))
    # Map t ∈ [0,1] to tolerance ∈ [1e-3, 1e-7]
    log_tol = -3 - 4*t  # -3 when t=0, -7 when t=1
    newton_rel_tol = T(10)^log_tol
end
```

### Option 2: Simplified Formula Fix
```julia
# Fix the current formula to achieve desired range
# Want: max_x → 1e-7, min_x → 1e-3
# log10(tolerance) = -3 - 4 * normalized_position
x_ratio = min(X₁[col], max_x) / min_x  # Range: [1, 1e4]
normalized_pos = (log10(x_ratio)) / 4   # Range: [0, 1]
rel_tol = T(10^(-3 - 4 * normalized_pos))
```

### Option 3: Power Law Scaling
```julia
# Use power law for smooth scaling
x_normalized = clamp(abs(X₁[col]) / max_x, min_x/max_x, T(1))
# Map [1e-4, 1] to [1e-3, 1e-7] using power law
alpha = log(1e-3/1e-7) / log(1e-4)  # ≈ 1
rel_tol = T(1e-7) * (x_normalized)^alpha
```

## Recommendation

I recommend **Option 1** (Direct Linear Interpolation) because:
1. **Clarity**: The logic is explicit and easy to understand
2. **Correctness**: Guaranteed to produce exact bounds (1e-7 to 1e-3)
3. **Robustness**: Handles edge cases clearly
4. **Performance**: Simple calculations with minimal overhead
5. **Debugging**: Easy to trace and verify behavior

The current formula can be fixed with Option 2, but Option 1 provides better maintainability and clarity for future developers.

## Implementation Notes

- The dynamic range of 1e4 is maintained across all approaches
- All options respect the bounds set by `min_rel_tol` and `max_rel_tol`
- The first iteration (i=0) continues to use uniform 10% tolerance
- Consider logging the tolerance values during debugging to verify behavior