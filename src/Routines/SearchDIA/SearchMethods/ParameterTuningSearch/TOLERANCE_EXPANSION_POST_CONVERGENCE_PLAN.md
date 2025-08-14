# Plan: Tolerance Expansion Post-Convergence Optimization

## Status: PENDING APPROVAL

## Simple Concept

When we achieve convergence:
1. We used some tolerance (e.g., 20 ppm) to collect PSMs
2. The fitted model from those PSMs might show tighter tolerance (e.g., 18 ppm)
3. **Key insight**: We should test if using 20 * 1.5 = 30 ppm to collect PSMs yields more PSMs
4. If yes, fit a new model from the larger PSM set

## The Problem

Currently, when we converge, we immediately accept the fitted model. But the fitted model is based on PSMs collected with the current tolerance. If we had used a wider tolerance for collection, we might have found more valid PSMs, leading to a better model fit.

## Simple Solution

After convergence, before finalizing:
1. Take the tolerance that was USED for PSM collection (not the fitted result)
2. Expand it by 50%
3. Collect PSMs again with this expanded tolerance
4. If we get more PSMs, fit a new model from this larger set
5. Use the better result (more PSMs = better)

## Example Walkthrough

```
Iteration N:
- Current collection tolerance: 20 ppm
- Collect PSMs with 20 ppm → 1000 PSMs
- Fit model from these PSMs → fitted tolerance: 18 ppm
- Check convergence → CONVERGED!

Post-Convergence Expansion:
- Take the 20 ppm used for collection (NOT the 18 ppm fitted)
- Expand to 20 * 1.5 = 30 ppm
- Collect PSMs with 30 ppm → 1200 PSMs (more!)
- Fit NEW model from 1200 PSMs → fitted tolerance: 22 ppm
- Use this better result (more PSMs, still reasonable tolerance)
```

## Core Implementation Idea

```julia
# After convergence is detected
if converged
    # Get the tolerance that was USED for collection
    current_collection_tol = getCurrentCollectionTolerance(search_context, ms_file_idx)
    
    # Try with expanded tolerance
    expanded_tol = current_collection_tol * 1.5
    expanded_psms = collect_with_tolerance(expanded_tol, ...)
    
    if size(expanded_psms, 1) > size(current_psms, 1)
        # More PSMs found! Fit new model
        new_model = fit_model(expanded_psms)
        use_this_model(new_model)
    else
        # No improvement, keep original
        use_original_model()
    end
end
```

## Key Points

1. **Use collection tolerance, not fitted tolerance**: The tolerance used to collect PSMs is what limits our discovery
2. **Simple 50% expansion**: Start with fixed 1.5x factor
3. **PSM count as decision metric**: More PSMs = better (initially)
4. **Refit from expanded set**: The new model is fitted from the larger PSM set

## Questions for Approval

1. Should we make the expansion factor (1.5x) configurable or keep it fixed?
2. Should we require a minimum improvement (e.g., 5% more PSMs) or accept any improvement?
3. Should this be optional (controlled by a parameter) or always performed?
4. Should we limit how wide the final fitted tolerance can be?

Please approve this simple concept, and I'll create a detailed implementation plan.