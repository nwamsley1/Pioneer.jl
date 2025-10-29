# Use Unified RT Model Fitting in FirstPassSearch

## Executive Summary

Replace FirstPassSearch's simple spline fitting with the robust `fit_irt_model()` function from ParameterTuningSearch. Add optional RT alignment plotting for FirstPassSearch (disabled by default).

## Current State

### FirstPassSearch RT Fitting (utils.jl:366-377)
```julia
irt_to_rt_spline = UniformSpline(best_rts, best_irts, 3, 5)
rt_to_irt_spline = UniformSpline(best_irts, best_rts, 3, 5)
```

**Issues:**
- Uses basic UniformSpline (no penalty, no RANSAC)
- Fixed parameters: degree=3, n_knots=5
- No outlier removal
- No adaptive knot selection based on PSM count
- Falls back to IdentityModel on any error

### ParameterTuningSearch RT Fitting (utils.jl:390-512)
```julia
fit_irt_model(params, psms) returns (model, rt, irt, mad)
```

**Benefits:**
- Conditional RANSAC for <1000 PSMs
- Smoothing penalty (configurable λ)
- Adaptive knots: 1 per 100 PSMs, minimum 3
- Outlier removal based on MAD threshold
- Linear model fallback only on error
- Consistent with ParameterTuning phase

---

## Proposed Changes

### 1. Add Configuration Parameters

**File**: `assets/example_config/defaultSearchParams.json`

Add to `first_search.irt_mapping` section:
```json
{
  "first_search": {
    "irt_mapping": {
      "max_prob_to_impute_irt": 0.75,
      "fwhm_nstd": 4,
      "irt_nstd": 4,
      "plot_rt_alignment": false,
      "use_robust_fitting": true
    }
  }
}
```

**New Parameters:**
- `plot_rt_alignment`: Generate RT alignment plots after FirstPass (default: false)
- `use_robust_fitting`: Use fit_irt_model() vs simple UniformSpline (default: true for consistency)

### 2. Add Fields to FirstPassSearchParameters

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

Add to `FirstPassSearchParameters` struct (after line 115):
```julia
plot_rt_alignment::Bool
use_robust_fitting::Bool
```

Extract in constructor from `params.first_search.irt_mapping`.

### 3. Create Wrapper Function for fit_irt_model

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

Add new function before `map_retention_times!`:
```julia
"""
    fit_irt_model_firstpass(psms::Arrow.Table, params::FirstPassSearchParameters)

Wrapper around ParameterTuningSearch.fit_irt_model() for FirstPassSearch.
Converts FirstPassSearchParameters to format expected by fit_irt_model().

# Arguments
- `psms`: Arrow.Table with columns: rt, irt_predicted
- `params`: FirstPassSearchParameters

# Returns
Tuple: (rt_model, valid_rt, valid_irt, irt_mad)
"""
function fit_irt_model_firstpass(
    psms::Arrow.Table,
    params::FirstPassSearchParameters
)
    # Convert to DataFrame for fit_irt_model
    df = DataFrame(
        rt = psms[:rt],
        irt_predicted = psms[:irt_predicted]
    )

    # Create minimal ParameterTuningSearchParameters-like object
    # Or call fit_irt_model directly with required parameters
    # Need to pass: lambda_penalty, ransac_threshold, min_psms, knots config

    # Call shared fitting function
    return fit_irt_model(params, df)
end
```

**Note**: This requires making `fit_irt_model()` accessible to FirstPassSearch, or extracting it to a shared utilities module.

### 4. Update map_retention_times!

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl` (lines 345-418)

Replace lines 363-381 with:
```julia
try
    # Filter for high-confidence PSMs
    best_hits = psms[:prob] .> params.min_prob_for_irt_mapping

    if params.use_robust_fitting
        # Use robust fitting with RANSAC/penalty
        best_psms_df = DataFrame(
            rt = psms[:rt][best_hits],
            irt_predicted = psms[:irt_predicted][best_hits]
        )

        rt_model, valid_rt, valid_irt, irt_mad = fit_irt_model_firstpass(best_psms_df, params)

        # Store RT to iRT model
        setRtIrtMap!(search_context, rt_model, ms_file_idx)

        # Fit inverse model (iRT to RT)
        irt_to_rt_df = DataFrame(
            rt = valid_irt,  # Swap: fit iRT as input
            irt_predicted = valid_rt  # Swap: fit RT as output
        )
        irt_model, _, _, _ = fit_irt_model_firstpass(irt_to_rt_df, params)
        setIrtRtMap!(search_context, irt_model, ms_file_idx)

        # Optionally generate plots
        if params.plot_rt_alignment
            plot_rt_alignment_firstpass(
                valid_rt, valid_irt, rt_model, ms_file_idx,
                getDataOutDir(search_context)
            )
        end
    else
        # Use simple UniformSpline (legacy behavior)
        best_rts = psms[:rt][best_hits]
        best_irts = psms[:irt_predicted][best_hits]

        irt_to_rt_spline = UniformSpline(best_rts, best_irts, 3, 5)
        rt_to_irt_spline = UniformSpline(best_irts, best_rts, 3, 5)

        setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
        setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)
    end
catch e
    # Existing error handling...
end
```

### 5. Add RT Alignment Plotting Function

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

Add new function:
```julia
"""
    plot_rt_alignment_firstpass(rt, irt, rt_model, ms_file_idx, output_dir)

Generate RT alignment diagnostic plot for FirstPassSearch.

# Arguments
- `rt`: Observed retention times
- `irt`: Predicted indexed retention times
- `rt_model`: Fitted RT conversion model
- `ms_file_idx`: MS file index
- `output_dir`: Output directory path
"""
function plot_rt_alignment_firstpass(
    rt::Vector{Float32},
    irt::Vector{Float32},
    rt_model::RtConversionModel,
    ms_file_idx::Int,
    output_dir::String
)
    # Create plot similar to ParameterTuning RT plots
    p = scatter(rt, irt,
        xlabel="Observed RT (min)",
        ylabel="Predicted iRT",
        title="FirstPass RT Alignment (File $ms_file_idx)",
        legend=false,
        markersize=2,
        alpha=0.5
    )

    # Add fitted model line
    rt_range = range(minimum(rt), maximum(rt), length=100)
    fitted_irt = [rt_model(r) for r in rt_range]
    plot!(p, rt_range, fitted_irt, color=:red, linewidth=2)

    # Save to FirstPass RT alignment folder
    rt_plot_folder = joinpath(output_dir, "qc_plots", "rt_alignment_plots", "firstpass")
    !isdir(rt_plot_folder) && mkpath(rt_plot_folder)

    savefig(p, joinpath(rt_plot_folder, "file_$(ms_file_idx)_rt_alignment.pdf"))
end
```

---

## Implementation Strategy

### Option A: Extract fit_irt_model to Shared Utilities (RECOMMENDED)

Move `fit_irt_model()` from ParameterTuningSearch/utils.jl to a shared location:
- Create `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
- Move `fit_irt_model()` there
- Both ParameterTuningSearch and FirstPassSearch can import it
- Requires minimal parameter passing (just RT/iRT data + config)

**Benefits:**
- Single source of truth for RT fitting
- Easier to maintain and test
- Can be reused by other search methods

**Implementation:**
1. Create new file `CommonSearchUtils/rt_alignment_utils.jl`
2. Move `fit_irt_model()` and helper functions (`fit_linear_irt_model`, etc.)
3. Make the function accept parameters directly instead of ParameterTuningSearchParameters:
   ```julia
   function fit_irt_model(
       psms::DataFrame;
       lambda_penalty::Float32 = Float32(0.1),
       ransac_threshold::Int = 1000,
       min_psms::Int = 10,
       spline_degree::Int = 3,
       max_knots::Int = 7,
       outlier_threshold::Float32 = Float32(5.0)
   )
   ```
4. Update both ParameterTuningSearch and FirstPassSearch to call it

### Option B: Duplicate fit_irt_model in FirstPassSearch

Copy the function to FirstPassSearch/utils.jl.

**Drawbacks:**
- Code duplication
- Maintenance burden (2 places to update)
- Not recommended

---

## Configuration Migration

### Breaking Changes
None - new parameters have defaults that maintain current behavior.

### Behavior Changes (if use_robust_fitting=true)
- RT fitting uses RANSAC for <1000 PSMs
- Adaptive knot selection instead of fixed 5 knots
- Outlier removal improves RT alignment quality
- Consistent methodology across ParameterTuning and FirstPass

---

## Testing Plan

1. **Test with use_robust_fitting=false**: Verify legacy behavior preserved
2. **Test with use_robust_fitting=true**: Verify improved RT alignment
3. **Test with plot_rt_alignment=true**: Verify plots generated
4. **Test with <1000 PSMs**: Verify RANSAC is used
5. **Test with >=1000 PSMs**: Verify standard fitting is used
6. **Compare RT alignment quality**: Before vs after on real dataset

---

## Timeline Estimate

| Task | Time | Cumulative |
|------|------|------------|
| Extract fit_irt_model to shared utils | 20 min | 20 min |
| Add configuration parameters | 10 min | 30 min |
| Update FirstPassSearchParameters | 10 min | 40 min |
| Create fit_irt_model_firstpass wrapper | 15 min | 55 min |
| Update map_retention_times! | 30 min | 85 min |
| Add plotting function | 20 min | 105 min |
| Testing | 30 min | 135 min |
| **Total** | **~2.25 hours** | |

---

## Commit Strategy

**Commit 1**: Extract fit_irt_model to shared utilities
```bash
git commit -m "Extract RT model fitting to shared utilities

- Move fit_irt_model() from ParameterTuningSearch to CommonSearchUtils
- Create rt_alignment_utils.jl for shared RT fitting functions
- Refactor to accept parameters directly instead of struct
- Update ParameterTuningSearch to use shared function
- Enables reuse in FirstPassSearch and future search methods

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

**Commit 2**: Use robust RT fitting in FirstPassSearch
```bash
git commit -m "Use unified RT model fitting in FirstPassSearch

- Replace simple UniformSpline with fit_irt_model()
- Add configuration: plot_rt_alignment, use_robust_fitting
- Implement optional RT alignment plotting for FirstPass
- Consistent methodology with ParameterTuning phase

Benefits:
- RANSAC for limited data (<1000 PSMs)
- Adaptive knot selection (1 per 100 PSMs)
- Outlier removal via MAD threshold
- Better RT alignment quality

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## Risks and Mitigation

### Risk 1: Different PSM filtering between phases
**Context**: ParameterTuning uses all PSMs, FirstPass filters by probability threshold.

**Mitigation**: fit_irt_model() already handles different PSM counts gracefully with adaptive knots and conditional RANSAC.

### Risk 2: Performance impact
**Context**: fit_irt_model() is more complex than simple UniformSpline.

**Mitigation**:
- Make it optional with `use_robust_fitting` flag
- Performance difference is negligible for RT fitting (not in critical path)
- Better RT alignment may improve downstream performance

### Risk 3: Breaking changes in RT models
**Context**: Different fitting may produce different RT models.

**Mitigation**:
- Keep `use_robust_fitting=false` option for legacy behavior
- Test thoroughly before making default=true
- Document changes in release notes

---

## Future Enhancements

1. **Unified RT plotting across all phases**: Standardize plot generation
2. **RT model diagnostics**: Add metrics for model quality (R², residual distribution)
3. **Cross-file RT model borrowing**: Use robust files to help failed files (like ParameterTuning does)
4. **RT model persistence**: Save/load RT models for faster reprocessing
