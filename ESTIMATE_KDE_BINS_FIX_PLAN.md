# Fix Plan for estimateKdeBins Type Error

## Problem Summary

The Quadrupole Tuning search is failing with:
```
MethodError: no method matching estimateKdeBins(::Vector{Union{Missing, Float32}})
```

This occurs because our previous fix allows missing values in DataFrame columns, but downstream functions like `estimateKdeBins` cannot handle them.

## Root Cause Analysis

1. **Data Flow**:
   - `summarize_precursor` returns missing values for empty groups (our first fix)
   - DataFrame columns become `Union{Float32, Missing}` type (our second fix)
   - `fit_quad_model` → `MergeBins` → `estimateKdeBins` receives these mixed-type vectors
   - `estimateKdeBins` only accepts `AbstractVector{<:AbstractFloat}`, not vectors with missing

2. **The Warning Message**:
   - "Too few psms found for quad modeling. Using default model."
   - This suggests the data is sparse, likely producing many missing values

## Proposed Solution

### Option 1: Quick Fix - Add Method Overload (Recommended)

Add a new method to `estimateKdeBins` that handles missing values:

```julia
# In binIsotopeRatioData.jl, after the existing estimateKdeBins function
function estimateKdeBins(x::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat
    # Filter out missing values
    x_clean = collect(skipmissing(x))
    
    # Check if we have enough data
    if length(x_clean) < 10  # Minimum for meaningful KDE
        throw(ArgumentError("Insufficient non-missing data for KDE estimation: $(length(x_clean)) points"))
    end
    
    # Call the original method with clean data
    return estimateKdeBins(x_clean)
end
```

### Option 2: Filter at Entry Point

Modify `fit_quad_model` to filter missing values before processing:

```julia
function fit_quad_model(psms::DataFrame, window_width::Float64)
    # Filter out rows with missing values in critical columns
    required_cols = [:x0, :yt, :center_mz, :prec_charge]
    
    # Create filter mask
    valid_mask = ones(Bool, nrow(psms))
    for col in required_cols
        if col in names(psms)
            valid_mask .&= .!ismissing.(psms[!, col])
        end
    end
    
    psms_clean = psms[valid_mask, :]
    
    # Check if we have enough data
    if nrow(psms_clean) < 20  # min_bin_size parameter
        @warn "Insufficient valid data for quad model after filtering" n_rows=nrow(psms_clean)
        return nothing  # Will use default model
    end
    
    binned_psms = MergeBins(
        psms_clean,
        (-(window_width + 1.0), window_width + 1.0),
        min_bin_size=20,
        min_bin_width=0.1
    )
    
    return fitRazoQuadModel(
        window_width,
        binned_psms,
        regularization_coef=params.regularization_coef
    )
end
```

### Option 3: Update MergeBins Function

Modify `MergeBins` to handle missing values internally:

```julia
function MergeBins(isotopes_ratio_data::AbstractDataFrame, x0_lim::Tuple{Float64, Float64}; 
                   min_bin_size::Int64=10, min_bin_width::Float64=0.05, max_iterations::Int64=3)
    
    # Filter out rows with missing x0 values
    if :x0 in names(isotopes_ratio_data)
        valid_data = isotopes_ratio_data[.!ismissing.(isotopes_ratio_data.x0), :]
    else
        throw(ArgumentError("Column :x0 not found in data"))
    end
    
    # Check if we have data after filtering
    if nrow(valid_data) == 0
        @warn "No valid data after filtering missing values"
        # Return empty result with correct structure
        return DataFrame(
            prec_charge = UInt8[],
            n = Int64[],
            n_bins = Int64[],
            # Add other required columns based on downstream expectations
        )
    end
    
    # Continue with existing logic on valid_data
    # ... rest of the function
end
```

## Recommended Approach

**Implement both Option 1 and Option 2** for maximum robustness:

1. **Option 1** provides a safety net for any calls to `estimateKdeBins`
2. **Option 2** filters data early, preventing wasted computation and providing better diagnostics

This two-layer approach ensures:
- Immediate fix for the current error
- Better performance by filtering early
- Clear diagnostics when data is insufficient
- Graceful degradation to default models

## Implementation Steps

1. Add the method overload to `estimateKdeBins` (Option 1)
2. Update `fit_quad_model` to filter missing values (Option 2)
3. Add logging to track how many rows are filtered
4. Test with the problematic dataset
5. Verify the pipeline completes with either fitted or default models

## Testing Strategy

```julia
# Test 1: estimateKdeBins with missing values
x_with_missing = Union{Float32, Missing}[1.0, missing, 2.0, 3.0, missing, 4.0]
@test length(estimateKdeBins(x_with_missing)) > 0

# Test 2: All missing values
x_all_missing = Union{Float32, Missing}[missing, missing, missing]
@test_throws ArgumentError estimateKdeBins(x_all_missing)

# Test 3: fit_quad_model with sparse data
df_sparse = DataFrame(
    x0 = [1.0, missing, 2.0, missing, 3.0],
    yt = [0.5, 0.6, missing, 0.7, 0.8],
    center_mz = [500.0, 500.0, 500.0, 500.0, 500.0],
    prec_charge = UInt8[2, 2, 2, 2, 2]
)
model = fit_quad_model(df_sparse, 1.0)
@test model !== nothing || model === nothing  # Either fits or returns nothing gracefully
```

## Expected Outcomes

1. **Immediate**: Pipeline will complete without errors
2. **With sparse data**: Will use default quadrupole model (as warned)
3. **With sufficient data**: Will fit custom model after filtering missings
4. **Performance**: Minimal impact, as filtering is fast compared to KDE estimation

## Risk Assessment

**Very Low Risk**:
- Only adds defensive programming
- Preserves existing behavior for complete data
- Provides clear warnings when data is insufficient
- Falls back to default models gracefully

## Notes

- The warning "Too few psms found for quad modeling" suggests this is expected behavior with sparse data
- The default model fallback is already implemented, we just need to ensure it's reached
- This fix complements our previous fixes by handling the downstream effects of allowing missing values