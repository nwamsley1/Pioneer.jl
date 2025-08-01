# Quadrupole Tuning Type Propagation Fix Plan

## Error Summary

**Error**: `MethodError: no method matching estimateKdeBins(::Vector{Union{Missing, Float32}})`

**Root Cause**: The `estimateKdeBins` function expects `AbstractVector{<:AbstractFloat}` but receives `Vector{Union{Missing, Float32}}` due to our previous fix that allows missing values in the DataFrame columns.

**Call Chain**:
1. `summarize_precursor` can return missing values
2. DataFrame columns become `Union{Float32, Missing}`
3. `fit_quad_model` → `MergeBins` → `estimateKdeBins` receives mixed type vectors
4. `estimateKdeBins` cannot handle missing values

## Comprehensive Fix Strategy

### 1. Immediate Fix: Filter Missing Values Before Processing

Add missing value filtering at the entry point of quadrupole model fitting:

```julia
function fit_quad_model(psms::DataFrame, window_width::Float64)
    # Filter out rows with missing values in critical columns
    psms_clean = psms[
        .!ismissing.(psms.x0) .& 
        .!ismissing.(psms.yt) .& 
        .!ismissing.(psms.center_mz) .& 
        .!ismissing.(psms.prec_charge), 
        :
    ]
    
    # Check if we have enough data after filtering
    if nrow(psms_clean) < 10  # Minimum required for meaningful model
        @warn "Insufficient data for quadrupole model fitting after filtering missing values" n_rows=nrow(psms_clean)
        return nothing  # Or return a default model
    end
    
    binned_psms = MergeBins(
        psms_clean,
        (-(window_width + 1.0), window_width + 1.0),
        min_bin_size=20,
        min_bin_width=0.1
    )
    # ... rest of function
```

### 2. Add Type Conversion Helper

Create a utility function to safely extract non-missing values:

```julia
"""
    skipmissing_vector(v::AbstractVector{Union{T, Missing}}) where T

Convert a vector with possible missing values to a vector of concrete type,
skipping all missing values.

# Returns
Vector{T} containing only non-missing elements
"""
function skipmissing_vector(v::AbstractVector{Union{T, Missing}}) where T
    return collect(skipmissing(v))
end
```

### 3. Update estimateKdeBins to Handle Missing Values

Add a method that handles Union types:

```julia
# Add this method to handle Union{Missing, T} vectors
function estimateKdeBins(x::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat
    x_clean = skipmissing_vector(x)
    if isempty(x_clean)
        throw(ArgumentError("No non-missing values in input vector"))
    end
    return estimateKdeBins(x_clean)
end
```

### 4. Update MergeBins Function

Modify MergeBins to handle missing values gracefully:

```julia
function MergeBins(isotopes_ratio_data::AbstractDataFrame, x0_lim::Tuple{Float64, Float64}; kwargs...)
    # Filter out rows with missing x0 values early
    valid_data = isotopes_ratio_data[.!ismissing.(isotopes_ratio_data.x0), :]
    
    if nrow(valid_data) == 0
        # Return empty result structure
        return DataFrame(
            prec_charge = UInt8[],
            n = Int64[],
            n_bins = Int64[],
            # ... other required columns
        )
    end
    
    # Continue with existing logic on valid_data
    # ...
end
```

### 5. Systematic Type Handling Throughout Pipeline

#### 5.1 Update DataFrame Operations

Wherever we perform operations on DataFrame columns that might contain missing:

```julia
# Before
mean_x0 = mean(df.x0)

# After
mean_x0 = mean(skipmissing(df.x0))
```

#### 5.2 Update Statistical Functions

Functions that compute statistics should handle missing values:

```julia
# Update calls to std, mean, median, etc.
std_yt = std(skipmissing(group.yt))
median_x0 = median(skipmissing(group.x0))
```

### 6. Add Validation at Key Points

Add validation functions to check data integrity:

```julia
"""
    validate_quad_data(df::DataFrame) -> DataFrame

Validate and clean quadrupole tuning data, removing rows with missing values
in critical columns and ensuring minimum data requirements.
"""
function validate_quad_data(df::DataFrame)
    required_cols = [:x0, :yt, :center_mz, :prec_charge]
    
    # Check columns exist
    missing_cols = setdiff(required_cols, names(df))
    if !isempty(missing_cols)
        throw(ArgumentError("Missing required columns: $missing_cols"))
    end
    
    # Filter rows with any missing values in required columns
    valid_mask = true
    for col in required_cols
        valid_mask = valid_mask .& .!ismissing.(df[!, col])
    end
    
    df_clean = df[valid_mask, :]
    
    # Validate minimum data
    if nrow(df_clean) < 10
        @warn "Insufficient valid data for quadrupole analysis" n_rows=nrow(df_clean)
    end
    
    return df_clean
end
```

## Implementation Priority

1. **High Priority - Immediate Fixes**:
   - Add missing value filtering in `fit_quad_model`
   - Add method to `estimateKdeBins` for Union types
   - These fix the immediate error

2. **Medium Priority - Robustness**:
   - Update `MergeBins` to handle missing values
   - Add validation functions
   - Ensures pipeline handles edge cases gracefully

3. **Low Priority - Systematic Updates**:
   - Update all statistical operations to use `skipmissing`
   - Add comprehensive tests
   - Long-term maintainability

## Testing Strategy

### Unit Tests
```julia
@testset "Quadrupole functions handle missing values" begin
    # Test estimateKdeBins with missing values
    x_with_missing = Union{Float32, Missing}[1.0, missing, 2.0, 3.0, missing, 4.0]
    bins = estimateKdeBins(x_with_missing)
    @test !isempty(bins)
    
    # Test with all missing
    x_all_missing = Union{Float32, Missing}[missing, missing, missing]
    @test_throws ArgumentError estimateKdeBins(x_all_missing)
    
    # Test fit_quad_model with mixed data
    df_mixed = DataFrame(
        x0 = [1.0, missing, 2.0, 3.0],
        yt = [0.5, 0.6, missing, 0.7],
        center_mz = [500.0, 500.0, 500.0, missing],
        prec_charge = UInt8[2, 2, 2, 2]
    )
    model = fit_quad_model(df_mixed, 1.0)
    @test model !== nothing
end
```

### Integration Tests
- Run full quadrupole tuning with datasets containing missing values
- Verify results are consistent with filtered data
- Check performance impact

## Risk Assessment

**Low to Medium Risk**:
- Changes are mostly defensive programming
- Preserves existing behavior for valid data
- Adds graceful handling of edge cases
- Some performance overhead from filtering operations

## Benefits

1. **Robustness**: Pipeline handles incomplete data gracefully
2. **Clarity**: Explicit about missing value handling
3. **Debugging**: Better error messages for data issues
4. **Compatibility**: Works with our DataFrame schema changes

## Alternative Approach

If filtering missing values causes too much data loss, consider imputation:
- Use median/mean for missing numerical values
- Use mode for categorical values
- Document imputation strategy clearly

However, filtering is preferred for scientific accuracy.