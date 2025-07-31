# Quadrupole Tuning DataFrame Type Mismatch Fix Plan

## Error Summary

**Error Type**: Type conversion error - "cannot convert a value to missing for assignment"

**Location**: `src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/utils.jl:288` in `append!` operation

**When it occurs**: During DataFrame append operation when combining PSM results

## Detailed Error Analysis

### Root Cause

The error occurs because of a type mismatch between DataFrames:

1. **Initial State**: `total_psms = DataFrame()` creates an empty DataFrame
2. **First Iteration**: When `append!(total_psms, processed_psms)` is first called:
   - `processed_psms` has columns with concrete types (e.g., `:center_mz => Vector{Float32}`)
   - But some rows from `summarize_precursor` return `missing` values
   - DataFrame creates columns that can accommodate both Float32 and Missing
3. **Type Conflict**: The DataFrame column types become incompatible:
   - Some columns are typed as `Vector{Missing}` (when first rows have missing values)
   - Later append operations try to add Float32 values to Missing-only columns
   - Julia cannot convert Float32 to Missing type

### Why This Happens Now

Our previous fix to handle empty groups returns all `missing` values. If the first group processed returns missing values, the DataFrame columns get typed as `Vector{Missing}` instead of `Vector{Union{Float32, Missing}}`.

## Proposed Solution

### Option 1: Specify Column Types at Initialization (Recommended)

Create `total_psms` with the correct schema upfront:

```julia
# Instead of: total_psms = DataFrame()
total_psms = DataFrame(
    scan_idx = Int64[],
    precursor_idx = UInt32[],
    center_mz = Union{Float32, Missing}[],
    δ = Union{Float32, Missing}[],
    yt = Union{Float32, Missing}[],
    x0 = Union{Float32, Missing}[],
    x1 = Union{Float32, Missing}[],
    prec_charge = Union{UInt8, Missing}[],
    half_width_mz = Float32[]
)
```

### Option 2: Use `promote=true` in append!

```julia
append!(total_psms, processed_psms; promote=true)
```

This tells DataFrame to automatically promote column types to accommodate new data.

### Option 3: Filter Out Missing Rows

Filter out rows with missing values before appending:

```julia
valid_rows = .!ismissing.(processed_psms.center_mz)
append!(total_psms, processed_psms[valid_rows, :])
```

## Recommended Approach

**Use Option 1** - Initialize with correct schema. This is the most robust solution because:
1. It's explicit about expected column types
2. Avoids runtime type promotion overhead
3. Makes the code's intent clear
4. Handles all cases where some values might be missing

## Implementation Steps

1. Find the `total_psms = DataFrame()` initialization
2. Replace with schema-aware initialization
3. Ensure all columns that can have missing values use `Union{T, Missing}` type
4. Test with datasets that produce missing values

## Testing Strategy

### Test Cases
1. **All Missing**: Test with data that produces all missing values
2. **Mixed**: Test with data that has both missing and non-missing values
3. **No Missing**: Test with data that has no missing values
4. **Empty Groups**: Ensure our previous fix still works correctly

### Integration Test
Run the full Quadrupole Tuning search on:
- The dataset that triggered the original error
- Datasets with sparse isotope patterns
- Datasets with varying charge states

## Risk Assessment

**Low Risk** - This change:
- Only affects DataFrame initialization
- Makes column types more permissive (allows both values and missing)
- Doesn't change any algorithm logic
- Follows Julia/DataFrames best practices

## Alternative Quick Fix

If we need an immediate workaround, we can use Option 2 (promote=true) as a minimal change:
```julia
append!(total_psms, processed_psms; promote=true)
```

However, Option 1 is preferred for long-term maintainability.

## Expected Outcome

After this fix:
- DataFrame append operations will succeed regardless of missing values
- The Quadrupole Tuning search will handle all data patterns correctly
- No type conversion errors will occur
- Performance will be slightly better with pre-defined schema