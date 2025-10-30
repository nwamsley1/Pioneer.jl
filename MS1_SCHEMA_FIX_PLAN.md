# MS1 Feature Schema Standardization Fix Plan

## Problem Statement

The search pipeline is failing during ScoringSearch with an `ArgumentError: mismatched schemas between different arrow batches`. This occurs when Arrow.jl attempts to load multiple PSM files that have inconsistent schemas for MS1-related features.

## Error Analysis

### Error Location
- **File**: `score_psms.jl:505`
- **Function**: `load_psms_for_lightgbm()`
- **Line**: `Arrow.Table(file_paths)` - attempts to load multiple arrow files simultaneously

### Root Cause

SecondPassSearch creates different Arrow schemas depending on whether MS1 PSMs were found for each file:

**Schema A: Files WITH MS1 PSMs** (from MS1 scoring + join)
```julia
# Actual MS1 scoring produces these columns with specific types:
:rt_ms1                         => Float32
:weight_ms1                     => Float32
:gof_ms1                        => Float16 (HALF precision)
:max_matched_residual_ms1       => Float16
:max_unmatched_residual_ms1     => Float16
:fitted_spectral_contrast_ms1   => Float16
:error_ms1                      => Float16
:m0_error_ms1                   => Float16
:n_iso_ms1                      => UInt8
:big_iso_ms1                    => UInt8
:rt_max_intensity_ms1           => Float32
:rt_diff_max_intensity_ms1      => Float32

# Plus additional columns from join:
:m0_ms1                         => Bool
:spectral_contrast_ms1          => Float16
:matched_ratio_ms1              => Float16
:ms_file_idx_ms1                => UInt32
:scan_idx_ms1                   => UInt32
:pair_idx_ms1                   => UInt32
```

**Schema B: Files WITHOUT MS1 PSMs** (manually created)
```julia
# Only 12 columns, all Float32:
:rt_ms1                         => Float32
:weight_ms1                     => Float32
:gof_ms1                        => Float32  # Should be Float16
:max_matched_residual_ms1       => Float32  # Should be Float16
:max_unmatched_residual_ms1     => Float32  # Should be Float16
:fitted_spectral_contrast_ms1   => Float32  # Should be Float16
:error_ms1                      => Float32  # Should be Float16
:m0_error_ms1                   => Float32  # Should be Float16
:n_iso_ms1                      => Float32  # Should be UInt8
:big_iso_ms1                    => Float32  # Should be UInt8
:rt_max_intensity_ms1           => Float32
:rt_diff_max_intensity_ms1      => Float32

# Missing these columns entirely:
:m0_ms1, :spectral_contrast_ms1, :matched_ratio_ms1,
:ms_file_idx_ms1, :scan_idx_ms1, :pair_idx_ms1
```

**Additional Issue: ms1_ms2_rt_diff**
- One schema shows `Int64` type for this field (should be Float32)
- This is calculated after the MS1 join at line 459
- Type inconsistency suggests a bug in the calculation or data type promotion

## Current Code Location

**File**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

```julia
# Lines 427-465
if size(ms1_psms, 1) > 0
    # Path A: Join with actual MS1 PSMs
    psms = leftjoin(psms, ms1_psms, on = :precursor_idx,
                   makeunique = true, renamecols = "" => "_ms1")
    ms1_cols = filter(col -> endswith(String(col), "_ms1"), names(psms))
    miss_mask = ismissing.(psms[!, ms1_cols[1]])
else
    # Path B: Manually create columns (INCOMPLETE SCHEMA)
    ms1_cols = [
        :rt_ms1, :weight_ms1, :gof_ms1, :max_matched_residual_ms1,
        :max_unmatched_residual_ms1, :fitted_spectral_contrast_ms1,
        :error_ms1, :m0_error_ms1, :n_iso_ms1, :big_iso_ms1,
        :rt_max_intensity_ms1, :rt_diff_max_intensity_ms1
    ]
    miss_mask = trues(size(psms, 1))
    for col in ms1_cols
        psms[!, col] = -1*ones(Float32, size(psms, 1))  # ALL Float32
    end
end

# Common processing
for col in ms1_cols
    psms[!, col] = coalesce.(psms[!, col], zero(nonmissingtype(eltype(psms[!, col]))))
    disallowmissing!(psms, col)
end

# Calculate MS1-MS2 RT difference (line 459)
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = ifelse.(psms[!,:rt_ms1] .== -1,
                      -1,
                      abs.(rt_to_irt_model.(psms[!,:rt]) .- rt_to_irt_model.(psms[!,:rt_ms1])))
```

## Investigation Questions

Before implementing the fix, we need to determine:

1. **What is the canonical MS1 schema?**
   - Examine MS1 scoring output to get definitive column list and types
   - Check `parseMs1Psms()` output schema
   - Verify which columns come from `leftjoin` operation

2. **Why does ms1_ms2_rt_diff become Int64?**
   - Check if certain RT model outputs return Int64
   - Verify type promotion rules in the calculation
   - May need explicit Float32 conversion

3. **Are all MS1 columns necessary for downstream analysis?**
   - Check model_config.jl feature sets
   - Determine which columns are actually used in scoring
   - Consider removing unused columns vs. standardizing all

4. **What are appropriate sentinel values?**
   - Current: -1 for Float32 columns
   - Need sentinel for UInt8 (0? UInt8(255)?)
   - Need sentinel for Bool (false?)
   - Need sentinel for UInt32 (UInt32(0)?)

## Proposed Solution

### Approach: Complete Schema Definition

Create a standardized MS1 schema that ALWAYS produces identical columns and types, regardless of whether MS1 PSMs were found.

#### Step 1: Define Canonical MS1 Schema

Create a constant or function defining the complete MS1 schema:

```julia
# Somewhere in SecondPassSearch.jl or utils.jl
const MS1_FEATURE_SCHEMA = [
    (:rt_ms1, Float32, Float32(-1)),
    (:weight_ms1, Float32, Float32(-1)),
    (:gof_ms1, Float16, Float16(-1)),
    (:max_matched_residual_ms1, Float16, Float16(-1)),
    (:max_unmatched_residual_ms1, Float16, Float16(-1)),
    (:fitted_spectral_contrast_ms1, Float16, Float16(-1)),
    (:error_ms1, Float16, Float16(-1)),
    (:m0_error_ms1, Float16, Float16(-1)),
    (:n_iso_ms1, UInt8, UInt8(0)),
    (:big_iso_ms1, UInt8, UInt8(0)),
    (:rt_max_intensity_ms1, Float32, Float32(-1)),
    (:rt_diff_max_intensity_ms1, Float32, Float32(-1)),
    (:m0_ms1, Bool, false),
    (:spectral_contrast_ms1, Float16, Float16(-1)),
    (:matched_ratio_ms1, Float16, Float16(-1)),
    (:ms_file_idx_ms1, UInt32, UInt32(0)),
    (:scan_idx_ms1, UInt32, UInt32(0)),
    (:pair_idx_ms1, UInt32, UInt32(0)),
]
```

#### Step 2: Modify SecondPassSearch.jl

Replace lines 427-455 with new logic:

```julia
# Initialize all MS1 columns with sentinel values FIRST
for (col_name, col_type, sentinel_value) in MS1_FEATURE_SCHEMA
    psms[!, col_name] = fill(sentinel_value, size(psms, 1))
end

# If we have MS1 data, update the values
if size(ms1_psms, 1) > 0
    # Join to get MS1 features
    joined_data = leftjoin(
        select(psms, :precursor_idx),
        ms1_psms,
        on = :precursor_idx,
        makeunique = true
    )

    # Update PSM columns with MS1 values where available
    for (col_name, col_type, sentinel_value) in MS1_FEATURE_SCHEMA
        joined_col = Symbol(string(col_name))
        if hasproperty(joined_data, joined_col)
            # Update non-missing values
            for i in 1:nrow(psms)
                if !ismissing(joined_data[i, joined_col])
                    psms[i, col_name] = col_type(joined_data[i, joined_col])
                end
            end
        end
    end

    miss_mask = .!ismissing.(joined_data[!, first(MS1_FEATURE_SCHEMA)[1]])
else
    miss_mask = falses(size(psms, 1))
end

# No need for coalesce/disallowmissing since we pre-filled with non-missing values
```

#### Step 3: Fix ms1_ms2_rt_diff Type

Ensure explicit Float32 conversion:

```julia
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== -1,
                      -1,
                      abs.(rt_to_irt_model.(psms[!,:rt]) .- rt_to_irt_model.(psms[!,:rt_ms1]))))
```

#### Step 4: Update ms1_features_missing Column

Ensure it's added correctly:

```julia
psms[!, :ms1_features_missing] = miss_mask
```

## Testing Plan

### Unit Testing
1. Test schema generation for files WITH MS1 PSMs
2. Test schema generation for files WITHOUT MS1 PSMs
3. Verify Arrow.Table() can load mixed files
4. Verify sentinel values don't affect downstream scoring

### Integration Testing
1. Run complete pipeline on the failing dataset
2. Verify no schema mismatch errors
3. Check that scoring models handle sentinel values correctly
4. Validate final results match expected output

### Validation Checks
```julia
# Check all files have same schema
files = readdir(second_pass_folder, join=true)
schemas = [Arrow.Schema(Arrow.Table(f)) for f in files]
@assert all(s == first(schemas) for s in schemas) "Schemas don't match!"

# Check MS1 column types
for col_name, col_type, _ in MS1_FEATURE_SCHEMA
    @assert eltype(psms[!, col_name]) == col_type "Type mismatch for $col_name"
end
```

## Implementation Order

1. **Investigation Phase** (estimate: 1 hour)
   - Read MS1 scoring code to understand output schema
   - Check actual Arrow files to see what columns/types exist
   - Verify feature usage in model_config.jl

2. **Schema Definition** (estimate: 30 min)
   - Create MS1_FEATURE_SCHEMA constant
   - Document rationale for each column and sentinel value

3. **Code Modification** (estimate: 1-2 hours)
   - Refactor SecondPassSearch.jl lines 427-465
   - Add explicit type conversions
   - Update comments

4. **Testing** (estimate: 1-2 hours)
   - Run on failing dataset
   - Verify schema consistency
   - Check downstream effects

5. **Documentation** (estimate: 30 min)
   - Update CLAUDE.md if needed
   - Add comments explaining schema standardization
   - Document sentinel value meanings

## Risks and Mitigation

### Risk 1: Unknown MS1 columns
**Mitigation**: Parse actual MS1 output before defining schema

### Risk 2: Sentinel values affect scoring
**Mitigation**: Review feature filtering in ScoringSearch to ensure sentinels are handled

### Risk 3: Performance degradation
**Mitigation**: Pre-filling columns is O(n), same as current approach

### Risk 4: Breaking changes to existing data
**Mitigation**: This only affects new runs; document schema version if needed

## Alternative Approaches Considered

### Alternative 1: Type Conversion Layer
Add conversion logic after file loading in ScoringSearch
- **Pros**: Minimal changes to SecondPassSearch
- **Cons**: Bandaid fix, doesn't prevent future inconsistencies

### Alternative 2: Schema-aware Arrow Writing
Use Arrow's schema parameter when writing files
- **Pros**: Enforces schema at write time
- **Cons**: Requires defining schema separately, more complex

### Alternative 3: Drop MS1 Features for Files Without MS1
Only include MS1 columns in files that have MS1 data
- **Pros**: Minimal memory usage
- **Cons**: Requires ScoringSearch to handle optional columns, complex

## Questions for Review

1. **Should we investigate the MS1 schema first, or proceed with assumptions?**
   - Option A: Examine actual MS1 output before proceeding
   - Option B: Use error message schema as ground truth

2. **What sentinel values are appropriate for each type?**
   - Current proposal: -1 for floats, 0 for unsigned ints, false for bool
   - Alternative: Use missing values and handle in ScoringSearch?

3. **Should we validate all files have same schema after writing?**
   - Add assertion/test after SecondPassSearch completes?
   - Add to pipeline validation step?

4. **Is this the right level of fix, or should we look deeper?**
   - Current: Fix SecondPassSearch output
   - Alternative: Fix Arrow.Table loading to handle schema mismatches?

## Next Steps

1. **User Review**: Review this plan and provide feedback
2. **Investigation**: Examine MS1 scoring output schema
3. **Implementation**: Execute fix based on approved plan
4. **Testing**: Validate fix on failing dataset
5. **Documentation**: Update relevant docs

---

**Created**: 2025-10-29
**Status**: Awaiting Review
**Assignee**: Claude
**Reviewer**: Nathan Wamsley
