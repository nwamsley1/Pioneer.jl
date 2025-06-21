# MaxLFQ FileReference & TransformPipeline Integration Plan (Simplified)

## Overview

This document outlines a **simplified** plan to integrate FileReferences and TransformPipeline into MaxLFQ search. The focus is on type safety and pipeline composition **without over-optimizing the working algorithms**.

**Key Principle**: Don't fix what isn't broken - the existing MaxLFQ implementation handles memory management well.

## Key Correction: Arrow.Table Memory Usage

**Important**: `DataFrame(Arrow.Table(file_path))` is **memory-efficient** - it creates a lazy wrapper that doesn't load data into memory. Only `DataFrame(Tables.columntable(Arrow.Table(file_path)))` loads everything into memory for modification operations.

## Analysis: Current MaxLFQ Implementation is Fine

### Memory Usage Assessment

The existing MaxLFQ implementation already handles memory well:

1. **Pre-allocation (Lines 354-375)**: The batch processing approach limits memory usage to reasonable chunks
2. **MaxLFQ Algorithm**: Matrix operations are per-protein-group and scale appropriately  
3. **Batch Processing**: The existing batch_size parameter provides sufficient control

**Conclusion**: The current memory management works. Focus on integration, not optimization.

## Simplified Integration Strategy

### Phase 1: FileReference Integration (Type Safety)

**Add ProteinQuantFileReference Type:** ✅ Already implemented
```julia
struct ProteinQuantFileReference <: FileReference
    file_path::String
    schema::FileSchema
    row_count::Int
    sorted_by::Tuple{Vararg{Symbol}}
    n_protein_groups::Int
    n_experiments::Int
end
```

**Input Validation:** ✅ Already implemented
```julia
function validate_maxlfq_input(ref::PSMFileReference)
    # Check required columns and sort order
    # Implemented in FileOperations.jl
end
```

### Phase 2: Simple Pipeline Preprocessing

**Current Manual Operations (Lines 338-342):**
```julia
subdf = subdf[(
    subdf[!,:pg_qval].<=q_value_threshold
).&(subdf[!,:global_qval_pg].<=q_value_threshold
).&(subdf[!,:use_for_protein_quant])
,:]
```

**Replace With Simple Pipeline:**
```julia
# Use existing operations - no custom functions needed
preprocessing_pipeline = TransformPipeline() |>
    filter_by_multiple_thresholds([
        (:pg_qval, q_value_threshold),
        (:global_qval_pg, q_value_threshold)
    ]) |>
    filter_rows(row -> row.use_for_protein_quant)

# Apply directly to DataFrame - keep existing approach
```

### Phase 3: Keep Existing LFQ Algorithm

**No Changes to Memory Management:**
- Keep existing batch_size parameter
- Keep existing pre-allocation patterns (they work fine)
- Keep existing matrix operations in getProtAbundance
- **Only add**: FileReference validation and basic pipeline preprocessing

## Implementation Details

### MaxLFQSearch.jl Changes (Simple Integration)

```julia
# Current approach (lines 181-189)
LFQ(DataFrame(Arrow.Table(precursors_long_path)), ...)

# New approach with FileReference validation only
function summarize_results!(results, params, search_context)
    # Create FileReference for validation
    precursor_ref = PSMFileReference(precursors_long_path)
    
    # Validate input
    validate_maxlfq_input(precursor_ref)
    
    # Call existing LFQ function - no changes to algorithm
    LFQ(
        DataFrame(Arrow.Table(precursors_long_path)),  # Keep this lazy approach
        protein_long_path,
        precursor_quant_col,
        collect(getFileIdToName(getMSData(search_context))),
        params.q_value_threshold,
        batch_size = params.batch_size,
        min_peptides = params.min_peptides
    )
    
    # Create output reference for metadata tracking
    protein_ref = ProteinQuantFileReference(protein_long_path)
end
```

### Optional: Pipeline Preprocessing in LFQ Function

**Only if desired**, replace manual filtering in LFQ function:
```julia
# In LFQ function, replace lines 338-342:
# OLD:
# subdf = subdf[(subdf[!,:pg_qval].<=q_value_threshold)...]

# NEW:
preprocessing_pipeline = TransformPipeline() |>
    filter_by_multiple_thresholds([
        (:pg_qval, q_value_threshold),
        (:global_qval_pg, q_value_threshold)
    ]) |>
    filter_rows(row -> row.use_for_protein_quant)

for (desc, op) in preprocessing_pipeline.operations
    subdf = op(subdf)
end
```

**This is optional and minimal - the manual filtering works fine.**

## Benefits of This Simplified Approach

### 1. Type Safety & Validation
- **FileReference integration** - Metadata tracking and validation
- **Input validation** - Catch configuration errors early
- **Schema checking** - Ensure required columns exist

### 2. Code Quality
- **Minimal changes** - Don't break working algorithms
- **Simple integration** - Use existing pipeline operations where beneficial
- **Clear interfaces** - FileReference abstractions

### 3. Maintainability
- **Consistent patterns** - Match ScoringSearch approach where appropriate
- **No over-engineering** - Keep working memory management
- **Easy testing** - Minimal surface area for integration bugs

## Success Metrics (Simplified)

### Integration Quality
- [ ] FileReference validation catches input errors
- [ ] No performance regression from existing implementation
- [ ] Identical output results

### Code Quality  
- [ ] Clean integration with existing MaxLFQ function
- [ ] Optional pipeline preprocessing works correctly
- [ ] Comprehensive input validation

## Migration Strategy

1. **Add FileReference validation** - Integrate validation functions
2. **Optional pipeline preprocessing** - Only if it adds value
3. **Test thoroughly** - Ensure no regressions
4. **Keep existing as fallback** - Don't remove working code

---

*Plan updated: January 2025*
*Key insight: Don't fix what isn't broken - the existing MaxLFQ implementation works well*
*Approach: Focus on type safety and validation, not optimization*