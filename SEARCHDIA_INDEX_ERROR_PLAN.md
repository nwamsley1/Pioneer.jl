# SearchDIA Index Error Analysis and Fix Plan

## Error Description

**Location**: `FirstPassSearch.jl:514` in `summarize_results!` function  
**Error Type**: `BoundsError: attempt to access 2344362-element Arrow.Primitive{Float32, Vector{Float32}} at index [0x003a9112]`  
**Index Value**: 0x003a9112 = 3,837,202 (decimal)  
**Array Size**: 2,344,362 elements

## Exact Error Code

```julia
# FirstPassSearch.jl lines 510-527
precursors = getPrecursors(getSpecLib(search_context))
i = 1
for (pid, val) in pairs(precursor_dict)
    i += 1
    setPredIrt!(search_context, pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])  # ← ERROR HERE (line 514)
    partner_pid = getPartnerPrecursorIdx(precursors)[pid]
    if ismissing(partner_pid)
        continue
    end
    
    # Match-between-runs logic continues...
    if !haskey(precursor_dict, partner_pid)
        insert!(precursor_dict, partner_pid, val)
        setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])  # ← ALSO FAILS (line 524)
    else
        setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[partner_pid])  # ← ALSO FAILS (line 526)
    end
end
```

## Root Cause Analysis

### The Problem
The error occurs when SearchDIA tries to use `pid` (precursor ID from `precursor_dict`) as a **direct array index** to access the retention time array via `getIrt(...)[pid]`.

### Data Mismatch
1. **SearchDIA expectation**: `pid` values should be contiguous integers from 1 to N (array indices)
2. **Actual library data**: `pid` values come from our modified `base_pep_id` which are non-contiguous
3. **Array access**: `getIrt()` returns an array with 2,344,362 elements (fragments, not precursors)
4. **Index value**: SearchDIA is trying to access index 3,837,202 which exceeds array bounds

### Library Structure Analysis

**From our BuildSpecLib output:**
- Total precursors: 84,184
- Unique base_pep_ids: 29,024 (non-contiguous)
- base_pep_id range: 1 to 31,606
- Multiple precursors share the same base_pep_id (for pairing)

**SearchDIA expects:**
- Precursor indices to be contiguous: 1, 2, 3, ..., N
- Direct array indexing: `array[precursor_idx]` should work
- IRT array sized to match number of precursors, not fragments

## Contributing Factors

### 1. Our Modification Impact
Our fix to make `base_pep_id` unique per modification variant:
```julia
# In add_mods function
current_base_pep_id += 1  # Creates potentially large, non-contiguous IDs
```

This created `base_pep_id` values that don't correspond to array positions.

### 2. Index vs ID Confusion
- **base_pep_id**: Should be used for target-decoy pairing logic
- **precursor_idx**: Should be used for array indexing (contiguous 1..N)
- SearchDIA is using `base_pep_id` values as array indices

### 3. Array Size Mismatch
- IRT array has 2,344,362 elements (likely fragments)
- Library has 84,184 precursors
- SearchDIA tries to access precursor-based indices in fragment-sized array

## Proposed Solutions

### Option A: Fix Library Generation (Recommended)
1. **Keep base_pep_id for pairing** - Our pairing fix works correctly
2. **Add separate precursor_idx column** - Contiguous 1..N for array indexing
3. **Update SearchDIA to use precursor_idx** - For array access operations

**Implementation:**
```julia
# In BuildSpecLib after pairing fix
function add_contiguous_precursor_idx!(df::DataFrame)
    df.precursor_idx = 1:nrow(df)  # Contiguous 1..N
    return df
end
```

### Option B: Fix SearchDIA Indexing
1. **Create index mapping** - Map base_pep_id → array_position
2. **Use mapping for array access** - Translate IDs before indexing
3. **Maintain pairing logic** - Keep using base_pep_id for target-decoy pairs

**Implementation:**
```julia
# In SearchDIA
base_pep_id_to_index = Dict(id => i for (i, id) in enumerate(unique_base_pep_ids))
irt_value = getIrt(precursors)[base_pep_id_to_index[pid]]
```

### Option C: Restructure Library Format
1. **Use contiguous base_pep_id** - Renumber to be 1..N contiguous
2. **Store pairing separately** - Use pair_id mapping instead of base_pep_id matching
3. **Update pairing logic** - Use separate pairing table

## Recommended Implementation Plan

### Phase 1: Immediate Fix (Option A)
1. **Add precursor_idx column** to BuildSpecLib output
2. **Update SearchDIA** to use precursor_idx for array access
3. **Keep base_pep_id** for target-decoy pairing logic
4. **Test with keap1 library**

### Phase 2: Validation
1. **Verify pairing still works** - pair_id should still show [2] for all pairs
2. **Verify SearchDIA indexing** - No more bounds errors
3. **Check search results** - Ensure PSMs are correctly identified

### Phase 3: Cleanup (Optional)
1. **Standardize naming** - Clarify base_pep_id vs precursor_idx usage
2. **Add documentation** - Document the difference between pairing IDs and array indices
3. **Update other methods** - Ensure consistent usage throughout SearchDIA

## Implementation Details

### BuildSpecLib Changes
```julia
# After creating final DataFrame in build_fasta_df
seq_df = add_charge_specific_partner_columns!(seq_df)
seq_df.precursor_idx = 1:nrow(seq_df)  # Add contiguous indexing
return seq_df
```

### SearchDIA Changes
```julia
# Replace all instances of using base_pep_id as array index
# OLD: getIrt(precursors)[pid]  
# NEW: getIrt(precursors)[getPrecursorIdx(precursors)[pid]]
```

## Expected Outcome
- **Pairing preserved**: All pairs remain 2 members each
- **SearchDIA works**: No bounds errors during search
- **Performance maintained**: Minimal overhead from additional column
- **Backward compatibility**: Other parts of pipeline unaffected

## Testing Strategy
1. **BuildSpecLib test**: Verify library has both base_pep_id and precursor_idx
2. **Pairing test**: Confirm pair_sizes = [2] still works
3. **SearchDIA test**: Run search without bounds errors
4. **Integration test**: Full pipeline with results validation

## Risk Assessment
- **Low risk**: Adding precursor_idx is additive, doesn't break existing logic
- **Targeted fix**: Only affects array indexing, not pairing logic
- **Reversible**: Can revert by removing precursor_idx column if issues arise