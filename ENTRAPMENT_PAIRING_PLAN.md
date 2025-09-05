# Entrapment Pairing System Implementation Plan

## Overview
This document details the implementation of a new entrapment pairing system that creates two new columns:
- `entrapment_pair_id`: Groups original target with all its entrapment variants
- `entrapment_target_idx`: Points back to the original target row for entrapment sequences

## Problem Statement
The existing system has target-decoy pairing via `add_charge_specific_partner_columns!()`, but lacks a mechanism to link entrapment sequences back to their original targets. This new system creates entrapment groups where:

1. **Original targets** (`entrapment_group_id == 0`): Get unique `entrapment_pair_id`, `entrapment_target_idx` points to self
2. **Entrapment sequences** (`entrapment_group_id > 0`): Share `entrapment_pair_id` with their target, `entrapment_target_idx` points to original target
3. **Decoys** (`is_decoy == true`): Both columns set to `missing`

## Key Design Decisions

### Pairing Key: `(base_pep_id, precursor_charge)`
- **Why not base_entrap_id**: We want entrapment sequences to link back to their original target, so we use the peptide-level ID
- **Why not base_seq_id**: We need modification-specific pairing, so peptide-level (post-modification) is correct
- **Why include charge**: Same peptide with different charges should be paired separately

### Column Types
- `entrapment_pair_id`: `Union{UInt32, Missing}` - missing for decoys
- `entrapment_target_idx`: `Union{UInt32, Missing}` - missing for decoys, points to row index for non-decoys

## Implementation Details

### 1. Function: `add_entrapment_partner_columns!(df::DataFrame)`

**Location**: `/src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`

**Algorithm**:
1. **Initialize columns** with missing values
2. **First pass**: Create lookup dictionary `(base_pep_id, precursor_charge) -> row_index` for original targets only
3. **Second pass**: 
   - Skip decoys (leave as missing)
   - Original targets: Assign new pair_id, target_idx = self
   - Entrapment sequences: Use target's pair_id, target_idx = target row

**Key Logic**:
```julia
if row.entrapment_group_id == 0
    # Original target
    entrapment_pair_ids[idx] = new_pair_id
    entrapment_target_idxs[idx] = UInt32(idx)  # Points to self
else
    # Entrapment sequence
    target_idx = target_lookup[(row.base_pep_id, row.precursor_charge)]
    entrapment_pair_ids[idx] = target_pair_id
    entrapment_target_idxs[idx] = target_idx  # Points to target
end
```

### 2. Pipeline Integration

**Location**: `/src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**Integration Point**: After `add_charge_specific_partner_columns!()` call

**Rationale**: 
- After all filtering is complete (ensures valid row indices)
- After target-decoy pairing (logical sequence)
- Before final output (ensures columns are included)

### 3. Expected Data Structure

**Example Output**:
```
Row | base_pep_id | charge | entrapment_group_id | decoy | entrapment_pair_id | entrapment_target_idx
----|-------------|--------|-------------------|-------|-------------------|---------------------
1   | 100         | 2      | 0                 | false | 1                 | 1                    (target)
2   | 100         | 2      | 1                 | false | 1                 | 1                    (entrap -> target)
3   | 100         | 2      | 2                 | false | 1                 | 1                    (entrap -> target)
4   | 100         | 2      | 0                 | true  | missing           | missing              (decoy)
5   | 101         | 3      | 0                 | false | 2                 | 5                    (target)
6   | 101         | 3      | 1                 | false | 2                 | 5                    (entrap -> target)
```

## Implementation Status

### ✅ Completed Tasks

1. **Function Implementation**
   - ✅ Created `add_entrapment_partner_columns!()` in `pair_decoys.jl`
   - ✅ Comprehensive documentation and logging
   - ✅ Proper handling of missing values for decoys
   - ✅ Statistics reporting for validation

2. **Pipeline Integration**
   - ✅ Added function call in `chronologer_prep.jl`
   - ✅ Positioned after target-decoy pairing
   - ✅ Columns automatically included in final DataFrame

3. **Testing**
   - ✅ Code compilation successful
   - ✅ No syntax errors or type issues
   - ✅ Ready for integration testing

### 🔄 Next Steps (Post-Implementation)

1. **Integration Testing**
   - Test with keap1 dataset
   - Verify column creation and values
   - Validate entrapment linking logic

2. **Validation Checks**
   - Ensure all original targets have self-referencing `entrapment_target_idx`
   - Confirm entrapment sequences link to correct targets
   - Verify decoys have missing values
   - Check entrapment groups share same `entrapment_pair_id`

## Technical Specifications

### Data Flow
```
Original Targets (entrapment_group_id=0, decoy=false)
↓
Create lookup: (base_pep_id, precursor_charge) → target_row_index
↓
Assign unique entrapment_pair_id
↓
Set entrapment_target_idx = self
↓
Entrapment Sequences (entrapment_group_id>0, decoy=false)
↓
Lookup target using (base_pep_id, precursor_charge)
↓
Use target's entrapment_pair_id
↓
Set entrapment_target_idx = target_row_index
↓
Decoys (decoy=true)
↓
Leave both columns as missing
```

### Error Handling
- **Missing targets**: Debug message logged, entrapment sequence gets missing values
- **Duplicate targets**: Uses first occurrence in lookup dictionary
- **Type safety**: All indices converted to UInt32, missing handled explicitly

### Performance Considerations
- **Two-pass algorithm**: O(n) first pass for lookup creation, O(n) second pass for assignment
- **Memory efficient**: Pre-allocate vectors, use efficient data structures
- **Dictionary lookup**: O(1) average case for target finding

## Files Modified

1. **`/src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`**
   - Added `add_entrapment_partner_columns!()` function
   - ~94 lines of new code including documentation

2. **`/src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`**
   - Added function call after target-decoy pairing
   - Updated pipeline integration comments

## Success Criteria

- [x] Function compiles without errors
- [x] Pipeline integration successful
- [ ] Original targets have self-referencing `entrapment_target_idx`
- [ ] Entrapment sequences correctly link to their targets  
- [ ] Decoys have `missing` values for both columns
- [ ] All entrapment groups share the same `entrapment_pair_id`
- [ ] Library builds successfully with new columns in output

## Usage Example

```julia
# The function is automatically called during library building
BuildSpecLib("/path/to/params.json")

# The resulting precursors_table.arrow will contain the new columns:
# - entrapment_pair_id
# - entrapment_target_idx
```

## Future Enhancements

1. **Validation Functions**: Create helper functions to verify entrapment pairing correctness
2. **Statistics Export**: Add detailed pairing statistics to build reports
3. **Configuration Options**: Allow disabling entrapment pairing if not needed
4. **Performance Optimization**: Consider batch processing for very large datasets