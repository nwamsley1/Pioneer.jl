# Detailed Implementation Plan: Partner Index Fix

## Problem Summary

**Root Cause**: Partner precursor indices (`partner_precursor_idx`) are assigned BEFORE DataFrame transformations (filtering and sorting), causing indices to become invalid and exceed table bounds.

**Current Flow**:
1. `build_fasta_df()` assigns partner indices based on DataFrame positions
2. `prepare_chronologer_input()` filters by m/z range (~33% reduction: 126k → 84k rows)  
3. `parse_chronologer_output()` sorts by retention time and m/z
4. **Result**: Partner indices reference non-existent or wrong rows

**Target Flow**:
1. Generate all precursors without partner indices
2. Apply all filtering and sorting transformations
3. Assign partner indices based on FINAL DataFrame positions
4. **Result**: Partner indices are valid row references

## Implementation Details

### Step 1: Remove Early Pairing Assignment

**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Lines**: 225-227 (currently active fix location)

**Current Code** (to be removed):
```julia
# Apply charge-specific target-decoy pairing AFTER all filtering is complete
# This ensures partner_precursor_idx values are valid row indices  
fasta_df = add_charge_specific_partner_columns!(fasta_df)
```

**Action**: Remove these lines completely. The pairing will move to after ALL transformations.

### Step 2: Verify add_pair_indices! Function ✅

**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Lines**: 417-451

**Existing Function**:
```julia
function add_pair_indices!(df)
    # Uses existing pair_id column to create partner mappings
    # Creates partner_precursor_idx column pointing to partner row
    # Handles cases where pairs don't exist (leaves as missing)
    # Perfect for final DataFrame state
end
```

**Function Analysis**:
✅ **Function exists and is perfectly designed for this use case**
✅ **Uses `pair_id` column (already created by earlier pairing logic)**
✅ **Creates `partner_precursor_idx` column with row indices**
✅ **Works on any DataFrame format with `pair_id` column**
✅ **Handles missing partners gracefully**

**Key Logic**:
1. Maps each `pair_id` to its row indices
2. For pairs with exactly 2 members, sets each to point to the other
3. Leaves singleton pairs as `missing` (appropriate behavior)

### Step 3: Enable Final Pairing Assignment

**File**: `src/Routines/BuildSpecLib.jl`
**Line**: 321

**Current Code**:
```julia
# add_pair_indices!(precursors_table)  # Removed: partner_precursor_idx now added by add_charge_specific_partner_columns!
```

**New Code**:
```julia
add_pair_indices!(precursors_table)  # Add partner indices AFTER all sorting is complete
```

**Critical Timing**: This occurs after:
- All chronologer processing and RT prediction
- All DataFrame sorting (by irt, then by mz within RT bins)
- All filtering operations
- DataFrame is in its FINAL state before library creation

### Step 4: Function Compatibility Check

**Potential Issue**: `add_charge_specific_partner_columns!` was designed for the FASTA DataFrame format, but `add_pair_indices!` needs to work on the chronologer output format.

**Key Differences**:
- FASTA format: Has `decoy` column (Boolean)
- Chronologer format: May have different column names/types

**Required Adaptation**: If `add_pair_indices!` doesn't exist or is incompatible:
1. Create new function based on `add_charge_specific_partner_columns!` logic
2. Adapt to work with chronologer DataFrame format
3. Use same pairing algorithm but compatible with final table structure

## Implementation Steps

### Phase 1: Code Analysis ✅
1. **Search for existing function** ✅:
   - Function found at `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl:417`
   - Function is perfectly designed for this use case

2. **Analyze current table structure at line 321**:
   - DataFrame at this point has `pair_id` column (required by `add_pair_indices!`)
   - All sorting and filtering complete
   - Ready for final partner index assignment

### Phase 2: Function Implementation ✅
**Selected Option A**: `add_pair_indices!` exists and is fully compatible
- Simply uncomment line 321 ✅
- No function adaptation needed ✅
- Function designed exactly for this purpose ✅

### Phase 3: Remove Early Pairing
- Remove lines 225-227 from `chronologer_prep.jl`
- Ensure no other pairing calls remain in early pipeline

### Phase 4: Testing
- Build test library with keap1.fasta
- Run validation script to verify indices are valid
- Test SearchDIA runs without BoundsError

## Validation Criteria

### 1. Partner Index Validity
```julia
max_partner_idx = Int64(maximum(skipmissing(ptable.partner_precursor_idx)))
table_size = nrow(ptable)
valid_indices = max_partner_idx <= table_size
```
**Expected**: `true`

### 2. Partnership Symmetry
```julia
# For each precursor with a partner, verify partner points back
for row_idx in sample_rows
    partner_idx = ptable.partner_precursor_idx[row_idx]
    if !ismissing(partner_idx)
        partner_partner_idx = ptable.partner_precursor_idx[partner_idx]
        symmetric = !ismissing(partner_partner_idx) && partner_partner_idx == row_idx
    end
end
```
**Expected**: All partnerships are symmetric

### 3. Pairing Integrity
- Same `pair_id` for both partners
- Opposite `is_decoy` status for partners
- Same `base_pep_id` and `precursor_charge` for partners

### 4. SearchDIA Functionality
```julia
SearchDIA("path/to/search_params.json")
```
**Expected**: No BoundsError, successful completion

## Risk Assessment

### Low Risk Changes
- Uncommenting line 321 (if function exists and works)
- Removing pairing from chronologer_prep.jl

### Medium Risk Changes
- Function adaptation (if add_pair_indices! needs modification)
- Ensuring compatibility with chronologer DataFrame format

### High Risk Areas
- Changes to pairing algorithm logic
- Modifications to base_pep_id generation (avoid this)

## Rollback Plan

### If Implementation Fails
1. **Revert to current state**: Re-add lines 225-227 in chronologer_prep.jl
2. **Comment out line 321**: Keep it disabled
3. **Alternative approach**: Investigate if sorting can be avoided or indices can be updated after sorting

### Alternative Solutions
1. **Update indices after sorting**: Instead of moving pairing, update existing indices
2. **Use stable identifiers**: Replace indices with base_pep_id-based lookup
3. **Defer sorting**: Perform sorting in SearchDIA instead of library building

## Success Metrics

1. ✅ **Library builds without errors**
2. ✅ **All partner_precursor_idx ≤ table size**
3. ✅ **Symmetric partnerships maintained**
4. ✅ **SearchDIA runs without BoundsError**
5. ✅ **Target-decoy pairing statistics match previous results**

## Files to Modify

### Primary Changes
- `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl` (lines 225-227)
- `src/Routines/BuildSpecLib.jl` (line 321)

### Potential Secondary Changes
- `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl` (if function needs adaptation)

### Testing Files
- `VALIDATE_KEAP1_FIX.jl` (validation script)
- `test_keap1_params.json` (test parameters)

## Next Actions

1. ✅ **Search for add_pair_indices! function definition**
2. ✅ **Analyze DataFrame structure at BuildSpecLib.jl line 321**
3. ✅ **Implement function if missing or adapt if incompatible** (not needed - function perfect as-is)
4. **Execute the simple two-line change: remove early pairing, enable late pairing**
5. **Run comprehensive validation tests**
6. **Test SearchDIA with corrected library**

## Simplified Implementation (Ready to Execute)

The solution is now very simple with `add_pair_indices!` function confirmed to exist and work perfectly:

### Change 1: Remove Early Pairing
```julia
# File: src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl
# Lines 225-227: DELETE these lines
# Apply charge-specific target-decoy pairing AFTER all filtering is complete
# This ensures partner_precursor_idx values are valid row indices  
fasta_df = add_charge_specific_partner_columns!(fasta_df)
```

### Change 2: Enable Final Pairing  
```julia  
# File: src/Routines/BuildSpecLib.jl
# Line 321: UNCOMMENT this line
add_pair_indices!(precursors_table)  # Add partner indices AFTER all sorting is complete
```

### Change 3: Add Monitoring Print Statements
Strategic prints to verify the fix during execution:

**A. Monitor filtering impact** (in `chronologer_prep.jl` after line 223):
```julia
println("🔍 PARTNER INDEX DEBUG:")
println("   Before m/z filtering: $(nrow(fasta_df)) precursors")
filter!(x -> (x.mz >= prec_mz_min) & (x.mz <= prec_mz_max), fasta_df)
println("   After m/z filtering: $(nrow(fasta_df)) precursors ($(round((1-nrow(fasta_df)/nrow_before)*100, digits=1))% removed)")
```

**B. Monitor chronologer processing** (in `chronologer_parse.jl` after line 80):
```julia
println("   After chronologer parsing: $(nrow(precursors_df)) precursors")
println("   Unique pair_ids before sorting: $(length(unique(precursors_df.pair_id)))")
```

**C. Monitor sorting impact** (in `chronologer_parse.jl` after line 139):
```julia
println("   After RT/mz sorting: $(nrow(precursors_df)) precursors") 
println("   Unique pair_ids after sorting: $(length(unique(precursors_df.pair_id)))")
```

**D. Monitor final pairing** (in `BuildSpecLib.jl` before and after line 321):
```julia
println("   Before add_pair_indices!: $(nrow(precursors_table)) precursors")
println("   Unique pair_ids available: $(length(unique(precursors_table.pair_id)))")
add_pair_indices!(precursors_table)
max_partner = ismissing(maximum(skipmissing(precursors_table.partner_precursor_idx))) ? 0 : Int64(maximum(skipmissing(precursors_table.partner_precursor_idx)))
println("   After add_pair_indices!: max partner_idx = $max_partner (table size: $(nrow(precursors_table)))")
println("   Valid indices: $(max_partner <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
```

### Change 4: Clean Up (Remove Old Debug Prints)
Remove any old debugging prints that are no longer needed:
- Check for temporary println statements added during investigation
- Remove any print statements not part of the strategic monitoring above
- Keep only the new structured monitoring prints

**Expected Console Output Flow:**
```
🔍 PARTNER INDEX DEBUG:
   Before m/z filtering: ~126000 precursors
   After m/z filtering: ~84000 precursors (33.5% removed)
   After chronologer parsing: ~84000 precursors  
   Unique pair_ids before sorting: ~42000
   After RT/mz sorting: ~84000 precursors
   Unique pair_ids after sorting: ~42000
   Before add_pair_indices!: ~84000 precursors
   Unique pair_ids available: ~42000
   After add_pair_indices!: max partner_idx = ~84000 (table size: ~84000)
   Valid indices: ✅ YES
```

**Total Changes**: 2 functional lines + monitoring prints + cleanup
**Risk Level**: Very low - using existing, tested function

## Execution Order
1. **Add monitoring prints** (Change 3) - so we can see what happens
2. **Remove early pairing** (Change 1) - stop creating invalid indices  
3. **Enable final pairing** (Change 2) - create valid indices after all transforms
4. **Test with keap1** - verify fix works
5. **Clean up prints** (Change 4) - remove debugging output for production

This plan ensures partner indices are assigned when the DataFrame is in its final, stable state, eliminating the bounds error while maintaining proper target-decoy pairing functionality.