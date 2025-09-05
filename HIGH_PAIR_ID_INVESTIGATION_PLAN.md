# Investigation Plan: High pair_id Values Issue

## Current Situation Summary

### What We've Fixed So Far:
1. ✅ **base_pep_id preservation**: Reduced from 29,024 to 11,076 unique IDs
2. ✅ **base_prec_id incrementation**: Added `current_base_prec_id += 1` in add_mods()  
3. ✅ **pair_id assignment**: Changed from base_pep_id to base_prec_id
4. ✅ **Partner indices**: Fixed bounds checking (max ≤ table size)

### Current Problem:
- **9,720 precursors still have very high pair_id values** (>1M, max = 4294967295 = UInt32 max)
- **All high pair_id precursors are decoys** (based on previous analysis)
- **Improvement**: Reduced from 12,817 to 9,720 (24% reduction)

## Current Data Flow Analysis

### 1. Library Building Pipeline
```
1. Digestion → base_pep_id = 1,2,3... (unique per original peptide)
2. Add entrapments → inherit base_pep_id from targets ✅
3. Add modifications → base_prec_id = 1,2,3... (unique per variant) ✅
4. Add decoys → inherit base_pep_id + base_prec_id from targets ✅  
5. Add charges → base_prec_id incremented further per charge ✅
6. Assign pair_id = base_prec_id ✅
7. Pairing function → assigns final pair_id values ❌ (PROBLEM HERE)
```

### 2. Pairing Function Flow (pair_decoys.jl)
```
Input: DataFrame with pair_id = base_prec_id
Process:
1. Create lookup: (base_pep_id, charge, is_target) → row_index
2. For each target: find decoy with same (base_pep_id, charge)
3. Assign shared pair_id to both target and decoy
4. Handle unpaired decoys with high counter values
Output: DataFrame with final pair_id values
```

## Problem Analysis

### ROOT CAUSE IDENTIFIED: base_prec_id Overwritten in add_charge()

**The Real Issue**: The `add_charge()` function is **overwriting** the carefully crafted base_prec_id values from `add_mods()` with its own incremental counter!

**Data Flow Problem**:
```
1. add_mods() → base_prec_id = 100, 101, 102... (unique per modification) ✅
2. add_charge() → base_prec_id = 1, 2, 3, 4... (overwrites previous!) ❌
3. pair_id assignment → uses corrupted base_prec_id values
4. Pairing function → can't match properly due to wrong base_prec_id
```

**Example of the Issue**:
```
After add_mods():
  PEPTIDE[none] → base_prec_id = 100
  PEPTIDE[Ox@3] → base_prec_id = 101

After add_charge() (WRONG BEHAVIOR):
  PEPTIDE[none]+2 → base_prec_id = 1  ← Should be 100!
  PEPTIDE[none]+3 → base_prec_id = 2  ← Should be 100!
  PEPTIDE[Ox@3]+2 → base_prec_id = 3  ← Should be 101!
  PEPTIDE[Ox@3]+3 → base_prec_id = 4  ← Should be 101!
```

### Evidence Supporting This Theory:
- UInt32 max value (4294967295) suggests counter overflow in pairing function
- add_charge() increments base_prec_id per charge state (line 395)
- This breaks the modification-specific pairing we intended
- Pairing function gets confused and assigns overflow values to unpaired decoys

## Detailed Investigation Steps

### Step 1: Confirm Pairing Function Behavior
**Check if pairing function overwrites base_prec_id-based pair_id:**

1. **Before pairing**: Log pair_id values going INTO the pairing function
2. **After pairing**: Log pair_id values coming OUT of the pairing function  
3. **Compare**: See if base_prec_id values are being replaced

### Step 2: Trace Counter Overflow Issue
**Investigate pair_id_counter in pairing function:**

1. **Initial value**: Check starting value of pair_id_counter
2. **Increment logic**: Verify counter doesn't overflow UInt32
3. **Unpaired handling**: Check if unpaired decoys get excessive counter values

### Step 3: Examine Pairing Logic vs. base_prec_id System
**Two possible approaches:**

**Option A**: Modify pairing function to PRESERVE base_prec_id-based pair_id
- Don't create new pair_id values
- Use existing base_prec_id values for pairing
- Only assign partner_precursor_idx

**Option B**: Fix counter overflow in current pairing system
- Keep current pairing logic
- Fix counter overflow/initialization issues
- Ensure counter stays within reasonable bounds

## Proposed Solution (UPDATED with Root Cause Fix)

### Primary Fix: Preserve base_prec_id in add_charge() Function
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

#### Change 1: Fix add_charge() to Preserve base_prec_id (IMPLEMENTED ✅)
**Lines 374-395:**

**OLD CODE**:
```julia
base_prec_id = one(UInt32)  # Creates new counter
for fasta_peptide in fasta_peptides
    for charge in range(UInt8(min_charge), UInt8(max_charge))
        # ... FastaEntry creation ...
        base_prec_id,  # Uses incremented counter
        # ...
    base_prec_id += one(UInt32)  # Overwrites modification-specific ID!
```

**NEW CODE** (FIXED):
```julia
# NOTE: Preserve base_prec_id from modification variants (don't increment per charge)
for fasta_peptide in fasta_peptides
    for charge in range(UInt8(min_charge), UInt8(max_charge))
        # ... FastaEntry creation ...
        get_base_prec_id(fasta_peptide),  # Preserve from modification variant
        # ...
    # NOTE: No increment - all charge states share same base_prec_id
```

### Expected Result from This Fix:
- **Modification variants will now have consistent base_prec_id across charge states** ✅
- **Different modifications will have different base_prec_id values** ✅  
- **Pairing function should work correctly with proper base_prec_id values** ❌ (Still needed secondary fix)
- **High pair_id count reduced but not eliminated** (9,720 → 8,276, 15% reduction)

### Secondary Fix (IMPLEMENTED ✅): Pairing Function Adjustments
**File**: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`

**Root Cause**: The pairing function was completely overwriting the carefully crafted base_prec_id-based pair_id values.

#### Implemented Solution: Preserve Existing pair_id Values
**OLD LOGIC** (lines 218, 243-247):
```julia
new_pair_ids = Vector{UInt32}(undef, n)    # Creates blank array
# ...
pair_id_counter += 1                        # Sequential counter
new_pair_ids[idx] = pair_id_counter        # Overwrites base_prec_id values!
new_pair_ids[decoy_idx] = pair_id_counter
```

**NEW LOGIC** (IMPLEMENTED):
```julia
new_pair_ids = copy(df.pair_id)            # Preserve existing base_prec_id values!
# ...
target_pair_id = df[idx, :pair_id]         # Use existing correct value
new_pair_ids[idx] = target_pair_id         # Keep base_prec_id-based pair_id
new_pair_ids[decoy_idx] = target_pair_id   # Share with decoy partner

# Unpaired decoys get safe sequential IDs from max_existing_pair_id + 1
max_existing_pair_id = maximum(df.pair_id)
unpaired_counter = max_existing_pair_id + 1
```

**Key Changes**:
1. **Initialize with existing values**: `copy(df.pair_id)` instead of blank array
2. **Preserve target pair_ids**: Use existing base_prec_id-based values
3. **Safe unpaired handling**: Start unpaired counter from max existing + 1
4. **No UInt32 overflow**: Eliminates counter overflow issues

### Expected Results:
1. **No more high pair_id values** - all stay within base_prec_id range
2. **Proper modification variant pairing** - each variant has unique pair_id
3. **Entrapment tracking preserved** - base_pep_id remains unchanged
4. **Target-decoy pairing works** - targets and decoys with same base_prec_id share pair_id

### Alternative Solution (Option B)

If Option A doesn't work, investigate and fix the counter overflow:

1. **Check counter initialization**: Ensure pair_id_counter starts at reasonable value
2. **Fix overflow**: Use UInt64 or check bounds before assignment
3. **Debug unpaired logic**: Verify why so many decoys are unpaired

## Testing Strategy (UPDATED)

### Phase 1: Test Primary Fix (add_charge() correction) - IMPLEMENTED ✅
1. **Rebuild keap1 library** with corrected add_charge() function
2. **Check pair_id distribution** - should eliminate most high values
3. **Verify base_prec_id consistency** across charge states

### Phase 2: Debug Remaining Issues (if any)
1. **If high pair_id count still >1000**: investigate pairing function issues
2. **Add debug logging** to pairing function if needed
3. **Check for edge cases** in unpaired decoy handling

### Phase 3: Implement Secondary Fix (if needed)
1. **Apply pairing function changes** if primary fix insufficient
2. **Test complete pipeline** with both fixes
3. **Verify all functionality** preserved

### Phase 4: Final Validation
1. **Check target-decoy pairing** still works correctly
2. **Verify base_pep_id preservation** for entrapments  
3. **Test SearchDIA** runs without errors
4. **Run comprehensive validation** suite
5. **Confirm all three objectives met**:
   - ✅ base_pep_id preserved for entrapment tracking
   - ✅ No high pair_id values (reasonable range)
   - ✅ Target-decoy pairing functional

## Files to Modify

### Primary Changes:
- `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl` (pairing logic)

### Debug/Testing:
- Add temporary logging to understand current behavior
- Modify validation script to check pair_id ranges

## Success Criteria

1. ✅ **All pair_id values reasonable**: Max pair_id < 100,000
2. ✅ **No UInt32 overflow**: Max pair_id ≠ 4294967295  
3. ✅ **High pair_id count near zero**: <100 instead of 9,720
4. ✅ **Target-decoy pairing preserved**: Same functionality as before
5. ✅ **base_pep_id unchanged**: Entrapment tracking still works
6. ✅ **Library builds and SearchDIA runs**: No regressions

## Risk Assessment

### Low Risk (Option A):
- Logical fix that aligns with our base_prec_id system
- Removes complex counter logic
- Preserves all pairing functionality

### Medium Risk:
- Need to ensure base_prec_id values are truly unique per variant
- Verify no edge cases in unpaired decoy handling

This plan will systematically resolve the remaining high pair_id issue while preserving all the fixes we've already implemented.