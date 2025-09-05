# Fix Plan: Make pair_id Independent from base_pep_id

## Problem Statement

After fixing `base_pep_id` to be preserved (not unique per modification variant), the `pair_id` assignment broke because it was based on `base_pep_id`. This causes:

1. **Multiple modification variants sharing same pair_id** - breaks target-decoy pairing
2. **Ridiculously high pair_id values** (like 4294959104) for unpaired decoys
3. **Pairing logic confusion** - can't distinguish between modification variants

## Current Issue Location

**File:** `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Line 567:**
```julia
_pair_id[n] = get_base_pep_id(peptide)  # PROBLEM: Now multiple variants share same pair_id
```

## Root Cause Analysis

### Before base_pep_id fix:
- `base_pep_id` was unique per modification variant (1, 2, 3, 4, 5...)
- `pair_id = base_pep_id` worked because each variant had unique ID
- Target-decoy pairing worked correctly

### After base_pep_id fix:
- `base_pep_id` is preserved from original target (same for all modification variants)
- `pair_id = base_pep_id` means multiple variants get same pair_id
- Target-decoy pairing breaks - can't distinguish modification variants

## Solution Options

### Option 1: Use base_prec_id for pair_id (RECOMMENDED)

**Concept:** `base_prec_id` already behaves like the old `base_pep_id` - it's unique per precursor variant.

**CONFIRMED from commit history analysis:**

**Before our base_pep_id changes (commit 49ce1223):**
```julia
_pair_id[n] = get_base_prec_id(peptide)  # ORIGINAL: Used base_prec_id for pairing
```

**After our base_pep_id changes (commit 22bef33d):**
```julia
_pair_id[n] = get_base_pep_id(peptide)   # BROKEN: Changed to base_pep_id
```

**Evidence from commit history:**
- **Original add_mods():** Used `get_base_pep_id(fasta_peptide)` (preserved original ID)
- **Original add_charge():** `base_prec_id += one(UInt32)` - incremented for each charge state
- **Original pair_id assignment:** Used `get_base_prec_id(peptide)` - not base_pep_id!
- **Our modification introduced the bug:** When we added `current_base_pep_id` increment logic, we mistakenly changed pair_id to use base_pep_id instead of base_prec_id

**Implementation:**
```julia
# Line 567: REVERT to original behavior:
# Change from:
_pair_id[n] = get_base_pep_id(peptide)   # CURRENT (broken)
# Back to:
_pair_id[n] = get_base_prec_id(peptide)  # ORIGINAL (correct)
```

**Advantages:**
- ✅ Minimal code change - just one line
- ✅ base_prec_id already tracks unique precursor variants
- ✅ Maintains existing pairing logic
- ✅ No new ID generation needed

**Potential Issues:**
- ⚠️ Need to verify base_prec_id is properly incremented in add_mods()
- ⚠️ May have gaps in base_prec_id sequence (less clean)

### Option 2: Create New Sequential pair_id

**Concept:** Generate completely new sequential IDs for pair_id, independent of both base_pep_id and base_prec_id.

**Implementation:**
```julia
# Add counter for unique pair_id generation
unique_pair_id = UInt32(0)
for (n, peptide) in enumerate(fasta_peptides)
    # ... existing code ...
    unique_pair_id += 1
    _pair_id[n] = unique_pair_id  # Sequential 1, 2, 3, 4...
end
```

**Advantages:**
- ✅ Clean sequential IDs (1, 2, 3, 4...)
- ✅ Completely independent from other ID systems
- ✅ Easy to understand and debug

**Disadvantages:**
- ❌ More code changes
- ❌ Creates third ID system to maintain

### Option 3: Hash-based pair_id

**Concept:** Create pair_id based on hash of (base_pep_id + modifications + charge).

**Implementation:**
```julia
# Create unique hash for each variant
mod_string = getModString(get_structural_mods(peptide))
charge = get_charge(peptide)
base_id = get_base_pep_id(peptide)
_pair_id[n] = hash((base_id, mod_string, charge)) % UInt32
```

**Advantages:**
- ✅ Deterministic - same input = same pair_id
- ✅ Captures all variant information

**Disadvantages:**
- ❌ Risk of hash collisions
- ❌ Non-sequential IDs (harder to debug)
- ❌ More complex logic

## Recommended Implementation (Option 1)

### Step 1: Verify base_prec_id Behavior

**Check these locations to ensure base_prec_id is properly incremented:**

1. **add_mods() function** - Verify it maintains unique base_prec_id per modification
2. **add_charge() function** - Verify it creates unique base_prec_id per charge state
3. **add_decoy_sequences()** - Verify decoys get proper base_prec_id

### Step 2: Update pair_id Assignment

**File:** `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Line 567:**

```julia
# Change from:
_pair_id[n] = get_base_pep_id(peptide)  # BROKEN: multiple variants share same ID

# Change to:
_pair_id[n] = get_base_prec_id(peptide)  # FIXED: unique per precursor variant
```

### Step 3: Test the Fix

1. **Build library** and check pair_id values
2. **Verify no ridiculously high pair_id values** (like 4294959104)
3. **Check target-decoy pairing** works correctly
4. **Verify modification variants** have different pair_id values

## Expected Results After Fix

### Before Fix (Current Problem):
```
Sequence: PEPTIDE, Mods: none,     Charge: 2 → base_pep_id: 100, pair_id: 100
Sequence: PEPTIDE, Mods: [Ox@3],   Charge: 2 → base_pep_id: 100, pair_id: 100  ← PROBLEM: Same pair_id
Sequence: PEPTIDE, Mods: [Ox@5],   Charge: 2 → base_pep_id: 100, pair_id: 100  ← PROBLEM: Same pair_id
```

### After Fix (Expected):
```
Sequence: PEPTIDE, Mods: none,     Charge: 2 → base_pep_id: 100, pair_id: 1001
Sequence: PEPTIDE, Mods: [Ox@3],   Charge: 2 → base_pep_id: 100, pair_id: 1002  ← FIXED: Unique pair_id
Sequence: PEPTIDE, Mods: [Ox@5],   Charge: 2 → base_pep_id: 100, pair_id: 1003  ← FIXED: Unique pair_id
```

## Verification Steps

1. **Count unique pair_ids** - should be close to total precursor count
2. **Check for missing partners** - reduced warnings about unpaired precursors
3. **Validate ID ranges** - pair_id values should be reasonable (not billions)
4. **Test pairing logic** - each modification variant should find correct decoy partner

## Rollback Plan

If Option 1 doesn't work:
1. Revert the one-line change
2. Implement Option 2 (sequential pair_id generation)
3. Test thoroughly before proceeding

## Files to Examine/Modify

### Primary Changes:
- `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl` (line 567)

### Verification Files:
- `src/Routines/BuildSpecLib/fasta/fasta_digest.jl` (base_prec_id creation)
- `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl` (add_charge function)
- `src/Routines/BuildSpecLib/fasta/fasta_utils.jl` (add_decoy_sequences)

### Testing:
- Build keap1 library
- Check pair_id values in final precursors_table.arrow
- Verify target-decoy pairing works

This plan maintains the fixed `base_pep_id` behavior while creating proper `pair_id` values for target-decoy pairing.