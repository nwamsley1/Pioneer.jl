# Plan: Move Modifications Before Entrapment and Handle Mods in Entrapment Sequences

## Problem Analysis

### Current Pipeline Order
```
1. combine_shared_peptides
2. assign_base_seq_ids
3. add_entrapment_sequences  ← Sets modifications to `missing`
4. assign_base_entrap_ids
5. add_mods                  ← Adds modifications after entrapment
6. assign_base_pep_ids
7. add_decoy_sequences       ← Properly handles modifications with adjust_mod_positions
8. add_charge
9. assign_base_prec_ids
```

### Issues
1. **Entrapment sequences lose modifications**: `add_entrapment_sequences` sets both structural_mods and isotopic_mods to `missing`
2. **Inconsistent modification handling**: `add_decoy_sequences` properly adjusts modification positions after shuffling, but `add_entrapment_sequences` doesn't
3. **Logical flow problem**: Modifications are added AFTER entrapment, so entrapment sequences never get the benefit of shuffled modifications

## Solution Overview

Move `add_mods` to BEFORE `add_entrapment_sequences` and update `add_entrapment_sequences` to handle modifications properly using the same logic as `add_decoy_sequences`.

### New Pipeline Order
```
1. combine_shared_peptides
2. assign_base_seq_ids
3. add_mods                  ← MOVED: Now before entrapment
4. assign_base_pep_ids       ← MOVED: After mods, before entrapment
5. add_entrapment_sequences  ← UPDATED: Will handle modification shuffling
6. assign_base_entrap_ids    ← After entrapment sequences
7. add_decoy_sequences
8. add_charge
9. assign_base_prec_ids
```

## Detailed Implementation Plan

### Phase 1: Update Pipeline Order in chronologer_prep.jl

#### File: `/src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**Current Steps to Move:**
- Move `add_mods` call from Step 4 to after Step 2 (assign_base_seq_ids)
- Move `assign_base_pep_ids` call from Step 5 to after add_mods
- Keep `add_entrapment_sequences` but update it to handle modifications
- Move `assign_base_entrap_ids` to after entrapment sequences

**Specific Changes:**
1. **Lines ~194-200**: Move `add_mods` call and logging to after `assign_base_seq_ids`
2. **Lines ~202-207**: Move `assign_base_pep_ids` call to after `add_mods`
3. **Lines ~183-192**: Update `add_entrapment_sequences` call position and logging
4. **Update step numbering** in all log messages

### Phase 2: Update add_entrapment_sequences Function

#### File: `/src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

**Current Issues in add_entrapment_sequences:**
- Lines 180-181: `missing, #structural_mods` and `missing, #isotopic_mods`
- No position tracking for modification adjustment

**Required Changes:**

1. **Add modification adjustment logic** (similar to add_decoy_sequences lines 518-529):
   ```julia
   # Adjust modification positions based on sequence shuffling
   adjusted_structural_mods = adjust_mod_positions(
       get_structural_mods(target_entry),
       shuffle_seq.new_positions,
       seq_length
   )
   
   adjusted_isotopic_mods = adjust_mod_positions(
       get_isotopic_mods(target_entry),
       shuffle_seq.new_positions,
       seq_length
   )
   ```

2. **Update FastaEntry constructor** to use adjusted modifications:
   ```julia
   entrapment_fasta_entries[n] = FastaEntry(
       # ... other fields ...
       adjusted_structural_mods,  # Instead of missing
       adjusted_isotopic_mods,    # Instead of missing
       # ... remaining fields ...
   )
   ```

3. **Add sequence length tracking**:
   ```julia
   seq_length = UInt8(length(get_sequence(target_entry)))
   ```

### Phase 3: Update ID Assignment Logic

Since we're changing the pipeline order, we need to ensure proper ID assignment:

1. **base_seq_id**: Assigned after combine_shared_peptides (unchanged)
2. **base_pep_id**: Now assigned after add_mods (includes modification information)
3. **base_entrap_id**: Assigned after add_entrapment_sequences (includes entrapment grouping)
4. **base_prec_id**: Assigned after add_charge (unchanged)

### Phase 4: Handle base_prec_id in Entrapment

**Current Issue**: `add_entrapment_sequences` increments base_prec_id starting from maximum existing value.

**Solution**: Since modifications are now added before entrapment, we need to ensure base_prec_id values don't conflict:
- Keep the current logic but ensure it works with modified peptides
- The increment logic (lines 145-151) should still work correctly

### Phase 5: Update Function Documentation

Update function docstrings to reflect new pipeline order:

1. **add_entrapment_sequences**: Update to mention it handles modifications
2. **add_mods**: Update to mention it's called before entrapment
3. **Pipeline documentation**: Update any comments about order

## Implementation Steps

### Step 1: Backup and Preparation
- Current code is already committed
- Create this plan document
- Test current pipeline to establish baseline

### Step 2: Update chronologer_prep.jl Pipeline Order
- Move add_mods to before add_entrapment_sequences
- Move assign_base_pep_ids to after add_mods
- Update step numbering and logging
- Test compilation

### Step 3: Update add_entrapment_sequences Function
- Add modification adjustment logic
- Update FastaEntry constructor calls
- Add sequence length tracking
- Test with simple modifications

### Step 4: Integration Testing
- Test full pipeline with keap1 test data
- Verify modification pairing works correctly
- Check that entrapment sequences have proper modifications

### Step 5: Validation
- Run modification pairing check script
- Ensure target-decoy pairs still match correctly
- Verify all ID assignment works properly

## Expected Outcomes

1. **Entrapment sequences will have shuffled modifications**: Just like decoy sequences
2. **Consistent modification handling**: Both entrapment and decoy sequences use same adjustment logic
3. **Proper ID hierarchy**: base_seq_id → base_pep_id (with mods) → base_entrap_id (with entrapment) → base_prec_id (with charge)
4. **Maintained target-decoy pairing**: Using (base_pep_id, precursor_charge) as updated

## Risk Assessment

**Low Risk Changes:**
- Moving pipeline steps (well-defined interfaces)
- Reusing existing modification adjustment logic

**Medium Risk Changes:**
- Updating add_entrapment_sequences (may affect base_prec_id assignment)
- ID assignment timing changes

**Mitigation:**
- Test with small datasets first
- Keep existing test data for comparison
- Commit each phase separately for easy rollback

## Files to Modify

1. **`/src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`** - Pipeline order
2. **`/src/Routines/BuildSpecLib/fasta/fasta_utils.jl`** - add_entrapment_sequences function
3. **This plan document** - For tracking progress

## Success Criteria

- [ ] Pipeline compiles without errors
- [ ] Entrapment sequences have non-missing modifications
- [ ] Target-decoy pairing still works (42,092 pairs expected)
- [ ] All ID fields populated correctly in output DataFrame
- [ ] Library builds successfully to completion