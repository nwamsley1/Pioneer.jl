# Fix Plan: Preserve base_pep_id for Entrapment Sequence Tracking

## Problem Statement  

Currently, `base_pep_id` is being reassigned during modification generation, breaking the link between:
- Original target peptides
- Their entrapment sequences  
- Their decoy sequences
- All modification variants

This prevents proper tracking of which entrapment sequences belong to which original targets.

## VERIFIED: base_pep_id Flow Through Pipeline

After thorough analysis, I've confirmed:

1. **base_pep_id IS preserved through the entire pipeline** - it flows from initial creation all the way to the final `precursors_table.arrow`
2. **The ONLY place it's incorrectly modified** is in `add_mods()` at lines 320-334 of `chronologer_prep.jl`
3. **Current result**: 29,024 unique base_pep_ids in the final table (should be much fewer - only original peptides)
4. **No other transformations modify base_pep_id**:
   - ✅ Preserved through chronologer processing
   - ✅ Preserved through sorting
   - ✅ Preserved through filtering
   - ✅ Present in final precursors_table.arrow

## Current Problematic Flow

```
1. Digestion → base_pep_id assigned (1,2,3...)
2. Add entrapment → entrapment inherits target's base_pep_id ✓
3. Add modifications → NEW base_pep_id assigned (breaks link!) ✗
4. Add decoys → decoys inherit modified base_pep_id
5. Add charges → charges inherit modified base_pep_id
```

## Root Cause

In `chronologer_prep.jl` lines 320-334:
```julia
current_base_pep_id += 1  # Creates NEW ID, breaking entrapment link
...
current_base_pep_id,  # Overwrites original base_pep_id
```

## Desired Behavior

`base_pep_id` should:
1. Be assigned once during initial peptide digestion
2. Be preserved through ALL transformations:
   - Entrapment sequences (same base_pep_id)
   - Modification variants (same base_pep_id)
   - Decoy sequences (same base_pep_id)
   - Charge states (same base_pep_id)
3. Link all variants back to the original target peptide

## Solution Design

### Option 1: Preserve Original base_pep_id (RECOMMENDED)

**Change in `add_mods()` function:**
```julia
# Line 334: Use original base_pep_id instead of new one
get_base_pep_id(fasta_peptide),  # Preserve original ID
```

**Impact:**
- All modification variants share same base_pep_id
- Entrapment sequences remain linked to targets
- Decoys still paired correctly via (base_pep_id, charge, mods)

**New Pairing Logic Required:**
- Modify `add_charge_specific_partner_columns!` to use:
  - `(base_pep_id, charge, modification_string)` for unique pairing
  - OR create a new `mod_variant_id` field for pairing

### Option 2: Add Separate Tracking Field

**Add new field `original_base_pep_id`:**
```julia
struct FastaEntry
    ...
    base_pep_id::UInt32          # Used for pairing
    original_base_pep_id::UInt32 # Links to original target
    ...
end
```

**Impact:**
- More complex but preserves current pairing logic
- Requires changes to FastaEntry struct
- Clear separation of concerns

### Option 3: Use entrapment_group_id for Linking

**Current state:**
- `entrapment_group_id` already tracks entrapment groups (1,2,3...)
- Could extend to track original peptide ID

**Proposed change:**
```julia
# In add_entrapment_sequences:
entrapment_group_id,  # Current: 1,2,3...
# Change to:
(get_base_pep_id(target_entry) << 8) | entrapment_group_id  # Encode both
```

**Impact:**
- No struct changes needed
- Preserves current base_pep_id logic
- Can decode original peptide from entrapment_group_id

## Pipeline Data Flow Verification

### Locations where base_pep_id is created/modified:
1. **fasta_digest.jl:166-200** - Initial creation (starts at 1, increments for each peptide)
2. **fasta_utils.jl:600-631** - Alternative creation path (similar increment logic)
3. **chronologer_prep.jl:320-334** - ❌ **PROBLEM: Reassigns new values for modifications**

### Locations where base_pep_id is preserved:
1. **add_entrapment_sequences()** - Line 183: `get_base_pep_id(target_entry)`
2. **add_decoy_sequences()** - Line 542: `get_base_pep_id(target_entry)`  
3. **add_charge()** - Line 388: `get_base_pep_id(fasta_peptide)`
4. **build_fasta_df()** - Line 567: `get_base_pep_id(peptide)`
5. **chronologer processing** - Passed through unchanged
6. **Final Arrow output** - Column preserved in precursors_table.arrow

## Simplified Implementation (ONLY Preserve base_pep_id)

### Single Goal: 
**Ensure each original target sequence has a unique base_pep_id that persists through the entire BuildSpecLib process and appears unchanged in the final precursors_table.arrow**

### Step 1: Fix add_mods() to Preserve base_pep_id

**File:** `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Lines:** 289-340

**Changes:**
```julia
# Remove line 292:
# current_base_pep_id = UInt32(0)

# Remove line 320:
# current_base_pep_id += 1

# Change line 334:
get_base_pep_id(fasta_peptide),  # Preserve original base_pep_id
```

**That's it!** No pairing logic changes needed.

### Expected Result:
- Original targets get base_pep_id during initial digestion (1, 2, 3...)
- Entrapment sequences inherit same base_pep_id as their target
- Decoy sequences inherit same base_pep_id as their target
- All modification variants keep same base_pep_id as their parent
- All charge states keep same base_pep_id as their parent
- Final table has ~10,000 unique base_pep_ids instead of 29,024

### Verification:

**Add validation to ensure:**
1. Original targets and entrapment sequences share base_pep_id
2. Each modification variant can still find its decoy partner
3. SearchDIA correctly groups entrapment sequences

## Testing Strategy

### Test 1: Verify base_pep_id Preservation
```julia
# After library build, check:
df = DataFrame(Arrow.Table("precursors_table.arrow"))
grouped = groupby(df, :base_pep_id)
for group in grouped
    # Should contain:
    # - Original target(s)
    # - Entrapment sequences (if entrapment_r > 0)
    # - Decoy sequences
    # - All modification variants
    @assert all(group.base_pep_id .== first(group.base_pep_id))
end
```

### Test 2: Verify Correct Pairing
```julia
# Check that each target has exactly one decoy partner
for group in groupby(df, [:base_pep_id, :prec_charge, :structural_mods])
    targets = sum(.!group.is_decoy)
    decoys = sum(group.is_decoy)
    @assert targets == decoys == 1
end
```

### Test 3: Verify Entrapment Grouping
```julia
# Check entrapment sequences share base_pep_id with target
targets = df[.!df.is_decoy .& (df.entrapment_group_id .== 0), :]
for target in eachrow(targets)
    entrapments = df[(df.base_pep_id .== target.base_pep_id) .& 
                     (df.entrapment_group_id .> 0), :]
    @assert nrow(entrapments) == entrapment_r
end
```

## Risk Assessment

### Low Risk
- Preserving base_pep_id is conceptually cleaner
- Entrapment tracking becomes more intuitive
- No struct changes needed

### Medium Risk  
- Pairing logic needs modification string comparison
- Performance impact of string operations in pairing
- Need to ensure modification string format is consistent

### Mitigation
- Cache modification strings during pairing
- Use hash of modifications instead of string comparison
- Comprehensive testing with multiple modification scenarios

## Alternative: Minimal Change Approach

If changing pairing logic is too risky, use Option 3:
- Encode original base_pep_id in upper bits of entrapment_group_id
- No changes to pairing logic needed
- Can extract original peptide ID when needed for analysis

## Rollback Plan

If issues arise:
1. Revert add_mods() to use incremental base_pep_id
2. Document that base_pep_id is NOT stable across modifications
3. Use alternative tracking method for entrapment analysis

## Success Criteria

1. ✅ Entrapment sequences share base_pep_id with original targets
2. ✅ All modification variants share same base_pep_id
3. ✅ Target-decoy pairing still works correctly
4. ✅ SearchDIA runs without errors
5. ✅ Entrapment analysis can group sequences correctly

## Implementation Priority

1. **High Priority:** Fix base_pep_id preservation (Option 1, Step 1)
2. **High Priority:** Update pairing logic (Option 1, Step 2)
3. **Medium Priority:** Add comprehensive tests
4. **Low Priority:** Performance optimizations if needed