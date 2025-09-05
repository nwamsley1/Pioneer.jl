# Base Prec ID Pairing Plan

## Current Problem Analysis

### The Issue
We have two different ID systems that are conflated:

1. **base_pep_id**: Should track **entrapment relationships** - shared between original targets and their entrapment sequences
2. **base_prec_id**: Should track **modification-specific pairing** - each modification variant needs unique ID for proper target-decoy pairing

Currently, we're using base_pep_id for both purposes, which causes:
- ❌ **Modification pairing issues**: Modified targets pair with unmodified decoys
- ❌ **Complex logic**: One ID trying to serve two different purposes

### Root Cause
The pairing function `add_charge_specific_partner_columns!` uses `(base_pep_id, charge)` for pairing, but base_pep_id is shared across modification variants. We need it to use `(base_prec_id, charge)` where base_prec_id is unique per modification variant.

## Proposed Solution: Dual ID System

### Core Concept
- **base_pep_id**: Entrapment tracking (shared across modification variants)
- **base_prec_id**: Modification-specific pairing (unique per modification variant)

### New Pipeline Order
```
1. combine_shared_peptides()           # Merge I/L equivalent sequences
2. reassign_base_pep_ids!()           # Clean base_pep_id assignment  
3. add_entrapment_sequences()         # Inherit base_pep_id for tracking
4. add_mods()                         # Create modification variants
5. reassign_base_prec_ids!()          # Unique base_prec_id per variant
6. add_decoys()                       # Inherit base_prec_id for pairing  
7. add_charge()                       # Preserve base_prec_id per charge
8. add_charge_specific_partner_columns!() # Pair by (base_prec_id, charge)
```

## Detailed Implementation Plan

### Step 1: Create `reassign_base_prec_ids!()` Function
**File**: `src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

**New Function**:
```julia
function reassign_base_prec_ids!(fasta_entries::Vector{FastaEntry})
    """
    Reassign sequential base_prec_id values starting from 1.
    Should be called after add_mods() to ensure each modification variant 
    gets unique base_prec_id for proper target-decoy pairing.
    
    Returns:
    - Int: Number of entries processed
    """
    
    # Simply assign sequential base_prec_id starting from 1
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        
        # Create new FastaEntry with sequential base_prec_id
        fasta_entries[i] = FastaEntry(
            get_id(entry),
            get_description(entry),
            get_gene(entry),
            get_protein(entry),
            get_organism(entry),
            get_proteome(entry),
            get_sequence(entry),
            get_start_idx(entry),
            get_structural_mods(entry),
            get_isotopic_mods(entry),
            get_charge(entry),
            get_base_pep_id(entry),  # ← Preserve base_pep_id for entrapment tracking
            UInt32(i),               # ← Sequential base_prec_id for modification pairing
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(fasta_entries)
end
```

### Step 2: Update Pipeline Order
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**Current Order** (incorrect):
```julia
# Combine shared peptides
fasta_entries = combine_shared_peptides(fasta_entries)
reassign_base_pep_ids!(fasta_entries)  # After combine_shared_peptides

# Add entrapment sequences 
fasta_entries = add_entrapment_sequences(fasta_entries, ...)

# Add modifications
fasta_entries = add_mods(fasta_entries, ...)

# PROBLEM: Using base_pep_id reassignment for modification variants 
reassign_base_pep_ids!(fasta_entries)  # Wrong timing and purpose

# Add decoys
fasta_entries = add_decoy_sequences(fasta_entries)

# Add charges
fasta_entries = add_charge(fasta_entries, ...)
```

**New Order** (correct):
```julia
# Combine shared peptides
fasta_entries = combine_shared_peptides(fasta_entries)

# FIX 1: Clean base_pep_id assignment for entrapment tracking
before_count = length(unique([get_base_pep_id(e) for e in fasta_entries]))
entries_processed = reassign_base_pep_ids!(fasta_entries)
after_count = length(unique([get_base_pep_id(e) for e in fasta_entries]))
@user_info "Base Pep ID Reassignment: $entries_processed entries processed"
@user_info "  Before: $before_count → After: $after_count unique base_pep_ids" 
@user_info "  Purpose: Clean base_pep_id assignment for entrapment tracking"

# Add entrapment sequences (inherit base_pep_id for tracking)
fasta_entries = add_entrapment_sequences(fasta_entries, ...)

# Add modifications (creates variants sharing base_pep_id)
fasta_entries = add_mods(fasta_entries, ...)

# FIX 2: Unique base_prec_id assignment for modification-specific pairing
before_prec_count = length(unique([get_base_prec_id(e) for e in fasta_entries]))
prec_entries_processed = reassign_base_prec_ids!(fasta_entries)
after_prec_count = length(unique([get_base_prec_id(e) for e in fasta_entries]))
@user_info "Base Prec ID Reassignment: $prec_entries_processed entries processed"
@user_info "  Before: $before_prec_count → After: $after_prec_count unique base_prec_ids"
@user_info "  Purpose: Unique base_prec_id per modification variant for proper pairing"

# Add decoys (inherit base_prec_id for pairing)
fasta_entries = add_decoy_sequences(fasta_entries)

# Add charges (preserve base_prec_id per charge - already fixed)
fasta_entries = add_charge(fasta_entries, ...)
```

### Step 3: Update Pairing Function (Already Correct!)
**File**: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`

The `add_charge_specific_partner_columns!()` function **already pairs by base_pep_id**, but we need to verify it should use base_prec_id instead:

**Current pairing key**:
```julia
key = (row.base_pep_id, row.precursor_charge, !row.decoy)
```

**Should remain the same** because:
- After our fixes, base_pep_id will actually contain the base_prec_id values due to our reassignment
- OR we need to update it to use base_prec_id explicitly

**Analysis needed**: Check if the pairing function needs to be updated to use base_prec_id instead of base_pep_id.

### Step 4: Update add_charge() Function (Already Fixed!)
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

The add_charge() function is already fixed to preserve base_prec_id:
```julia
get_base_prec_id(fasta_peptide),  # Preserve from modification variant
```

No additional changes needed here.

## Expected Results

### Before Fix
```
Pipeline: combine_shared_peptides → reassign_base_pep_ids → entrapments → mods → reassign_base_pep_ids → decoys → charges → pairing

Issues:
- base_pep_id used for both entrapment tracking AND modification pairing
- Modified targets pair with unmodified decoys  
- Complex ID management with dual purposes
```

### After Fix  
```
Pipeline: combine_shared_peptides → reassign_base_pep_ids → entrapments → mods → reassign_base_prec_ids → decoys → charges → pairing

Benefits:
- base_pep_id: Clean entrapment tracking (shared across variants)
- base_prec_id: Clean modification pairing (unique per variant)
- Modified targets pair with modified decoys correctly
- Clear separation of concerns
```

### Validation Metrics
1. **Entrapment Coverage**: Should remain 100% (base_pep_id preserved)
2. **Modification Pairing**: 0% mismatched modifications (base_prec_id ensures matching)
3. **Pairing Success**: 100% targets paired with decoys
4. **ID Uniqueness**: 
   - base_pep_id: One per base peptide sequence (shared across variants)
   - base_prec_id: One per modification variant (unique pairing)

## Implementation Steps

### Phase 1: Create reassign_base_prec_ids!() Function ✅
- Add function to fasta_utils.jl
- Mirror reassign_base_pep_ids!() but for base_prec_id field
- Include diagnostic logging

### Phase 2: Update Pipeline Order
- Reorder function calls in chronologer_prep.jl
- Add diagnostic logging to track ID assignments
- Update comments to clarify purpose of each step

### Phase 3: Verify Pairing Function
- Check if add_charge_specific_partner_columns needs to use base_prec_id
- Update pairing key if necessary
- Test pairing logic with new ID system

### Phase 4: Test and Validate
- Rebuild keap1 library with new pipeline
- Verify entrapment tracking preserved
- Confirm modification pairing works correctly  
- Run comprehensive validation suite

## Key Design Principles

### Separation of Concerns
- **base_pep_id**: Entrapment relationships only
- **base_prec_id**: Modification pairing only
- Clear, single-purpose ID system

### Timing is Critical
- base_pep_id assigned after combine_shared_peptides, before entrapments
- base_prec_id assigned after add_mods, before decoys
- Each ID serves its purpose at the right pipeline stage

### Preservation Strategy
- Entrapments inherit base_pep_id (for tracking)
- Decoys inherit base_prec_id (for pairing)  
- Charges preserve base_prec_id (for pairing)

## Risk Assessment

### Low Risk
- reassign_base_prec_ids!() is identical to working reassign_base_pep_ids!()
- Pipeline reordering uses same functions in different sequence
- Clear separation improves maintainability

### Medium Risk  
- Need to verify pairing function uses correct ID field
- Extensive testing required to ensure no regressions
- Must preserve all existing functionality

### Success Criteria
1. ✅ **Entrapment tracking preserved**: 100% coverage maintained
2. ✅ **Modification pairing perfect**: 0% mismatched modifications
3. ✅ **All targets paired**: No unpaired decoy warnings
4. ✅ **Clean ID system**: Single purpose per ID type
5. ✅ **No regressions**: All validation tests pass

This plan provides a clear path to proper modification-specific pairing while preserving entrapment tracking functionality.