# Base Pep ID Duplication Fix Plan

## Current Problem Analysis

### Issue Discovered
Multiple completely different sequences are sharing the same `base_pep_id`, violating the fundamental design principle where each unique original peptide sequence should have its own identifier.

**Example from base_pep_id = 1270:**
- `HSPGGGHGGRSQGSQR` 
- `SGPSHRHQQGSGGGGR`
- `GGSHPGHRGGQGSQSR`
- `GRHGSGSGQSPGHGQR`  
- `IYEIAYR`
- `YIYEIAR`
- `AIYEIYR`
- `YAIYEIR`

These 8 completely different sequences (not I/L equivalent) all share `base_pep_id = 1270`.

### Root Cause Analysis

**Pipeline Flow:**
1. **FASTA Digestion** (`fasta_digest.jl`) → Each peptide gets unique sequential `base_pep_id`
2. **Combine Shared Peptides** (`fasta_utils.jl`) → I/L equivalent sequences merged (**BUG LOCATION**)
3. **Add Entrapments** → Inherit `base_pep_id` from targets ✅
4. **Add Decoys** → Inherit `base_pep_id` from targets ✅
5. **Add Modifications** → Preserve `base_pep_id` ✅
6. **Add Charges** → Preserve `base_pep_id` ✅

**Bug in `combine_shared_peptides` Function:**
The function correctly identifies I/L equivalent sequences but has logic errors in `base_pep_id` assignment that cause non-equivalent sequences to share the same ID.

### Impact Assessment
- ❌ **Entrapment tracking broken**: Different peptides grouped incorrectly
- ❌ **Data integrity compromised**: Fundamental ID system violated
- ❌ **Analysis validity questioned**: Statistics based on base_pep_id counts are wrong
- ✅ **High pair_id fix still works**: This is a separate issue from the recent fix

## Proposed Solution

### Strategy: Post-Processing Base Pep ID Reassignment

**Key Insight**: Instead of fixing the complex `combine_shared_peptides` logic, we'll implement a **clean-up step** that reassigns `base_pep_id` values after shared peptides are combined.

### Implementation Plan

#### Step 1: Add Reassignment Function
**File**: `src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

**New Function**: `reassign_base_pep_ids!(fasta_entries::Vector{FastaEntry})`

```julia
function reassign_base_pep_ids!(fasta_entries::Vector{FastaEntry})
    """
    Reassign base_pep_id values to ensure each unique sequence has its own ID.
    This fixes any duplication issues from combine_shared_peptides.
    
    Algorithm:
    1. Group entries by I/L-equivalent sequence
    2. Assign same base_pep_id to all entries in each group
    3. Use sequential IDs across groups
    """
    
    # Create sequence -> base_pep_id mapping
    seq_to_base_pep_id = Dict{String, UInt32}()
    current_base_pep_id = one(UInt32)
    
    for entry in fasta_entries
        sequence = get_sequence(entry)
        il_equivalent = replace(sequence, 'I' => 'L')
        
        if !haskey(seq_to_base_pep_id, il_equivalent)
            seq_to_base_pep_id[il_equivalent] = current_base_pep_id
            current_base_pep_id += one(UInt32)
        end
    end
    
    # Update all entries with correct base_pep_id
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        sequence = get_sequence(entry)
        il_equivalent = replace(sequence, 'I' => 'L')
        correct_base_pep_id = seq_to_base_pep_id[il_equivalent]
        
        # Create new FastaEntry with corrected base_pep_id
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
            correct_base_pep_id,  # ← Fixed base_pep_id
            get_base_prec_id(entry),
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(seq_to_base_pep_id)  # Return number of unique sequences
end
```

#### Step 2: Integrate Into Pipeline
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**Location**: After `combine_shared_peptides` call (line 174)

```julia
# Combine shared peptides
fasta_entries = combine_shared_peptides(fasta_entries)

# FIX: Reassign base_pep_id values to ensure uniqueness
unique_sequences = reassign_base_pep_ids!(fasta_entries)
@user_info "Reassigned base_pep_id values: $unique_sequences unique sequences"

# Add entrapment sequences if specified
fasta_entries = add_entrapment_sequences(
    fasta_entries,
    UInt8(_params.fasta_digest_params["entrapment_r"])
)
```

#### Step 3: Add Diagnostic Logging
Add logging to verify the fix:

```julia
# Before reassignment
before_count = length(unique([get_base_pep_id(e) for e in fasta_entries]))
before_sequences = length(unique([get_sequence(e) for e in fasta_entries]))

# After reassignment
unique_sequences = reassign_base_pep_ids!(fasta_entries)
after_count = length(unique([get_base_pep_id(e) for e in fasta_entries]))

@user_info "Base Pep ID Reassignment Results:"
@user_info "  Before: $before_count unique base_pep_ids for $before_sequences sequences"
@user_info "  After:  $after_count unique base_pep_ids for $unique_sequences sequences"
@user_info "  Status: $(after_count == unique_sequences ? "✅ FIXED" : "❌ STILL BROKEN")"
```

## Expected Results

### Before Fix
- 8 different sequences sharing `base_pep_id = 1270`
- Unpredictable duplication patterns throughout the library
- Entrapment tracking compromised

### After Fix
- Each unique sequence gets its own `base_pep_id`
- I/L equivalent sequences share the same `base_pep_id` (intended behavior)
- Clean 1:1 mapping between unique sequences and base_pep_ids
- Entrapment tracking restored

### Validation Metrics
1. **Uniqueness Check**: `length(unique(base_pep_ids)) == length(unique(il_equivalent_sequences))`
2. **Entrapment Preservation**: Same entrapment coverage as before
3. **Library Size**: Same total precursor count
4. **I/L Equivalence**: Only I/L equivalent sequences share base_pep_ids

## Testing Strategy

### Phase 1: Unit Test the Function
Create test with known sequences:
```julia
test_entries = [
    entry_with_sequence("PEPTIDE"),    # Should get base_pep_id = 1
    entry_with_sequence("PEPTLDE"),    # Should get base_pep_id = 1 (I/L equivalent)
    entry_with_sequence("DIFFERENT"),  # Should get base_pep_id = 2
]
reassign_base_pep_ids!(test_entries)
# Verify correct assignments
```

### Phase 2: Integration Test
1. Rebuild keap1 library with fix
2. Verify `base_pep_id = 1270` issue resolved
3. Check random sampling of other base_pep_ids
4. Ensure entrapment coverage maintained

### Phase 3: Validation Suite
Run existing `TEST_KEAP1_LIBRARY.jl` and verify:
- No sequences share base_pep_id unless I/L equivalent
- Entrapment-target coverage remains 100%
- All other functionality preserved

## Implementation Notes

### Why This Approach?
1. **Clean Separation**: Fixes the issue without modifying complex existing logic
2. **Minimal Risk**: Only touches base_pep_id assignment, preserves all other data
3. **Clear Intent**: Makes the fix obvious and maintainable
4. **Defensive**: Protects against future similar bugs

### Alternative Approaches Considered
1. **Fix `combine_shared_peptides` Logic**: More complex, higher risk of breaking shared peptide functionality
2. **Remove `combine_shared_peptides`**: Would break I/L equivalence handling
3. **Fix Earlier in Pipeline**: Would require changes to multiple functions

### Potential Edge Cases
1. **Empty Sequences**: Handle gracefully (shouldn't occur in practice)
2. **Very Large Libraries**: Memory usage should be minimal (just ID reassignment)
3. **Special Characters**: I/L replacement should handle standard amino acids

### Backward Compatibility
- No changes to file formats or APIs
- No changes to downstream processing
- Only internal base_pep_id values change (user never sees these directly)

## Files to Modify

### Primary Changes
- `src/Routines/BuildSpecLib/fasta/fasta_utils.jl` - Add reassignment function
- `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl` - Integrate into pipeline

### Testing
- Update `TEST_KEAP1_LIBRARY.jl` - Add specific base_pep_id uniqueness checks
- Create unit test for `reassign_base_pep_ids!` function

## Success Criteria

1. ✅ **Each unique sequence has unique base_pep_id**
2. ✅ **I/L equivalent sequences share base_pep_id**  
3. ✅ **No regression in entrapment tracking**
4. ✅ **Library builds successfully**
5. ✅ **All validation tests pass**
6. ✅ **Performance impact minimal (<5% build time increase)**

This plan provides a surgical fix that addresses the root cause while minimizing risk to the complex existing pipeline logic.