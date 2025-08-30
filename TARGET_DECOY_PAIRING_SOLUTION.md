# Target-Decoy Pairing Solution for Pioneer.jl BuildSpecLib

## Executive Summary

The target-decoy pairing in Pioneer.jl is broken because each precursor gets a unique `pair_id` instead of sharing IDs with their target/decoy partner. The solution involves using the commented-out `add_precursor_partner_columns!` function from `pair_decoys.jl` which properly creates charge-specific target-decoy pairs.

## The Multi-Layered Problem

### Layer 1: Current Bug
**File**: `chronologer_prep.jl`, line 554
```julia
_pair_id[n] = get_base_prec_id(peptide)  # Each precursor gets unique ID
```

### Layer 2: Naive Fix Would Create New Problem
Simply using `base_pep_id` would group ALL charge states together:
```
Peptide A, charge 2, target   → pair_id = 10
Peptide A, charge 2, decoy    → pair_id = 10
Peptide A, charge 3, target   → pair_id = 10  (WRONG!)
Peptide A, charge 3, decoy    → pair_id = 10  (WRONG!)
```
This creates groups of 4+ precursors per pair_id instead of proper target-decoy pairs.

### Layer 3: Correct Requirement
Each target-decoy pair should be **charge-specific**:
```
Peptide A, charge 2, target   → pair_id = 10
Peptide A, charge 2, decoy    → pair_id = 10
Peptide A, charge 3, target   → pair_id = 11  (Different pair!)
Peptide A, charge 3, decoy    → pair_id = 11
```

## The Complete Solution

### Option 1: Use Existing Function (Recommended)

#### Step 1: Uncomment the Pairing Function
**File**: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`
**Lines**: 97-177

The `add_precursor_partner_columns!` function is already implemented but commented out. It:
1. Takes the complete precursor DataFrame after all processing
2. Creates lookup dictionary for all precursors
3. For each target, finds its decoy by sequence reversal
4. Matches based on (proteome, sequence, mods, **charge**)
5. Assigns shared `pair_id` to matched pairs
6. Sets `partner_precursor_idx` for direct partner lookup

#### How the Matching Works
```julia
# For each target precursor:
key = (row.proteome_identifiers, row.sequence, mods, row.prec_charge)
lookup[key] = row_index

# Find its decoy partner:
partner_seq, partner_mods = reverseSequence(row.sequence, row.structural_mods)
partner_key = (row.proteome_identifiers, partner_seq, partner_mods, row.prec_charge)

if haskey(lookup, partner_key)
    partner_idx = lookup[partner_key]
    # Both get same pair_id
    pair_ids[idx] = next_pair_id
    pair_ids[partner_idx] = next_pair_id
    
    # Both get partner references
    partner_indices[idx] = partner_idx
    partner_indices[partner_idx] = idx
end
```

#### Step 2: Modify chronologer_prep.jl
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**Change line 554** from:
```julia
_pair_id[n] = get_base_prec_id(peptide)
```
**To**:
```julia
_pair_id[n] = zero(UInt32)  # Temporary placeholder
```

**Add after DataFrame creation** (around line 571):
```julia
# Create initial DataFrame
seq_df = DataFrame(
    (proteome_identifiers = _proteome_identifiers[1:n],
     sequence = _sequence[1:n],
     # ... other columns ...
     pair_id = _pair_id[1:n])
)

# Apply proper target-decoy pairing
seq_df_with_pairs, lookup = add_precursor_partner_columns!(seq_df)

# Replace the original DataFrame
seq_df = seq_df_with_pairs
```

#### Step 3: Remove Conflicting Function
**Remove or comment out the call to `add_pair_indices!`** since it assumes each `pair_id` appears exactly twice, which won't work with unique pair_ids.

### Option 2: Custom Implementation (Alternative)

If uncommenting the existing function is too risky, implement a charge-specific pairing system:

```julia
function create_charge_specific_pair_ids(fasta_peptides::Vector{FastaEntry})
    pair_counter = Dict{Tuple{UInt32, UInt8}, UInt32}()  # (base_pep_id, charge) -> pair_id
    pair_ids = Vector{UInt32}(undef, length(fasta_peptides))
    next_pair_id = UInt32(1)
    
    for (n, peptide) in enumerate(fasta_peptides)
        pep_id = get_base_pep_id(peptide)
        charge = get_charge(peptide)
        
        # Key ignores target/decoy status - both share same pair_id
        key = (pep_id, charge)
        
        if !haskey(pair_counter, key)
            pair_counter[key] = next_pair_id
            next_pair_id += 1
        end
        
        pair_ids[n] = pair_counter[key]
    end
    
    return pair_ids
end

# Use in prepare_chronologer_input:
_pair_id = create_charge_specific_pair_ids(fasta_peptides)
```

## Implementation Details

### Sequence Reversal Logic
The `reverseSequence` function in `pair_decoys.jl` properly handles:
- **Sequence reversal**: `PEPTIDE` → `EDITPEP` (preserves C-terminal residue)
- **Modification mapping**: Adjusts modification positions for reversed sequence
- **Missing modifications**: Handles cases with no modifications

### Lookup Key Structure
The pairing uses a comprehensive key:
```julia
key = (
    row.proteome_identifiers,  # Protein accession(s)
    row.sequence,              # Peptide sequence
    mods,                      # Set of modifications
    row.prec_charge            # Charge state
)
```

This ensures:
- Same peptide at different charges = different pairs
- Different modifications = different pairs
- Different proteins = different pairs

### Partner Index Benefits
Setting `partner_precursor_idx` provides:
1. **O(1) partner lookup**: Direct array access instead of search
2. **MBR efficiency**: Instant partner identification
3. **Validation**: Easy to verify pairing integrity

## Expected Results

### Before Fix
```julia
julia> gptable = groupby(ptable, :pair_id)
julia> unique([size(subdf,1) for (key, subdf) in pairs(gptable)])
1-element Vector{Int64}: 1  # Each pair_id has only 1 member
```

### After Fix
```julia
julia> gptable = groupby(ptable, :pair_id)
julia> unique([size(subdf,1) for (key, subdf) in pairs(gptable)])
1-element Vector{Int64}: 2  # Each pair_id has exactly 2 members
```

## Validation Testing

The `pair_decoys.jl` file includes validation tests (lines 183-191):
```julia
# Test proper pairing
value_counts(df, col) = combine(groupby(df, col), nrow)
value_counts(ptable, :pair_id)
pair_id_counts = unique(value_counts(ptable, :pair_id)[!,:nrow])
@test length(pair_id_counts) == 1     # Only one group size
@test pair_id_counts[1] == 2          # Each group has exactly 2 members
```

## Impact on Downstream Systems

### MBR (Match-Between-Runs)
**Before**: Cannot find partners, fails to add decoys
```julia
partner_pid = getPartnerPrecursorIdx(precursors)[pid]  # Returns missing
```

**After**: Direct partner lookup works
```julia
partner_pid = getPartnerPrecursorIdx(precursors)[pid]  # Returns valid index
if !haskey(precursor_dict, partner_pid)
    insert!(precursor_dict, partner_pid, val)  # Successfully adds partner
end
```

### CV Fold Assignment
**Before**: Cannot ensure pairs stay together
```julia
# Each precursor has unique pair_id, so no pairing constraints
```

**After**: Can implement pair-aware assignment
```julia
# Assign pairs to same CV fold based on shared pair_id
for pair_id in unique_pair_ids
    fold = rand(cv_folds)
    # All precursors with this pair_id go to same fold
end
```

### Target/Decoy Ratios
**Before**: Imbalanced ratios (e.g., 4.5:1)
```
Targets: 133,252
Decoys: 29,408
```

**After**: Balanced ratios (close to 1:1 with MBR)
```
Targets: ~81,000
Decoys: ~81,000
```

## Migration Strategy

### Phase 1: Immediate Fix
1. Uncomment `add_precursor_partner_columns!` function
2. Modify `chronologer_prep.jl` to use placeholder pair_ids
3. Apply proper pairing after DataFrame creation

### Phase 2: Validation
1. Run test builds with new pairing system
2. Validate pair_id distributions
3. Test MBR functionality
4. Verify CV fold assignment works

### Phase 3: Integration
1. Update CV fold assignment to use proper pairs
2. Enhance MBR diagnostics to show pairing success
3. Add pairing validation to build process

## Risk Assessment

### Low Risk
- Function already exists and was designed for this purpose
- Modification is late in pipeline (doesn't affect earlier processing)
- Easy to revert if issues arise

### Medium Risk
- Changes final precursor table structure
- May affect downstream code expecting current format
- Need to verify all systems can handle partner_precursor_idx column

### Mitigation
- Extensive testing on small datasets first
- Maintain backward compatibility where possible
- Document all changes for easier troubleshooting

## Conclusion

The target-decoy pairing bug is fixable using existing code that was designed for this exact purpose. The `add_precursor_partner_columns!` function properly creates charge-specific pairs and provides efficient partner lookup, addressing both the immediate pairing problem and enabling proper MBR functionality.

The fix is surgical - it doesn't change the overall workflow, just corrects the final pair_id assignment and adds partner indexing. This should restore the intended 1:1 target/decoy ratio and enable proper CV fold assignment without data leakage.