# Target-Decoy Pairing Bug in Pioneer.jl BuildSpecLib

## Executive Summary

The target-decoy pairing mechanism in Pioneer.jl is fundamentally broken. Each precursor has a unique `pair_id` instead of sharing IDs with their target/decoy partner. This breaks Match-Between-Runs (MBR) functionality and prevents achieving the expected 1:1 target-to-decoy ratio.

## The Bug

**Location**: `chronologer_prep.jl`, line 554
```julia
_pair_id[n] = get_base_prec_id(peptide)  # BUG: Uses unique ID instead of shared pair ID
```

The `pair_id` column in the final Arrow table is incorrectly set to `base_prec_id`, which becomes unique per precursor after charge state expansion.

## Data Flow and Column Mapping

### Final Arrow Table Columns
- **`pair_id`** ← `get_base_prec_id(peptide)` (UNIQUE per precursor - **THIS IS THE BUG**)
- **`base_pep_id`** ← `get_base_pep_id(peptide)` (shared between target/decoy pairs)
- **`entrapment_group_id`** ← `get_entrapment_pair_id(peptide)` (always 0 for regular sequences)
- **`partner_precursor_idx`** ← Set by `add_pair_indices!` function (attempts to find pairs but fails due to unique pair_ids)

## The Workflow

### 1. Protein Digestion (`fasta_digest.jl`)
```julia
# Each target peptide gets unique sequential IDs
base_pep_id = 1, 2, 3, ...
base_prec_id = 1, 2, 3, ...
entrapment_pair_id = 0  # Always 0 for regular sequences
```

### 2. Decoy Generation (`add_decoy_sequences` in `fasta_utils.jl`)
```julia
# Decoys INHERIT the same IDs from their target
decoy_fasta_entries[n] = FastaEntry(
    ...
    get_base_pep_id(target_entry),      # Shared with target ✓
    get_base_prec_id(target_entry),     # Shared with target ✓
    get_entrapment_pair_id(target_entry), # Shared with target (0)
    true  # is_decoy flag
)
```
**Result**: Target and decoy share `base_pep_id` and `base_prec_id`

### 3. Entrapment Sequences (`add_entrapment_sequences` in `fasta_utils.jl`)
- Creates shuffled versions of target sequences for FDR validation
- Each entrapment group gets a unique `entrapment_pair_id` (1, 2, 3, ...)
- Does NOT interfere with regular target-decoy pairing
- Only affects sequences marked with non-zero `entrapment_group_id`

### 4. Charge State Expansion (`add_charge` in `chronologer_prep.jl`)
```julia
base_prec_id = one(UInt32)  # RESETS to 1!
for fasta_peptide in fasta_peptides
    for charge in range(min_charge, max_charge)
        push!(fasta_peptides_wcharge,
            FastaEntry(
                ...
                get_base_pep_id(fasta_peptide),  # Preserved ✓
                base_prec_id,                     # REASSIGNED! ✗
                get_entrapment_pair_id(fasta_peptide), # Preserved
                is_decoy(fasta_peptide)
            ))
        base_prec_id += one(UInt32)  # Sequential increment
    end
end
```
**Problem**: `base_prec_id` is reassigned sequentially, breaking the shared IDs between targets and decoys!

### 5. Final Table Creation (`prepare_chronologer_input` in `chronologer_prep.jl`)
```julia
# Line 554 - THE BUG:
_pair_id[n] = get_base_prec_id(peptide)  # Now unique per precursor!

# Creates DataFrame with:
DataFrame(
    ...
    base_pep_id = _base_pep_id[1:n],      # Still shared
    pair_id = _pair_id[1:n]                # Unique per precursor!
)
```

The `prepare_chronologer_input` function orchestrates the entire workflow:
1. Calls digestion, decoy generation, and entrapment addition
2. Applies modifications
3. Expands charge states (where the bug occurs)
4. Creates the final precursor table with incorrect pair_ids

### 6. Partner Index Assignment (`add_pair_indices!` in `chronologer_prep.jl`)
```julia
# Attempts to find pairs based on pair_id
# But since each precursor has unique pair_id, no pairs are found!
if length(rows) == 2  # Never true because each pair_id is unique
    other_idx = (rows[1] == i) ? rows[2] : rows[1]
    precursor_pair_idx[i] = other_idx
end
```

## Impact

1. **No functional target-decoy pairs**: Each precursor appears unpaired
2. **MBR fails**: Cannot add partner precursors because partners aren't identified
3. **Incorrect ratios**: Instead of 1:1 target:decoy ratio after MBR, we see 4.5:1 or worse
4. **CV fold splitting**: Cannot properly assign pairs to same fold
5. **Data leakage**: When pairs are split across CV folds (which happens randomly)

## Evidence

```julia
julia> ptable = DataFrame(Tables.columntable(Arrow.Table("precursors.arrow")))
julia> gptable = groupby(ptable, :pair_id)
julia> unique([size(subdf,1) for (key, subdf) in pairs(gptable)])
1-element Vector{Int64}: 1  # Each pair_id has only 1 member!
```

## Solution Options

### Option 1: Use `base_pep_id` for pairing
```julia
# In prepare_chronologer_input, line 554:
_pair_id[n] = get_base_pep_id(peptide)  # Use peptide ID instead
```

### Option 2: Preserve `base_prec_id` during charge expansion
```julia
# In add_charge function:
# Don't reassign base_prec_id, use the original value
base_prec_id = get_base_prec_id(fasta_peptide)
```

### Option 3: Implement proper pair ID assignment
- Create unique pair IDs during initial digestion
- Ensure both target and decoy get the same pair ID
- Preserve through all transformations

### Option 4: Use the commented-out pairing function
- Uncomment `add_precursor_partner_columns!` in `pair_decoys.jl`
- This function properly matches targets to decoys and assigns shared pair_ids

## Recommendation

**Immediate fix**: Change line 554 in `chronologer_prep.jl`:
```julia
# Current (broken):
_pair_id[n] = get_base_prec_id(peptide)

# Fixed:
_pair_id[n] = get_base_pep_id(peptide)
```

This would ensure target-decoy pairs share the same `pair_id` since `base_pep_id` is correctly preserved through all transformations.

**Long-term fix**: Implement a robust pairing system that:
1. Assigns unique pair IDs during peptide generation
2. Maintains pairing through all transformations
3. Validates pairing integrity before writing final tables
4. Properly handles entrapment sequences separately from regular target-decoy pairs