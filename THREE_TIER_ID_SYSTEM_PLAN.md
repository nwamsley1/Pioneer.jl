# Three-Tier ID System Implementation Plan

## Overview
Implement a three-tier ID system for proper tracking and pairing of peptides through the library building pipeline:
- **base_seq_id**: Tracks base sequences (assigned after combine_shared_peptides)
- **base_pep_id**: Tracks peptides with modifications (assigned after add_mods)
- **base_prec_id**: Tracks precursors with charge states (assigned after add_charge)

## Current Problem
The current two-tier system (base_pep_id and base_prec_id) is insufficient because:
1. We need to track sequences at three different levels of the pipeline
2. Entrapment sequences need to inherit from base sequences
3. Modifications create peptide variants
4. Charge states create precursor variants
5. Pairing needs to consider both peptide and precursor levels

## Proposed Solution: Three-Tier ID System

### 1. Update FastaEntry Structure
**File**: `src/Structs/FastaEntry.jl`

Add new field to FastaEntry struct:
```julia
struct FastaEntry
    id::String
    description::String
    gene::String
    protein::String 
    organism::String
    proteome::String
    sequence::String
    start_idx::UInt32
    structural_mods::StructuralModifications
    isotopic_mods::IsotopicModifications
    charge::UInt8
    base_seq_id::UInt32      # NEW: Base sequence ID
    base_pep_id::UInt32      # EXISTING: Peptide ID (with mods)
    base_prec_id::UInt32     # EXISTING: Precursor ID (with charge)
    entrapment_pair_id::UInt8
    decoy::Bool
end
```

Update all getter functions:
- Add `get_base_seq_id(entry::FastaEntry)`
- Keep existing `get_base_pep_id(entry::FastaEntry)`
- Keep existing `get_base_prec_id(entry::FastaEntry)`

### 2. Create ID Assignment Functions
**File**: `src/Routines/BuildSpecLib/fasta/fasta_utils.jl`

#### 2.1 Create `assign_base_seq_ids!` function
```julia
function assign_base_seq_ids!(fasta_entries::Vector{FastaEntry})
    """
    Assign sequential base_seq_id values starting from 1.
    Called after combine_shared_peptides to identify unique base sequences.
    """
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        fasta_entries[i] = FastaEntry(
            # ... all existing fields ...
            UInt32(i),               # base_seq_id - sequential
            get_base_pep_id(entry),  # preserve existing
            get_base_prec_id(entry), # preserve existing
            # ... remaining fields ...
        )
    end
    return length(fasta_entries)
end
```

#### 2.2 Update `reassign_base_pep_ids!` function
```julia
function assign_base_pep_ids!(fasta_entries::Vector{FastaEntry})
    """
    Assign sequential base_pep_id values starting from 1.
    Called after add_mods to identify unique peptides (sequence + mods).
    """
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        fasta_entries[i] = FastaEntry(
            # ... all existing fields ...
            get_base_seq_id(entry),  # preserve base_seq_id
            UInt32(i),               # base_pep_id - sequential
            get_base_prec_id(entry), # preserve existing
            # ... remaining fields ...
        )
    end
    return length(fasta_entries)
end
```

#### 2.3 Create `assign_base_prec_ids!` function
```julia
function assign_base_prec_ids!(fasta_entries::Vector{FastaEntry})
    """
    Assign sequential base_prec_id values starting from 1.
    Called after add_charge to identify unique precursors (peptide + charge).
    """
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        fasta_entries[i] = FastaEntry(
            # ... all existing fields ...
            get_base_seq_id(entry),  # preserve base_seq_id
            get_base_pep_id(entry),  # preserve base_pep_id
            UInt32(i),               # base_prec_id - sequential
            # ... remaining fields ...
        )
    end
    return length(fasta_entries)
end
```

### 3. Update Pipeline Order
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

New pipeline order:
```julia
# Step 1: Combine shared peptides (I/L equivalence)
fasta_entries = combine_shared_peptides(fasta_entries)

# Step 2: Assign base_seq_id for sequence tracking
seq_count = assign_base_seq_ids!(fasta_entries)
@user_info "Base Seq ID Assignment: $seq_count unique base sequences"

# Step 3: Add entrapment sequences (inherit base_seq_id)
fasta_entries = add_entrapment_sequences(fasta_entries, ...)

# Step 4: Add modifications (creates peptide variants)
fasta_entries = add_mods(fasta_entries, ...)

# Step 5: Assign base_pep_id for peptide tracking
pep_count = assign_base_pep_ids!(fasta_entries)
@user_info "Base Pep ID Assignment: $pep_count unique peptides"

# Step 6: Add decoy sequences (inherit base_pep_id)
fasta_entries = add_decoy_sequences(fasta_entries)

# Step 7: Add charge states (creates precursor variants)
fasta_entries = add_charge(fasta_entries, ...)

# Step 8: Assign base_prec_id for precursor tracking
prec_count = assign_base_prec_ids!(fasta_entries)
@user_info "Base Prec ID Assignment: $prec_count unique precursors"

# Step 9: Build DataFrame with all IDs
seq_df = build_fasta_df(fasta_entries, ...)

# Step 10: Pair using combination of base_pep_id and base_prec_id
seq_df = add_charge_specific_partner_columns!(seq_df)
```

### 4. Update DataFrame Creation
**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

In `build_fasta_df` function, add base_seq_id to vectors and DataFrame:
```julia
# Add to vector allocations
_base_seq_id = Vector{UInt32}(undef, prec_alloc_size)

# In the loop
_base_seq_id[n] = get_base_seq_id(peptide)

# In DataFrame creation
seq_df = DataFrame(
    # ... existing columns ...
    base_seq_id = _base_seq_id[1:n],
    base_pep_id = _base_pep_id[1:n],
    base_prec_id = _base_prec_id[1:n],
    pair_id = _pair_id[1:n]  # Will be assigned by pairing function
)
```

### 5. Update Pairing Function
**File**: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`

Update `add_charge_specific_partner_columns!` to use combination of base_pep_id and base_prec_id:
```julia
function add_charge_specific_partner_columns!(df::DataFrame)
    """
    Creates target-decoy pairs using (base_pep_id, base_prec_id) combination.
    This ensures proper pairing at both peptide and precursor levels.
    """
    n = nrow(df)
    
    # Create lookup: (base_pep_id, base_prec_id, is_target) -> row_index
    lookup = Dict{Tuple{UInt32, UInt32, Bool}, Int}()
    for (idx, row) in enumerate(eachrow(df))
        key = (row.base_pep_id, row.base_prec_id, !row.decoy)
        lookup[key] = idx
    end
    
    # Initialize pair_id vector
    pair_ids = Vector{Union{UInt32, Missing}}(missing, n)
    next_pair_id = UInt32(1)
    
    # Process targets to find their decoy partners
    for idx in 1:n
        row = df[idx, :]
        if row.decoy
            continue
        end
        
        # Look for decoy partner with same base_pep_id and base_prec_id
        decoy_key = (row.base_pep_id, row.base_prec_id, false)
        
        if haskey(lookup, decoy_key)
            decoy_idx = lookup[decoy_key]
            # Assign same pair_id to both
            pair_ids[idx] = next_pair_id
            pair_ids[decoy_idx] = next_pair_id
            next_pair_id += 1
        else
            # Unpaired target
            pair_ids[idx] = next_pair_id
            next_pair_id += 1
        end
    end
    
    # Assign unique pair_ids to unpaired decoys
    for idx in 1:n
        if ismissing(pair_ids[idx])
            pair_ids[idx] = next_pair_id
            next_pair_id += 1
        end
    end
    
    # Add pair_id column to DataFrame
    df.pair_id = pair_ids
    
    return df
end
```

### 6. Update Functions to Preserve IDs

#### 6.1 Update `add_entrapment_sequences`
Ensure it preserves base_seq_id from the target sequence.

#### 6.2 Update `add_decoy_sequences`  
Ensure it preserves both base_seq_id and base_pep_id from the target.

#### 6.3 Update `add_charge`
Ensure it preserves base_seq_id and base_pep_id.

## Implementation Steps

### Phase 1: Update FastaEntry Structure
1. Modify FastaEntry struct to include base_seq_id
2. Update constructor and getter functions
3. Update all existing FastaEntry creation sites

### Phase 2: Create ID Assignment Functions
1. Implement assign_base_seq_ids!
2. Rename/update reassign_base_pep_ids! to assign_base_pep_ids!
3. Rename/update reassign_base_prec_ids! to assign_base_prec_ids!

### Phase 3: Update Pipeline
1. Modify chronologer_prep.jl pipeline order
2. Add logging for each ID assignment step
3. Update build_fasta_df to include all three ID columns

### Phase 4: Update Pairing Logic
1. Modify add_charge_specific_partner_columns!
2. Use (base_pep_id, base_prec_id) combination for pairing
3. Test pairing accuracy

### Phase 5: Validation
1. Verify all IDs are properly preserved through pipeline
2. Check that pairing correctly matches targets with decoys
3. Ensure modifications are properly handled
4. Validate that entrapment tracking still works

## Expected Outcomes

### ID Hierarchy
- **base_seq_id**: Groups all variants of the same base sequence
  - Shared across: entrapments, modifications, decoys, charges
  
- **base_pep_id**: Groups all charge variants of the same peptide
  - Shared across: decoys, charges
  - Unique per: modification variant
  
- **base_prec_id**: Unique per precursor
  - Unique per: charge state
  
- **pair_id**: Links target-decoy pairs
  - Based on: (base_pep_id, base_prec_id) combination

### Final Output
The precursors_table.arrow will contain:
- `base_seq_id`: Base sequence identifier
- `base_pep_id`: Peptide identifier (with modifications)
- `base_prec_id`: Precursor identifier (with charge)
- `pair_id`: Target-decoy pair identifier

## Testing Strategy

### Unit Tests
1. Test each ID assignment function independently
2. Verify ID preservation through each pipeline step
3. Test pairing logic with known examples

### Integration Tests
1. Build test library with known sequences
2. Verify correct ID assignments at each level
3. Check that modification pairing works correctly
4. Validate entrapment tracking

### Validation Metrics
1. No sequences with different modifications share base_pep_id
2. No peptides with different charges share base_prec_id
3. All target-decoy pairs have matching modifications
4. Entrapment sequences properly tracked via base_seq_id

## Risk Assessment

### Low Risk
- Adding new field to FastaEntry is straightforward
- ID assignment functions are simple sequential assignments
- Pipeline order is logical and clear

### Medium Risk
- Ensuring all FastaEntry creation sites are updated
- Verifying ID preservation through all functions
- Testing comprehensive pairing scenarios

### Mitigation
- Careful testing at each implementation phase
- Preserve existing functionality while adding new
- Clear logging to track ID assignments

This three-tier system provides clear separation of concerns and enables proper tracking at each level of the peptide library building process.