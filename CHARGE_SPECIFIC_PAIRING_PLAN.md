# Charge-Specific Target-Decoy Pairing Plan

## Executive Summary

Instead of fixing the broken `pair_id` assignment directly, we'll use the correctly preserved `base_pep_id` as a shared identifier across charge states, then create charge-specific pairs by matching on `(base_pep_id, charge, target/decoy)` combinations. This avoids relying on sequence reversal and works with any decoy generation method.

## The Elegant Solution

### Current State (What Works)
- `base_pep_id` is correctly shared between target and decoy versions of the same peptide
- `base_pep_id` is preserved through all processing steps (digestion → decoy generation → charge expansion)
- Each charge state of the same peptide shares the same `base_pep_id`

### Current Problem (What's Broken)
- `pair_id` becomes unique per precursor instead of shared between target/decoy pairs
- Can't rely on sequence reversal because decoy generation method may have changed

### The Solution Strategy
1. **Keep** `base_pep_id` as the shared peptide identifier
2. **Use** `base_pep_id` initially for `pair_id` (groups all charge states together)
3. **Create** charge-specific pairs using `(base_pep_id, charge)` matching
4. **Assign** unique `pair_id` to each charge-specific target-decoy pair

## Detailed Implementation Plan

### Phase 1: Fix Initial pair_id Assignment

**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
**Line**: 554

**Change from**:
```julia
_pair_id[n] = get_base_prec_id(peptide)  # Unique per precursor
```

**Change to**:
```julia
_pair_id[n] = get_base_pep_id(peptide)   # Shared across charge states
```

**Result**: All charge states of the same peptide (both target and decoy) now share the same initial `pair_id`.

### Phase 2: Create Modified Pairing Function

**File**: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`

Create a new function `add_charge_specific_partner_columns!` based on the existing commented function:

```julia
function add_charge_specific_partner_columns!(df::DataFrame)
    """
    Creates charge-specific target-decoy pairs using (base_pep_id, charge) matching.
    Each target finds its decoy partner with same base_pep_id and charge.
    Assigns unique pair_id to each charge-specific pair.
    """
    result_df = copy(df)
    n = nrow(df)
    
    # Create lookup dictionary: (base_pep_id, charge, is_target) -> row_index
    lookup = Dict{Tuple{UInt32, UInt8, Bool}, Int}()
    for (idx, row) in enumerate(eachrow(df))
        key = (row.base_pep_id, row.precursor_charge, !row.decoy)  # !decoy = is_target
        lookup[key] = idx
    end
    
    # Initialize output vectors
    partner_indices = Vector{Union{UInt32, Missing}}(missing, n)
    new_pair_ids = Vector{UInt32}(undef, n)
    
    # Process each target to find its decoy partner
    pair_id_counter = UInt32(0)
    processed_pairs = Set{Tuple{Int, Int}}()
    
    for idx in 1:n
        row = df[idx, :]
        
        # Skip if already processed or if it's a decoy (process targets first)
        if row.decoy
            continue
        end
        
        # Look for decoy partner with same base_pep_id and charge
        target_key = (row.base_pep_id, row.precursor_charge, true)   # is_target = true
        decoy_key = (row.base_pep_id, row.precursor_charge, false)   # is_target = false
        
        if haskey(lookup, decoy_key)
            decoy_idx = lookup[decoy_key]
            
            # Avoid double processing
            pair_tuple = (min(idx, decoy_idx), max(idx, decoy_idx))
            if pair_tuple ∉ processed_pairs
                pair_id_counter += 1
                
                # Assign same pair_id to both
                new_pair_ids[idx] = pair_id_counter
                new_pair_ids[decoy_idx] = pair_id_counter
                
                # Set partner references
                partner_indices[idx] = UInt32(decoy_idx)
                partner_indices[decoy_idx] = UInt32(idx)
                
                push!(processed_pairs, pair_tuple)
            end
        else
            # No decoy partner found - assign unique pair_id
            pair_id_counter += 1
            new_pair_ids[idx] = pair_id_counter
        end
    end
    
    # Handle any remaining decoys without targets
    for idx in 1:n
        if df[idx, :].decoy && ismissing(partner_indices[idx])
            pair_id_counter += 1
            new_pair_ids[idx] = pair_id_counter
        end
    end
    
    # Update the DataFrame
    result_df.pair_id = new_pair_ids
    result_df.partner_precursor_idx = partner_indices
    
    return result_df
end
```

### Phase 3: Integrate into Processing Pipeline

**File**: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`

**After DataFrame creation** (around line 571):
```julia
# Create initial DataFrame with base_pep_id as pair_id
seq_df = DataFrame(
    (proteome_identifiers = _proteome_identifiers[1:n],
     sequence = _sequence[1:n],
     # ... other columns ...
     base_pep_id = _base_pep_id[1:n],
     pair_id = _pair_id[1:n],     # Initially same as base_pep_id
     decoy = _decoy[1:n],
     precursor_charge = _precursor_charge[1:n])
)

# Apply charge-specific pairing
seq_df = add_charge_specific_partner_columns!(seq_df)

# Remove the old add_pair_indices! call if present
# add_pair_indices!(seq_df)  # REMOVE THIS LINE
```

## How the Matching Logic Works

### Input State
After Phase 1, precursors look like this:
```
Row  | base_pep_id | charge | decoy | pair_id | sequence
-----|-------------|--------|-------|---------|----------
1    | 10          | 2      | false | 10      | PEPTIDE
2    | 10          | 2      | true  | 10      | EDITPEP  
3    | 10          | 3      | false | 10      | PEPTIDE
4    | 10          | 3      | true  | 10      | EDITPEP
5    | 11          | 2      | false | 11      | ANOTHER
6    | 11          | 2      | true  | 11      | REHTONA
```

### Matching Process
1. **Row 1** (target, base_pep_id=10, charge=2):
   - Look for: (base_pep_id=10, charge=2, is_target=false)
   - **Found**: Row 2 (decoy, base_pep_id=10, charge=2)
   - **Assign**: Both get new_pair_id = 1
   - **Link**: partner_indices[1] = 2, partner_indices[2] = 1

2. **Row 3** (target, base_pep_id=10, charge=3):
   - Look for: (base_pep_id=10, charge=3, is_target=false)
   - **Found**: Row 4 (decoy, base_pep_id=10, charge=3)
   - **Assign**: Both get new_pair_id = 2
   - **Link**: partner_indices[3] = 4, partner_indices[4] = 3

3. **Row 5** (target, base_pep_id=11, charge=2):
   - Look for: (base_pep_id=11, charge=2, is_target=false)
   - **Found**: Row 6 (decoy, base_pep_id=11, charge=2)
   - **Assign**: Both get new_pair_id = 3

### Final Output State
```
Row  | base_pep_id | charge | decoy | pair_id | partner_idx | sequence
-----|-------------|--------|-------|---------|-------------|----------
1    | 10          | 2      | false | 1       | 2           | PEPTIDE
2    | 10          | 2      | true  | 1       | 1           | EDITPEP  
3    | 10          | 3      | false | 2       | 4           | PEPTIDE
4    | 10          | 3      | true  | 2       | 3           | EDITPEP
5    | 11          | 2      | false | 3       | 6           | ANOTHER
6    | 11          | 2      | true  | 3       | 5           | REHTONA
```

## Advantages of This Approach

### 1. Robust to Decoy Generation Changes
- Doesn't rely on sequence reversal
- Works with any decoy generation method
- Only requires that target/decoy share `base_pep_id`

### 2. Charge-State Specific
- Each charge state creates independent pairs
- No confusion between different charge states
- Proper 2-member pairs for each pair_id

### 3. Direct Partner Lookup
- `partner_precursor_idx` provides O(1) access
- Perfect for MBR implementation
- Easy validation and debugging

### 4. Backward Compatible
- Uses existing `base_pep_id` infrastructure
- Doesn't change earlier processing steps
- Easy to validate and revert if needed

## Validation and Testing

### Expected Results
```julia
# Test 1: Each pair_id should have exactly 2 members
gptable = groupby(ptable, :pair_id)
pair_sizes = [nrow(subdf) for (key, subdf) in pairs(gptable)]
@test all(pair_sizes .== 2)  # All pairs have exactly 2 members

# Test 2: Partners should reference each other
for i in 1:nrow(ptable)
    partner_idx = ptable[i, :partner_precursor_idx]
    if !ismissing(partner_idx)
        @test ptable[partner_idx, :partner_precursor_idx] == i
        @test ptable[i, :pair_id] == ptable[partner_idx, :pair_id]
    end
end

# Test 3: Target/decoy pairs should have same base_pep_id and charge
for (pair_id, subdf) in pairs(groupby(ptable, :pair_id))
    if nrow(subdf) == 2
        @test subdf[1, :base_pep_id] == subdf[2, :base_pep_id]
        @test subdf[1, :precursor_charge] == subdf[2, :precursor_charge]
        @test subdf[1, :decoy] != subdf[2, :decoy]  # One target, one decoy
    end
end
```

### MBR Impact Verification
```julia
# Test that MBR can now find partners
precursor_dict = get_precursor_dict_after_first_pass()
partner_found_count = 0
partner_added_count = 0

for pid in keys(precursor_dict)
    partner_pid = getPartnerPrecursorIdx(precursors)[pid]
    if !ismissing(partner_pid)
        partner_found_count += 1
        if !haskey(precursor_dict, partner_pid)
            # MBR should successfully add this partner
            partner_added_count += 1
        end
    end
end

@test partner_found_count > 0  # Should find many partners
# After MBR, should have close to 1:1 ratio
```

## Implementation Timeline

### Phase 1: Setup and Initial Commit
1. Update this plan with testing strategy
2. Commit current state as baseline
3. Create BuildSpecLib JSON for keap1.fasta testing

### Phase 2: Core Implementation
1. Modify line 554 in `chronologer_prep.jl`
2. Create `add_charge_specific_partner_columns!` function
3. Integrate into processing pipeline
4. Commit after each major change

### Phase 3: Iterative Testing with Revise.jl
1. Use `] activate .` then `using Revise, Pioneer`
2. Test BuildSpecLib("test_params.json") after each change
3. Fix issues iteratively without restarting Julia
4. Commit working versions

### Phase 4: Validation
1. Inspect final precursors.arrow for proper pairing
2. Verify each pair_id has exactly 2 members
3. Check partner_precursor_idx references

## Testing Strategy

### Test Dataset
- **FASTA**: `/data/fasta/keap1.fasta` 
- **Method**: Altimeter prediction models
- **Output**: Library will be created in specified output directory

### BuildSpecLib JSON Configuration
```json
{
    "output_dir": "/tmp/keap1_test_library",
    "lib_name": "keap1_test",
    "fasta_dir": "data/fasta/keap1.fasta",
    "prediction_model": "altimeter",
    "nce_params": {
        "NCE": 25.0,
        "charge": 2,
        "dynamic": true
    }
}
```

### Testing Workflow
1. Make code change
2. Julia already running with `using Revise, Pioneer`
3. Run `BuildSpecLib("test_keap1_params.json")`
4. Check for errors, fix, repeat
5. When successful, inspect output library files
6. Commit working version

### Expected Library Output Location
- **Main directory**: `/tmp/keap1_test_library/keap1_test.poin/`
- **Key files**:
  - `precursors.arrow` (main file to inspect for pairing)
  - `config.json` (build parameters)
  - `build_log.txt` (error messages if any)

### Validation Commands
```julia
using DataFrames, Arrow, Tables

# Load and inspect the library
ptable = DataFrame(Tables.columntable(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors.arrow")))

# Check pairing
gptable = groupby(ptable, :pair_id)
pair_sizes = [nrow(subdf) for (key, subdf) in pairs(gptable)]
println("Pair sizes: ", unique(pair_sizes))  # Should be [2] after fix

# Check partner references
missing_partners = sum(ismissing.(ptable.partner_precursor_idx))
println("Missing partners: $missing_partners / $(nrow(ptable))")

# Check target/decoy balance
n_targets = sum(.!ptable.decoy)
n_decoys = sum(ptable.decoy)
println("Targets: $n_targets, Decoys: $n_decoys, Ratio: $(n_targets/n_decoys)")
```

## Risk Mitigation

### Low Risk Changes
- Only modifies final pair_id assignment
- Uses existing, well-tested `base_pep_id` values
- Easy to revert by changing one line

### Validation Strategy
- Extensive automated testing
- Manual inspection of small datasets
- A/B testing with old vs new pairing

### Fallback Plan
- Can immediately revert line 554 change
- Old behavior preserved until new function is called
- Gradual rollout possible

## Expected Impact

### MBR Functionality
- **Before**: Cannot find partners, poor target/decoy ratios
- **After**: Direct partner lookup, balanced ratios close to 1:1

### CV Fold Assignment
- **Before**: Random assignment, potential data leakage
- **After**: Can ensure pairs stay together, no leakage

### Performance
- **Before**: MBR searches fail, wasted computation
- **After**: O(1) partner lookup, efficient MBR processing

This approach is elegant, robust, and solves the core problem without relying on specific decoy generation methods.