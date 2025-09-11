# Plan: Handling Modified Peptides in Decoy Generation

## Problem Statement

Currently, the decoy generation system treats each modified variant of a peptide as a separate entity when checking for duplicates. This causes issues:

1. **False Duplicates**: When generating a decoy for an oxidized peptide (e.g., with Met oxidation), the reversed/shuffled sequence might match the unmodified version, causing unnecessary fallback to shuffle.

2. **Inconsistent Decoys**: Different modification states of the same peptide might get different decoy sequences, which is incorrect. All modification variants should share the same decoy sequence.

3. **Example Case**:
   - Original: `PEPTMIDER` (unmodified)
   - Original: `PEPTM[+16]IDER` (oxidized Met)
   - Reversed: `REDIMTPEP` 
   - Both should get the SAME reversed sequence, with modifications adjusted accordingly

## Current Flow (Problematic)

```
1. Process PEPTMIDER → Reverse to REDIMTPEP → Check if exists → Add
2. Process PEPTM[+16]IDER → Reverse to REDIM[+16]TPEP → Check if exists
   → Might match REDIMTPEP (without considering mods) → False duplicate!
```

## Solution Design

### Key Principle
**Decoy generation should be based on the BASE PEPTIDE SEQUENCE, not the modified sequence.**

### Implementation Strategy

#### Phase 1: Group by Base Peptide
1. Before generating decoys, group all peptides by their base sequence (ignoring modifications)
2. Use the existing `base_pep_id` which already tracks peptide identity (sequence + mods)
3. Create a new grouping based on sequence alone

#### Phase 2: Generate One Decoy Per Base Sequence
1. For each unique base sequence:
   - Generate ONE decoy sequence (reverse or shuffle)
   - Store this mapping: `base_sequence → decoy_sequence`
2. Check for duplicates ONLY against other base sequences, not modification variants

#### Phase 3: Apply Decoy to All Variants
1. For each peptide with modifications:
   - Look up the decoy sequence for its base sequence
   - Apply the same decoy sequence
   - Adjust modification positions according to the transformation

### Detailed Algorithm

```julia
function add_decoy_sequences_grouped(target_fasta_entries; method="reverse")
    # Step 1: Group by base sequence
    sequence_groups = Dict{String, Vector{FastaEntry}}()
    for entry in target_fasta_entries
        base_seq = get_sequence(entry)
        if !haskey(sequence_groups, base_seq)
            sequence_groups[base_seq] = Vector{FastaEntry}()
        end
        push!(sequence_groups[base_seq], entry)
    end
    
    # Step 2: Generate decoys for unique base sequences
    base_to_decoy = Dict{String, String}()
    used_decoys = Set{Tuple{String, UInt8}}()  # (sequence, charge)
    
    for (base_seq, entries) in sequence_groups
        # Get representative charge (might vary, handle all charges)
        charges = unique([get_charge(e) for e in entries])
        
        # Generate decoy for this base sequence
        decoy_seq = generate_decoy(base_seq, method)
        
        # Check if this decoy conflicts with other BASE sequences
        for charge in charges
            if (decoy_seq, charge) in used_decoys
                # Fall back to shuffle
                decoy_seq = generate_decoy(base_seq, "shuffle")
                # Keep trying until unique
            end
        end
        
        base_to_decoy[base_seq] = decoy_seq
        for charge in charges
            push!(used_decoys, (decoy_seq, charge))
        end
    end
    
    # Step 3: Create decoy entries for all variants
    decoy_entries = Vector{FastaEntry}()
    for entry in target_fasta_entries
        base_seq = get_sequence(entry)
        decoy_seq = base_to_decoy[base_seq]
        
        # Adjust modifications for the decoy sequence
        adjusted_mods = adjust_mod_positions(entry, decoy_seq)
        
        # Create decoy entry
        push!(decoy_entries, create_decoy_entry(entry, decoy_seq, adjusted_mods))
    end
    
    return decoy_entries
end
```

### Special Considerations

#### 1. Deterministic Shuffling
For the shuffle method to produce consistent results for the same base sequence:
- Use a seeded random number generator based on the sequence hash
- Or store the shuffle result and reuse it

#### 2. Position Mapping
- The existing `shuffle_seq.new_positions` tracking already handles this
- Need to ensure it's preserved and reused for all modification variants

#### 3. Performance
- Grouping adds minimal overhead
- Actually IMPROVES performance by reducing redundant decoy generation

#### 4. Edge Cases
- **Different charges**: Same sequence with different charges should get the same decoy sequence
- **Multiple modification sites**: All combinations should map to the same decoy
- **Entrapment sequences**: Apply the same logic

## Entrapment Sequence Handling

### The Entrapment Challenge
Entrapment sequences have the same fundamental issue as decoys, but with additional complexity:
- They need to be grouped by base sequence BEFORE generating entrapments
- Multiple entrapment sequences (if entrapment_r > 1) need to be different from each other
- But all modification variants of the same peptide should get the same set of entrapment sequences

### Entrapment Algorithm

```julia
function add_entrapment_sequences_grouped(target_fasta_entries; method="shuffle", entrapment_r=1)
    # Step 1: Group by base sequence (same as decoys)
    sequence_groups = Dict{String, Vector{FastaEntry}}()
    for entry in target_fasta_entries
        base_seq = get_sequence(entry)
        if !haskey(sequence_groups, base_seq)
            sequence_groups[base_seq] = Vector{FastaEntry}()
        end
        push!(sequence_groups[base_seq], entry)
    end
    
    # Step 2: Generate entrapments for unique base sequences
    base_to_entrapments = Dict{String, Vector{String}}()
    used_sequences = Set{Tuple{String, UInt8}}()  # Track all used sequences
    
    for (base_seq, entries) in sequence_groups
        charges = unique([get_charge(e) for e in entries])
        entrapments = Vector{String}()
        
        # Generate entrapment_r unique entrapments for this base
        for i in 1:entrapment_r
            attempts = 0
            while attempts < max_attempts
                entrap_seq = generate_entrapment(base_seq, method)
                
                # Check uniqueness against all charges
                is_unique = true
                for charge in charges
                    if (entrap_seq, charge) in used_sequences
                        is_unique = false
                        break
                    end
                end
                
                if is_unique
                    push!(entrapments, entrap_seq)
                    for charge in charges
                        push!(used_sequences, (entrap_seq, charge))
                    end
                    break
                end
                attempts += 1
            end
        end
        
        base_to_entrapments[base_seq] = entrapments
    end
    
    # Step 3: Create entrapment entries for all variants
    entrapment_entries = Vector{FastaEntry}()
    for entry in target_fasta_entries
        base_seq = get_sequence(entry)
        entrapments = base_to_entrapments[base_seq]
        
        for (i, entrap_seq) in enumerate(entrapments)
            # Adjust modifications for the entrapment sequence
            adjusted_mods = adjust_mod_positions(entry, entrap_seq)
            
            # Create entrapment entry with group_id = i
            push!(entrapment_entries, create_entrapment_entry(
                entry, entrap_seq, adjusted_mods, entrapment_group_id=i
            ))
        end
    end
    
    return entrapment_entries
end
```

## Key Insight: Using Existing IDs

Looking at the chronologer_prep.jl workflow:
1. **Step 2**: Modifications are added
2. **Step 3**: `base_pep_id` is assigned (tracks sequence + modifications)
3. **Step 4**: Entrapments are added
4. **Step 5**: `base_target_id` is assigned
5. **Step 6**: Decoys are added

The issue is that `base_pep_id` includes modifications, so it doesn't group modification variants together. We need a new ID that tracks just the base sequence.

### Proposed New ID: `base_sequence_id`
- Assigned after modifications are added but before entrapments/decoys
- Groups all entries with the same amino acid sequence
- Ignores modification state
- Used as the grouping key for consistent decoy/entrapment generation

## Implementation Steps

### Step 1: Refactor Data Structure
- Create a `PeptideGroup` structure to hold base sequence and all variants
- Track the decoy sequence once per group

### Step 2: Modify Duplicate Detection
- Only check for duplicates against OTHER base peptides
- Allow "duplicates" within the same peptide's modification variants

### Step 3: Update Functions
1. `add_decoy_sequences()` - Main refactor needed
2. `add_entrapment_sequences()` - Apply same logic
3. Consider creating shared helper functions

### Step 4: Testing
- Test with peptides having multiple oxidation states
- Test with multiple modification types
- Verify all variants get the same decoy
- Verify position adjustments work correctly

## Benefits

1. **Correctness**: All modification variants properly share the same decoy
2. **Performance**: Generate each decoy only once per unique sequence
3. **Reduced Fallbacks**: No false duplicates from modification variants
4. **Consistency**: Deterministic behavior for all variants

## Potential Issues to Watch

1. **Memory**: Storing mappings might increase memory usage slightly
2. **Complexity**: Code becomes more complex but more correct
3. **Backward Compatibility**: Need to ensure existing libraries still work

## Alternative Simpler Approach

If full refactoring is too complex, a simpler fix:
1. When checking for duplicates, compare only the base sequence (ignore mods)
2. But track the full sequence+mods for uniqueness
3. This is less elegant but might be easier to implement

## Recommendation

Implement the full solution (Phase 1-3) because:
- It's the correct approach scientifically
- It will prevent future issues
- It actually simplifies the logic once implemented
- Performance will be better

The key insight is that **decoy generation is a property of the peptide sequence, not the modification state**.