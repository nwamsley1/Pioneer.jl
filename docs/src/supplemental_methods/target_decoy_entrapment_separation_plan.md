# Plan: Separate Target/Decoy and Entrapment Groups Before Inference

## Current Problem

The current implementation handles target/decoy and entrapment group separation **within** the inference algorithm:
- Lines 220-228: Group proteins by `(is_target, entrap_id)` when merging indistinguishable proteins
- Lines 240-246: Filter peptides by matching target/decoy and entrap_id when assigning
- Lines 269-298: Similar grouping logic in Case F (no unique peptides)
- Lines 323-332: Group by `(pep_set_hash, is_target, entrap_id)` in merge-first greedy

**Issues:**
1. Complexity scattered throughout algorithm
2. Easy to miss edge cases
3. Inefficient - checking metadata repeatedly
4. Harder to reason about correctness

## Proposed Solution

**Pre-partition input by (is_target, entrap_id), run inference separately, merge results**

### Architecture

```
infer_proteins(proteins, peptides)
  ├─ Step 1: Partition by (is_target, entrap_id)
  │    └─ Create subgroups: {(target, entrap1), (decoy, entrap1), (target, entrap2), ...}
  │
  ├─ Step 2: Run inference on each subgroup
  │    └─ infer_proteins_single_group(subgroup_proteins, subgroup_peptides)
  │         ├─ Build mappings
  │         ├─ Find connected components
  │         ├─ Phase 1: Select unique proteins
  │         ├─ Phase 2: Greedy with merge-first
  │         └─ Assign peptides
  │
  └─ Step 3: Merge results from all subgroups
       └─ Combine peptide_to_protein and use_for_quant dictionaries
```

### Benefits

1. **Simpler core algorithm**: No metadata checks during inference
2. **Cleaner separation of concerns**: Partitioning vs. inference logic
3. **Easier to test**: Can test single-group inference independently
4. **More efficient**: Metadata checks happen once upfront
5. **Easier to parallelize**: Each subgroup could be processed in parallel
6. **Better correctness**: Impossible to accidentally mix target/decoy or entrapments

---

## Detailed Implementation Plan

### Step 1: Create Helper Function `partition_by_metadata`

**Purpose**: Split input into subgroups by (is_target, entrap_id)

**Input:**
```julia
proteins::Vector{ProteinKey}
peptides::Vector{PeptideKey}
```

**Output:**
```julia
subgroups::Dictionary{Tuple{Bool, UInt8}, Tuple{Vector{ProteinKey}, Vector{PeptideKey}}}
# Key: (is_target, entrap_id)
# Value: (proteins_in_group, peptides_in_group)
```

**Algorithm:**
```julia
function partition_by_metadata(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)
    subgroups = Dictionary{Tuple{Bool, UInt8}, Tuple{Vector{ProteinKey}, Vector{PeptideKey}}}()

    for i in 1:length(peptides)
        peptide = peptides[i]
        protein = proteins[i]

        # Key based on peptide metadata (since peptide determines the group)
        key = (peptide.is_target, peptide.entrap_id)

        if !haskey(subgroups, key)
            insert!(subgroups, key, (ProteinKey[], PeptideKey[]))
        end

        # Add to appropriate subgroup
        push!(subgroups[key][1], protein)
        push!(subgroups[key][2], peptide)
    end

    return subgroups
end
```

**Key Decision**: Use peptide's metadata to determine group, not protein's
- Rationale: A peptide can only belong to one group, but a protein (in "A;B" format) might contain proteins from different groups
- When we split protein names by ";", we'll filter to only proteins matching the peptide's metadata

---

### Step 2: Extract Core Inference into `infer_proteins_single_group`

**Purpose**: Run inference on a single (is_target, entrap_id) group

**Input:**
```julia
proteins::Vector{ProteinKey}
peptides::Vector{PeptideKey}
# All entries have same is_target and entrap_id
```

**Output:**
```julia
InferenceResult(peptide_to_protein, use_for_quant)
```

**Changes to make:**
1. **Remove metadata grouping logic** (lines 220-228, 269-298, 323-332)
   - No need to group by `(is_target, entrap_id)` anymore
   - All proteins/peptides in the input already have matching metadata

2. **Simplify protein splitting** (lines 130-138)
   ```julia
   # OLD: Split and keep all proteins
   for protein_part in split(protein_key.name, ";")
       individual_protein = ProteinKey(protein_part, protein_key.is_target, protein_key.entrap_id)
       push!(...)
   end

   # NEW: Split and filter by metadata
   for protein_part in split(protein_key.name, ";")
       individual_protein = ProteinKey(protein_part, protein_key.is_target, protein_key.entrap_id)
       # Only include if metadata matches peptide
       if (individual_protein.is_target == peptide_key.is_target) &&
          (individual_protein.entrap_id == peptide_key.entrap_id)
           push!(...)
       end
   end
   ```

3. **Simplify merge logic** (lines 323-332)
   ```julia
   # OLD: Group by (pep_set_hash, is_target, entrap_id)
   group_key = (pep_set_hash, protein.is_target, protein.entrap_id)

   # NEW: Group only by pep_set_hash (metadata guaranteed to match)
   group_key = pep_set_hash
   ```

4. **Remove conditional peptide assignment** (lines 240-246)
   ```julia
   # OLD: Assign only if metadata matches
   if (is_target == peptide_key.is_target) && (entrap_id == peptide_key.entrap_id)
       insert!(peptide_to_protein, peptide_key, final_protein)
   end

   # NEW: Always assign (metadata guaranteed to match)
   insert!(peptide_to_protein, peptide_key, final_protein)
   ```

**Function signature:**
```julia
function infer_proteins_single_group(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::InferenceResult
    # ... existing logic with simplifications above
end
```

---

### Step 3: Refactor Main `infer_proteins` as Orchestrator

**Purpose**: Partition, delegate to single-group inference, merge results

**New structure:**
```julia
function infer_proteins(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::InferenceResult
    # Validate input
    if length(proteins) != length(peptides)
        throw(ArgumentError("proteins and peptides vectors must have the same length"))
    end

    # Handle empty input
    if isempty(proteins)
        return InferenceResult(
            Dictionary{PeptideKey, ProteinKey}(),
            Dictionary{PeptideKey, Bool}()
        )
    end

    # Step 1: Partition by (is_target, entrap_id)
    subgroups = partition_by_metadata(proteins, peptides)

    # Step 2: Run inference on each subgroup
    subgroup_results = Dictionary{Tuple{Bool, UInt8}, InferenceResult}()

    for (metadata_key, (subgroup_proteins, subgroup_peptides)) in pairs(subgroups)
        result = infer_proteins_single_group(subgroup_proteins, subgroup_peptides)
        insert!(subgroup_results, metadata_key, result)
    end

    # Step 3: Merge results
    merged_peptide_to_protein = Dictionary{PeptideKey, ProteinKey}()
    merged_use_for_quant = Dictionary{PeptideKey, Bool}()

    for (metadata_key, result) in pairs(subgroup_results)
        for (peptide, protein) in pairs(result.peptide_to_protein)
            insert!(merged_peptide_to_protein, peptide, protein)
        end
        for (peptide, use_quant) in pairs(result.use_for_quant)
            insert!(merged_use_for_quant, peptide, use_quant)
        end
    end

    return InferenceResult(merged_peptide_to_protein, merged_use_for_quant)
end
```

---

## Edge Cases and Considerations

### Edge Case 1: Protein Group Spanning Multiple Metadata Groups

**Scenario:**
```julia
proteins = [ProteinKey("A;B", true, UInt8(1))]
peptides = [PeptideKey("pep1", true, UInt8(1))]
```

Where:
- Protein A is target, entrap 1
- Protein B is decoy, entrap 1
- Combined as "A;B"

**Current behavior**: When splitting "A;B", both A and B are created with `(true, 1)` metadata from the ProteinKey

**New behavior**: When splitting "A;B", we filter to only proteins matching peptide metadata
- If peptide is `(true, 1)`, only protein A is kept
- If peptide is `(false, 1)`, only protein B is kept

**Action**: Update protein splitting logic to filter by metadata (Step 2, point 2 above)

### Edge Case 2: Empty Subgroups

**Scenario**: After partitioning, some metadata groups might be empty

**Solution**: The loop in Step 3 naturally handles this - empty groups produce empty results

### Edge Case 3: Single Metadata Group

**Scenario**: All peptides have the same `(is_target, entrap_id)`

**Solution**: Only one subgroup is created, inference runs once, no overhead

### Edge Case 4: Peptide with No Matching Proteins After Filtering

**Scenario**: Protein "A;B" but peptide metadata doesn't match either A or B after filtering

**Solution**: This should be impossible if input is well-formed (protein-peptide edges must be valid)
- Add assertion/validation in `partition_by_metadata` to catch malformed input

---

## Testing Strategy

### Unit Tests for `partition_by_metadata`

**Test 1: Single metadata group**
```julia
proteins = [ProteinKey("A"), ProteinKey("B")]
peptides = [PeptideKey("pep1", true, 1), PeptideKey("pep2", true, 1)]
# Expect: 1 subgroup with 2 proteins, 2 peptides
```

**Test 2: Multiple metadata groups**
```julia
proteins = [ProteinKey("A"), ProteinKey("B"), ProteinKey("C")]
peptides = [PeptideKey("pep1", true, 1), PeptideKey("pep2", false, 1), PeptideKey("pep3", true, 2)]
# Expect: 3 subgroups
```

**Test 3: Empty input**
```julia
proteins = []
peptides = []
# Expect: empty subgroups
```

### Integration Tests for `infer_proteins`

**Test 1: Target/decoy separation**
```julia
proteins = [
    ProteinKey("A", true, 1),
    ProteinKey("B", false, 1)
]
peptides = [
    PeptideKey("pep1", true, 1),
    PeptideKey("pep2", false, 1)
]
# Expect: pep1 → A, pep2 → B, never mixed
```

**Test 2: Entrapment group separation**
```julia
proteins = [
    ProteinKey("A", true, 1),
    ProteinKey("A", true, 2)
]
peptides = [
    PeptideKey("same_seq", true, 1),
    PeptideKey("same_seq", true, 2)
]
# Expect: Two separate protein groups, never mixed
```

**Test 3: Combined protein groups with metadata filtering**
```julia
proteins = [
    ProteinKey("A;B", true, 1),
    ProteinKey("A;B", false, 1)
]
peptides = [
    PeptideKey("pep1", true, 1),
    PeptideKey("pep2", false, 1)
]
# Expect: After splitting and filtering:
#   pep1 only sees protein A (target)
#   pep2 only sees protein B (decoy)
```

### Regression Tests

**Use existing test suite** from `test/UnitTests/test_protein_inference.jl`
- All existing tests should still pass
- Behavior should be identical to current implementation

---

## Migration Strategy

### Phase 1: Create New Functions Without Breaking Old Code

1. Add `partition_by_metadata` as a new helper function
2. Copy `infer_proteins` to `infer_proteins_single_group` (internal function)
3. Simplify `infer_proteins_single_group` to remove metadata logic
4. Don't modify public `infer_proteins` yet

### Phase 2: Test New Implementation

1. Create test file `test_metadata_separation.jl`
2. Write unit tests for `partition_by_metadata`
3. Write integration tests comparing old vs. new `infer_proteins`
4. Ensure identical behavior on all test cases

### Phase 3: Replace Implementation

1. Refactor public `infer_proteins` to use new architecture
2. Run full test suite
3. Fix any regressions
4. Remove old code

### Phase 4: Update Documentation

1. Update docstrings to reflect new architecture
2. Update `protein_inference.tex` if needed
3. Add examples showing target/decoy and entrapment separation

---

## Performance Considerations

### Potential Overhead

1. **Partitioning cost**: O(n) where n = number of peptides
   - Negligible compared to inference cost

2. **Multiple inference calls**: One per metadata group
   - In practice, typically 2-4 groups (target/decoy × 1-2 entrapments)
   - Each group is smaller, so total work may be similar or less

3. **Merging results**: O(n) to combine dictionaries
   - Negligible

### Potential Improvements

1. **Parallelization**: Each subgroup can be processed independently
   ```julia
   using Base.Threads

   subgroup_results = Dictionary{Tuple{Bool, UInt8}, InferenceResult}()
   @threads for (metadata_key, (subgroup_proteins, subgroup_peptides)) in collect(pairs(subgroups))
       result = infer_proteins_single_group(subgroup_proteins, subgroup_peptides)
       insert!(subgroup_results, metadata_key, result)
   end
   ```

2. **Reduced connected component complexity**: Smaller graphs per subgroup
   - Target and decoy proteins never in same component
   - Faster DFS traversal

---

## Code Structure Summary

### Before Refactoring
```
infer_proteins(proteins, peptides)
  ├─ Build mappings (with metadata checks)
  ├─ Find components (mixed target/decoy/entrapments)
  └─ For each component:
       ├─ Check indistinguishable (group by metadata)
       ├─ Phase 1: Unique proteins
       ├─ Phase 2: Greedy (group by metadata in merge)
       └─ Assign peptides (filter by metadata)
```

### After Refactoring
```
infer_proteins(proteins, peptides) [PUBLIC API]
  ├─ partition_by_metadata(proteins, peptides)
  ├─ For each (is_target, entrap_id) subgroup:
  │    └─ infer_proteins_single_group(subgroup_proteins, subgroup_peptides)
  │         ├─ Build mappings (no metadata checks needed)
  │         ├─ Find components (homogeneous metadata)
  │         └─ For each component:
  │              ├─ Check indistinguishable (no metadata grouping)
  │              ├─ Phase 1: Unique proteins
  │              ├─ Phase 2: Greedy (no metadata checks in merge)
  │              └─ Assign peptides (no metadata filtering)
  └─ Merge subgroup results
```

---

## Recommendation

**This refactoring should be done in a separate commit/PR after the merge-first fix**

**Rationale:**
1. The merge-first fix (current commit) already solves the immediate bug
2. This refactoring is a larger architectural change
3. Easier to review and test separately
4. If issues arise, easier to bisect which change caused them

**Next Steps:**
1. Merge current branch with merge-first fix
2. Create new branch `refactor/metadata-separation`
3. Implement this plan in phases
4. Thoroughly test before merging

This approach provides a cleaner, more maintainable codebase with better separation of concerns and potential for future optimizations like parallelization.
