# Protein Inference Wrapper Implementation Plan

**Date**: 2025-10-13
**Branch**: `supplemental_methods`
**Status**: AWAITING APPROVAL

## Objective

Refactor the protein inference API to use automatic grouping by entrapment and target/decoy status, ensuring that inference is performed separately for each unique combination of `(is_target, entrap_id)`.

## Motivation

Currently, `infer_proteins()` accepts two parallel vectors and performs inference on the entire dataset. However, protein inference should be performed independently for each combination of:
- **Target/Decoy status** (`is_target: Bool`)
- **Entrapment group** (`entrap_id: UInt8`)

This ensures:
1. Target and decoy proteins are never merged together
2. Different entrapment groups are kept separate (important for FDR calibration)
3. Each subset gets correct parsimony-based inference

## Proposed API Refactoring

### Rename and Change Signatures

**Step 1**: Rename current function to internal implementation:
```julia
# Current function becomes internal implementation
function _infer_proteins_single_group(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::InferenceResult
    # ... existing implementation (no logic changes)
end
```

**Step 2**: Create new public API with automatic grouping:
```julia
"""
    infer_proteins(proteins::Vector{ProteinKey}, peptides::Vector{PeptideKey})::InferenceResult

Perform protein inference with automatic grouping by entrapment and target/decoy status.

This function groups input protein-peptide pairs by (is_target, entrap_id) and
performs protein inference separately for each group. This ensures:
- Target and decoy proteins are never merged together
- Different entrapment groups are kept independent
- Correct parsimony-based inference within each population

# Arguments
- `proteins`: Vector of ProteinKey representing protein assignments
- `peptides`: Vector of PeptideKey representing peptide assignments (must be same length)

# Returns
- `InferenceResult`: Combined results from all groups, mapping unique peptides to proteins

# Example
```julia
proteins = [
    ProteinKey("P1", true, 0x00),
    ProteinKey("P2", true, 0x00),
    ProteinKey("REV_P1", false, 0x00),
]
peptides = [
    PeptideKey("PEP1", true, 0x00),
    PeptideKey("PEP2", true, 0x00),
    PeptideKey("REV_PEP1", false, 0x00),
]
result = infer_proteins(proteins, peptides)
```

# Implementation Note
Internally groups by (protein.is_target, protein.entrap_id) and calls
`_infer_proteins_single_group()` for each group.
"""
function infer_proteins(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::InferenceResult
```

### Relationship Between Functions

- **New Internal**: `_infer_proteins_single_group(proteins, peptides)`
  - Low-level implementation (existing algorithm)
  - Assumes all inputs belong to the same group
  - No grouping logic
  - Not exported (internal use only)

- **New Public**: `infer_proteins(proteins, peptides)`
  - High-level API (same signature as current)
  - Automatically groups by (is_target, entrap_id)
  - Calls `_infer_proteins_single_group()` for each group
  - Combines results
  - This is the function users call

## Implementation Plan

### Step 1: Rename Existing Function

Rename the current `infer_proteins()` to `_infer_proteins_single_group()`:
- This becomes an internal implementation detail
- No logic changes, just a rename
- Not exported from the module

### Step 2: Create Grouping Helper Function

```julia
function _group_by_population(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::Dictionary{Tuple{Bool, UInt8}, Tuple{Vector{ProteinKey}, Vector{PeptideKey}}}

    groups = Dictionary{Tuple{Bool, UInt8}, Tuple{Vector{ProteinKey}, Vector{PeptideKey}}}()

    for i in 1:length(proteins)
        prot = proteins[i]
        pep = peptides[i]

        # Group key is (protein.is_target, protein.entrap_id)
        # We use protein's status since that determines the grouping
        group_key = (prot.is_target, prot.entrap_id)

        if !haskey(groups, group_key)
            insert!(groups, group_key, (ProteinKey[], PeptideKey[]))
        end

        push!(groups[group_key][1], prot)
        push!(groups[group_key][2], pep)
    end

    return groups
end
```

**Design Decision**: Use `protein.is_target` and `protein.entrap_id` as the grouping key because:
- Proteins define the inference groups
- Peptides inherit their group from their parent protein
- In valid data, `protein.is_target == peptide.is_target` (could add validation)

### Step 3: Implement New Public API

```julia
function infer_proteins(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)::InferenceResult
    # Validate input lengths match
    if length(proteins) != length(peptides)
        throw(ArgumentError("proteins and peptides vectors must have the same length"))
    end

    # Step 1: Group by (is_target, entrap_id)
    groups = _group_by_population(proteins, peptides)

    # Step 2: Perform inference on each group
    combined_result = Dictionary{PeptideKey, ProteinKey}()

    for ((is_target, entrap_id), (group_proteins, group_peptides)) in pairs(groups)
        # Perform inference on this group using internal function
        group_result = _infer_proteins_single_group(group_proteins, group_peptides)

        # Merge results into combined dictionary
        for (pep, prot) in pairs(group_result.peptide_to_protein)
            # Sanity check: ensure no conflicts (same peptide assigned to different proteins)
            if haskey(combined_result, pep)
                if combined_result[pep] != prot
                    @warn "Peptide $(pep.sequence) assigned to multiple proteins across groups"
                end
            else
                insert!(combined_result, pep, prot)
            end
        end
    end

    return InferenceResult(combined_result)
end
```

### Step 4: Add Validation (Optional)

```julia
function _validate_grouped_pairs(
    proteins::Vector{ProteinKey},
    peptides::Vector{PeptideKey}
)
    """Validate that protein and peptide status fields match."""
    for i in 1:length(proteins)
        prot = proteins[i]
        pep = peptides[i]

        if prot.is_target != pep.is_target
            @warn "Mismatched target status: protein $(prot.name) ($(prot.is_target)) " *
                  "paired with peptide $(pep.sequence) ($(pep.is_target))"
        end
        if prot.entrap_id != pep.entrap_id
            @warn "Mismatched entrap_id: protein $(prot.name) ($(prot.entrap_id)) " *
                  "paired with peptide $(pep.sequence) ($(pep.entrap_id))"
        end
    end
end
```

### Step 5: Update Documentation

Update `src/utils/proteinInference.jl` docstring to reflect new behavior:

```julia
"""
    infer_proteins(proteins::Vector{ProteinKey}, peptides::Vector{PeptideKey})::InferenceResult

Perform protein inference with automatic grouping by entrapment and target/decoy status.

Automatically groups inputs by (protein.is_target, protein.entrap_id) and performs
inference separately for each group. This ensures target/decoy separation and
independent processing of entrapment groups.

# Implementation Note
Internally calls `_infer_proteins_single_group()` for each group and combines results.
"""
```

### Step 6: Add Tests

Add comprehensive tests in `test/UnitTests/test_protein_inference.jl`:

```julia
@testset "infer_proteins automatic grouping" begin
    @testset "Single group (backward compatibility)" begin
        # Test with only one group - should work as before
        # All proteins have same (is_target, entrap_id)
    end

    @testset "Target and decoy separation" begin
        # Test that targets and decoys are processed separately
        # Mix of is_target=true and is_target=false
    end

    @testset "Multiple entrapment groups" begin
        # Test with entrap_id 0, 1, 2
        # Verify each group is processed independently
    end

    @testset "Combined target/decoy and entrapment groups" begin
        # Test all combinations: (true, 0), (false, 0), (true, 1), (false, 1)
    end

    @testset "Cross-group peptide handling" begin
        # Same peptide sequence in different groups
        # Should be allowed (different PeptideKey instances)
    end

    @testset "Empty input" begin
        # Test with empty vectors
    end
end
```

## File Locations

- **Implementation**: `src/utils/proteinInference.jl`
  - Rename existing `infer_proteins()` to `_infer_proteins_single_group()`
  - Add `_group_by_population()` helper
  - Add new `infer_proteins()` wrapper
- **Tests**: `test/UnitTests/test_protein_inference.jl`
  - Existing tests continue to work (single group = no change in behavior)
  - Add new testset for multi-group scenarios
- **Documentation**: Update docstrings in `proteinInference.jl`

## Edge Cases to Handle

1. **Empty input**: `proteins = [], peptides = []` → return empty InferenceResult
2. **Single group**: Should behave identically to current `infer_proteins()` behavior
   - All existing tests should pass without modification
3. **Cross-group peptides**: Same peptide sequence appearing in different groups
   - Expected: Each group infers independently, peptide may map to different proteins
   - Validation: Warn if same PeptideKey (with identical status/entrap) assigned to different proteins
4. **Mismatched status**: Protein and peptide have different `is_target` or `entrap_id`
   - Decision: Use protein's status for grouping (peptide inherits)
   - Add validation warnings (optional)
5. **Large datasets**: With many groups, consider memory usage
   - Current approach processes sequentially (low memory)
   - Could parallelize if needed (future optimization)

## Performance Considerations

- **Time complexity**: O(G × N) where G is number of groups, N is total pairs
  - Grouping: O(N) single pass
  - Inference per group: O(n_g × log n_g) for group of size n_g
  - Total: Sum of group costs ≈ O(N log N) if groups are balanced

- **Space complexity**: O(N)
  - Groups dictionary: O(N) total across all groups
  - Combined result: O(unique peptides)

- **Memory profile**: Sequential processing keeps peak memory low

## Migration Strategy

This is a **transparent enhancement**:
- External API remains unchanged: `infer_proteins(proteins, peptides)`
- Internal implementation changes to add automatic grouping
- Existing code continues to work without modification
- Single-group inputs behave identically to before

**Backward Compatibility**:
- All existing tests should pass without modification
- Existing callers get automatic grouping for free
- No breaking changes to API

## Testing Strategy

1. **Existing unit tests**: Should pass without modification
   - All current tests use single group (same is_target, entrap_id)
   - Behavior should be identical
2. **New unit tests**: Test multi-group scenarios
   - Target/decoy separation
   - Multiple entrapment groups
   - Cross-group peptides
3. **Integration test**: Run full SearchDIA pipeline, verify:
   - Target and decoy proteins are separate
   - Entrapment groups are independent
   - Results match expected protein groups

## Open Questions

1. **Validation strictness**: Should mismatched status throw error or just warn?
   - Recommendation: Warn only (permissive), since protein status is authoritative
   - Could make this optional for performance

2. **Conflict handling**: What if same PeptideKey assigned to different proteins in different groups?
   - Current plan: Warn but allow (each group is independent)
   - This shouldn't happen with proper data (same PeptideKey implies same status/entrap)
   - Alternative: Use first assignment, ignore others silently

3. **Return type**: Should we return metadata about which groups were processed?
   - Current: Return standard InferenceResult (no breaking changes)
   - Alternative: Add `groups_processed::Set{Tuple{Bool, UInt8}}` field
   - Decision: Keep simple for now, can enhance later if needed

4. **Internal function naming**: Is `_infer_proteins_single_group()` clear?
   - Alternative: `_infer_proteins_impl()`
   - Alternative: `_infer_proteins_core()`
   - Current choice emphasizes single-group assumption

## Success Criteria

Implementation is complete when:
- [ ] Current `infer_proteins()` renamed to `_infer_proteins_single_group()`
- [ ] `_group_by_population()` helper function implemented
- [ ] New `infer_proteins()` wrapper implemented
- [ ] Grouping logic correctly separates by (is_target, entrap_id)
- [ ] All existing unit tests pass without modification
- [ ] New multi-group unit tests added and passing (≥5 test cases)
- [ ] Documentation updated with new behavior
- [ ] Integration test with real data succeeds

## Timeline Estimate

- Implementation: 30-45 minutes
- Testing: 30-45 minutes
- Documentation: 15 minutes
- **Total**: ~1.5-2 hours

---

## Approval Checklist

Please review and approve:
- [ ] API refactoring approach (rename existing, create new wrapper with same name)
- [ ] Function signatures (parallel vectors in/out, same as current)
- [ ] Grouping strategy (use protein's is_target/entrap_id)
- [ ] Conflict handling (warn but allow cross-group peptide differences)
- [ ] Validation approach (optional warnings, not errors)
- [ ] Testing strategy (existing tests must pass)
- [ ] Backward compatibility (transparent enhancement, no breaking changes)

**Next Steps After Approval**:
1. Rename `infer_proteins()` to `_infer_proteins_single_group()`
2. Implement `_group_by_population()` helper
3. Implement new `infer_proteins()` wrapper
4. Verify existing tests pass
5. Add new multi-group tests
6. Update documentation
7. Run integration test
8. Commit changes
