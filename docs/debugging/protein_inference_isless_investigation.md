# Protein Inference `isless` Error Investigation

**Date**: 2025-10-11
**Branch**: `supplemental_methods`
**Commit**: `f9e6f08e` (refactor: remove use_for_quant dictionary)

## Problem Statement

During integration testing, the protein inference algorithm fails with:

```
MethodError: no method matching isless(::Pioneer.PeptideKey, ::Pioneer.PeptideKey)
```

**Location**: `src/utils/proteinInference.jl:328`

```julia
pep_set_hash = hash(sort(collect(remaining_peps)))
```

**Critical Mystery**: Unit tests pass (86/86) with the exact same code, but integration test fails immediately.

## Diagnostic Output from Integration Test

### Original Function Inputs
- **Total input pairs**: 69,914
- **Unique proteins**: 11,699
- **Unique peptides**: 69,914
- **Multi-protein assignments**: 3,880 (peptides with `;` in protein name)

### Sample Input (First 10 pairs)
```
[1] Protein: 'P50995' (target=true, entrap=0)
    Peptide: 'LLISLSQGNR' (target=true, entrap=0)
[2] Protein: 'P99999' (target=true, entrap=0)
    Peptide: 'GIIWGEDTLMEYLENPK' (target=true, entrap=0)
[3] Protein: 'P09874' (target=true, entrap=0)
    Peptide: 'TTNFAGILSQGLR' (target=true, entrap=0)
...
```

### Component State When Error Occurs
- **Candidate proteins**: 18
- **Remaining peptides**: 7

**First 3 candidate proteins**:
- O60814
- P33778
- P02293

**Remaining peptides**:
- KESYSIYVYK
- ESYSIYVYK
- ETYSSYIYK
- KESYSVYVYK
- AMGIMNSFVNDIFER

### Specific Failure Point

**Protein**: O60814
**Intersection size**: 3 peptides
**Failing peptides**: KESYSVYVYK, ESYSVYVYK, and one other

**Error**:
```
MethodError(isless, (Pioneer.PeptideKey("KESYSVYVYK", true, 0x00),
                     Pioneer.PeptideKey("ESYSVYVYK", true, 0x00)),
            0x00000000000069bf)
```

## Stack Trace Analysis

The error occurs at **TWO different line numbers** in the same function:

1. **Line 332** (original error in try block):
   ```
   infer_proteins(...) at proteinInference.jl:332
   ```

2. **Line 453** (error in catch block diagnostics):
   ```
   infer_proteins(...) at proteinInference.jl:453
   ```

This indicates the error occurred twice:
1. First during normal execution at line 328/332
2. Again when the diagnostic code tried to re-sort the same peptides at line 449

### Secondary Error: Scoping Bug

After the primary `isless` error, a secondary error occurred:

```
UndefVarError: `remaining_peps` not defined in `Pioneer`
```

**Location**: Line 459 (after the inner catch block)

**Cause**: The variable `remaining_peps` is declared inside the inner try block (line 438) but referenced outside it (line 455-459) when trying to print diagnostic information.

## Root Cause Analysis

### Primary Issue: Missing Comparison Methods

The `PeptideKey` struct lacks these required methods for sorting:

```julia
struct PeptideKey
    sequence::String
    is_target::Bool
    entrap_id::UInt8
end
```

**Missing methods**:
- `Base.isless(::PeptideKey, ::PeptideKey)` - Required for `sort()`
- `Base.:(==)(::PeptideKey, ::PeptideKey)` - For equality comparison
- `Base.hash(::PeptideKey, ::UInt)` - For Dictionary keys (already works implicitly)

Similarly, `ProteinKey` lacks these methods (used in other parts of the code at line 346).

### Why This Code Path Exists

The merge-first greedy set cover algorithm (lines 319-390) groups proteins by their remaining peptide sets using:

```julia
pep_set_hash = hash(sort(collect(remaining_peps)))
```

The `sort()` is needed to ensure identical peptide sets produce identical hashes regardless of set iteration order. This line was introduced in commit `86275bfd` ("fix: implement merge-first approach for indistinguishable proteins").

## The Unit Test vs Integration Mystery

### Unit Test Status
- **Result**: 86/86 passing
- **Commit**: f9e6f08e (same code with `hash(sort(...))`)
- **Scale**: Typically 4-8 input pairs, 2-4 unique proteins

### Why Tests Pass

Several possibilities:

1. **Different Code Path**: Unit tests might not enter the greedy set cover loop if all peptides are resolved in Phase 1 (unique peptide selection)

2. **Julia Sort Optimization**: For very small collections (2-3 elements), Julia's `sort()` might use a different algorithm that doesn't require `isless` comparisons, or the compiler optimizes it away

3. **Set Iteration Order**: With small sets, the hash of the unsorted collection might coincidentally match across runs, avoiding the need for sorting

4. **Test Data Characteristics**: The carefully constructed test cases might not trigger the specific component topology that requires merging

### Verification Experiment

Testing with a minimal struct:
```julia
struct TestKey; x::String; end
s = Set([TestKey("a"), TestKey("b")])
sort(collect(s))  # MethodError
```

This confirms that structs without `isless` CANNOT be sorted, regardless of size.

### Conclusion

**The unit tests DO exercise the merge-first code** (Cases I, J, K test this explicitly), but somehow avoid triggering the sort operation. This suggests the tests take a different code path through the algorithm, possibly due to:
- Specific component structure
- Number of candidate proteins
- Distribution of remaining peptides

The integration test with 11,699 proteins and complex protein-peptide topology hits a case that absolutely requires sorting.

## Solution

### Fix 1: Add Comparison Methods to PeptideKey

In `src/structs/protein_inference_types.jl`, after the `PeptideKey` struct definition (~line 57):

```julia
# Comparison methods for PeptideKey (required for sorting)
Base.isless(a::PeptideKey, b::PeptideKey) =
    (a.sequence, a.is_target, a.entrap_id) < (b.sequence, b.is_target, b.entrap_id)

Base.:(==)(a::PeptideKey, b::PeptideKey) =
    a.sequence == b.sequence && a.is_target == b.is_target && a.entrap_id == b.entrap_id

Base.hash(k::PeptideKey, h::UInt) =
    hash((k.sequence, k.is_target, k.entrap_id), h)
```

### Fix 2: Add Comparison Methods to ProteinKey

In `src/structs/protein_inference_types.jl`, after the `ProteinKey` struct definition (~line 41):

```julia
# Comparison methods for ProteinKey (required for sorting)
Base.isless(a::ProteinKey, b::ProteinKey) =
    (a.name, a.is_target, a.entrap_id) < (b.name, b.is_target, b.entrap_id)

Base.:(==)(a::ProteinKey, b::ProteinKey) =
    a.name == b.name && a.is_target == b.is_target && a.entrap_id == b.entrap_id

Base.hash(k::ProteinKey, h::UInt) =
    hash((k.name, k.is_target, k.entrap_id), h)
```

**Note**: `ProteinKey` comparison is used at line 346: `protein_names = sort([p.name for p in proteins])`

### Fix 3: Fix Scoping Bug in Diagnostic Code

In `src/utils/proteinInference.jl`, around line 435-465, declare `remaining_peps` outside the inner try block:

```julia
# Before (broken):
for (idx, protein) in enumerate(candidate_proteins)
    println("\nTesting protein $idx: $(protein.name)")
    try
        remaining_peps = intersect(...)  # Scoped inside try
        ...
    catch inner_e
        println("  ✗ ERROR at this protein!")
        println("  Error: ", inner_e)
        println("\n  Peptides for this protein:")
        for (pi, pep) in enumerate(collect(remaining_peps))  # ERROR: out of scope!
            ...
        end
    end
end

# After (fixed):
for (idx, protein) in enumerate(candidate_proteins)
    println("\nTesting protein $idx: $(protein.name)")
    remaining_peps = Set{PeptideKey}()  # Declare outside
    try
        remaining_peps = intersect(...)  # Assignment, not declaration
        ...
    catch inner_e
        println("  ✗ ERROR at this protein!")
        println("  Error: ", inner_e)
        println("\n  Peptides for this protein:")
        for (pi, pep) in enumerate(collect(remaining_peps))  # Now in scope!
            ...
        end
    end
end
```

## Testing Plan

1. **Add comparison methods** to both structs
2. **Fix scoping bug** in diagnostic code
3. **Run unit tests**: Verify 86/86 still pass
4. **Run integration test**: Should progress past protein inference
5. **Verify diagnostics**: Remove comparison methods temporarily to test diagnostic output works

## Lessons Learned

1. **Struct Comparison Requirements**: Any struct used with `sort()` must implement `Base.isless`
2. **Test Coverage Limitations**: Unit tests can pass even when integration fails due to:
   - Different code paths through complex algorithms
   - Scale differences (small vs large data)
   - Specific data characteristics
3. **Julia Scoping**: Variables declared inside try blocks are not accessible in catch blocks
4. **Diagnostic Code Quality**: Error diagnostic code can itself contain bugs that mask the real issue

## References

- Commit introducing merge-first: `86275bfd`
- Commit removing use_for_quant: `f9e6f08e`
- Julia Base.Order documentation: https://docs.julialang.org/en/v1/base/sort/#Base.Order.Ordering
- Test file: `test/UnitTests/test_protein_inference.jl`
