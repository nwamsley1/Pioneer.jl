# Fix Plan: Invalid Characters in Altimeter Peptide Sequences

## Problem Summary

The Altimeter API is rejecting peptide sequences containing dash characters (`-`), causing the library build to fail:
```
ValueError: Invalid AA in requested peptide. AA:- Sequence:IL-ASSRLN- Peptide_index:952
```

Investigation revealed 4 sequences with dashes in the precursors table:
- `IL-ASSRLN-`
- `I-LRASLNS-`
- `ILSASNRL--`
- `R-NALLSIS-`

## Root Cause Analysis

### 1. Incomplete Character Filtering in FASTA Digestion
**Location:** `src/Routines/BuildSpecLib/fasta/fasta_digest.jl:177`

**Current Implementation:**
```julia
if (occursin("[H", peptide)) | (occursin("U", peptide)) | (occursin("O", peptide)) |
   (occursin("X", peptide)) | occursin("Z", peptide) | occursin("B", peptide)
    continue
end
```

**Problems:**
- Only checks for specific non-standard amino acids
- Does NOT filter dashes (`-`) or other unexpected characters
- Not robust to future edge cases
- Multiple `occursin` calls are inefficient

### 2. Missing Strip Operation for SplineCoefficientModel
**Location:** `src/Routines/BuildSpecLib/koina/koina_batch_prep.jl:149`

**Issue:**
- `RetentionTimeModel` uses `strip.()` to clean sequences (line 198)
- `SplineCoefficientModel` (Altimeter) does NOT use `strip.()`
- Whitespace or special characters are sent directly to API

## Proposed Solution

### Change 1: Robust Amino Acid Validation (Fast Implementation)

Replace the current character-by-character checking with a **Set-based validation** that only allows the 20 standard amino acids.

**File:** `src/Routines/BuildSpecLib/fasta/fasta_digest.jl`

**Implementation Strategy:**

1. **Define a global constant** (at module level, before the function):
```julia
# Set of valid amino acid characters for fast lookup (O(1) per character)
const VALID_AAS = Set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
```

2. **Replace line 177-179** with:
```julia
# Skip peptides containing non-standard amino acids
if !all(aa -> aa ∈ VALID_AAS, peptide)
    continue
end
```

**Why This Approach is Optimal:**

- **Fast:** Set lookup is O(1) per character
- **Memory efficient:** Set is created once at compile time
- **Robust:** Automatically rejects ANY unexpected character (-, X, U, O, Z, B, *, gaps, etc.)
- **Future-proof:** No need to maintain a list of invalid characters
- **Clean:** Single, readable condition

**Performance Comparison:**
- Old: 6 separate `occursin()` calls = O(6n) where n is peptide length
- New: Single iteration with Set lookups = O(n) where n is peptide length
- Net result: ~6x faster for filtering logic

### Change 2: Add Strip to SplineCoefficientModel

**File:** `src/Routines/BuildSpecLib/koina/koina_batch_prep.jl`
**Line:** 149

**Current:**
```julia
Dict(
    "name" => "peptide_sequences",
    "shape" => [nrow(batch_data), 1],
    "datatype" => "BYTES",
    "data" => batch_data.koina_sequence
)
```

**Updated:**
```julia
Dict(
    "name" => "peptide_sequences",
    "shape" => [nrow(batch_data), 1],
    "datatype" => "BYTES",
    "data" => strip.(batch_data.koina_sequence)
)
```

**Rationale:**
- Consistent with `RetentionTimeModel` implementation
- Defensive programming against whitespace issues
- Negligible performance cost

## Testing Strategy

### 1. Unit Test for Amino Acid Validation
Create test cases in `test/Routines/BuildSpecLib/fasta/`:
```julia
@testset "Invalid Amino Acid Filtering" begin
    # Test with valid sequence
    valid_peptides, _ = digest_sequence("PEPTIDE", r"[KR][^P]", 20, 5, 0)
    @test length(valid_peptides) > 0

    # Test with dash
    invalid_peptides, _ = digest_sequence("PEP-TIDE", r"[KR][^P]", 20, 5, 0)
    # After filtering, should be empty

    # Test with other invalid characters
    test_cases = ["PEPTUXDE", "PEP*TIDE", "PEPZIDE", "PEP TIDE"]
    for seq in test_cases
        peptides, _ = digest_sequence(seq, r"[KR][^P]", 20, 5, 0)
        # Verify filtered out
    end
end
```

### 2. Integration Test
```bash
# Re-run the failing library build command
julia --project=. build_spectral_library.jl params.json
```

**Expected Result:**
- No sequences with invalid characters reach Altimeter
- Build completes without `ValueError` from API

### 3. Performance Validation
```julia
# Benchmark the new vs old implementation
using BenchmarkTools
peptide = "PEPTIDEPEPTIDE"
@btime occursin("-", peptide)  # Old approach (one check)
@btime all(aa -> aa ∈ VALID_AAS, peptide)  # New approach (all checks)
```

## Implementation Order

1. **Add `VALID_AAS` constant** to `fasta_digest.jl` (after imports, before functions)
2. **Update filtering logic** in `digest_fasta()` function (line 177)
3. **Add `strip.()`** to SplineCoefficientModel batch prep (line 149)
4. **Add unit tests** for the new validation
5. **Run integration test** with problematic FASTA file

## Benefits

✅ **Robustness:** Automatically handles all current and future invalid characters
✅ **Performance:** ~6x faster than multiple `occursin()` calls
✅ **Maintainability:** No need to update invalid character list
✅ **Clarity:** Single, self-documenting validation condition
✅ **Consistency:** Aligns with Koina API expectations (20 standard AAs only)

## Risks and Mitigations

**Risk:** Some users might have FASTA files with non-standard amino acids they want to keep
**Mitigation:** Add logging to report filtered peptides:
```julia
if !all(aa -> aa ∈ VALID_AAS, peptide)
    @debug_l2 "Filtered peptide with non-standard amino acids: $peptide"
    continue
end
```

## Documentation Updates

After implementation, update:
1. FASTA digestion function docstring to mention filtering
2. Build log to show count of filtered peptides
3. User documentation about supported amino acids
