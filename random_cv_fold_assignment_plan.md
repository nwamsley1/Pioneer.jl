# Plan: Replace CV Fold Assignment with Simple Random Assignment

## Overview

This plan describes the minimal intervention to replace the current peptide-based CV fold assignment with simple random assignment per unique protein group.

**Goal**: Each unique protein name gets randomly assigned to a CV fold (0 or 1), with the assignment being consistent across all instances of that protein.

---

## Current Implementation Summary

The current system has 3 key components:

1. **`build_protein_cv_fold_mapping()`** (Lines 1700-1756)
   - Scans PSM files to find highest-scoring peptide per protein
   - Assigns protein to CV fold of that peptide
   - Returns: `Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}`

2. **`assign_protein_group_cv_folds!()`** (Lines 1771-1794)
   - Takes the mapping and applies it to protein groups DataFrame
   - Adds `:cv_fold` column to DataFrame

3. **`perform_probit_analysis_multifold()`** (Lines 1888-2014)
   - Main orchestrator
   - Calls the above functions if `protein_to_cv_fold` is not provided
   - Uses the mapping for training/testing

---

## Proposed Changes

### Strategy: Replace `build_protein_cv_fold_mapping()` with Random Version

**Minimal Intervention**:
1. Comment out the existing `build_protein_cv_fold_mapping()` function
2. Create new `build_random_protein_cv_fold_mapping()` function
3. Update the single call site in `perform_probit_analysis_multifold()`

This approach:
- ✅ Minimal code changes (only 2 locations modified)
- ✅ Preserves all other logic (assignment, model training, etc.)
- ✅ Easy to revert (just uncomment and change function name back)
- ✅ No changes to function signatures or data structures

---

## Detailed Implementation Plan

### Step 1: Comment Out Original Function (Lines 1682-1756)

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`

**Action**: Wrap the entire function in a block comment

**Lines to modify**: 1682-1756

```julia
#=
"""
    build_protein_cv_fold_mapping(psm_paths::Vector{String}, precursors::LibraryPrecursors)
    -> Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}

Build a mapping from protein names to their CV fold assignments based on highest-scoring peptides.

[ORIGINAL IMPLEMENTATION - TEMPORARILY DISABLED]
[Full docstring and function code here...]
"""
function build_protein_cv_fold_mapping(
    psm_paths::Vector{String},
    precursors::LibraryPrecursors
)
    # ... original implementation ...
end
=#
```

**Note**: Use `#= ... =#` block comment to preserve the entire function including docstring

---

### Step 2: Create New Random Assignment Function (Insert after Line 1756)

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`

**Location**: Insert immediately after the commented-out function (after line 1756)

**New Function**:

```julia
"""
    build_protein_cv_fold_mapping(psm_paths::Vector{String}, precursors::LibraryPrecursors)
    -> Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}

Build a mapping from protein names to CV fold assignments using RANDOM assignment.

**TEMPORARY IMPLEMENTATION**: Replaces peptide-based assignment for testing.
Each unique protein is randomly assigned to fold 0 or 1 with equal probability.

# Arguments
- `psm_paths`: Vector of paths to PSM files (used to discover proteins)
- `precursors`: Library precursors (not used in random version, kept for signature compatibility)

# Returns
- Dictionary mapping protein_name to named tuple with:
  - `best_score`: Set to 1.0f0 (dummy value, not used in CV)
  - `cv_fold`: Randomly assigned to 0 or 1

# Process
1. Scans PSM files to collect all unique protein names
2. Randomly assigns each protein to fold 0 or 1
3. Uses fixed random seed (1234) for reproducibility

# Notes
- **This is a simplified replacement for the peptide-based assignment**
- All instances of the same protein get the same CV fold
- Random seed ensures reproducible results across runs
"""
function build_protein_cv_fold_mapping(
    psm_paths::Vector{String},
    precursors::LibraryPrecursors
)
    using Random

    # Fixed seed for reproducibility
    rng = MersenneTwister(1234)

    # Create mapping: protein_name -> (best_score=1.0, cv_fold=random)
    protein_to_cv_fold = Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}()

    # Collect all unique protein names across all PSM files
    all_proteins = Set{String}()

    for psm_path in psm_paths
        # Skip if PSM file doesn't exist
        if !isfile(psm_path)
            @user_warn "PSM file not found: $psm_path"
            continue
        end

        # Load PSM data
        psms = DataFrame(Arrow.Table(psm_path))

        # Verify required column exists
        if !hasproperty(psms, :inferred_protein_group)
            @user_warn "PSM file missing :inferred_protein_group column: $psm_path"
            continue
        end

        # Filter for valid PSMs with protein group assignments
        psms = filter(row -> !ismissing(row.inferred_protein_group), psms)

        # Add unique proteins to the set
        for protein_name in psms.inferred_protein_group
            push!(all_proteins, protein_name)
        end
    end

    # Convert to sorted vector for deterministic iteration
    protein_names = sort(collect(all_proteins))

    # Randomly assign each protein to fold 0 or 1
    n_proteins = length(protein_names)

    @user_info "Random CV fold assignment: assigning $n_proteins unique proteins to 2 folds"

    for protein_name in protein_names
        # Randomly choose fold 0 or 1
        cv_fold = rand(rng, UInt8[0, 1])

        # Create entry with dummy best_score
        value = (best_score = 1.0f0, cv_fold = cv_fold)
        insert!(protein_to_cv_fold, protein_name, value)
    end

    # Report fold distribution
    fold_0_count = count(p -> p.cv_fold == 0, values(protein_to_cv_fold))
    fold_1_count = count(p -> p.cv_fold == 1, values(protein_to_cv_fold))

    @user_info "Random CV fold distribution: Fold 0: $fold_0_count proteins, Fold 1: $fold_1_count proteins"

    return protein_to_cv_fold
end
```

---

### Step 3: No Changes Needed to Other Functions

**Important**: The following functions require **NO CHANGES**:

1. ✅ **`assign_protein_group_cv_folds!()` (Lines 1771-1794)**
   - Works identically with random mapping
   - Still applies mapping to DataFrame

2. ✅ **`apply_probit_scores_multifold!()` (Lines 1810-1852)**
   - Still uses `protein_to_cv_fold` dictionary
   - No changes needed

3. ✅ **`perform_probit_analysis_multifold()` (Lines 1888-2014)**
   - Still calls `build_protein_cv_fold_mapping()` at line 1909
   - Function name unchanged, so no modifications needed
   - The new random implementation will be called automatically

---

## Code Changes Summary

### File: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`

#### Change 1: Comment Out Original (Lines 1682-1756)
```julia
# Add opening block comment at line 1682
#=

# Add closing block comment at line 1756
=#
```

#### Change 2: Insert New Function (After Line 1756)
```julia
# Insert the complete new function shown in Step 2 above
```

**Total lines modified**: ~2 (add block comments)
**Total lines added**: ~85 (new function with documentation)
**Total lines deleted**: 0 (original function preserved in comments)

---

## Testing the Changes

### Before Running

1. **Verify the changes**:
   ```julia
   # Check that original function is commented out
   # Check that new function is present
   # Check function signature matches exactly
   ```

2. **Expected behavior**:
   - Each unique protein gets random fold assignment (0 or 1)
   - Same protein always gets same fold across all files
   - ~50% of proteins in each fold (random distribution)
   - All downstream code works identically

### Validation

Run the ScoringSearch pipeline and verify:

1. ✅ Log messages show random assignment:
   ```
   Random CV fold assignment: assigning N unique proteins to 2 folds
   Random CV fold distribution: Fold 0: X proteins, Fold 1: Y proteins
   ```

2. ✅ Probit training proceeds normally for both folds

3. ✅ Protein group files updated correctly

4. ✅ Final results have pg_score, pg_qval, etc.

---

## Key Properties Preserved

1. **Consistent Assignment**: Each protein name → same CV fold globally
2. **Two-Fold CV**: Still uses folds 0 and 1
3. **Dictionary Structure**: Same data structure maintained
4. **Downstream Compatibility**: All other functions work unchanged

---

## Differences from Original

| Aspect | Original (Peptide-Based) | New (Random) |
|--------|-------------------------|--------------|
| Assignment basis | Highest-scoring peptide's CV fold | Random (50/50) |
| `best_score` field | Actual peptide probability | Dummy value (1.0) |
| Fold distribution | Depends on peptide scores | ~50% each fold |
| Determinism | Based on data | Based on random seed |
| Library dependency | Uses `precursors.pid_to_cv_fold` | Ignores precursors |

---

## Reverting the Changes

To restore original behavior:

1. Remove the new random function
2. Uncomment the original function (remove `#=` and `=#`)
3. Done - no other changes needed

---

## Alternative Approaches (Not Recommended)

### Alternative 1: Modify Function In-Place
**Why not**: Loses original implementation, harder to revert

### Alternative 2: Add Parameter to Choose Method
**Why not**: More complex, unnecessary for temporary change

### Alternative 3: Create Separate Function Name
**Why not**: Requires changing call site, more invasive

---

## Implementation Checklist

- [ ] Backup current file
- [ ] Add opening block comment at line 1682: `#=`
- [ ] Add closing block comment at line 1756: `=#`
- [ ] Insert new function after line 1756
- [ ] Verify function signature matches exactly
- [ ] Check imports (Random, MersenneTwister)
- [ ] Test compilation
- [ ] Run integration test
- [ ] Verify log messages
- [ ] Check fold distribution is ~50/50
- [ ] Validate results

---

## Notes and Caveats

1. **Random Seed**: Using fixed seed (1234) ensures reproducibility
   - Same protein names → same fold assignments across runs
   - Change seed to get different random assignments

2. **Best Score Field**: Set to dummy value 1.0
   - Not used in CV training/scoring
   - Only used for determining fold in original implementation
   - Safe to ignore in random version

3. **Library Precursors**: Parameter kept for signature compatibility
   - Not used in random implementation
   - Could be removed but would require changing call site

4. **Performance**: Random version is faster
   - No need to scan precursor indices
   - No need to track best scores
   - Only needs to collect unique protein names

5. **Statistical Properties**:
   - Original: Fold assignment correlated with peptide detectability
   - Random: Fold assignment independent of any protein properties
   - Random version may give **better** CV estimates (less correlation)

---

**Document Created**: 2025-01-18
**Purpose**: Temporary modification for testing random CV fold assignment
**Reversibility**: High (simple uncomment to restore)
