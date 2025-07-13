# Implementation Plan: Target/Decoy Handling in Protein EFDR Analysis

## Overview

Currently, the protein entrapment pairing groups only by protein name, which incorrectly pairs target and decoy proteins together. This needs to be updated to consider the target/decoy status when forming entrapment pairs.

## Current Issues

1. **Incorrect Pairing**: Target and decoy proteins with the same name are being paired together
2. **No Decoy Filtering**: Decoy proteins should be filtered out before EFDR calculation
3. **Grouping Logic**: The groupby operation needs to include target/decoy status

## Proposed Changes

### 1. Update `assign_protein_entrapment_pairs!` Function

**Location**: `src/core/protein_entrapment_pairing.jl`

**Current grouping**:
```julia
grouped = groupby(df, :protein)
```

**Proposed grouping**:
```julia
# Check for target column
if !hasproperty(df, :target)
    error("DataFrame must have :target column to distinguish target/decoy proteins")
end

# Group by protein name AND target status
grouped = groupby(df, [:protein, :target])
```

### 2. Update `run_protein_efdr_analysis` Function

**Location**: `src/api.jl`

**Add decoy filtering** (lines ~475):
```julia
# Check for target column and filter
original_rows = nrow(protein_results)
if hasproperty(protein_results, :target)
    non_target_rows = sum(.!protein_results.target)
    if non_target_rows > 0
        verbose && @warn "Filtering out $non_target_rows non-target (decoy) rows before EFDR calculation"
        filter!(x -> x.target, protein_results)
    end
else
    verbose && @warn "No 'target' column found. Assuming all rows are targets."
end
```

Note: This filtering already exists! Just need to ensure it's working properly.

### 3. Update `add_original_target_protein_scores!` Function

**Location**: `src/core/protein_scoring_v2.jl`

**Current key generation**:
```julia
key = (row[file_col], row.protein)
```

**Proposed key generation**:
```julia
# Include target status in the key if available
key = if hasproperty(protein_results, :target)
    (row[file_col], row.protein, row.target)
else
    (row[file_col], row.protein)
end
```

### 4. Update Documentation

#### Required Column Updates

Add `:target` to required columns documentation:
- In `PROTEIN_ANALYSIS_WORKFLOW.md`
- In function docstrings
- In examples

#### Pairing Logic Documentation

Update to reflect that pairing considers both protein name and target/decoy status.

## Implementation Steps

### Phase 1: Core Function Updates
1. Update `assign_protein_entrapment_pairs!` to group by `[:protein, :target]`
2. Add validation to ensure target column exists
3. Update tests to include target/decoy scenarios

### Phase 2: Score Propagation Updates
1. Modify `add_original_target_protein_scores!` to consider target status
2. Ensure dictionary keys include target status when available
3. Update score lookup logic

### Phase 3: Documentation Updates
1. Update `PROTEIN_ANALYSIS_WORKFLOW.md` with new requirements
2. Update function docstrings
3. Add examples showing target/decoy handling

### Phase 4: Testing
1. Create test cases with mixed target/decoy proteins
2. Verify decoys are filtered before EFDR calculation
3. Ensure targets and decoys are never paired together

## Code Changes Summary

### 1. `protein_entrapment_pairing.jl`

```julia
function assign_protein_entrapment_pairs!(df::DataFrame)
    # Check required columns
    required_cols = [:protein, :entrap_id, :target]  # Added :target
    missing_cols = [col for col in required_cols if !hasproperty(df, col)]
    if !isempty(missing_cols)
        error("DataFrame missing required columns: $missing_cols")
    end
    
    # ... initialization code ...
    
    # Group by protein name AND target status
    grouped = groupby(df, [:protein, :target])
    
    # ... rest of the function remains the same ...
end
```

### 2. `protein_scoring_v2.jl`

```julia
function add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)
    # ... validation code ...
    
    # Build dictionary with target status in key
    if hasproperty(protein_results, :target)
        if file_col == :ms_file_idx
            protein_to_target = Dictionary{Tuple{Int, String, Bool}, Float32}()
        else
            protein_to_target = Dictionary{Tuple{String, String, Bool}, Float32}()
        end
    else
        # Fallback to original behavior if no target column
        if file_col == :ms_file_idx
            protein_to_target = Dictionary{Tuple{Int, String}, Float32}()
        else
            protein_to_target = Dictionary{Tuple{String, String}, Float32}()
        end
    end
    
    # ... rest of function with updated key generation ...
end
```

## Expected Behavior After Implementation

1. **Decoy Filtering**: All decoy proteins removed before EFDR calculation
2. **Separate Pairing**: Target proteins only paired with other target proteins
3. **Clean EFDR**: No contamination from decoy proteins in EFDR estimates
4. **Backward Compatibility**: Works with data lacking target column (with warning)

## Validation Criteria

1. No decoy proteins in final EFDR results
2. Target and decoy proteins with same name have different `entrap_pair_id` (or decoys have missing pair_id)
3. EFDR calculations only consider target proteins
4. All existing tests pass with updated logic

## Risk Assessment

- **Low Risk**: Changes are localized to entrapment module
- **Backward Compatibility**: Maintained through conditional checks
- **Data Integrity**: No changes to underlying data structures