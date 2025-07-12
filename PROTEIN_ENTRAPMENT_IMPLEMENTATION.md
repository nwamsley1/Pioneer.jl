# Protein-Level Entrapment Analysis Implementation Checklist

## Status: COMPLETED ✅

All phases of the protein-level entrapment analysis have been successfully implemented:
- ✅ Core protein pairing functions (`protein_entrapment_pairing.jl`)
- ✅ Protein scoring functions (`protein_scoring.jl`)
- ✅ API integration (`run_protein_efdr_analysis` in `api.jl`)
- ✅ Comprehensive tests (`test_protein_entrapment.jl`)
- ✅ Usage examples (`protein_analysis.jl`)
- ✅ **BUG FIX**: Created protein-specific EFDR functions (`protein_efdr.jl`) to avoid precursor_idx requirement

## Overview
Implementation checklist for protein-level entrapment analysis in EntrapmentAnalysis module. This checklist is written for Claude to track progress and implementation details.

## Reference Data
**Example protein data file**: `/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow`

### First Task: Investigate Protein Data Structure
- [x] Load and examine the protein data file using Arrow.jl
- [x] Document the exact column names and types
- [x] Identify patterns in protein naming
- [x] Check how entrap_id values are distributed
- [x] Determine if ms_file_idx exists or if file_name is used
- [x] Note any unexpected columns or data patterns

#### Findings from Data Investigation:
Based on the Arrow schema provided:
- **File identification**: Uses `file_name` column (String), NOT `ms_file_idx`
- **Entrapment**: `entrap_id` column (UInt8) where 0 = original, >0 = entrapment versions
- **Protein names**: Full protein names in `protein` column (e.g., "sp|P00330|ADH1_YEAST")
- **Score columns**: 
  - `pg_score` (Float32) - per-file protein group score
  - `global_pg_score` (Float32) - global protein group score
  - `pg_pep` (Float32) - posterior error probability
- **Q-value columns**:
  - `qval` (Float32) - per-file q-value
  - `global_qval` (Float32) - global q-value
- **Other important columns**:
  - `target` (Bool) - target vs decoy
  - `species` (String) - organism
  - `peptides` (Array of UInt32) - peptide indices
  - `n_peptides` (UInt32) - peptide count
  - `abundance` (Float32) - protein abundance

## Phase 1: Core Protein Pairing Functions
**File to create**: `src/core/protein_entrapment_pairing.jl`

### 1. `assign_protein_entrapment_pairs!(df::DataFrame)`

#### Implementation steps:
- [x] Start with examining the actual protein data to understand structure
- [x] Create function with proper error checking for required columns
- [x] Implementation logic:
  ```julia
  # Key insight: Unlike precursors, proteins pair by name alone
  # Group by: protein (name only)
  # Within each protein name group:
  #   - Find all rows where entrap_id == 0 (originals)
  #   - Find all rows where entrap_id > 0 (entrapments)
  #   - Use same round-robin pairing as precursor version
  ```
- [x] Handle edge cases:
  - [x] Proteins with no entrapment versions
  - [x] Proteins with only entrapment versions (no original)
  - [x] Multiple files having same protein
- [x] Test with actual data from the reference file

### 2. `add_protein_entrap_pair_ids!(protein_results::DataFrame, protein_library::DataFrame)`

#### Implementation considerations:
- [x] First check if this function is even needed - the protein results might already have entrap_id
- [x] If needed, map by protein name (not by index like precursors)
- [x] Consider if protein_library is necessary or if all info is in results

## Phase 2: Protein Scoring Functions
**File to create**: `src/core/protein_scoring.jl`

### 3. `add_original_target_protein_scores!(protein_results::DataFrame; score_col=:pg_score)`

#### Key differences from precursor version:
- [x] No need for library_precursors parameter (entrap_id is in results)
- [x] Simpler lookup - just use entrap_id directly
- [x] Check actual data for:
  - [x] Does ms_file_idx exist? Or use file_name?
  - [x] What score columns exist? (pg_score, global_pg_score, pg_pep)
- [x] Implementation approach:
  ```julia
  # Build dictionary: (ms_file_idx/file_name, pair_id) -> target_score
  # For each row:
  #   if entrap_id == 0: use own score
  #   else: lookup target score from same pair
  ```

### 4. `create_global_protein_results_df(protein_results::DataFrame; score_col=:global_pg_score)`

#### Investigation needed:
- [x] Check if protein results already have global scores calculated
- [x] Determine grouping key: just protein name? Or protein + file?
- [x] Verify score column names in actual data

## Phase 3: API Integration
**Decision**: Extend `src/api.jl` or create new `src/protein_api.jl`

### 5. `run_protein_efdr_analysis(protein_results_path::String; kwargs...)`

#### Pre-implementation checks:
- [x] Load example data and check structure
- [x] Determine if protein data has all needed info (no separate library needed?)
- [x] Check default score/qval pairs from actual data
- [x] Verify if target column exists

#### Implementation strategy:
- [x] Auto-detect protein vs precursor data by columns
- [x] Simplified flow if entrap_id already in results
- [x] Reuse existing EFDR calculation functions
- [x] Handle file-specific vs global scores appropriately

## Phase 4: Testing Strategy

### Create `test/test_protein_entrapment.jl`:
- [x] First create test data that mimics real structure
- [x] Use actual column names from reference file
- [x] Test scenarios:
  ```julia
  # Scenario 1: Single protein with one entrapment
  # Scenario 2: Single protein with multiple entrapments  
  # Scenario 3: Multiple files, same proteins
  # Scenario 4: Missing scores handling
  ```

## Phase 5: Example Creation

### Create `examples/protein_analysis.jl`:
- [x] Use actual file path for real example
- [x] Show both automated and manual analysis
- [x] Demonstrate difference between per-file and global analysis

## Implementation Order

1. [x] **First**: Load and thoroughly examine `/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow`
2. [x] **Second**: Create protein_entrapment_pairing.jl based on actual data structure
3. [x] **Third**: Create protein_scoring.jl adapted to actual score columns
4. [x] **Fourth**: Extend API with auto-detection
5. [x] **Fifth**: Create tests using realistic data
6. [x] **Sixth**: Create examples with real file

## Key Questions to Answer with Data Investigation

1. [x] Does the protein data have `ms_file_idx` or just `file_name`? - **Answer: file_name only**
2. [x] Are all entrap_id values present in results, or do we need a library? - **Answer: entrap_id is in results**
3. [x] What are the exact score column names? (pg_score vs protein_score, etc.) - **Answer: pg_score, global_pg_score, pg_pep**
4. [x] How are protein names formatted? Any special parsing needed? - **Answer: Full UniProt format like "sp|P00330|ADH1_YEAST"**
5. [x] Is there already a global score, or do we calculate it? - **Answer: global_pg_score already exists**
6. [x] What's the typical size of protein data? Performance considerations? - **Answer: Much smaller than precursor data**

## Commit Messages

1. "Investigate protein data structure and add implementation plan"
2. "Add protein entrapment pairing functions"
3. "Add protein scoring functions for EFDR analysis"
4. "Extend API to support protein-level EFDR analysis"
5. "Add comprehensive tests for protein entrapment"
6. "Add protein analysis examples and complete implementation"

## Notes for Implementation

- Keep similar API to precursor analysis for consistency
- Protein analysis should be simpler (no mods, charge states)
- Consider performance with large datasets
- Maintain backwards compatibility
- Document differences from precursor analysis clearly

## Bug Fix: Protein EFDR Analysis Error

### Problem
The initial implementation reused the precursor-level `add_efdr_columns!` function which required a `precursor_idx` column. Protein data doesn't have precursor indices, causing the error:
```
ERROR: DataFrame must have :precursor_idx column
```

### Solution
Created protein-specific EFDR functions in `src/core/protein_efdr.jl`:
- `add_protein_efdr_columns!` - Works directly with protein data using `entrap_id`
- `compare_protein_efdr_methods` - Protein-specific comparison without precursor_idx
- `calculate_protein_efdr_calibration_error` - Protein-specific calibration

These functions:
1. Extract `entrap_id` directly from the protein DataFrame
2. Don't require precursor_idx or library DataFrames
3. Reuse the existing EFDR calculation methods (CombinedEFDR, PairedEFDR)
4. Follow the same interface patterns for consistency

The API was updated to use these protein-specific functions, eliminating the need for dummy DataFrames and temporary columns.

### Parameter Ordering Bug Fix
**Problem**: InexactError when trying to convert float q-values to integers:
```
ERROR: InexactError: Int64(9.5229028374888e-5)
```

**Cause**: Parameters were passed to EFDR constructors in wrong order - `qvals` and `entrap_labels` were swapped.

**Solution**: Fixed parameter order in `protein_efdr.jl` line 86:
```julia
# Wrong: method_type(scores, original_target_scores, qvals, entrap_labels; r=r)
# Correct:
method = method_type(scores, original_target_scores, entrap_labels, qvals; r=r)
```

## Next Immediate Action
The protein-level entrapment analysis is now fully functional and ready for use.