# Fix Plan: Protein Group File Writing Issue in ScoringSearch

## Problem Statement

When ScoringSearch skips probit analysis due to insufficient data (e.g., targets: 2224, decoys: 40), the final `protein_groups_long.arrow` and `protein_groups_wide.arrow` files are not created in the main output directory, causing MaxLFQSearch to fail.

## Current Data Flow

### Normal Flow (Sufficient Data for Probit)
1. **ScoringSearch Step 12**: Creates protein groups in `temp_data/passing_proteins/`
2. **ScoringSearch Step 14**: Performs probit regression (updates scores)
3. **ScoringSearch Step 15-23**: Calculates global scores, q-values, updates PSMs
4. **IntegrateChromatogramsSearch**: Integrates peaks, updates PSMs
5. **MaxLFQSearch**: 
   - Calls `LFQ()` which reads PSM files and creates `protein_groups_long.arrow`
   - Calls `writeProteinGroupsCSV()` which creates `protein_groups_wide.arrow`

### Broken Flow (Insufficient Data - Probit Skipped)
1. **ScoringSearch Step 12**: Creates protein groups in `temp_data/passing_proteins/`
2. **ScoringSearch Step 14**: SKIPS probit regression (prints warning)
3. **ScoringSearch Step 15-23**: Still executes but protein groups remain in temp folder
4. **IntegrateChromatogramsSearch**: Integrates peaks, updates PSMs
5. **MaxLFQSearch**: 
   - Calls `LFQ()` which expects to read PSM files with protein info
   - **FAILS** because `protein_groups_long.arrow` doesn't exist

## Root Cause Analysis (UPDATED)

After detailed code analysis, the exact issue has been identified:

### What Happens During Normal Flow (Probit Runs)
1. **Step 14** (`perform_protein_probit_regression`): 
   - Loads protein groups from `temp_data/passing_proteins/`
   - Performs probit regression to refine scores
   - **CRITICAL**: Overwrites `pg_score` column with probit scores (lines ~1765-1783 in utils.jl)
   - **CRITICAL**: Writes files back with `writeArrow(file_path(ref), df)`

2. **Step 23** (`update_psms_with_probit_scores_refs`):
   - Reads the updated protein group files
   - Adds columns to PSMs: `pg_score`, `global_pg_score`, `pg_qval`, `global_pg_qval`, `pg_pep`

3. **MaxLFQ**: 
   - Reads PSMs expecting: `pg_qval`, `global_pg_qval`, `use_for_protein_quant`, `inferred_protein_group`

### What Happens When Probit is Skipped
1. **Step 14**: 
   - Prints "Skipping Probit analysis: insufficient data"
   - **PROBLEM**: Does NOT write protein group files back to disk
   - **PROBLEM**: Protein groups retain original log-sum scores, not refined

2. **Step 23**: 
   - Still runs, tries to read protein group files
   - May fail or use stale/missing data

3. **MaxLFQ**: 
   - Fails because PSMs lack required columns or have incorrect values

## Investigation Findings

### Where Files Are Written
- **protein_groups_long.arrow**: Created by `LFQ()` function in `src/utils/maxLFQ.jl`
- **protein_groups_wide.arrow**: Created by `writeProteinGroupsCSV()` in MaxLFQSearch
- Both depend on PSM files having proper protein group assignments

### Critical Discovery
The `LFQ()` function (line 377 in maxLFQ.jl) reads from PSM files, not from protein group files:
```julia
prot = DataFrame(Arrow.Table(file_path(prot_ref)))  # prot_ref is a PSM file reference
```

This means the protein group information must be properly embedded in the PSM files by ScoringSearch.

## Exact Code Locations Requiring Fix

### Critical Code in `perform_protein_probit_regression` (utils.jl)

**Lines ~1740-1750**: Check for sufficient data
```julia
if n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000
    # Performs probit analysis
    perform_probit_analysis_multifold(...)
else
    @info "Skipping Probit analysis: insufficient data"
    # MISSING: No file write-back happens here!
end
```

**Lines ~1765-1783** (Inside `perform_probit_analysis_multifold`): 
```julia
# This ONLY runs when probit is performed:
df[!, :pg_score] = Float32.(prob_scores)
sort!(df, [:pg_score, :target], rev = [true, true])
writeArrow(file_path(ref), df)  # Critical write-back
```

## Proposed Solution

### Option 1: Minimal Fix - Ensure Files Are Written (SELECTED FOR IMPLEMENTATION)

Modify `perform_protein_probit_regression` to always write protein group files back, even when probit is skipped:

### Option 2: Create Protein Files Directly

Add a new step in ScoringSearch to write protein group files directly when probit is skipped:

1. After Step 23, check if probit was skipped
2. If skipped, merge protein group files from temp folder
3. Write to main output directory as `protein_groups_long.arrow`

### Option 3: Modify MaxLFQ to Handle Missing Files

Update MaxLFQSearch to handle the case where protein files don't exist:

1. Check if `protein_groups_long.arrow` exists
2. If not, create it from PSM files using protein inference
3. Continue with normal processing

## Recommended Implementation (Option 1)

### Implementation Details

In `perform_protein_probit_regression` (utils.jl), modify the else block:

```julia
if n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000
    # Existing probit analysis code
    perform_probit_analysis_multifold(...)
else
    @info "Skipping Probit analysis: insufficient data (targets: $n_targets, decoys: $n_decoys)"
    
    # NEW CODE: Write protein group files back even without probit
    for ref in pg_refs
        if exists(ref)
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            # Keep original pg_score values (log-sum scores)
            # But ensure file is properly sorted for downstream processing
            sort!(df, [:pg_score, :target], rev = [true, true])
            writeArrow(file_path(ref), df)
        end
    end
    @info "Wrote protein group files with original scores"
end
```

This ensures:
1. Protein group files are always written back to disk
2. Files are properly sorted for downstream processing
3. Step 23 can read the files and update PSMs correctly
4. MaxLFQ receives PSMs with all required columns

### Additional Fix Required

In `src/utils/maxLFQ.jl` line 383, fix the typo:
```julia
# Change from:
(:qlobal_pg_qval, q_value_threshold)
# To:
(:global_pg_qval, q_value_threshold)
```

## Testing Plan

1. **Create Test Case**: 
   - Run with data that has < 50 decoys to trigger skip
   - Verify warning message appears
   - Check that protein files are created

2. **Verify No Regression**:
   - Run with normal data (sufficient targets/decoys)
   - Ensure protein files are created only once
   - Verify scores are correct

3. **Check File Contents**:
   - Compare protein files from skipped vs normal runs
   - Ensure essential columns are present
   - Verify PSM files have `inferred_protein_group` column

## Implementation Priority

**High Priority**:
- Ensure PSMs always have `inferred_protein_group` column
- Make sure basic protein scores are available

**Medium Priority**:
- Add proper fallback scoring when probit is skipped
- Improve warning messages to indicate what's being done

**Low Priority**:
- Optimize to avoid duplicate file writes
- Add configuration option to force probit even with low numbers

## Files to Modify

1. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`
   - Add probit skip tracking
   - Ensure PSM updates always happen

2. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`
   - Modify `perform_protein_probit_regression` to handle skip case
   - Add `update_psms_with_basic_protein_scores_refs` function

3. `src/utils/maxLFQ.jl` (optional)
   - Add better error messages when expected columns are missing

## Success Criteria

1. Pipeline completes without errors when probit is skipped
2. `protein_groups_long.arrow` and `protein_groups_wide.arrow` are created
3. No duplicate file writes in normal case
4. PSM files contain necessary protein information
5. MaxLFQSearch runs successfully