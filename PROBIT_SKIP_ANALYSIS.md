# Probit Analysis Skip: Detailed Workflow Comparison and Error Analysis

## Overview

This document details the differences in the SearchDIA pipeline when probit analysis is performed versus when it is skipped due to insufficient data (e.g., < 50 targets or < 50 decoys).

## Workflow Comparison

### When Probit Analysis RUNS (Normal Flow)

1. **Step 12: Protein Inference** (`perform_protein_inference_pipeline`)
   - Each PSM file is processed independently
   - Protein groups are created **per file**
   - Results: 3 protein group files in `temp_data/passing_proteins/`
   - Each file contains only the protein groups for its corresponding PSMs

2. **Step 14: Probit Regression** (`perform_protein_probit_regression`)
   - **CRITICAL MERGE**: All protein group files are loaded into a single DataFrame
     ```julia
     all_protein_groups = DataFrame()
     for pg_path in passing_pg_paths
         append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
     end
     ```
   - Probit model is trained on the merged data
   - Scores are recalculated using probit features
   - **Files are written back** via `apply_probit_scores_multifold!`
   - Each file still contains its original protein groups, but with updated scores

3. **Step 23: Update PSMs** (`update_psms_with_probit_scores_refs`)
   - For each PSM file, loads its **paired** protein group file
   - Builds lookup dictionary from that single file
   - Updates PSMs with protein scores
   - PSMs get columns: `pg_score`, `global_pg_score`, `pg_qval`, `qlobal_pg_qval`

4. **MaxLFQ Processing**
   - Merges all PSM files
   - Filters based on `pg_qval` and `qlobal_pg_qval` thresholds
   - Runs LFQ algorithm
   - Creates `protein_groups_long.arrow`
   - Calls `writeProteinGroupsCSV` to create `protein_groups_wide.arrow`

### When Probit Analysis is SKIPPED (Problem Flow - Before Fix)

1. **Step 12: Protein Inference** (Same as above)
   - Protein groups created per file
   - Fragmented across 3 files

2. **Step 14: Skipped**
   - Original implementation: Just wrote files back unchanged
   - **MISSING**: No merge step
   - **PROBLEM**: Protein groups remain fragmented

3. **Step 23: Update PSMs** 
   - Loads paired protein group file
   - **CRITICAL ISSUE**: Protein group lookup is incomplete
   - Many PSMs reference protein groups in OTHER files
   - Result: Most PSMs get `missing` values for protein scores

4. **MaxLFQ Processing**
   - PSMs with `missing` scores are filtered out
   - Result: 0 rows in `protein_groups_long.arrow`
   - `writeProteinGroupsCSV` fails or produces empty output

### After Our Fix (Current State)

We modified Step 14 (when probit is skipped) to:
1. Merge all protein groups (mimicking probit behavior)
2. Sort by original scores
3. **Write ALL protein groups to EACH file**
   - This ensures complete lookup in Step 23
   - Should prevent missing scores

## Current Error Analysis

Despite our fix, we're getting:
```
SystemError: opening file ".../protein_groups_wide.arrow": No such file or directory
```

### Theory: Why protein_groups_wide.arrow is Still Missing

**Hypothesis 1: Duplicate Protein Groups**
- By writing ALL protein groups to EACH file, we've created duplicates
- When Step 23 processes multiple files, it might be creating duplicate entries
- This could cause issues in the LFQ algorithm or writeProteinGroupsCSV

**Hypothesis 2: Memory or Processing Issue**
- Writing all protein groups to each file increases file size 3x
- This might cause memory issues or processing failures in LFQ
- The function might be failing before creating the output file

**Hypothesis 3: Mismatch in Protein Group Structure**
- The merged protein groups might have a different structure than expected
- Missing columns or incorrect data types could cause silent failures
- `writeProteinGroupsCSV` might be failing to process the data

**Hypothesis 4: File Path or Permissions Issue**
- The output directory might not be writable
- The file path construction might be incorrect
- Previous partial writes might be blocking new file creation

## Recommended Investigation Steps

1. **Check if protein_groups_long.arrow is created**
   - If yes: Problem is in `writeProteinGroupsCSV`
   - If no: Problem is in LFQ processing

2. **Add logging to track duplicate protein groups**
   - Count unique vs total protein groups after merging
   - Check if duplicates are causing issues

3. **Verify the merge creates valid data**
   - Check column names and types
   - Ensure no missing required columns

4. **Consider Alternative Fix**
   - Instead of writing ALL groups to EACH file
   - Build a GLOBAL lookup in Step 23 from ALL files
   - This avoids duplicates while ensuring complete coverage

## Alternative Solution Proposal

Instead of duplicating all protein groups across files, modify Step 23 to:

```julia
# Build global lookup from ALL protein group files
global_pg_lookup = Dict{ProteinKey, Tuple{Float32, Float32}}()
for ref in all_pg_refs  # ALL refs, not just paired
    pg_table = Arrow.Table(file_path(ref))
    # Add to global lookup...
end

# Then process each PSM file using the global lookup
```

This would:
- Avoid data duplication
- Ensure complete protein group coverage
- Maintain original file structure
- Reduce memory usage

## Conclusion

The current fix successfully addresses the lookup issue but may have introduced a new problem with duplicate protein groups. The missing `protein_groups_wide.arrow` suggests that either:
1. The LFQ processing is failing due to duplicates
2. The writeProteinGroupsCSV function cannot handle the modified data structure
3. There's a file system or path issue

The most robust solution would be to implement the alternative approach: building a global lookup in Step 23 without modifying the protein group files.