# Empty Protein Groups Root Cause Analysis

## Problem Found

The protein group files contain 0 rows because `apply_probit_scores_multifold!` is being called on empty/non-existent files when `skip_scoring = true`.

## The Bug Flow

1. **Step 12 (protein inference)**: Creates protein group files ONLY if `nrow(protein_groups_df) > 0`
2. **Step 14 (probit)**: With `skip_scoring = true`, calls `apply_probit_scores_multifold!` on `pg_refs`
3. **The Issue**: `pg_refs` includes references to files that don't exist (weren't created because they had 0 rows)
4. **Result**: `transform_and_write!` loads non-existent files (gets empty DataFrame) and writes them back as empty

## Why This Happens

Looking at the logging output:
- Building protein group lookup finds `n_protein_groups = 0`
- Available keys are empty: `available_keys_sample = Pioneer.ProteinKey[]`

This suggests that in Step 12, NO protein groups passed the filtering criteria.

## The Real Root Cause

The protein groups are being filtered out by `min_peptides` requirement. With `min_peptides = 2` (or higher), if proteins don't have enough peptides, they get filtered out completely.

## Solution Options

### Option 1: Skip Processing Non-Existent Files
In `apply_probit_scores_multifold!`, check if file exists before processing:
```julia
for ref in pg_refs
    if !exists(ref)
        continue
    end
    transform_and_write!(ref) do df
        # ...
    end
end
```

### Option 2: Lower min_peptides Threshold
Allow proteins with just 1 peptide during low-data scenarios

### Option 3: Handle Empty pg_refs
Don't call `apply_probit_scores_multifold!` if pg_refs is empty

## Recommended Fix

Implement Option 1 - check file existence in `apply_probit_scores_multifold!` to avoid creating empty files from non-existent sources.

## Testing Plan

1. Run test with new logging to see how many protein groups are created initially
2. Check if files are being created in Step 12
3. Verify that Step 14 doesn't create empty files from non-existent ones