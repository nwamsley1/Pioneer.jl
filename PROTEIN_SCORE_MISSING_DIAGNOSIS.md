# Protein Score Missing Issue Diagnosis

## Current State Analysis

### What We Know:
1. **Protein groups ARE created** (2264 groups across 3 files in Step 12)
2. **Skip_scoring implementation works** (processes files, keeps original scores)
3. **Step 23 fails to find protein groups** (all PSMs get missing scores)
4. **MaxLFQ receives 0 protein groups** after filtering

### The Data Flow:
```
Step 12: Create protein groups (2264 groups) ✓
    ↓
Step 14: Process with skip_scoring=true ✓
    ↓
Step 23: Update PSMs with protein scores ✗ (all missing)
    ↓
MaxLFQ: Filter by pg_qval thresholds ✗ (0 groups remain)
```

## Key Issue: Step 23 Lookup Failure

The protein group lookup in Step 23 builds a key like this:
```julia
key = ProteinKey(
    psms_df[i, :inferred_protein_group],  # From PSM
    psms_df[i, :target],                   # From PSM
    psms_df[i, :entrapment_group_id]       # From PSM
)
```

And tries to find it in:
```julia
pg_score_lookup[key]  # Built from paired protein group file only
```

## Hypothesis: Column Name Mismatch

The protein group files use `:entrap_id` but PSMs use `:entrapment_group_id`. This could cause the key mismatch.

## Added Logging Will Show:

1. **In LFQ function:**
   - Total input rows and unique proteins
   - Number of missing pg_qval and qlobal_pg_qval values
   - Number of rows filtered at each pipeline step

2. **In update_psms_with_probit_scores_refs:**
   - Number of protein groups loaded from each file
   - First 5 missing keys and available keys when lookup fails

3. **After pipeline operations:**
   - Warning when batch is filtered to 0 rows

## Expected Output from Logging:

When you run the test again, we should see:
1. How many PSMs have missing pg_qval values
2. The exact keys being looked up vs what's available
3. Where the filtering happens in MaxLFQ

## Next Steps:

1. **Run the test** with the new logging
2. **Analyze the output** to see:
   - Are the column names matching?
   - Are the protein names matching?
   - Are entrapment IDs matching?
3. **Fix the root cause** based on findings

## Potential Fixes:

### Fix 1: Column Name Alignment
Ensure consistent column naming between protein groups (`:entrap_id`) and PSMs (`:entrapment_group_id`)

### Fix 2: Global Lookup
Build lookup from ALL protein group files, not just paired file:
```julia
# Build global lookup
for all_pg_refs...
    # Add to lookup
end
```

### Fix 3: Skip Missing Check
For rows with `use_for_protein_quant = false`, don't require protein scores

### Fix 4: Handle Missing Values
In MaxLFQ, handle missing pg_qval values gracefully instead of filtering them out