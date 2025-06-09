# ML Protein Scoring FDR Integration

## Overview

This document describes how ML-based protein scoring (probit regression) is integrated with FDR filtering in the ScoringSearch module.

## Problem Statement

Previously, the ML protein scoring calculated `probit_score` and `probit_qval` columns but these were not used for final FDR filtering. The system continued to use the traditional `pg_score` for filtering, effectively ignoring the ML improvements.

## Solution Implemented

Instead of maintaining separate columns for ML scores, the system now **overwrites** the existing scoring columns when ML scoring is performed:

1. `pg_score` → Overwritten with ML-based probit scores
2. `pg_qval` → Overwritten with ML-based q-values
3. `global_pg_score` → Recalculated using new ML scores
4. `global_pg_qval` → Recalculated using new global scores

## Implementation Details

### Out-of-Memory (OOM) Probit Analysis

In `perform_probit_analysis_oom`:
```julia
# Original pg_scores saved for comparison reporting
original_pg_scores = copy(df.pg_score)

# Overwrite with ML scores
df[!, :pg_score] = Float32.(prob_scores)
df[!, :pg_qval] = Float32.(probit_qvalues)

# After all files processed, update global scores
```

### In-Memory Probit Analysis

In `perform_probit_analysis`:
```julia
# Calculate ML scores
prob_scores = calculate_probit_scores(X, β_fitted, X_mean, X_std)
probit_qvalues = calculate_qvalues_from_scores(prob_scores, y)

# Overwrite existing columns
all_protein_groups[!, :pg_score] = Float32.(prob_scores)
all_protein_groups[!, :pg_qval] = Float32.(probit_qvalues)

# Recalculate global scores
# Write back to individual files
```

### Global Score Update

After ML scoring, a final pass updates global scores:
1. Collect maximum scores across all files for each protein
2. Update `global_pg_score` with new maximums
3. Recalculate `global_pg_qval` based on updated global scores

## Benefits

1. **No downstream changes required**: The existing FDR filtering code continues to work unchanged
2. **Seamless integration**: ML scoring automatically improves results when performed
3. **Transparent**: Users see improved results without configuration changes
4. **Backwards compatible**: If ML scoring fails or is skipped, original scores remain

## Usage

The system automatically uses ML scores when available:

1. If probit analysis succeeds → `pg_score` contains ML scores
2. If probit analysis is skipped → `pg_score` contains traditional scores
3. FDR filtering always uses `pg_score` and `pg_qval` columns

## Verification

To verify ML scoring is being used:
1. Check log messages for "ML-based protein scoring complete"
2. Compare reported improvements in probit analysis output
3. Final filtered results should match probit-reported numbers at each FDR threshold

## Technical Notes

- Original scores are preserved temporarily for comparison reporting
- All score columns maintain Float32 precision
- Global scores are recalculated to maintain consistency
- File writing uses `writeArrow` for Windows compatibility