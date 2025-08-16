# Plan: Exact Probit Workflow Mimicry Without Scoring

## Problem Statement

`protein_groups_long.arrow` still has 0 rows even after our fix. We need to exactly mimic what probit does, minus the actual score transformation.

## Current Probit Workflow Analysis

### What Probit Does (from perform_probit_analysis_multifold)

```julia
# Lines 1662-1669: Load and merge all protein groups
@info "Using in-memory probit regression"
all_protein_groups = DataFrame()
for pg_path in passing_pg_paths
    if isfile(pg_path) && endswith(pg_path, ".arrow")
        append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
    end
end

# Lines 1671-1675: Check targets/decoys
n_targets = sum(all_protein_groups.target)
n_decoys = sum(.!all_protein_groups.target)

# Lines 1756-1860: Core probit processing
# 1. Build CV fold mapping
# 2. Assign CV folds to protein groups
# 3. Train models per fold
# 4. Apply scores per fold
# 5. Calculate improvement metrics

# Lines 1853-1856: CRITICAL - Update individual files
if !isempty(pg_refs)
    @info "Updating individual protein group files with probit scores"
    apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names)
end
```

### The Key Function: apply_probit_scores_multifold!

```julia
# Lines 1675-1715
function apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names)
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Assign CV folds
            cv_folds = Vector{UInt8}(undef, nrow(df))
            for (i, protein_name) in enumerate(df.protein_name)
                if haskey(protein_to_cv_fold, protein_name)
                    cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
                else
                    cv_folds[i] = UInt8(0)
                end
            end
            df[!, :cv_fold] = cv_folds
            
            # Save original scores
            df[!, :old_pg_score] = copy(df.pg_score)
            
            # Apply model to each fold (THIS IS THE ONLY SCORING PART)
            for (fold, model) in models
                mask = df.cv_fold .== fold
                if sum(mask) > 0
                    X = Matrix{Float64}(df[mask, feature_names])
                    df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
                end
            end
            
            # Sort by pg_score and target
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Remove temporary cv_fold column
            select!(df, Not(:cv_fold))
            
            return df
        end
    end
end
```

## Current "Skip Probit" Implementation Issues

Our current implementation writes ALL protein groups to EACH file:
```julia
# This is WRONG - creates duplicates
for ref in pg_refs
    writeArrow(file_path(ref), all_protein_groups)  # Writing ALL to EACH
end
```

The problem: This creates duplicate protein groups across files, which breaks downstream processing.

## The Correct Solution

We need to mimic `apply_probit_scores_multifold!` WITHOUT the scoring:

```julia
function apply_original_scores_like_probit!(pg_refs::Vector{ProteinGroupFileReference})
    # Process each file individually (like probit does)
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Skip CV fold assignment (not needed without scoring)
            
            # Skip saving old scores (not needed)
            
            # Skip applying models (THE ONLY PART WE ACTUALLY SKIP)
            
            # CRITICAL: Still sort by pg_score and target (original scores)
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            return df
        end
    end
end
```

Wait, this is too simple. Let me analyze what's REALLY different...

## The Real Issue: Data Flow

### When Probit Runs:
1. Loads all groups into memory (merge)
2. Processes them together
3. **Writes back to ORIGINAL files** (each file keeps its own groups)

### When Probit is Skipped (Our Current Fix):
1. Loads all groups into memory (merge) ✓
2. Skips processing ✓
3. **Writes ALL groups to EACH file** ✗ (This is wrong!)

## The Correct Fix

We need to:
1. Load all protein groups (for consistency checking)
2. Process each file individually
3. Write each file back with ONLY its original groups

```julia
function mimic_probit_without_scoring!(pg_refs, passing_pg_paths)
    # Step 1: Load all (like probit does) - for validation only
    all_protein_groups = DataFrame()
    for pg_path in passing_pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
        end
    end
    
    # Log statistics (like probit does)
    n_targets = sum(all_protein_groups.target)
    n_decoys = sum(.!all_protein_groups.target)
    @info "Processing protein groups without probit" total=nrow(all_protein_groups) targets=n_targets decoys=n_decoys
    
    # Step 2: Process each file INDIVIDUALLY (like apply_probit_scores_multifold!)
    for ref in pg_refs
        if exists(ref)
            # Load this file's data
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            
            # Sort by pg_score and target (mimicking probit behavior)
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Write back THIS FILE'S data only
            writeArrow(file_path(ref), df)
            
            @info "Processed protein group file" path=file_path(ref) n_groups=nrow(df)
        end
    end
end
```

## Implementation Plan

1. **Remove the current fix** that writes all groups to each file
2. **Implement the corrected function** that:
   - Loads all groups (for consistency with probit path)
   - Processes each file individually
   - Sorts each file by score
   - Writes back only that file's groups
3. **Test** that protein groups are properly found in Step 23

## Why This Should Work

- Each file maintains its original protein groups (no duplicates)
- Files are sorted consistently (like probit does)
- The merge step ensures we see all data (for logging/validation)
- Step 23's lookup will find the groups in the paired files
- No structural changes that could break downstream processing

## Key Insight

The probit workflow doesn't redistribute protein groups across files - it only updates scores within each file. Our fix was creating duplicates by writing all groups to each file, which breaks the 1:1 pairing assumption in Step 23.