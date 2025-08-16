# Plan: Add Skip-Scoring Option to Probit Functions

## Overview

Instead of having a separate code path when probit is skipped, we'll modify the existing probit functions to accept a `skip_scoring` parameter. This ensures we follow the EXACT same workflow, just without the actual score transformation.

## Implementation Strategy

### 1. Modify perform_probit_analysis_multifold

**Current signature:**
```julia
function perform_probit_analysis_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{...}} = nothing,
    show_improvement = true
)
```

**New signature:**
```julia
function perform_probit_analysis_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{...}} = nothing,
    show_improvement = true,
    skip_scoring = false  # NEW PARAMETER
)
```

**Changes needed:**
1. Skip model training when `skip_scoring = true`
2. Skip score calculation
3. Still perform all other operations (CV fold assignment, sorting, file writing)

### 2. Modify apply_probit_scores_multifold!

**Current function:**
```julia
function apply_probit_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{...},
    models::Dict{UInt8, Vector{Float64}},
    feature_names::Vector{Symbol}
)
```

**New signature:**
```julia
function apply_probit_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{...},
    models::Dict{UInt8, Vector{Float64}},
    feature_names::Vector{Symbol};
    skip_scoring = false  # NEW PARAMETER
)
```

**Changes in the function:**
```julia
for ref in pg_refs
    transform_and_write!(ref) do df
        # Assign CV folds (STILL DO THIS)
        cv_folds = Vector{UInt8}(undef, nrow(df))
        for (i, protein_name) in enumerate(df.protein_name)
            if haskey(protein_to_cv_fold, protein_name)
                cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
            else
                cv_folds[i] = UInt8(0)
            end
        end
        df[!, :cv_fold] = cv_folds
        
        # Save original scores (STILL DO THIS)
        df[!, :old_pg_score] = copy(df.pg_score)
        
        # Apply models (SKIP THIS IF skip_scoring = true)
        if !skip_scoring
            for (fold, model) in models
                mask = df.cv_fold .== fold
                if sum(mask) > 0
                    X = Matrix{Float64}(df[mask, feature_names])
                    df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
                end
            end
        end
        # If skip_scoring, pg_score remains unchanged
        
        # Sort (STILL DO THIS)
        sort!(df, [:pg_score, :target], rev = [true, true])
        
        # Remove temp column (STILL DO THIS)
        select!(df, Not(:cv_fold))
        
        return df
    end
end
```

### 3. Modify perform_probit_analysis_oom (Out-of-Memory version)

**Add skip_scoring parameter:**
```julia
function perform_probit_analysis_oom(
    pg_refs::Vector{ProteinGroupFileReference},
    total_protein_groups::Int,
    max_protein_groups_in_memory::Int,
    qc_folder::String;
    skip_scoring = false  # NEW PARAMETER
)
```

**Changes:**
1. Skip model training when `skip_scoring = true`
2. In the file processing loop, skip score calculation but still sort and write

### 4. Update perform_protein_probit_regression

**Current check:**
```julia
if n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000
    perform_probit_analysis_multifold(...)
else
    @info "Skipping Probit analysis: insufficient data"
    # Current workaround code
end
```

**New approach:**
```julia
if n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000
    # Run normal probit
    perform_probit_analysis_multifold(
        all_protein_groups,
        qc_folder,
        pg_refs,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        skip_scoring = false  # Normal operation
    )
else
    @info "Insufficient data for probit scoring - using original scores"
    # Run SAME function but skip scoring
    perform_probit_analysis_multifold(
        all_protein_groups,
        qc_folder,
        pg_refs,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        skip_scoring = true  # Skip scoring only
    )
end
```

## Benefits of This Approach

1. **Exact same code path**: We follow the identical workflow regardless of data sufficiency
2. **All side effects preserved**: Any file operations, state changes, or data transformations happen consistently
3. **Easier maintenance**: Only one code path to maintain
4. **Debugging clarity**: Can easily toggle skip_scoring to test
5. **Minimal changes**: Just adding conditional checks around actual scoring

## Implementation Steps

1. Add `skip_scoring` parameter to all relevant functions
2. Add conditional checks around:
   - Model training
   - Score calculation
   - Improvement metrics (when skip_scoring = true)
3. Update the call sites to pass `skip_scoring = true` when appropriate
4. Test that both paths produce valid output

## Code Changes Summary

### Files to modify:
1. `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl`
   - `perform_probit_analysis_multifold` (lines ~1719-1860)
   - `apply_probit_scores_multifold!` (lines ~1675-1715)
   - `perform_probit_analysis_oom` (lines ~850-950)
   - `perform_protein_probit_regression` (lines ~630-720)

### Key principle:
- When `skip_scoring = true`:
  - DO: All file operations, sorting, CV fold assignment, data structure updates
  - SKIP: Model training, score calculation, improvement metrics

This approach guarantees we follow the exact same workflow, ensuring any hidden side effects or state changes are preserved.