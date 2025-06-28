# Multi-Fold Cross-Validation Probit Analysis Implementation Plan

## Overview
This document outlines the implementation of `perform_probit_analysis_multifold`, which automatically detects the number of CV folds from the LibraryPrecursors and performs cross-validation probit analysis accordingly.

## Key Requirements
1. Automatically detect the number of CV folds from LibraryPrecursors
2. Determine protein group CV fold from the highest-scoring constituent peptide
3. Handle cases where peptides within a protein group have different CV folds
4. Support any number of folds that exist in the precursor library

## Current System Understanding
- LibraryPrecursors contains `pid_to_cv_fold::Vector{UInt8}` with CV fold for each precursor
- CV folds are assigned during library building based on protein groups
- `getCvFold(precursors, precursor_idx)` returns the cv_fold for a precursor
- PSMs inherit cv_fold from their precursors
- Protein groups may contain peptides with different cv_folds

## Implementation Steps

### Step 1: Efficient Fold Detection from Library

```julia
# Get unique CV folds directly from LibraryPrecursors
function detect_unique_cv_folds(precursors::LibraryPrecursors)
    return sort(unique(precursors.pid_to_cv_fold))
end
```

### Step 2: Protein Group CV Fold Assignment

For each protein group:
1. Find all its constituent peptides from PSM data
2. Look up each peptide's precursor cv_fold using the library
3. Select the cv_fold of the highest-scoring peptide
4. Add cv_fold column to protein groups DataFrame

```julia
function assign_protein_group_cv_folds!(
    all_protein_groups::DataFrame,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors
)
    # Create mapping: protein_name -> (best_score, cv_fold)
    protein_to_cv_fold = Dict{String, Tuple{Float32, UInt8}}()
    
    # Process PSM files to find best peptide per protein
    for ref in pg_refs
        psm_path = get_corresponding_psm_path(ref)
        psms = DataFrame(Arrow.Table(psm_path))
        
        # Filter for PSMs that are used for protein quantification
        psms = filter(row -> row.use_for_protein_quant, psms)
        
        # Group by inferred_protein_group
        for group in groupby(psms, :inferred_protein_group)
            protein_name = first(group.inferred_protein_group)
            
            # Find highest scoring PSM
            best_idx = argmax(group.prob)
            best_score = group.prob[best_idx]
            precursor_idx = group.precursor_idx[best_idx]
            
            # Get cv_fold from library (more reliable than PSM file)
            cv_fold = getCvFold(precursors, precursor_idx)
            
            # Update if this is the best score for this protein
            if !haskey(protein_to_cv_fold, protein_name) || 
               best_score > protein_to_cv_fold[protein_name][1]
                protein_to_cv_fold[protein_name] = (best_score, cv_fold)
            end
        end
    end
    
    # Assign cv_fold to protein groups
    all_protein_groups.cv_fold = [
        get(protein_to_cv_fold, name, (0.0f0, UInt8(0)))[2]
        for name in all_protein_groups.protein_name
    ]
end
```

### Step 3: Main Function Structure

```julia
function perform_probit_analysis_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    show_improvement = true
)
    # 1. Detect unique CV folds from library
    unique_cv_folds = detect_unique_cv_folds(precursors)
    n_folds = length(unique_cv_folds)
    
    @info "Detected $n_folds CV folds in the library: $unique_cv_folds"
    
    # 2. Assign CV folds to protein groups based on best peptide
    assign_protein_group_cv_folds!(all_protein_groups, pg_refs, precursors)
    
    # 3. Check distribution
    fold_counts = countmap(all_protein_groups.cv_fold)
    @info "Protein group distribution across folds: $fold_counts"
    
    # 4. Train probit model for each fold
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides]
    models = Dict{UInt8, Vector{Float64}}()
    
    for test_fold in unique_cv_folds
        # Get training data (all folds except test_fold)
        train_mask = all_protein_groups.cv_fold .!= test_fold
        
        # Skip if insufficient training data
        n_train_targets = sum(all_protein_groups[train_mask, :target])
        n_train_decoys = sum(.!all_protein_groups[train_mask, :target])
        
        if n_train_targets < 10 || n_train_decoys < 10
            @warn "Insufficient training data for fold $test_fold (targets: $n_train_targets, decoys: $n_train_decoys)"
            continue
        end
        
        X_train = Matrix{Float64}(all_protein_groups[train_mask, feature_names])
        y_train = all_protein_groups[train_mask, :target]
        
        # Fit model
        β_fitted = fit_probit_model(X_train, y_train)
        models[test_fold] = β_fitted
        
        @info "Fitted model for fold $test_fold using $(sum(train_mask)) training samples"
    end
    
    # 5. Apply models to their respective test folds
    all_protein_groups.old_pg_score = all_protein_groups.pg_score  # Save for comparison
    
    for test_fold in unique_cv_folds
        if !haskey(models, test_fold)
            continue  # Skip if no model was trained
        end
        
        test_mask = all_protein_groups.cv_fold .== test_fold
        if sum(test_mask) > 0
            X_test = Matrix{Float64}(all_protein_groups[test_mask, feature_names])
            prob_scores = calculate_probit_scores(X_test, models[test_fold])
            all_protein_groups[test_mask, :pg_score] = Float32.(prob_scores)
        end
    end
    
    # 6. Report improvement if requested
    if show_improvement
        report_cv_improvement(all_protein_groups)
    end
    
    # 7. Update protein group files if provided
    if !isempty(pg_refs)
        apply_probit_scores_multifold!(pg_refs, models, feature_names, precursors)
    end
end
```

### Step 4: Helper Functions

1. `get_corresponding_psm_path(pg_ref)` - Map protein group file to PSM file
   ```julia
   function get_corresponding_psm_path(pg_ref::ProteinGroupFileReference)
       pg_path = file_path(pg_ref)
       # Replace "passing_proteins" with "scored_PSMs" in path
       return replace(pg_path, "passing_proteins" => "scored_PSMs")
   end
   ```

2. `apply_probit_scores_multifold!` - Update protein group files with CV-aware scoring
   ```julia
   function apply_probit_scores_multifold!(
       pg_refs::Vector{ProteinGroupFileReference},
       models::Dict{UInt8, Vector{Float64}},
       feature_names::Vector{Symbol},
       precursors::LibraryPrecursors
   )
       for ref in pg_refs
           # Load file
           df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
           
           # Assign CV folds
           assign_cv_folds_to_df!(df, ref, precursors)
           
           # Apply appropriate model to each fold
           for (fold, model) in models
               mask = df.cv_fold .== fold
               if sum(mask) > 0
                   X = Matrix{Float64}(df[mask, feature_names])
                   df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
               end
           end
           
           # Write back
           writeArrow(file_path(ref), df)
       end
   end
   ```

### Step 5: Testing Strategy

Create mock LibraryPrecursors with different cv_fold configurations:
- 2-fold: pid_to_cv_fold with values in {0, 1}
- 3-fold: pid_to_cv_fold with values in {0, 1, 2}
- 5-fold: pid_to_cv_fold with values in {0, 1, 2, 3, 4}

Test cases:
1. Automatic detection of fold count from library
2. Correct CV fold assignment to protein groups based on best peptide
3. Model training with proper hold-out
4. Handling of imbalanced folds
5. Edge cases (empty folds, single protein)

## Advantages of This Approach

1. **Efficiency**: Uses LibraryPrecursors directly instead of scanning PSM files
2. **Reliability**: Gets cv_fold from the source of truth (library)
3. **Flexibility**: Automatically adapts to any number of folds
4. **Consistency**: Maintains cv_fold assignments from library building
5. **Robustness**: Handles peptides with different cv_folds within protein groups