# Detailed Step-by-Step Comparison: Probit Run vs Skip

## Context
When analyzing protein groups in ScoringSearch Step 14, the system decides whether to run probit analysis based on data sufficiency (>50 targets, >50 decoys, >1000 total groups).

## Path 1: When Probit Analysis RUNS (Sufficient Data)

### Step 14A: Initial Setup (lines 643-653)
```julia
passing_pg_paths = [file_path(ref) for ref in pg_refs]
total_protein_groups = 0
for ref in pg_refs
    if exists(ref)
        total_protein_groups += row_count(ref)
    end
end
max_protein_groups_in_memory_limit = 5 * max_psms_in_memory
```

### Step 14B: In-Memory Path Selected (lines 662-669)
```julia
@info "Using in-memory probit regression"
all_protein_groups = DataFrame()
for pg_path in passing_pg_paths
    if isfile(pg_path) && endswith(pg_path, ".arrow")
        append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
    end
end
```
**STATE**: `all_protein_groups` contains merged data from all files

### Step 14C: Data Validation (lines 671-685)
```julia
n_targets = sum(all_protein_groups.target)
n_decoys = sum(.!all_protein_groups.target)
if n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000
    # Proceed with probit
    perform_probit_analysis_multifold(...)
```

### Step 14D: Inside perform_probit_analysis_multifold (lines 1719-1860)

#### 14D.1: CV Fold Detection (lines 1719-1722)
```julia
unique_cv_folds = detect_unique_cv_folds(precursors)
n_folds = length(unique_cv_folds)
@info "Multi-fold probit regression analysis" n_folds=n_folds
```

#### 14D.2: Build Protein-to-CV-Fold Mapping (lines 1756-1762)
```julia
if protein_to_cv_fold === nothing
    psm_paths = [get_corresponding_psm_path(ref) for ref in pg_refs]
    protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, precursors)
end
```
**STATE**: `protein_to_cv_fold` dictionary created

#### 14D.3: Assign CV Folds to Protein Groups (line 1765)
```julia
assign_protein_group_cv_folds!(all_protein_groups, protein_to_cv_fold)
```
**STATE**: `all_protein_groups` now has `:cv_fold` column

#### 14D.4: Train Probit Models (lines 1777-1808)
```julia
models = Dict{UInt8, Vector{Float64}}()
for test_fold in unique_cv_folds
    # Train model for each fold
    β_fitted = fit_probit_model(X_train, y_train)
    models[test_fold] = β_fitted
end
```
**STATE**: `models` dictionary contains trained models per fold

#### 14D.5: Apply Models to All Groups (lines 1810-1849)
```julia
# Apply models to the merged dataframe
for fold in unique(all_protein_groups.cv_fold)
    fold_mask = all_protein_groups.cv_fold .== fold
    if haskey(models, fold) && sum(fold_mask) > 0
        X_fold = Matrix{Float64}(all_protein_groups[fold_mask, feature_names])
        all_protein_groups[fold_mask, :pg_score] = Float32.(calculate_probit_scores(X_fold, models[fold]))
    end
end
```
**STATE**: `all_protein_groups` has updated `:pg_score` values

#### 14D.6: CRITICAL - Update Individual Files (lines 1853-1856)
```julia
if !isempty(pg_refs)
    @info "Updating individual protein group files with probit scores"
    apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names)
end
```

### Step 14E: Inside apply_probit_scores_multifold! (lines 1675-1715)
```julia
for ref in pg_refs
    transform_and_write!(ref) do df
        # Add CV folds
        cv_folds = Vector{UInt8}(undef, nrow(df))
        for (i, protein_name) in enumerate(df.protein_name)
            if haskey(protein_to_cv_fold, protein_name)
                cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
            else
                cv_folds[i] = UInt8(0)
            end
        end
        df[!, :cv_fold] = cv_folds
        
        # Save old scores
        df[!, :old_pg_score] = copy(df.pg_score)
        
        # Apply models
        for (fold, model) in models
            mask = df.cv_fold .== fold
            if sum(mask) > 0
                X = Matrix{Float64}(df[mask, feature_names])
                df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
            end
        end
        
        # Sort
        sort!(df, [:pg_score, :target], rev = [true, true])
        
        # Remove temp column
        select!(df, Not(:cv_fold))
        
        return df
    end
end
```
**RESULT**: Each file updated with new scores, sorted, written back

## Path 2: When Probit Analysis is SKIPPED (Insufficient Data)

### Step 14A-C: Same initial setup and validation

### Step 14D: Skip Branch Executed (lines 687-714, current implementation)
```julia
@info "Skipping Probit analysis: insufficient data (targets: $n_targets, decoys: $n_decoys)"

# Process each file individually
total_pg_processed = 0
for ref in pg_refs
    if exists(ref)
        # Load this specific file's data
        df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
        n_groups = nrow(df)
        
        # Sort by pg_score and target
        sort!(df, [:pg_score, :target], rev = [true, true])
        
        # Write back THIS file's data only
        writeArrow(file_path(ref), df)
        
        total_pg_processed += n_groups
        @info "Processed protein group file" path=file_path(ref) n_groups=n_groups
    end
end
```
**RESULT**: Each file sorted and written back with original scores

## Key Differences Identified

### 1. **Data Merging**
- **Probit**: Merges all files into `all_protein_groups` DataFrame
- **Skip**: Processes files individually, no merge

### 2. **CV Fold Mapping**
- **Probit**: Creates `protein_to_cv_fold` dictionary
- **Skip**: No CV fold mapping created

### 3. **Protein Group Visibility**
- **Probit**: All protein groups are seen together in memory
- **Skip**: Each file processed in isolation

### 4. **Column Additions**
- **Probit**: Temporarily adds `:cv_fold` and `:old_pg_score` columns
- **Skip**: No column additions

### 5. **Score Updates**
- **Probit**: Updates `:pg_score` with model predictions
- **Skip**: Keeps original `:pg_score` values

### 6. **State Changes**
- **Probit**: Creates and uses `models` dictionary, `protein_to_cv_fold` mapping
- **Skip**: No persistent state created

## Step 23: Update PSMs with Scores

### How update_psms_with_probit_scores_refs Works (lines 719-842)

```julia
for paired_ref in paired_files
    # Load protein groups from PAIRED file only
    pg_table = Arrow.Table(file_path(pg_ref))
    
    # Build lookup from THIS file
    pg_score_lookup = Dict{ProteinKey, Tuple{Float32, Float32}}()
    for i in 1:n_pg_rows
        key = ProteinKey(pg_table[:protein_name][i], pg_table[:target][i], pg_table[:entrap_id][i])
        pg_score_lookup[key] = (pg_table[:pg_score][i], pep_val)
    end
    
    # Process PSMs
    transform_and_write!(psm_ref) do psms_df
        for i in 1:n_psms
            key = ProteinKey(psms_df[i, :inferred_protein_group], ...)
            if !haskey(pg_score_lookup, key)
                # Set to missing - THIS IS THE PROBLEM
                probit_pg_scores[i] = missing
            end
        end
    end
end
```

## The Critical Discovery

### The Problem is NOT in Step 14!

Both paths:
1. Process files individually
2. Sort by score
3. Write back to original files

### The Real Issue: File Pairing Assumption

Step 23 assumes protein groups for a PSM file are in the PAIRED protein group file. But:
- Protein inference (Step 12) creates protein groups PER FILE
- A protein appearing in multiple MS files gets separate entries
- Step 23 only looks in the paired file

### Why Probit "Works"

It's not that probit redistributes protein groups - it's that when probit runs:
1. The merged DataFrame operation might be deduplicating protein groups
2. The CV fold mapping might be creating a shared reference
3. The transform_and_write! might be doing something different

### Next Investigation Steps

1. Check if protein groups are duplicated across files
2. Verify if transform_and_write! behaves differently than our manual approach
3. Examine if the CV fold mapping affects protein group distribution
4. Check if there's hidden state in SearchContext

## Hypothesis

The issue might be that we're not using `transform_and_write!` in the skip path, which might have side effects or special handling that our simple `writeArrow` doesn't have.