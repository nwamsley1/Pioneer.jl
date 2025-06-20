# Protein Inference Pipeline API Plan

## Status: 🚧 IN PROGRESS (2025-01)

## Overview

This document outlines the plan to simplify `perform_protein_inference` using a composable Pipeline API while preserving the core `infer_proteins` algorithm from `src/utils/proteinInference.jl`.

## Current Architecture

1. **Core Algorithm**: `infer_proteins` implements the parsimony principle for protein inference
2. **File Processing**: `perform_file_protein_inference` wraps `infer_proteins` for single files
3. **Orchestration**: `perform_protein_inference` processes multiple files and tracks mappings

## Proposed Decomposition

### 1. Pre-Inference Pipeline (Prepare PSMs)

```julia
pre_inference_pipeline = TransformPipeline() |>
    add_peptide_metadata(precursors) |>        # Add sequence, protein info from precursor_idx
    deduplicate_peptides() |>                  # Remove duplicate peptide-protein pairs
    validate_peptide_data()                    # Ensure required columns exist
```

### 2. Inference Operation (Wrap existing algorithm)

```julia
# Apply infer_proteins to prepared data
apply_protein_inference_algorithm(precursors) 
```

### 3. Post-Inference Pipeline (Process results)

```julia
post_inference_pipeline = TransformPipeline() |>
    add_inferred_protein_column() |>           # Add protein assignments to PSMs
    add_quantification_flag() |>               # Add use_for_quant based on inference
    filter_by_min_peptides(min_peptides) |>    # Apply minimum peptide filter
    calculate_protein_scores() |>              # Calculate initial protein scores
    add_protein_features(protein_catalog)      # Add coverage, n_peptides, etc.
```

### 4. New Composable Interface

```julia
function perform_protein_inference_pipeline(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    precursors::LibraryPrecursors,
    protein_catalog;
    min_peptides::Int = 2
)
    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    
    for (idx, psm_ref) in enumerate(psm_refs)
        # Step 1: Prepare PSMs
        prepared_psms = apply_pipeline(psm_ref, pre_inference_pipeline)
        
        # Step 2: Run inference (uses existing infer_proteins internally)
        inference_result = apply_inference_to_dataframe(
            load_dataframe(prepared_psms), 
            precursors
        )
        
        # Step 3: Create protein groups
        pg_path = joinpath(output_folder, "protein_groups_$(idx).arrow")
        pg_ref = create_protein_groups_from_inference(
            prepared_psms,
            inference_result,
            pg_path,
            post_inference_pipeline
        )
        
        # Step 4: Update PSMs with protein assignments
        update_psms_with_inference!(psm_ref, inference_result)
        
        # Track mappings
        push!(pg_refs, pg_ref)
        psm_to_pg_mapping[file_path(psm_ref)] = pg_path
    end
    
    return pg_refs, psm_to_pg_mapping
end
```

## Key Design Elements

### 1. Preserve Core Algorithm
- `infer_proteins` remains completely unchanged
- We only change how data flows to and from it

### 2. New Helper Functions

```julia
# Wrapper that calls infer_proteins with proper data format
function apply_inference_to_dataframe(df::DataFrame, precursors)
    # Extract and format data for infer_proteins
    proteins_vec, peptides_vec = extract_for_inference(df, precursors)
    
    # Call existing algorithm
    raw_result = infer_proteins(proteins_vec, peptides_vec)
    
    # Convert to structured result
    return convert_inference_result(raw_result)
end

# Create protein groups from inference results
function create_protein_groups_from_inference(
    psms_ref::PSMFileReference,
    inference_result::InferenceResult,
    output_path::String,
    post_pipeline::TransformPipeline
)
    # Group PSMs by inferred protein
    protein_groups_df = group_psms_by_protein(
        load_dataframe(psms_ref),
        inference_result
    )
    
    # Apply post-processing pipeline
    processed_groups = apply_pipeline_to_df(protein_groups_df, post_pipeline)
    
    # Write and return reference
    writeArrow(output_path, processed_groups)
    return ProteinGroupFileReference(output_path)
end
```

### 3. Pipeline Operations

Each operation is a pure function that transforms data:

- **`add_peptide_metadata`**: Adds sequence and protein info from library
- **`deduplicate_peptides`**: Removes duplicate peptide entries
- **`add_inferred_protein_column`**: Adds inference results to PSMs
- **`calculate_protein_scores`**: Computes log-sum scores
- **`add_protein_features`**: Adds peptide coverage, etc.

## Benefits

1. **Unchanged Core**: `infer_proteins` algorithm remains exactly as is
2. **Testability**: Each pipeline operation can be unit tested
3. **Clarity**: Data transformations are explicit
4. **Flexibility**: Easy to add/modify pre/post processing steps
5. **Performance**: Still processes files individually
6. **Debugging**: Can inspect data at each pipeline stage

## Implementation Steps

1. Create `protein_inference_pipeline.jl` with new operations
2. Implement pipeline operations as pure functions
3. Create wrappers for `infer_proteins` integration
4. Update ScoringSearch to use new interface
5. Add comprehensive tests for each operation
6. Deprecate old `perform_protein_inference` with migration guide

## Example Usage in ScoringSearch

```julia
# Build protein catalog once
protein_catalog = count_protein_peptides(getPrecursors(getSpecLib(search_context)))

# Perform inference with clear pipeline
pg_refs, mappings = perform_protein_inference_pipeline(
    passing_refs,
    passing_proteins_folder,
    getPrecursors(getSpecLib(search_context)),
    protein_catalog,
    min_peptides = params.min_peptides
)

# Continue with rest of workflow...
```

## Key Files to Modify

1. **New file**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/protein_inference_pipeline.jl`
2. **Update**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`
3. **Deprecate**: Current `perform_protein_inference` in `utils.jl`

## Testing Strategy

1. Unit tests for each pipeline operation
2. Integration test comparing old vs new implementation
3. Performance benchmarks
4. Edge case testing (empty files, single peptide, etc.)

This approach maintains the existing algorithmic behavior while making the data flow and transformations explicit and testable.