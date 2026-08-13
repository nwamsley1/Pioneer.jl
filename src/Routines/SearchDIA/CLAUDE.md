# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the SearchDIA module.

## SearchDIA Overview

SearchDIA is the core search engine of Pioneer.jl that performs data-independent acquisition (DIA) proteomics analysis. It implements a multi-stage search strategy with 7 active search methods, each building upon previous results to progressively refine peptide identification and quantification.

## Architecture

### Search Pipeline Flow

```
SearchDIA.jl (main entry)
├── 1. ParameterTuningSearch    - Establishes fragment mass tolerance
├── 2. NceTuningSearch          - Calibrates collision energy models
├── 3. QuadTuningSearch         - Models quadrupole transmission
├── 4. FragmentIndexSearch      - Builds fragment index for candidate selection
├── 5. MainSearch               - PSM identification, LightGBM scoring, global FDR, fold-split output
├── 6. ScoringSearch            - Rescoring with additional features, FDR control, protein grouping
├── 7. IntegrateChromatogramSearch - Peak integration and quantification
└── 8. MaxLFQSearch             - Label-free quantification
```

**Note**: The old SecondPassSearch is fully replaced by MainSearch as of 2025-03. MainSearch
performs global prescore aggregation, filters to passing precursors, and writes CV fold-split
Arrow files directly to `temp_data/main_search_psms/`. ScoringSearch reads these files and
trains a second LightGBM with additional features (20 features vs MainSearch's 15).

### Key Abstractions

Each search method implements this interface:
```julia
SearchMethod (abstract type)
├── getSearchMethodName()     # Returns method identifier
├── getParameters()           # Returns SearchParameters subtype
├── getResults()              # Returns SearchResults subtype
├── performSearch!()          # Main search logic
├── summarizeResults!()       # Post-processing and output
└── processChunk!()           # Parallel chunk processing
```

### Core Data Structures

**SearchContext** - Central state container:
- `spec_lib`: Spectral library
- `ms_data`: Mass spec data iterator
- `params`: Search parameters
- `models`: Calibrated models (RT, mass error, quad transmission)
- `results`: Accumulated results from all search methods

**PSM Types**:
- `SimpleScoredPSM` - Basic PSM with score
- `ComplexScoredPSM` - PSM with detailed features for ML scoring
- `Ms1ScoredPSM` - PSM with MS1-level quantification

## Common Development Tasks

### Adding a New Search Method

1. Create new folder in `SearchMethods/`
2. Define types:
```julia
struct MySearchMethod <: SearchMethod
    search_context::SearchContext
end

struct MySearchParameters <: SearchParameters
    # parameters...
end

struct MySearchResults <: SearchResults
    # results...
end
```

3. Implement required interface methods
4. Add to `SearchMethods.jl` imports
5. Add to `routines` array in `SearchDIA.jl`

### Modifying Protein Group Scoring

Protein group handling spans multiple modules:

**Core Protein Inference** (`src/utils/proteinInference.jl`):
- `infer_proteins()` - Type-safe protein inference using ProteinKey and PeptideKey
- Implements greedy set cover algorithm for minimal protein sets
- Handles complex cases: distinct, differentiable, indistinguishable proteins

**ScoringSearch Integration** (`ScoringSearch/protein_inference_pipeline.jl`):
- `perform_protein_inference_pipeline()` - Composable pipeline for protein inference
- `apply_inference_to_dataframe()` - Wrapper around core infer_proteins algorithm
- `group_psms_by_protein()` - Aggregates PSMs into protein groups

**Legacy Functions** (`ScoringSearch/utils.jl`):
- `getProteinGroupsDict()` - Creates protein groups from PSMs (legacy approach)
- `merge_sorted_protein_groups()` - Memory-efficient merging using heap
- `writeProteinGroups()` - Outputs final protein groups

**Type-Safe Structures**:
```julia
struct ProteinKey
    name::String
    is_target::Bool
    entrap_id::UInt8
end

struct PeptideKey
    sequence::String
    is_target::Bool
    entrap_id::UInt8
end

struct InferenceResult
    peptide_to_protein::Dictionary{PeptideKey, ProteinKey}
    use_for_quant::Dictionary{PeptideKey, Bool}
end
```

The `merge_sorted_protein_groups` function is critical for handling large datasets:
```julia
# Uses heap-based merging to combine sorted protein group tables
# Preserves sort order while minimizing memory usage
# Handles entrapment groups correctly
```

### Working with ML Protein Scoring

ML scoring integration in `utils_protein_ml.jl` (memory-efficient OOM approach):
1. **Memory-Efficient Architecture**: Uses out-of-memory (OOM) processing to handle large datasets
2. **Sample-Train-Apply Workflow**:
   - Samples protein groups proportionally from files for training
   - Trains LightGBM models on the sample
   - Applies models file-by-file to avoid loading all data into memory
3. **CV Fold Consistency**: Maintains cross-validation fold assignment based on constituent peptides
4. **Feature Engineering**: Extracts top-N precursor scores plus protein-level statistics
5. **Scalability**: Designed for experiments with hundreds to thousands of files

**Key Functions**:
- `apply_ml_protein_scoring_oom!()` - Main OOM workflow orchestrator (in `proteinGroupScoringOOM.jl`)
- `get_protein_groups_with_ml()` - Integration function that calls OOM approach when enabled
- Memory usage is constant regardless of dataset size

**Performance Benefits**:
- Constant memory usage regardless of dataset size
- Handles experiments with thousands of MS files
- Maintains statistical power through strategic sampling
- Follows established OOM patterns from PSM scoring in `percolatorSortOf.jl`

### Debugging Search Results

Key output files to examine:
- `precursors_{method}.arrow` - PSMs from each method
- `passing_proteins_{method}.arrow` - Protein groups
- `models.jld2` - Calibrated models
- `psm_plots/` - Diagnostic plots

## Performance Considerations

### Thread Task Partitioning
`partitionThreadTasks.jl` divides work based on:
- MS2 scan density
- Available threads
- Memory constraints

### Memory-Efficient Operations
- Arrow file streaming for large datasets
- Heap-based merging in `merge_sorted_protein_groups`
- Sparse matrix operations in quantification

### RT Indexing
`buildRTIndex.jl` creates efficient lookup structures:
- Binary search for RT windows
- Pre-computed scan ranges
- Minimizes chromatogram extraction overhead

## Testing SearchDIA

### Integration Test
```julia
# Full pipeline test (search-only, runs against the committed ecoli library)
SearchDIA("./test/integration/search_ecoli.json")
```

### Unit Testing Individual Methods
```julia
# Test specific search method
include("test/UnitTests/SearchMethods/test_scoring_search.jl")
```

### Debugging Tips
1. Enable verbose logging in parameters
2. Check intermediate Arrow files
3. Use `@show` macros in `processChunk!` methods
4. Examine `models.jld2` for calibration issues

## Common Pitfalls

1. **Memory Usage**: Watch for loading entire datasets into memory
2. **Thread Safety**: SearchContext is shared - avoid mutations in parallel code
3. **Sort Order**: Many operations assume sorted data (by scan, RT, or score)
4. **Entrapment Groups**: Always consider entrapment when grouping proteins

## Key Files Reference

- `SearchDIA.jl` - Main orchestration logic (top-level pipeline driver)
- `LibrarySearch.jl` - Fragment-index search wrapper around `searchFragmentIndexPartitionMajorHinted`
- `process_scans_fused.jl` - Fused per-precursor scan loop dispatch
- `PSMs/` - PSM data structures and scoring metrics
- `CommonSearchUtils/fusedMatch.jl` - Fused match+score+build (`run_fused!`); replaces classic
  `selectTransitions! + matchPeaks! + buildDesignMatrix! + sortSparse! + ScoreFragmentMatches!`
- `CommonSearchUtils/fusedScan.jl` - `FusedScratch`, `SparseArrayFused`, `finalize_column!`
- `WriteOutputs/` - Result formatting and plotting

## Recent Development Updates

### MainSearch Consolidation (2025-03)
The old FirstPassSearch and SecondPassSearch have been consolidated into MainSearch.
The pipeline no longer performs a second fragment-index search. MainSearch handles the
complete flow:

1. **Fragment matching + spectral deconvolution** (`deconvolve_spectra` / `deconvolve_scans!`) → all (scan, precursor) PSMs
2. **Per-file LightGBM scoring** (15 features) → `temp_data/prescore_scores/{file}.arrow`
3. **Global prescore aggregation** — PEP-calibrate each file's scores, log-odds combine
   top-sqrt(n) across files, compute global q-values
4. **Filter to passing precursors** at `search.global_prescore_qvalue_threshold` (default 2%)
5. **Write fold-split Arrow files** to `temp_data/main_search_psms/` with additional
   columns (irt_diff, prec_mz, pair_id, etc.) via `add_search_columns!`
6. **Filter fragment index** to passing precursors for downstream chromatogram integration

ScoringSearch then reads the fold-split files and trains a second LightGBM with 20 features
(adds log2_intensity_explained, prec_mz, irt_pred, irt_diff, longest_y, b_count, charge).

**Key functions**:
- `aggregate_prescore_globally!()` — `MainSearch/utils.jl`
- `calibrate_file_scores()` / `combine_scores()` — `MainSearch/prescore_aggregation.jl`
- `_logodds_combine()` — log-odds averaging with floor clamping

### Protein Inference Refactoring (2025-01)
- Replaced complex NamedTuples with `ProteinKey` and `PeptideKey` types
- Unified API with single `infer_proteins()` function
- ~350 lines of code reduction

## Module-Specific Documentation

- SearchMethods: See `SearchMethods/CLAUDE.md` for detailed guidance on implementing and modifying search methods
- Data Structures: See `DataStructures_CLAUDE.md` for PSM types, scoring systems, and SearchContext
- Common Utilities: See `CommonSearchUtils/CLAUDE.md` for core algorithms and shared functionality
