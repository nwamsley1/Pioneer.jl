# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the SearchDIA module.

## SearchDIA Overview

SearchDIA is the core search engine of Pioneer.jl that performs data-independent acquisition (DIA) proteomics analysis. It implements a multi-pass search strategy with 9 sequential search methods, each building upon previous results to progressively refine peptide identification and quantification.

## Architecture

### Search Pipeline Flow

```
SearchDIA.jl (main entry)
├── 1. ParameterTuningSearch - Establishes fragment mass tolerance
├── 2. NceTuningSearch - Calibrates collision energy models
├── 3. QuadTuningSearch - Models quadrupole transmission
├── 4. FirstPassSearch - Initial PSM identification with RT calibration
├── 5. HuberTuningSearch - Optimizes Huber loss parameters
├── 6. SecondPassSearch - Refined search with calibrated parameters
├── 7. ScoringSearch - FDR control and protein group scoring
├── 8. IntegrateChromatogramSearch - Peak integration and quantification
└── 9. MaxLFQSearch - Label-free quantification
```

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

Protein group handling is in `ScoringSearch/utils.jl`:

**Key Functions**:
- `getProteinGroupsDict()` - Creates protein groups from PSMs
- `merge_sorted_protein_groups()` - Memory-efficient merging using heap
- `writeProteinGroups()` - Outputs final protein groups

**Protein Group Key Structure**:
```julia
(protein_name::String, is_target::Bool, entrap_id::UInt8)
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
   - Trains XGBoost random forest models on the sample
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
# Full pipeline test
SearchDIA("./data/ecoli_test/ecoli_test_params.json")
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

- `searchRAW.jl` - Main orchestration logic
- `PSMs/` - PSM data structures and scoring metrics
- `CommonSearchUtils/matchPeaks.jl` - Core spectral matching
- `CommonSearchUtils/selectTransitions/` - Transition selection strategies
- `WriteOutputs/` - Result formatting and plotting

## Module-Specific Documentation

- SearchMethods: See `SearchMethods/CLAUDE.md` for detailed guidance on implementing and modifying search methods
- Data Structures: See `DataStructures_CLAUDE.md` for PSM types, scoring systems, and SearchContext
- Common Utilities: See `CommonSearchUtils/CLAUDE.md` for core algorithms and shared functionality
- Transition Selection: See `CommonSearchUtils/selectTransitions/CLAUDE.md` for transition selection framework