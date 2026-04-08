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
├── 7. PrecursorScoringSearch - Precursor rescoring and FDR control
├── 8. IntegrateChromatogramSearch - Peak integration and quantification
├── 9. ProteinInferenceSearch - Protein-group assignment on integrated passing precursors
├── 10. ProteinScoringSearch - Protein-group scoring and q-values
└── 11. MaxLFQSearch - Label-free quantification
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

Protein group handling spans multiple modules:

**Core Protein Inference** (`src/utils/proteinInference.jl`):
- `infer_proteins()` - Type-safe protein inference using ProteinKey and PeptideKey
- Implements greedy set cover algorithm for minimal protein sets
- Handles complex cases: distinct, differentiable, indistinguishable proteins

**Protein Inference and Scoring Integration** (`ProteinInferenceSearch/*`, `ProteinScoringSearch/*`):
- `run_protein_inference!()` - Annotates passing precursors with inferred protein groups and quant flags
- `apply_inference_to_dataframe()` - Wrapper around core infer_proteins algorithm
- `build_protein_group_tables()` - Aggregates annotated passing precursors into per-run protein-group tables

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

Protein handling now runs in dedicated post-integration steps:
1. **Separate Steps**:
   - `ProteinInferenceSearch` annotates integrated passing precursors with inferred protein groups
   - `ProteinScoringSearch` builds protein-group tables, fits protein probit models, and computes protein-level q-values
2. **Multifold Probit Architecture**:
   - Builds protein CV folds from constituent precursor folds
   - Fits fold-specific probit models on protein-level features
   - Applies held-out models back to each fold
3. **Optional Diagnostics**:
   - `proteinScoring.write_qc_plots`
   - `proteinScoring.log_feature_importance`
4. **Shared Protein-Group Builder**: `build_protein_group_tables()` lives in `ProteinScoringSearch/protein_inference_pipeline.jl`
- **Memory Budget**: In-memory row limits are derived from `optimization.machine_learning.max_in_memory_table_mb`

**Performance Benefits**:
- Keeps protein rescoring bounded by the shared in-memory table budget
- Handles multifold protein rescoring without mixing training and held-out folds
- Reuses the same protein inference outputs for both protein-group tables and protein q-values

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

## Recent Development Updates (2025-01)

### Protein Inference Refactoring
The protein inference system was recently refactored to use type-safe structures:

**Completed Changes**:
- Replaced complex NamedTuples with `ProteinKey` and `PeptideKey` types
- Unified API with single `infer_proteins()` function
- Added comprehensive test suite (89 passing tests)
- Removed legacy functions and conversion utilities
- ~350 lines of code reduction

**Migration Guide**:
- Old: `infer_proteins_typed()` → New: `infer_proteins()`
- Old: Complex NamedTuple keys → New: `ProteinKey` and `PeptideKey` structs
- Old: Multiple function variants → New: Single type-safe function

**Testing**:
- Integration: `SearchDIA("./data/ecoli_test/ecoli_test_params.json")`
- Unit tests: `include("test/UnitTests/test_protein_inference.jl")`

### SearchMethods Refactoring Status
All planned SearchMethods refactoring has been completed:
- ✅ FileReference type hierarchy with PSMFileReference and ProteinGroupFileReference
- ✅ Algorithm wrappers for protein inference and PSM score updates  
- ✅ PrecursorScoringSearch interface using file references
- ✅ Generic heap-based merge supporting N sort keys
- ✅ MaxLFQSearch simplified to use MSData directly

## Module-Specific Documentation

- SearchMethods: See `SearchMethods/CLAUDE.md` for detailed guidance on implementing and modifying search methods
- Data Structures: See `DataStructures_CLAUDE.md` for PSM types, scoring systems, and SearchContext
- Common Utilities: See `CommonSearchUtils/CLAUDE.md` for core algorithms and shared functionality
- Transition Selection: See `CommonSearchUtils/selectTransitions/CLAUDE.md` for transition selection framework
