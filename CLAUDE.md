# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Pioneer.jl is a high-performance data-independent acquisition (DIA) proteomics analysis tool that:
- Builds in silico spectral libraries from protein sequences using Koina API
- Searches mass spectrometry data to identify and quantify peptides/proteins
- Implements advanced algorithms for collision energy calibration, isotope correction, and robust quantification

Version: 0.1.12

## Essential Commands

### Development Setup
```bash
# Start Julia with optimal performance (adjust threads based on system)
julia --threads 15 --gcthreads 7,1

# In Julia REPL
pkg> activate .
pkg> develop ./
```

### Running Tests
```bash
# Run all tests
pkg> test

# Run specific test file
julia> include("test/UnitTests/buildDesignMatrix.jl")

# Run integration test
julia> SearchDIA("./data/ecoli_test/ecoli_test_params.json")
```

### Building Documentation
```bash
cd docs
julia --project=. make.jl
# Output in docs/build/
```

### Core Functionality
```julia
using Pioneer

# Build spectral library
params = GetBuildLibParams(output_dir, lib_name, fasta_dir)
BuildSpecLib(params)

# Search DIA data
params = GetSearchParams(lib_path, ms_data_path, output_dir)
SearchDIA(params)

# Parse empirical library
ParseSpecLib("path/to/params.json")

# Convert mzML files
convertMzML("path/to/params.json")
```

## Architecture Overview

### Key Module Structure
- `src/Routines/BuildSpecLib/` - Library construction pipeline
  - `fasta/` - Protein digestion and peptide generation
  - `koina/` - API integration for spectrum prediction
  - `fragments/` - Fragment indexing and processing
  - `chronologer/` - Retention time prediction

- `src/Routines/SearchDIA/` - DIA search pipeline
  - `SearchMethods/` - Multiple search strategies (FirstPass, Scoring, MaxLFQ, etc.)
  - `PSMs/` - Peptide-spectrum match handling
  - `CommonSearchUtils/` - Shared utilities for RT indexing, peak matching

### Critical Data Flow
1. FASTA → Enzymatic digestion → Peptides with modifications
2. Peptides → Koina models (Altimeter/Prosit/UniSpec) → Predicted spectra
3. Raw MS data → PioneerConverter → Arrow format files
4. Arrow data + Library → Multi-pass search → PSMs → FDR control → Protein quantification

## SearchDIA Deep Dive

### Search Pipeline Architecture

SearchDIA executes a series of search methods in a specific order, each building on the results of previous steps:

1. **ParameterTuningSearch** - Establishes fragment mass tolerance
2. **NceTuningSearch** - Calibrates collision energy models
3. **QuadTuningSearch** - Models quadrupole transmission
4. **FirstPassSearch** - Initial PSM identification with RT calibration
5. **HuberTuningSearch** - Optimizes Huber loss parameters
6. **SecondPassSearch** - Refined search with calibrated parameters
7. **ScoringSearch** - FDR control and protein group scoring
8. **IntegrateChromatogramSearch** - Peak integration and quantification
9. **MaxLFQSearch** - Label-free quantification

### SearchMethod Interface

All search methods implement this interface:
- `get_parameters(::SearchMethod, params::PioneerParameters)` - Extract relevant parameters
- `init_search_results(::SearchParameters, ::SearchContext)` - Initialize result containers
- `process_file!(::SearchResults, ::SearchParameters, ::SearchContext, ms_file_idx, spectra)` - Process single file
- `process_search_results!(...)` - Post-process results (optional)
- `reset_results!(::SearchResults)` - Clean up between files
- `summarize_results!(...)` - Aggregate results across all files

### Key Data Structures

**SearchContext** - Central state container:
- Spectral library reference
- Thread-local search data structures
- Model storage (RT, mass error, NCE, quadrupole)
- Output paths and QC folders
- Precursor dictionary for match-between-runs

**PSM Types**:
- `SimpleScoredPSM` - Basic scored match (first/second pass)
- `ComplexScoredPSM` - Full feature PSM (scoring search)
- `Ms1ScoredPSM` - MS1-based quantification

**Search Results Types**:
- Each SearchMethod defines its own results structure
- Common pattern: store file paths, models, and aggregated statistics

### Protein Group Handling

**Protein Inference** (`utils/proteinInference.jl`):
- Implements parsimony principle via greedy set cover
- Handles indistinguishable proteins and shared peptides
- Returns dictionary mapping peptides to protein groups with retention flags
- **Entrapment Groups**: Proteins from different entrapment groups are treated as separate entities
  - Each protein group key includes `(protein_name, target, entrap_id)`
  - Ensures proper FDR estimation across different sample preparations

**Standard Scoring** (`getProteinGroupsDict` in utils.jl`):
- Log-sum scoring: `score = -sum(log(1 - prob))` for all peptides
- Filters by minimum peptide count (default: 2)
- Tracks peptide sets per protein group
- Respects entrapment group boundaries in scoring

**ML-Enhanced Scoring** (`utils_protein_ml.jl`):
- Extracts top-N precursor scores as features
- Adds protein-level statistics (peptide count, score distribution)
- Uses XGBoost random forests with cross-validation
- Maintains CV fold consistency with precursor models

**Merging Functions**:
- `merge_sorted_psms_scores` - Heap-based merging of PSM files
- `merge_sorted_protein_groups` - Similar for protein group files
- Both handle large datasets via batched processing

### Common Utilities

**Fragment Matching** (`matchPeaks.jl`):
- Core spectral matching algorithm
- Handles isotope patterns and mass errors
- Returns detailed match statistics

**RT Indexing** (`buildRTIndex.jl`):
- Creates retention time indices for fast chromatogram extraction
- Supports isotope trace selection
- Enables efficient XIC generation

**Design Matrix** (`buildDesignMatrix.jl`):
- Sparse matrix construction for deconvolution
- Handles overlapping isotope patterns
- Critical for accurate quantification

**Transition Selection**:
- Multiple strategies based on search phase
- `StandardTransitionSelection` - Fragment-based
- `RTIndexTransitionSelection` - Chromatogram-based
- `QuadEstimationSelection` - For transmission modeling

### Testing Patterns

**Integration Test**:
```julia
SearchDIA("./data/ecoli_test/ecoli_test_params.json")
```
- Tests full pipeline with real data
- Validates output files and folder structure
- Checks protein group generation

**Unit Test Examples**:
- `buildDesignMatrix.jl` - Tests sparse matrix construction
- `matchPeaks.jl` - Validates fragment matching
- `proteinInference.jl` - Tests all protein grouping cases

### Output Structure

SearchDIA creates:
```
results_dir/
├── pioneer_search_log.txt
├── config.json (copy of parameters)
├── temp_data/ (intermediate files)
├── qc_plots/
│   ├── rt_alignment_plots/
│   ├── mass_error_plots/
│   └── (other QC plots)
├── precursors_long.[arrow|tsv]
├── precursors_wide.[arrow|tsv]
├── protein_groups_long.[arrow|tsv]
└── protein_groups_wide.[arrow|tsv]
```

### Debugging SearchDIA

1. **Check search progress**: Monitor `pioneer_search_log.txt`
2. **Verify intermediate files**: Look in `temp_data/` folders
3. **Inspect models**: Check SearchContext model dictionaries
4. **Trace PSM flow**: Follow PSM files through each search phase
5. **Protein issues**: Check protein inference dictionary and merge operations

### Key Technical Patterns
- Heavy use of Arrow files for data interchange
- Custom fragment indexing for fast lookup
- Multi-threaded processing with task partitioning
- Spline-based models for RT conversion and quadrupole transmission
- XGBoost integration for PSM rescoring and protein group ML scoring
- Huber loss optimization for robust parameter estimation

### Performance Considerations
- Fragment index uses binary search with pre-computed bounds
- RT indexing enables efficient chromatogram extraction
- Design matrix construction optimized for sparse operations
- Thread task partitioning based on data characteristics

## Development Notes

### Parameter Files
All major functions use JSON parameter files. Examples in `data/example_config/`:
- `defaultBuildLibParams.json` - Library building parameters
- `defaultSearchParams.json` - DIA search parameters
- Key parameters control tolerances, scoring, and output formats

### Testing Strategy
- Unit tests in `test/UnitTests/` cover individual components
- Integration test uses E. coli dataset in `data/ecoli_test/`
- CI runs on Julia 1.11 via GitHub Actions

### Common Debugging Entry Points
- `src/Routines/SearchDIA/searchRAW.jl` - Main search orchestration
- `src/Routines/BuildSpecLib/build/buildPioneerLib.jl` - Library building pipeline
- `src/Routines/SearchDIA/CommonSearchUtils/matchPeaks.jl` - Core peak matching logic
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils.jl` - Protein group scoring (getProteinGroupsDict, writeProteinGroups)

### External Dependencies
- PioneerConverter (.NET) - Required for Thermo RAW file conversion
- Koina API - External service for spectrum prediction (requires internet)
- XGBoost - Machine learning for PSM rescoring

### Supported Prediction Models
- **unispec** - Instrument-specific models for QE, QEHFX, LUMOS, ELITE, VELOS
- **altimeter** - Spline coefficient model (requires local deployment at port 8000)
- **prosit_2020_hcd** - Instrument-agnostic HCD fragmentation model
- **AlphaPeptDeep** - Supports QE, LUMOS, TIMSTOF, SCIEXTOF instruments

### Current Development Focus
- Refactoring ScoringSearch and MaxLFQSearch for better encapsulation
- See src/Routines/SearchDIA/SearchMethods/REFACTORING_PLAN_V2.md for detailed plan
- Adding FileReference abstraction layer for type-safe file operations
- Protein group ML scoring integration using XGBoost random forests
- See PROTEIN_ML_SCORING_INTEGRATION.md for implementation details
- Features top-N precursor scores with cross-validation consistency

### Protein Group ML Scoring

Pioneer now supports ML-enhanced protein group scoring that uses top N precursor scores as features:

**Features**:
- XGBoost random forest models trained per CV fold
- Uses top N precursor scores plus protein-level statistics as features
- Out-of-memory (OOM) support for large experiments
- Automatic fallback to standard scoring if ML fails
- Maintains CV fold consistency with precursor scoring

**Configuration** (in parameters JSON):
```json
"machine_learning": {
    "use_ml_protein_scoring": true,
    "n_top_precursors": 5,
    "num_protein_trees": 100
}
```

**Debugging**:
```julia
# Test protein scoring on a single file
include("src/utils/ML/test_protein_scoring.jl")
test_protein_scoring_pipeline(psm_path, pg_path, precursors)

# Diagnose issues across multiple files
diagnose_protein_scoring_issue(passing_psms_paths, passing_pg_paths, precursors)
```

**Implementation Files**:
- Core ML: `src/utils/ML/proteinGroupScoring.jl`
- Integration: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/utils_protein_ml.jl`
- Testing: `src/utils/ML/test_protein_scoring.jl`

### Module-Specific Documentation
- SearchDIA module: See `src/Routines/SearchDIA/CLAUDE.md` for deep dive into search pipeline
- SearchMethods: See `src/Routines/SearchDIA/SearchMethods/CLAUDE.md` for implementing search methods
- SearchDIA Data Structures: See `src/Routines/SearchDIA/DataStructures_CLAUDE.md` for PSM types and scoring systems
- Common Search Utilities: See `src/Routines/SearchDIA/CommonSearchUtils/CLAUDE.md` for core algorithms
- Transition Selection: See `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/CLAUDE.md` for transition selection

## Memories
- remember me as memory update

## Current Development: SearchMethods Refactoring (2025-01)

### Overview
Refactoring ScoringSearch and MaxLFQSearch to improve encapsulation using file references instead of direct file access. See REFACTORING_PLAN_V2.md for detailed plan.

### Completed Work
1. **Phase 1**: Created abstract FileReference type hierarchy
   - FileReferences.jl: Abstract base type with PSMFileReference and ProteinGroupFileReference
   - SearchResultReferences.jl: Type-safe result management
   - FileOperations.jl: Updated to use abstract types

2. **Phase 2**: Added algorithm wrappers
   - apply_protein_inference: Wraps getProteinGroupsDict from utils
   - update_psms_with_scores: Streaming PSM score updates

3. **Phase 3**: Created ScoringSearch interface
   - scoring_interface.jl: Reference-only functions for ScoringSearch
   - All file access through FileOperations layer
   - Comprehensive unit tests

4. **Phase 4.1**: Added method results storage to SearchContext
   - Added method_results field as Dict{Type{<:SearchMethod}, Any}
   - Added store_results!, get_results, has_results accessor functions
   - Updated SearchContext constructor

5. **Phase 4.3**: Implemented generic heap-based merge in FileOperations
   - Generic heap-based merge supporting arbitrary number of sort keys
   - Dynamic heap type construction based on column types
   - Memory-efficient batch processing
   - Specialized merge functions: merge_psm_scores, merge_protein_groups_by_score

### Test Files Created
- test/UnitTests/test_file_references.jl
- test/UnitTests/test_result_references.jl  
- test/UnitTests/test_scoring_interface.jl

### In Progress
- Phase 4.2: Update ScoringSearch to use new reference-based approach
- Phase 4.4: Update MaxLFQSearch to use references from SearchContext

### Next Steps
1. Integrate reference creation in ScoringSearch after protein inference
2. Store ScoringSearchResultRefs in SearchContext
3. Update MaxLFQSearch to retrieve and use references
4. Clean up scoring_reference_wrapper.jl (untracked file)

### Important Notes
- The protein inference and MaxLFQ algorithms already exist and are NOT being reimplemented
- We are only wrapping them with better abstractions and safety checks
- All tests passing for completed phases
- Maintains backward compatibility while adding type safety
- All tests passing for completed phases