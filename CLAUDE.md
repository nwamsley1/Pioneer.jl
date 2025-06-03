# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Pioneer.jl is a high-performance data-independent acquisition (DIA) proteomics analysis tool that:
- Builds in silico spectral libraries from protein sequences using Koina API
- Searches mass spectrometry data to identify and quantify peptides/proteins
- Implements advanced algorithms for collision energy calibration, isotope correction, and robust quantification

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

### Key Technical Patterns
- Heavy use of Arrow files for data interchange
- Custom fragment indexing for fast lookup
- Multi-threaded processing with task partitioning
- Spline-based models for RT conversion and quadrupole transmission
- XGBoost integration for PSM rescoring
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

### External Dependencies
- PioneerConverter (.NET) - Required for Thermo RAW file conversion
- Koina API - External service for spectrum prediction (requires internet)
- XGBoost - Machine learning for PSM rescoring