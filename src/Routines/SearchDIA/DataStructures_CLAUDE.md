# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with SearchDIA data structures.

## Data Structure Overview

SearchDIA uses a sophisticated type hierarchy designed for high-performance proteomics analysis. Data structures are distributed across several locations but form a cohesive system for representing PSMs (Peptide-Spectrum Matches), search infrastructure, and intermediate computations.

## PSM Type Hierarchy

### Abstract Base Types
```julia
abstract type PSM end
abstract type UnscoredPSM{T<:AbstractFloat} <: PSM end
abstract type ScoredPSM{H,L<:AbstractFloat} <: PSM end
```

### PSM Complexity Levels

**1. SimpleUnscoredPSM/SimpleScoredPSM** - Basic fragment statistics
- **Purpose**: Fast initial screening with minimal features
- **Features**: Fragment counts (b/y/p/i), basic ranks, intensity, error
- **Use Cases**: FirstPassSearch, simple filtering

**2. ComplexUnscoredPSM/ComplexScoredPSM** - Detailed fragmentation analysis  
- **Purpose**: Comprehensive scoring with isotope and residual analysis
- **Features**: Isotope tracking, longest ion series, deconvolution metrics
- **Use Cases**: SecondPassSearch, ML training, detailed scoring

**3. Ms1UnscoredPSM/Ms1ScoredPSM** - Precursor-level features
- **Purpose**: MS1-based identification and quantification
- **Features**: Isotope envelope analysis, precursor mass accuracy
- **Use Cases**: Precursor-only searches, MS1 quantification

### Precision Stratification

Scored PSMs use dual precision to optimize memory:
- **H (High Precision)**: Critical values requiring accuracy (errors, weights) 
- **L (Low Precision)**: Features for ML training (typically Float16)

Example:
```julia
struct SimpleScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    error::H                    # High precision for accurate error
    spectral_contrast::L        # Low precision sufficient for ML features
    # ...
end
```

## Match Ion System

### Ion Type Hierarchy
```julia
abstract type MatchIon{T<:AbstractFloat} <: Ion{T} end
```

**FragmentMatch** - Fragment ion matches
- Links predicted fragments to observed peaks
- Contains intensity, mass error, isotope information
- Includes ranking information for scoring

**PrecursorMatch** - Precursor ion matches  
- Isotope envelope matching
- Mass accuracy for MS1 features

### Key Design Patterns

**1. Compact Encoding**: Uses UInt8 for counts, ranks, types to minimize memory
**2. Pre-allocated Arrays**: Working arrays reused across scans to avoid allocations
**3. Sparse Representations**: SparseArray for efficient spectral data storage

## Search Infrastructure

### SearchContext - Central State Manager
```julia
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    spec_lib::L                 # Spectral library
    temp_structures::Vector     # Per-thread working data
    mass_spec_data_reference::M # MS data access
    
    # Calibrated models (stored per MS file)
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
    
    # Global parameters
    huber_delta::Ref{Float32}
    n_threads::Int64
    # ...
end
```

### Thread-Safe Data Structures

**SimpleLibrarySearch** - Per-thread working memory
```julia
mutable struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    # Pre-allocated match arrays
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    
    # PSM storage for all three types
    scored_psms::Vector{SimpleScoredPSM{Float32, Float16}}
    complex_scored_psms::Vector{ComplexScoredPSM{Float32, Float16}}
    ms1_scored_psms::Vector{Ms1ScoredPSM{Float32, Float16}}
    
    # Working arrays for computations
    Hs::SparseArray                    # Sparse spectral representation
    precursor_weights::Vector{Float32} # Deconvolution weights
    residuals::Vector{Float32}         # Spectral residuals
    # ...
end
```

## Spectral Scoring System

### Score Type Hierarchy
```julia
abstract type SpectralScores{T<:AbstractFloat} end
```

**SpectralScoresSimple** - Basic spectral similarity
- scribe, city_block, spectral_contrast
- entropy_score, matched_ratio
- Used with SimpleScoredPSM

**SpectralScoresComplex** - Advanced deconvolution metrics
- fitted_spectral_contrast, goodness-of-fit
- max_matched_residual, fitted_manhattan_distance  
- Used with ComplexScoredPSM

**SpectralScoresMs1** - Precursor-level scoring
- Similar to complex but for MS1 analysis
- Used with Ms1ScoredPSM

### Scoring Philosophy

1. **Iterative Improvement**: Scores computed by iteratively removing worst-fitting peaks
2. **Multiple Metrics**: Different distance metrics capture different aspects of similarity
3. **Deconvolution**: Complex scoring includes spectral deconvolution and residual analysis

## Data Access Patterns

### Mass Spec Data Reference
```julia
abstract type MassSpecDataReference end

struct ArrowTableReference{N} <: MassSpecDataReference
    file_paths::NTuple{N, String}           # MS data files
    first_pass_psms::Vector{String}         # PSM output paths
    second_pass_psms::Vector{String}
    passing_psms::Vector{String}
    passing_proteins::Vector{String}
    rt_index_paths::Vector{String}          # RT indexing
end
```

### File Organization Pattern
- **Input**: Raw MS data converted to Arrow format
- **Intermediate**: PSMs/proteins from each search method  
- **Output**: Final filtered results and protein groups

## Performance Optimizations

### Memory Management
1. **Pre-allocation**: Working arrays sized once per thread
2. **Block Growth**: PSM vectors grow in blocks to minimize allocations
3. **Sparse Storage**: SparseArray for efficient spectral representation
4. **Precision Stratification**: Mixed precision to reduce memory usage

### Threading Patterns
```julia
function processChunk!(method::SearchMethod, batch_id::Int, thread_id::Int)
    # Get thread-specific data structures
    search_data = getSearchData(method.search_context)[thread_id]
    
    # Process assigned scans using pre-allocated arrays
    for scan_idx in getAssignedScans(batch_id)
        # Reuse arrays from search_data
        ScoreFragmentMatches!(search_data.unscored_psms, ...)
    end
end
```

### Custom Collections
- **ArrayDict**: Fast integer key → value mapping
- **Counter**: Efficient counting with integer overflow protection
- **SparseArray**: Memory-efficient spectral storage

## Development Patterns

### Adding New PSM Types
1. **Define unscored version** inheriting from UnscoredPSM{T}
2. **Define scored version** inheriting from ScoredPSM{H,L}  
3. **Implement ModifyFeatures!** for accumulating features
4. **Implement Score!** for conversion to scored PSM
5. **Add to SimpleLibrarySearch** working arrays
6. **Update interface methods** (getters/setters)

### Extending Spectral Scoring
1. **Define new SpectralScores subtype** with required metrics
2. **Implement scoring function** following iterative improvement pattern
3. **Add to corresponding PSM Score! method**
4. **Update getDistanceMetrics** if needed

### Testing Data Structures
```julia
# Test PSM type conversion
unscored = SimpleUnscoredPSM{Float32}()
unscored = ModifyFeatures!(unscored, prec_id, match, error_model, rank)

# Test scoring pipeline
scored_psms = Vector{SimpleScoredPSM{Float32,Float16}}()
Score!(scored_psms, unscored_psms, spectral_scores, ...)

# Validate structure integrity
@test length(scored_psms) == expected_count
@test all(psm -> psm.precursor_idx > 0, scored_psms)
```

## Common Debugging Patterns

### PSM Validation
```julia
# Check for valid precursor indices
@assert all(psm -> psm.precursor_idx > 0, psms)

# Validate scoring pipeline  
@assert length(unscored_psms) == length(spectral_scores)

# Check thread safety
@assert length(search_data.ion_matches) >= max_expected_matches
```

### Memory Usage Analysis
```julia
# Monitor PSM vector growth
println("PSM vector length: $(length(scored_psms))")
println("PSM vector capacity: $(length(scored_psms.ref.mem))")

# Check working array sizes
println("Ion matches allocated: $(length(search_data.ion_matches))")
println("Spectral scores allocated: $(length(search_data.spectral_scores))")
```

## File Locations

- **PSMs/**: PSM.jl, UnscoredPSMs.jl, ScoredPSMs.jl, spectralDistanceMetrics.jl
- **SearchMethods/SearchTypes.jl**: SearchContext, search infrastructure
- **src/structs/**: Shared ion types (MatchIon.jl, LibraryIon.jl, etc.)
- **src/structs/**: Specialized collections (ArrayDict.jl, SparseArray.jl, etc.)

## Integration Points

### With SearchMethods
- SearchContext coordinates all data structures
- Each search method processes PSMs through common interface
- Results flow between methods via SearchContext state

### With CommonUtils  
- matchPeaks.jl generates FragmentMatch/PrecursorMatch objects
- buildDesignMatrix.jl consumes ScoredPSM data
- Fragment indexing works with LibraryIon structures

### With ML Pipeline
- ComplexScoredPSM features feed EvoTrees training
- Protein group scoring uses aggregated PSM scores
- Cross-validation assignments preserved across data structures