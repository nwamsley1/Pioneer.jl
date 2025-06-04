# SearchDIA Data Structures Analysis

## Overview

The SearchDIA module in Pioneer.jl implements a sophisticated data-independent acquisition (DIA) proteomics search engine. This analysis examines the core data structures, their relationships, type hierarchies, and design patterns.

## Core Architectural Patterns

### 1. Abstract Type Hierarchies

```julia
# Base interfaces
abstract type PSM end                          # Peptide-spectrum matches
abstract type SearchMethod end                 # Different search algorithms  
abstract type SearchResults end               # Algorithm-specific results
abstract type SearchParameters end            # Algorithm-specific parameters
abstract type SearchDataStructures end        # Working data containers
abstract type MassSpecDataReference end       # MS data file references
abstract type SpectralLibrary end            # Library implementations
abstract type Ion{T<:AbstractFloat} end      # All ion types
abstract type MatchIon{T<:AbstractFloat} <: Ion{T} end  # Matched ions
abstract type LibraryIon{T<:AbstractFloat} <: Ion{T} end # Library ions
```

### 2. Precision-Aware Design

Most data structures are parameterized by floating-point precision (`T<:AbstractFloat`), with separate high-precision (`H`) and low-precision (`L`) fields where memory efficiency matters:

```julia
struct ComplexScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    # High precision for critical measurements
    error::H
    weight::H
    
    # Low precision for scores (Float16 typically)
    poisson::L
    spectral_contrast::L
    # ... other scores
end
```

### 3. Performance-Optimized Collections

Custom data structures prioritize cache efficiency and computational performance:

- **SparseArray**: Column-sparse matrix optimized for spectral fitting
- **ArrayDict**: Fast key-value mapping with preallocated arrays
- **Counter**: Efficient counting with minimal allocations

## Data Structure Categories

### 1. PSM (Peptide-Spectrum Match) Types

#### Abstract PSM Hierarchy

```julia
abstract type PSM end
  ├── abstract type UnscoredPSM{T<:AbstractFloat} <: PSM
  │   ├── SimpleUnscoredPSM{T}     # Basic fragment statistics
  │   ├── ComplexUnscoredPSM{T}    # Detailed fragmentation features
  │   └── Ms1UnscoredPSM{T}        # Precursor isotope features
  └── abstract type ScoredPSM{H,L<:AbstractFloat} <: PSM  
      ├── SimpleScoredPSM{H,L}     # Post-scoring simple features
      ├── ComplexScoredPSM{H,L}    # Post-scoring complex features  
      └── Ms1ScoredPSM{H,L}        # Post-scoring precursor features
```

#### Key PSM Features

**SimpleUnscoredPSM** - Basic fragment matching:
```julia
struct SimpleUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8        # Best fragment rank observed
    topn::UInt8             # Number of top-N fragments matched
    b_count::UInt8          # b-ion count
    y_count::UInt8          # y-ion count  
    p_count::UInt8          # Precursor ion count
    i_count::UInt8          # Isotope count
    intensity::T            # Total matched intensity
    error::T                # Mass error accumulator
    precursor_idx::UInt32   # Link to precursor
end
```

**ComplexUnscoredPSM** - Advanced fragment analysis:
```julia
struct ComplexUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    # Rankings for different ion types
    best_rank::UInt8
    best_rank_iso::UInt8
    topn::UInt8
    topn_iso::UInt8
    
    # Fragmentation coverage
    longest_y::UInt8         # Longest y-ion series
    longest_b::UInt8         # Longest b-ion series
    isotope_count::UInt8
    
    # Ion-specific intensities
    b_int::T                 # Total b-ion intensity
    y_int::T                 # Total y-ion intensity
    
    # Additional counters
    non_cannonical_count::UInt8
    error::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end
```

### 2. Spectral Scoring Types

#### Spectral Distance Metrics

```julia
abstract type SpectralScores{T<:AbstractFloat} end

# For simple searches - basic similarity metrics
struct SpectralScoresSimple{T<:AbstractFloat} <: SpectralScores{T}
    scribe::T                        # Scribe score (normalized spectral angle)
    city_block::T                    # City block distance  
    spectral_contrast::T             # Cosine similarity
    matched_ratio::T                 # Log ratio of matched/unmatched
    entropy_score::T                 # Spectral entropy
    percent_theoretical_ignored::T   # Fraction of theory ignored
end

# For complex searches - deconvolution-based metrics  
struct SpectralScoresComplex{T<:AbstractFloat} <: SpectralScores{T}
    spectral_contrast::T             # Raw spectral contrast
    fitted_spectral_contrast::T     # After deconvolution fitting
    gof::T                          # Goodness of fit (-log2 of residuals)
    max_matched_residual::T         # Worst matched peak residual
    max_unmatched_residual::T       # Worst unmatched peak residual
    fitted_manhattan_distance::T    # Manhattan distance after fitting
    matched_ratio::T                # Log ratio matched/unmatched
end
```

### 3. Ion Types and Fragment Representation

#### Ion Type Hierarchy

```julia
abstract type Ion{T<:AbstractFloat} end
  ├── abstract type LibraryIon{T} <: Ion{T}
  │   ├── LibraryPrecursorIon{T}      # Complete precursor information
  │   └── abstract type LibraryFragmentIon{T} <: LibraryIon{T}
  │       ├── SimpleFrag{T}            # Basic fragment with intensity
  │       ├── abstract type AltimeterFragment{T} <: LibraryFragmentIon{T}
  │       │   ├── DetailedFrag{T}      # Full fragment annotation
  │       │   └── SplineDetailedFrag{N,T}  # Spline-based intensity
  │       ├── PrecursorBinFragment{T}  # For indexing
  │       └── IndexFragment{T}         # Fragment index entry
  └── abstract type MatchIon{T} <: Ion{T}
      ├── FragmentMatch{T}             # Observed fragment match
      └── PrecursorMatch{T}            # Observed precursor match
```

#### Fragment Representations

**DetailedFrag** - Complete fragment annotation:
```julia
struct DetailedFrag{T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32          # Parent precursor ID
    mz::T                    # Fragment m/z
    intensity::Float16       # Predicted intensity
    
    # Ion classification
    ion_type::UInt16         # Numeric ion type identifier
    is_y::Bool              # y-ion flag
    is_b::Bool              # b-ion flag  
    is_p::Bool              # Precursor flag
    is_isotope::Bool        # Isotope flag
    
    # Physical properties
    frag_charge::UInt8      # Fragment charge
    ion_position::UInt8     # Position in peptide
    prec_charge::UInt8      # Precursor charge
    rank::UInt8             # Intensity rank
    sulfur_count::UInt8     # Sulfur atoms (for isotope correction)
end
```

**SplineDetailedFrag** - NCE-dependent intensity modeling:
```julia
struct SplineDetailedFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}
    # Same fields as DetailedFrag except:
    intensity::NTuple{N, T}  # Spline coefficients instead of single value
    # ... other fields identical
end
```

### 4. Search Infrastructure Types

#### SearchContext - Central Coordination

```julia
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    # Core components
    spec_lib::L                              # Spectral library
    temp_structures::AbstractVector{<:SearchDataStructures}  # Per-thread working data
    mass_spec_data_reference::M              # MS file references
    
    # Output management
    data_out_dir::Base.Ref{String}          # Main output directory
    qc_plot_folder::Base.Ref{String}        # QC plots location
    rt_alignment_plot_folder::Base.Ref{String}
    mass_err_plot_folder::Base.Ref{String}
    
    # Calibration models (per MS file)
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    ms1_mass_error_model::Dict{Int64, MassErrorModel}
    nce_model::Dict{Int64, NceModel}
    
    # RT conversion models
    irt_rt_map::Dict{Int64, RtConversionModel}    # iRT to RT
    rt_irt_map::Dict{Int64, RtConversionModel}    # RT to iRT
    
    # Search state
    precursor_dict::Base.Ref{Dictionary}    # Best precursor results
    irt_errors::Dict{Int64, Float32}        # RT calibration errors
    irt_obs::Dict{UInt32, Float32}          # Observed iRT values
    
    # Configuration
    n_threads::Int64
    n_precursors::Int64  
    buffer_size::Int64
end
```

#### SimpleLibrarySearch - Working Data Container

```julia
mutable struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    # Match data (reused across spectra)
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    mass_err_matches::Vector{FragmentMatch{Float32}}
    
    # Fast lookup structures
    id_to_col::ArrayDict{UInt32, UInt16}    # Precursor ID to column mapping
    prec_count::Counter{UInt32, UInt8}      # Precursor match counting
    ion_templates::Vector{DetailedFrag{Float32}}  # Template fragments
    iso_splines::I                          # Isotope correction model
    
    # PSM collections for different search modes
    scored_psms::Vector{SimpleScoredPSM{Float32, Float16}}
    unscored_psms::Vector{SimpleUnscoredPSM{Float32}}
    spectral_scores::Vector{SpectralScoresSimple{Float16}}
    
    complex_scored_psms::Vector{ComplexScoredPSM{Float32, Float16}}
    complex_unscored_psms::Vector{ComplexUnscoredPSM{Float32}}
    complex_spectral_scores::Vector{SpectralScoresComplex{Float16}}
    
    ms1_scored_psms::Vector{Ms1ScoredPSM{Float32, Float16}}
    ms1_unscored_psms::Vector{Ms1UnscoredPSM{Float32}}
    ms1_spectral_scores::Vector{SpectralScoresMs1{Float16}}
    
    # Working arrays for deconvolution
    Hs::SparseArray                         # Design matrix
    prec_ids::Vector{UInt32}               # Precursor IDs
    precursor_weights::Vector{Float32}      # Fitted weights
    temp_weights::Vector{Float32}          # Temporary weights
    residuals::Vector{Float32}             # Fit residuals
    isotopes::Vector{Float32}              # Isotope corrections
    precursor_transmission::Vector{Float32} # Quad transmission
end
```

### 5. Mass Spectrometry Data Types

#### MS Data Hierarchy

```julia
abstract type MassSpecData end
  └── abstract type NonIonMobilityData <: MassSpecData
      ├── BasicNonIonMobilityMassSpecData{T}   # Single-file data
      └── BatchNonIonMobilityMassSpecData{T}   # Multi-file batched data
```

#### ArrowTableReference - File Management

```julia
struct ArrowTableReference{N} <: MassSpecDataReference
    file_paths::NTuple{N, String}           # Paths to Arrow files
    file_id_to_name::NTuple{N, String}      # Parsed file names
    
    # Result storage paths (per file)
    first_pass_psms::Vector{String}         # First pass results
    second_pass_psms::Vector{String}        # Second pass results  
    passing_psms::Vector{String}           # Final passing PSMs
    passing_proteins::Vector{String}       # Final protein results
    rt_index_paths::Vector{String}         # RT index paths
    failed_search_indicator::Vector{Bool}  # Search failure flags
end
```

### 6. Performance-Critical Data Structures

#### SparseArray - Optimized Spectral Fitting

```julia
mutable struct SparseArray{Ti<:Integer,T<:AbstractFloat}
    n_vals::Int64              # Number of non-zero values
    m::Int64                   # Number of rows (peaks)
    n::Int64                   # Number of columns (precursors)
    
    # Sparse matrix in CSC format
    rowval::Vector{Ti}         # Row indices
    colval::Vector{UInt16}     # Column values (for sorting)
    nzval::Vector{T}           # Predicted intensities
    
    # Additional arrays for spectral analysis
    matched::Vector{Bool}      # Which entries are matched
    isotope::Vector{UInt8}     # Isotope annotations
    x::Vector{T}              # Observed intensities
    colptr::Vector{Ti}        # Column pointers (CSC format)
end
```

#### ArrayDict - Fast Lookup

```julia
mutable struct ArrayDict{I<:Unsigned,C<:Real}
    keys::Vector{I}            # Key storage
    vals::Vector{C}            # Value storage  
    size::Int64               # Current size
end
```

#### Counter - Efficient Counting

```julia
mutable struct Counter{I,C<:Unsigned}
    ids::Vector{I}            # Entity IDs
    counts::Vector{C}         # Count values
    size::Int64              # Number of active entries
    matches::Int64           # Number of matches found
end
```

## Data Flow and Relationships

### 1. Search Execution Flow

```
MS Data → SearchContext → SearchMethod → SearchResults
    ↓           ↓              ↓              ↓
Arrow Files  Calibration   Algorithm     Output Files
             Models        Parameters
```

### 2. PSM Processing Pipeline

```
UnscoredPSM → Spectral Analysis → ScoredPSM → FDR Control → Final Results
     ↓              ↓                  ↓            ↓             ↓
Raw Features   Distance Metrics   Full Scores   Q-values    Proteins
```

### 3. Fragment Matching Process

```
LibraryFragments → FragmentIndex → FragmentMatch → PSM Features
       ↓               ↓               ↓              ↓
   Theoretical     Fast Lookup     Observed      Statistics
   Predictions                     Matches
```

## Design Patterns and Performance Considerations

### 1. Memory Efficiency Patterns

- **Precision Stratification**: Different precision for different data types
- **Pre-allocation**: Large working arrays allocated once and reused
- **Compact Representations**: UInt8/UInt16 for categorical data
- **Sparse Storage**: Only store non-zero values in spectral fitting

### 2. Cache-Friendly Design

- **Structure of Arrays**: Related data stored in separate arrays
- **Contiguous Memory**: Vector-based storage for hot loops
- **Minimal Indirection**: Direct indexing where possible

### 3. Type Stability

- **Parameterized Types**: Generic over precision types
- **Abstract Interfaces**: Well-defined method contracts
- **Concrete Implementations**: Avoid Union types in hot paths

### 4. Computational Efficiency

- **Vectorized Operations**: SIMD-friendly data layouts
- **Minimal Allocations**: Reuse working memory
- **Specialized Algorithms**: Custom sorting for sparse data

## Testing and Validation Patterns

### 1. Type-Safe Constructors

Many types provide default constructors for testing:
```julia
SimpleUnscoredPSM{Float32}() = SimpleUnscoredPSM(...)  # Zero-initialized
FragmentMatch{Float32}() = FragmentMatch(...)          # Empty match
```

### 2. Interface Compliance

Abstract types define required methods:
```julia
# PSM interface requirements
ScoreFragmentMatches!(results, matches, ...)
ModifyFeatures!(score, transition, mass, intensity)
Score!(PSMs_dict, unscored_PSMs)
```

### 3. Validation Methods

- Getter methods ensure type safety
- Reset methods for memory reuse
- Range checks in critical paths

## Key Insights for Development

### 1. Extensibility Points

- **SearchMethod**: Add new search algorithms
- **PSM Types**: Define new feature sets  
- **SpectralScores**: Implement new similarity metrics
- **SearchParameters**: Configure algorithm behavior

### 2. Performance Hotspots

- **Fragment Matching**: Core of spectral search
- **Spectral Fitting**: Deconvolution via least squares
- **PSM Scoring**: Feature calculation and machine learning
- **RT Alignment**: Cross-run normalization

### 3. Memory Management

- **Working Memory**: Reused across spectra
- **Result Storage**: Accumulated across files
- **Intermediate Data**: Carefully managed lifecycle

This data structure analysis reveals a sophisticated, performance-oriented design that balances flexibility with computational efficiency. The type system enables both extensibility and optimization, while the data flow patterns support complex multi-pass search strategies typical of modern proteomics workflows.