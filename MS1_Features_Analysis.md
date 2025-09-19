# Pioneer.jl MS1 Feature Scoring and Mass Tolerance Estimation Analysis

## Overview

This document provides a comprehensive analysis of how MS1 (precursor-level) scoring features work in Pioneer.jl and how MS1 mass tolerance is estimated. Pioneer.jl implements a sophisticated MS1 feature scoring system that operates alongside MS2 fragment-based scoring to provide enhanced peptide identification confidence.

## MS1 Feature Architecture

### Core Data Structures

Pioneer.jl implements MS1 features through a specialized type hierarchy:

**Ms1UnscoredPSM** (`src/Routines/SearchDIA/PSMs/UnscoredPSMs.jl:57-73`)
```julia
struct Ms1UnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    m0::Bool                           # M0 isotope peak detected
    n_iso::UInt8                       # Number of isotope peaks observed
    big_iso::UInt8                     # Highest isotope index observed
    m0_error::Union{Missing,T}         # Mass error of M0 peak (ppm)
    error::T                           # Cumulative mass error
    precursor_idx::UInt32              # Precursor identifier
end
```

**Ms1ScoredPSM** (`src/Routines/SearchDIA/PSMs/ScoredPSMs.jl:93-118`)
```julia
struct Ms1ScoredPSM{H,L<:AbstractFloat} <: ScoredPSM{H,L}
    # Isotope Pattern Features
    m0::Bool                           # M0 peak presence
    n_iso::UInt8                       # Isotope peak count
    big_iso::UInt8                     # Maximum isotope index
    m0_error::Union{Missing,L}         # M0 mass error
    error::L                           # Log-transformed cumulative error

    # Spectral Similarity Features
    spectral_contrast::L               # Raw spectral contrast
    fitted_spectral_contrast::L        # Fitted spectral contrast
    gof::L                            # Goodness of fit
    max_matched_residual::L           # Maximum matched residual
    max_unmatched_residual::L         # Maximum unmatched residual
    fitted_manhattan_distance::L      # Fitted Manhattan distance
    matched_ratio::L                  # Matched/unmatched ratio
    weight::H                         # Precursor weight

    # Identifiers
    precursor_idx::UInt32
    ms_file_idx::UInt32
    scan_idx::UInt32
end
```

**SpectralScoresMs1** (`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:42-51`)
```julia
struct SpectralScoresMs1{T<:AbstractFloat} <: SpectralScores{T}
    spectral_contrast::T
    fitted_spectral_contrast::T
    gof::T                            # Goodness of fit (-log)
    max_matched_residual::T           # -log normalized residual
    max_unmatched_residual::T         # -log normalized residual
    fitted_manhattan_distance::T      # -log normalized distance
    matched_ratio::T                  # Log ratio matched/unmatched
end
```

## MS1 Feature Calculation Pipeline

### 1. Precursor Matching (`src/structs/MatchIon.jl:100-122`)

MS1 features begin with `PrecursorMatch` objects that capture matches between theoretical precursor masses and observed MS1 peaks:

```julia
struct PrecursorMatch{T<:AbstractFloat} <: MatchIon{T}
    predicted_intensity::T            # Theoretical intensity
    intensity::T                      # Observed intensity
    theoretical_mz::T                 # Expected m/z
    observed_mz::T                    # Measured m/z
    iso_idx::UInt8                    # Isotope index (1=M0, 2=M+1, etc.)
    peak_ind::Int64                   # Peak index in spectrum
    prec_id::UInt32                   # Precursor identifier
end
```

### 2. Feature Extraction (`src/Routines/SearchDIA/PSMs/UnscoredPSMs.jl:219-248`)

The `ModifyFeatures!` function for `Ms1UnscoredPSM` accumulates isotope pattern information:

```julia
function ModifyFeatures!(score::Ms1UnscoredPSM{T}, prec_id::UInt32,
                        match::PrecursorMatch{T}, errdist::MassErrorModel,
                        m_rank::Int64) where {T<:Real}

    n_iso += one(UInt8)  # Count each isotope peak

    if match.iso_idx > big_iso
        big_iso = match.iso_idx  # Track highest isotope
    end

    # Calculate PPM error for this isotope peak
    ppm_err = (getMZ(match) - match.observed_mz)/(getMZ(match)/1e6)
    error += abs(ppm_err)

    # Special handling for M0 peak
    if getIsoIdx(match) == one(UInt8)  # M0 isotope
        m0 = true
        m0_error = abs(ppm_err)
    end

    return updated_Ms1UnscoredPSM
end
```

**Key Features Calculated:**
- **m0**: Boolean indicating M0 peak detection
- **n_iso**: Count of detected isotope peaks (≥2 required for scoring)
- **big_iso**: Highest isotope index (indicates envelope completeness)
- **m0_error**: Specific mass error for the monoisotopic peak
- **error**: Cumulative absolute mass error across all isotopes

### 3. Spectral Scoring (`src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:433-536`)

The `getDistanceMetrics` function for MS1 calculates sophisticated similarity metrics:

```julia
function getDistanceMetrics(w::Vector{T}, r::Vector{T}, H::SparseArray{Ti,T},
                           spectral_scores::Vector{SpectralScoresMs1{U}})
```

**Algorithm Overview:**
1. **Residual Calculation**: Computes fit residuals between theoretical and observed isotope patterns
2. **Spectral Contrast**: Measures correlation between theoretical and observed patterns
3. **Fitted Spectral Contrast**: Uses deconvoluted intensities for interference correction
4. **Goodness of Fit**: Evaluates overall pattern matching quality
5. **Residual Analysis**: Identifies and quantifies peak mismatches

**Key Calculations:**
```julia
# Spectral contrast (raw correlation)
spectral_contrast = dot_product/(sqrt(h2_sum)*sqrt(x2_sum))

# Fitted spectral contrast (interference-corrected)
fitted_spectral_contrast = fitted_dotp/(sqrt(fitted_dotp_norm1)*sqrt(sum_of_fitted_peaks_squared))

# Goodness of fit (log-normalized residual sum)
gof = -log2(sum_of_residuals/sum_of_fitted_peaks + 1e-10)

# Matched vs unmatched ratio
matched_ratio = log2(matched_sum/unmatched_sum)
```

### 4. Final Scoring (`src/Routines/SearchDIA/PSMs/ScoredPSMs.jl:329-387`)

The `Score!` function converts unscored PSMs to scored PSMs with quality filters:

```julia
function Score!(scored_psms::Vector{Ms1ScoredPSM{H, L}},
                unscored_PSMs::Vector{Ms1UnscoredPSM{H}},
                spectral_scores::Vector{SpectralScoresMs1{L}},
                weight::Vector{H}, IDtoCOL::ArrayDict{UInt32, UInt16},
                cycle_idx::Int64, last_val::Int64, n_vals::Int64, scan_idx::Int64)
```

**Quality Filters Applied:**
- **M0 Detection Required**: `unscored_PSMs[i].m0 == true`
- **Minimum Isotopes**: `unscored_PSMs[i].n_iso >= 2`
- **Weight Threshold**: `weight[i] >= 1e-6`

**Feature Transformations:**
- **Log Error**: `Float16(log2(unscored_PSMs[i].error + 1e-6))`
- **Clamped Ratio**: `min(spectral_scores[scores_idx].matched_ratio, Float16(10.0))`

## MS1 Mass Tolerance Estimation

### Current Implementation

Pioneer.jl maintains separate mass error models for MS1 and MS2 data through the `SearchContext`:

```julia
# From src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:216,264,469-476
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    mass_error_model::Dict{Int64, MassErrorModel}      # MS2 fragment tolerance
    ms1_mass_error_model::Dict{Int64, MassErrorModel}  # MS1 precursor tolerance
end

function getMs1MassErrorModel(s::SearchContext, index::I) where {I<:Integer}
   if haskey(s.ms1_mass_error_model, index)
       return s.ms1_mass_error_model[index]
   else
       @user_warn "Mass error model not found for ms_file_idx $index. Returning default +/- 30ppm"
       return MassErrorModel(zero(Float32), (30.0f0, 30.0f0))
   end
end
```

### Mass Error Model Structure

The `MassErrorModel` encapsulates both bias and tolerance:

```julia
struct MassErrorModel
    offset::Float32        # Systematic bias in ppm
    tolerances::Tuple{Float32, Float32}  # (left_tol, right_tol) in ppm
end
```

### Isotope Error Bounds

MS1 tolerance is further refined through isotope error bounds that account for mass defects:

```julia
# From selectTransitions implementations
isotope_err_bounds::Tuple{I, I} = (3, 1)  # (negative, positive) isotope offsets

# Applied in precursor selection:
mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge
```

Where `NEUTRON = 1.003355` Da represents the neutron mass difference.

### Parameter Tuning Integration

MS1 mass tolerance estimation appears to be integrated with the existing `ParameterTuningSearch` framework, which optimizes mass tolerances based on PSM count and quality metrics. The system likely:

1. **Collects Precursor Matches**: Uses initial tolerance estimates to gather MS1 precursor matches
2. **Analyzes Mass Errors**: Calculates PPM errors for precursor isotope patterns
3. **Fits Error Distribution**: Estimates systematic bias and tolerance bounds
4. **Validates Model**: Tests tolerance expansion for improved PSM recovery

## Feature Integration in Machine Learning

### Feature Sets

MS1 features are incorporated into the machine learning pipeline through feature sets defined in `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl`:

```julia
# MS1-related features likely included in:
const ADVANCED_FEATURE_SET = [
    # ... other features ...
    :m0_error,
    :spectral_contrast,
    :fitted_spectral_contrast,
    :gof,
    :matched_ratio,
    # ... isotope pattern features ...
]
```

### Scoring Pipeline Integration

MS1 features are processed alongside MS2 features in the SecondPassSearch:

```julia
# From src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl
getDistanceMetrics(weights, residuals, Hs, getMs1SpectralScores(search_data))

ScoreFragmentMatches!(
    getMs1UnscoredPsms(search_data),
    getIdToCol(search_data),
    precursor_matches,
    n_matches,
    mass_error_model,
    1  # rank parameter
)
```

## Performance Characteristics

### Computational Complexity

- **Isotope Pattern Matching**: O(n_peaks × n_precursors × n_isotopes)
- **Spectral Scoring**: O(n_precursors × n_isotopes²) for correlation calculations
- **Feature Extraction**: O(n_matches) linear in number of precursor matches

### Memory Usage

MS1 feature structures are designed for memory efficiency:
- **Float16** for low-precision spectral scores
- **UInt8** for counts and indices
- **Boolean** for binary features
- Sparse representation for isotope patterns

### Threading Considerations

MS1 features are computed in thread-local `SearchDataStructures`:

```julia
# From src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:183-185
ms1_scored_psms::Vector{Ms1ScoredPSM{Float32, Float16}}
ms1_unscored_psms::Vector{Ms1UnscoredPSM{Float32}}
ms1_spectral_scores::Vector{SpectralScoresMs1{Float16}}
```

This enables parallel processing across different MS files and scan ranges.

## Quality Control and Validation

### Feature Quality Metrics

1. **Isotope Pattern Completeness**: Measured by `n_iso` and `big_iso`
2. **Mass Accuracy**: Quantified through `m0_error` and cumulative `error`
3. **Spectral Similarity**: Assessed via multiple correlation metrics
4. **Signal Quality**: Evaluated through matched/unmatched ratios

### Filtering Criteria

The pipeline applies stringent quality filters:
- **M0 Peak Required**: Ensures monoisotopic peak detection
- **Minimum Isotope Count**: Requires ≥2 isotopes for reliable pattern analysis
- **Weight Threshold**: Filters low-confidence precursors
- **Error Bounds**: Excludes precursors with excessive mass errors

## Future Enhancement Opportunities

### 1. Adaptive Tolerance Estimation

The existing `feature/ms1-mass-tolerance-estimation` branch (referenced in git history) suggests ongoing work on:
- Dynamic tolerance adjustment based on instrument characteristics
- File-specific mass error model calibration
- Integration with existing ParameterTuningSearch framework

### 2. Enhanced Isotope Pattern Analysis

Potential improvements include:
- **Theoretical Isotope Modeling**: More sophisticated isotope distribution predictions
- **Pattern Shape Analysis**: Beyond just peak presence/absence
- **Charge State Validation**: Using isotope spacing for charge determination

### 3. Machine Learning Integration

Advanced ML approaches could leverage:
- **Deep Learning**: Neural networks for isotope pattern recognition
- **Feature Engineering**: Automated discovery of informative MS1 features
- **Transfer Learning**: Models trained on high-quality data applied to challenging datasets

## Conclusion

Pioneer.jl implements a comprehensive MS1 feature scoring system that significantly enhances peptide identification confidence. The system combines:

1. **Robust Isotope Pattern Analysis**: Through sophisticated pattern matching algorithms
2. **Spectral Similarity Scoring**: Using multiple correlation and goodness-of-fit metrics
3. **Mass Accuracy Assessment**: Via detailed error analysis and tolerance estimation
4. **Quality-Based Filtering**: Ensuring only high-confidence features contribute to scoring

The modular design allows for independent optimization of MS1 and MS2 features while maintaining computational efficiency through careful data structure design and parallel processing capabilities. The ongoing development of adaptive mass tolerance estimation promises further improvements in sensitivity and specificity across diverse instrument platforms and sample types.

## Code References

Key files for MS1 feature implementation:

- `src/Routines/SearchDIA/PSMs/UnscoredPSMs.jl:57-248` - Core MS1 PSM structures
- `src/Routines/SearchDIA/PSMs/ScoredPSMs.jl:93-387` - MS1 scoring pipeline
- `src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:42-536` - Spectral scoring algorithms
- `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:216,469-479` - Mass error model management
- `src/structs/MatchIon.jl:100-122` - Precursor match data structures
- `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl` - Integration with search pipeline