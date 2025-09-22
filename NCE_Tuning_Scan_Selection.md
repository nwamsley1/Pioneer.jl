# NCE Tuning Search: Scan Selection and PSM Collection

## Overview

The NCE (Normalized Collision Energy) tuning search in Pioneer.jl employs an intelligent scan selection strategy combined with a robust PSM (Peptide-Spectrum Match) collection pipeline. This document details the mechanisms used to select representative scans and collect PSMs for parameter optimization.

## Scan Selection Strategy

### 1. FilteredMassSpecData Creation

The system creates a filtered view of mass spectrometry data with configurable parameters:

- **Initial scan count**: 500 scans (configurable via `initial_scan_count`)
- **Maximum scans**: 80,000 (configurable via `max_parameter_tuning_scans`)
- **Scan scaling factor**: 10x between attempts (configurable via `scan_scale_factor`)
- **Target MS order**: MS2 scans only (UInt8(2))

### 2. Intelligent Scan Prioritization

The scan selection uses a three-step RT-binning strategy to ensure representative sampling:

#### Step 1: RT Binning
- Divides the retention time range into 15 bins (default)
- Assigns each scan to its corresponding RT bin based on retention time
- Ensures coverage across the entire chromatographic run
- Handles edge cases (empty data, single RT value)

#### Step 2: Peak Density Sorting
- Within each RT bin, sorts scans by TIC (Total Ion Current)
- TIC serves as a proxy for peak density and information content
- Higher TIC scans are prioritized as they typically contain more informative spectra
- Sorting is performed in descending order (highest TIC first)

#### Step 3: Interleaved Sampling
- Takes scans in round-robin fashion from all bins
- First scan from bin 1, then first from bin 2, etc.
- Ensures temporal distribution of selected scans
- Prevents bias toward any particular retention time region

### 3. Dynamic Scan Scaling

When convergence fails, the system automatically scales up:

- Multiplies scan count by 10x (e.g., 500 → 5,000 → 50,000)
- Uses `append!` method to add more scans from the pre-computed priority order
- Maintains the original prioritization strategy
- Stops at `max_parameter_tuning_scans` limit

## PSM Collection Pipeline

### 1. Library Search Pipeline

The PSM collection follows a multi-step process:

```
Fragment Index Search → PSM Generation → Column Addition → Scoring → Filtering
```

### 2. Three-Phase Strategy with Bias Exploration

The search explores different mass bias scenarios to handle instrument calibration drift:

- **Phase 1: Zero Bias**
  - Searches with no bias shift
  - Explores the nominal mass range
  - Most likely to succeed for well-calibrated instruments

- **Phase 2: Positive Bias**
  - Applies positive bias shift equal to maximum theoretical tolerance
  - Explores systematic positive mass errors
  - Handles instruments with positive calibration drift

- **Phase 3: Negative Bias**
  - Applies negative bias shift equal to maximum theoretical tolerance
  - Explores systematic negative mass errors
  - Handles instruments with negative calibration drift

### 3. Multi-Score Threshold Search

Within each phase, the system tries multiple score thresholds sequentially:

- **Default thresholds**: [22, 17, 15] (configurable array in `min_score`)
- Starts with the strictest threshold (highest quality PSMs)
- Progressively relaxes to lower thresholds if needed
- Returns immediately upon achieving convergence with any score
- Clear logging indicates which phase and score achieved convergence

### 4. Iterative Mass Error Refinement

For each score threshold, the system performs iterative refinement:

1. **Initial PSM Collection**
   - Collects PSMs at base tolerance (5 ppm default)
   - Applies phase-specific bias shift
   - Filters by current score threshold

2. **Tolerance Expansion** (per iteration)
   - Multiplies tolerance by 2x (configurable via `mass_tolerance_scale_factor`)
   - Allows capture of PSMs with larger mass errors

3. **Bias Detection**
   - Collects PSMs with expanded tolerance
   - Calculates median PPM error from fragment matches
   - Identifies systematic mass offset

4. **Bias Adjustment**
   - Updates mass error model with detected bias
   - Keeps expanded tolerance from previous step
   - Prepares for refined PSM collection

5. **Final Collection**
   - Collects PSMs with bias-adjusted parameters
   - Should yield more/better PSMs if bias was significant

6. **Convergence Check**
   - Minimum 100 PSMs required (configurable via `min_psms`)
   - Mass offset must be < initial_tolerance/4
   - Fitted tolerance must be reasonable

## Key Collection Parameters

### Fragment Matching

- **Selection strategy**: Top N fragments per precursor for mass tolerance calibration
- **Max fragments**: Configurable via `max_frags_for_mass_err_estimation`
- **Intensity filtering**: Removes low-intensity matches for robust error estimation
- **Error calculation**: PPM error = (observed - theoretical) / (theoretical/1e6)

### PSM Scoring and Filtering

- **Scoring features**: 8 features used in probit regression
  - entropy_score
  - city_block
  - scribe
  - spectral_contrast
  - y_count
  - error
  - TIC
  - intercept

- **FDR control**: q-value threshold ≤ 0.01
- **Target filtering**: Only target PSMs used for parameter fitting (decoys removed)
- **Best PSM selection**: Currently commented out, all passing PSMs retained

### Convergence Criteria

- **Minimum PSM count**: ≥ 100 PSMs (configurable)
- **Mass offset stability**: |mass_offset| < init_tolerance/4
- **Reasonable tolerance**: Fitted tolerance within expected bounds
- **Multiple attempts**: Up to 3 iterations per score, per phase

## Post-Convergence Optimization

### Tolerance Expansion Test

After achieving convergence, the system tests tolerance expansion:

1. Expands collection tolerance by 50% (1.5x)
2. Collects PSMs with expanded tolerance
3. If PSM count increases by >10%, refits models with expanded set
4. Maximizes PSM yield while maintaining accuracy

### Best Attempt Tracking

The system maintains fallback options:

- Tracks best parameters across all attempts
- Records phase, score, iteration, and PSM count
- Falls back to best attempt if convergence fails
- Ensures usable parameters for downstream analysis

## Implementation Details

### Memory Efficiency

- **FilteredMassSpecData**: Stores only sampled scans in memory
- **TopN filtering**: Reduces memory footprint by keeping only most intense peaks
- **Thread-local storage**: Minimizes memory contention in parallel processing
- **Incremental append**: Adds scans as needed rather than pre-allocating maximum

### Performance Optimizations

- **Parallel processing**: PSM collection distributed across threads
- **Binary search**: Efficient fragment index searching
- **Sparse operations**: Optimized matrix operations for large datasets
- **Pre-computed priority**: Scan order calculated once, reused for scaling

### Thread Safety

- Each thread has independent SearchDataStructures
- Results collected in thread-local storage
- Final merging performed in single-threaded context
- No shared state mutations during parallel execution

## Configuration Parameters

### Key Parameters (from JSON configuration)

```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "min_score": [22, 17, 15],  // Score thresholds to try
      "min_count": 7,              // Minimum fragment count
      "min_spectral_contrast": 0.5 // Spectral quality threshold
    },
    "iteration_settings": {
      "init_mass_tol_ppm": 5.0,           // Initial mass tolerance
      "mass_tolerance_scale_factor": 2.0,  // Expansion per iteration
      "iterations_per_phase": 3,           // Max iterations per phase
      "scan_scale_factor": 10.0           // Scan count scaling
    },
    "search_settings": {
      "initial_scan_count": 500,           // Starting scan count
      "max_parameter_tuning_scans": 80000  // Maximum scan limit
    }
  }
}
```

## Failure Handling

### Fallback Strategies

When convergence fails after all attempts:

1. **Best Attempt Recovery**
   - Uses parameters from iteration with most PSMs
   - Applies best mass error model and RT alignment
   - Tests 1.5x tolerance expansion on best parameters

2. **Parameter Borrowing**
   - Attempts to borrow from successfully tuned neighboring files
   - Prefers files processed earlier in the run
   - Falls back to conservative defaults if unavailable

3. **Conservative Defaults**
   - Mass tolerance: ±50 ppm
   - RT model: Identity transformation
   - IRT error: typemax(Float32) (effectively no RT filtering)
   - Logs warnings for user awareness

## Quality Control

### Generated Outputs

- RT alignment plots (per file and combined PDF)
- Mass error distribution plots (per file and combined PDF)
- Diagnostic summary in logs
- Parameter history tracking

### SearchContext Updates

Upon successful completion:
- Stores calibrated MassErrorModel
- Stores RT-to-IRT conversion model (SplineRtConversionModel)
- Updates IRT error estimates
- Makes parameters available to all downstream search methods

## Recent Enhancements (2025-01)

1. **Multi-Score Search**: `min_score` parameter now accepts arrays for flexible threshold exploration
2. **Robust Convergence**: 3-phase approach with systematic bias exploration
3. **Dynamic Scaling**: Configurable scan scaling with intelligent prioritization
4. **Post-Convergence Optimization**: Automatic tolerance expansion testing
5. **Improved Fallback**: Best attempt tracking with parameter borrowing

---

## Comparison with July 2024 Implementation

### Major Architectural Changes

The parameter tuning search has undergone significant architectural improvements between July 2024 and the current implementation:

#### July 2024 Implementation (`src/Routines/LibrarySearch/parameterTuningSearch.jl`)

**Structure:**
- Simple linear script without modular organization
- Direct iteration over MS files with basic while loop
- No systematic scan selection strategy
- Fixed sampling rate (`sample_rate` parameter)

**Key Characteristics:**
1. **Random Sampling**: Used `sample_rate` parameter for random scan sampling
2. **Single-Phase Search**: No systematic bias exploration
3. **Fixed Iterations**: Limited to `max_presearch_iters` attempts
4. **Basic Convergence**: Only checked PSM count at fixed FDR threshold
5. **No Fallback Strategy**: Warning message but no parameter recovery

#### Current Implementation (2025-01)

**Structure:**
- Modular, object-oriented design with SearchMethod interface
- Sophisticated state management with `IterationState`
- Intelligent scan selection with RT binning
- Dynamic scaling with configurable factors

**Major Improvements:**

### 1. Scan Selection Evolution

**July 2024:**
```julia
# Random sampling with fixed rate
sample_rate = params_[:presearch_params]["sample_rate"]
```

**Current:**
```julia
# Intelligent RT-binning with peak density prioritization
- RT binning (15 bins across chromatogram)
- TIC-based sorting within bins
- Interleaved round-robin sampling
- Dynamic scaling (500 → 5,000 → 50,000 scans)
```

### 2. Search Strategy Enhancement

**July 2024:**
```julia
while n <= params_[:presearch_params]["max_presearch_iters"]
    # Single attempt with fixed parameters
    mass_err_model = MassErrorModel(0.0f0, (frag_tol_ppm, frag_tol_ppm))
```

**Current:**
```julia
# Three-phase bias exploration
for phase in 1:3
    # Multiple score thresholds per phase
    for min_score in [22, 17, 15]
        # Iterative refinement with tolerance expansion
        for iter in 1:iterations_per_phase
            # Expand → Detect bias → Adjust → Refit cycle
```

### 3. Convergence Criteria

**July 2024:**
```julia
# Simple PSM count threshold
if sum(rtPSMs[!,:q_value].<=max_qval) >= min_samples
    break
end
```

**Current:**
```julia
# Multi-factor convergence check
- PSM count ≥ min_psms
- |mass_offset| < init_tolerance/4
- Reasonable fitted tolerance
- Post-convergence expansion test
```

### 4. Error Recovery

**July 2024:**
```julia
if n >= max_presearch_iters
    @warn "Presearch did not find $min_samples precursors"
    # Continue with whatever was found
```

**Current:**
```julia
# Sophisticated fallback strategy
- Best attempt tracking throughout search
- Parameter borrowing from successful files
- Conservative defaults as last resort
- Detailed logging of fallback reason
```

### 5. Performance Optimizations

**July 2024:**
- Sequential processing
- Full data loading
- No memory optimization

**Current:**
- Parallel PSM collection
- FilteredMassSpecData for memory efficiency
- Thread-local storage to minimize contention
- Pre-computed scan priorities

### 6. Quality Control

**July 2024:**
- Basic RT alignment plots
- Simple mass error plots
- Manual PDF merging

**Current:**
- Comprehensive diagnostic tracking
- Parameter history recording
- Iteration state visualization
- Automated multi-page PDF generation

### Key Improvements Summary

1. **Robustness**: From single-attempt to multi-phase exploration with fallback
2. **Efficiency**: From random sampling to intelligent prioritization
3. **Scalability**: From fixed sampling to dynamic scaling
4. **Reliability**: From basic convergence to sophisticated validation
5. **Diagnostics**: From minimal logging to comprehensive tracking
6. **Modularity**: From script to object-oriented SearchMethod interface

The current implementation represents a complete redesign focused on robustness, efficiency, and reliability, particularly for challenging datasets where the July 2024 version would have failed to converge.