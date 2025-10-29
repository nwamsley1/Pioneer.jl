# CLAUDE.md - ParameterTuningSearch Implementation Guide

This file provides comprehensive documentation of the ParameterTuningSearch module in Pioneer.jl, detailing the current robust 3-phase convergence strategy with scan scaling.

## Purpose

ParameterTuningSearch is the first search method in the Pioneer.jl pipeline, responsible for establishing optimal fragment mass tolerance and retention time alignment models. It serves as the foundation for all subsequent search methods by calibrating critical parameters based on the actual data characteristics.

## Core Algorithm Flow

### High-Level Control Flow Schematic

```
╔═══════════════════════════════════════════════════════════════════╗
║                     PARAMETERTUNINGSEARCH                        ║
║                  Mass Error & RT Calibration                     ║
╚═══════════════════════════════════════════════════════════════════╝
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Initialize Models   │
                    │ - Mass Error: ±5ppm │
                    │ - Quad: Gaussian    │
                    └─────────────────────┘
                               │
                               ▼
        ┌──────────────────────────────────────────────┐
        │     SCAN SCALING LOOP (Outer Loop)          │
        │  Initial: 500 scans → Scale by 10x → Max: 80k│
        └──────────────────────────────────────────────┘
                               │
                  ┌────────────┴────────────┐
                  ▼                          ▼
        ┌─────────────────┐       ┌──────────────────┐
        │ Create/Expand   │       │ If Failed:       │
        │ FilteredSpectra │       │ Scale scans 10x  │
        └─────────────────┘       └──────────────────┘
                  │
                  ▼
    ┌──────────────────────────────────────────┐
    │        3-PHASE LOOP (Middle Loop)        │
    │   Phase 1: Bias = 0 ppm                  │
    │   Phase 2: Bias = +max_tol ppm           │
    │   Phase 3: Bias = -max_tol ppm           │
    └──────────────────────────────────────────┘
                  │
         ┌────────┴────────┬────────────┐
         ▼                 ▼            ▼
    [Phase 1]         [Phase 2]    [Phase 3]
    Zero Bias      Positive Bias  Negative Bias
         │                 │            │
         └────────┬────────┴────────────┘
                  ▼
    ┌──────────────────────────────────────────┐
    │    SCORE THRESHOLD LOOP (per Phase)      │
    │    Tries each min_score in sequence      │
    │    e.g., [22, 17, 15]                    │
    └──────────────────────────────────────────┘
                  │
         ┌────────┼────────┬────────────┐
         ▼        ▼        ▼            ▼
    [Score 22] [Score 17] [Score 15] [Score N]
         │        │        │            │
         └────────┴────────┴────────────┘
                  ▼
    ┌──────────────────────────────────────────┐
    │   ITERATION LOOP (Inner Loop per Score)  │
    │      Up to 3 iterations per score        │
    └──────────────────────────────────────────┘
                  │
                  ▼
         [See Mass Error Estimation 
          Process Flow Below]
```

### Mass Error Estimation Process Flow

```
╔═══════════════════════════════════════════════════════════════════╗
║              MASS ERROR ESTIMATION WITHIN EACH ITERATION          ║
╚═══════════════════════════════════════════════════════════════════╝

Step 0: Initial PSM Collection (Before Iteration Loop)
┌──────────────────────────────────────────────────────────────────┐
│  1. Fragment Index Search with current tolerance + phase bias    │
│  2. Score PSMs (spectral contrast, matched ratio)                │
│  3. Filter by FDR (q-value < 0.01)                              │
│  4. Attempt convergence check                                    │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
                    [If not converged,
                     enter iteration loop]
                            │
                            ▼
Step 1: Expand Mass Tolerance
┌──────────────────────────────────────────────────────────────────┐
│  tolerance *= mass_tolerance_scale_factor (default: 2.0)         │
│  Apply to both left and right tolerances                         │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 2: Collect PSMs for Bias Detection
┌──────────────────────────────────────────────────────────────────┐
│  Library Search Pipeline:                                        │
│  1. Fragment index search with expanded tolerance                │
│  2. Match fragments to spectra                                   │
│  3. Calculate match scores                                       │
│  4. Return PSM DataFrame                                         │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 3: Fit Mass Error Model for Bias
┌──────────────────────────────────────────────────────────────────┐
│  Fragment-Level Mass Error Calculation:                          │
│  1. Extract matched fragments from PSMs                          │
│  2. Calculate PPM error for each fragment:                       │
│     ppm_err = (observed_mz - theoretical_mz) / (theoretical/1e6) │
│  3. Filter by intensity (remove low-quality matches)             │
│  4. Calculate median PPM error → systematic bias                 │
│  5. Subtract bias, calculate MAD for tolerance bounds            │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 4: Adjust Bias in Mass Error Model
┌──────────────────────────────────────────────────────────────────┐
│  Update MassErrorModel with detected bias:                       │
│  - Keep expanded tolerance from Step 1                           │
│  - Replace bias offset with newly calculated value               │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 5: Collect PSMs with Adjusted Bias
┌──────────────────────────────────────────────────────────────────┐
│  Repeat library search with bias-corrected model                 │
│  Should yield more/better PSMs if bias was significant           │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 6: Fit Final Models
┌──────────────────────────────────────────────────────────────────┐
│  1. Mass Error Model:                                            │
│     - Refit using all fragments from Step 5 PSMs                 │
│     - Calculate final bias and tolerance bounds                  │
│  2. RT Alignment Model:                                          │
│     - Fit spline: empirical_RT → library_iRT                     │
│     - Remove outliers iteratively                                │
└──────────────────────────────────────────────────────────────────┘
                            │
                            ▼
Step 7: Convergence Check & Optional Expansion
┌──────────────────────────────────────────────────────────────────┐
│  Convergence Criteria:                                           │
│  - PSM count ≥ min_psms (default: 100)                          │
│  - |mass_offset| < init_tolerance/4                              │
│  - Fitted tolerance reasonable                                   │
│                                                                  │
│  If Converged:                                                   │
│  - Test 1.5x tolerance expansion                                 │
│  - If >10% more PSMs, refit with expanded set                    │
│  - Store final models in SearchContext                           │
└──────────────────────────────────────────────────────────────────┘
```

## Current Implementation: 3-Phase Convergence with Multi-Score Search

### Overview

The module implements a sophisticated convergence strategy that combines:
1. **3-phase search** with different bias shifts to explore the parameter space
2. **Multi-score thresholds** within each phase to handle varying data quality
3. **Scan count scaling** between attempts to gather more data when needed
4. **Mass tolerance expansion** within iterations to progressively widen the search
5. **Post-convergence optimization** to capture additional PSMs

### Key Configuration Parameters

```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "min_score": [22, 17, 15],  // Can be single value or array
      "min_count": 7,
      "min_spectral_contrast": 0.5
    },
    "iteration_settings": {
      "init_mass_tol_ppm": [20.0, 30.0],  // Explicit tolerances per iteration
      "ms1_tol_ppm": 20.0,
      "scan_counts": [500, 5000, 50000, 80000]  // Explicit scan counts to try
    }
  }
}
```

### Multi-Score Search Feature (2025-01)

The `min_score` parameter now accepts an array of thresholds to try sequentially:
- **Single value**: `"min_score": 22` - Backwards compatible
- **Array**: `"min_score": [22, 17, 15]` - Tries each threshold in order

Each phase independently explores all score thresholds with its specific bias setting:
- Tries stricter thresholds first (higher scores = better quality)
- Relaxes to lower thresholds only if needed
- Returns immediately upon convergence with any score
- Logs clearly which phase and score achieved convergence

### Phase Strategy

**Phase 1: Zero Bias**
- Starts with no bias shift
- Explores the nominal mass range
- Most likely to succeed for well-calibrated instruments

**Phase 2: Positive Bias Shift**
- Applies positive bias shift equal to maximum theoretical tolerance
- Explores systematic positive mass errors
- Handles instruments with positive calibration drift

**Phase 3: Negative Bias Shift**
- Applies negative bias shift equal to maximum theoretical tolerance
- Explores systematic negative mass errors
- Handles instruments with negative calibration drift

### Post-Convergence Tolerance Expansion

After achieving convergence, the algorithm tests whether expanding the collection tolerance by 50% yields significantly more PSMs:

```julia
function test_tolerance_expansion!(...)
    expanded_tolerance = collection_tolerance * 1.5
    expanded_psms = collect_with_expanded_tolerance()
    
    if psm_count_increased
        refit_model_with_expanded_psms()
        return expanded_results
    else
        return original_results
    end
end
```

## Core Components

### PSM Collection Pipeline

1. **Library Search** (`library_search`)
   - Fragment index search for candidate precursors
   - Detailed PSM generation with spectral matching
   - Returns DataFrame of PSMs with scores

2. **Column Addition** (`add_tuning_search_columns!`)
   - Adds RT, iRT, charge, TIC columns
   - Calculates matched ratios
   - Prepares for scoring

3. **Scoring and Filtering** (`filter_and_score_psms!`)
   - Probit regression scoring
   - FDR control with q-values
   - Best PSM per precursor selection

### Model Fitting

**RT Model** (`fit_irt_model`):
- Spline fitting with outlier removal
- Returns SplineRtConversionModel
- Stores RT/iRT pairs for QC plots

**Mass Error Model** (`fit_mass_err_model`):
- Fragment-based mass error estimation
- Intensity filtering for robust fitting
- Returns MassErrorModel with offset and tolerances

### Convergence Criteria

Convergence is achieved when:
1. Sufficient PSMs collected (≥ `min_psms`)
2. Mass offset is reasonable (< init_tolerance/4)
3. Fitted tolerance is not excessive
4. Mass error distribution is stable

## Fallback Strategies

When convergence fails after all attempts:

1. **Parameter Borrowing** (`get_fallback_parameters`):
   - Attempts to borrow from successfully tuned neighboring files
   - Prefers files processed earlier in the run
   - Falls back to conservative defaults if no files available

2. **Conservative Defaults**:
   - Mass tolerance: ±50 ppm
   - RT model: Identity transformation
   - Logs warning for user awareness

## Output and QC

### Generated Files
- RT alignment plots per file
- Mass error distribution plots
- Combined PDFs for all files
- Diagnostic summary in logs

### SearchContext Updates
- Stores calibrated MassErrorModel
- Stores RT-to-iRT conversion model
- Updates IRT error estimates
- Makes parameters available to downstream methods

## Key Implementation Details

### Thread Safety
- Each thread has independent SearchDataStructures
- Results collected in thread-local storage
- Final merging in single-threaded context

### Memory Management
- FilteredMassSpecData created once and reused
- Scans appended incrementally
- Efficient sparse matrix operations

### Performance Optimizations
- Binary search in fragment index
- Exponential bound expansion
- TopN peak filtering for large tolerances
- Parallel processing where possible

## Common Issues and Solutions

### Issue: Phase 2/3 Failing with 0 PSMs
**Cause**: Bias shift calculation was using uninitialized max tolerance
**Solution**: Calculate bias shift from iteration settings correctly

### Issue: Scan Scaling Within Phases
**Cause**: Scan count was doubling within phase iterations
**Solution**: Separated scan scaling (between attempts) from mass tolerance scaling (within phases)

### Issue: Parameter Extraction Errors
**Cause**: Looking for parameters in wrong JSON hierarchy
**Solution**: Properly extract from `iteration_settings` under `tuning_params`

### Issue: Convergence with Conservative Tolerances
**Cause**: Not exploring tolerance boundaries adequately
**Solution**: Added post-convergence tolerance expansion test

## Integration with Downstream Methods

The calibrated parameters are critical for:
- **FirstPassSearch**: Uses mass error model and RT conversion
- **SecondPassSearch**: Refined search with calibrated parameters
- **All quantification methods**: Rely on RT alignment for XIC extraction

## Recent Changes (2025-01)

1. **Multi-Score Search Enhancement** (Latest):
   - `min_score` parameter now accepts arrays: `[22, 17, 15]`
   - Each phase independently tries all score thresholds
   - Maintains backwards compatibility with single values
   - Other search methods automatically use first value from array
   - Clear logging shows which score achieved convergence

2. **Removed deprecated parameters**:
   - `tol_ppm` → `init_mass_tol_ppm` in iteration_settings
   - `max_tolerance_ppm` → calculated dynamically

3. **Implemented robust convergence**:
   - 3-phase approach with bias shifts
   - Scan scaling with configurable factor
   - Post-convergence tolerance expansion
   - Score threshold exploration within phases

4. **Code cleanup**:
   - Removed unused functions from old implementation
   - Commented out verbose @info logging
   - Kept only essential @warn statements

5. **Configuration improvements**:
   - Centralized iteration settings
   - Made scan scaling configurable
   - Better parameter organization
   - Flexible min_score configuration

## Performance Characteristics

- **Memory Usage**: Moderate - stores PSMs and filtered spectra
- **CPU Usage**: High during fragment matching and model fitting
- **I/O**: Minimal - works with in-memory data
- **Typical Runtime**: 30-120 seconds per file depending on scan count
- **Scaling**: Linear with scan count, quadratic with tolerance

## Future Improvements

1. **Adaptive strategies**: Dynamic adjustment of scale factors
2. **Machine learning**: Predict optimal parameters from data characteristics
3. **Parallel phase execution**: Run phases concurrently when possible
4. **Better diagnostics**: Track convergence metrics for analysis