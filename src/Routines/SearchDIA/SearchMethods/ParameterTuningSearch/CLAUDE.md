# CLAUDE.md - ParameterTuningSearch Implementation Guide

This file provides comprehensive documentation of the ParameterTuningSearch module in Pioneer.jl, detailing the current robust 3-phase convergence strategy with scan scaling.

## Purpose

ParameterTuningSearch is the first search method in the Pioneer.jl pipeline, responsible for establishing optimal fragment mass tolerance and retention time alignment models. It serves as the foundation for all subsequent search methods by calibrating critical parameters based on the actual data characteristics.

## Current Implementation: 3-Phase Convergence with Scan Scaling

### Overview

The module implements a sophisticated convergence strategy that combines:
1. **3-phase search** with different bias shifts to explore the parameter space
2. **Scan count scaling** between attempts to gather more data when needed
3. **Mass tolerance expansion** within phases to progressively widen the search
4. **Post-convergence optimization** to capture additional PSMs

### Key Configuration Parameters

```json
{
  "parameter_tuning": {
    "iteration_settings": {
      "init_mass_tol_ppm": 5.0,
      "mass_tolerance_scale_factor": 2.0,
      "iterations_per_phase": 3,
      "scan_scale_factor": 10.0
    },
    "search_settings": {
      "initial_scan_count": 500,
      "max_parameter_tuning_scans": 80000
    }
  }
}
```

### Main Convergence Algorithm (`process_file!`)

The algorithm operates with nested loops:
1. **Outer loop**: Scan scaling attempts
2. **Middle loop**: 3 phases with different bias shifts
3. **Inner loop**: Iterations within each phase

```julia
# Pseudocode structure
for scan_attempt in 1:max_attempts
    filtered_spectra = create_or_expand_spectra(scan_count)
    
    for phase in 1:3
        bias_shift = get_phase_bias(phase)  # 0, +max_tol, -max_tol
        
        # Initial attempt at base tolerance
        if check_convergence(initial_psms)
            return success
        end
        
        # Iteration loop with tolerance expansion
        for iteration in 1:iterations_per_phase
            expand_tolerance()
            collect_psms_for_bias()
            adjust_bias()
            collect_psms_with_adjusted_bias()
            if check_convergence()
                return success
            end
        end
    end
    
    # Scale up scan count for next attempt
    scan_count *= scan_scale_factor
end
```

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

### Within-Phase Iteration Strategy

Each phase executes the following steps:

1. **Initial Attempt** (before iteration loop):
   - Collect PSMs at base tolerance with phase bias
   - Fit models and check convergence
   - If converged, test tolerance expansion

2. **Iteration Loop** (up to `iterations_per_phase`):
   - **Step 1**: Expand tolerance by `mass_tolerance_scale_factor`
   - **Step 2**: Collect PSMs with expanded tolerance
   - **Step 3**: Fit model to determine optimal bias
   - **Step 4**: Adjust bias based on fitted model
   - **Step 5**: Collect PSMs with adjusted bias
   - **Step 6**: Fit final models
   - **Step 7**: Check convergence and optionally expand

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

1. **Removed deprecated parameters**:
   - `tol_ppm` → `init_mass_tol_ppm` in iteration_settings
   - `max_tolerance_ppm` → calculated dynamically

2. **Implemented robust convergence**:
   - 3-phase approach with bias shifts
   - Scan scaling with configurable factor
   - Post-convergence tolerance expansion

3. **Code cleanup**:
   - Removed unused functions from old implementation
   - Commented out verbose @info logging
   - Kept only essential @warn statements

4. **Configuration improvements**:
   - Centralized iteration settings
   - Made scan scaling configurable
   - Better parameter organization

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