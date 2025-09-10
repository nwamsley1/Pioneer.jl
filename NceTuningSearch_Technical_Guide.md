# NceTuningSearch Technical Guide

## Overview

The **NceTuningSearch** module is the second search method in Pioneer.jl's 9-stage DIA pipeline. Its primary purpose is to calibrate normalized collision energy (NCE) parameters by building instrument-specific models that predict optimal NCE values based on precursor mass-to-charge ratio (m/z) and charge state.

## Purpose and Context

### Why NCE Calibration Matters
- **Fragmentation Efficiency**: Different precursor ions require different collision energies for optimal fragmentation
- **Instrument Variability**: Each mass spectrometer has unique NCE characteristics that must be learned from data
- **Mass Dependence**: Heavier precursors typically require higher NCE values
- **Charge State Effects**: Higher charge states often need different NCE adjustments

### Pipeline Position
NceTuningSearch runs after **ParameterTuningSearch** (which establishes mass tolerances) and before **QuadTuningSearch** (which models quadrupole transmission). It uses the calibrated mass error models from ParameterTuningSearch to perform accurate spectral matching.

## Algorithm Architecture

### High-Level Flow

```
┌─────────────────────────────────────────────────────────────┐
│                    NceTuningSearch Flow                     │
└─────────────────────────────────────────────────────────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │ Check Library Type  │
                    │ (Skip if basic      │
                    │ FragmentIndexLib)   │
                    └─────────────────────┘
                               │
                               ▼
        ┌──────────────────────────────────────────────┐
        │        PSM COLLECTION LOOP                   │
        │     Up to 10 iterations to gather            │
        │     sufficient PSMs (min_samples)            │
        └──────────────────────────────────────────────┘
                               │
                  ┌────────────┴────────────┐
                  ▼                          ▼
        ┌─────────────────┐       ┌──────────────────┐
        │ Library Search  │       │ Process & Filter │
        │ with NCE grid   │       │ PSMs (FDR, etc.) │
        └─────────────────┘       └──────────────────┘
                  │                          │
                  └────────┬───────────────┘
                           ▼
    ┌──────────────────────────────────────────────────┐
    │          PIECEWISE NCE MODEL FITTING             │
    │                                                  │
    │  f(mz, charge) = {                               │
    │    left_slope * mz + left_intercept + charge_slope * charge  (mz ≤ 500) │
    │    right_value + charge_slope * charge           (mz > 500)   │
    │  }                                               │
    └──────────────────────────────────────────────────┘
                           │
                           ▼
    ┌──────────────────────────────────────────────────┐
    │              VISUALIZATION                       │
    │  • Scatter plot of actual PSM data points        │
    │  • Fitted model lines per charge state           │
    │  • PSM counts and model parameters               │
    └──────────────────────────────────────────────────┘
                           │
                           ▼
    ┌──────────────────────────────────────────────────┐
    │          STORE RESULTS                           │
    │  • Save NCE model to SearchContext               │
    │  • Generate combined PDF plots                   │
    │  • Log diagnostic information                    │
    └──────────────────────────────────────────────────┘
```

## Detailed Implementation

### 1. PSM Collection Strategy

The module uses an iterative approach to collect sufficient PSMs for reliable NCE modeling:

#### Grid Search Process
```julia
# NCE grid: 30 values from 21.0 to 40.0
nce_grid = LinRange{Float32}(21.0f0, 40.0f0, 30)

# Iterative collection (up to 10 attempts)
for i in 1:10
    psms = library_search(spectra, search_context, params, ms_file_idx)
    processed_psms = process_psms!(psms, spectra, search_context, params)
    
    if nrow(processed_psms) >= min_samples
        break  # Sufficient data collected
    end
end
```

#### PSM Processing Pipeline
1. **Library Search**: Fragment index search across NCE grid values
2. **Column Addition**: RT, iRT, charge, TIC information
3. **Scoring**: Probit regression scoring using presearch model
4. **FDR Control**: Q-value calculation and filtering (q ≤ 0.01)
5. **Best PSM Selection**: Highest SCRIBE score per precursor
6. **Target Filtering**: Remove decoy matches

### 2. Piecewise NCE Model

#### Mathematical Formulation

The NCE model uses a **piecewise linear function** with a breakpoint at 500 m/z:

```
f(mz, charge) = {
    left_slope × mz + left_intercept + charge_slope × charge,  if mz ≤ 500
    right_value + charge_slope × charge,                       if mz > 500
}
```

Where:
- **left_slope**: Slope of the linear region (low m/z)
- **left_intercept**: Y-intercept of the linear region
- **right_value**: Constant NCE value for high m/z (enforced continuity: `left_slope × 500 + left_intercept`)
- **charge_slope**: Linear charge dependence (applied uniformly)

#### Model Fitting Process

1. **Objective Function**: Minimize sum of squared residuals
   ```julia
   SSR = Σ(predicted_NCE - observed_NCE)²
   ```

2. **Optimization Method**: LBFGS algorithm from Optim.jl

3. **Continuity Constraint**: Right value automatically calculated to ensure smooth transition at breakpoint

4. **Initial Guess**: Linear regression on left-hand side data (mz ≤ 500)

### 3. Quality Control and Diagnostics

#### Diagnostic Logging
The enhanced version now provides comprehensive diagnostics:

```julia
@info "NceTuningSearch: Data summary:" *
      "\n  - Total PSMs: $(nrow(processed_psms))" *
      "\n  - Precursor m/z range: $(min_mz)-$(max_mz)" *
      "\n  - NCE range: $(min_nce)-$(max_nce)" *
      "\n  - Charge states: $(unique_charges)" *
      "\n  - PSMs per charge: $(psms_per_charge)"

@info "NceTuningSearch: Fitted NCE model parameters:" *
      "\n  - Breakpoint: $(nce_model.breakpoint) m/z" *
      "\n  - Left slope: $(nce_model.left_slope)" *
      "\n  - Left intercept: $(nce_model.left_intercept)" *
      "\n  - Right value: $(nce_model.right_value)" *
      "\n  - Charge slope: $(nce_model.charge_slope)"
```

#### Visualization Enhancements
- **Scatter Plot**: All PSM data points with transparency (shows data density)
- **Fitted Lines**: Model predictions for each charge state
- **Legend**: Charge state with PSM count (e.g., "+2 (n=1234)")
- **Annotations**: NCE values at maximum m/z for each charge
- **Axes**: Labeled precursor m/z vs NCE

## Parameter Configuration

### Key Parameters

```json
{
  "parameter_tuning": {
    "fragment_settings": {
      "min_score": 22,              // Minimum fragment match score
      "min_count": 7,               // Minimum fragment count
      "min_spectral_contrast": 0.5, // Minimum spectral contrast
      "min_log2_ratio": -3.0        // Minimum log2 matched ratio
    },
    "search_settings": {
      "min_samples": 1000           // Minimum PSMs for NCE modeling
    }
  }
}
```

### Fixed Internal Parameters
- **NCE Grid**: 30 values from 21.0 to 40.0
- **Breakpoint**: 500.0 m/z (defined by `NCE_MODEL_BREAKPOINT` constant)
- **Max Q-value**: 0.01 (1% FDR threshold)
- **Max Iterations**: 10 attempts to collect sufficient PSMs

## Integration with Pipeline

### Input Dependencies
- **SearchContext**: Contains calibrated mass error models from ParameterTuningSearch
- **Spectral Library**: Must be advanced library (not basic FragmentIndexLibrary)
- **Mass Spec Data**: DIA spectra for NCE calibration

### Output Products
1. **NCE Models**: Stored in SearchContext for each MS file
2. **QC Plots**: Combined PDF showing NCE calibration curves
3. **Diagnostic Logs**: Model parameters and data statistics

### Downstream Usage
The calibrated NCE models are used by:
- **QuadTuningSearch**: For accurate spectral matching during quad modeling
- **FirstPassSearch**: Initial PSM discovery with optimized NCE
- **SecondPassSearch**: Comprehensive search with calibrated parameters
- **All subsequent methods**: Benefit from improved spectral matching

## Common Issues and Troubleshooting

### Insufficient PSMs Warning
```
"Could not collect enough PSMs for NCE alignment. In 10 iterations collected N samples"
```
**Causes:**
- Low-complexity sample
- Restrictive filtering parameters
- Poor spectral library match to sample
- Instrument-specific issues

**Solutions:**
- Lower `min_samples` threshold
- Relax `min_score` or other filtering parameters
- Check mass calibration from ParameterTuningSearch
- Verify spectral library compatibility

### Model Fitting Failures
**Symptoms:** Missing NCE models in SearchContext
**Common Causes:**
- Insufficient data for optimization
- Extreme outliers in NCE values
- Numerical instability in LBFGS

**Debugging Steps:**
1. Check diagnostic logs for data distribution
2. Examine scatter plots for outliers
3. Verify charge state distribution
4. Consider manual parameter bounds

### Poor Model Quality
**Indicators:**
- Wide scatter around fitted lines
- Unrealistic parameter values
- Large residuals

**Improvements:**
- Increase `min_samples` for more robust fitting
- Filter extreme NCE outliers
- Check for systematic biases in data collection

## Performance Characteristics

### Computational Complexity
- **Memory Usage**: Moderate (stores PSMs and model parameters)
- **CPU Usage**: Moderate (LBFGS optimization per file)
- **I/O**: Low (works with in-memory data)
- **Scaling**: Linear with number of files, sublinear with PSM count

### Typical Runtime
- **Per File**: 10-60 seconds depending on PSM collection
- **Bottlenecks**: Library search iterations, not model fitting
- **Memory**: ~100MB for typical datasets

## Future Enhancements

### Potential Improvements
1. **Adaptive Grid Search**: Dynamic NCE range based on data
2. **Robust Regression**: M-estimators instead of least squares
3. **Multi-breakpoint Models**: Multiple piecewise segments
4. **Cross-Validation**: Model selection and validation
5. **Bayesian Approach**: Uncertainty quantification in predictions

### Advanced Features
- **Transfer Learning**: Use models from similar instruments
- **Temporal Calibration**: Account for instrument drift over time
- **Sample-Specific Adjustment**: Fine-tune for different sample types

## Conclusion

NceTuningSearch provides crucial NCE calibration that significantly improves spectral matching accuracy throughout the Pioneer.jl pipeline. Its piecewise linear model captures the nonlinear relationship between precursor properties and optimal collision energy while maintaining computational efficiency and interpretability.

The enhanced diagnostic logging and visualization make it easier to assess calibration quality and troubleshoot issues, ensuring robust performance across diverse experimental conditions.