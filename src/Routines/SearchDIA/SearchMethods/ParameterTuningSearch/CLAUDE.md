# CLAUDE.md - ParameterTuningSearch Deep Dive

This file provides a comprehensive understanding of how ParameterTuningSearch works in Pioneer.jl, documenting the intricate details of the parameter tuning process.

## Purpose

ParameterTuningSearch is the first search method in the Pioneer.jl pipeline, responsible for establishing optimal fragment mass tolerance and retention time alignment models. It serves as the foundation for all subsequent search methods by calibrating critical parameters based on the actual data characteristics.

## Core Workflow

### 1. Iterative PSM Collection (`collect_psms`)

The search performs up to `max_presearch_iters` iterations to collect sufficient high-quality PSMs:

```julia
for i in 1:getMaxPresearchIters(params)
    new_psms = library_search(spectra, search_context, params, ms_file_idx)
    # Add columns and append to existing PSMs
    # Check if we have enough PSMs passing FDR threshold
    if filter_and_score_psms!(psms, params, search_context) >= getMinPsms(params)
        break
    end
end
```

**Key Points:**
- Uses a specialized `library_search` that employs the presearch fragment index (coarser binning)
- Each iteration adds more PSMs to the pool
- Stops when sufficient PSMs pass the FDR threshold (default: q-value ≤ 0.01)

### 2. Library Search Process

The `library_search` function orchestrates a two-phase search:

#### Phase 1: Fragment Index Search (`searchFragmentIndex`)
- Scans through spectra and queries the fragment index
- Uses initial mass tolerance (default: 20 ppm) with mass offset from previous iteration
- Applies sampling rate (default: 0.02) to reduce computational load
- Scores precursors based on matched fragment count
- Filters by minimum score (default: 22)

#### Phase 2: Detailed PSM Generation (`getPSMS`)
- Takes high-scoring precursors from Phase 1
- Performs detailed fragment matching using `matchPeaks`
- Builds design matrices for quantification
- Calculates spectral scores (entropy, spectral contrast, etc.)
- Generates PSM records with all features

### 3. PSM Scoring and Filtering (`filter_and_score_psms!`)

**Scoring Process:**
- Uses probit regression with features: entropy_score, city_block, scribe, spectral_contrast, y_count, error, TIC, intercept
- Calculates q-values with FDR control
- Filters PSMs by q-value threshold (default: 0.01)
- Selects best PSM per precursor based on probability score
- **Important:** Only uses target PSMs for parameter fitting (excludes decoys)

### 4. RT Model Fitting (`fit_irt_model`)

**Two-stage spline fitting:**
1. Initial fit: `UniformSpline(irt_predicted, rt, degree=3, n_knots=5)`
2. Outlier removal: Remove points beyond `MAD * outlier_threshold`
3. Final fit: Refit spline on cleaned data

**Returns:**
- SplineRtConversionModel
- Valid RT/iRT pairs used for fitting
- Median Absolute Deviation (MAD) for error estimation

### 5. Mass Error Model Fitting

**Process:**
1. Collect matched fragments using current RT model
2. Use specialized `mass_error_search` that:
   - Only considers top 5 fragments per precursor
   - Applies RT filtering
   - Collects PPM errors for all matches
3. Fit mass error model:
   - Calculate median mass offset
   - Determine left/right tolerances using quantiles (default: 1% and 99%)

### 6. Convergence Loop in `process_file!`

The current implementation has a convergence loop with up to 5 attempts:

```julia
while n_attempts < 5
    # Set current mass error model
    # Collect PSMs
    # Fit RT model
    # Fit mass error model
    
    if abs(getMassOffset(mass_err_model)) > (getFragTolPpm(params)/4)
        # Mass offset too large - retry with updated offset
    elseif (getLeftTol(mass_err_model) + getRightTol(mass_err_model)) > init_mass_tol
        # Tolerance too wide - increase initial tolerance by 1.5x
    else
        # Converged - break
    end
end
```

## Key Data Structures

### ParameterTuningSearchParameters
- `frag_tol_ppm`: Initial fragment tolerance (default: 20.0)
- `min_psms`: Minimum PSMs required (default: 3500)
- `sample_rate`: Fraction of spectra to sample (default: 0.02)
- `max_presearch_iters`: Maximum collection iterations (default: 10)
- `min_index_search_score`: Minimum fragment matches (default: 22)

### ParameterTuningSearchResults
- `mass_err_model`: MassErrorModel with offset and tolerances
- `rt_to_irt_model`: SplineRtConversionModel
- `irt`, `rt`: Vectors of RT/iRT pairs used for model
- `ppm_errs`: Fragment mass errors for QC plots

## Critical Implementation Details

### Fragment Index Search Strategy

The fragment index uses a hierarchical structure:
1. RT bins (grouped by retention time)
2. Fragment bins (grouped by m/z within RT bins)
3. Individual fragments

The search uses exponential bound expansion for efficiency:
- Start with narrow m/z window
- Exponentially expand bounds until all potential matches found
- Use branchless binary search within bounds

### Mass Error Handling

Two types of mass error corrections:
1. **Forward correction** (theoretical → empirical): Used in peak matching
   - `getMzBounds(mem, theoretical_mz)` returns bounds for empirical peaks
2. **Reverse correction** (empirical → theoretical): Used in fragment index
   - `getMzBoundsReverse(mem, corrected_empirical_mz)` returns bounds for theoretical fragments

### Thread Safety

- Each thread has its own `SearchDataStructures` instance
- Results are collected independently per thread
- Final merging happens in single-threaded context

## Common Issues and Edge Cases

### 1. Mass Bias Exceeding Initial Tolerance
**Symptom:** Actual mass bias is larger than the search window
**Current Handling:** Iterative refinement with updated offset
**Problem:** May require multiple iterations, reducing efficiency

### 2. Inadequate Edge Sampling
**Symptom:** True tolerance is close to initial guess
**Current Handling:** Increase tolerance by 1.5x if fitted tolerance ≥ initial
**Problem:** May not adequately sample distribution tails

### 3. Large Isolation Windows
**Symptom:** Too many precursor candidates per spectrum
**Current Handling:** None - relies on minimum score filtering
**Problem:** Slow performance, many false positives

### 4. Low PSM Yield
**Symptom:** Cannot collect enough PSMs even after max iterations
**Current Handling:** Proceeds with available PSMs or fails
**Problem:** Poor model fits, unreliable parameters

## QC Output

Two main QC plots are generated:
1. **RT Alignment Plot**: Scatter plot of RT vs iRT with fitted spline
2. **Mass Error Plot**: Histogram of mass errors with model boundaries

These are saved per file and then merged into combined PDFs.

## Integration with Downstream Methods

The calibrated parameters are stored in SearchContext and used by:
- **FirstPassSearch**: Uses mass error model and RT conversion
- **SecondPassSearch**: Refined search with calibrated parameters
- **All subsequent methods**: Rely on RT alignment for chromatogram extraction

## Performance Characteristics

- **Memory Usage**: Moderate - stores PSMs and fragment matches
- **CPU Usage**: High during fragment matching and model fitting
- **I/O**: Minimal - works with in-memory data
- **Typical Runtime**: 30-60 seconds per file depending on sampling rate

## Future Improvements Needed

1. **Adaptive Sampling**: Dynamically adjust sample rate based on PSM yield
2. **TopN Peak Filtering**: For large tolerance scenarios
3. **Better Convergence Criteria**: More sophisticated than current simple checks
4. **Diagnostic Reporting**: Track and report convergence metrics
5. **Fallback Strategies**: Handle files that fail to converge gracefully