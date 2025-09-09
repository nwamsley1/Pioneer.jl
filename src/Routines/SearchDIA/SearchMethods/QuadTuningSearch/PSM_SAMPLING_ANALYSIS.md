# PSM Sampling Strategy in QuadTuningSearch

## Overview
QuadTuningSearch uses an iterative sampling strategy to collect PSMs for quadrupole transmission modeling. Unlike ParameterTuningSearch which samples scans, QuadTuningSearch processes **all scans** but collects PSMs iteratively until sufficient data is gathered.

## Current Implementation (Fixed - No Longer Uses Loop)

### PSM Collection Process

The `collect_psms` function in `utils.jl` has been simplified to run a single collection pass:

```julia
# Single library search processes all scans
psms = library_search(spectra, search_context, params, ms_file_idx)

# Process and filter PSMs
processed_psms = process_initial_psms(psms, spectra, search_context)
quad_psms = perform_quad_transmission_search(...)
# Apply quality filters and return final PSMs
```

### Key Characteristics

1. **No Scan Sampling**: QuadTuningSearch processes ALL MS2 scans in a single pass
2. **Single Collection**: One library search followed by processing and filtering
3. **Warning System**: Warns if insufficient PSMs are collected rather than attempting multiple iterations
4. **Comprehensive Logging**: Detailed progress logging at each processing step

### PSM Filtering Criteria

The `filter_quad_psms` function applies strict criteria:
- **Isotope Constraint**: Only M0 and M1 isotopes (`iso_idx < 3`)
- **Charge State**: Only 2+ charges (`charge == 2`)
- **Fragment Matches**: Minimum number of matched fragments (`n_matches >= min_quad_tuning_fragments`)
- **Non-zero Weight**: Must have positive abundance (`weight > 0`)

### library_search Function

The generic `library_search` function used by QuadTuningSearch:
1. Uses the full fragment index (not presearch index)
2. Processes all scans according to thread partitioning
3. No scan sampling - processes entire dataset

## Comparison with v0.2.1

### What's the Same
- **10-iteration loop** structure unchanged
- **Unique precursor tracking** unchanged
- **PSM filtering criteria** unchanged (M0/M1, charge 2+, min fragments)
- **No scan sampling** - both versions process all scans

### What Changed

#### 1. Loop Removal (Major Fix)
- **v0.2.1**: Ineffective 10-iteration loop that only ran once due to logic flaw
- **Current**: Single-pass collection with proper logic and comprehensive logging

#### 2. Default Thresholds
- **v0.2.1**: `min_quad_tuning_psms = 1000` (in example docs)
- **Current**: `min_quad_tuning_psms = 5000` (default in code)
- Both allow JSON override, but default increased 5x

#### 3. Error Handling
- **Current**: Better error handling, continues on failure, warns about insufficient PSMs
- **v0.2.1**: Would throw errors and stop

#### 4. Logging
- **Current**: Comprehensive logging at each step showing PSM counts
- **v0.2.1**: Minimal logging

## Comparison with ParameterTuningSearch

| Aspect | QuadTuningSearch | ParameterTuningSearch |
|--------|-----------------|----------------------|
| **Scan Sampling** | NO - processes all scans | YES - samples subset |
| **Sampling Strategy** | Single-pass PSM collection | Convergence-based scan sampling |
| **Data Structure** | Regular MassSpecData | FilteredMassSpecData (scan subset) |
| **Iterations** | Single pass | Dynamic until convergence |
| **Early Stop** | Warning if insufficient PSMs | When score converges |
| **Precursor Filtering** | All precursors (no duplicates) | All precursors |

## Why No Scan Sampling?

QuadTuningSearch doesn't sample scans because:

1. **Need Isotope Patterns**: Requires M0/M1 isotope pairs which are less common
2. **Charge State Restriction**: Only uses 2+ charges, further limiting data
3. **Strict Quality Filters**: Multiple criteria reduce available PSMs significantly
4. **Statistical Power**: Needs sufficient data for reliable transmission modeling

## PSM Count Requirements

The minimum PSM requirement ensures:
- Sufficient isotope ratio measurements
- Coverage across m/z range
- Statistical power for spline fitting
- Robustness to outliers

Current default of 5000 PSMs provides better coverage than v0.2.1's 1000, especially for:
- Wide isolation windows
- Complex samples
- Instruments with variable transmission

## Implementation Details

### Thread Partitioning
Uses `partition_scans` to distribute work across threads:
- Balances scan density
- Maintains data locality
- No scan sampling at this level

### Memory Management
- Processes all scans but tracks unique precursors
- Iterative approach prevents memory explosion
- Reuses search data structures across iterations

## Recommendations

1. **Keep Current Approach**: No scan sampling is appropriate for quad tuning
2. **Configurable Threshold**: The 5000 PSM default is reasonable but should remain configurable
3. **Consider Adaptive Stopping**: Could stop earlier if transmission model converges
4. **Add Progress Reporting**: Report PSM collection progress to user

## Summary

QuadTuningSearch's PSM sampling strategy has been corrected and is fundamentally different from ParameterTuningSearch:
- **No scan sampling** - needs all available data
- **Single-pass collection** (fixed from ineffective loop)
- **Strict filtering** for high-quality isotope measurements
- **Higher PSM requirements** (5000 vs 1000 in v0.2.1) for better model fitting
- **Comprehensive logging** at each processing step

### Key Fix Applied
The previous 10-iteration loop was logically flawed because:
1. `library_search` processes all scans deterministically, returning the same PSMs each time
2. After iteration 1, all precursors are in `unique_precursors` set
3. Iteration 2+ filter out all PSMs as duplicates, causing immediate return
4. **Result**: Loop effectively only ran once anyway

The fix removes this unnecessary complexity and makes the single-pass behavior explicit with proper logging and warnings.

This approach is well-suited for quadrupole transmission modeling which requires:
- Isotope ratio measurements (M0/M1)
- Specific charge states (2+) 
- High data quality for accurate spline fitting
- Coverage across the full m/z range