# MS1 vs MS2 Mass Tolerance Estimation in Pioneer.jl

## Overview

Pioneer.jl employs two distinct approaches for mass tolerance estimation: one for MS2 fragments (in ParameterTuningSearch) and one for MS1 precursors (in FirstPassSearch). This analysis examines the key differences, sampling strategies, and potential improvements.

## MS2 Mass Tolerance Estimation (ParameterTuningSearch)

### Sampling Strategy

**Data Structure: FilteredMassSpecData**
- **Location**: `src/structs/FilteredMassSpecData.jl:24-73`
- **Purpose**: Intelligent scan selection based on RT binning and peak density
- **Key Features**:
  - RT-based stratified sampling across time dimension
  - Peak density-based prioritization for representative scans
  - Configurable scan limits (500 → 10,000 → 80,000 scan scaling)
  - TopN peak filtering to reduce noise

**Scan Selection Algorithm** (`src/structs/FilteredMassSpecData.jl:83-100`):
```julia
function compute_rt_bins(spectra::MassSpecData, n_bins::Int = 15)
    rt_values = getRetentionTimes(spectra)
    rt_min, rt_max = extrema(rt_values)
    bin_width = (rt_max - rt_min) / n_bins
    # Assigns each scan to RT bin for stratified sampling
end
```

**Multi-Phase Convergence Strategy** (`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl:300-350`):
1. **Phase 1**: Zero bias exploration
2. **Phase 2**: Positive bias shift (+max_tolerance)
3. **Phase 3**: Negative bias shift (-max_tolerance)
4. **Multi-Score Thresholds**: `[22, 17, 15]` tried sequentially within each phase
5. **Scan Scaling**: 10x scaling between failed attempts (500 → 5,000 → 50,000 scans)

**Mass Error Model Fitting** (`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:fit_mass_err_model`):
- **Fragment-level analysis**: Collects PPM errors from matched fragments
- **Exponential distribution fitting**: Separate modeling for positive/negative tails
- **Intensity filtering**: Removes low-quality fragment matches
- **Robust estimation**: Uses quantile-based fallbacks if exponential fitting fails

## MS1 Mass Tolerance Estimation (FirstPassSearch)

### Sampling Strategy

**Data Structure: Raw MassSpecData**
- **Location**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:474-481`
- **Sampling Approach**: Simple intensity-based selection
- **Sample Size**: Fixed 3,000 most intense PSMs with q-value ≤ 0.001

**Sample Selection** (`src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:474-481`):
```julia
temp_psms = temp_psms[temp_psms[!,:q_value].<=0.001,:]
most_intense = sortperm(temp_psms[!,:log2_summed_intensity], rev = true)
ms1_errs = vcat(
    mass_error_search(
        spectra,
        temp_psms[most_intense[1:(min(3000, length(most_intense)))],:scan_idx],
        temp_psms[most_intense[1:(min(3000, length(most_intense)))],:precursor_idx],
        # ... other parameters
    )
)
```

**Mass Error Extraction** (`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:mass_error_search`):
- **MS1 chromatogram extraction**: Finds nearest MS1 scans to MS2 identifications
- **Precursor isotope matching**: Matches theoretical isotope patterns to observed peaks
- **PPM error calculation**: Direct mass error calculation for each isotope peak

**Model Fitting** (`src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:mass_err_ms1`):
```julia
function mass_err_ms1(ppm_errs::Vector{Float32}, params::P) where {P<:FragmentIndexSearchParameters}
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err

    # Separate positive and negative errors for asymmetric modeling
    r_errs = abs.(ppm_errs[ppm_errs .> 0.0f0])
    l_errs = abs.(ppm_errs[pmp_errs .< 0.0f0])

    # Fit exponential distribution to each tail
    r_bound = quantile(Exponential(...), 1 - frag_err_quantile)
    l_bound = quantile(Exponential(...), 1 - frag_err_quantile)

    return MassErrorModel(Float32(mass_err), (Float32(abs(l_bound)), Float32(abs(r_bound))))
end
```

## Key Differences

### 1. Sampling Sophistication

**MS2 (ParameterTuningSearch)**:
- ✅ **Intelligent FilteredMassSpecData**: RT-stratified sampling ensures representative coverage
- ✅ **Adaptive scaling**: Progressive scan count increase (500 → 80,000)
- ✅ **Multi-phase convergence**: Systematic bias exploration
- ✅ **Peak density optimization**: Priority-based scan selection

**MS1 (FirstPassSearch)**:
- ❌ **Simple intensity ranking**: No RT stratification or temporal coverage
- ❌ **Fixed sample size**: Always 3,000 scans regardless of data size
- ❌ **No convergence testing**: Single-shot estimation
- ❌ **No bias exploration**: No systematic bias shift testing

### 2. Convergence Strategy

**MS2 Features**:
- **Multi-phase approach**: Tests different bias scenarios systematically
- **Tolerance expansion**: Progressive widening when needed
- **Quality thresholds**: Multiple score cutoffs within each phase
- **Post-convergence optimization**: 50% tolerance expansion test

**MS1 Limitations**:
- **Single-shot estimation**: No iterative refinement
- **No bias testing**: Assumes zero systematic bias
- **No expansion testing**: Fixed tolerance estimation
- **No convergence criteria**: Accepts whatever data is available

### 3. Data Quality and Coverage

**MS2 Advantages**:
- **Fragment intensity filtering**: Removes low-quality matches
- **Multiple scoring thresholds**: Adapts to data quality
- **Scan density consideration**: Balances coverage vs quality
- **Temporal stratification**: Ensures instrument drift capture

**MS1 Issues**:
- **Intensity bias**: Favors high-intensity precursors only
- **Temporal clustering**: May miss RT-dependent systematic errors
- **Limited sample diversity**: Fixed 3,000 scan limit regardless of dataset size
- **No quality filtering**: Includes potentially poor MS1 matches

### 4. Statistical Robustness

**MS2 Strengths**:
- **Robust exponential fitting**: Handles asymmetric error distributions
- **Quantile-based fallbacks**: Graceful degradation when fitting fails
- **Large sample sizes**: Can utilize up to 80,000 scans for robust estimation
- **Outlier removal**: Multiple filtering stages

**MS1 Weaknesses**:
- **Simple median estimation**: Limited robustness to outliers
- **Small sample constraint**: Fixed 3,000 scan limit
- **No outlier detection**: MAD-based filtering only after initial estimation
- **Single estimation attempt**: No fallback strategies

## Potential Improvements for MS1 Mass Tolerance Estimation

### 1. Adopt FilteredMassSpecData for MS1

**Benefits**:
- **RT stratification**: Ensure representative sampling across chromatographic time
- **Adaptive scaling**: Start with smaller samples and expand as needed
- **Peak density optimization**: Select scans with optimal MS1 peak characteristics

**Implementation**:
```julia
# In FirstPassSearch
filtered_ms1_spectra = FilteredMassSpecData(
    spectra,
    max_scans = initial_scan_count,
    target_ms_order = 1,  # MS1 scans only
    topn = nothing,       # Keep all MS1 peaks
    min_intensity = 0.0,
    n_rt_bins = 15
)
```

### 2. Multi-Phase MS1 Convergence

**Proposed Algorithm**:
```julia
function estimate_ms1_tolerance_robust(spectra, psms, params)
    for phase_bias in [0.0, +initial_tol, -initial_tol]
        for score_threshold in [high_confidence, medium_confidence, low_confidence]
            for scan_count in [1000, 5000, 15000]
                filtered_data = create_filtered_ms1_data(spectra, scan_count)
                ms1_errors = collect_ms1_errors(filtered_data, psms, phase_bias)

                if check_ms1_convergence(ms1_errors, min_sample_size)
                    return fit_ms1_model(ms1_errors, phase_bias)
                end
            end
        end
    end
    return fallback_ms1_model(params)
end
```

### 3. Enhanced MS1 Error Collection

**Current Issues**:
- Uses MS2-based PSM identifications to find MS1 scans
- Limited to 3,000 highest intensity PSMs
- No systematic RT coverage

**Proposed Improvements**:
```julia
function collect_ms1_errors_improved(filtered_spectra, target_precursors, ms1_tolerance)
    # Stratified sampling across RT bins
    rt_bins = compute_rt_bins(filtered_spectra, n_bins=20)

    ms1_errors = Float32[]
    for bin in rt_bins
        bin_scans = get_ms1_scans_in_bin(filtered_spectra, bin)
        for scan_idx in bin_scans
            # Extract all precursor candidates in scan
            candidates = find_precursor_candidates(filtered_spectra, scan_idx, target_precursors)

            # Match to theoretical isotope patterns
            for candidate in candidates
                isotope_errors = match_isotope_pattern(candidate, ms1_tolerance)
                append!(ms1_errors, isotope_errors)
            end
        end
    end

    return ms1_errors
end
```

### 4. MS1-Specific Convergence Criteria

**Proposed Criteria**:
```julia
function check_ms1_convergence(ms1_errors, min_precursors=500)
    return (
        length(ms1_errors) >= min_precursors &&
        abs(median(ms1_errors)) < initial_tolerance/4 &&
        mad(ms1_errors) < initial_tolerance &&
        length(unique_precursors(ms1_errors)) >= min_precursors/2
    )
end
```

### 5. Integration with Existing Infrastructure

**Advantages of Using FilteredMassSpecData**:
- **Code reuse**: Leverages existing intelligent sampling logic
- **Consistent architecture**: Same pattern as ParameterTuningSearch
- **Memory efficiency**: In-memory filtered view rather than full data loading
- **Scalability**: Handles datasets of any size gracefully

**Modified FirstPassSearch Flow**:
```julia
function process_chunk!(method::FirstPassSearch, batch_id::Int, thread_id::Int)
    # 1. Create FilteredMassSpecData for MS1 and MS2 separately
    filtered_ms1 = FilteredMassSpecData(spectra, target_ms_order=1, max_scans=scan_count)
    filtered_ms2 = FilteredMassSpecData(spectra, target_ms_order=2, max_scans=scan_count)

    # 2. Perform MS2-based PSM identification (existing logic)
    psms = library_search(filtered_ms2, search_context, params, ms_file_idx)

    # 3. Robust MS1 tolerance estimation using filtered MS1 data
    ms1_model = estimate_ms1_tolerance_robust(filtered_ms1, psms, params)

    # 4. Store results
    setMs1MassErrorModel!(search_context, ms_file_idx, ms1_model)
end
```

## Performance Considerations

### Memory Usage
- **Current MS1**: Minimal - processes only 3,000 scans
- **Proposed MS1**: Moderate - FilteredMassSpecData scales linearly with scan count
- **Trade-off**: Slightly higher memory for significantly better tolerance estimation

### Computational Cost
- **Current MS1**: Very fast - single-shot estimation
- **Proposed MS1**: Moderate - multi-phase convergence with early termination
- **Benefit**: More robust estimation prevents downstream search failures

### Implementation Complexity
- **FilteredMassSpecData adoption**: Low - existing infrastructure
- **Multi-phase logic**: Medium - adaptation of ParameterTuningSearch patterns
- **Integration**: Low - fits existing FirstPassSearch architecture

## Conclusion

The current MS1 mass tolerance estimation in FirstPassSearch is significantly less sophisticated than the MS2 approach in ParameterTuningSearch. The main limitations are:

1. **Simplistic sampling**: Intensity-based selection without RT stratification
2. **No convergence testing**: Single-shot estimation without robustness checks
3. **Fixed sample size**: Doesn't adapt to dataset characteristics
4. **Limited bias exploration**: Assumes zero systematic error

**Adopting FilteredMassSpecData for MS1 tolerance estimation would provide**:
- ✅ Representative sampling across chromatographic time
- ✅ Adaptive scaling for different dataset sizes
- ✅ Consistent architecture with existing ParameterTuningSearch
- ✅ Improved robustness through multi-phase convergence
- ✅ Better handling of systematic mass errors

The implementation would be straightforward given the existing infrastructure and would significantly improve MS1 tolerance estimation robustness without major architectural changes.

## Code References

**MS2 Implementation**:
- `src/structs/FilteredMassSpecData.jl:24-400` - Intelligent sampling infrastructure
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl:200-400` - Multi-phase convergence
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:fit_mass_err_model` - Robust model fitting

**MS1 Implementation**:
- `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:474-481` - Simple sampling
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:mass_err_ms1` - Basic model fitting
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl:mass_error_search` - MS1 error collection