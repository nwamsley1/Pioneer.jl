# NCE Tuning Sampling Optimization Plan

## Overview

This document outlines a comprehensive plan to optimize NCE tuning performance through intelligent scan sampling while maintaining calibration quality. The plan includes data point tracking and adaptive sampling strategies based on retention time distribution.

## Current Implementation Analysis

### Current NCE Tuning Flow

**Location**: `src/Routines/SearchDIA/SearchMethods/NceTuningSearch/NceTuningSearch.jl:207`

```julia
# Current approach - processes ALL MS2 scans
psms = library_search(spectra, search_context, params, ms_file_idx)
processed_psms = process_psms!(psms, spectra, search_context, params)
```

**Issues with Current Approach**:
1. **No data point visibility**: No tracking of how many PSMs/scans are used for NCE fitting
2. **Performance bottleneck**: Processing all scans can be slow for large files
3. **Potential oversampling**: More data than needed for reliable NCE calibration
4. **No adaptive strategy**: Fixed approach regardless of file size or data density

### Current Warning System

**Location**: `src/Routines/SearchDIA/SearchMethods/NceTuningSearch/NceTuningSearch.jl:211-214`

```julia
# Warn if insufficient PSMs for reliable NCE modeling
if nrow(processed_psms) < 100
    @user_warn "Low PSM count ($(nrow(processed_psms))) may result in unreliable NCE calibration. Consider lowering filtering thresholds."
end
```

## Proposed Solution: Simple Adaptive Scan Sampling

### Strategy Overview

**Core Concept**: Create a lightweight data structure similar to FilteredMassSpec that pre-orders scans by retention time bins and progressively samples them until sufficient data quality is achieved.

**Key Benefits**:
1. **Progressive Sampling**: Start with 25% of scans, escalate to 50%, then 100% if needed
2. **RT Distribution**: Ensure even distribution across retention time via round-robin from RT bins
3. **Quality Thresholds**: Define minimum PSM counts for reliable NCE calibration
4. **Data Point Tracking**: Log sample sizes and PSM counts at each level
5. **Memory Efficient**: Only stores scan indices, not actual spectral data
6. **Reuses Existing Logic**: Leverages proven RT binning from FilteredMassSpec

### Simple Implementation Design

**Core Data Structure**:
```julia
# New lightweight struct for NCE sampling
struct NceSamplingData
    # Pre-computed scan ordering
    scan_priority_order::Vector{UInt32}     # All MS2 scans ordered by RT bins
    rt_bin_assignments::Vector{Int8}        # Which RT bin each scan belongs to
    total_ms2_scans::Int                    # Total available scans

    # Sampling state
    current_sample_size::Int                # How many scans currently sampled
    target_sample_fractions::Vector{Float64} # [0.25, 0.5, 1.0] sampling levels
    current_level::Int                      # Which level we're on (1, 2, or 3)

    # Quality thresholds
    min_psms_excellent::Int                 # 500+ PSMs = excellent
    min_psms_good::Int                      # 200+ PSMs = good
    min_psms_adequate::Int                  # 100+ PSMs = adequate
end
```

**Algorithm Flow**:
1. **Initialization**: Create NceSamplingData using existing FilteredMassSpec RT binning logic
2. **Progressive Sampling**: Start with 25% of scans (round-robin from RT bins)
3. **Quality Assessment**: Run NCE search and evaluate PSM count/quality
4. **Adaptive Escalation**: If insufficient quality, escalate to 50%, then 100%
5. **Early Termination**: Stop as soon as "GOOD" quality achieved (≥200 PSMs)

### Implementation Plan

#### Phase 1: Core Sampling Infrastructure

**New Functions in `utils.jl`**:

```julia
"""
Create lightweight sampling data structure using existing FilteredMassSpec logic.
Reuses proven RT binning and scan prioritization algorithms.
"""
function create_nce_sampling_data(spectra::MassSpecData, n_rt_bins::Int = 15)
    # Reuse existing FilteredMassSpec RT binning logic
    rt_bin_assignments, _, _, _ = compute_rt_bins(spectra, n_rt_bins)

    # Get all MS2 scans and sort within bins by TIC (proxy for data quality)
    all_sorted_scans, bin_starts, bin_ends = sort_scans_by_peak_density(
        spectra, UInt8(2), rt_bin_assignments, n_rt_bins
    )

    # Create interleaved priority order (round-robin from RT bins)
    scan_priority_order = create_priority_order(all_sorted_scans, bin_starts, bin_ends)

    return NceSamplingData(
        UInt32.(scan_priority_order),
        rt_bin_assignments,
        length(scan_priority_order),
        0,                          # current_sample_size
        [0.25, 0.5, 1.0],          # target_sample_fractions
        1,                          # current_level
        500, 200, 100              # quality thresholds
    )
end

"""
Get the next batch of scans for current sampling level.
Returns empty vector if no more levels to try.
"""
function get_next_sample_scans(sampling_data::NceSamplingData)
    if sampling_data.current_level > length(sampling_data.target_sample_fractions)
        return UInt32[]  # No more levels to try
    end

    # Calculate target size for current level
    target_fraction = sampling_data.target_sample_fractions[sampling_data.current_level]
    target_size = round(Int, sampling_data.total_ms2_scans * target_fraction)

    # Get scans from priority order up to target size
    if target_size > sampling_data.current_sample_size
        new_scans = sampling_data.scan_priority_order[
            (sampling_data.current_sample_size + 1):target_size
        ]
        return new_scans
    else
        return UInt32[]  # Already have enough scans for this level
    end
end

"""
Advance to next sampling level and update state.
"""
function advance_sampling_level!(sampling_data::NceSamplingData, actual_scans_used::Int)
    sampling_data.current_sample_size = actual_scans_used
    sampling_data.current_level += 1
end
```

#### Phase 2: Main Adaptive Sampling Function

**Replace in `process_file!`**:

```julia
# OLD:
psms = library_search(spectra, search_context, params, ms_file_idx)
processed_psms = process_psms!(psms, spectra, search_context, params)

# NEW:
processed_psms, sampling_stats = adaptive_nce_sampling_simple(
    spectra, search_context, params, ms_file_idx, file_name
)
```

**New Function in `utils.jl`**:

```julia
"""
Simple adaptive sampling strategy for NCE tuning that progressively
samples 25% → 50% → 100% of scans until sufficient data quality.

# Algorithm
1. Create sampling plan using existing FilteredMassSpec RT binning
2. Progressive sampling loop: try 25%, 50%, 100%
3. Quality assessment after each level
4. Early termination when GOOD quality achieved (≥200 PSMs)
5. RT-distributed sampling via round-robin from bins

# Returns
(processed_psms::DataFrame, sampling_stats::NamedTuple)
"""
function adaptive_nce_sampling_simple(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NceTuningSearchParameters,
    ms_file_idx::Int64,
    file_name::String
)

    # Step 1: Create sampling plan
    sampling_data = create_nce_sampling_data(spectra)

    if sampling_data.total_ms2_scans == 0
        @warn "No MS2 scans found in $file_name for NCE tuning"
        return DataFrame(), (total_scans=0, sampled_scans=0, final_psms=0, level=:none)
    end

    @info "NCE Tuning $file_name: $(sampling_data.total_ms2_scans) total MS2 scans available"

    # Step 2: Progressive sampling loop
    best_psms = DataFrame()
    final_stats = (total_scans=0, sampled_scans=0, final_psms=0, level=:none)

    while sampling_data.current_level <= 3
        # Get next batch of scans to try
        additional_scans = get_next_sample_scans(sampling_data)

        if isempty(additional_scans)
            break  # No more scans to add at this level
        end

        # Current level info
        level_names = [:quarter, :half, :full]
        level_descriptions = ["25%", "50%", "100%"]
        current_level_name = level_names[sampling_data.current_level]
        current_description = level_descriptions[sampling_data.current_level]

        @info "NCE Tuning $file_name: Trying $current_description sampling ($(length(additional_scans)) additional scans)"

        # Step 3: Perform library search on current scan set
        # Create scan indices for library search (all scans up to current level)
        current_scan_count = sampling_data.current_sample_size + length(additional_scans)
        all_scans_to_search = sampling_data.scan_priority_order[1:current_scan_count]

        # Use existing library search but with scan filtering
        psms = library_search_with_scan_filter(
            spectra, search_context, params, ms_file_idx, all_scans_to_search
        )
        processed_psms = process_psms!(psms, spectra, search_context, params)

        psm_count = nrow(processed_psms)

        # Step 4: Evaluate quality
        quality_level = if psm_count >= sampling_data.min_psms_excellent
            "EXCELLENT"
        elseif psm_count >= sampling_data.min_psms_good
            "GOOD"
        elseif psm_count >= sampling_data.min_psms_adequate
            "ADEQUATE"
        else
            "POOR"
        end

        @info "NCE Tuning $file_name: $current_description sampling → $psm_count PSMs ($quality_level)"

        # Step 5: Check if we should continue or stop
        best_psms = processed_psms  # Always keep the latest (more data)
        final_stats = (
            total_scans = sampling_data.total_ms2_scans,
            sampled_scans = current_scan_count,
            final_psms = psm_count,
            level = current_level_name
        )

        # Stop if we have good quality or this is the final level
        if psm_count >= sampling_data.min_psms_good || sampling_data.current_level >= 3
            if psm_count >= sampling_data.min_psms_good
                @info "NCE Tuning $file_name: Sufficient quality achieved with $current_description sampling"
            end
            break
        else
            @info "NCE Tuning $file_name: Need more data ($psm_count < $(sampling_data.min_psms_good)), escalating to next level"
        end

        # Advance to next level
        advance_sampling_level!(sampling_data, current_scan_count)
    end

    # Step 6: Final quality warning if needed
    final_psm_count = nrow(best_psms)
    if final_psm_count < sampling_data.min_psms_adequate
        @warn "NCE Tuning $file_name: LOW PSM COUNT ($final_psm_count) may result in unreliable NCE calibration"
    end

    return best_psms, final_stats
end
```

#### Phase 3: Library Search Integration

**Add scan filtering to existing library search**:

```julia
"""
Enhanced library search that accepts a scan filter.
Integrates with existing NCE tuning pipeline seamlessly.
"""
function library_search_with_scan_filter(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NceTuningSearchParameters,
    ms_file_idx::Int64,
    scan_indices::Vector{UInt32}
)
    # Create a view/filter of spectra containing only the specified scans
    # This could be implemented as:
    # 1. A filtered wrapper around the original spectra
    # 2. Or modify existing library search to accept scan filter
    # 3. Or create temporary FilteredMassSpecData with only these scans

    # For now, integrate with existing library_search by:
    # - Adding scan_filter parameter to LibrarySearchNceTuning
    # - Or using FilteredMassSpecData as a lightweight wrapper

    filtered_spectra = create_scan_filtered_view(spectra, scan_indices)

    return library_search(filtered_spectra, search_context, params, ms_file_idx)
end

"""
Create a lightweight view of spectra containing only specified scans.
Could reuse FilteredMassSpecData or create simpler wrapper.
"""
function create_scan_filtered_view(spectra::MassSpecData, scan_indices::Vector{UInt32})
    # Implementation could:
    # 1. Use existing FilteredMassSpecData constructor with specific scan list
    # 2. Create new lightweight wrapper type
    # 3. Modify existing MassSpecData interface to accept scan filters

    # Simplest approach: use existing FilteredMassSpecData
    return FilteredMassSpecData(
        spectra,
        max_scans = length(scan_indices),
        target_ms_order = UInt8(2),
        # Override the scan selection to use our specific indices
        scan_selection_override = scan_indices
    )
end
```

## Performance Benefits

### Expected Improvements

1. **Speed**: 50-75% reduction in NCE tuning time for large files
   - Most files will achieve good quality with 25% sampling
   - Only files with poor data quality require full sampling
   - RT-distributed sampling ensures representative coverage

2. **Scalability**: Better handling of very large DIA files
   - Memory efficient: only stores scan indices, not spectral data
   - Constant memory overhead regardless of file size
   - Progressive scaling avoids loading unnecessary data

3. **Quality**: Maintained calibration quality with intelligent sampling
   - Round-robin RT bin sampling ensures temporal coverage
   - TIC-based prioritization within bins selects high-quality scans
   - Adaptive escalation guarantees sufficient data for reliable models

4. **Visibility**: Clear reporting of data usage and quality
   - Logs sampling level used for each file
   - Reports PSM counts and quality assessments
   - Tracks total vs. sampled scan counts

### Simple Implementation Advantages

1. **Reuses Existing Logic**: Leverages proven RT binning from FilteredMassSpec
2. **Minimal Code Changes**: Small modifications to existing NCE tuning pipeline
3. **Memory Efficient**: Lightweight data structure with only scan indices
4. **Backwards Compatible**: Fallback to 100% sampling maintains current behavior
5. **Easy Integration**: Uses existing library search with filtered scan views
6. **Performance Predictable**: 3-level sampling with clear quality thresholds

### Fallback Safety

- Always escalates to full scan sampling if needed
- Maintains current quality thresholds and warnings
- Provides clear logging of sampling decisions
- No change in NCE model fitting algorithms
- Graceful degradation for poor quality data

## Implementation Timeline

### Phase 1: Core Infrastructure (Week 1)
- Implement `NceSamplingData` struct and utilities
- Add `create_nce_sampling_data` function
- Reuse existing FilteredMassSpec RT binning logic
- Unit tests for sampling data structure

### Phase 2: Adaptive Sampling (Week 2)
- Implement `adaptive_nce_sampling_simple` function
- Add progressive sampling loop (25% → 50% → 100%)
- Integrate quality assessment and early termination
- Test with representative datasets

### Phase 3: Library Search Integration (Week 3)
- Implement `library_search_with_scan_filter`
- Create scan-filtered view mechanism
- Integration with existing NCE tuning pipeline
- Performance benchmarking

### Phase 4: Testing and Validation (Week 4)
- Compare NCE models from sampled vs full data
- Validate RT distribution coverage
- Performance improvement measurement
- Documentation and final integration

## Testing Strategy

### Unit Tests
- `NceSamplingData` creation and state management
- `get_next_sample_scans` with various sampling levels
- `advance_sampling_level!` state transitions
- RT bin distribution validation
- Quality threshold evaluation

### Integration Tests
- `adaptive_nce_sampling_simple` with representative datasets
- Compare NCE models from sampled vs full data
- Validate RT distribution coverage across sampling levels
- Library search integration with scan filtering
- Performance benchmarking on large files

### Validation Criteria
- NCE model correlation > 0.95 between sampling levels
- Maintained charge state coverage across RT range
- 50-75% performance improvement for large files
- Early termination rate >70% at 25% sampling for good quality data
- No regression in final NCE model quality

### Test Datasets
- **Small files** (< 1000 scans): Should use 25% sampling and achieve good quality
- **Medium files** (1000-5000 scans): May need 50% sampling
- **Large files** (> 5000 scans): Should show significant performance gains
- **Poor quality files**: Should escalate to 100% sampling with appropriate warnings

## Code References

**Files to Modify**:
1. `src/Routines/SearchDIA/SearchMethods/NceTuningSearch/NceTuningSearch.jl:207` - Replace PSM collection call
2. `src/Routines/SearchDIA/SearchMethods/NceTuningSearch/utils.jl` - Add new sampling functions
3. `src/structs/FilteredMassSpecData.jl` - Reuse existing RT binning functions

**Key Functions to Add**:
- `create_nce_sampling_data()` - Initialize sampling data structure
- `get_next_sample_scans()` - Get scans for current sampling level
- `advance_sampling_level!()` - Progress to next sampling level
- `adaptive_nce_sampling_simple()` - Main adaptive sampling function
- `library_search_with_scan_filter()` - Library search with scan filtering

**Key Functions to Reuse**:
- `compute_rt_bins()` (FilteredMassSpecData.jl:83) - RT binning logic
- `sort_scans_by_peak_density()` (FilteredMassSpecData.jl:115) - Scan prioritization
- `create_priority_order()` (FilteredMassSpecData.jl:189) - Round-robin ordering

**Integration Points**:
- `process_file!` (NceTuningSearch.jl:207) - Replace library_search call
- Add `NceSamplingData` struct definition
- Minimal changes to existing pipeline

## Summary

This simplified plan provides comprehensive NCE tuning optimization through:

1. **Intelligent Sampling**: Progressive 25% → 50% → 100% sampling with RT distribution
2. **Proven Algorithms**: Reuses FilteredMassSpec RT binning and scan prioritization
3. **Simple Integration**: Minimal changes to existing codebase
4. **Performance Gains**: 50-75% reduction in processing time for most files
5. **Quality Assurance**: Maintains calibration quality with adaptive escalation
6. **Memory Efficiency**: Lightweight data structure with scan indices only

The approach strikes an optimal balance between implementation simplicity, performance gains, and maintenance of data quality for reliable NCE calibration.