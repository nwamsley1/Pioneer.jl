# Scan Sampling Improvement Plan for ParameterTuningSearch

## Overview
Replace the current sample_rate approach with explicit scan count sampling to provide better control over PSM collection and improve convergence in challenging datasets.

## Current Issues
1. Sample rate (0.02) is abstract - unclear how many scans are actually being processed
2. No visibility into scan counts vs PSM yield
3. Inefficient retry strategy when insufficient PSMs are found
4. No explicit control over computational cost

## Proposed Changes

### 1. Modify FilteredMassSpecData Constructor

**Current:**
```julia
FilteredMassSpecData(
    spectra,
    params.sample_rate,  # e.g., 0.02
    1000,  # topn_peaks
    target_ms_order = UInt8(2)
)
```

**Proposed:**
```julia
FilteredMassSpecData(
    spectra,
    max_scans = 2500,  # Explicit scan count
    topn_peaks = 1000,
    target_ms_order = UInt8(2)
)
```

**Implementation Details:**
- Change constructor to accept `max_scans::Int` instead of `sample_rate::Float32`
- Randomly sample up to `max_scans` MS2 spectra
- Store the actual number of sampled scans for reporting

### 2. Update ParameterTuningSearchParameters

**Add new parameters:**
```julia
struct ParameterTuningSearchParameters
    # ... existing fields ...
    initial_scan_count::Int64  # Default: 2500
    expanded_scan_count::Int64  # Default: 5000
    # Remove or deprecate: sample_rate
end
```

### 3. Modify process_file! Logic

**New flow:**
```julia
function process_file!(...)
    # ... initialization ...
    
    # Create initial filtered data with 2500 scans
    filtered_spectra = FilteredMassSpecData(
        spectra,
        max_scans = params.initial_scan_count,  # 2500
        topn_peaks = 1000,
        target_ms_order = UInt8(2)
    )
    
    @info "Initial sampling: $(length(filtered_spectra)) MS2 scans from $(count_ms2_scans(spectra)) total"
    
    # Main convergence loop
    while n_attempts < 5
        # ... set models ...
        
        # Collect PSMs
        psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
        final_psm_count = size(psms, 1)
        
        @info "Iteration $(n_attempts + 1): Collected $final_psm_count PSMs from $(length(filtered_spectra)) scans"
        
        # Check if we have enough PSMs
        if final_psm_count < min_psms_for_fitting
            @warn "Insufficient PSMs: $final_psm_count < $min_psms_for_fitting"
            
            # First retry: Expand to 5000 scans
            if n_attempts == 0 && length(filtered_spectra) < params.expanded_scan_count
                additional_scans = params.expanded_scan_count - length(filtered_spectra)
                n_added = append!(filtered_spectra, max_additional_scans = additional_scans)
                @info "Expanded sampling: added $n_added scans (now $(length(filtered_spectra)) total)"
                n_attempts += 1
                continue
            end
            
            # Second retry: Increase mass tolerance
            if n_attempts == 1
                current_tol = getLeftTol(results.mass_err_model[])
                new_tol = min(current_tol * 1.5f0, 50.0f0)
                results.mass_err_model[] = MassErrorModel(
                    getMassOffset(results.mass_err_model[]),
                    (new_tol, new_tol)
                )
                @info "Expanding mass tolerance from $current_tol to $new_tol ppm"
                n_attempts += 1
                continue
            end
            
            # Give up after 2 attempts
            break
        end
        
        # ... rest of convergence logic ...
    end
end
```

### 4. Update FilteredMassSpecData Implementation

**Key changes needed:**
```julia
mutable struct FilteredMassSpecData{T<:AbstractFloat} <: MassSpecData
    # ... existing fields ...
    total_ms2_scans::Int  # Track total MS2 scans in original data
    max_scans_sampled::Int  # Track how many we tried to sample
end

function FilteredMassSpecData(
    spectra::MassSpecData,
    max_scans::Int = 2500,
    topn_peaks::Union{Nothing, Int} = nothing;
    target_ms_order::UInt8 = UInt8(2)
)
    # Count total MS2 scans
    ms2_indices = findall(i -> getMsOrder(spectra, i) == target_ms_order, 1:length(spectra))
    total_ms2 = length(ms2_indices)
    
    # Randomly sample up to max_scans
    n_to_sample = min(max_scans, total_ms2)
    sampled_indices = sample(ms2_indices, n_to_sample, replace=false)
    sort!(sampled_indices)  # Maintain scan order
    
    # ... rest of construction ...
    
    return FilteredMassSpecData(
        # ... fields ...
        total_ms2_scans = total_ms2,
        max_scans_sampled = n_to_sample
    )
end

# Add helper to get unsampled scan count
function getUnsampledCount(fms::FilteredMassSpecData)
    return fms.total_ms2_scans - length(fms)
end

# Update append! to work with scan counts
function append!(
    fms::FilteredMassSpecData,
    max_additional_scans::Int = 2500
)
    # Sample from remaining unsampled scans
    unsampled = setdiff(1:fms.total_ms2_scans, fms.original_scan_indices)
    n_to_add = min(max_additional_scans, length(unsampled))
    
    if n_to_add == 0
        return 0
    end
    
    new_indices = sample(unsampled, n_to_add, replace=false)
    # ... append logic ...
    
    return n_to_add
end
```

### 5. Update collect_psms Function

**Changes needed:**
- Remove internal sampling logic (now handled by FilteredMassSpecData)
- Remove the iteration loop for adding more scans (now handled in process_file!)
- Simplify to just perform search and scoring

```julia
function collect_psms(
    filtered_spectra::FilteredMassSpecData,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:ParameterTuningSearchParameters}
    
    # Single search on provided filtered_spectra
    psms = library_search(filtered_spectra, search_context, params, ms_file_idx)
    
    if !iszero(size(psms, 1))
        # Map indices and add columns
        psms[!, :filtered_scan_idx] = psms[!, :scan_idx]
        psms[!, :scan_idx] = [
            getOriginalScanIndex(filtered_spectra, idx) 
            for idx in psms[!, :filtered_scan_idx]
        ]
        
        add_tuning_search_columns!(psms, spectra, 
            getPrecursors(getSpecLib(search_context)), params)
        
        # Score and filter
        filter_and_score_psms!(psms, params, search_context)
    end
    
    return psms
end
```

## Implementation Steps

1. **Update FilteredMassSpecData.jl**
   - Modify constructor to use max_scans
   - Add scan counting utilities
   - Update append! for scan-based sampling

2. **Update types.jl**
   - Add initial_scan_count and expanded_scan_count parameters
   - Update parameter extraction

3. **Simplify collect_psms**
   - Remove iteration loop
   - Single-pass search only

4. **Update process_file!**
   - Add scan count logging
   - Implement two-stage retry strategy:
     - First: Expand scans to 5000
     - Second: Increase mass tolerance
   - Better progress reporting

5. **Update parameter files**
   - Add new scan count parameters to default configs
   - Maintain backward compatibility if possible

## Benefits

1. **Transparency**: Clear visibility into scan counts and PSM yields
2. **Control**: Explicit control over computational cost
3. **Efficiency**: Better retry strategy with targeted expansions
4. **Debugging**: Easier to diagnose low PSM issues
5. **Scalability**: Can adjust defaults based on dataset size

## Testing Strategy

1. Test with standard E. coli dataset
2. Test with artificially reduced data (simulate low PSM scenario)
3. Verify scan expansion works correctly
4. Ensure backward compatibility with existing parameters
5. Check memory usage with 5000 scan sampling

## Migration Notes

- The sample_rate parameter should be deprecated but kept for backward compatibility
- If sample_rate is provided, convert to scan count: `max_scans = Int(round(total_ms2_scans * sample_rate))`
- Default values should be tuned based on typical DIA experiments

## Expected Behavior

### Scenario 1: Good Data
- Sample 2500 scans initially
- Get sufficient PSMs (>1000)
- Proceed with parameter fitting
- Total time: ~15-20 seconds

### Scenario 2: Challenging Data
- Sample 2500 scans initially
- Get insufficient PSMs (<1000)
- Expand to 5000 scans
- Retry and get sufficient PSMs
- Total time: ~25-30 seconds

### Scenario 3: Very Difficult Data
- Sample 2500 scans initially
- Get insufficient PSMs
- Expand to 5000 scans
- Still insufficient PSMs
- Increase mass tolerance by 1.5x
- Either succeed or fail gracefully
- Total time: ~35-40 seconds

## Configuration Example

```json
{
  "parameter_tuning": {
    "initial_scan_count": 2500,
    "expanded_scan_count": 5000,
    "min_psms": 1000,
    "frag_tol_ppm": 20.0,
    "topn_peaks": 1000
  }
}
```