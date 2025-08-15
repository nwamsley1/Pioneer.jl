# Intelligent Scan Selection for FilteredMassSpecData

## Overview

Replace random scan selection with intelligent prioritization based on:
1. **Retention time coverage**: Divide RT space into equal regions
2. **Peak density**: Prioritize scans with more peaks
3. **Interleaved sampling**: Ensure coverage across all RT regions

## Data Structure Changes

### Updated FilteredMassSpecData Type

```julia
struct FilteredMassSpecData
    # Existing fields
    original_spectra::MassSpecData
    filtered_indices::Vector{Int32}  # Indices of scans included
    topn::Int64
    target_ms_order::UInt8
    
    # New fields for intelligent selection
    scan_priority_order::Vector{Int32}    # Pre-computed scan ordering
    n_scans_sampled::Int32                # How many scans from priority order are included
    rt_bin_assignments::Vector{Int8}      # Which RT bin each scan belongs to
end
```

## Algorithm for Computing Scan Priority Order

### Step 1: RT Binning

```julia
function compute_rt_bins(spectra::MassSpecData, n_bins::Int = 15)
    # Get RT range
    rt_values = getRetentionTimes(spectra)
    rt_min, rt_max = extrema(rt_values)
    bin_width = (rt_max - rt_min) / n_bins
    
    # Assign each scan to a bin
    rt_bin_assignments = Vector{Int8}(undef, length(rt_values))
    for (idx, rt) in enumerate(rt_values)
        bin_idx = min(ceil(Int, (rt - rt_min) / bin_width), n_bins)
        rt_bin_assignments[idx] = bin_idx
    end
    
    return rt_bin_assignments, rt_min, rt_max, bin_width
end
```

### Step 2: Sort Scans Within Bins (Optimized)

```julia
function sort_scans_by_peak_density(spectra::MassSpecData, 
                                   target_ms_order::UInt8,
                                   rt_bin_assignments::Vector{Int8},
                                   n_bins::Int = 15)
    # Separate scans by MS order
    ms_orders = getMsOrders(spectra)
    target_scan_indices = findall(==(target_ms_order), ms_orders)
    n_target_scans = length(target_scan_indices)
    
    # Pre-allocate single array for all scans
    all_sorted_scans = Vector{Int32}(undef, n_target_scans)
    
    # Count scans per bin for pre-allocation
    bin_counts = zeros(Int32, n_bins)
    for scan_idx in target_scan_indices
        bin_idx = rt_bin_assignments[scan_idx]
        bin_counts[bin_idx] += 1
    end
    
    # Calculate bin boundaries in the output array
    bin_starts = Vector{Int32}(undef, n_bins)
    bin_ends = Vector{Int32}(undef, n_bins)
    cumsum = 0
    for i in 1:n_bins
        bin_starts[i] = cumsum + 1
        cumsum += bin_counts[i]
        bin_ends[i] = cumsum
    end
    
    # Reset bin_counts to use as write positions
    fill!(bin_counts, 0)
    
    # Write scans to their bin positions
    for scan_idx in target_scan_indices
        bin_idx = rt_bin_assignments[scan_idx]
        write_pos = bin_starts[bin_idx] + bin_counts[bin_idx]
        all_sorted_scans[write_pos] = scan_idx
        bin_counts[bin_idx] += 1
    end
    
    # In-place sort each bin by peak count (descending)
    peak_counts = getPeakCounts(spectra)  # Or compute from intensity arrays
    for i in 1:n_bins
        if bin_starts[i] <= bin_ends[i]
            # Sort this bin's slice in-place
            bin_slice = view(all_sorted_scans, bin_starts[i]:bin_ends[i])
            sort!(bin_slice, by=idx -> peak_counts[idx], rev=true)
        end
    end
    
    # Return the sorted array with bin boundaries for interleaving
    return all_sorted_scans, bin_starts, bin_ends
end
```

### Step 3: Interleave Bins for Priority Order (Optimized)

```julia
function create_priority_order(all_sorted_scans::Vector{Int32},
                              bin_starts::Vector{Int32},
                              bin_ends::Vector{Int32})
    n_bins = length(bin_starts)
    n_total_scans = length(all_sorted_scans)
    
    # Pre-allocate priority order array
    priority_order = Vector{Int32}(undef, n_total_scans)
    
    # Track current position in each bin
    bin_positions = copy(bin_starts)
    
    # Interleave: take one scan from each bin in round-robin fashion
    write_idx = 1
    scans_remaining = n_total_scans
    
    while scans_remaining > 0
        for bin_idx in 1:n_bins
            if bin_positions[bin_idx] <= bin_ends[bin_idx]
                # Write next scan from this bin
                priority_order[write_idx] = all_sorted_scans[bin_positions[bin_idx]]
                bin_positions[bin_idx] += 1
                write_idx += 1
                scans_remaining -= 1
                
                if scans_remaining == 0
                    break
                end
            end
        end
    end
    
    return priority_order
end
```

## Constructor Implementation

```julia
function FilteredMassSpecData(
    spectra::MassSpecData;
    max_scans::Int = 500,
    topn::Int = 200,
    target_ms_order::UInt8 = UInt8(2),
    n_rt_bins::Int = 15
)
    # Step 1: Compute RT bins
    rt_bin_assignments, rt_min, rt_max, bin_width = compute_rt_bins(spectra, n_rt_bins)
    
    # Step 2: Sort scans within bins by peak density (returns sorted array + boundaries)
    all_sorted_scans, bin_starts, bin_ends = sort_scans_by_peak_density(
        spectra, target_ms_order, rt_bin_assignments, n_rt_bins
    )
    
    # Step 3: Create interleaved priority order
    scan_priority_order = create_priority_order(all_sorted_scans, bin_starts, bin_ends)
    
    # Step 4: Sample initial scans following priority order
    n_to_sample = min(max_scans, length(scan_priority_order))
    filtered_indices = scan_priority_order[1:n_to_sample]
    
    # Step 5: Apply topN filtering to selected scans
    # ... existing topN logic ...
    
    return FilteredMassSpecData(
        spectra,
        filtered_indices,
        topn,
        target_ms_order,
        scan_priority_order,
        n_to_sample,
        rt_bin_assignments
    )
end
```

## Append Implementation

```julia
function Base.append!(
    filtered::FilteredMassSpecData;
    max_additional_scans::Int
)
    # Determine how many more scans we can add
    n_available = length(filtered.scan_priority_order) - filtered.n_scans_sampled
    n_to_add = min(max_additional_scans, n_available)
    
    if n_to_add <= 0
        @warn "No more scans available to add"
        return filtered
    end
    
    # Get next scans from priority order
    start_idx = filtered.n_scans_sampled + 1
    end_idx = filtered.n_scans_sampled + n_to_add
    new_scan_indices = filtered.scan_priority_order[start_idx:end_idx]
    
    # Add to filtered indices
    append!(filtered.filtered_indices, new_scan_indices)
    
    # Update counter
    filtered.n_scans_sampled += n_to_add
    
    # Apply topN filtering to new scans
    # ... existing topN logic ...
    
    @info "Added $n_to_add scans (now have $(length(filtered.filtered_indices)) total)"
    
    return filtered
end
```

## Example Usage

```julia
# Initial creation with 500 scans
filtered_spectra = FilteredMassSpecData(
    spectra,
    max_scans = 500,
    topn = 200,
    target_ms_order = UInt8(2)
)

# Scans are selected in this order:
# 1. Best scan from RT bin 1
# 2. Best scan from RT bin 2
# 3. Best scan from RT bin 3
# ... through bin 15
# 16. Second-best scan from RT bin 1
# 17. Second-best scan from RT bin 2
# ... and so on

# Later, append more scans following the same priority
append!(filtered_spectra, max_additional_scans = 500)
# This adds the next 500 scans from the priority order
```

## Benefits

1. **Better RT Coverage**: Ensures all retention time regions are represented early
2. **Quality First**: Within each region, highest-quality scans (most peaks) are selected first
3. **Deterministic**: Same input always produces same selection order
4. **Efficient Scaling**: When adding more scans, we continue from where we left off
5. **Memory Efficient**: Pre-computed order avoids repeated sorting

## Implementation Notes

1. **Peak Count Metric**: Can use either:
   - Number of peaks above intensity threshold
   - Total number of peaks
   - Sum of peak intensities (TIC)

2. **RT Bin Count**: 15 bins is suggested, but could be configurable:
   - More bins = finer RT resolution
   - Fewer bins = more scans per bin for better quality sorting

3. **MS Order Handling**: 
   - Separate priority orders for MS1 and MS2
   - Switch between them based on target_ms_order

4. **Edge Cases**:
   - Empty bins: Skip in interleaving
   - Few scans: All scans included regardless of priority
   - Uneven distribution: Some bins may exhaust before others

## Testing Strategy

```julia
# Test 1: Verify RT coverage
filtered = FilteredMassSpecData(spectra, max_scans=15)
# Should have exactly one scan from each RT bin

# Test 2: Verify peak density ordering
for bin_idx in 1:15
    bin_scans = filter(i -> rt_bin_assignments[i] == bin_idx, filtered.filtered_indices)
    peak_counts = [getPeakCount(spectra, idx) for idx in bin_scans]
    @test issorted(peak_counts, rev=true)
end

# Test 3: Verify append continuation
initial_count = filtered.n_scans_sampled
append!(filtered, max_additional_scans=100)
@test filtered.n_scans_sampled == initial_count + 100
```

## Questions for Approval

1. Should the number of RT bins (15) be configurable or fixed?
2. What metric should we use for "number of peaks"?
   - Raw peak count?
   - Peaks above intensity threshold?
   - Total ion current (TIC)?
3. Should we maintain separate priority orders for MS1 and MS2, or compute on-the-fly based on target_ms_order?
4. Should the priority order be recomputed if topN changes, or keep the original order?