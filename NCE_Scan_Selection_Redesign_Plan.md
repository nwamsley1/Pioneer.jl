# NCE Tuning Search - Scan Selection Redesign Plan

## Overview

Replace the current `FilteredMassSpecData` approach in NCE tuning search with a memory-efficient scan indexing strategy that:
- Never loads scan data into memory until needed
- Maintains intelligent RT-based prioritization
- Uses progressive sampling **without replacement**
- Provides detailed timing and progress information

## Key Concept: Sampling Without Replacement

Given a prioritized vector of scan indices:
- First iteration: Sample indices [1 to 10] (10% of 100 scans)
- Second iteration: Sample indices [11 to 30] (next 20% of 100 scans)
- Third iteration: Sample indices [31 to 70] (next 40% of 100 scans)
- Continue until convergence or all scans used

**Critical**: We never re-process scans from previous iterations.

## Implementation Architecture

### Phase 1: Scan Index Building

```julia
struct NCEScanPriorityIndex
    scan_indices::Vector{Int32}        # Ordered scan indices by priority
    rt_values::Vector{Float32}         # RT values for each scan
    tic_values::Vector{Float32}        # TIC values for each scan
    ms_orders::Vector{UInt8}           # MS order for filtering
    total_ms2_count::Int32              # Total MS2 scans
    n_rt_bins::Int32                    # Number of RT bins used
end
```

### Phase 2: Priority Ordering Algorithm

```
1. Extract metadata (RT, TIC, MS order) for all scans - NO peak data
2. Filter for MS2 scans only
3. Divide RT range into N bins (default: 15)
4. Within each bin, sort scans by TIC (descending)
5. Create priority vector by round-robin from bins:
   - Take highest TIC scan from bin 1
   - Take highest TIC scan from bin 2
   - ...
   - Take highest TIC scan from bin N
   - Take second highest from bin 1
   - Continue until all scans ordered
```

### Phase 3: Progressive Sampling Without Replacement

```
Initial sample: First a% of priority vector (e.g., scans 1-10 if a=10%)
If more PSMs needed:
    Sample NEXT 2*a% (e.g., scans 11-30)
If more PSMs needed:
    Sample NEXT 4*a% (e.g., scans 31-70)
Continue until convergence or exhaustion
```

## Detailed Implementation

### Step 1: Build NCE Scan Priority Index

```julia
function build_nce_scan_priority_index(
    spectra::MassSpecData;
    n_rt_bins::Int = 15,
    target_ms_order::UInt8 = UInt8(2),
    verbose::Bool = true
)::NCEScanPriorityIndex

    start_time = time()

    if verbose
        println("┌─ Building NCE scan priority index...")
        println("│  Total scans in file: $(length(spectra))")
    end

    # Step 1.1: Extract metadata (NO peak data loaded)
    metadata_start = time()
    rt_values = getRetentionTimes(spectra)
    tic_values = getTICs(spectra)
    ms_orders = getMsOrders(spectra)
    metadata_time = time() - metadata_start

    if verbose
        println("│  ✓ Metadata extraction: $(round(metadata_time, digits=3))s")
    end

    # Step 1.2: Filter for MS2 scans
    filter_start = time()
    ms2_mask = ms_orders .== target_ms_order
    ms2_indices = findall(ms2_mask)
    n_ms2 = length(ms2_indices)
    filter_time = time() - filter_start

    if verbose
        println("│  ✓ MS2 filtering: $(round(filter_time, digits=3))s")
        println("│    Found $n_ms2 MS2 scans")
    end

    # Step 1.3: Create RT bins
    binning_start = time()
    ms2_rt = rt_values[ms2_indices]
    ms2_tic = tic_values[ms2_indices]

    rt_min, rt_max = extrema(ms2_rt)
    bin_width = (rt_max - rt_min) / n_rt_bins

    # Assign each MS2 scan to a bin
    bin_assignments = zeros(Int, n_ms2)
    for (i, rt) in enumerate(ms2_rt)
        bin_idx = min(max(1, ceil(Int, (rt - rt_min) / bin_width)), n_rt_bins)
        bin_assignments[i] = bin_idx
    end
    binning_time = time() - binning_start

    if verbose
        println("│  ✓ RT binning: $(round(binning_time, digits=3))s")
        println("│    RT range: $(round(rt_min, digits=2)) - $(round(rt_max, digits=2)) min")
        println("│    Bins: $n_rt_bins × $(round(bin_width, digits=2)) min")
    end

    # Step 1.4: Sort within bins by TIC and create priority order
    sorting_start = time()

    # Group scans by bin
    bins = [Int32[] for _ in 1:n_rt_bins]
    for (i, bin_idx) in enumerate(bin_assignments)
        push!(bins[bin_idx], ms2_indices[i])
    end

    # Sort each bin by TIC (descending)
    for bin in bins
        if !isempty(bin)
            sort!(bin, by=idx -> tic_values[idx], rev=true)
        end
    end

    # Round-robin selection from bins
    priority_order = Int32[]
    max_bin_size = maximum(length.(bins))

    for round in 1:max_bin_size
        for bin_idx in 1:n_rt_bins
            if round <= length(bins[bin_idx])
                push!(priority_order, bins[bin_idx][round])
            end
        end
    end

    sorting_time = time() - sorting_start

    if verbose
        println("│  ✓ Priority ordering: $(round(sorting_time, digits=3))s")

        # Show distribution across bins
        bin_counts = [length(bin) for bin in bins]
        println("│    Scans per bin: min=$(minimum(bin_counts)), max=$(maximum(bin_counts)), mean=$(round(mean(bin_counts), digits=1))")
    end

    total_time = time() - start_time
    if verbose
        println("└─ Index building complete: $(round(total_time, digits=3))s total")
        println("   Priority vector contains $n_ms2 scans")
    end

    return NCEScanPriorityIndex(
        priority_order,
        rt_values[priority_order],
        tic_values[priority_order],
        ms_orders[priority_order],
        Int32(n_ms2),
        Int32(n_rt_bins)
    )
end
```

### Step 2: Progressive Sampling Without Replacement

```julia
function progressive_nce_psm_collection!(
    scan_index::NCEScanPriorityIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NCETuningSearchParameters,
    ms_file_idx::Int64;
    initial_percent::Float64 = 10.0,
    min_psms_required::Int = 100,
    verbose::Bool = true
)
    total_scans = length(scan_index.scan_indices)
    scans_processed = 0  # Track position in priority vector
    iteration = 0
    converged = false
    all_psms = DataFrame()  # Accumulate PSMs across iterations

    # Calculate sampling schedule (percentages of total)
    sample_schedule = Float64[]
    remaining = 100.0
    current_chunk = initial_percent
    while remaining > 0
        chunk = min(current_chunk, remaining)
        push!(sample_schedule, chunk)
        remaining -= chunk
        current_chunk *= 2  # Double the chunk size each time
    end

    if verbose
        println("\n┌─ NCE Progressive PSM Collection")
        println("│  Total MS2 scans: $total_scans")
        println("│  Target PSMs: $min_psms_required")
        println("│  Sampling schedule: $(sample_schedule)%")
        println("│")
    end

    # Progressive sampling loop
    for chunk_percent in sample_schedule
        iteration += 1

        # Calculate range for THIS iteration (no overlap with previous)
        n_scans_this_chunk = ceil(Int, total_scans * chunk_percent / 100)
        start_idx = scans_processed + 1
        end_idx = min(scans_processed + n_scans_this_chunk, total_scans)

        if start_idx > total_scans
            if verbose
                println("│  ⚠ All scans exhausted")
            end
            break
        end

        # Extract scan indices for this chunk
        chunk_scan_indices = scan_index.scan_indices[start_idx:end_idx]
        n_chunk_scans = length(chunk_scan_indices)

        if verbose
            println("│  ┌─ Iteration $iteration")
            println("│  │  Processing scans $start_idx to $end_idx ($n_chunk_scans scans)")
            println("│  │  Cumulative: $end_idx / $total_scans ($(round(100*end_idx/total_scans, digits=1))%)")
        end

        # Collect PSMs for this chunk
        chunk_start = time()
        chunk_psms = collect_nce_psms_for_scans(
            chunk_scan_indices,
            spectra,
            search_context,
            params,
            ms_file_idx,
            verbose
        )
        chunk_time = time() - chunk_start

        # Combine with previous PSMs
        if !isempty(chunk_psms)
            all_psms = isempty(all_psms) ? chunk_psms : vcat(all_psms, chunk_psms)
        end

        total_psms = nrow(all_psms)

        if verbose
            println("│  │  Time: $(round(chunk_time, digits=3))s")
            println("│  │  PSMs this chunk: $(nrow(chunk_psms))")
            println("│  │  Total PSMs: $total_psms")
        end

        # Check convergence
        if total_psms >= min_psms_required
            converged = true
            if verbose
                println("│  └─ ✓ CONVERGED with $total_psms PSMs")
            end
            break
        else
            if verbose
                println("│  └─ ✗ Need more PSMs ($total_psms < $min_psms_required)")
            end
        end

        # Update position in priority vector
        scans_processed = end_idx
    end

    if verbose
        println("└─ Collection complete")
        println("   Iterations: $iteration")
        println("   Scans used: $scans_processed / $total_scans ($(round(100*scans_processed/total_scans, digits=1))%)")
        println("   PSMs found: $(nrow(all_psms))")
        println("   Converged: $converged")
    end

    return all_psms, converged, scans_processed
end
```

### Step 3: PSM Collection for Specific Scans

```julia
function collect_nce_psms_for_scans(
    scan_indices::Vector{Int32},
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NCETuningSearchParameters,
    ms_file_idx::Int64,
    verbose::Bool
)::DataFrame

    if verbose
        println("│  │  ├─ Loading and processing $(length(scan_indices)) scans...")
    end

    # This is where we FIRST access actual scan data
    # Library search only processes the specified scan indices
    psms = library_search_with_indices(
        spectra,
        scan_indices,  # Only these specific scans
        search_context,
        params,
        ms_file_idx
    )

    # Score and filter PSMs
    if !isempty(psms)
        # Add necessary columns
        add_nce_search_columns!(psms, spectra, search_context)

        # Score PSMs
        score_nce_psms!(psms, params)

        # Apply FDR filtering
        filter!(row -> row.q_value <= params.max_q_value, psms)
    end

    return psms
end
```

### Step 4: Integration with NCE Tuning Search

```julia
function process_file!(
    results::NCETuningSearchResults,
    params::NCETuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    verbose = params.verbose_logging

    if verbose
        println("\n" * "═"^60)
        println(" NCE Tuning Search - File $ms_file_idx")
        println("═"^60)
    end

    # Build scan priority index (metadata only, no peak data)
    index_start = time()
    scan_index = build_nce_scan_priority_index(
        spectra;
        n_rt_bins = params.n_rt_bins,
        verbose = verbose
    )
    index_time = time() - index_start

    # Log memory usage after index building
    if verbose
        log_memory_usage("After index building")
    end

    # Progressive PSM collection with sampling without replacement
    collection_start = time()
    all_psms, converged, scans_used = progressive_nce_psm_collection!(
        scan_index,
        spectra,
        search_context,
        params,
        ms_file_idx;
        initial_percent = params.initial_sample_percent,
        min_psms_required = params.min_psms_for_nce,
        verbose = verbose
    )
    collection_time = time() - collection_start

    # Fit NCE model from collected PSMs
    if converged && !isempty(all_psms)
        model_start = time()
        nce_model = fit_nce_model(all_psms, params)
        store_nce_model!(search_context, ms_file_idx, nce_model)
        model_time = time() - model_start

        if verbose
            println("\n✓ NCE model fitted successfully")
            println("  Model fitting time: $(round(model_time, digits=3))s")
        end
    else
        if verbose
            println("\n✗ Failed to collect sufficient PSMs for NCE modeling")
            println("  Using default NCE model")
        end
        # Use default or borrowed model
        apply_fallback_nce_model!(search_context, ms_file_idx)
    end

    # Log final statistics
    if verbose
        total_time = index_time + collection_time
        println("\n" * "─"^40)
        println("Timing Summary:")
        println("  Index building: $(round(index_time, digits=3))s")
        println("  PSM collection: $(round(collection_time, digits=3))s")
        println("  Total time: $(round(total_time, digits=3))s")
        println("  Scans used: $scans_used / $(length(scan_index.scan_indices))")
        println("  Memory efficiency: No scan data preloaded")
        log_memory_usage("Final")
    end

    return results
end
```

## Memory Management

### Key Principle: Lazy Loading

```julia
# NEVER do this:
all_scan_data = [getMzArray(spectra, i) for i in 1:length(spectra)]

# ALWAYS do this:
for scan_idx in selected_scans
    scan_data = getMzArray(spectra, scan_idx)  # Load only when needed
    process(scan_data)
    # scan_data goes out of scope, eligible for GC
end
```

### Memory Monitoring

```julia
function log_memory_usage(checkpoint::String)
    GC.gc()  # Force collection for accurate reading
    mem_used = Base.gc_live_bytes() / 1e9
    println("Memory [$checkpoint]: $(round(mem_used, digits=2)) GB")
end
```

## Configuration Parameters

### Add to NCE Tuning Search Config

```json
{
  "nce_tuning": {
    "progressive_sampling": {
      "initial_sample_percent": 10.0,
      "min_psms_for_nce": 100,
      "n_rt_bins": 15,
      "verbose_logging": true
    }
  }
}
```

## Key Differences from Current Implementation

### Memory Usage
- **Current**: FilteredMassSpecData loads all selected scan data upfront
- **New**: Only metadata loaded initially, scan data loaded on-demand

### Sampling Strategy
- **Current**: All scans loaded, then processed
- **New**: Progressive sampling without replacement

### Example Memory Comparison

For 10,000 MS2 scans with 1000 peaks each:
- **Current**: ~400 MB upfront (all scan data)
- **New**: ~400 KB upfront (metadata only)
- **New during processing**: ~40 MB per 1000 scans processed

## Implementation Timeline

### Phase 1: Core Implementation (Day 1)
1. Create `NCEScanPriorityIndex` struct
2. Implement `build_nce_scan_priority_index`
3. Add comprehensive logging

### Phase 2: Progressive Sampling (Day 2)
1. Implement `progressive_nce_psm_collection!`
2. Ensure sampling without replacement
3. Add memory monitoring

### Phase 3: Integration (Day 3)
1. Modify NCE tuning search to use new approach
2. Update library search interface
3. Test with real data

### Phase 4: Optimization (Day 4)
1. Profile and identify bottlenecks
2. Fine-tune default parameters
3. Add unit tests

## Success Metrics

- ✓ Zero scan data loaded until needed
- ✓ Memory usage reduced by >90% during index building
- ✓ Clear logging of every operation with timing
- ✓ Sampling without replacement working correctly
- ✓ NCE model quality maintained or improved

## Testing Checklist

- [ ] Verify no scan data loaded during index building
- [ ] Confirm sampling without replacement
- [ ] Check memory usage at each checkpoint
- [ ] Validate NCE model quality
- [ ] Test with files of varying sizes
- [ ] Verify convergence behavior
- [ ] Check timing for each component