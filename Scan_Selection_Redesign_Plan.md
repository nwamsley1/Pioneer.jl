# Scan Selection Redesign Plan

## Overview

Replace the current `FilteredMassSpecData` approach with a memory-efficient scan indexing strategy that:
- Never loads scan data into memory until needed
- Maintains intelligent RT-based prioritization
- Uses progressive sampling with exponential expansion
- Provides detailed timing and progress information

## Core Design Principles

1. **Memory Efficiency**: Work with scan indices only, no data preloading
2. **RT Coverage**: Ensure even sampling across retention time range
3. **Quality Priority**: Select highest TIC scans within each RT bin
4. **Progressive Sampling**: Start small, expand geometrically if needed
5. **Observability**: Detailed logging at each step with timing information

## Implementation Architecture

### Phase 1: Scan Index Building

```julia
struct ScanPriorityIndex
    scan_indices::Vector{Int32}        # Ordered scan indices
    rt_bin_boundaries::Vector{Float32} # RT bin edges
    tic_values::Vector{Float32}        # TIC for priority sorting
    rt_values::Vector{Float32}         # RT for binning
    ms_orders::Vector{UInt8}           # MS order filtering
    total_ms2_count::Int32              # Total MS2 scans available
    n_rt_bins::Int32                    # Number of RT bins
end
```

### Phase 2: Priority Ordering Algorithm

```
1. Extract metadata (RT, TIC, MS order) for all scans
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

### Phase 3: Progressive Sampling Strategy

```
Initial sample: a% of priority vector (default a=5%)
If convergence fails:
    Sample next 2*a% (total 3*a%)
If convergence fails:
    Sample next 4*a% (total 7*a%)
Continue: 8*a%, 16*a%, 32*a%, 64*a%, 100%
```

## Detailed Implementation Steps

### Step 1: Build Scan Priority Index

```julia
function build_scan_priority_index(
    spectra::MassSpecData;
    n_rt_bins::Int = 15,
    target_ms_order::UInt8 = UInt8(2),
    verbose::Bool = true
)::ScanPriorityIndex

    start_time = time()

    # Step 1.1: Extract metadata vectors (no peak data)
    if verbose
        println("┌─ Building scan priority index...")
        println("│  Extracting metadata from $(length(spectra)) scans...")
    end

    metadata_start = time()
    rt_values = getRetentionTimes(spectra)
    tic_values = getTICs(spectra)
    ms_orders = getMsOrders(spectra)
    metadata_time = time() - metadata_start

    if verbose
        println("│  ✓ Metadata extraction: $(round(metadata_time, digits=3))s")
    end

    # Step 1.2: Filter for target MS order
    filter_start = time()
    ms2_mask = ms_orders .== target_ms_order
    ms2_indices = findall(ms2_mask)
    n_ms2 = length(ms2_indices)
    filter_time = time() - filter_start

    if verbose
        println("│  ✓ MS2 filtering: $(round(filter_time, digits=3))s")
        println("│    Found $n_ms2 MS2 scans out of $(length(spectra)) total")
    end

    # Step 1.3: Create RT bins
    binning_start = time()
    rt_min, rt_max = extrema(rt_values[ms2_indices])
    bin_width = (rt_max - rt_min) / n_rt_bins
    rt_bin_boundaries = collect(rt_min:bin_width:rt_max)
    binning_time = time() - binning_start

    if verbose
        println("│  ✓ RT binning setup: $(round(binning_time, digits=3))s")
        println("│    RT range: $rt_min - $rt_max minutes")
        println("│    Bin width: $(round(bin_width, digits=2)) minutes")
    end

    # Step 1.4: Assign scans to bins and sort by TIC
    sorting_start = time()
    scan_priority_order = create_priority_ordering(
        ms2_indices, rt_values, tic_values,
        rt_bin_boundaries, n_rt_bins, verbose
    )
    sorting_time = time() - sorting_start

    if verbose
        println("│  ✓ Priority ordering: $(round(sorting_time, digits=3))s")
    end

    total_time = time() - start_time
    if verbose
        println("└─ Index building complete: $(round(total_time, digits=3))s total")
    end

    return ScanPriorityIndex(
        scan_priority_order,
        rt_bin_boundaries,
        tic_values[ms2_indices],
        rt_values[ms2_indices],
        ms_orders[ms2_indices],
        Int32(n_ms2),
        Int32(n_rt_bins)
    )
end
```

### Step 2: Progressive Sampling Implementation

```julia
function progressive_psm_collection!(
    results::ParameterTuningSearchResults,
    scan_index::ScanPriorityIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters;
    initial_sample_percent::Float64 = 5.0,
    verbose::Bool = true
)

    total_scans = length(scan_index.scan_indices)
    scans_processed = 0
    sample_iteration = 0
    converged = false

    # Calculate sampling schedule
    sample_percentages = [initial_sample_percent]
    current_percent = initial_sample_percent
    while sum(sample_percentages) < 100.0
        next_percent = min(current_percent * 2, 100.0 - sum(sample_percentages))
        push!(sample_percentages, next_percent)
        current_percent = next_percent
    end

    if verbose
        println("\n┌─ Progressive PSM Collection")
        println("│  Total MS2 scans available: $total_scans")
        println("│  Sampling schedule: $(sample_percentages)%")
        println("│  Min PSMs for convergence: $(getMinPsms(params))")
    end

    # Progressive sampling loop
    for (iteration, sample_percent) in enumerate(sample_percentages)
        sample_iteration = iteration

        # Calculate scan range for this iteration
        start_idx = scans_processed + 1
        n_to_sample = ceil(Int, total_scans * sample_percent / 100)
        end_idx = min(scans_processed + n_to_sample, total_scans)
        scans_to_process = scan_index.scan_indices[start_idx:end_idx]

        if verbose
            println("│")
            println("│  ┌─ Iteration $iteration: Sampling $(length(scans_to_process)) scans")
            println("│  │  Scan range: $start_idx - $end_idx")
            println("│  │  Cumulative: $(end_idx) / $total_scans scans")
        end

        # Collect PSMs for this batch
        batch_start = time()
        psms = collect_psms_for_scans(
            scans_to_process,
            spectra,
            search_context,
            params,
            verbose
        )
        batch_time = time() - batch_start

        n_psms = size(psms, 1)
        if verbose
            println("│  │  ✓ PSM collection: $(round(batch_time, digits=3))s")
            println("│  │  Found $n_psms PSMs")
        end

        # Check convergence
        if n_psms >= getMinPsms(params)
            converged = true
            if verbose
                println("│  └─ ✓ CONVERGED with $n_psms PSMs")
            end
            break
        else
            if verbose
                println("│  └─ ✗ Not converged ($n_psms < $(getMinPsms(params)))")
            end
        end

        scans_processed = end_idx

        # Early exit if we've processed all scans
        if scans_processed >= total_scans
            if verbose
                println("│  └─ ⚠ All scans processed without convergence")
            end
            break
        end
    end

    if verbose
        println("└─ Progressive sampling complete")
        println("   Iterations: $sample_iteration")
        println("   Scans used: $scans_processed / $total_scans")
        println("   Converged: $converged")
    end

    return converged, scans_processed
end
```

### Step 3: PSM Collection for Scan Batch

```julia
function collect_psms_for_scans(
    scan_indices::Vector{Int32},
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    verbose::Bool = false
)::DataFrame

    if verbose
        println("│  │  ├─ Processing $(length(scan_indices)) scans...")
    end

    # Perform library search only on specified scans
    # This is where we actually access the scan data for the first time
    psms = library_search_indexed(
        spectra,
        scan_indices,  # Only these scans
        search_context,
        params
    )

    # Add metadata and score
    if !isempty(psms)
        add_tuning_search_columns!(psms, spectra, search_context)
        filter_and_score_psms!(psms, params, search_context)
    end

    return psms
end
```

## Integration Points

### Modified ParameterTuningSearch Flow

```julia
function process_file!(
    results::ParameterTuningSearchResults,
    params::ParameterTuningSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    verbose = true  # Enable detailed logging

    # Build scan priority index (no data loading)
    println("\n═══ File $ms_file_idx Processing ═══")
    scan_index = build_scan_priority_index(spectra; verbose=verbose)

    # Initialize models
    initialize_models!(search_context, ms_file_idx, params)

    # Run 3-phase search with progressive sampling
    for phase in 1:3
        println("\n─── Phase $phase ───")

        # Apply phase bias
        apply_phase_bias!(search_context, ms_file_idx, phase, params)

        # Try each score threshold
        for min_score in getMinIndexSearchScores(params)
            println("  Testing min_score = $min_score")

            # Progressive sampling
            converged, scans_used = progressive_psm_collection!(
                results, scan_index, spectra,
                search_context, params;
                verbose=verbose
            )

            if converged
                println("✓ CONVERGED in Phase $phase with score $min_score")
                return results
            end
        end
    end

    println("✗ Failed to converge after all phases")
    apply_fallback_strategy!(results, search_context, ms_file_idx)
    return results
end
```

## Performance Monitoring

### Timing Checkpoints

```julia
struct TimingStats
    index_building::Float64
    metadata_extraction::Float64
    rt_binning::Float64
    priority_ordering::Float64
    psm_collection::Vector{Float64}  # Per iteration
    total_time::Float64
end

function log_timing_stats(stats::TimingStats, file_path::String)
    open(file_path, "a") do io
        println(io, "Timing Statistics:")
        println(io, "  Index Building: $(stats.index_building)s")
        println(io, "    - Metadata: $(stats.metadata_extraction)s")
        println(io, "    - RT Binning: $(stats.rt_binning)s")
        println(io, "    - Priority Ordering: $(stats.priority_ordering)s")
        println(io, "  PSM Collection:")
        for (i, time) in enumerate(stats.psm_collection)
            println(io, "    - Iteration $i: $(time)s")
        end
        println(io, "  Total Time: $(stats.total_time)s")
    end
end
```

### Memory Usage Tracking

```julia
function log_memory_usage(label::String)
    gc_stats = Base.gc_num()
    mem_used = Sys.total_memory() - Sys.free_memory()
    println("Memory [$label]: $(round(mem_used / 1e9, digits=2)) GB")
    println("  GC collections: $(gc_stats.collect)")
    println("  GC time: $(round(gc_stats.total_time / 1e9, digits=3))s")
end
```

## Configuration

### New Parameters to Add

```json
{
  "parameter_tuning": {
    "progressive_sampling": {
      "initial_sample_percent": 5.0,
      "expansion_factor": 2.0,
      "n_rt_bins": 15,
      "verbose_logging": true,
      "timing_log_path": "timing_stats.log"
    }
  }
}
```

## Benefits Over Current Approaches

### vs FilteredMassSpecData
- **Memory**: O(n) metadata only vs O(n*m) with peak data
- **Flexibility**: Can adjust sampling without rebuilding
- **Speed**: No upfront data copying

### vs July 2024 Random Sampling
- **Coverage**: Guaranteed RT coverage
- **Quality**: Prioritizes high-TIC scans
- **Efficiency**: Progressive expansion vs fixed rate

## Implementation Phases

### Phase 1: Core Implementation
1. Implement `ScanPriorityIndex` struct
2. Create `build_scan_priority_index` function
3. Add verbose logging throughout

### Phase 2: Progressive Sampling
1. Implement `progressive_psm_collection!`
2. Integrate with existing PSM collection
3. Add timing statistics

### Phase 3: Integration
1. Modify `process_file!` to use new approach
2. Update configuration parameters
3. Remove `FilteredMassSpecData` dependency

### Phase 4: Optimization
1. Profile and optimize bottlenecks
2. Add parallel processing where beneficial
3. Fine-tune default parameters

## Testing Strategy

### Unit Tests
- Test scan priority ordering
- Verify RT bin assignment
- Check progressive sampling percentages

### Integration Tests
- Compare PSM yields with current method
- Verify memory usage improvements
- Benchmark timing vs FilteredMassSpecData

### Validation
- Run on test dataset
- Compare convergence rates
- Verify parameter estimation quality

## Questions for Implementation

1. **RT Binning**: Should bins be equal width or equal scan count?
2. **TIC vs Peak Count**: Which metric better indicates scan quality?
3. **Expansion Factor**: Is 2x geometric growth optimal, or should we use 1.5x or 3x?
4. **Minimum Sample**: Should we enforce a minimum number of scans per iteration?
5. **Parallelization**: Should scan indexing be parallelized across threads?
6. **Caching**: Should we cache the scan index for reuse across phases?

## Risk Mitigation

### Potential Issues
1. **Sparse RT regions**: Some bins may have few scans
   - Solution: Merge adjacent sparse bins

2. **Memory spikes**: PSM collection may create temporary peaks
   - Solution: Process in smaller chunks with GC hints

3. **Slow convergence**: May need many iterations
   - Solution: Adaptive expansion factor based on PSM yield

## Success Metrics

- Memory usage reduced by >50%
- Convergence rate similar or better
- Total runtime within 10% of current
- Clear logging of all operations
- No loss in parameter estimation quality