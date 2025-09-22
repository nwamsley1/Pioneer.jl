# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Memory-efficient scan selection for quadrupole tuning search.
Uses intelligent isolation window m/z binning with TIC-based prioritization.
"""

"""
    QuadScanPriorityIndex

Contains prioritized scan indices for memory-efficient quadrupole tuning.
Only stores metadata, never loads actual scan data until needed.
"""
struct QuadScanPriorityIndex
    scan_indices::Vector{Int32}        # Ordered scan indices by priority
    center_mz_values::Vector{Float32}  # Center m/z values for each scan
    tic_values::Vector{Float32}        # TIC values for each scan
    ms_orders::Vector{UInt8}           # MS order for verification
    total_ms2_count::Int32             # Total MS2 scans
    n_mz_bins::Int32                   # Number of m/z bins used
end

"""
    build_quad_scan_priority_index(spectra::MassSpecData;
                                   n_mz_bins::Int = 20,
                                   target_ms_order::UInt8 = UInt8(2),
                                   verbose::Bool = true)::QuadScanPriorityIndex

Build prioritized scan index for quadrupole tuning without loading scan data.

# Process
1. Extract metadata (center m/z, TIC, MS order) for all scans - NO peak data
2. Filter for MS2 scans only
3. Divide center m/z range into N bins for even coverage across isolation windows
4. Within each bin, sort scans by TIC (descending)
5. Create priority vector by round-robin from bins

# Arguments
- `spectra`: MassSpecData object
- `n_mz_bins`: Number of m/z bins for even coverage (default: 20)
- `target_ms_order`: MS order to filter for (default: 2)
- `verbose`: Enable detailed logging with timing

# Returns
QuadScanPriorityIndex containing ordered scan indices and metadata

# Notes
- Binning by center m/z ensures good coverage across the isolation window range
- This is critical for quadrupole transmission modeling which needs data across m/z
- Higher TIC scans within each bin provide better isotope patterns for modeling
"""
function build_quad_scan_priority_index(
    spectra::MassSpecData;
    n_mz_bins::Int = 20,
    target_ms_order::UInt8 = UInt8(2),
    verbose::Bool = true
)::QuadScanPriorityIndex

    start_time = time()

    if verbose
        println("┌─ Building Quad scan priority index...")
        println("│  Total scans in file: $(length(spectra))")
    end

    # Step 1: Extract metadata (NO peak data loaded)
    metadata_start = time()
    if verbose
        println("│  ├─ Extracting metadata (center m/z, TIC, MS order)...")
    end

    center_mz_values = getCenterMzs(spectra)
    tic_values = getTICs(spectra)
    ms_orders = getMsOrders(spectra)

    metadata_time = time() - metadata_start
    if verbose
        println("│  │  ✓ Metadata extraction: $(round(metadata_time, digits=3))s")
    end

    # Step 2: Filter for MS2 scans
    filter_start = time()
    if verbose
        println("│  ├─ Filtering for MS2 scans...")
    end

    ms2_mask = ms_orders .== target_ms_order
    ms2_indices = findall(ms2_mask)
    n_ms2 = length(ms2_indices)

    filter_time = time() - filter_start
    if verbose
        println("│  │  ✓ MS2 filtering: $(round(filter_time, digits=3))s")
        println("│  │  Found $n_ms2 MS2 scans")
    end

    if n_ms2 == 0
        if verbose
            println("│  └─ ⚠ No MS2 scans found, returning empty index")
        end
        return QuadScanPriorityIndex(
            Int32[], Float32[], Float32[], UInt8[],
            Int32(0), Int32(n_mz_bins)
        )
    end

    # Step 3: Create center m/z bins
    binning_start = time()
    if verbose
        println("│  ├─ Creating center m/z bins...")
    end

    ms2_center_mz = center_mz_values[ms2_indices]
    ms2_tic = tic_values[ms2_indices]

    # Filter out missing values
    valid_mask = .!ismissing.(ms2_center_mz)
    if !all(valid_mask)
        if verbose
            println("│  │  ⚠ Filtering out $(count(.!valid_mask)) scans with missing center m/z")
        end
        ms2_indices = ms2_indices[valid_mask]
        ms2_center_mz = ms2_center_mz[valid_mask]
        ms2_tic = ms2_tic[valid_mask]
        n_ms2 = length(ms2_indices)
    end

    if n_ms2 == 0
        if verbose
            println("│  └─ ⚠ No valid MS2 scans found after filtering, returning empty index")
        end
        return QuadScanPriorityIndex(
            Int32[], Float32[], Float32[], UInt8[],
            Int32(0), Int32(n_mz_bins)
        )
    end

    mz_min, mz_max = extrema(ms2_center_mz)

    # Handle edge case of single m/z value
    if mz_min == mz_max
        if verbose
            println("│  │  ⚠ All scans have same center m/z ($(mz_min)), using single bin")
        end
        priority_order = Int32.(ms2_indices[sortperm(ms2_tic, rev=true)])
        binning_time = time() - binning_start
        total_time = time() - start_time

        if verbose
            println("│  │  ✓ Single-bin ordering: $(round(binning_time, digits=3))s")
            println("└─ Index building complete: $(round(total_time, digits=3))s total")
        end

        return QuadScanPriorityIndex(
            priority_order,
            center_mz_values[priority_order],
            tic_values[priority_order],
            ms_orders[priority_order],
            Int32(n_ms2),
            Int32(1)
        )
    end

    bin_width = (mz_max - mz_min) / n_mz_bins

    # Assign each MS2 scan to a bin
    bin_assignments = zeros(Int, n_ms2)
    for (i, center_mz) in enumerate(ms2_center_mz)
        bin_idx = min(max(1, ceil(Int, (center_mz - mz_min) / bin_width)), n_mz_bins)
        bin_assignments[i] = bin_idx
    end

    binning_time = time() - binning_start
    if verbose
        println("│  │  ✓ Center m/z binning: $(round(binning_time, digits=3))s")
        println("│  │  m/z range: $(round(mz_min, digits=2)) - $(round(mz_max, digits=2)) m/z")
        println("│  │  Bins: $n_mz_bins × $(round(bin_width, digits=2)) m/z")
    end

    # Step 4: Sort within bins by TIC and create priority order
    sorting_start = time()
    if verbose
        println("│  ├─ Sorting scans within bins by TIC...")
    end

    # Group scans by bin
    bins = [Int32[] for _ in 1:n_mz_bins]
    for (i, bin_idx) in enumerate(bin_assignments)
        push!(bins[bin_idx], ms2_indices[i])
    end

    # Sort each bin by TIC (descending - highest TIC first)
    for bin in bins
        if !isempty(bin)
            sort!(bin, by=idx -> tic_values[idx], rev=true)
        end
    end

    # Round-robin selection from bins
    priority_order = Int32[]
    max_bin_size = maximum(length.(bins))

    for round in 1:max_bin_size
        for bin_idx in 1:n_mz_bins
            if round <= length(bins[bin_idx])
                push!(priority_order, bins[bin_idx][round])
            end
        end
    end

    sorting_time = time() - sorting_start
    if verbose
        println("│  │  ✓ Priority ordering: $(round(sorting_time, digits=3))s")

        # Show distribution across bins
        bin_counts = [length(bin) for bin in bins]
        non_empty_bins = count(x -> x > 0, bin_counts)
        if non_empty_bins > 0
            println("│  │  Scans per bin: min=$(minimum(bin_counts[bin_counts .> 0])), max=$(maximum(bin_counts)), mean=$(round(mean(bin_counts), digits=1))")
            println("│  │  Non-empty bins: $non_empty_bins / $n_mz_bins")
        end
    end

    total_time = time() - start_time
    if verbose
        println("└─ Index building complete: $(round(total_time, digits=3))s total")
        println("   Priority vector contains $n_ms2 scans")
        println("   Memory: Metadata only, no scan data loaded")
    end

    return QuadScanPriorityIndex(
        priority_order,
        center_mz_values[priority_order],
        tic_values[priority_order],
        ms_orders[priority_order],
        Int32(n_ms2),
        Int32(n_mz_bins)
    )
end

"""
    progressive_quad_psm_collection!(scan_index::QuadScanPriorityIndex,
                                    spectra::MassSpecData,
                                    search_context::SearchContext,
                                    params::QuadTuningSearchParameters,
                                    ms_file_idx::Int64;
                                    initial_percent::Float64 = 10.0,
                                    min_psms_required::Int = 1000,
                                    verbose::Bool = true)

Progressively collect PSMs using sampling without replacement for quadrupole tuning.

# Process
1. Start with initial_percent of prioritized scans
2. If insufficient PSMs, sample next 2*initial_percent scans
3. Continue doubling until min_psms_required or all scans used
4. Never re-process scans from previous iterations

# Arguments
- `scan_index`: QuadScanPriorityIndex with prioritized scan order
- `spectra`: MassSpecData object
- `search_context`: SearchContext for library search
- `params`: Quad tuning search parameters
- `ms_file_idx`: File index for context
- `initial_percent`: Starting percentage of scans to sample
- `min_psms_required`: Minimum PSMs needed for quad modeling
- `verbose`: Enable detailed logging

# Returns
(all_psms::DataFrame, converged::Bool, scans_processed::Int)

# Notes
- Uses higher initial sampling (10%) than NCE since quad needs more data
- Targets different PSM count requirement than NCE tuning
"""
function progressive_quad_psm_collection!(
    scan_index::QuadScanPriorityIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64;
    initial_percent::Float64 = 10.0,  # Higher than NCE since quad needs more data
    min_psms_required::Int = 1000,    # Different requirement than NCE
    verbose::Bool = true
)
    total_scans = length(scan_index.scan_indices)
    scans_processed = 0  # Track position in priority vector
    iteration = 0
    converged = false
    all_psms = DataFrame()  # Accumulate PSMs across iterations

    if total_scans == 0
        if verbose
            println("│  ⚠ No scans available for processing")
        end
        return all_psms, false, 0
    end

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
        println("\n┌─ Quad Progressive PSM Collection")
        println("│  Total MS2 scans: $total_scans")
        println("│  Target PSMs: $min_psms_required")
        println("│  Initial sample: $(initial_percent)%")
        println("│  Sampling schedule: $(sample_schedule)%")
        println("│")
    end

    # Progressive sampling loop - sampling WITHOUT replacement
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

        # Extract scan indices for this chunk (NEVER re-process previous scans)
        chunk_scan_indices = scan_index.scan_indices[start_idx:end_idx]
        n_chunk_scans = length(chunk_scan_indices)

        if verbose
            println("│  ┌─ Iteration $iteration")
            println("│  │  Processing scans $start_idx to $end_idx ($n_chunk_scans scans)")
            println("│  │  Cumulative: $end_idx / $total_scans ($(round(100*end_idx/total_scans, digits=1))%)")
        end

        # Collect PSMs for this chunk - this is where scan data is FIRST loaded
        chunk_start = time()
        chunk_psms = collect_quad_psms_for_scans(
            chunk_scan_indices,
            spectra,
            search_context,
            params,
            ms_file_idx,
            verbose
        )
        chunk_time = time() - chunk_start

        # Combine with previous PSMs (accumulate across iterations)
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

        # Update position in priority vector (for next iteration)
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

"""
    collect_quad_psms_for_scans(scan_indices::Vector{Int32},
                               spectra::MassSpecData,
                               search_context::SearchContext,
                               params::QuadTuningSearchParameters,
                               ms_file_idx::Int64,
                               verbose::Bool)::DataFrame

Collect PSMs from specific scan indices for quadrupole tuning.
This is where scan data is FIRST accessed after index building.

# Arguments
- `scan_indices`: Specific scans to process (from priority vector)
- `spectra`: MassSpecData object
- `search_context`: SearchContext for library search
- `params`: Quad tuning search parameters
- `ms_file_idx`: File index for context
- `verbose`: Enable logging

# Returns
DataFrame of processed PSMs from specified scans

# Notes
- Uses IndexedMassSpecData to only process selected scans
- Eliminates post-filtering waste by filtering BEFORE library search
- Maps virtual scan indices back to actual scan numbers
"""
function collect_quad_psms_for_scans(
    scan_indices::Vector{Int32},
    spectra::MassSpecData,
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64,
    verbose::Bool
)::DataFrame

    if verbose
        println("│  │  ├─ Creating indexed view of $(length(scan_indices)) scans...")
    end

    # Create indexed view that only exposes the target scans
    indexed_start = time()
    indexed_spectra = IndexedMassSpecData(spectra, scan_indices)
    indexed_time = time() - indexed_start

    if verbose
        println("│  │  │  Indexed view creation: $(round(indexed_time, digits=3))s")
        println("│  │  ├─ Running library search on indexed scans...")
    end

    # Library search now only processes the selected scans (massive performance improvement!)
    search_start = time()
    psms = library_search(indexed_spectra, search_context, params, ms_file_idx)
    search_time = time() - search_start

    if verbose
        println("│  │  │  Library search: $(round(search_time, digits=3))s")
        println("│  │  │  PSMs found: $(nrow(psms))")
    end

    # Process PSMs using indexed spectra (no filtering needed!)
    if !isempty(psms)
        process_start = time()
        processed_psms = process_initial_psms(psms, indexed_spectra, search_context)
        process_time = time() - process_start

        if verbose
            println("│  │  │  PSM processing: $(round(process_time, digits=3))s")
            println("│  │  │  Final PSMs: $(nrow(processed_psms))")
        end

        # CRITICAL: Map virtual scan indices back to actual scan indices
        # The library search used indices 1, 2, 3... but we need the actual scan indices
        if !isempty(processed_psms) && haskey(processed_psms, :scan_idx)
            # Create mapping from virtual to actual scan indices
            scan_mapping = create_scan_mapping(indexed_spectra)

            # Map all scan_idx values from virtual to actual
            processed_psms[!, :scan_idx] = [scan_mapping[Int32(virtual_idx)] for virtual_idx in processed_psms[!, :scan_idx]]

            if verbose
                println("│  │  │  Scan indices mapped back to original numbering")
            end
        end

        return processed_psms
    end

    return DataFrame()
end

"""
    log_memory_usage(checkpoint::String)

Log current memory usage for monitoring.
"""
function log_memory_usage(checkpoint::String)
    GC.gc()  # Force collection for accurate reading
    gc_stats = Base.gc_num()
    mem_used = Sys.total_memory() - Sys.free_memory()
    println("Memory [$checkpoint]: $(round(mem_used / 1e9, digits=2)) GB")
    println("  GC collections: $(gc_stats.collect)")
    println("  GC time: $(round(gc_stats.total_time / 1e9, digits=3))s")
end