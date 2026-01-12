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
Memory-efficient scan selection for NCE tuning search.
Uses intelligent RT-binning with TIC-based prioritization.
"""

"""
    NCEScanPriorityIndex

Contains prioritized scan indices for memory-efficient NCE tuning.
Only stores metadata, never loads actual scan data until needed.
"""
struct NCEScanPriorityIndex
    scan_indices::Vector{Int32}        # Ordered scan indices by priority
    rt_values::Vector{Float32}         # RT values for each scan
    tic_values::Vector{Float32}        # TIC values for each scan
    ms_orders::Vector{UInt8}           # MS order for verification
    total_ms2_count::Int32             # Total MS2 scans
    n_rt_bins::Int32                   # Number of RT bins used
end

"""
    build_nce_scan_priority_index(spectra::MassSpecData;
                                  n_rt_bins::Int = 15,
                                  target_ms_order::UInt8 = UInt8(2),
                                  verbose::Bool = true)::NCEScanPriorityIndex

Build prioritized scan index for NCE tuning without loading scan data.

# Process
1. Extract metadata (RT, TIC, MS order) for all scans - NO peak data
2. Filter for MS2 scans only
3. Divide RT range into N bins
4. Within each bin, sort scans by TIC (descending)
5. Create priority vector by round-robin from bins

# Arguments
- `spectra`: MassSpecData object
- `n_rt_bins`: Number of RT bins for even coverage (default: 15)
- `target_ms_order`: MS order to filter for (default: 2)

# Returns
NCEScanPriorityIndex containing ordered scan indices and metadata
"""
function build_nce_scan_priority_index(
    spectra::MassSpecData;
    n_rt_bins::Int = 15,
    target_ms_order::UInt8 = UInt8(2)
)::NCEScanPriorityIndex

    # Step 1: Extract metadata (NO peak data loaded)
    rt_values = getRetentionTimes(spectra)
    tic_values = getTICs(spectra)
    ms_orders = getMsOrders(spectra)

    # Step 2: Filter for MS2 scans
    ms2_mask = ms_orders .== target_ms_order
    ms2_indices = findall(ms2_mask)
    n_ms2 = length(ms2_indices)

    if n_ms2 == 0
        return NCEScanPriorityIndex(
            Int32[], Float32[], Float32[], UInt8[],
            Int32(0), Int32(n_rt_bins)
        )
    end

    # Step 3: Create RT bins
    ms2_rt = rt_values[ms2_indices]
    ms2_tic = tic_values[ms2_indices]

    rt_min, rt_max = extrema(ms2_rt)

    # Handle edge case of single RT value
    if rt_min == rt_max
        priority_order = Int32.(ms2_indices[sortperm(ms2_tic, rev=true)])
        return NCEScanPriorityIndex(
            priority_order,
            rt_values[priority_order],
            tic_values[priority_order],
            ms_orders[priority_order],
            Int32(n_ms2),
            Int32(1)
        )
    end

    bin_width = (rt_max - rt_min) / n_rt_bins

    # Assign each MS2 scan to a bin
    bin_assignments = zeros(Int, n_ms2)
    for (i, rt) in enumerate(ms2_rt)
        bin_idx = min(max(1, ceil(Int, (rt - rt_min) / bin_width)), n_rt_bins)
        bin_assignments[i] = bin_idx
    end

    # Step 4: Sort within bins by TIC and create priority order
    # Group scans by bin
    bins = [Int32[] for _ in 1:n_rt_bins]
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
        for bin_idx in 1:n_rt_bins
            if round <= length(bins[bin_idx])
                push!(priority_order, bins[bin_idx][round])
            end
        end
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

"""
    progressive_nce_psm_collection!(scan_index::NCEScanPriorityIndex,
                                   spectra::MassSpecData,
                                   search_context::SearchContext,
                                   params::NceTuningSearchParameters,
                                   ms_file_idx::Int64;
                                   )

Progressively collect PSMs using sampling without replacement.

# Process
1. Start with initial_percent of prioritized scans
2. If insufficient PSMs, sample next 2*initial_percent scans
3. Continue doubling until min_psms_required or all scans used
4. Never re-process scans from previous iterations

# Arguments
- `scan_index`: NCEScanPriorityIndex with prioritized scan order
- `spectra`: MassSpecData object
- `search_context`: SearchContext for library search
- `params`: NCE tuning search parameters
- `ms_file_idx`: File index for context

# Returns
(all_psms::DataFrame, converged::Bool, scans_processed::Int)
"""
function progressive_nce_psm_collection!(
    scan_index::NCEScanPriorityIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NceTuningSearchParameters,
    ms_file_idx::Int64
)
    total_scans = length(scan_index.scan_indices)
    scans_processed = 0  # Track position in priority vector
    iteration = 0
    converged = false
    all_psms = DataFrame()  # Accumulate PSMs across iterations

    if total_scans == 0
        return all_psms, false, 0
    end

    # Calculate sampling parameters from NCE tuning configuration
    min_psms_required = params.min_psms
    initial_percent = max(params.initial_percent, 100.0 * params.min_initial_scans / total_scans)

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

    # Progressive sampling loop - sampling WITHOUT replacement
    for chunk_percent in sample_schedule
        iteration += 1

        # Calculate range for THIS iteration (no overlap with previous)
        n_scans_this_chunk = ceil(Int, total_scans * chunk_percent / 100)
        start_idx = scans_processed + 1
        end_idx = min(scans_processed + n_scans_this_chunk, total_scans)

        if start_idx > total_scans
            break
        end

        # Extract scan indices for this chunk (NEVER re-process previous scans)
        chunk_scan_indices = scan_index.scan_indices[start_idx:end_idx]

        # Collect PSMs for this chunk - this is where scan data is FIRST loaded
        chunk_psms = collect_nce_psms_for_scans(
            chunk_scan_indices,
            spectra,
            search_context,
            params,
            ms_file_idx
        )

        # Combine with previous PSMs (accumulate across iterations)
        if !isempty(chunk_psms)
            all_psms = isempty(all_psms) ? chunk_psms : vcat(all_psms, chunk_psms)
        end

        total_psms = nrow(all_psms)

        # Check convergence
        if total_psms >= min_psms_required
            converged = true
            break
        end

        # Update position in priority vector (for next iteration)
        scans_processed = end_idx
    end


    return all_psms, converged, scans_processed
end

"""
    collect_nce_psms_for_scans(scan_indices::Vector{Int32},
                              spectra::MassSpecData,
                              search_context::SearchContext,
                              params::NceTuningSearchParameters,
                              ms_file_idx::Int64,
                              )::DataFrame

Collect PSMs from specific scan indices for NCE tuning.
This is where scan data is FIRST accessed after index building.

# Arguments
- `scan_indices`: Specific scans to process (from priority vector)
- `spectra`: MassSpecData object
- `search_context`: SearchContext for library search
- `params`: NCE tuning search parameters
- `ms_file_idx`: File index for context

# Returns
DataFrame of processed PSMs from specified scans
"""
function collect_nce_psms_for_scans(
    scan_indices::Vector{Int32},
    spectra::MassSpecData,
    search_context::SearchContext,
    params::NceTuningSearchParameters,
    ms_file_idx::Int64
)::DataFrame

    # Create indexed view that only exposes the target scans
    indexed_spectra = IndexedMassSpecData(spectra, scan_indices)

    # Library search now only processes the selected scans (massive performance improvement!)
    psms = library_search(indexed_spectra, search_context, params, ms_file_idx)

    # Process PSMs using indexed spectra (no filtering needed!)
    if !isempty(psms)
        psms = process_psms!(psms, indexed_spectra, search_context, params)

        # CRITICAL: Map virtual scan indices back to actual scan indices
        # The library search used indices 1, 2, 3... but we need the actual scan indices
        if !isempty(psms) && "scan_idx" in names(psms)
            # Create mapping from virtual to actual scan indices
            scan_mapping = create_scan_mapping(indexed_spectra)

            # Map all scan_idx values from virtual to actual
            psms[!, :scan_idx] = [scan_mapping[Int32(virtual_idx)] for virtual_idx in psms[!, :scan_idx]]
        end
    end

    return psms
end