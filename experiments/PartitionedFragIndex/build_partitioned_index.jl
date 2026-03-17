"""
    build_partitioned_index(nfi::NativeFragmentIndex{T}; partition_width=5.0f0)

Reorganize a NativeFragmentIndex into a PartitionedFragmentIndex by splitting fragments
according to their precursor m/z into partitions of width `partition_width` Da.

Every partition preserves the **full** frag bin and RT bin structure from the original
index. Frag bins that have no fragments for a given partition get empty ranges
(first > last), so the exponential + binary search navigates the identical "landscape"
as the baseline. RT bins are copied verbatim since they index into frag bins at the
same positions.
"""
function build_partitioned_index(nfi::Pioneer.NativeFragmentIndex{T};
                                  partition_width::T = T(5.0)) where {T<:AbstractFloat}
    rt_bins_orig   = Pioneer.getRTBins(nfi)
    frag_bins_orig = Pioneer.getFragBins(nfi)
    fragments_orig = Pioneer.getFragments(nfi)
    n_frag_bins    = length(frag_bins_orig)

    # Step 1: Find prec_mz range across all fragments → compute n_partitions
    min_prec_mz = typemax(T)
    max_prec_mz = typemin(T)
    for frag in fragments_orig
        pmz = Pioneer.getPrecMZ(frag)
        min_prec_mz = min(min_prec_mz, pmz)
        max_prec_mz = max(max_prec_mz, pmz)
    end
    n_partitions = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / partition_width))
    println("  Prec m/z range: [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))]")
    println("  Partition width: $partition_width Da → $n_partitions partitions")

    # Step 2: Distribute fragments preserving full frag bin structure.
    # Walk each frag bin in the original and distribute its fragments across partitions.
    # Each partition gets the same number of frag bins (n_frag_bins) in the same order,
    # but with partition-local first_bin/last_bin ranges.
    partition_fragments = [Pioneer.IndexFragment{T}[] for _ in 1:n_partitions]
    partition_fb_first  = [Vector{UInt32}(undef, n_frag_bins) for _ in 1:n_partitions]
    partition_fb_last   = [Vector{UInt32}(undef, n_frag_bins) for _ in 1:n_partitions]
    starts = Vector{UInt32}(undef, n_partitions)  # reusable per-frag-bin buffer

    for j in 1:n_frag_bins
        fbin = frag_bins_orig[j]
        frag_range = Pioneer.getSubBinRange(fbin)

        # Snapshot current end of each partition's fragments vector
        for k in 1:n_partitions
            starts[k] = UInt32(length(partition_fragments[k]) + 1)
        end

        # Distribute this frag bin's fragments to their respective partitions
        for fi in frag_range
            frag = fragments_orig[fi]
            pmz = Pioneer.getPrecMZ(frag)
            k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1,
                       1, n_partitions)
            push!(partition_fragments[k], frag)
        end

        # Record first/last for each partition's copy of frag bin j
        # When no fragments were added, new_end < starts[k] → empty UnitRange
        for k in 1:n_partitions
            new_end = UInt32(length(partition_fragments[k]))
            partition_fb_first[k][j] = starts[k]
            partition_fb_last[k][j]  = new_end
        end
    end

    # Step 3: Assemble NativeFragmentIndex per partition.
    # Each partition gets the full frag bin array (same lb/ub, partition-local ranges)
    # and a verbatim copy of the original RT bins (indices into frag bins are unchanged).
    partitions = Vector{Pioneer.NativeFragmentIndex{T}}(undef, n_partitions)
    for k in 1:n_partitions
        p_frag_bins = Vector{Pioneer.FragIndexBin{T}}(undef, n_frag_bins)
        for j in 1:n_frag_bins
            orig = frag_bins_orig[j]
            p_frag_bins[j] = Pioneer.FragIndexBin{T}(
                Pioneer.getLow(orig), Pioneer.getHigh(orig),
                partition_fb_first[k][j], partition_fb_last[k][j])
        end
        partitions[k] = Pioneer.NativeFragmentIndex{T}(
            p_frag_bins, copy(rt_bins_orig), partition_fragments[k])
    end

    return PartitionedFragmentIndex{T}(partitions, min_prec_mz, partition_width, n_partitions)
end

"""
    build_partitioned_index_from_lib(spec_lib; partition_width=5.0f0,
        frag_bin_tol_ppm=10.0f0, rt_bin_tol=1.0f0,
        rank_to_score=UInt8[8,4,4,2,2,1,1],
        y_start_index=UInt8(4), b_start_index=UInt8(3))

Build a LocalPartitionedFragmentIndex from scratch using the spectral library's
DetailedFrag data. Each partition gets its own independently-constructed index
with LocalFragment entries (UInt16 local_id + UInt8 score = 4 bytes).

Precursor IDs are remapped to partition-local UInt16 values (1..N, N ≤ 65535).
Partitions that would exceed 65535 unique precursors are automatically split.
"""
function build_partitioned_index_from_lib(
    spec_lib::Pioneer.SpectralLibrary;
    partition_width::Float32 = 5.0f0,
    frag_bin_tol_ppm::Float32 = 10.0f0,
    rt_bin_tol::Float32 = 1.0f0,
    rank_to_score::Vector{UInt8} = UInt8[8, 4, 4, 2, 2, 1, 1],
    y_start_index::UInt8 = UInt8(4),
    b_start_index::UInt8 = UInt8(3),
)
    precursors = Pioneer.getPrecursors(spec_lib)
    frag_lookup = Pioneer.getFragmentLookupTable(spec_lib)
    detailed_frags = Pioneer.getFragments(frag_lookup)
    prec_mzs = Pioneer.getMz(precursors)
    prec_irts = Pioneer.getIrt(precursors)
    n_precursors = length(prec_mzs)
    max_rank = length(rank_to_score)

    # ── Step 1: Assign precursors to initial partitions by prec_mz ───────────
    min_prec_mz = Float32(Inf)
    max_prec_mz = Float32(-Inf)
    for i in 1:n_precursors
        pmz = prec_mzs[i]
        min_prec_mz = min(min_prec_mz, pmz)
        max_prec_mz = max(max_prec_mz, pmz)
    end
    n_initial = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / partition_width))
    println("  Prec m/z range: [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))]")
    println("  Partition width: $partition_width Da → $n_initial initial partitions")

    # Collect global precursor IDs per initial partition (sorted by prec_mz within)
    initial_partition_pids = [UInt32[] for _ in 1:n_initial]
    for pid in UInt32(1):UInt32(n_precursors)
        pmz = prec_mzs[pid]
        k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1, 1, n_initial)
        push!(initial_partition_pids[k], pid)
    end

    # ── Step 2: Split partitions exceeding MAX_LOCAL_PRECS ───────────────────
    # Build final partition specs: each is a Vector{UInt32} of global pids
    final_partition_pids = Vector{UInt32}[]
    for pids in initial_partition_pids
        if length(pids) <= MAX_LOCAL_PRECS
            push!(final_partition_pids, pids)
        else
            # Sort by prec_mz and split into chunks of MAX_LOCAL_PRECS
            sort!(pids, by = pid -> prec_mzs[pid])
            for i in 1:MAX_LOCAL_PRECS:length(pids)
                chunk_end = min(i + MAX_LOCAL_PRECS - 1, length(pids))
                push!(final_partition_pids, pids[i:chunk_end])
            end
        end
    end
    n_partitions = length(final_partition_pids)
    if n_partitions != n_initial
        println("  Split to $n_partitions partitions ($(n_partitions - n_initial) extra from UInt16 limit)")
    end

    # ── Step 3: Build SimpleFrags + local ID mapping per partition ────────────
    partition_frags = [Pioneer.SimpleFrag{Float32}[] for _ in 1:n_partitions]
    partition_local_to_global = [UInt32[] for _ in 1:n_partitions]
    # global_to_local: reusable buffer (cleared per partition)
    global_to_local = Dict{UInt32, UInt16}()

    for k in 1:n_partitions
        pids = final_partition_pids[k]
        empty!(global_to_local)
        local_to_global = zeros(UInt32, length(pids))
        for (i, pid) in enumerate(pids)
            lid = UInt16(i)
            global_to_local[pid] = lid
            local_to_global[i] = pid
        end

        for pid in pids
            pmz = prec_mzs[pid]
            pirt = prec_irts[pid]
            frag_range = Pioneer.getPrecFragRange(frag_lookup, pid)
            lid = global_to_local[pid]

            rank = 0
            for fi in frag_range
                dfrag = detailed_frags[fi]

                # Apply fragment index ion-type filters
                if dfrag.is_y
                    dfrag.ion_position < y_start_index && continue
                elseif dfrag.is_b
                    dfrag.ion_position < b_start_index && continue
                elseif dfrag.is_p
                    continue
                end
                dfrag.is_isotope && continue

                rank += 1
                rank > max_rank && break

                # Store with local ID in the SimpleFrag's prec_id field
                # (we'll extract it during index build)
                push!(partition_frags[k], Pioneer.SimpleFrag{Float32}(
                    Pioneer.getMz(dfrag),
                    UInt32(lid),  # local ID stored as UInt32 temporarily
                    pmz,
                    pirt,
                    UInt8(0),     # charge unused
                    rank_to_score[rank],
                ))
            end
        end

        partition_local_to_global[k] = local_to_global
    end

    total_frags = sum(length, partition_frags)
    println("  Total SimpleFrags across partitions: $total_frags")
    max_local = maximum(length, final_partition_pids)
    println("  Max precursors in a partition: $max_local")

    # ── Step 4: Build LocalPartition per partition ────────────────────────────
    partitions = Vector{LocalPartition{Float32}}(undef, n_partitions)

    for k in 1:n_partitions
        frags_k = partition_frags[k]
        l2g = partition_local_to_global[k]
        n_local = UInt16(length(l2g))

        if isempty(frags_k)
            partitions[k] = LocalPartition{Float32}(
                Pioneer.FragIndexBin{Float32}[],
                Pioneer.FragIndexBin{Float32}[],
                LocalFragment[],
                l2g,
                n_local,
            )
            continue
        end

        partitions[k] = _build_local_partition(frags_k, l2g, n_local,
                                                frag_bin_tol_ppm, rt_bin_tol)
    end

    part_stats = [(length(getFragBins(p)), length(getRTBins(p)),
                   length(getFragments(p))) for p in partitions]
    total_fb = sum(first, part_stats)
    total_rt = sum(x -> x[2], part_stats)
    total_fr = sum(last, part_stats)
    println("  Total frag_bins: $total_fb  rt_bins: $total_rt  fragments: $total_fr")
    frag_mem = sum(sizeof(getFragments(p)) for p in partitions)
    bin_mem = sum(sizeof(getFragBins(p)) + sizeof(getRTBins(p)) for p in partitions)
    l2g_mem = sum(sizeof(p.local_to_global) for p in partitions)
    println("  Memory: fragments=$(round(frag_mem/1024^2, digits=1))MB  bins=$(round(bin_mem/1024^2, digits=1))MB  l2g=$(round(l2g_mem/1024^2, digits=1))MB")

    # Compute per-partition prec_mz bounds (sorted by prec_mz_min)
    partition_bounds = Vector{Tuple{Float32, Float32}}(undef, n_partitions)
    for k in 1:n_partitions
        pids = final_partition_pids[k]
        if isempty(pids)
            partition_bounds[k] = (Float32(Inf), Float32(-Inf))
        else
            pmin = Float32(Inf)
            pmax = Float32(-Inf)
            for pid in pids
                pmz = prec_mzs[pid]
                pmin = min(pmin, pmz)
                pmax = max(pmax, pmz)
            end
            partition_bounds[k] = (pmin, pmax)
        end
    end

    return LocalPartitionedFragmentIndex{Float32}(partitions, partition_bounds, n_partitions)
end

"""
Build a CompactPartition from a vector of SimpleFrags using the standard
hierarchical binning algorithm (sort by iRT → bin by RT → sort by frag m/z →
bin by m/z), producing CompactFragment entries (prec_id + score only).
"""
function _build_compact_partition(
    frag_ions::Vector{Pioneer.SimpleFrag{Float32}},
    frag_bin_tol_ppm::Float32,
    rt_bin_tol::Float32,
)
    sort!(frag_ions, by = x -> Pioneer.getIRT(x))

    n = length(frag_ions)
    compact_fragments = Vector{CompactFragment}(undef, n)
    rt_bins = Vector{Pioneer.FragIndexBin{Float32}}(undef, n)   # upper bound
    frag_bins = Vector{Pioneer.FragIndexBin{Float32}}(undef, n) # upper bound
    rt_bin_idx = 0
    frag_bin_idx = 0

    start_idx = 1
    start_irt = Pioneer.getIRT(frag_ions[1])

    for i in 1:n
        stop_irt = Pioneer.getIRT(frag_ions[i])
        if (stop_irt - start_irt > rt_bin_tol) && (i > start_idx)
            stop_idx = i - 1
            stop_irt_val = Pioneer.getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> Pioneer.getMZ(x))
            first_fb = frag_bin_idx + 1
            frag_bin_idx = _build_compact_frag_bins!(compact_fragments, frag_bins,
                frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
            rt_bin_idx += 1
            rt_bins[rt_bin_idx] = Pioneer.FragIndexBin{Float32}(
                start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))
            start_idx = i
            start_irt = Pioneer.getIRT(frag_ions[i])
        end
    end

    # Last RT bin
    stop_idx = n
    stop_irt_val = Pioneer.getIRT(frag_ions[stop_idx])
    sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> Pioneer.getMZ(x))
    first_fb = frag_bin_idx + 1
    frag_bin_idx = _build_compact_frag_bins!(compact_fragments, frag_bins,
        frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
    rt_bin_idx += 1
    rt_bins[rt_bin_idx] = Pioneer.FragIndexBin{Float32}(
        start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))

    return CompactPartition{Float32}(
        frag_bins[1:frag_bin_idx],
        rt_bins[1:rt_bin_idx],
        compact_fragments,
    )
end

"""
Build fragment m/z bins within an RT bin, producing CompactFragment entries.
Sorts each bin's fragments by prec_mz (preserves the standard binning algorithm).
Returns updated frag_bin_idx.
"""
function _build_compact_frag_bins!(
    compact_fragments::Vector{CompactFragment},
    frag_bins::Vector{Pioneer.FragIndexBin{Float32}},
    frag_bin_idx::Int,
    frag_ions::Vector{Pioneer.SimpleFrag{Float32}},
    start::Int, stop::Int,
    frag_bin_tol_ppm::Float32,
)
    start_idx = start
    start_mz = Pioneer.getMZ(frag_ions[start])

    for i in start:stop
        stop_mz = Pioneer.getMZ(frag_ions[i])
        diff_mz = stop_mz - start_mz
        mean_mz = (stop_mz + start_mz) / 2
        if (diff_mz / (mean_mz / 1.0f6) > frag_bin_tol_ppm) && (i > start_idx)
            bin_stop = i - 1
            bin_stop_mz = Pioneer.getMZ(frag_ions[bin_stop])
            sort!(@view(frag_ions[start_idx:bin_stop]), by = x -> Pioneer.getPrecMZ(x))
            frag_bin_idx += 1
            frag_bins[frag_bin_idx] = Pioneer.FragIndexBin{Float32}(
                start_mz, bin_stop_mz, UInt32(start_idx), UInt32(bin_stop))
            for idx in start_idx:bin_stop
                sf = frag_ions[idx]
                compact_fragments[idx] = CompactFragment(
                    Pioneer.getPrecID(sf), Pioneer.getScore(sf))
            end
            start_idx = i
            start_mz = Pioneer.getMZ(frag_ions[i])
        end
    end

    # Last frag bin
    stop_mz = Pioneer.getMZ(frag_ions[stop])
    sort!(@view(frag_ions[start_idx:stop]), by = x -> Pioneer.getPrecMZ(x))
    frag_bin_idx += 1
    frag_bins[frag_bin_idx] = Pioneer.FragIndexBin{Float32}(
        start_mz, stop_mz, UInt32(start_idx), UInt32(stop))
    for idx in start_idx:stop
        sf = frag_ions[idx]
        compact_fragments[idx] = CompactFragment(
            Pioneer.getPrecID(sf), Pioneer.getScore(sf))
    end

    return frag_bin_idx
end

# ── LocalPartition build functions ───────────────────────────────────────────

"""
Build a LocalPartition from SimpleFrags whose prec_id field already contains
local UInt16 IDs (stored as UInt32). Produces LocalFragment entries.
"""
function _build_local_partition(
    frag_ions::Vector{Pioneer.SimpleFrag{Float32}},
    local_to_global::Vector{UInt32},
    n_local::UInt16,
    frag_bin_tol_ppm::Float32,
    rt_bin_tol::Float32,
)
    sort!(frag_ions, by = x -> Pioneer.getIRT(x))

    n = length(frag_ions)
    local_fragments = Vector{LocalFragment}(undef, n)
    rt_bins = Vector{Pioneer.FragIndexBin{Float32}}(undef, n)
    frag_bins = Vector{Pioneer.FragIndexBin{Float32}}(undef, n)
    rt_bin_idx = 0
    frag_bin_idx = 0

    start_idx = 1
    start_irt = Pioneer.getIRT(frag_ions[1])

    for i in 1:n
        stop_irt = Pioneer.getIRT(frag_ions[i])
        if (stop_irt - start_irt > rt_bin_tol) && (i > start_idx)
            stop_idx = i - 1
            stop_irt_val = Pioneer.getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> Pioneer.getMZ(x))
            first_fb = frag_bin_idx + 1
            frag_bin_idx = _build_local_frag_bins!(local_fragments, frag_bins,
                frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
            rt_bin_idx += 1
            rt_bins[rt_bin_idx] = Pioneer.FragIndexBin{Float32}(
                start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))
            start_idx = i
            start_irt = Pioneer.getIRT(frag_ions[i])
        end
    end

    # Last RT bin
    stop_idx = n
    stop_irt_val = Pioneer.getIRT(frag_ions[stop_idx])
    sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> Pioneer.getMZ(x))
    first_fb = frag_bin_idx + 1
    frag_bin_idx = _build_local_frag_bins!(local_fragments, frag_bins,
        frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
    rt_bin_idx += 1
    rt_bins[rt_bin_idx] = Pioneer.FragIndexBin{Float32}(
        start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))

    return LocalPartition{Float32}(
        frag_bins[1:frag_bin_idx],
        rt_bins[1:rt_bin_idx],
        local_fragments,
        local_to_global,
        n_local,
    )
end

"""
Build fragment m/z bins producing LocalFragment entries (UInt16 local IDs).
The SimpleFrag prec_id field already contains the local ID as UInt32.
"""
function _build_local_frag_bins!(
    local_fragments::Vector{LocalFragment},
    frag_bins::Vector{Pioneer.FragIndexBin{Float32}},
    frag_bin_idx::Int,
    frag_ions::Vector{Pioneer.SimpleFrag{Float32}},
    start::Int, stop::Int,
    frag_bin_tol_ppm::Float32,
)
    start_idx = start
    start_mz = Pioneer.getMZ(frag_ions[start])

    for i in start:stop
        stop_mz = Pioneer.getMZ(frag_ions[i])
        diff_mz = stop_mz - start_mz
        mean_mz = (stop_mz + start_mz) / 2
        if (diff_mz / (mean_mz / 1.0f6) > frag_bin_tol_ppm) && (i > start_idx)
            bin_stop = i - 1
            bin_stop_mz = Pioneer.getMZ(frag_ions[bin_stop])
            sort!(@view(frag_ions[start_idx:bin_stop]), by = x -> Pioneer.getPrecMZ(x))
            frag_bin_idx += 1
            frag_bins[frag_bin_idx] = Pioneer.FragIndexBin{Float32}(
                start_mz, bin_stop_mz, UInt32(start_idx), UInt32(bin_stop))
            for idx in start_idx:bin_stop
                sf = frag_ions[idx]
                local_fragments[idx] = LocalFragment(
                    UInt16(Pioneer.getPrecID(sf)), Pioneer.getScore(sf))
            end
            start_idx = i
            start_mz = Pioneer.getMZ(frag_ions[i])
        end
    end

    # Last frag bin
    stop_mz = Pioneer.getMZ(frag_ions[stop])
    sort!(@view(frag_ions[start_idx:stop]), by = x -> Pioneer.getPrecMZ(x))
    frag_bin_idx += 1
    frag_bins[frag_bin_idx] = Pioneer.FragIndexBin{Float32}(
        start_mz, stop_mz, UInt32(start_idx), UInt32(stop))
    for idx in start_idx:stop
        sf = frag_ions[idx]
        local_fragments[idx] = LocalFragment(
            UInt16(Pioneer.getPrecID(sf)), Pioneer.getScore(sf))
    end

    return frag_bin_idx
end
