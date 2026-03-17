#
# Build a FixedBinPartitionedIndex from the spectral library.
#
# Reuses the existing partition assignment (prec_mz → partition, local IDs)
# from build_partitioned_index.jl, then restructures each partition into
# fixed-width m/z bins with CSR-style pointers.
#

"""
    build_fixed_bin_index(spec_lib; partition_width=5.0f0, rt_bin_tol=1.0f0,
        rank_to_score=UInt8[8,4,4,2,2,1,1], y_start_index=UInt8(4), b_start_index=UInt8(3))

Build a FixedBinPartitionedIndex with fixed 0.005 Da m/z bins per RT bin.
"""
function build_fixed_bin_index(
    spec_lib::Pioneer.SpectralLibrary;
    partition_width::Float32 = 5.0f0,
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

    # ── Step 1: Partition assignment (same as build_partitioned_index_from_lib) ─
    min_prec_mz = Float32(Inf)
    max_prec_mz = Float32(-Inf)
    for i in 1:n_precursors
        pmz = prec_mzs[i]
        min_prec_mz = min(min_prec_mz, pmz)
        max_prec_mz = max(max_prec_mz, pmz)
    end
    n_initial = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / partition_width))

    initial_partition_pids = [UInt32[] for _ in 1:n_initial]
    for pid in UInt32(1):UInt32(n_precursors)
        pmz = prec_mzs[pid]
        k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1, 1, n_initial)
        push!(initial_partition_pids[k], pid)
    end

    # Split partitions exceeding UInt16 limit
    final_partition_pids = Vector{UInt32}[]
    for pids in initial_partition_pids
        if length(pids) <= MAX_LOCAL_PRECS
            push!(final_partition_pids, pids)
        else
            sort!(pids, by = pid -> prec_mzs[pid])
            for i in 1:MAX_LOCAL_PRECS:length(pids)
                chunk_end = min(i + MAX_LOCAL_PRECS - 1, length(pids))
                push!(final_partition_pids, pids[i:chunk_end])
            end
        end
    end
    n_partitions = length(final_partition_pids)
    println("  Prec m/z range: [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))]")
    println("  $n_partitions partitions ($(n_partitions - n_initial) extra from UInt16 limit)")

    # ── Step 2: For each partition, collect SimpleFrags with local IDs ────────
    # Group by local_id so we can sort by iRT then by frag m/z

    # We'll collect (frag_mz, local_id, prec_irt, score) tuples per partition
    # and then bin them.

    PartFrag = @NamedTuple{frag_mz::Float32, local_id::UInt16, prec_irt::Float32, score::UInt8}
    partition_frags_raw = [PartFrag[] for _ in 1:n_partitions]
    partition_local_to_global = [UInt32[] for _ in 1:n_partitions]

    global_to_local = Dict{UInt32, UInt16}()

    for k in 1:n_partitions
        pids = final_partition_pids[k]
        empty!(global_to_local)
        l2g = zeros(UInt32, length(pids))
        for (i, pid) in enumerate(pids)
            global_to_local[pid] = UInt16(i)
            l2g[i] = pid
        end
        partition_local_to_global[k] = l2g

        for pid in pids
            pmz = prec_mzs[pid]
            pirt = prec_irts[pid]
            lid = global_to_local[pid]
            frag_range = Pioneer.getPrecFragRange(frag_lookup, pid)

            rank = 0
            for fi in frag_range
                dfrag = detailed_frags[fi]
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

                push!(partition_frags_raw[k], (
                    frag_mz = Pioneer.getMz(dfrag),
                    local_id = lid,
                    prec_irt = pirt,
                    score = rank_to_score[rank],
                ))
            end
        end
    end

    total_frags = sum(length, partition_frags_raw)
    println("  Total fragments: $total_frags")

    # ── Step 3: Build FixedBinPartition for each partition ────────────────────
    partitions = Vector{FixedBinPartition}(undef, n_partitions)
    partition_bounds = Vector{Tuple{Float32, Float32}}(undef, n_partitions)

    total_mem_bins = 0
    total_mem_frags = 0

    for k in 1:n_partitions
        pids = final_partition_pids[k]
        frags = partition_frags_raw[k]
        l2g = partition_local_to_global[k]
        n_local = UInt16(length(l2g))

        if isempty(pids)
            partition_bounds[k] = (Float32(Inf), Float32(-Inf))
            partitions[k] = FixedBinPartition(
                FixedBinRTData[], l2g, n_local, 0.0f0, Int32(0))
            continue
        end

        # Compute prec_mz bounds
        pmin = minimum(pid -> prec_mzs[pid], pids)
        pmax = maximum(pid -> prec_mzs[pid], pids)
        partition_bounds[k] = (pmin, pmax)

        if isempty(frags)
            partitions[k] = FixedBinPartition(
                FixedBinRTData[], l2g, n_local, 0.0f0, Int32(0))
            continue
        end

        # Find fragment m/z range for this partition
        frag_mz_min = minimum(f -> f.frag_mz, frags)
        frag_mz_max = maximum(f -> f.frag_mz, frags)
        # Add small margin
        mz_min = floor(frag_mz_min / FIXED_BIN_WIDTH) * FIXED_BIN_WIDTH
        n_mz_bins = ceil(Int32, (frag_mz_max - mz_min) / FIXED_BIN_WIDTH) + Int32(1)

        # Sort fragments by iRT
        sort!(frags, by = f -> f.prec_irt)

        # Build RT bins by walking sorted fragments
        rt_data_list = FixedBinRTData[]
        n = length(frags)
        start_idx = 1
        start_irt = frags[1].prec_irt

        function _build_rt_bin!(frag_slice, irt_lo, irt_hi)
            # Sort this RT bin's fragments by frag_mz
            sort!(frag_slice, by = f -> f.frag_mz)

            # Build CSR pointers: frag_ptrs[bin] = index of first fragment in bin
            frag_ptrs = Vector{UInt16}(undef, n_mz_bins + 1)
            local_fragments = Vector{LocalFragment}(undef, length(frag_slice))

            frag_out_idx = UInt16(1)
            frag_in_idx = 1

            for bin in Int32(1):n_mz_bins
                frag_ptrs[bin] = frag_out_idx
                bin_mz_hi = mz_min + bin * FIXED_BIN_WIDTH

                while frag_in_idx <= length(frag_slice) && frag_slice[frag_in_idx].frag_mz < bin_mz_hi
                    f = frag_slice[frag_in_idx]
                    local_fragments[frag_out_idx] = LocalFragment(f.local_id, f.score)
                    frag_out_idx += UInt16(1)
                    frag_in_idx += 1
                end
            end
            frag_ptrs[n_mz_bins + 1] = frag_out_idx

            push!(rt_data_list, FixedBinRTData(
                frag_ptrs, local_fragments[1:frag_out_idx-1], irt_lo, irt_hi))
        end

        for i in 1:n
            stop_irt = frags[i].prec_irt
            if (stop_irt - start_irt > rt_bin_tol) && (i > start_idx)
                stop_idx = i - 1
                _build_rt_bin!(frags[start_idx:stop_idx], start_irt, frags[stop_idx].prec_irt)
                start_idx = i
                start_irt = frags[i].prec_irt
            end
        end
        # Last RT bin
        _build_rt_bin!(frags[start_idx:n], start_irt, frags[n].prec_irt)

        total_mem_bins += sum(sizeof(rd.frag_ptrs) for rd in rt_data_list)
        total_mem_frags += sum(sizeof(rd.fragments) for rd in rt_data_list)

        partitions[k] = FixedBinPartition(
            rt_data_list, l2g, n_local, mz_min, n_mz_bins)
    end

    println("  Memory: frag_ptrs=$(round(total_mem_bins/1024^2, digits=1))MB  fragments=$(round(total_mem_frags/1024^2, digits=1))MB  l2g=$(round(sum(sizeof(p.local_to_global) for p in partitions)/1024^2, digits=1))MB")
    println("  Total: $(round((total_mem_bins + total_mem_frags + sum(sizeof(p.local_to_global) for p in partitions))/1024^2, digits=1))MB")

    return FixedBinPartitionedIndex(partitions, partition_bounds, n_partitions)
end
