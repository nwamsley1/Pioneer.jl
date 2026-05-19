using Pioneer

const DEFAULT_INDEX_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted.poin/partitioned_fragment_index.jls"

index_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INDEX_PATH

function print_examples(label, d; n = 10)
    println(label, ": ", length(d))
    shown = 0
    for key in sort!(collect(keys(d)))
        println("  ", key, " => ", d[key])
        shown += 1
        shown >= n && break
    end
end

function main(index_path)
    pfi = Pioneer.deserialize_from_jls(index_path)

    precursor_locations = Dict{UInt32, Vector{Tuple{Int, Int}}}()
    precursor_partitions = Dict{UInt32, Set{Int}}()
    empty_partitions = Int[]
    rt_bins_with_duplicate_precursor_fragments = 0
    total_rt_bins = 0
    total_fragment_entries = 0

    for (partition_idx, partition) in enumerate(Pioneer.getPartitions(pfi))
        isempty(partition.local_to_global) && push!(empty_partitions, partition_idx)

        fragments = Pioneer.getFragments(partition)
        rt_bins = Pioneer.getRTBins(partition)
        total_fragment_entries += length(fragments)
        total_rt_bins += length(rt_bins)

        for (rt_bin_idx, rt_bin) in enumerate(rt_bins)
            seen_in_bin = Set{UInt32}()
            n_fragment_entries_in_rt_bin = 0
            first_frag_bin = Int(rt_bin.first_bin)
            last_frag_bin = Int(rt_bin.last_bin)
            frag_bins = Pioneer.getFragBins(partition)

            if first_frag_bin <= last_frag_bin
                for frag_bin_idx in first_frag_bin:last_frag_bin
                    first_frag = Int(frag_bins.first_bins[frag_bin_idx])
                    last_frag = Int(frag_bins.last_bins[frag_bin_idx])
                    n_fragment_entries_in_rt_bin += max(0, last_frag - first_frag + 1)

                    for frag_idx in first_frag:last_frag
                        local_id = Int(fragments[frag_idx].local_id)
                        global_id = partition.local_to_global[local_id]
                        push!(seen_in_bin, global_id)
                    end
                end
            end

            rt_bins_with_duplicate_precursor_fragments +=
                length(seen_in_bin) < n_fragment_entries_in_rt_bin ? 1 : 0

            for global_id in seen_in_bin
                push!(get!(precursor_locations, global_id, Tuple{Int, Int}[]),
                      (partition_idx, rt_bin_idx))
                push!(get!(precursor_partitions, global_id, Set{Int}()), partition_idx)
            end
        end
    end

    multi_location = Dict{UInt32, Vector{Tuple{Int, Int}}}()
    multi_partition = Dict{UInt32, Vector{Int}}()

    for (precursor_id, locations) in precursor_locations
        unique_locations = unique(locations)
        if length(unique_locations) > 1
            multi_location[precursor_id] = unique_locations
        end
    end

    for (precursor_id, partitions) in precursor_partitions
        if length(partitions) > 1
            multi_partition[precursor_id] = sort!(collect(partitions))
        end
    end

    println("index_path: ", index_path)
    println("n_partitions: ", pfi.n_partitions)
    println("total_rt_bins: ", total_rt_bins)
    println("total_fragment_entries: ", total_fragment_entries)
    println("unique_precursors_seen: ", length(precursor_locations))
    println("empty_partitions: ", length(empty_partitions))
    if !isempty(empty_partitions)
        println("empty_partition_examples: ", empty_partitions[1:min(end, 20)])
    end
    println("rt_bins_with_duplicate_precursor_fragments: ", rt_bins_with_duplicate_precursor_fragments)
    print_examples("precursors_in_multiple_partitions", multi_partition)
    print_examples("precursors_in_multiple_partition_rt_bin_locations", multi_location)
end

main(index_path)
