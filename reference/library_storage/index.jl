"""Disk-only partition. Search continues to use Pioneer's native partition type."""
struct StoredPartition{I<:Unsigned,D<:Unsigned}
    lows::Vector{Float32}
    high_deltas::Vector{D}
    ends::Vector{UInt32}
    padding::UInt8
    rt_bins::Vector{P.FragIndexBin{Float32}}
    precursor_ids::Vector{I}
    rank_masks::Vector{UInt8}
    local_to_global::Vector{UInt32}
    n_local_precursors::I
    skip_hints::Vector{UInt16}
end

struct StoredIndex{I<:Unsigned}
    partitions::Vector{StoredPartition{I}}
    partition_bounds::Vector{Tuple{Float32,Float32}}
end

function encode_partition(partition, ::Type{I}) where {I<:Unsigned}
    bins = partition.fragment_bins
    n = length(bins.lows)
    require(length(bins.first_bins) == length(bins.last_bins) == n,
        "Fragment-bin columns must have equal lengths")
    require(length(bins.highs) >= n, "Missing high bounds")
    require(all(x -> isequal(x, Inf32), @view bins.highs[n + 1:end]), "SIMD padding must contain positive infinities")
    require(length(bins.highs) - n <= typemax(UInt8), "Too much SIMD padding")
    deltas = Vector{UInt32}(undef, n)
    for i in 1:n
        low, high = bins.lows[i], bins.highs[i]
        require(isfinite(low) && isfinite(high) && 0 < low <= high,
            "Fragment bounds must be positive, finite, and ordered")
        first = i == 1 ? UInt32(1) : bins.last_bins[i - 1] + UInt32(1)
        require(bins.first_bins[i] == first, "Fragment-bin entry ranges must be contiguous")
        deltas[i] = reinterpret(UInt32, high) - reinterpret(UInt32, low)
    end
    D = unsigned_type(maximum(deltas; init = UInt32(0)))
    return StoredPartition(
        bins.lows, D.(deltas), bins.last_bins, UInt8(length(bins.highs) - n),
        partition.rt_bins, I[P.getPrecID(f) for f in partition.fragments],
        UInt8[P.getScore(f) for f in partition.fragments], partition.local_to_global,
        partition.n_local_precs, partition.skip_hints
    )
end

"""
    encode_index(index) -> StoredIndex

Separate IDs and full eight-bit rank masks, omit derivable range starts, and
store exact high-bound bit deltas. Delta width is selected per partition.
"""
function encode_index(index::P.AbstractLocalPartitionedFragmentIndex{Float32})
    I = P.local_id_type(index)
    partitions = StoredPartition{I}[encode_partition(p, I) for p in index.partitions]
    return StoredIndex{I}(partitions, index.partition_bounds)
end

function decode_partition(stored::StoredPartition{I}) where {I}
    n = length(stored.lows)
    require(length(stored.high_deltas) == length(stored.ends) == n,
        "Stored fragment-bin columns must have equal lengths")
    require(length(stored.precursor_ids) == length(stored.rank_masks),
        "Stored ID and mask columns must have equal lengths")
    highs = fill(Inf32, n + Int(stored.padding))
    firsts = Vector{UInt32}(undef, n)
    for i in 1:n
        highs[i] = reinterpret(Float32,
            reinterpret(UInt32, stored.lows[i]) + UInt32(stored.high_deltas[i]))
        firsts[i] = i == 1 ? UInt32(1) : stored.ends[i - 1] + UInt32(1)
    end
    F = P.local_fragment_type(I)
    entries = F[F(id, mask) for (id, mask) in zip(stored.precursor_ids, stored.rank_masks)]
    Q = P.local_partition_type(I){Float32}
    return Q(
        P.SoAFragBins(stored.lows, highs, firsts, stored.ends), stored.rt_bins,
        entries, stored.local_to_global, stored.n_local_precursors, stored.skip_hints
    )
end

"""Reconstruct the native index arrays once; add no indirection during search."""
function decode_index(stored::StoredIndex{I}) where {I}
    Q = P.local_partition_type(I){Float32}
    partitions = Q[decode_partition(p) for p in stored.partitions]
    return P.local_index_type(I){Float32}(
        partitions, stored.partition_bounds, length(partitions)
    )
end

write_index(path, index) = P.serialize_to_jls(path, encode_index(index))
read_index(path) = decode_index(P.deserialize_from_jls(path))
