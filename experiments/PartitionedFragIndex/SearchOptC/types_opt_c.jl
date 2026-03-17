#
# Fixed-width bin types for O(1) fragment lookup
#

const FIXED_BIN_WIDTH = 0.005f0
const INV_FIXED_BIN_WIDTH = 1.0f0 / FIXED_BIN_WIDTH

"""
    FixedBinRTData

Per-RT-bin data within a FixedBinPartition. Contains the frag bin pointers
and the fragment array for this RT bin.

`frag_ptrs` has length `n_mz_bins + 1`. `frag_ptrs[i]` is the index of the
first fragment in m/z bin `i`. `frag_ptrs[end]` is one past the last fragment.
Fragments in bin `i` are `fragments[frag_ptrs[i] : frag_ptrs[i+1]-1]`.
Empty bins have `frag_ptrs[i] == frag_ptrs[i+1]`.
"""
struct FixedBinRTData
    frag_ptrs::Vector{UInt16}        # length n_mz_bins + 1, CSR-style pointers (max 65535 frags per RT bin)
    fragments::Vector{LocalFragment}  # fragments for this RT bin, sorted by mz bin
    irt_lo::Float32
    irt_hi::Float32
end

"""
    FixedBinPartition

A partition using fixed-width m/z bins for O(1) fragment lookup.
Each RT bin stores a CSR-style pointer array mapping m/z bucket → fragment range.

Given peak m/z `x`:
  bucket = floor(Int, (x - mz_min) * inv_bin_width) + 1
  frag_range = frag_ptrs[bucket] : frag_ptrs[bucket+1] - 1
  score all fragments[frag_range]

No exponential search, no binary search. One multiply + two loads.
"""
struct FixedBinPartition
    rt_data::Vector{FixedBinRTData}   # one per RT bin
    local_to_global::Vector{UInt32}
    n_local_precs::UInt16
    mz_min::Float32                   # lower edge of first m/z bin
    n_mz_bins::Int32                  # number of m/z bins
end

"""
    FixedBinPartitionedIndex

Top-level partitioned index using fixed-width m/z bins.
"""
struct FixedBinPartitionedIndex
    partitions::Vector{FixedBinPartition}
    partition_bounds::Vector{Tuple{Float32, Float32}}  # prec_mz bounds per partition
    n_partitions::Int
end

getPartitions(pfi::FixedBinPartitionedIndex) = pfi.partitions
getPartition(pfi::FixedBinPartitionedIndex, k::Int) = pfi.partitions[k]
getNPartitions(pfi::FixedBinPartitionedIndex) = pfi.n_partitions

@inline function get_partition_range(pfi::FixedBinPartitionedIndex, query_min::Float32, query_max::Float32)
    bounds = pfi.partition_bounds
    n = pfi.n_partitions
    first_k = n + 1
    lo, hi = 1, n
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if bounds[mid][2] >= query_min
            first_k = mid
            hi = mid - 1
        else
            lo = mid + 1
        end
    end
    last_k = first_k - 1
    lo, hi = first_k, n
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if bounds[mid][1] <= query_max
            last_k = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return first_k, last_k
end

"""
Compute the m/z bin index for a given m/z value.
One multiply + one truncation. No division.
"""
@inline function mz_to_bin(mz::Float32, mz_min::Float32)
    return floor(Int32, (mz - mz_min) * INV_FIXED_BIN_WIDTH) + Int32(1)
end
