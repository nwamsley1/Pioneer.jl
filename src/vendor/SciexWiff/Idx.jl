# The `SampleSubtree/Sample{N}/Idx` stream: one 54-byte record per (cycle, experiment) slot.
# Layout and evidence: notes/format.md §2.

const IDX_HEADER = 32
const IDX_RECORD = 54
const SCAN_FILE_HEADER = 44   # Idx offsets are relative to this byte of `.wiff.scan`

"""
    ScanIndex

Columns of the `Idx` stream, one entry per record (empty records included).
`block_offset` is absolute in `.wiff.scan` (the `ffffffff` marker); `meta_start` is where the block's
protobuf metadata begins (the end of the previous non-empty block).
"""
struct ScanIndex
    block_offset::Vector{Int64}
    block_size::Vector{Int32}
    meta_start::Vector{Int64}
    rt_ms::Vector{Float64}
    tic::Vector{Float64}
    base_peak_intensity::Vector{Float64}
    base_peak_bin::Vector{Float64}
end

Base.length(ix::ScanIndex) = length(ix.block_size)

@inline _rd(::Type{T}, b::AbstractVector{UInt8}, o::Integer) where {T} =
    GC.@preserve b ltoh(unsafe_load(Ptr{T}(pointer(b, o + 1))))

function ScanIndex(idx::Vector{UInt8})
    n, rem = divrem(length(idx) - IDX_HEADER, IDX_RECORD)
    (length(idx) >= IDX_HEADER && rem == 0) ||
        throw(ArgumentError("Idx stream of $(length(idx)) bytes is not 32 + 54·n"))
    off = Vector{Int64}(undef, n); sz = Vector{Int32}(undef, n); ms = Vector{Int64}(undef, n)
    rt = Vector{Float64}(undef, n); tic = similar(rt); bpi = similar(rt); bpb = similar(rt)
    prev_end = Int64(SCAN_FILE_HEADER)
    for r in 1:n
        o = IDX_HEADER + IDX_RECORD * (r - 1)
        # offsets are u32 and wrap in a .wiff.scan over 4 GiB; blocks are contiguous and in order, so a block
        # starts after the previous one ends: add 2^32 until it does
        off[r] = SCAN_FILE_HEADER + Int64(_rd(UInt32, idx, o))
        while off[r] < prev_end
            off[r] += Int64(1) << 32
        end
        sz[r] = Int32(_rd(UInt32, idx, o + 4))
        rt[r] = _rd(Float64, idx, o + 8)
        tic[r] = _rd(Float64, idx, o + 18)
        bpi[r] = _rd(Float64, idx, o + 26)
        bpb[r] = _rd(Float64, idx, o + 42)
        ms[r] = prev_end
        sz[r] > 0 && (prev_end = off[r] + sz[r])
    end
    ScanIndex(off, sz, ms, rt, tic, bpi, bpb)
end
