# Per-scan peak blocks: (bin delta, intensity) UInt32 words -> four byte planes -> zstd.
#
# Byte-identical to TimsSlices.jl's per-slice `.tdfs` blocks (docs/format.md there), so a reader that already
# decodes those (`TimsSlices.decode_slice!`) decodes these: 2·n words, per peak `bin − previous bin` (the
# accumulator starts at 0xFFFFFFFF, so the first delta is bin + 1) then the intensity; plane p holds byte p of
# every word; one zstd frame per scan.

using Zstd_jll: libzstd

mutable struct ZstdCtx
    cctx::Ptr{Cvoid}
    dctx::Ptr{Cvoid}
    function ZstdCtx()
        c = ccall((:ZSTD_createCCtx, libzstd), Ptr{Cvoid}, ())
        d = ccall((:ZSTD_createDCtx, libzstd), Ptr{Cvoid}, ())
        (c == C_NULL || d == C_NULL) && error("ZSTD context allocation failed")
        z = new(c, d)
        finalizer(z) do z
            ccall((:ZSTD_freeCCtx, libzstd), Csize_t, (Ptr{Cvoid},), z.cctx)
            ccall((:ZSTD_freeDCtx, libzstd), Csize_t, (Ptr{Cvoid},), z.dctx)
        end
        z
    end
end

_zstd_check(r::Csize_t, what) = ccall((:ZSTD_isError, libzstd), Cuint, (Csize_t,), r) != 0 &&
    error("zstd $what: ", unsafe_string(ccall((:ZSTD_getErrorName, libzstd), Cstring, (Csize_t,), r)))

"Per-thread scratch for encoding / decoding blocks."
mutable struct BlockCodec
    words::Vector{UInt32}
    planes::Vector{UInt8}
    zbuf::Vector{UInt8}
    zstd::ZstdCtx
end
BlockCodec() = BlockCodec(UInt32[], UInt8[], UInt8[], ZstdCtx())

function _transpose!(planes::Vector{UInt8}, words::Vector{UInt32}, n::Int)
    length(planes) < 4n && resize!(planes, 4n)
    @inbounds for pl in 0:3
        sh = 8pl; base = pl * n
        @simd for i in 1:n
            planes[base+i] = (words[i] >> sh) % UInt8
        end
    end
    planes
end

function _untranspose!(words::Vector{UInt32}, planes::AbstractVector{UInt8}, n::Int)
    length(words) < n && resize!(words, n)
    @inbounds @simd for i in 1:n
        words[i] = UInt32(planes[i]) | UInt32(planes[n+i]) << 8 | UInt32(planes[2n+i]) << 16 | UInt32(planes[3n+i]) << 24
    end
    words
end

"""
    encode_block!(codec, bins, ints, n, level) -> nbytes

Encode peaks `bins[1:n]` (strictly increasing) / `ints[1:n]` into `codec.zbuf[1:nbytes]`. `n == 0` gives 0 bytes.
"""
function encode_block!(c::BlockCodec, bins::AbstractVector{UInt32}, ints::AbstractVector{UInt32}, n::Int, level::Integer)
    n == 0 && return 0
    w = c.words
    length(w) < 2n && resize!(w, 2n)
    acc = typemax(UInt32)
    @inbounds for k in 1:n
        w[2k-1] = bins[k] - acc; w[2k] = ints[k]; acc = bins[k]
    end
    _transpose!(c.planes, w, 2n)
    bound = Int(ccall((:ZSTD_compressBound, libzstd), Csize_t, (Csize_t,), 8n))
    length(c.zbuf) < bound && resize!(c.zbuf, bound)
    r = GC.@preserve c ccall((:ZSTD_compressCCtx, libzstd), Csize_t,
        (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint),
        c.zstd.cctx, pointer(c.zbuf), length(c.zbuf), pointer(c.planes), 8n, level)
    _zstd_check(r, "compress")
    Int(r)
end

"""
    decode_block!(bins, ints, codec, payload, n) -> (bins, ints)

Inverse of `encode_block!`: fills `bins[1:n]`, `ints[1:n]` (grown if needed) from one scan's zstd block.
"""
function decode_block!(bins::Vector{UInt32}, ints::Vector{UInt32}, c::BlockCodec, payload::AbstractVector{UInt8}, n::Int)
    n == 0 && return bins, ints
    length(c.planes) < 8n && resize!(c.planes, 8n)
    r = GC.@preserve c payload ccall((:ZSTD_decompressDCtx, libzstd), Csize_t,
        (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t),
        c.zstd.dctx, pointer(c.planes), length(c.planes), pointer(payload), length(payload))
    _zstd_check(r, "decompress")
    Int(r) == 8n || error("block decompressed to $(Int(r)) bytes, expected $(8n)")
    _untranspose!(c.words, c.planes, 2n)
    length(bins) < n && (resize!(bins, n); resize!(ints, n))
    acc = typemax(UInt32); w = c.words
    @inbounds for k in 1:n
        acc += w[2k-1]; bins[k] = acc; ints[k] = w[2k]
    end
    bins, ints
end
