# ©2026 Washington University
#
# This file is part of TimsSlices.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Thin wrapper over libzstd with per-thread contexts and preallocated buffers (no allocation per call).
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

zstd_is_error(code::Csize_t) = ccall((:ZSTD_isError, libzstd), Cuint, (Csize_t,), code) != 0
zstd_error_name(code::Csize_t) = unsafe_string(ccall((:ZSTD_getErrorName, libzstd), Cstring, (Csize_t,), code))
zstd_compress_bound(n::Integer) = Int(ccall((:ZSTD_compressBound, libzstd), Csize_t, (Csize_t,), n))

"""
    zstd_compress!(dst, ctx, src, n, level) -> nbytes

Compress `src[1:n]` into `dst` (grown if needed), returning the compressed length.
"""
function zstd_compress!(dst::Vector{UInt8}, ctx::ZstdCtx, src::Vector{UInt8}, n::Integer, level::Integer)
    bound = zstd_compress_bound(n)
    length(dst) < bound && resize!(dst, bound)
    r = GC.@preserve dst src ccall((:ZSTD_compressCCtx, libzstd), Csize_t,
        (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Cint),
        ctx.cctx, pointer(dst), length(dst), pointer(src), n, level)
    zstd_is_error(r) && error("zstd compress: ", zstd_error_name(r))
    Int(r)
end

"""
    zstd_decompress!(dst, ctx, src, expected) -> nbytes

Decompress the zstd frame in `src` (a byte vector or view) into `dst`, which is resized to at least
`expected` bytes. Returns the decompressed length (checked against `expected` when it is > 0).
"""
function zstd_decompress!(dst::Vector{UInt8}, ctx::ZstdCtx, src::AbstractVector{UInt8}, expected::Integer)
    cap = expected > 0 ? expected : zstd_frame_content_size(src)
    length(dst) < cap && resize!(dst, cap)
    r = GC.@preserve dst src ccall((:ZSTD_decompressDCtx, libzstd), Csize_t,
        (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t),
        ctx.dctx, pointer(dst), length(dst), pointer(src), length(src))
    zstd_is_error(r) && error("zstd decompress: ", zstd_error_name(r))
    expected > 0 && Int(r) != expected && error("zstd decompress: got $(Int(r)) bytes, expected $expected")
    Int(r)
end

function zstd_frame_content_size(src::AbstractVector{UInt8})
    r = GC.@preserve src ccall((:ZSTD_getFrameContentSize, libzstd), Culonglong, (Ptr{UInt8}, Csize_t), pointer(src), length(src))
    (r == typemax(Culonglong) || r == typemax(Culonglong) - 1) && error("zstd: frame content size unknown")
    Int(r)
end
