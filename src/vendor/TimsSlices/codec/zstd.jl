# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
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
