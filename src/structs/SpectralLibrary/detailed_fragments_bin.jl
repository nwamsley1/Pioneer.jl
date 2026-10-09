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

# detailed_fragments.bin: a library's fragments and per-precursor fragment ranges as raw arrays, read by
# memory-mapping (the search pages in the fragments it touches instead of loading the whole table).
#
# Layout (little-endian, the in-memory layout of the element type):
#   header (64 bytes): magic "PIONDFRG", version::UInt32, kind::UInt32, elsize::UInt32, n_coef::UInt32,
#                      n_frags::UInt64, n_ranges::UInt64, frags_offset::UInt64, ranges_offset::UInt64, 0::UInt64
#   frags  at frags_offset  (page-aligned): n_frags elements
#   ranges at ranges_offset (page-aligned): n_ranges UInt64 = prec_frag_ranges (precursor i: ranges[i]:ranges[i+1]-1)
# Padding bytes (between sections and inside padded element types) are zero, so a file's bytes depend only on its
# contents.

using Mmap: Mmap

const DETAILED_FRAGS_MAGIC = b"PIONDFRG"
const DETAILED_FRAGS_VERSION = UInt32(1)
const DETAILED_FRAGS_HEADER_BYTES = 64
const DETAILED_FRAGS_ALIGN = 4096

"Element-type code and spline coefficient count stored in the header."
detailed_frags_kind(::Type{SplineCompactFrag{N, Float32}}) where {N} = (UInt32(1), UInt32(N))
detailed_frags_kind(::Type{CompactFrag{Float32}}) = (UInt32(2), UInt32(0))
function detailed_frags_type(kind::UInt32, n_coef::UInt32)
    kind == 1 && return SplineCompactFrag{Int(n_coef), Float32}
    kind == 2 && return CompactFrag{Float32}
    error("detailed_fragments.bin: unknown fragment kind $kind")
end

_aligned(n::Integer) = cld(n, DETAILED_FRAGS_ALIGN) * DETAILED_FRAGS_ALIGN
"Whether a struct type has padding bytes (primitive types such as Float32 have no fields and no padding)."
_has_padding(::Type{T}) where {T} =
    fieldcount(T) > 0 && sizeof(T) != sum(i -> sizeof(fieldtype(T, i)), 1:fieldcount(T); init = 0)

"""
    DetailedFragsWriter{T}(path)

Streams fragments of type `T` to a `detailed_fragments.bin` at `path`: `append_frags!` per batch, then
`finish_detailed_frags!(w, ranges)` writes the ranges and the header.
"""
mutable struct DetailedFragsWriter{T}
    io::IOStream
    n_frags::Int
    scratch::Vector{UInt8}
end

function DetailedFragsWriter{T}(path::AbstractString) where {T}
    isbitstype(T) || error("detailed_fragments.bin needs an isbits fragment type, got $T")
    io = open(path, "w")
    write(io, zeros(UInt8, DETAILED_FRAGS_ALIGN))                  # header placeholder + padding to the frags
    return DetailedFragsWriter{T}(io, 0, UInt8[])
end

function append_frags!(w::DetailedFragsWriter{T}, frags::AbstractVector{T}) where {T}
    if _has_padding(T)
        _write_packed!(w.io, w.scratch, frags)
    elseif frags isa Vector{T}
        write(w.io, frags)
    else
        write(w.io, collect(frags))
    end
    w.n_frags += length(frags)
    return w
end

"Write padded elements field by field into a zeroed buffer (the padding bytes of a value are undefined)."
function _write_packed!(io::IOStream, buf::Vector{UInt8}, frags::AbstractVector{T}) where {T}
    chunk = 1 << 16
    for lo in 1:chunk:length(frags)
        hi = min(lo + chunk - 1, length(frags))
        n = (hi - lo + 1) * sizeof(T)
        resize!(buf, n); fill!(buf, 0x00)
        GC.@preserve buf begin
            p = pointer(buf)
            for (j, i) in enumerate(lo:hi)
                _store_fields!(p + (j - 1) * sizeof(T), frags[i])
            end
        end
        write(io, view(buf, 1:n))
    end
    return nothing
end
@generated function _store_fields!(p::Ptr{UInt8}, x::T) where {T}
    stores = [:(unsafe_store!(Ptr{$(fieldtype(T, i))}(p + $(Int(fieldoffset(T, i)))), getfield(x, $i))) for i in 1:fieldcount(T)]
    return Expr(:block, stores..., :(return nothing))
end

function finish_detailed_frags!(w::DetailedFragsWriter{T}, ranges::AbstractVector{UInt64}) where {T}
    io = w.io
    frags_end = DETAILED_FRAGS_ALIGN + w.n_frags * sizeof(T)
    ranges_offset = _aligned(frags_end)
    write(io, zeros(UInt8, ranges_offset - frags_end))
    write(io, ranges isa Vector{UInt64} ? ranges : collect(UInt64, ranges))
    seek(io, 0)
    kind, n_coef = detailed_frags_kind(T)
    write(io, DETAILED_FRAGS_MAGIC, DETAILED_FRAGS_VERSION, kind, UInt32(sizeof(T)), n_coef,
          UInt64(w.n_frags), UInt64(length(ranges)), UInt64(DETAILED_FRAGS_ALIGN), UInt64(ranges_offset), UInt64(0))
    close(io)
    return nothing
end

"Write `frags` and `ranges` (prec_frag_ranges) to a detailed_fragments.bin at `path`."
function write_detailed_frags(path::AbstractString, frags::AbstractVector{T}, ranges::AbstractVector{UInt64}) where {T}
    w = DetailedFragsWriter{T}(path)
    try
        append_frags!(w, frags)
    catch
        close(w.io); rethrow()
    end
    finish_detailed_frags!(w, ranges)
    return path
end

"""
    mmap_detailed_frags(path) -> (frags, ranges)

The fragments and per-precursor ranges of a detailed_fragments.bin as read-only memory-mapped `Vector`s.
"""
function mmap_detailed_frags(path::AbstractString)
    open(path, "r") do io
        magic = read(io, 8)
        magic == DETAILED_FRAGS_MAGIC || error("$path is not a detailed_fragments.bin file")
        version = read(io, UInt32)
        version == DETAILED_FRAGS_VERSION || error("$path: format version $version, this Pioneer reads $DETAILED_FRAGS_VERSION")
        kind = read(io, UInt32); elsize = read(io, UInt32); n_coef = read(io, UInt32)
        T = detailed_frags_type(kind, n_coef)
        elsize == sizeof(T) || error("$path: element size $elsize, expected $(sizeof(T)) for $T")
        n_frags = read(io, UInt64); n_ranges = read(io, UInt64)
        frags_offset = read(io, UInt64); ranges_offset = read(io, UInt64)
        return _mmap_sections(io, T, n_frags, n_ranges, frags_offset, ranges_offset)
    end
end
function _mmap_sections(io::IOStream, ::Type{T}, n_frags::UInt64, n_ranges::UInt64,
                        frags_offset::UInt64, ranges_offset::UInt64) where {T}
    frags = Mmap.mmap(io, Vector{T}, Int(n_frags), Int64(frags_offset); grow = false)   # read-only: io opened "r"
    ranges = Mmap.mmap(io, Vector{UInt64}, Int(n_ranges), Int64(ranges_offset); grow = false)
    return frags, ranges
end

"""
    load_detailed_frags_and_ranges(lib_dir) -> (frags, prec_frag_ranges)

A library's fragments and per-precursor ranges: memory-mapped from detailed_fragments.bin, or loaded from the legacy
detailed_fragments.jls / .jld2 + precursor_to_fragment_indices.jls / .jld2 of older libraries.
"""
function load_detailed_frags_and_ranges(lib_dir::AbstractString)
    bin = joinpath(lib_dir, "detailed_fragments.bin")
    isfile(bin) && return mmap_detailed_frags(bin)
    frag_path = joinpath(lib_dir, "detailed_fragments")
    frags = if isfile(frag_path * ".jls")
        deserialize_from_jls(frag_path * ".jls")
    elseif isfile(frag_path * ".jld2")
        @user_warn "Loading legacy JLD2 format for detailed_fragments. Consider rebuilding library."
        jldopen(frag_path * ".jld2", "r") do file
            read(file, "data")
        end
    else
        error("Fragment file not found: $(frag_path).bin, .jls or .jld2")
    end
    ranges_path = joinpath(lib_dir, "precursor_to_fragment_indices")
    ranges = if isfile(ranges_path * ".jls")
        deserialize_from_jls(ranges_path * ".jls")
    elseif isfile(ranges_path * ".jld2")
        @user_warn "Loading legacy JLD2 format for precursor_to_fragment_indices. Consider rebuilding library."
        load(ranges_path * ".jld2")["pid_to_fid"]
    else
        error("precursor_to_fragment_indices file not found in $lib_dir")
    end
    return frags, ranges
end
