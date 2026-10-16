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

"""
Serialization utilities for Pioneer spectral library data.

Uses Julia's base Serialization with Zlib compression (written via CodecZlib, read via LibDeflate).
Files use the `.jls` extension.
"""

"""
    serialize_to_jls(filepath::String, data; level::Int=6)

Serialize data to a file with Zlib compression.

# Arguments
- `filepath`: Path to output file (should end with .jls)
- `data`: Julia object to serialize
- `level`: Compression level 1-9 (1=fastest, 6=balanced, 9=best compression)

# Example
```julia
serialize_to_jls("fragments.jls", fragment_vector)
serialize_to_jls("fragments.jls", fragment_vector; level=9)  # Max compression
```
"""
function serialize_to_jls(filepath::String, data; level::Int=6)
    mkpath(dirname(filepath))
    open(filepath, "w") do file_io
        stream = ZlibCompressorStream(file_io; level=level)
        try
            serialize(stream, data)
        finally
            close(stream)
        end
    end
    return nothing
end

"""
    deserialize_from_jls(filepath::String)

Deserialize data from a Zlib-compressed .jls file.

# Arguments
- `filepath`: Path to input file

# Returns
- The deserialized Julia object

# Example
```julia
fragments = deserialize_from_jls("fragments.jls")
```
"""
function deserialize_from_jls(filepath::String)
    # Inflate the whole file at once with libdeflate (~2.5x faster than streaming through zlib), then
    # deserialize from memory. zlib does not record the inflated size, so the output buffer is an
    # oversized anonymous mapping: pages never written are never committed, and the buffer is unmapped
    # as soon as the object is read.
    compressed = Mmap.mmap(filepath)
    n_out = max(4 * length(compressed), 1 << 20)
    while true
        inflated = Mmap.mmap(Vector{UInt8}, n_out)
        try
            n = LibDeflate.zlib_decompress!(LibDeflate.Decompressor(), inflated, compressed)
            if n isa LibDeflate.LibDeflateError
                n === LibDeflate.LibDeflateErrors.insufficient_output_space || error("Could not inflate $filepath: $n")
                n_out *= 2
                continue
            end
            return GC.@preserve inflated deserialize(IOBuffer(unsafe_wrap(Array, pointer(inflated), n)))
        finally
            Base.finalize(inflated)
        end
    end
end
