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

Uses Julia's base Serialization with Zlib compression (via CodecZlib).
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
    open(filepath, "r") do file_io
        stream = ZlibDecompressorStream(file_io)
        try
            return deserialize(stream)
        finally
            close(stream)
        end
    end
end
