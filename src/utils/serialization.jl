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

const DETAILED_FRAGMENTS_FORMAT_TAG = "pioneer_detailed_fragments_v1"

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

function serialize_library_detailed_frags(filepath::String, frags; level::Int=6)
    frag_type = eltype(frags)
    frag_kind = String(nameof(frag_type))

    if frag_kind ∉ ("DetailedFrag", "SplineDetailedFrag")
        throw(ArgumentError("Unsupported detailed fragment type $(frag_type) for library serialization"))
    end

    payload = Dict{String, Any}(
        "__pioneer_format__" => DETAILED_FRAGMENTS_FORMAT_TAG,
        "kind" => frag_kind,
        "prec_id" => getfield.(frags, :prec_id),
        "mz" => getfield.(frags, :mz),
        "intensity" => getfield.(frags, :intensity),
        "ion_type" => getfield.(frags, :ion_type),
        "is_y" => getfield.(frags, :is_y),
        "is_b" => getfield.(frags, :is_b),
        "is_p" => getfield.(frags, :is_p),
        "is_isotope" => getfield.(frags, :is_isotope),
        "frag_charge" => getfield.(frags, :frag_charge),
        "ion_position" => getfield.(frags, :ion_position),
        "prec_charge" => getfield.(frags, :prec_charge),
        "rank" => getfield.(frags, :rank),
        "sulfur_count" => getfield.(frags, :sulfur_count),
    )

    serialize_to_jls(filepath, payload; level=level)
    return nothing
end

function deserialize_library_detailed_frags(filepath::String)
    payload = deserialize_from_jls(filepath)

    if !(payload isa AbstractDict) || get(payload, "__pioneer_format__", nothing) != DETAILED_FRAGMENTS_FORMAT_TAG
        return payload
    end

    n = length(payload["prec_id"])
    kind = payload["kind"]

    if kind == "DetailedFrag"
        return [
            DetailedFrag(
                payload["prec_id"][i],
                payload["mz"][i],
                payload["intensity"][i],
                payload["ion_type"][i],
                payload["is_y"][i],
                payload["is_b"][i],
                payload["is_p"][i],
                payload["is_isotope"][i],
                payload["frag_charge"][i],
                payload["ion_position"][i],
                payload["prec_charge"][i],
                payload["rank"][i],
                payload["sulfur_count"][i],
            ) for i in 1:n
        ]
    elseif kind == "SplineDetailedFrag"
        return [
            SplineDetailedFrag(
                payload["prec_id"][i],
                payload["mz"][i],
                payload["intensity"][i],
                payload["ion_type"][i],
                payload["is_y"][i],
                payload["is_b"][i],
                payload["is_p"][i],
                payload["is_isotope"][i],
                payload["frag_charge"][i],
                payload["ion_position"][i],
                payload["prec_charge"][i],
                payload["rank"][i],
                payload["sulfur_count"][i],
            ) for i in 1:n
        ]
    end

    throw(ArgumentError("Unsupported detailed fragment payload kind: $(kind)"))
end
