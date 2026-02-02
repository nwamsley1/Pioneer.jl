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
Save detailed fragment information to file using Julia Serialization with Zlib compression.
"""
function save_detailed_frags(filepath::String, fragments::Vector{T}) where T <: AbstractKoinaFragment
    # Ensure .jls extension
    if !endswith(filepath, ".jls")
        filepath = replace(filepath, ".jld2" => ".jls")
        if !endswith(filepath, ".jls")
            filepath = filepath * ".jls"
        end
    end

    # Ensure directory exists
    mkpath(dirname(filepath))

    serialize_to_jls(filepath, fragments)
end

"""
Read fragments from a saved file. Supports both new .jls and legacy .jld2 formats.
"""
function read_detailed_frags(filepath::String)
    # Try .jls first
    jls_path = replace(filepath, ".jld2" => ".jls")
    if isfile(jls_path)
        return deserialize_from_jls(jls_path)
    end

    # Fall back to legacy .jld2
    if isfile(filepath) && endswith(filepath, ".jld2")
        @warn "Loading legacy JLD2 format. Consider rebuilding library."
        data = load(filepath)
        # Handle both key names used historically
        return haskey(data, "fragments") ? data["fragments"] : data["data"]
    end

    error("Fragment file not found: $jls_path or $filepath")
end

"""
Save/append fragment data to Arrow format.
"""
function save_fragments_arrow(filepath::String, fragments::DataFrame; mode::Symbol=:write)
    if mode == :write
        Arrow.write(filepath, fragments)
    elseif mode == :append
        Arrow.append(filepath, fragments)
    else
        throw(ArgumentError("Invalid mode: $mode. Must be :write or :append"))
    end
end