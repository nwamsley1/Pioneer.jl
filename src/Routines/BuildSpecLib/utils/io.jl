# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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
        @user_warn "Loading legacy JLD2 format. Consider rebuilding library."
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
