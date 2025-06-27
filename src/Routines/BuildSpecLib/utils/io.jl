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
Save detailed fragment information to file.
"""
function save_detailed_frags(filepath::String, fragments::Vector{T}) where T <: AbstractKoinaFragment
    if !endswith(filepath, ".jld2")
        throw(ArgumentError("Output file must have .jld2 extension"))
    end
    
    # Ensure directory exists
    mkpath(dirname(filepath))
    
    jldsave(filepath; fragments)
end

"""
Read fragments from a saved file.
"""
function read_detailed_frags(filepath::String)::Vector{AbstractKoinaFragment}
    if !isfile(filepath)
        error("Fragment file not found: $filepath")
    end
    
    data = load(filepath)
    return data["fragments"]
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