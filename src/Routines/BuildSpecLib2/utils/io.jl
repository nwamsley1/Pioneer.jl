"""
Save detailed fragment information to file.
"""
function save_detailed_frags(filepath::String, fragments::Vector{T}) where T <: AbstractFragment
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
function read_detailed_frags(filepath::String)::Vector{AbstractFragment}
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