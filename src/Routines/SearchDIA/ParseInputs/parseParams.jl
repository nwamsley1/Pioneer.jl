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

# Parameter default functions are loaded from paramDefaults.jl

struct PioneerParameters
    global_settings::NamedTuple
    search::NamedTuple
    acquisition::NamedTuple
    optimization::NamedTuple
    protein_scoring::NamedTuple
    maxLFQ::NamedTuple
    output::NamedTuple
    logging::NamedTuple
    paths::NamedTuple
end

"""
    params_to_dict(params::PioneerParameters)

Convert PioneerParameters struct back to nested Dict for JSON serialization.
Used to write complete merged configuration (user + defaults) to output directory.

# Returns
Dict with all parameter sections, ready for JSON.json() conversion
"""
function params_to_dict(params::PioneerParameters)
    # Recursively convert NamedTuples to Dicts
    function namedtuple_to_dict(nt::NamedTuple)
        result = Dict{String, Any}()
        for (k, v) in pairs(nt)
            key_str = string(k)
            result[key_str] = if v isa NamedTuple
                namedtuple_to_dict(v)
            else
                v
            end
        end
        return result
    end

    return Dict{String, Any}(
        "global" => namedtuple_to_dict(params.global_settings),
        "search" => namedtuple_to_dict(params.search),
        "acquisition" => namedtuple_to_dict(params.acquisition),
        "optimization" => namedtuple_to_dict(params.optimization),
        "proteinScoring" => namedtuple_to_dict(params.protein_scoring),
        "maxLFQ" => namedtuple_to_dict(params.maxLFQ),
        "output" => namedtuple_to_dict(params.output),
        "logging" => namedtuple_to_dict(params.logging),
        "paths" => namedtuple_to_dict(params.paths)
    )
end

# Helper function to merge settings with global defaults
function merge_with_globals(specific_settings::NamedTuple, global_settings::NamedTuple)
    merged = Dict{Symbol,Any}()
    for (k, v) in pairs(global_settings)
        merged[k] = v
    end
    for (k, v) in pairs(specific_settings)
        merged[k] = v
    end
    return NamedTuple(merged)
end


function parse_pioneer_parameters(json_path::String; apply_defaults::Bool = true)
    # Read user JSON
    user_params = JSON.parsefile(json_path, dicttype=Dict{String,Any})
    
    # Apply defaults if requested
    if apply_defaults
        # Determine if we should use simplified defaults based on parameter presence
        is_simplified = !haskey(user_params, "search")
        
        # Get appropriate defaults
        defaults = get_default_parameters(is_simplified)
        
        # Merge user params over defaults
        params = merge_with_defaults(user_params, defaults)
    else
        params = user_params
    end
    
    # Convert nested dictionaries to NamedTuples
    function dict_to_namedtuple(d::AbstractDict)
        # Create pairs for the NamedTuple constructor
        symbol_pairs = []

        for (k, v) in pairs(d)
            symbol_key = Symbol(k)
            if v isa AbstractDict
                push!(symbol_pairs, symbol_key => dict_to_namedtuple(v))
            else
                push!(symbol_pairs, symbol_key => v)
            end
        end
        
        return NamedTuple(symbol_pairs)
    end

    function expand_user_paths(nt::NamedTuple)
        vals_expanded = map(nt) do val
            expanduser(val)
        end
        return vals_expanded
    end

    # Parse each section
    global_settings = dict_to_namedtuple(params["global"])
    search = dict_to_namedtuple(params["search"])
    acquisition = dict_to_namedtuple(params["acquisition"])
    optimization = dict_to_namedtuple(params["optimization"])
    protein_scoring = dict_to_namedtuple(params["proteinScoring"])
    maxLFQ = dict_to_namedtuple(params["maxLFQ"])
    output = dict_to_namedtuple(params["output"])

    # Parse logging section (defaults already applied if apply_defaults=true)
    logging = dict_to_namedtuple(params["logging"])

    paths = expand_user_paths(dict_to_namedtuple(params["paths"]))

    return PioneerParameters(
        global_settings,
        search,
        acquisition,
        optimization,
        protein_scoring,
        maxLFQ,
        output,
        logging,
        paths
    )
end