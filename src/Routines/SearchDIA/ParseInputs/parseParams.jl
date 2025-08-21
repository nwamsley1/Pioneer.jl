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

struct PioneerParameters
    global_settings::NamedTuple
    parameter_tuning::NamedTuple
    first_search::NamedTuple
    quant_search::NamedTuple
    acquisition::NamedTuple
    rt_alignment::NamedTuple
    optimization::NamedTuple
    protein_inference::NamedTuple
    maxLFQ::NamedTuple
    output::NamedTuple
    logging::NamedTuple
    paths::NamedTuple
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


function parse_pioneer_parameters(json_path::String)
    # Read and parse JSON
    params = JSON.parsefile(json_path)
    
    # Convert nested dictionaries to NamedTuples
    function dict_to_namedtuple(d::Dict)
        # Create pairs for the NamedTuple constructor
        symbol_pairs = []
        
        for (k, v) in pairs(d)
            symbol_key = Symbol(k)
            if v isa Dict
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
    parameter_tuning = dict_to_namedtuple(params["parameter_tuning"])
    first_search = dict_to_namedtuple(params["first_search"])
    quant_search = dict_to_namedtuple(params["quant_search"])
    acquisition = dict_to_namedtuple(params["acquisition"])
    rt_alignment = dict_to_namedtuple(params["rt_alignment"])
    optimization = dict_to_namedtuple(params["optimization"])
    protein_inference = dict_to_namedtuple(params["proteinInference"])
    maxLFQ = dict_to_namedtuple(params["maxLFQ"])
    output = dict_to_namedtuple(params["output"])
    
    # Parse logging section with defaults if not present
    logging = if haskey(params, "logging")
        dict_to_namedtuple(params["logging"])
    else
        # Default logging configuration
        NamedTuple((
            console_level = "normal",
            enable_progress = true,
            enable_warnings = true,
            rotation_size_mb = 100.0,
            max_warnings = 10000
        ))
    end
    
    paths = expand_user_paths(dict_to_namedtuple(params["paths"]))

    return PioneerParameters(
        global_settings,
        parameter_tuning,
        first_search,
        quant_search,
        acquisition,
        rt_alignment,
        optimization,
        protein_inference,
        maxLFQ,
        output,
        logging,
        paths
    )
end