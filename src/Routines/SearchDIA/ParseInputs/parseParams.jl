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
        paths
    )
end