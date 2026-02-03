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

# Build parameter default functions are loaded from buildParamDefaults.jl

struct InvalidParametersError <: Exception
    message::String
    params::Dict{String, Any}
end

function check_params_bsp(json_string::String)
    # Parse user parameters
    user_params = JSON.parse(json_string)
    
    # Always use full defaults - user values will be merged over them
    defaults = get_build_default_parameters()

    # Merge user params over defaults
    params = merge_with_build_defaults(user_params, defaults)

    # Backward compatibility: Move calibration_raw_file from library_params to top level if needed
    if !haskey(params, "calibration_raw_file") && haskey(params, "library_params") && haskey(params["library_params"], "calibration_raw_file")
        params["calibration_raw_file"] = params["library_params"]["calibration_raw_file"]
        delete!(params["library_params"], "calibration_raw_file")
    end

    # Helper function to check if a key exists and has the correct type
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end

    # Check required top-level parameters
    required_params = ["fasta_digest_params", "library_params", "variable_mods", "fixed_mods", "isotope_mod_groups", "max_koina_requests", "max_koina_batch", "match_lib_build_batch", "fasta_paths", "fasta_names", "library_path", "predict_fragments", "include_contaminants"]
    for param in required_params
        check_param(params, param, param in ["predict_fragments", "include_contaminants"] ? Bool : Any)
    end

    # Check if calibration_raw_file is provided (optional but recommended)
    if !haskey(params, "calibration_raw_file") || isempty(params["calibration_raw_file"])
        @user_warn "No calibration_raw_file provided. Fragment m/z bounds will not be auto-detected from raw data."
        params["calibration_raw_file"] = ""  # Set to empty string for consistency
    end
    
    # nce_params defaults now come from JSON file

    # Check fasta_digest_params
    fasta_digest_params = params["fasta_digest_params"]
    check_param(fasta_digest_params, "min_length", Integer)
    check_param(fasta_digest_params, "max_length", Integer)
    check_param(fasta_digest_params, "min_charge", Integer)
    check_param(fasta_digest_params, "max_charge", Integer)
    check_param(fasta_digest_params, "cleavage_regex", String)
    check_param(fasta_digest_params, "missed_cleavages", Integer)
    check_param(fasta_digest_params, "max_var_mods", Integer)
    check_param(fasta_digest_params, "add_decoys", Bool)
    check_param(fasta_digest_params, "entrapment_r", Real)

    if haskey(fasta_digest_params, "specificity")
        check_param(fasta_digest_params, "specificity", String)
    end
    if !haskey(fasta_digest_params, "specificity") && haskey(fasta_digest_params, "enzymaticity")
        fasta_digest_params["specificity"] = fasta_digest_params["enzymaticity"]
    end

    if !haskey(fasta_digest_params, "specificity")
        fasta_digest_params["specificity"] = "full"
    else
        specificity = lowercase(fasta_digest_params["specificity"])
        if !(specificity in ["full", "semi-n", "semi-c", "semi"])
            error("specificity must be one of \"full\", \"semi-n\", \"semi-c\", or \"semi\"; got: $(fasta_digest_params["specificity"])")
        end
        fasta_digest_params["specificity"] = specificity
    end
    
    # Check decoy_method with default value
    if !haskey(fasta_digest_params, "decoy_method")
        fasta_digest_params["decoy_method"] = "shuffle"
    else
        decoy_method = fasta_digest_params["decoy_method"]
        if !(decoy_method in ["shuffle", "reverse"])
            error("decoy_method must be either 'shuffle' or 'reverse', got: $decoy_method")
        end
    end
    
    # Check entrapment_method with default value and validation rules
    if !haskey(fasta_digest_params, "entrapment_method")
        fasta_digest_params["entrapment_method"] = "shuffle"
    else
        entrapment_method = fasta_digest_params["entrapment_method"]
        decoy_method = fasta_digest_params["decoy_method"]
        
        if !(entrapment_method in ["shuffle", "reverse"])
            error("entrapment_method must be either 'shuffle' or 'reverse', got: $entrapment_method")
        end
        
        # Validate the combination rules
        if decoy_method == "reverse" && entrapment_method == "reverse"
            error("When decoy_method is 'reverse', entrapment_method must be 'shuffle'")
        end
    end

    # Check nce_params
    nce_params = params["nce_params"]
    check_param(nce_params, "nce", Real)
    check_param(nce_params, "default_charge", Integer)
    check_param(nce_params, "dynamic_nce", Bool)

    # Check library_params (defaults now come from JSON file)
    library_params = params["library_params"]
    check_param(library_params, "rt_bin_tol", Real)
    check_param(library_params, "frag_bin_tol_ppm", Real)
    check_param(library_params, "rank_to_score", Vector)
    check_param(library_params, "y_start_index", Integer)
    check_param(library_params, "b_start_index", Integer)
    check_param(library_params, "y_start", Integer)
    check_param(library_params, "b_start", Integer)
    check_param(library_params, "include_p_index", Bool)
    check_param(library_params, "include_p", Bool)
    check_param(library_params, "auto_detect_frag_bounds", Bool)
    # calibration_raw_file is optional - only check if it exists
    # check_param(library_params, "calibration_raw_file", String)  # Commented out - optional parameter
    check_param(library_params, "frag_mz_min", Real)
    check_param(library_params, "frag_mz_max", Real)
    check_param(library_params, "prec_mz_min", Real)
    check_param(library_params, "prec_mz_max", Real)
    check_param(library_params, "max_frag_charge", Integer)
    check_param(library_params, "max_frag_rank", Integer)
    check_param(library_params, "length_to_frag_count_multiple", Real)
    check_param(library_params, "min_frag_intensity", Real)
    check_param(library_params, "include_isotope", Bool)
    check_param(library_params, "include_internal", Bool)
    check_param(library_params, "include_immonium", Bool)
    check_param(library_params, "include_neutral_diff", Bool)
    check_param(library_params, "instrument_type", String)
    check_param(library_params, "prediction_model", String)

    # Check variable_mods and fixed_mods
    for mod_type in ["variable_mods", "fixed_mods"]
        mods = params[mod_type]
        check_param(mods, "pattern", Vector)
        check_param(mods, "mass", Vector)
        check_param(mods, "name", Vector)
    end

    # Check isotope_mod_groups
    isotope_mod_groups = params["isotope_mod_groups"]
    if !(isotope_mod_groups isa Vector)
        throw(InvalidParametersError("isotope_mod_groups should be an array", params))
    end
    for group in isotope_mod_groups
        check_param(group, "pattern", Vector)
        check_param(group, "mass", Vector)
        check_param(group, "mod_names", Vector)
        check_param(group, "label_name", String)
    end

    # expand any home directories "~"
    params["library_path"] = expanduser(params["library_path"])
    if !isempty(params["calibration_raw_file"])
        params["calibration_raw_file"] = expanduser(params["calibration_raw_file"])
    end

    # Handle .poin extension - add if not present, don't duplicate if already there
    if !endswith(params["library_path"], ".poin")
        params["_lib_dir"] = params["library_path"] * ".poin"
        params["_lib_name"] = basename(params["library_path"])
    else
        params["_lib_dir"] = params["library_path"]
        params["_lib_name"] = basename(params["library_path"][1:end-5])  # Strip .poin
    end

    for i in range(1,length(params["fasta_paths"]))
        params["fasta_paths"][i] = expanduser(params["fasta_paths"][i])
    end

    # Validate per-FASTA regex arrays if provided
    regex_keys = [
        "fasta_header_regex_accessions",
        "fasta_header_regex_genes",
        "fasta_header_regex_proteins",
        "fasta_header_regex_organisms",
    ]
    for key in regex_keys
        if haskey(params, key)
            if !(params[key] isa Vector)
                throw(InvalidParametersError("$key must be an array", params))
            end
            if length(params[key]) != length(params["fasta_paths"])
                throw(InvalidParametersError("Length of $key must match fasta_paths", params))
            end
        end
    end

    # Add default contaminants fasta
    if params["include_contaminants"]
        contam_path = asset_path("contaminants.fasta.gz")
        push!(params["fasta_paths"], contam_path)
        push!(params["fasta_names"], "CONTAM")
        push!(get!(params, "fasta_header_regex_accessions", String[]), "^\\w+\\|(\\w+(?:-\\d+)?)\\|")
        push!(get!(params, "fasta_header_regex_genes", String[]), " GN=(\\S+)")
        push!(get!(params, "fasta_header_regex_proteins", String[]), "^\\w+\\|\\w+?:-\\d+?\\|[^ ]+ (.*?) [^ ]+=")
        push!(get!(params, "fasta_header_regex_organisms", String[]), " OS=([^ ]+.*?) [^ ]+=")
    end
    
    # If all checks pass, return the validated parameters
    return params
end


"""
    checkParseSpecLibParams(json_path::String)

Validate the parameters for ParseSpecLib from a JSON configuration file.

Parameters:
- json_path: Path to JSON configuration file

Returns:
- Parsed and validated parameters dictionary
"""
function checkParseSpecLibParams(json_path::String)
    params = JSON.parsefile(json_path)
    
    # Helper function to check if a key exists and has the correct type
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end
    
    # Check required sections
    required_sections = ["library_params", "fixed_mods"]
    for section in required_sections
        check_param(params, section, Dict)
    end

    # Check required sections
    required_sections = ["isotope_mod_groups"]
    for section in required_sections
        check_param(params, section, Vector{Any})
    end
    
    # Validate library parameters
    lib_params = params["library_params"]
    required_lib_params = [
        ("input_lib_path", String),
        ("output_lib_path", String),
        ("rt_bin_tol", Number),
        ("frag_bin_tol_ppm", Number),
        ("rank_to_score", Vector),
        ("y_start_index", Number),
        ("b_start_index", Number),
        ("y_start", Number),
        ("b_start", Number),
        ("include_p_index", Bool),
        ("include_p", Bool),
        ("include_isotope", Bool),
        ("include_immonium", Bool),
        ("include_internal", Bool),
        ("include_neutral_diff", Bool),
        ("max_frag_charge", Number),
        ("max_frag_rank", Number),
        ("length_to_frag_count_multiple", Number),
        ("min_frag_intensity", Number),
        ("generate_decoys", Bool),
        ("generate_entrapment", Bool),
        ("entrapment_groups", Number),
        ("instrument_type", String)
    ]
    
    for (param, type) in required_lib_params
        check_param(lib_params, param, type)
    end
    
    # Validate fixed mods
    fixed_mods = params["fixed_mods"]
    for key in ["mass", "name"]
        check_param(fixed_mods, key, Vector)
    end
    
    # Validate isotope mod groups
    iso_groups = params["isotope_mod_groups"]
    for group in iso_groups
        check_param(group, "name", String)
        check_param(group, "channels", Vector)
        
        for channel in group["channels"]
            check_param(channel, "channel", String)
            check_param(channel, "mass", Number)
        end
    end
    
    # Optional: validate sulfur mod groups if present
    if haskey(params, "sulfur_mod_groups")
        for group in params["sulfur_mod_groups"]
            check_param(group, "name", String)
            check_param(group, "sulfur_count", Number)
        end
    end

    # expand any home directories "~"
    params["library_params"]["input_lib_path"] = expanduser(params["library_params"]["input_lib_path"])
    params["library_params"]["output_lib_path"] = expanduser(params["library_params"]["output_lib_path"])
    
    return params
end
