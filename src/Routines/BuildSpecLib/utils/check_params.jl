struct InvalidParametersError <: Exception
    message::String
    params::Dict{String, Any}
end

function check_params_bsp(json_string::String)
    params = JSON.parse(json_string)
    
    # Helper function to check if a key exists and has the correct type
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end

    # Check top-level parameters
    top_level_params = ["fasta_digest_params", "nce_params", "library_params", "variable_mods", "fixed_mods", "isotope_mod_groups", "max_koina_requests", "max_koina_batch", "match_lib_build_batch", "fasta_paths", "fasta_names", "out_dir", "lib_name", "new_lib_name", "out_name", "predict_fragments"]
    for param in top_level_params
        check_param(params, param, param in ["predict_fragments"] ? Bool : Any)
    end

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

    # Check nce_params
    nce_params = params["nce_params"]
    check_param(nce_params, "nce", Real)
    check_param(nce_params, "default_charge", Integer)
    check_param(nce_params, "dynamic_nce", Bool)

    # Check library_params
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
    check_param(library_params, "calibration_raw_file", String)
    check_param(library_params, "frag_mz_min", Real)
    check_param(library_params, "frag_mz_max", Real)
    check_param(library_params, "prec_mz_min", Real)
    check_param(library_params, "prec_mz_max", Real)
    check_param(library_params, "max_frag_charge", Integer)
    check_param(library_params, "max_frag_rank", Integer)
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

    # If all checks pass, return the validated parameters
    return params
end