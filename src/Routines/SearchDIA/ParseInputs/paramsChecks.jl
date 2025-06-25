using JSON
#see src/Routines/BuildspecKib/utils
#struct InvalidParametersError <: Exception
#    message::String
#    params::Dict{String, Any}
#end

function checkParams(json_path::String)
    params = JSON.parsefile(json_path)
    # Helper function to check if a key exists and has the correct type
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end

    # Check top-level sections
    required_sections = [
        "global", "parameter_tuning", "first_search", "quant_search",
        "acquisition", "rt_alignment", "optimization", "output", "paths"
    ]
    
    for section in required_sections
        check_param(params, section, Dict)
    end

    # Validate global parameters
    global_params = params["global"]
    check_param(global_params, "isotope_settings", Dict)
    check_param(global_params["isotope_settings"], "err_bounds_first_pass", Vector)
    check_param(global_params["isotope_settings"], "err_bounds_quant_search", Vector)
    check_param(global_params["isotope_settings"], "combine_traces", Bool)
    check_param(global_params["isotope_settings"], "partial_capture", Bool)
    check_param(global_params["isotope_settings"], "min_fraction_transmitted", Real)
    check_param(global_params, "scoring", Dict)
    check_param(global_params["scoring"], "q_value_threshold", Real)

    check_param(global_params["huber_override"], "override_huber_delta_fit", Bool)
    check_param(global_params["huber_override"], "huber_delta", Real)
    check_param(global_params, "ms1_scoring", Bool)
    check_param(global_params, "ms1_quant", Bool)
    check_param(global_params, "match_between_runs", Bool)

    # Validate parameter tuning parameters
    tuning_params = params["parameter_tuning"]
    check_param(tuning_params, "fragment_settings", Dict)
    frag_settings = tuning_params["fragment_settings"]
    check_param(frag_settings, "min_count", Integer)
    check_param(frag_settings, "max_rank", Integer)
    check_param(frag_settings, "tol_ppm", Real)
    check_param(frag_settings, "min_score", Integer)
    check_param(frag_settings, "min_spectral_contrast", Real)
    check_param(frag_settings, "relative_improvement_threshold", Real)
    check_param(frag_settings, "min_log2_ratio", Real)
    check_param(frag_settings, "min_top_n", Vector)
    check_param(frag_settings, "n_isotopes", Integer)

    check_param(tuning_params, "search_settings", Dict)
    search_settings = tuning_params["search_settings"]
    check_param(search_settings, "sample_rate", Real)
    check_param(search_settings, "min_samples", Integer)
    check_param(search_settings, "max_presearch_iters", Integer)
    check_param(search_settings, "frag_err_quantile", Real)

    # Validate first search parameters
    first_search = params["first_search"]
    check_param(first_search, "fragment_settings", Dict)
    first_frag = first_search["fragment_settings"]
    check_param(first_frag, "min_count", Integer)
    check_param(first_frag, "max_rank", Integer)
    check_param(first_frag, "min_score", Integer)
    check_param(first_frag, "min_spectral_contrast", Real)
    check_param(first_frag, "relative_improvement_threshold", Real)
    check_param(first_frag, "min_log2_ratio", Real)
    check_param(first_frag, "min_top_n", Vector)
    check_param(first_frag, "n_isotopes", Integer)

    check_param(first_search, "scoring_settings", Dict)
    score_settings = first_search["scoring_settings"]
    check_param(score_settings, "n_train_rounds", Integer)
    check_param(score_settings, "max_iterations", Integer)
    check_param(score_settings, "max_q_value_probit_rescore", Real)
    check_param(score_settings, "max_local_fdr", Real)

    score_settings = first_search["irt_mapping"]
    check_param(score_settings, "max_prob_to_impute_irt", Real)
    check_param(score_settings, "fwhm_nstd", Real)
    check_param(score_settings, "irt_nstd", Real)
    # Validate quant search parameters
    quant_search = params["quant_search"]
    check_param(quant_search, "fragment_settings", Dict)
    check_param(quant_search, "chromatogram", Dict)
    
    quant_frag = quant_search["fragment_settings"]
    check_param(quant_frag, "min_count", Integer)
    check_param(quant_frag, "min_y_count", Integer)
    check_param(quant_frag, "max_rank", Integer)
    check_param(quant_frag, "min_spectral_contrast", Real)
    check_param(quant_frag, "min_log2_ratio", Real)
    check_param(quant_frag, "min_top_n", Vector)
    check_param(quant_frag, "n_isotopes", Integer)

    chrom_settings = quant_search["chromatogram"]
    check_param(chrom_settings, "smoothing_strength", Real)
    check_param(chrom_settings, "padding", Integer)
    check_param(chrom_settings, "max_apex_offset", Integer)

    # Validate acquisition parameters
    acq_params = params["acquisition"]
    check_param(acq_params, "nce", Integer)
    check_param(acq_params, "quad_transmission", Dict)
    quad_trans = acq_params["quad_transmission"]
    check_param(quad_trans, "fit_from_data", Bool)
    check_param(quad_trans, "overhang", Real)
    check_param(quad_trans, "smoothness", Real)

    # Validate RT alignment parameters
    rt_params = params["rt_alignment"]
    check_param(rt_params, "n_bins", Integer)
    check_param(rt_params, "bandwidth", Real)
    check_param(rt_params, "sigma_tolerance", Integer)
    check_param(rt_params, "min_probability", Real)

    # Validate optimization parameters
    opt_params = params["optimization"]
    check_param(opt_params, "deconvolution", Dict)
    check_param(opt_params, "machine_learning", Dict)

    deconv = opt_params["deconvolution"]
    check_param(deconv, "lambda", Real)
    check_param(deconv, "reg_type", String)
    check_param(deconv, "huber_delta", Real)
    check_param(deconv, "huber_exp", Real)
    check_param(deconv, "huber_iters", Integer)
    check_param(deconv, "newton_iters", Integer)
    check_param(deconv, "newton_accuracy", Real)
    check_param(deconv, "max_diff", Real)

    ml_params = opt_params["machine_learning"]
    check_param(ml_params, "max_psms_in_memory", Integer)
    check_param(ml_params, "min_trace_prob", Real)
    check_param(ml_params, "max_q_value_xgboost_rescore", Real)
    check_param(ml_params, "max_q_value_xgboost_mbr_rescore", Real)
    check_param(ml_params, "spline_points", Integer)
    check_param(ml_params, "interpolation_points", Integer)

    # Validate Protein Inference parameters
    output = params["proteinInference"]
    check_param(output, "min_peptides", Integer)

    # Validate MaxLFQ parameters
    output = params["maxLFQ"]
    check_param(output, "run_to_run_normalization", Bool)

    # Validate output parameters
    output = params["output"]
    check_param(output, "write_csv", Bool)
    check_param(output, "delete_temp", Bool)
    check_param(output, "write_decoys", Bool)
    check_param(output, "plots_per_page", Integer)

    # Validate path parameters
    paths = params["paths"]
    check_param(paths, "library", String)
    check_param(paths, "ms_data", String)
    check_param(paths, "results", String)

    return params
end