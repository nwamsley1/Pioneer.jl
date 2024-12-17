using JSON

struct InvalidParametersError <: Exception
    message::String
    params::Dict{String, Any}
end

function checkParams(json_string::String)
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
    check_param(params, "isotope_err_bounds", Vector)
    check_param(params, "choose_most_intense", Bool)
    check_param(params, "q_value", Real)

    # Check nested parameter groups
    nested_params = [
        "presearch_params", "first_search_params", "summarize_first_search_params",
        "quant_search_params", "irt_mapping_params", "deconvolution_params",
        "qc_plot_params", "normalization_params", "xgboost_params",
        "quad_transmission", "benchmark_params", "output_params"
    ]

    for group in nested_params
        check_param(params, group, Dict)
    end

    # Validate presearch_params
    presearch = params["presearch_params"]
    check_param(presearch, "min_index_search_score", Real)
    check_param(presearch, "n_frag_isotopes", Integer)
    check_param(presearch, "min_frag_count", Integer)
    check_param(presearch, "min_log2_matched_ratio", Real)
    check_param(presearch, "min_spectral_contrast", Real)
    check_param(presearch, "min_topn_of_m", Vector)
    check_param(presearch, "max_best_rank", Integer)
    check_param(presearch, "sample_rate", Real)
    check_param(presearch, "frag_tol_ppm", Real)
    check_param(presearch, "max_qval", Real)
    check_param(presearch, "min_samples", Integer)
    check_param(presearch, "frag_err_quantile", Real)
    check_param(presearch, "max_presearch_iters", Integer)

    # Validate first_search_params
    first_search = params["first_search_params"]
    check_param(first_search, "min_index_search_score", Real)
    check_param(first_search, "min_frag_count", Integer)
    check_param(first_search, "min_topn_of_m", Vector)
    check_param(first_search, "n_frag_isotopes", Integer)
    check_param(first_search, "min_log2_matched_ratio", Real)
    check_param(first_search, "min_spectral_contrast", Real)
    check_param(first_search, "max_best_rank", Integer)
    check_param(first_search, "n_train_rounds_probit", Integer)
    check_param(first_search, "max_iter_probit", Integer)
    check_param(first_search, "max_q_value_probit_rescore", Real)
    check_param(first_search, "max_precursors_passing", Integer)

    # Validate summarize_first_search_params
    summarize = params["summarize_first_search_params"]
    check_param(summarize, "max_precursors", Integer)
    check_param(summarize, "scan_count_range", Vector)
    check_param(summarize, "max_q_val_for_irt", Real)
    check_param(summarize, "max_prob_to_impute", Real)
    check_param(summarize, "min_inference_points", Integer)
    check_param(summarize, "fwhm_nstd", Real)
    check_param(summarize, "irt_nstd", Real)
    check_param(summarize, "default_irt_width", Real)
    check_param(summarize, "peak_width_quantile", Real)
    check_param(summarize, "max_irt_bin_size", Real)

    # Validate quant_search_params
    quant_search = params["quant_search_params"]
    check_param(quant_search, "WH_smoothing_strength", Real)
    check_param(quant_search, "min_frag_count", Integer)
    check_param(quant_search, "min_y_count", Integer)
    check_param(quant_search, "min_log2_matched_ratio", Real)
    check_param(quant_search, "min_spectral_contrast", Real)
    check_param(quant_search, "min_topn_of_m", Vector)
    check_param(quant_search, "n_frag_isotopes", Integer)
    check_param(quant_search, "max_best_rank", Integer)
    check_param(quant_search, "n_pad", Integer)
    check_param(quant_search, "max_apex_offset", Integer)

    # Validate irt_mapping_params
    irt_mapping = params["irt_mapping_params"]
    check_param(irt_mapping, "n_bins", Integer)
    check_param(irt_mapping, "bandwidth", Real)
    check_param(irt_mapping, "n_sigma_tol", Real)
    check_param(irt_mapping, "min_prob", Real)

    # Validate deconvolution_params
    deconvolution = params["deconvolution_params"]
    check_param(deconvolution, "lambda", Real)
    check_param(deconvolution, "huber_delta", Real)
    check_param(deconvolution, "huber_delta0", Real)
    check_param(deconvolution, "huber_delta_exp", Real)
    check_param(deconvolution, "huber_delta_iters", Integer)
    check_param(deconvolution, "max_iter_newton", Integer)
    check_param(deconvolution, "max_iter_bisection", Integer)
    check_param(deconvolution, "max_iter_outer", Integer)
    check_param(deconvolution, "accuracy_newton", Real)
    check_param(deconvolution, "accuracy_bisection", Real)
    check_param(deconvolution, "max_diff", Real)

    # Validate qc_plot_params
    qc_plot = params["qc_plot_params"]
    check_param(qc_plot, "n_files_per_plot", Integer)

    # Validate normalization_params
    normalization = params["normalization_params"]
    check_param(normalization, "n_rt_bins", Integer)
    check_param(normalization, "spline_n_knots", Integer)

    # Validate xgboost_params
    xgboost = params["xgboost_params"]
    check_param(xgboost, "max_n_samples", Integer)
    check_param(xgboost, "min_best_trace_prob", Real)
    check_param(xgboost, "precursor_prob_spline_points_per_bin", Integer)
    check_param(xgboost, "precursor_q_value_interpolation_points_per_bin", Integer)
    check_param(xgboost, "pg_prob_spline_points_per_bin", Integer)
    check_param(xgboost, "pg_q_value_interpolation_points_per_bin", Integer)

    # Validate quad_transmission
    quad_transmission = params["quad_transmission"]
    check_param(quad_transmission, "fit_from_data", Bool)
    check_param(quad_transmission, "overhang", Real)
    check_param(quad_transmission, "smoothness", Real)

    # Validate benchmark_params
    benchmark = params["benchmark_params"]
    check_param(benchmark, "results_folder", String)

    # Validate output_params
    output = params["output_params"]
    check_param(output, "write_csv", Bool)
    check_param(output, "delete_temp", Bool)

    # Check remaining top-level parameters
    check_param(params, "library_folder", String)
    check_param(params, "ms_data_dir", String)

    # If all checks pass, return the validated parameters
    return params
end
