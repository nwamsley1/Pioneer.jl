module ParamDefaults

export get_default_parameters, merge_with_defaults

"""
    get_default_parameters()

Returns the complete default parameter structure for SearchDIA.
All parameters have sensible defaults that work for most experiments.
Users only need to override the values they want to change.
"""
function get_default_parameters()
    return Dict(
        "logging" => Dict(
            "debug_console_level" => 0
        ),
        "global" => Dict(
            "isotope_settings" => Dict(
                "err_bounds_first_pass" => [1, 0],
                "err_bounds_quant_search" => [3, 0],
                "combine_traces" => true,
                "partial_capture" => true,
                "min_fraction_transmitted" => 0.25
            ),
            "scoring" => Dict(
                "q_value_threshold" => 0.01
            ),
            "normalization" => Dict(
                "n_rt_bins" => 100,
                "spline_n_knots" => 7
            ),
            "huber_override" => Dict(
                "override_huber_delta_fit" => false,
                "huber_delta" => 1055
            ),
            "match_between_runs" => true,
            "ms1_quant" => false,
            "ms1_scoring" => true
        ),
        "parameter_tuning" => Dict(
            "fragment_settings" => Dict(
                "min_count" => 7,
                "max_rank" => 25,
                "min_score" => [22, 17, 15],
                "min_spectral_contrast" => 0.5,
                "relative_improvement_threshold" => 1.5,
                "min_log2_ratio" => 1.5,
                "min_top_n" => [3, 3],
                "n_isotopes" => 1,
                "intensity_filter_quantile" => 0.50
            ),
            "search_settings" => Dict(
                "sample_rate" => 0.02,
                "initial_scan_count" => 10000,
                "max_parameter_tuning_scans" => 40000,
                "max_q_value" => 0.01,
                "min_samples" => 1000,
                "topn_peaks" => 200,
                "max_frags_for_mass_err_estimation" => 5,
                "min_quad_tuning_psms" => 5000,
                "min_quad_tuning_fragments" => 3,
                "max_presearch_iters" => 10,
                "frag_err_quantile" => 0.0025
            ),
            "iteration_settings" => Dict(
                "init_mass_tol_ppm" => 10.0,
                "mass_tolerance_scale_factor" => 1.5,
                "iterations_per_phase" => 3,
                "scan_scale_factor" => 4
            )
        ),
        "first_search" => Dict(
            "fragment_settings" => Dict(
                "min_count" => 4,
                "max_rank" => 25,
                "min_score" => 12,
                "min_spectral_contrast" => 0.5,
                "relative_improvement_threshold" => 1.25,
                "min_log2_ratio" => 0.0,
                "min_top_n" => [2, 3],
                "n_isotopes" => 1
            ),
            "scoring_settings" => Dict(
                "n_train_rounds" => 2,
                "max_iterations" => 20,
                "max_q_value_probit_rescore" => 0.05,
                "max_PEP" => 0.9
            ),
            "irt_mapping" => Dict(
                "max_prob_to_impute_irt" => 0.75,
                "fwhm_nstd" => 4,
                "irt_nstd" => 4
            )
        ),
        "quant_search" => Dict(
            "fragment_settings" => Dict(
                "min_count" => 3,
                "min_y_count" => 2,
                "max_rank" => 255,
                "min_spectral_contrast" => 0.0,
                "min_log2_ratio" => -1.7,
                "min_top_n" => [2, 3],
                "n_isotopes" => 2
            ),
            "chromatogram" => Dict(
                "smoothing_strength" => 1e-6,
                "padding" => 0,
                "max_apex_offset" => 2
            )
        ),
        "acquisition" => Dict(
            "nce" => 25,
            "quad_transmission" => Dict(
                "fit_from_data" => false,
                "overhang" => 0.25,
                "smoothness" => 5.0
            )
        ),
        "rt_alignment" => Dict(
            "n_bins" => 200,
            "bandwidth" => 0.25,
            "sigma_tolerance" => 4,
            "min_probability" => 0.95
        ),
        "optimization" => Dict(
            "deconvolution" => Dict(
                "lambda" => 0.0,
                "reg_type" => "none",
                "huber_delta" => 300,
                "huber_exp" => 1.5,
                "huber_iters" => 15,
                "newton_iters" => 50,
                "bisection_iters" => 100,
                "outer_iters" => 1000,
                "newton_accuracy" => 10,
                "max_diff" => 0.01
            ),
            "machine_learning" => Dict(
                "max_psms_in_memory" => 50000000,
                "min_trace_prob" => 0.75,
                "max_q_value_xgboost_mbr_rescore" => 0.20,
                "min_PEP_neg_threshold_xgboost_rescore" => 0.90,
                "spline_points" => 500,
                "interpolation_points" => 10
            )
        ),
        "proteinInference" => Dict(
            "min_peptides" => 1
        ),
        "maxLFQ" => Dict(
            "run_to_run_normalization" => false
        ),
        "output" => Dict(
            "write_csv" => true,
            "write_decoys" => false,
            "delete_temp" => true,
            "plots_per_page" => 12
        ),
        # paths section will be provided by user and not defaulted
        "paths" => Dict()
    )
end

"""
    merge_with_defaults(user_params::Dict, defaults::Dict)

Recursively merges user parameters over default parameters.
User values override defaults at any nesting level.
Missing sections or parameters are filled from defaults.
"""
function merge_with_defaults(user_params::Dict, defaults::Dict)
    result = deepcopy(defaults)
    
    function recursive_merge!(target::Dict, source::Dict)
        for (key, value) in source
            if haskey(target, key) && isa(target[key], Dict) && isa(value, Dict)
                # Both are dicts, merge recursively
                recursive_merge!(target[key], value)
            else
                # Override with user value
                target[key] = value
            end
        end
    end
    
    recursive_merge!(result, user_params)
    return result
end

end # module