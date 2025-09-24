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

using JSON
#see src/Routines/BuildspecKib/utils
#struct InvalidParametersError <: Exception
#    message::String
#    params::Dict{String, Any}
#end

# Parameter default functions are loaded from paramDefaults.jl

function checkParams(json_path::String)
    # Read user params
    user_params = JSON.parsefile(json_path)
    
    # Apply defaults before validation
    defaults = get_default_parameters()
    params = merge_with_defaults(user_params, defaults)
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
    # tol_ppm moved to iteration_settings as init_mass_tol_ppm
    # min_score can be either Integer or Vector of Integers
    if !haskey(frag_settings, "min_score")
        throw(InvalidParametersError("Missing parameter: min_score", frag_settings))
    elseif !(frag_settings["min_score"] isa Integer || frag_settings["min_score"] isa Vector)
        throw(InvalidParametersError("Invalid type for parameter min_score: expected Integer or Vector, got $(typeof(frag_settings["min_score"]))", frag_settings))
    elseif frag_settings["min_score"] isa Vector
        # If it's a vector, make sure all elements are integers
        for (i, val) in enumerate(frag_settings["min_score"])
            if !(val isa Integer)
                throw(InvalidParametersError("Invalid type for min_score[$i]: expected Integer, got $(typeof(val))", frag_settings))
            end
        end
    end
    check_param(frag_settings, "min_spectral_contrast", Real)
    check_param(frag_settings, "relative_improvement_threshold", Real)
    check_param(frag_settings, "min_log2_ratio", Real)
    check_param(frag_settings, "min_top_n", Vector)
    check_param(frag_settings, "n_isotopes", Integer)
    
    # Check optional intensity filter quantile
    if haskey(frag_settings, "intensity_filter_quantile")
        check_param(frag_settings, "intensity_filter_quantile", Real)
        # Validate it's between 0 and 1
        if frag_settings["intensity_filter_quantile"] < 0 || frag_settings["intensity_filter_quantile"] >= 1
            error("parameter_tuning.fragment_settings.intensity_filter_quantile must be between 0 and 1 (got $(frag_settings["intensity_filter_quantile"]))")
        end
    end

    check_param(tuning_params, "search_settings", Dict)
    search_settings = tuning_params["search_settings"]
    check_param(search_settings, "min_samples", Integer)
    check_param(search_settings, "max_presearch_iters", Integer)
    check_param(search_settings, "frag_err_quantile", Real)
    
    # Check optional parameters if present
    if haskey(search_settings, "topn_peaks")
        check_param(search_settings, "topn_peaks", Integer)
    end
    # max_tolerance_ppm removed - calculated dynamically from iteration_settings
    
    # Check iteration_settings (all fields are now required)
    if haskey(tuning_params, "iteration_settings")
        iter_settings = tuning_params["iteration_settings"]
        
        # All iteration_settings fields are required
        check_param(iter_settings, "init_mass_tol_ppm", Real)
        check_param(iter_settings, "mass_tolerance_scale_factor", Real)
        check_param(iter_settings, "iterations_per_phase", Integer)
        check_param(iter_settings, "scan_scale_factor", Real)
        
        # Validate ranges
        if iter_settings["mass_tolerance_scale_factor"] <= 1.0
            error("iteration_settings.mass_tolerance_scale_factor must be greater than 1.0")
        end
        if iter_settings["scan_scale_factor"] < 1.0
            error("iteration_settings.scan_scale_factor must be greater than or equal to 1.0")
        end
        if iter_settings["init_mass_tol_ppm"] <= 0.0
            error("iteration_settings.init_mass_tol_ppm must be positive")
        end
        if iter_settings["iterations_per_phase"] <= 0
            error("iteration_settings.iterations_per_phase must be positive")
        end
    else
        error("parameter_tuning.iteration_settings is required")
    end
    if haskey(search_settings, "initial_scan_count")
        check_param(search_settings, "initial_scan_count", Integer)
    end
    # Check for both old and new parameter names for backward compatibility
    if haskey(search_settings, "max_parameter_tuning_scans")
        check_param(search_settings, "max_parameter_tuning_scans", Integer)
    elseif haskey(search_settings, "expanded_scan_count")
        check_param(search_settings, "expanded_scan_count", Integer)
    end
    if haskey(search_settings, "max_frags_for_mass_err_estimation")
        check_param(search_settings, "max_frags_for_mass_err_estimation", Integer)
    end
    if haskey(search_settings, "max_q_value")
        check_param(search_settings, "max_q_value", Real)
        # Validate it's between 0 and 1
        if search_settings["max_q_value"] <= 0 || search_settings["max_q_value"] > 1
            error("parameter_tuning.search_settings.max_q_value must be between 0 and 1 (got $(search_settings["max_q_value"]))")
        end
    end

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
    check_param(score_settings, "max_PEP", Real)

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
    check_param(ml_params, "max_q_value_mbr_rescore", Real)
    check_param(ml_params, "min_PEP_neg_threshold_rescore", Real)
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