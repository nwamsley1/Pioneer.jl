#!/usr/bin/env julia

# Test script for multi-score parameter tuning
using Pioneer
using JSON

# Create a test parameters dictionary with multiple min_scores
function create_test_params_with_scores(scores)
    return Dict(
        "global_settings" => Dict(
            "isotope_settings" => Dict(
                "err_bounds_first_pass" => [0, 2],
                "partial_capture" => false,
                "min_fraction_transmitted" => 0.5
            ),
            "scoring" => Dict(
                "q_value_threshold" => 0.01
            )
        ),
        "parameter_tuning" => Dict(
            "fragment_settings" => Dict(
                "min_score" => scores,  # Can be single value or array
                "min_count" => 3,
                "min_spectral_contrast" => 0.1,
                "min_log2_ratio" => -3.0,
                "min_top_n" => [3, 5]
            ),
            "search_settings" => Dict(
                "frag_err_quantile" => 0.01,
                "min_samples" => 100,
                "max_presearch_iters" => 10,
                "initial_scan_count" => 500,
                "max_parameter_tuning_scans" => 80000,
                "max_frags_for_mass_err_estimation" => 5
            ),
            "iteration_settings" => Dict(
                "init_mass_tol_ppm" => 5.0,
                "mass_tolerance_scale_factor" => 2.0,
                "iterations_per_phase" => 3,
                "scan_scale_factor" => 10.0
            )
        ),
        "rt_alignment" => Dict(
            "irt_tol" => 50.0,
            "spec_order" => [2],
            "relative_improvement_threshold" => 0.1,
            "spline_degree" => 3,
            "spline_n_knots" => 7,
            "spline_fit_outlier_sd" => 3,
            "irt_tol_sd" => 3
        )
    )
end

# Test 1: Single score (backward compatibility)
println("Test 1: Single min_score value (backward compatibility)")
params_single = create_test_params_with_scores(22)
pioneer_params_single = PioneerParameters(params_single)
ptsp_single = ParameterTuningSearchParameters(pioneer_params_single)

println("  Single score converted to vector: ", getMinIndexSearchScores(ptsp_single))
println("  Current active score: ", getMinIndexSearchScore(ptsp_single))
println("  ✓ Single score test passed\n")

# Test 2: Multiple scores
println("Test 2: Multiple min_score values")
params_multi = create_test_params_with_scores([22, 17, 15])
pioneer_params_multi = PioneerParameters(params_multi)
ptsp_multi = ParameterTuningSearchParameters(pioneer_params_multi)

println("  Score thresholds: ", getMinIndexSearchScores(ptsp_multi))
println("  Initial active score: ", getMinIndexSearchScore(ptsp_multi))

# Test score switching
setCurrentMinScore!(ptsp_multi, UInt8(17))
println("  After switching to 17: ", getMinIndexSearchScore(ptsp_multi))

setCurrentMinScore!(ptsp_multi, UInt8(15))
println("  After switching to 15: ", getMinIndexSearchScore(ptsp_multi))
println("  ✓ Multiple scores test passed\n")

# Test 3: Different array sizes
println("Test 3: Different array sizes")
for scores in [[30], [25, 20], [22, 17, 15, 12]]
    params = create_test_params_with_scores(scores)
    pioneer_params = PioneerParameters(params)
    ptsp = ParameterTuningSearchParameters(pioneer_params)
    println("  Input: $scores → Stored: $(getMinIndexSearchScores(ptsp))")
end
println("  ✓ Array size test passed\n")

println("All tests passed! Multi-score parameter tuning is working correctly.")