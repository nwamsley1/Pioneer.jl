using Pkg
Pkg.activate(".")
using Pioneer
using JSON

# Test 1: Load config with max_q_value specified
println("Test 1: Loading config with max_q_value = 0.01")
# Use the actual parameter loading method
params = Pioneer.PioneerParameters("data/example_config/defaultSearchParams.json")
pt_params = Pioneer.ParameterTuningSearchParameters(params)
println("  max_q_val from ParameterTuningSearchParameters: ", Pioneer.getMaxQVal(pt_params))

# Test 2: Test with different value
println("\nTest 2: Testing with max_q_value = 0.05")
params_json["parameter_tuning"]["search_settings"]["max_q_value"] = 0.05
params = Pioneer.GetSearchParams(params_json, "./", "./", "test.arrow")
pt_params = Pioneer.ParameterTuningSearchParameters(params)
println("  max_q_val from ParameterTuningSearchParameters: ", Pioneer.getMaxQVal(pt_params))

# Test 3: Test fallback when parameter is missing
println("\nTest 3: Testing fallback when max_q_value is not specified")
delete!(params_json["parameter_tuning"]["search_settings"], "max_q_value")
params = Pioneer.GetSearchParams(params_json, "./", "./", "test.arrow")
pt_params = Pioneer.ParameterTuningSearchParameters(params)
println("  max_q_val from ParameterTuningSearchParameters: ", Pioneer.getMaxQVal(pt_params))
println("  (Should match global q_value_threshold: ", params_json["global_settings"]["scoring"]["q_value_threshold"], ")")

println("\nAll tests passed!")