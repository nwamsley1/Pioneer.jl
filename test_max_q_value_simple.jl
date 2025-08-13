using Pkg
Pkg.activate(".")
using Pioneer

println("Testing max_q_value parameter extraction...")

# Test 1: Load config with max_q_value = 0.01 (as set in defaultSearchParams.json)
println("\nTest 1: Loading default config with max_q_value = 0.01")
params = Pioneer.PioneerParameters("data/example_config/defaultSearchParams.json")
pt_params = Pioneer.ParameterTuningSearchParameters(params)
max_q_val = Pioneer.getMaxQVal(pt_params)
println("  max_q_val extracted: ", max_q_val)
@assert max_q_val == 0.01f0 "Expected 0.01, got $max_q_val"
println("  ✓ Test 1 passed!")

# Test 2: Verify global q_value_threshold is different (it's 0.001 in the config)
println("\nTest 2: Checking that parameter tuning q-value is independent of global")
global_q_val = params.global_settings.scoring.q_value_threshold
println("  Global q_value_threshold: ", global_q_val)
println("  Parameter tuning max_q_val: ", max_q_val)
@assert global_q_val == 0.001 "Expected global threshold to be 0.001"
@assert max_q_val != global_q_val "Parameter tuning should have independent q-value"
println("  ✓ Test 2 passed!")

println("\n✅ All tests passed! The max_q_value parameter is working correctly.")