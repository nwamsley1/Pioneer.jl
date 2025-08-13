using Pkg
Pkg.activate(".")
using Pioneer
using JSON

println("Testing max_q_value parameter extraction inline...")

# Parse the JSON directly
json_data = JSON.parsefile("data/example_config/defaultSearchParams.json")

# Check that the parameter is in the JSON
max_q_from_json = json_data["parameter_tuning"]["search_settings"]["max_q_value"]
println("max_q_value in JSON: ", max_q_from_json)
@assert max_q_from_json == 0.01

# Now test that parameter extraction works by loading parameters
# We'll use the internal function
params = Pioneer.parse_pioneer_parameters("data/example_config/defaultSearchParams.json")

# Create ParameterTuningSearchParameters
pt_params = Pioneer.ParameterTuningSearchParameters(params)

# Get the max_q_val
max_q_val = Pioneer.getMaxQVal(pt_params)
println("max_q_val extracted: ", max_q_val)

# Verify it matches
@assert max_q_val == Float32(0.01) "Expected 0.01, got $max_q_val"

# Also verify it's different from global
global_q = params.global_settings.scoring.q_value_threshold
println("Global q_value_threshold: ", global_q)
@assert global_q == 0.001 "Global should be 0.001"
@assert max_q_val != global_q "Should be independent values"

println("\n✅ Success! max_q_value parameter is correctly extracted and used.")