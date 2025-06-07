using Pioneer
using JSON

# Load test parameters and modify to only run Huber tuning
params_file = "./data/ecoli_test/ecoli_test_params.json"
params = JSON.parsefile(params_file)

# Modify to only run up to Huber tuning
params["search_method_names"] = ["ParameterTuning", "NCETuning", "QuadTuning", "FirstPass", "HuberTuning"]

# Write modified params
test_params_file = "./data/ecoli_test/test_huber_params.json"
open(test_params_file, "w") do f
    JSON.print(f, params, 4)
end

println("Running search up to Huber tuning...")
try
    SearchDIA(test_params_file)
catch e
    println("Error during search: $e")
    rethrow(e)
end