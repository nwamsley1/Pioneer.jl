#!/usr/bin/env julia

# Test script to compare probit regression vs XGBoost for PSM scoring

using Pioneer
using DataFrames
using Arrow
using Statistics

println("Testing Probit Regression vs XGBoost on E. coli dataset")
println("=" ^ 60)

# Function to run SearchDIA and get protein count
function run_search_and_count(params_file::String, description::String)
    println("\n$description")
    println("-" ^ 40)
    
    # Run SearchDIA
    try
        SearchDIA(params_file)
        
        # Check output files
        output_dir = dirname(params_file) * "/results"
        
        # Count proteins in protein groups file
        pg_file = "$output_dir/protein_groups_long.arrow"
        if isfile(pg_file)
            pg_df = Arrow.Table(pg_file) |> DataFrame
            n_proteins = length(unique(pg_df.protein_name))
            n_passing = sum(pg_df.global_qval .<= 0.01)
            println("Total protein groups: $n_proteins")
            println("Passing protein groups (q≤0.01): $n_passing")
            
            # Count PSMs
            psm_file = "$output_dir/precursors_long.arrow"
            if isfile(psm_file)
                psm_df = Arrow.Table(psm_file) |> DataFrame
                n_psms = nrow(psm_df)
                n_passing_psms = sum(psm_df.q_value .<= 0.01)
                println("Total PSMs: $n_psms")
                println("Passing PSMs (q≤0.01): $n_passing_psms")
            end
        else
            println("Warning: Protein groups file not found!")
        end
        
    catch e
        println("Error running SearchDIA: $e")
    end
end

# Path to ecoli test parameters
params_file = "./data/ecoli_test/ecoli_test_params.json"

if !isfile(params_file)
    println("Error: E. coli test parameters file not found at $params_file")
    println("Please ensure you're running from the Pioneer.jl root directory")
    exit(1)
end

# First run with XGBoost (default)
println("\n" * "=" ^ 60)
println("Running with XGBoost (default)")
run_search_and_count(params_file, "XGBoost Results")

# Now we need to modify the code to use probit regression
# This requires uncommenting the probit code and commenting out XGBoost
println("\n" * "=" ^ 60)
println("\nTo test with Probit Regression:")
println("1. Edit src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl")
println("2. Comment out lines 479-499 (XGBoost call)")
println("3. Uncomment lines 468-477 (Probit call)")
println("4. Rerun this script")

println("\n" * "=" ^ 60)
println("\nNote: The current implementation uses XGBoost by default.")
println("To switch to probit regression, manual code changes are needed.")