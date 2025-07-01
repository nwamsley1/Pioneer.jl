# Test script to verify the protein scoring fix

push!(LOAD_PATH, dirname(@__DIR__))
using EntrapmentAnalysis
using DataFrames

# Test 1: Basic functionality with the fixed code
println("Test 1: Basic protein scoring functionality")
test_df = DataFrame(
    protein = ["PROT_A", "PROT_A", "PROT_B", "PROT_B"],
    entrap_id = UInt8[0, 1, 0, 1],
    file_name = ["Rep1", "Rep1", "Rep1", "Rep1"],
    pg_score = Float32[100.0, 90.0, 80.0, 70.0]
)

println("Input data:")
println(test_df)

try
    add_original_target_protein_scores!(test_df; score_col=:pg_score)
    println("\nSuccess! Original target scores added:")
    println(select(test_df, :protein, :entrap_id, :file_name, :pg_score, :pg_score_original_target))
catch e
    println("\nError: $e")
end

# Test 2: Duplicate target protein (should throw error)
println("\n\nTest 2: Duplicate target protein detection")
test_df2 = DataFrame(
    protein = ["PROT_A", "PROT_A", "PROT_A"],
    entrap_id = UInt8[0, 0, 1],  # Two targets for same protein
    file_name = ["Rep1", "Rep1", "Rep1"],
    pg_score = Float32[100.0, 95.0, 90.0]
)

println("Input data with duplicate target:")
println(test_df2)

try
    add_original_target_protein_scores!(test_df2; score_col=:pg_score)
    println("\nUnexpected success - should have thrown error!")
catch e
    println("\nExpected error caught: $e")
end

# Test 3: Real-world scenario with multiple files
println("\n\nTest 3: Multiple files scenario")
test_df3 = DataFrame(
    protein = ["PROT_A", "PROT_A", "PROT_A", "PROT_A"],
    entrap_id = UInt8[0, 1, 0, 1],
    file_name = ["Rep1", "Rep1", "Rep2", "Rep2"],
    pg_score = Float32[100.0, 90.0, 105.0, 95.0]
)

println("Input data with multiple files:")
println(test_df3)

try
    add_original_target_protein_scores!(test_df3; score_col=:pg_score)
    println("\nSuccess! Original target scores added:")
    println(select(test_df3, :protein, :entrap_id, :file_name, :pg_score, :pg_score_original_target))
catch e
    println("\nError: $e")
end

println("\n\nAll tests completed!")