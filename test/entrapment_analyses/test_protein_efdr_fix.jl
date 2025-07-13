"""
Test script to verify the protein EFDR fix works correctly
"""

push!(LOAD_PATH, dirname(dirname(@__DIR__)))
using EntrapmentAnalysis
using DataFrames
using Arrow

println("Testing protein EFDR fix...")

# Test 1: Basic functionality with mock data
println("\n=== Test 1: Mock Data ===")
mock_data = DataFrame(
    file_name = ["file1", "file1", "file2", "file2"],
    target = fill(true, 4),
    entrap_id = UInt8[0, 1, 0, 1],
    protein = ["PROT1", "PROT1", "PROT1", "PROT1"],
    pg_score = Float32[100.0, 90.0, 95.0, 85.0],
    global_pg_score = Float32[100.0, 90.0, 100.0, 90.0],
    qval = Float32[0.01, 0.02, 0.01, 0.02],
    global_qval = Float32[0.01, 0.02, 0.01, 0.02]
)

# Save and test
temp_path = tempname() * ".arrow"
Arrow.write(temp_path, mock_data)

try
    results = run_protein_efdr_analysis(
        temp_path;
        output_dir=tempname(),
        verbose=true,
        plot_formats=[]
    )
    println("\n✅ Mock data test passed!")
catch e
    println("\n❌ Mock data test failed: $e")
    rethrow(e)
finally
    rm(temp_path, force=true)
end

# Test 2: Real data if available
real_data_path = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow"

if isfile(real_data_path)
    println("\n=== Test 2: Real Data ===")
    try
        results = run_protein_efdr_analysis(
            real_data_path;
            output_dir="protein_efdr_test_output",
            verbose=true,
            plot_formats=[:png]
        )
        println("\n✅ Real data test passed!")
        println("Check output in: protein_efdr_test_output/")
    catch e
        println("\n❌ Real data test failed: $e")
        rethrow(e)
    end
else
    println("\n⚠️  Real data file not found, skipping real data test")
end

println("\n=== All tests completed ===")

# Test 3: Verify parameter ordering fix
println("\n=== Test 3: Parameter Ordering Fix ===")
test_df = DataFrame(
    entrap_id = UInt8[0, 1, 0, 1],
    pg_score = Float32[100.0, 90.0, 80.0, 70.0],
    pg_score_original_target = Float32[100.0, 100.0, 80.0, 80.0],
    qval = Float32[0.0001, 0.0002, 0.0001, 0.0002]  # Very small q-values that would cause the error
)

try
    add_protein_efdr_columns!(test_df; score_qval_pairs=[(:pg_score, :qval)])
    println("✅ Parameter ordering fix verified - no InexactError!")
catch e
    println("❌ Parameter ordering still has issues: $e")
    rethrow(e)
end