# Test script for refactored EFDR implementation
using Test
using DataFrames
using Random

# Include the necessary files
include("entrapment_helper_funcs.jl")
include("efdr_funcs.jl")

# Create test data
Random.seed!(42)
n = 100

# Create mock data
score = rand(n)
original_target_score = rand(n)
entrapment_label = [i % 3 == 0 ? 0 : (i % 3) for i in 1:n]  # Mix of 0, 1, 2
qval = sort(rand(n) * 0.1)  # Q-values between 0 and 0.1, sorted
r = 1.5

println("Testing refactored EFDR implementation...")

# Test 1: Create method instances directly
println("\n1. Testing direct method instantiation...")
combined_method = CombinedEFDR(score, original_target_score, entrapment_label, qval; r=r)
paired_method = PairedEFDR(score, original_target_score, entrapment_label, qval; r=r)

@test combined_method.r == r
@test length(combined_method.score) == n
@test combined_method.score == score

# Test 2: Calculate EFDR using the new methods
println("2. Testing calculate_efdr with method instances...")
combined_efdr = calculate_efdr(combined_method)
paired_efdr = calculate_efdr(paired_method)

@test length(combined_efdr) == n
@test length(paired_efdr) == n

# Debug the values
println("Combined EFDR range: ", minimum(combined_efdr), " to ", maximum(combined_efdr))
println("Paired EFDR range: ", minimum(paired_efdr), " to ", maximum(paired_efdr))

# Test EFDR properties
@test !any(isnan.(combined_efdr))
@test !any(isnan.(paired_efdr))
@test all(0 .<= combined_efdr .<= 1)
@test all(0 .<= paired_efdr .<= 1)

# Test 3: Test backward compatibility functions
println("3. Testing backward compatibility...")
combined_efdr_old = get_combined_efdr(score, original_target_score, entrapment_label, qval, r)
paired_efdr_old = get_paired_efdr(score, original_target_score, entrapment_label, qval, r)

@test combined_efdr == combined_efdr_old
@test paired_efdr == paired_efdr_old

# Test 4: Test with DataFrames
println("4. Testing with DataFrames...")
library_precursors = DataFrame(
    entrapment_group_id = entrapment_label,
    entrap_pair_id = [UInt32(i ÷ 3) for i in 1:n]
)

prec_results = DataFrame(
    precursor_idx = 1:n,
    global_prob = score,
    global_qval = qval,
    prec_prob = score .* 0.95,
    qval = qval .* 1.1
)

# Add EFDR columns
add_efdr_columns!(prec_results, library_precursors; 
                 method_types=[CombinedEFDR, PairedEFDR],
                 score_qval_pairs=[(:global_prob, :global_qval)],
                 r=r)

@test hasproperty(prec_results, :global_prob_combined_efdr)
@test hasproperty(prec_results, :global_prob_paired_efdr)
@test hasproperty(prec_results, :global_prob_original_target)

# Test 5: Test error handling
println("5. Testing error handling...")
@test_throws ErrorException CombinedEFDR(score[1:50], original_target_score, entrapment_label, qval; r=r)

println("\nAll tests passed! ✅")