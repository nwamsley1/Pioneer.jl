# Example usage of the refactored EFDR implementation
using DataFrames
include("entrapment_helper_funcs.jl")
include("efdr_funcs.jl")

# Create example data
n = 1000
score = rand(n)
original_target_score = rand(n)
entrapment_label = [i % 4 == 0 ? 0 : 1 for i in 1:n]  # 25% targets, 75% entrapments
qval = sort(rand(n) * 0.1)  # Q-values between 0 and 0.1

# Method 1: Direct instantiation and calculation
println("=== Method 1: Direct API ===")
combined_method = CombinedEFDR(score, original_target_score, entrapment_label, qval; r=1.0)
paired_method = PairedEFDR(score, original_target_score, entrapment_label, qval; r=1.0)

combined_efdr = calculate_efdr(combined_method)
paired_efdr = calculate_efdr(paired_method)

println("Combined EFDR at FDR 0.01: ", combined_efdr[findfirst(qval .>= 0.01)])
println("Paired EFDR at FDR 0.01: ", paired_efdr[findfirst(qval .>= 0.01)])

# Method 2: Using DataFrames with the wrapper function
println("\n=== Method 2: DataFrame API ===")
library_precursors = DataFrame(
    entrapment_group_id = entrapment_label,
    entrap_pair_id = [UInt32(i ÷ 2) for i in 1:n]
)

prec_results = DataFrame(
    precursor_idx = 1:n,
    global_prob = score,
    global_qval = qval
)

# Add EFDR columns automatically
add_efdr_columns!(prec_results, library_precursors;
                 method_types = [CombinedEFDR, PairedEFDR],
                 score_qval_pairs = [(:global_prob, :global_qval)],
                 r = 1.0)

# Show the new columns
println("\nNew columns added: ", names(prec_results)[end-3:end])
println("\nFirst 5 rows:")
show(first(prec_results[:, [
    :global_qval, 
    :global_prob_combined_efdr, 
    :global_prob_paired_efdr
]], 5), truncate=false)

# Method 3: Creating custom methods with different parameters
println("\n\n=== Method 3: Custom Parameters ===")
custom_method = CombinedEFDR(score, original_target_score, entrapment_label, qval; r=2.0)
custom_efdr = calculate_efdr(custom_method)
println("EFDR with r=2.0 at FDR 0.01: ", custom_efdr[findfirst(qval .>= 0.01)])