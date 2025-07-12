# Test script for global score analysis implementation

push!(LOAD_PATH, dirname(dirname(@__DIR__)))

using EntrapmentAnalysis
using DataFrames
using Test

println("Testing global score analysis implementation...")

# Create test data
function create_test_data()
    # Create test results with multiple files
    prec_results = DataFrame(
        precursor_idx = UInt32[1, 1, 2, 2, 3, 3, 4, 4],
        ms_file_idx = [1, 2, 1, 2, 1, 2, 1, 2],
        global_prob = Float32[0.9, 0.8, 0.7, 0.95, 0.6, 0.5, 0.85, 0.4],
        prec_prob = Float32[0.9, 0.8, 0.7, 0.95, 0.6, 0.5, 0.85, 0.4],
        global_qval = Float32[0.01, 0.02, 0.03, 0.005, 0.04, 0.05, 0.015, 0.06],
        qval = Float32[0.01, 0.02, 0.03, 0.005, 0.04, 0.05, 0.015, 0.06],
        target = true
    )
    
    # Create test library
    library_precursors = DataFrame(
        precursor_idx = UInt32[1, 2, 3, 4],
        sequence = ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "PEPTIDE4"],
        entrapment_group_id = UInt8[0, 1, 0, 1],
        structural_mods = ["", "", "", ""],
        mod_key = ["", "", "", ""]
    )
    
    return prec_results, library_precursors
end

# Test create_global_results_df function
println("\nTest 1: create_global_results_df function")
prec_results, library_precursors = create_test_data()

global_df = create_global_results_df(prec_results; score_col=:global_prob)

@test nrow(global_df) == 4  # Should have 4 unique precursors
@test all(global_df.ms_file_idx .== 0)  # All should be 0
@test global_df[global_df.precursor_idx .== 2, :global_prob][1] ≈ 0.95  # Should keep max score

println("✓ create_global_results_df works correctly")

# Test that original dataframe is unchanged
@test nrow(prec_results) == 8
@test unique(prec_results.ms_file_idx) == [1, 2]
println("✓ Original dataframe unchanged")

# Test entrapment pair assignment
println("\nTest 2: Entrapment pair processing")
assign_entrapment_pairs!(library_precursors)
add_entrap_pair_ids!(prec_results, library_precursors)
add_entrap_pair_ids!(global_df, library_precursors)

# Test original target scores
println("\nTest 3: Original target scores")
add_original_target_scores!(prec_results, library_precursors, [:prec_prob])
add_original_target_scores!(global_df, library_precursors, [:global_prob])

@test hasproperty(prec_results, :prec_prob_original_target)
@test hasproperty(global_df, :global_prob_original_target)
println("✓ Original target scores added correctly")

# Test EFDR calculation
println("\nTest 4: EFDR calculation")
add_efdr_columns!(prec_results, library_precursors;
                 method_types = [CombinedEFDR],
                 score_qval_pairs = [(:prec_prob, :qval)],
                 r = 1.0)

add_efdr_columns!(global_df, library_precursors;
                 method_types = [CombinedEFDR],
                 score_qval_pairs = [(:global_prob, :global_qval)],
                 r = 1.0)

@test hasproperty(prec_results, :prec_prob_combined_efdr)
@test hasproperty(global_df, :global_prob_combined_efdr)
println("✓ EFDR columns added correctly")

println("\nAll tests passed! ✓")
println("\nSummary:")
println("- Original dataframe: $(nrow(prec_results)) rows")
println("- Global dataframe: $(nrow(global_df)) rows")
println("- Per-file analysis uses ms_file_idx: $(unique(prec_results.ms_file_idx))")
println("- Global analysis uses ms_file_idx: $(unique(global_df.ms_file_idx))")