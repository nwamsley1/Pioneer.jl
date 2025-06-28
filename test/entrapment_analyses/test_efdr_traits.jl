# Test for trait-based EFDR system
using Test, DataFrames

# Create test data
function create_test_data()
    # Create library_precursors with entrapment groups
    library_precursors = DataFrame(
        entrapment_group_id = [0, 1, 0, 1, 0, 2],
        entrap_pair_id = [UInt32(1), UInt32(1), UInt32(2), UInt32(2), UInt32(3), UInt32(3)]
    )
    
    # Create prec_results with scores and q-values
    prec_results = DataFrame(
        precursor_idx = [1, 2, 3, 4, 5, 6],
        global_prob = [0.9, 0.8, 0.95, 0.7, 0.85, 0.6],
        global_qval = [0.01, 0.02, 0.005, 0.03, 0.015, 0.04],
        prec_prob = [0.88, 0.82, 0.92, 0.75, 0.87, 0.65],
        qval = [0.012, 0.018, 0.008, 0.025, 0.013, 0.035]
    )
    
    return prec_results, library_precursors
end

@testset "EFDR Trait System Tests" begin
    prec_results, library_precursors = create_test_data()
    
    # Test adding EFDR columns
    add_efdr_columns!(prec_results, library_precursors)
    
    # Check that all expected columns were created
    expected_columns = [
        :global_prob_combined_efdr,
        :global_prob_paired_efdr,
        :prec_prob_combined_efdr,
        :prec_prob_paired_efdr
    ]
    
    for col in expected_columns
        @test hasproperty(prec_results, col)
        @test all(!ismissing(prec_results[!, col]))
        @test all(0 .<= prec_results[!, col] .<= 1)  # EFDR should be between 0 and 1
    end
    
    # Test with custom parameters
    prec_results2, library_precursors2 = create_test_data()
    add_efdr_columns!(prec_results2, library_precursors2;
                     methods = [CombinedEFDR()],  # Only combined
                     score_qval_pairs = [(:global_prob, :global_qval)],  # Only one pair
                     r = 2.0)
    
    @test hasproperty(prec_results2, :global_prob_combined_efdr)
    @test !hasproperty(prec_results2, :global_prob_paired_efdr)
    @test !hasproperty(prec_results2, :prec_prob_combined_efdr)
    
    # Test direct trait dispatch
    prec_results3, library_precursors3 = create_test_data()
    add_entrap_pair_ids!(prec_results3, library_precursors3)
    add_original_target_scores!(prec_results3, library_precursors3; score_col=:global_prob)
    
    entrapment_labels = [library_precursors3.entrapment_group_id[pid] for pid in prec_results3.precursor_idx]
    
    combined_efdr = calculate_efdr(CombinedEFDR(), 
                                  prec_results3.global_prob,
                                  prec_results3.global_prob_original_target,
                                  entrapment_labels,
                                  prec_results3.global_qval)
    
    paired_efdr = calculate_efdr(PairedEFDR(),
                                prec_results3.global_prob,
                                prec_results3.global_prob_original_target,
                                entrapment_labels,
                                prec_results3.global_qval)
    
    @test length(combined_efdr) == nrow(prec_results3)
    @test length(paired_efdr) == nrow(prec_results3)
    @test all(0 .<= combined_efdr .<= 1)
    @test all(0 .<= paired_efdr .<= 1)
    
    println("All EFDR trait system tests passed!")
end