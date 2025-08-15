using Pioneer
using Test

# Simple test for intelligent scan selection
@testset "Intelligent Scan Selection" begin
    # Load test data
    test_data_path = joinpath(@__DIR__, "..", "data", "ecoli_test", "raw")
    ms_data_path = joinpath(test_data_path, "ecoli_filtered_01.arrow")
    
    if isfile(ms_data_path)
        # Load mass spec data
        ms_data = Pioneer.BasicMassSpecData(ms_data_path)
        
        # Test 1: Create FilteredMassSpecData with intelligent selection
        filtered = Pioneer.FilteredMassSpecData(
            ms_data,
            max_scans = 15,  # Request just 15 scans (one per RT bin)
            target_ms_order = UInt8(2),
            n_rt_bins = 15
        )
        
        @test length(filtered) == 15
        @test length(filtered.scan_priority_order) > 0
        @test filtered.n_scans_sampled == 15
        @test length(filtered.rt_bin_assignments) == length(ms_data)
        
        # Test 2: Check RT coverage - should have scans distributed across RT range
        rts = [Pioneer.getRetentionTime(filtered, i) for i in 1:length(filtered)]
        rt_range = maximum(rts) - minimum(rts)
        original_rt_range = maximum(Pioneer.getRetentionTimes(ms_data)) - minimum(Pioneer.getRetentionTimes(ms_data))
        
        # RT coverage should be at least 50% of original range
        @test rt_range > 0.5 * original_rt_range
        
        # Test 3: Test append functionality
        initial_count = length(filtered)
        n_added = append!(filtered, max_additional_scans = 10)
        
        @test n_added == 10
        @test length(filtered) == initial_count + 10
        @test filtered.n_scans_sampled == 25
        
        # Test 4: Verify deterministic ordering
        filtered2 = Pioneer.FilteredMassSpecData(
            ms_data,
            max_scans = 15,
            target_ms_order = UInt8(2),
            n_rt_bins = 15
        )
        
        # Should get same scan indices in same order
        @test filtered2.original_scan_indices[1:15] == filtered.original_scan_indices[1:15]
        
        println("✓ All intelligent scan selection tests passed!")
    else
        @warn "Test data file not found at $ms_data_path, skipping tests"
    end
end