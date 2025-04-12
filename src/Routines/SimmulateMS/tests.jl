"""
Test the peak merging functionality
"""
function test_peak_merging()
    # Test case 1: Peaks within tolerance
    mz1 = [100.0, 200.0, 400.0, 400.05]
    intensity1 = [1000.0, 500.0, 800.0, 600.0]
    
    # 400.05 should merge with 400.0 (125 ppm difference)
    merged_mz, merged_intensity = merge_peaks(mz1, intensity1, 150.0)
    @test length(merged_mz) == 3
    @test isapprox(merged_mz[3], 400.03, atol=0.01)  # Weighted average
    @test isapprox(merged_intensity[3], 1400.0)  # Sum of intensities
    
    # Test case 2: No merging needed
    mz2 = [100.0, 200.0, 300.0, 400.0]
    intensity2 = [1000.0, 500.0, 800.0, 600.0]
    merged_mz, merged_intensity = merge_peaks(mz2, intensity2, 5.0)
    @test length(merged_mz) == 4
    
    # Test case 3: Multiple merges
    mz3 = [100.0, 100.001, 100.002, 200.0, 200.001]
    intensity3 = [1000.0, 500.0, 200.0, 800.0, 400.0]
    merged_mz, merged_intensity = merge_peaks(mz3, intensity3, 20.0)
    @test length(merged_mz) == 2
    @test isapprox(merged_mz[1], 100.001, atol=0.0005)  # Weighted average
    @test isapprox(merged_intensity[1], 1700.0)  # Sum of intensities
    
    println("Peak merging tests passed!")
end
