using Test
using Pioneer
using Random
using Arrow
using DataFrames

@testset "FilteredMassSpecData Tests" begin
    # Create test data path
    test_data_path = joinpath(@__DIR__, "..", "..", "data", "ecoli_test", "file_0_2000_2100.arrow")
    
    # Check if test data exists, if not skip tests
    if !isfile(test_data_path)
        @warn "Test data not found at $test_data_path, skipping FilteredMassSpecData tests"
        return
    end
    
    original = BasicMassSpecData(test_data_path)
    
    @testset "Basic Construction" begin
        # Test with sampling only
        filtered = FilteredMassSpecData(original, 0.1, nothing)
        @test length(filtered) < length(original)
        @test length(filtered) > 0
        
        # Test interface compliance
        @test length(getMzArray(filtered, 1)) > 0
        @test getRetentionTime(filtered, 1) > 0
        @test getScanNumber(filtered, 1) > 0
        @test getTIC(filtered, 1) > 0
    end
    
    @testset "TopN Filtering" begin
        # Test with topN filtering
        filtered = FilteredMassSpecData(original, 1.0, 100, target_ms_order=UInt8(2))
        
        for i in 1:min(10, length(filtered))
            mz_array = getMzArray(filtered, i)
            @test length(mz_array) <= 100
            # Verify m/z ordering
            @test issorted(mz_array)
        end
    end
    
    @testset "Index Mapping" begin
        filtered = FilteredMassSpecData(original, 0.5, nothing, seed=42)
        
        for i in 1:min(10, length(filtered))
            orig_idx = getOriginalScanIndex(filtered, i)
            # Verify metadata matches
            @test getScanNumber(filtered, i) == getScanNumber(original, orig_idx)
            @test getRetentionTime(filtered, i) ≈ getRetentionTime(original, orig_idx)
            @test getMsOrder(filtered, i) == getMsOrder(original, orig_idx)
        end
    end
    
    @testset "MS Order Filtering" begin
        # Test MS2 only filtering
        filtered_ms2 = FilteredMassSpecData(original, 1.0, nothing, target_ms_order=UInt8(2))
        for i in 1:length(filtered_ms2)
            @test getMsOrder(filtered_ms2, i) == UInt8(2)
        end
        
        # Test no MS order filtering
        filtered_all = FilteredMassSpecData(original, 1.0, nothing, target_ms_order=nothing)
        # Should have all scans
        @test length(filtered_all) == length(original)
    end
    
    @testset "Incremental Append" begin
        # Start with small sample
        filtered = FilteredMassSpecData(original, 0.1, nothing, seed=42)
        initial_length = length(filtered)
        initial_scans = copy(getOriginalScanIndices(filtered))
        
        # Append more
        n_added = append!(filtered, 0.1)
        @test length(filtered) == initial_length + n_added
        
        # Verify no duplicates
        all_indices = getOriginalScanIndices(filtered)
        @test length(unique(all_indices)) == length(all_indices)
        
        # Verify original scans unchanged
        @test all_indices[1:initial_length] == initial_scans
        
        # Verify new scans have correct data
        for i in (initial_length + 1):length(filtered)
            orig_idx = getOriginalScanIndex(filtered, i)
            @test getScanNumber(filtered, i) == getScanNumber(original, orig_idx)
        end
    end
    
    @testset "Zero Overhead Mode" begin
        # Test with no filtering or sampling
        filtered = FilteredMassSpecData(original, 1.0, nothing, target_ms_order=UInt8(2))
        
        # Should have same MS2 scans
        ms2_count = count(i -> getMsOrder(original, i) == UInt8(2), 1:length(original))
        @test length(filtered) == ms2_count
        
        # Verify data is copied correctly
        for i in 1:min(10, length(filtered))
            orig_idx = getOriginalScanIndex(filtered, i)
            orig_mz = getMzArray(original, orig_idx)
            filt_mz = getMzArray(filtered, i)
            @test length(orig_mz) == length(filt_mz)
            @test all(abs.(orig_mz .- filt_mz) .< 1e-6)
        end
    end
    
    @testset "Edge Cases" begin
        # Empty result
        filtered = FilteredMassSpecData(original, 0.0, nothing)
        @test length(filtered) == 0
        
        # TopN larger than peaks
        filtered = FilteredMassSpecData(original, 1.0, 10000)
        @test length(filtered) > 0
        for i in 1:min(5, length(filtered))
            orig_idx = getOriginalScanIndex(filtered, i)
            @test length(getMzArray(filtered, i)) == length(getMzArray(original, orig_idx))
        end
        
        # Append when fully sampled
        filtered = FilteredMassSpecData(original, 1.0, nothing, target_ms_order=UInt8(2))
        n_added = append!(filtered, 1.0)
        @test n_added == 0
    end
    
    @testset "Intensity Filtering" begin
        # Test with minimum intensity threshold
        min_intensity = 1000.0
        filtered = FilteredMassSpecData(original, 1.0, 50, min_intensity=min_intensity)
        
        for i in 1:min(5, length(filtered))
            intensities = getIntensityArray(filtered, i)
            @test all(intensities .>= min_intensity)
        end
    end
    
    @testset "Reproducibility" begin
        # Test that same seed gives same results
        filtered1 = FilteredMassSpecData(original, 0.5, 100, seed=123)
        filtered2 = FilteredMassSpecData(original, 0.5, 100, seed=123)
        
        @test length(filtered1) == length(filtered2)
        @test getOriginalScanIndices(filtered1) == getOriginalScanIndices(filtered2)
        
        # Different seed should give different results
        filtered3 = FilteredMassSpecData(original, 0.5, 100, seed=456)
        @test getOriginalScanIndices(filtered1) != getOriginalScanIndices(filtered3)
    end
end