using Test
using Pioneer

@testset "summarize_precursor tests" begin
    @testset "Empty input handling" begin
        # Test that empty inputs return all missing values
        result = Pioneer.summarize_precursor(
            UInt8[],      # iso_idx
            Float32[],    # center_mz
            Float32[],    # iso_mz
            UInt8[],      # prec_charge
            Float32[],    # weight
            Float32[]     # δ
        )
        
        @test ismissing(result.center_mz)
        @test ismissing(result.δ)
        @test ismissing(result.yt)
        @test ismissing(result.x0)
        @test ismissing(result.x1)
        @test ismissing(result.prec_charge)
        
        # Verify it's the expected NamedTuple structure
        @test isa(result, NamedTuple)
        @test keys(result) == (:center_mz, :δ, :yt, :x0, :x1, :prec_charge)
    end
    
    @testset "Single M0 isotope" begin
        # Test with only M0 isotope (iso_idx = 1)
        result = Pioneer.summarize_precursor(
            UInt8[1],         # iso_idx (M0)
            Float32[500.0],   # center_mz
            Float32[500.1],   # iso_mz
            UInt8[2],         # prec_charge
            Float32[1000.0],  # weight
            Float32[0.9]      # δ
        )
        
        @test result.center_mz == 500.0f0
        @test result.δ == 0.9f0
        @test result.yt == 10.0f0  # Expected extreme ratio
        @test result.x0 ≈ 0.1f0
        @test result.prec_charge == 2
        # x1 should be m1_mz - center_mz where m1_mz = iso_mz + NEUTRON/charge
    end
    
    @testset "Single M1 isotope" begin
        # Test with only M1 isotope (iso_idx = 2)
        result = Pioneer.summarize_precursor(
            UInt8[2],         # iso_idx (M1)
            Float32[500.0],   # center_mz
            Float32[500.5],   # iso_mz
            UInt8[2],         # prec_charge
            Float32[800.0],   # weight
            Float32[0.8]      # δ
        )
        
        @test result.center_mz == 500.0f0
        @test result.δ == 0.8f0
        @test result.yt == -10.0f0  # Expected extreme ratio
        @test result.x1 ≈ 0.5f0
        @test result.prec_charge == 2
    end
    
    @testset "M0/M1 pair" begin
        # Test with both M0 and M1 isotopes
        result = Pioneer.summarize_precursor(
            UInt8[1, 2],              # iso_idx (M0, M1)
            Float32[500.0, 500.0],    # center_mz
            Float32[500.1, 500.6],    # iso_mz
            UInt8[2, 2],              # prec_charge
            Float32[1000.0, 100.0],   # weight
            Float32[0.9, 0.1]         # δ
        )
        
        @test result.center_mz == 500.0f0
        @test result.δ == 0.9f0
        @test result.prec_charge == 2
        @test result.x0 ≈ 0.1f0
        @test result.x1 ≈ 0.6f0
        # yt should be log(weight[m0]/(weight[m1]*δ[m0]))
        expected_yt = log(1000.0f0 / (100.0f0 * 0.9f0))
        @test result.yt ≈ expected_yt
    end
    
    @testset "Invalid patterns return missing" begin
        # Test with 3 isotopes (should return all missing)
        result = Pioneer.summarize_precursor(
            UInt8[1, 2, 3],                    # iso_idx
            Float32[500.0, 500.0, 500.0],      # center_mz
            Float32[500.1, 500.6, 501.1],      # iso_mz
            UInt8[2, 2, 2],                    # prec_charge
            Float32[1000.0, 100.0, 10.0],      # weight
            Float32[0.9, 0.08, 0.02]           # δ
        )
        
        @test ismissing(result.center_mz)
        @test ismissing(result.δ)
        @test ismissing(result.yt)
        @test ismissing(result.x0)
        @test ismissing(result.x1)
        @test ismissing(result.prec_charge)
    end
end