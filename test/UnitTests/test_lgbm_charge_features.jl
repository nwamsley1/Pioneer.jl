using Test
using DataFrames
using Pioneer

@testset "LightGBM charge-state features" begin
    @test :charge in Pioneer.PRESCORE_FEATURES
    @test :charge in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(DataFrame(charge = UInt8[2, 3, 4]), [:charge])
    @test X[:, 1] == Float32[2, 3, 4]
end

@testset "LightGBM spectral contrast feature" begin
    @test :spectral_contrast in Pioneer.PRESCORE_FEATURES
    @test :spectral_contrast in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(
        DataFrame(spectral_contrast = Float32[0.0, 0.5, 1.0]),
        [:spectral_contrast],
    )
    @test X[:, 1] == Float32[0.0, 0.5, 1.0]
end
