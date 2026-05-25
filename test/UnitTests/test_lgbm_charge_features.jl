using Test
using DataFrames
using Pioneer

@testset "LightGBM charge-state features" begin
    @test :charge in Pioneer.PRESCORE_FEATURES
    @test :charge in Pioneer.ADVANCED_FEATURE_SET

    X = Pioneer.feature_matrix(DataFrame(charge = UInt8[2, 3, 4]), [:charge])
    @test X[:, 1] == Float32[2, 3, 4]
end
