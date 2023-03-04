using Titus
using Test
using Tables, Arrow
println(pwd())
include("../src/integrate.jl")
#table = Arrow.Table("../data/parquet/GAPDH_VGVNGFGR.arrow")
#Avoid compiliation time from first call
#@time table = Arrow.Table("../data/parquet/GAPDH_VGVNGFGR.arrow")

@testset "integrate.jl" begin

    test_indices_a = Vector{Int64}([1, 2, 10, 11, 12, 13, 60, 61, 62, 63, 64, 65, 66, 67, 68, 100])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (7, 15)
    test_indices_a = Vector{Int64}([1, 2, 3, 12, 13])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 5)
    test_indices_a = Vector{Int64}([1, 2, 3, 12, 100])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 3)
    test_indices_a = Vector{Int64}([1, 2, 3, 12])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 3)
    test_indices_a = Vector{Int64}([1, 2, 3, 12, 100, 101, 102, 103, 104])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (5, 9)
    test_indices_a = Vector{Int64}([1])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 1)
    test_indices_a = Vector{Int64}([1,2])
    @test getIntegrationBounds(test_indices_a, max_gap_size = 10) == (1, 2)

end