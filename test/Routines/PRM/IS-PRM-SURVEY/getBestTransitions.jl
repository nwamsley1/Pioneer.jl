@testset "getBestTramsitions.jl" begin

test_best_psm = (rt = 10,
    scan_idx = 10,
    name = ["y3+2","y3+1", "y4+1"],
    mz = Float64[200, 100, 200],
    intensity = Float64[1000, 2000, 500]
    )

@test getBestTransitions(test_best_psm, 1) == [2]
@test getBestTransitions(test_best_psm, 2) == [2, 3]


test_best_psm = (rt = 10,
    scan_idx = 10,
    name = ["y3+2", "y3+1", "y4+2", "y4+1", "y5+1"],
    mz = Float64[200, 100, 200],
    intensity = Float64[1000, 2000, 500, 500, 500]
    )

    @test getBestTransitions(test_best_psm, 3) == [2, 4, 5]
    @test getBestTransitions(test_best_psm, 2) == [2, 4]
    @test getBestTransitions(test_best_psm, 1000) == [2, 4, 5]


end