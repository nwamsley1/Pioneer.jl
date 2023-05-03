@testset "getIntegrationBounds.jl" begin
    heavy_scans = Float64[1, 2, 3, 4, 5]
    light_scans = Float64[0, 1.4, 1.7, 2.4, 4.1, 5.5]
    
    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(1.0)) == ([1, 2, 2, 4, 5], [2, 3, 4, 5, 6])

    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(0.01)) == (Float64[], Float64[])

    light_scans = Float64[0, 1.4, 1.7, 2.0, 2.4, 4.1, 5.5]

    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(1.0)) == ([1, 2, 2, 2, 4, 5], [2, 3, 4, 5, 6, 7])

    heavy_scans = Float64[1.0, 1.1, 1.3, 1.5, 2, 3, 4, 5]
    light_scans = Float64[1.4, 1.7, 2.2, 3.2, 4.2, 5.3]

    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(1.0)) == ([4, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6])
    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(0.25)) == ([4, 4, 5, 6, 7], [1, 2, 3, 4, 5])
    @test getScanPairs(heavy_scans, light_scans, max_diff = Float64(0.2)) == ([4, 4], [1, 2])


end