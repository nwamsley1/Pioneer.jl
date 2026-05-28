using Test
using Pioneer

@testset "natural isotope mass spacing" begin
    expected = 13.00335483507 - 12.0
    @test Pioneer.C13_C12_MASS_DIFF ≈ expected atol=1e-12
    @test Pioneer.C13_C12_MASS_DIFF_F32 === Float32(expected)
end
