using Test
using Pioneer

@testset "Intensity mass-error spline optimizer stability" begin
    n = 220_000
    log2_intensity = collect(Float64.(range(5, stop = 25, length = n)))
    residual_da = 1e-4 .* sin.(log2_intensity)

    spline, label = Pioneer.fit_best_monotone_bias(log2_intensity, residual_da;
        n_knots = 10,
        lr = 0.0001,
        max_iter = 5000)

    @test spline !== nothing
    @test label in ("monotone-increasing", "monotone-decreasing")
    @test all(isfinite, spline.(log2_intensity[1:1000:end]))
end
