using Test
using DataFrames

using Pioneer: Chromatogram, integrate_chrom

@testset "IntegrateChrom basic" begin
    # Construct a tiny synthetic chromatogram with a single peak
    N = 21
    rt = collect(Float32, 1:N)
    # Triangular peak centered at apex
    apex = Int(ceil(N/2))
    intensity = Float32.(max.(0, (N÷2) .- abs.((1:N) .- apex)))
    frac = ones(Float32, N)  # precursor_fraction_transmitted
    scan_idx = collect(Int32, 1:N)

    chrom = DataFrame(rt=rt, intensity=intensity, precursor_fraction_transmitted=frac, scan_idx=scan_idx)

    # Buffers/state
    b = zeros(Float32, N)
    u2 = zeros(Float32, N)
    state = Chromatogram(zeros(Float32, N), zeros(Float32, N), 0)

    avg_cycle_time = 1.0f0
    λ = 1.0f0
    area, apex_out, n_points = integrate_chrom(@view chrom[:, :], apex, b, u2, state, avg_cycle_time, λ;
                                              max_apex_offset=2, n_pad=0, isplot=false)

    @test area > 0
    @test n_points > 0
    @test isa(apex_out, Int32) || isa(apex_out, Int)
end

