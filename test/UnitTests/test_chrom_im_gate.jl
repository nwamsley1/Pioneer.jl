# Tests for the ion-mobility gate helpers used by chromatogram extraction on packet data.

using Test
using Pioneer
using Pioneer: im_lines_by_charge, chrom_im_tol_sigma, chrom_rt_tol_override, CHROM_IM_TOL_SIGMA

@testset "im_lines_by_charge — charge-indexed lines with pooled fill" begin
    @test isempty(im_lines_by_charge(Dict{Int, NTuple{3, Float32}}()))
    model = Dict{Int, NTuple{3, Float32}}(
        0 => (1.47f0, -0.00093f0, 0.015f0),
        2 => (1.49f0, -0.00095f0, 0.013f0),
        3 => (1.36f0, -0.00077f0, 0.024f0),
    )
    lines = im_lines_by_charge(model)
    @test length(lines) == 8
    @test lines[2] == model[2] && lines[3] == model[3]
    @test lines[1] == model[0] && lines[4] == model[0] && lines[8] == model[0]
    # a charge above 8 extends the vector
    model[10] = (1.2f0, -0.0006f0, 0.04f0)
    @test length(im_lines_by_charge(model)) == 10
    # no pooled entry: the first line stands in
    only3 = Dict{Int, NTuple{3, Float32}}(3 => (1.36f0, -0.00077f0, 0.024f0))
    @test all(l -> l == only3[3], im_lines_by_charge(only3))
end

@testset "chromatogram IM / RT tolerance env overrides" begin
    withenv("PIONEER_CHROM_IM_SIGMA" => nothing, "PIONEER_CHROM_RT_TOL" => nothing) do
        @test chrom_im_tol_sigma() == CHROM_IM_TOL_SIGMA
        @test chrom_rt_tol_override() == 0f0
    end
    withenv("PIONEER_CHROM_IM_SIGMA" => "4.5", "PIONEER_CHROM_RT_TOL" => "0.15") do
        @test chrom_im_tol_sigma() == 4.5f0
        @test chrom_rt_tol_override() == 0.15f0
    end
    withenv("PIONEER_CHROM_IM_SIGMA" => "abc") do
        @test chrom_im_tol_sigma() == CHROM_IM_TOL_SIGMA
    end
end
