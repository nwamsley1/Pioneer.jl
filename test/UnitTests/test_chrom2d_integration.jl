# integrate_chrom_2d: the 2D (cycle x ion-mobility) integrator used on slice data.
using Test, Statistics
using Pioneer
using Pioneer: integrate_chrom_2d, Chrom2DScratch, WHWorkspace, Chromatogram

"Build the flat (rt, scan_idx, cycle, im_scan, intensity) columns for a synthetic grid."
function grid_cols(nrow, ncol; stride = 8, peak = (r, c) -> 0.0f0, base = 0.0f0)
    rt = Float32[]; sidx = UInt32[]; cyc = UInt32[]; im = UInt16[]; w = Float32[]
    s = UInt32(0)
    for r in 1:nrow, c in 1:ncol          # row-major: all mobility scans of one cycle together
        s += UInt32(1)
        push!(rt, Float32(0.1 * r)); push!(sidx, s); push!(cyc, UInt32(r))
        push!(im, UInt16(100 + stride * c)); push!(w, Float32(peak(r, c) + base))
    end
    rt, sidx, cyc, im, w
end

run2d(cols, apex_row_pos, half) = begin
    sc = Chrom2DScratch(); ws = WHWorkspace(256)
    st = Chromatogram(zeros(Float32, 256), zeros(Float32, 256), 0)
    integrate_chrom_2d(cols[1], cols[2], cols[3], cols[4], cols[5], apex_row_pos, half, sc, ws, st, 1.0f-6)
end

@testset "integrate_chrom_2d" begin
    nrow, ncol = 15, 11
    gauss(r, c) = 1000 * exp(-((r - 8)^2) / (2 * 1.4^2) - ((c - 6)^2) / (2 * 1.2^2))

    @testset "recovers a clean peak and ignores the empty rim" begin
        cols = grid_cols(nrow, ncol; peak = gauss)
        apex_pos = findfirst(i -> cols[3][i] == 8 && cols[4][i] == UInt16(100 + 8 * 6), 1:length(cols[1]))
        area, apex_sidx, npts, lo, hi, withheld = run2d(cols, apex_pos, 24)
        total = sum(cols[5])
        @test !withheld
        @test area > 0.75 * total          # the band holds the great majority of a compact peak
        @test area <= total                 # and never more than was there
        @test cols[4][findfirst(==(apex_sidx), cols[2])] == UInt16(100 + 8 * 6)   # apex on the right slice
        @test npts >= 3 && lo < hi
    end

    @testset "flat baseline is removed" begin
        base = 20.0f0
        plain = grid_cols(nrow, ncol; peak = gauss)
        withbase = grid_cols(nrow, ncol; peak = gauss, base = base)
        apex_pos = findfirst(i -> plain[3][i] == 8 && plain[4][i] == UInt16(100 + 8 * 6), 1:length(plain[1]))
        a0 = run2d(plain, apex_pos, 24)[1]
        a1 = run2d(withbase, apex_pos, 24)[1]
        # a constant pedestal everywhere must be subtracted, not integrated
        @test isapprox(a1, a0; rtol = 0.05)
    end

    @testset "the band is in IM-scan units, not columns" begin
        cols = grid_cols(nrow, ncol; peak = gauss)          # slices 8 scans apart
        apex_pos = findfirst(i -> cols[3][i] == 8 && cols[4][i] == UInt16(100 + 8 * 6), 1:length(cols[1]))
        narrow = run2d(cols, apex_pos, 8)[1]                # +/- 1 slice
        wide   = run2d(cols, apex_pos, 32)[1]               # +/- 4 slices
        @test narrow < wide                                  # a wider band collects more of the peak
        @test narrow > 0
        # the same grid with slices 16 scans apart: the same 1/K0 band must now span FEWER slices
        cols16 = grid_cols(nrow, ncol; stride = 16, peak = gauss)
        apex16 = findfirst(i -> cols16[3][i] == 8 && cols16[4][i] == UInt16(100 + 16 * 6), 1:length(cols16[1]))
        @test run2d(cols16, apex16, 32)[1] < wide
    end

    @testset "apex is hill-climbed away from a bad seed" begin
        cols = grid_cols(nrow, ncol; peak = gauss)
        bad = findfirst(i -> cols[3][i] == 3 && cols[4][i] == UInt16(100 + 8 * 2), 1:length(cols[1]))
        area, apex_sidx, _, _, _, _ = run2d(cols, bad, 24)
        @test cols[4][findfirst(==(apex_sidx), cols[2])] == UInt16(100 + 8 * 6)
        @test cols[3][findfirst(==(apex_sidx), cols[2])] == UInt32(8)
        @test area > 0
    end

    @testset "degenerate input" begin
        empty_cols = (Float32[], UInt32[], UInt32[], UInt16[], Float32[])
        @test run2d(empty_cols, 1, 24)[1] == 0.0f0
        zeros_cols = grid_cols(nrow, ncol)                  # all zero
        @test run2d(zeros_cols, 1, 24)[1] == 0.0f0
        tiny = grid_cols(2, 3; peak = (r, c) -> 5.0f0)      # fewer than 3 cycles
        @test run2d(tiny, 1, 24)[1] == 0.0f0
    end
end
