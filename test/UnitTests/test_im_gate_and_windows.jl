# Tests for the ion-mobility helpers: the fragment-index gate (im_gate_lines / passes_im_gate) and the 1/K0 ->
# IM-scan conversion of the chromatogram window and integration band (im_half_width_scans).

using Test
using Pioneer
using Pioneer: ImGate, passes_im_gate, im_gate_lines, im_gate_sigma_mult, IM_GATE_TOL_SIGMA,
               im_half_width_scans, CHROM_IM_WINDOW_K0, CHROM_IM_BAND_K0

@testset "im_gate_lines — the z2 line for every charge, sigma scaled per charge" begin
    z2 = (1.49f0, -0.00095f0, 0.013f0)
    lines = im_gate_lines(z2)
    @test length(lines) == 9                                  # charges 0-8, lines[z + 1]
    @test lines[3] == z2                                      # z2 multiplier is 1
    @test lines[4] == (z2[1], z2[2], z2[3] * 1.8f0)           # z3
    @test lines[5] == (z2[1], z2[2], z2[3] * 2.2f0)           # z4 and above
    @test all(l -> l[1:2] == z2[1:2], lines)                  # same line, only sigma differs
    @test im_gate_sigma_mult(1) == im_gate_sigma_mult(6) == 2.2f0
end

@testset "passes_im_gate — tolerance edge, charge fallback, no gate" begin
    a, b, s = 1.5f0, -0.001f0, 0.01f0
    im_lib = Float32[1.0f0, 1.0f0, 1.0f0]
    charges = UInt8[2, 3, 12]                                 # 12 is beyond the line vector
    g = ImGate(im_gate_lines((a, b, s)), im_lib, charges, IM_GATE_TOL_SIGMA)
    scan_on_line = 500f0                                      # a + b * 500 = 1.0: exactly on the line
    @test all(pid -> passes_im_gate(g, UInt32(pid), scan_on_line), 1:3)
    # z2: tolerance is 4 x 0.01 = 0.04 1/K0 = 40 scans either side
    @test passes_im_gate(g, UInt32(1), 530f0)
    @test !passes_im_gate(g, UInt32(1), 550f0)
    # z3: 1.8x wider (0.072 = 72 scans), so 60 scans off passes for z3 but not z2
    @test passes_im_gate(g, UInt32(2), 560f0) && !passes_im_gate(g, UInt32(1), 560f0)
    # charge beyond the vector falls back to lines[1] (z0 entry, 2.2x sigma)
    @test passes_im_gate(g, UInt32(3), 580f0) && !passes_im_gate(g, UInt32(3), 600f0)
    @test passes_im_gate(nothing, UInt32(1), 0f0)
end

@testset "im_half_width_scans — 1/K0 half-width to whole IM scans, rounded up" begin
    # the two instrument slopes of the HYE runs (Ultra 2 50 ng, Ultra 250 pg)
    @test im_half_width_scans(CHROM_IM_WINDOW_K0, 0.0008653846f0) == 64   # 63.56 -> 64
    @test im_half_width_scans(CHROM_IM_WINDOW_K0, 0.0008499475f0) == 65   # 64.71 -> 65
    @test im_half_width_scans(CHROM_IM_BAND_K0, 0.0008653846f0) == 25     # 24.27 -> 25
    @test im_half_width_scans(CHROM_IM_BAND_K0, 0.0008499475f0) == 25     # 24.71 -> 25
    # a narrower mobility range packs more scans per 1/K0: 0.85-1.30 over 930 scans
    @test im_half_width_scans(CHROM_IM_WINDOW_K0, Float32(0.45 / 930)) == 114
    @test im_half_width_scans(CHROM_IM_WINDOW_K0, 0f0) == 0               # no slope -> off
end
