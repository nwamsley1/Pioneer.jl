using Test
using Random
using Pioneer

struct ScanBoundTestModel <: Pioneer.AbstractMassErrorModel end

function Pioneer.getCorrectedMzAndBounds(::ScanBoundTestModel, mz::Float32,
                                        intensity::Float32, ::Float32)
    return mz, prevfloat(mz - intensity), nextfloat(mz + intensity)
end

@testset "Scan-derived lookup covers individual intervals" begin
    corrected = Float32[]
    lows = Float32[]
    highs = Float32[]
    width = Ref(0f0)
    masses = Union{Missing, Float32}[100, 200, 300, missing]
    intensities = Union{Missing, Float32}[0.001, 0.020, 0.010, 1]
    n = Pioneer.prepare_scan_peaks!(corrected, lows, highs, ScanBoundTestModel(),
                                    masses, intensities, 0f0, width)
    @test n == 4
    @test isfinite(width[])
    @test corrected[4] == Inf32
    for j in 1:3, target in (lows[j], highs[j], corrected[j])
        lo, hi = Pioneer.scan_match_window(target, width[])
        @test lo <= corrected[j] <= hi
    end

    target = 200.018f0
    @test abs(corrected[2] - target) > 0.015f0
    lo, hi = Pioneer.scan_match_window(target, width[])
    start = Pioneer.bsearch_hybrid(corrected, lo, 1, 3)
    @test Pioneer.scan_for_nearest_in_window(corrected, lows, highs,
                                            start, 3, target, hi)[1] == 2

    # Omitting the reduction preserves the calibration preparation API and arrays.
    other_c, other_l, other_h = Float32[], Float32[], Float32[]
    @test Pioneer.prepare_scan_peaks!(other_c, other_l, other_h, ScanBoundTestModel(),
                                     masses, intensities, 0f0) == n
    @test isequal((other_c, other_l, other_h), (corrected, lows, highs))
    @test Pioneer.prepare_scan_peaks!(corrected, lows, highs, ScanBoundTestModel(),
                                     Float32[], Float32[], 0f0, width) == 0
    @test width[] == 0f0

    rng = MersenneTwister(9382)
    masses = sort!(Float32.(150 .+ 1850 .* rand(rng, 2000)))
    intensities = Float32.(0.0001 .+ 0.04 .* rand(rng, 2000))
    n = Pioneer.prepare_scan_peaks!(corrected, lows, highs, ScanBoundTestModel(),
                                    masses, intensities, 0f0, width)
    targets = vcat(copy(lows[1:n]), copy(highs[1:n]), masses)
    for target in targets
        eligible = findall(j -> lows[j] <= target <= highs[j], 1:n)
        expected = isempty(eligible) ? 0 : eligible[argmin(abs.(corrected[eligible] .- target))]
        lo, hi = Pioneer.scan_match_window(target, width[])
        start = Pioneer.bsearch_hybrid(corrected, lo, 1, n)
        found = Pioneer.scan_for_nearest_in_window(corrected, lows, highs,
                                                   start, n, target, hi)[1]
        @test found == expected
    end
end
