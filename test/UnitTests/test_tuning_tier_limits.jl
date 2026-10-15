# ParameterTuningSearch collection limits: the per-tier scan cap and the low-tier peak limit.

using Test
using Pioneer: tier_scan_cap, tier_max_peaks, projected_short, TUNING_TIER_SCAN_CAP_MIN,
               TUNING_TIER_SCAN_CAP_FRACTION, TUNING_LOW_TIER_MAX_SCORE, TUNING_PROJECTION_STOP_FRACTION

@testset "tuning tier limits" begin
    @testset "tier_scan_cap" begin
        @test tier_scan_cap(5_000) == 5_000                                   # never more than the file has
        @test tier_scan_cap(46_662) == TUNING_TIER_SCAN_CAP_MIN               # floor for small files
        @test tier_scan_cap(835_971) == ceil(Int, TUNING_TIER_SCAN_CAP_FRACTION * 835_971)   # fraction for large ones
        @test tier_scan_cap(0) == 0
    end
    @testset "tier_max_peaks" begin
        hi = Int(TUNING_LOW_TIER_MAX_SCORE) + 1; lo = Int(TUNING_LOW_TIER_MAX_SCORE)
        @test tier_max_peaks(hi, 0, 200) == 0          # high tiers keep every peak
        @test tier_max_peaks(lo, 0, 200) == 200        # low tiers limited
        @test tier_max_peaks(lo - 1, 0, 200) == 200
        @test tier_max_peaks(lo, 100, 200) == 100      # never loosens a tighter existing limit
        @test tier_max_peaks(lo, 300, 200) == 200
        @test tier_max_peaks(lo, 0, 0) == 0            # disabled
    end
    @testset "projected_short" begin
        cap = 83_598; target = 1_200
        @test projected_short(3, 5_000, cap, target)          # single-cell score>=8: ~50 projected
        @test !projected_short(450, 5_000, cap, target)       # healthy yield
        @test !projected_short(0, 0, 0, 0)                    # nothing to project, no target
        # the threshold itself: just below / at TUNING_PROJECTION_STOP_FRACTION of the target
        n = 10_000; f = TUNING_PROJECTION_STOP_FRACTION
        @test projected_short(floor(Int, f * target * n / cap) - 1, n, cap, target)
        @test !projected_short(ceil(Int, f * target * n / cap) + 1, n, cap, target)
    end
end
