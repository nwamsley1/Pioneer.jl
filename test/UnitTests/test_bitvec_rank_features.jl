using Test
using Pioneer

@testset "bitvec excess rank table ranks fragment masks" begin
    target_counts = fill(10, 256)
    decoy_counts = fill(10, 256)
    target_counts[0x05 + 1] = 1000
    decoy_counts[0x05 + 1] = 1
    target_counts[0x03 + 1] = 800
    decoy_counts[0x03 + 1] = 1

    ranks = Pioneer._bitvec_excess_rank_table(target_counts, decoy_counts, 1.0)

    @test length(ranks) == 256
    @test ranks[0x05 + 1] == UInt16(1)
    @test ranks[0x03 + 1] == UInt16(2)
    @test Pioneer._bitvec_pattern_rank(ranks, UInt16(0x05)) == UInt16(1)
    @test Pioneer._bitvec_pattern_rank(nothing, UInt16(0x05)) == UInt16(0)
end

@testset "bitvec calibration keeps MS2 priority order" begin
    priority = Int32[9, 4, 7, 2]

    @test Pioneer._bitvec_calibration_scan_order(priority) == priority
end
