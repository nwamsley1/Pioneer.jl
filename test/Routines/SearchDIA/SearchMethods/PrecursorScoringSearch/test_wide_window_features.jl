using Test
using Pioneer
using DataFrames

const TEST_FLANKING_WINDOW_FEATURES = [
    :flanking_ms1_m0_candidate_fraction,
    :flanking_frag_candidate_fraction,
    :flanking_ms1_frag_sum_corr,
    :flanking_frag_corr_mean,
    :flanking_n_correlated_fragments,
    :flanking_n_correlated_fragments_bitvec_rank,
    :flanking_frag_corr_strength,
    :flanking_frag_corr_effective_n,
    :flanking_frag_corr_best_m0,
    :flanking_signal_support,
    :flanking_core_n_scans,
]

@testset "flanking-window features are cross-run only" begin
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, TEST_FLANKING_WINDOW_FEATURES)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), TEST_FLANKING_WINDOW_FEATURES)
    @test collect(Pioneer.FLANKING_WINDOW_FEATURES) == TEST_FLANKING_WINDOW_FEATURES
end

@testset "flanking-core bounds use contiguous passing cycles around the best score" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[101, 102, 103, 104],
        cycle_idx = UInt32[1, 2, 3, 4],
        lgbm_score = Float32[0.1, 0.9, 0.8, 0.2],
        weight = ones(Float32, 4),
        irt_obs = Float32[1, 2, 3, 4],
        ms1_m0_intensity = Float32[1, 10, 20, 40],
        frag1_int = Float32[0, 3, 5, 0],
        frag2_int = Float32[0, 1, 2, 0],
        frag3_int = zeros(Float32, 4),
        frag4_int = zeros(Float32, 4),
        frag5_int = zeros(Float32, 4),
        frag6_int = zeros(Float32, 4),
        frag7_int = zeros(Float32, 4),
        frag8_int = zeros(Float32, 4),
    )
    best = Pioneer.select_best_per_precursor!(psms, :lgbm_score)
    rank_table = fill(UInt16(99), 256)
    rank_table[0x00 + 1] = UInt16(1)

    Pioneer.add_trace_and_fragment_features!(
        best,
        psms,
        Bool[false, true, true, false];
        bitvec_rank_table = rank_table,
    )

    @test hasproperty(best, :flanking_core_scan_min)
    @test hasproperty(best, :flanking_core_scan_max)
    @test hasproperty(best, :flanking_core_n_scans)
    @test hasproperty(best, :flanking_core_ms1_m0_signal)
    @test hasproperty(best, :flanking_core_frag_signal)
    @test best.flanking_core_scan_min == UInt32[102]
    @test best.flanking_core_scan_max == UInt32[103]
    @test best.flanking_core_n_scans == UInt16[2]
    @test best.flanking_core_ms1_m0_signal == Float32[30]
    @test best.flanking_core_frag_signal == Float32[11]
end

@testset "flanking-window scans exclude core scans" begin
    key = (Int32(50000), Int32(200))
    scans = Int32[100, 105, 110, 115, 120]
    scan_to_window_key = fill((Int32(0), Int32(0)), 120)
    scan_to_window_pos = zeros(Int32, 120)
    for (pos, scan) in pairs(scans)
        scan_to_window_key[Int(scan)] = key
        scan_to_window_pos[Int(scan)] = Int32(pos)
    end
    scan_index = (
        window_scans = Dict(key => scans),
        scan_to_window_key = scan_to_window_key,
        scan_to_window_pos = scan_to_window_pos,
    )

    @test Pioneer._wide_flank_window_scans(Int32[105, 110], scan_index) == Int32[100, 115, 120]
end

@testset "flanking-window feature values combine core signal with flank support" begin
    flank_ms1_m0 = Float32[4, 1]
    flank_fragments = zeros(Float32, 2, 8)
    flank_fragments[:, 1] .= Float32[3, 2]
    flank_fragments[:, 2] .= Float32[1, 0.5]
    rank_table = fill(UInt16(99), 256)
    rank_table[0x03 + 1] = UInt16(7)

    features = Pioneer._wide_flank_feature_values(
        40f0,
        28f0,
        flank_ms1_m0,
        flank_fragments;
        bitvec_rank_table = rank_table,
    )

    @test features.flanking_ms1_m0_candidate_fraction ≈ Float32(40 / 45)
    @test features.flanking_frag_candidate_fraction ≈ Float32(28 / 34.5)
    @test features.flanking_ms1_frag_sum_corr > 0.9f0
    @test features.flanking_frag_corr_mean > 0.99f0
    @test features.flanking_n_correlated_fragments == UInt8(2)
    @test features.flanking_n_correlated_fragments_bitvec_rank == UInt16(7)
    @test features.flanking_frag_corr_strength > 1.9f0
    @test features.flanking_frag_corr_effective_n > 1.99f0
    @test features.flanking_frag_corr_best_m0 > 0.9f0
    @test features.flanking_signal_support == 1f0
end
