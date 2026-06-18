using Test
using Pioneer
using DataFrames

const TEST_WIDE_WINDOW_FEATURES = [
    :wide_ms1_m0_candidate_fraction,
    :wide_frag_candidate_fraction,
    :wide_ms1_frag_sum_corr,
    :wide_frag_corr_mean,
    :wide_n_correlated_fragments,
    :wide_n_correlated_fragments_bitvec_rank,
    :wide_frag_corr_strength,
    :wide_frag_corr_effective_n,
    :wide_frag_corr_best_m0,
    :wide_signal_support,
    :wide_core_n_scans,
]

@testset "wide-window features are cross-run only" begin
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, TEST_WIDE_WINDOW_FEATURES)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), TEST_WIDE_WINDOW_FEATURES)
    @test collect(Pioneer.WIDE_WINDOW_FEATURES) == TEST_WIDE_WINDOW_FEATURES
end

@testset "wide-core bounds use contiguous passing cycles around the best score" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10, 10],
        scan_idx = UInt32[101, 102, 103, 104],
        cycle_idx = UInt32[1, 2, 3, 4],
        lgbm_score = Float32[0.1, 0.9, 0.8, 0.2],
        weight = ones(Float32, 4),
        irt_obs = Float32[1, 2, 3, 4],
        frag1_int = zeros(Float32, 4),
        frag2_int = zeros(Float32, 4),
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

    @test hasproperty(best, :wide_core_scan_min)
    @test hasproperty(best, :wide_core_scan_max)
    @test hasproperty(best, :wide_core_n_scans)
    @test best.wide_core_scan_min == UInt32[102]
    @test best.wide_core_scan_max == UInt32[103]
    @test best.wide_core_n_scans == UInt16[2]
end

@testset "wide-window feature values summarize raw candidate support" begin
    ms1_m0 = Float32[4, 10, 20, 10, 1]
    fragments = zeros(Float32, 5, 8)
    fragments[:, 1] .= Float32[3, 5, 10, 5, 2]
    fragments[:, 2] .= Float32[1, 2, 4, 2, 1]
    candidate_mask = Bool[false, true, true, true, false]
    rank_table = fill(UInt16(99), 256)
    rank_table[0x03 + 1] = UInt16(7)

    features = Pioneer._wide_window_feature_values(
        ms1_m0,
        fragments,
        candidate_mask;
        bitvec_rank_table = rank_table,
    )

    @test features.wide_ms1_m0_candidate_fraction ≈ Float32(40 / 45)
    @test features.wide_frag_candidate_fraction ≈ Float32(28 / 35)
    @test features.wide_ms1_frag_sum_corr > 0.9f0
    @test features.wide_frag_corr_mean > 0.99f0
    @test features.wide_n_correlated_fragments == UInt8(2)
    @test features.wide_n_correlated_fragments_bitvec_rank == UInt16(7)
    @test features.wide_frag_corr_strength > 1.9f0
    @test features.wide_frag_corr_effective_n > 1.99f0
    @test features.wide_frag_corr_best_m0 > 0.9f0
    @test features.wide_signal_support == 1f0
end
