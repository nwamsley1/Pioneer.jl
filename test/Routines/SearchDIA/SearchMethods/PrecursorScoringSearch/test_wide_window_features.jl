@testset "wide-window cross-run feature values" begin
    ms1_m0 = Float32[4, 10, 20, 10, 1]
    fragments = Float32[
        3 1 0 0 0 0
        5 2 0 0 0 0
        10 4 0 0 0 0
        5 2 0 0 0 0
        2 1 0 0 0 0
    ]
    candidate_mask = Bool[false, true, true, true, false]

    features = Pioneer._wide_window_feature_values(ms1_m0, fragments, candidate_mask)

    @test features.wide_ms1_m0_candidate_fraction ≈ Float32(40 / 45)
    @test features.wide_frag_candidate_fraction ≈ Float32(28 / 35)
    @test features.wide_ms1_frag_sum_corr > 0.9f0
    @test features.wide_frag_corr_mean > 0.99f0
    @test features.wide_n_correlated_fragments == UInt8(2)
    @test features.wide_frag_corr_best_m0 > 0.9f0
    @test features.wide_signal_support == 1f0
end

@testset "wide-window features handle missing raw evidence" begin
    features = Pioneer._wide_window_feature_values(
        zeros(Float32, 4),
        zeros(Float32, 4, 6),
        Bool[true, false, true, false],
    )

    for value in values(features)
        @test isfinite(Float32(value))
        @test value == zero(value)
    end
end

@testset "wide-window features are cross-run only" begin
    @test all(feature -> feature in Pioneer.ADVANCED_FEATURE_SET, Pioneer.WIDE_WINDOW_FEATURES)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), Pioneer.WIDE_WINDOW_FEATURES)
end
