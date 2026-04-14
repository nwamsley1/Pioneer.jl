using Test
using Pioneer

@testset "PrecursorScoringSearch model configurations" begin
    configs = Pioneer.create_model_configurations()

    # Expect five models after addition of simplified probit
    @test length(configs) == 5

    names = getfield.(configs, :name)
    @test "ProbitRegression" in names
    @test "ProbitRegressionSimple" in names

    # Intercept should appear only in probit model feature sets
    model_by_name = Dict(name => cfg for (name, cfg) in zip(names, configs))
    @test :intercept in model_by_name["ProbitRegression"].features
    @test :intercept in model_by_name["ProbitRegressionSimple"].features
    @test :intercept ∉ model_by_name["SimpleLightGBM"].features
    @test :intercept ∉ model_by_name["AdvancedLightGBM"].features
    @test :intercept ∉ model_by_name["SuperSimplified"].features

    complementary_features = [
        :complementary_pair_fraction,
        :cleavage_coverage,
        :hypergeom_match_score
    ]
    for feature in complementary_features
        @test feature in model_by_name["SimpleLightGBM"].features
        @test feature in model_by_name["AdvancedLightGBM"].features
        @test feature in model_by_name["ProbitRegression"].features
    end
    for name in ["SimpleLightGBM", "AdvancedLightGBM", "ProbitRegression", "ProbitRegressionSimple", "SuperSimplified"]
        @test :target_protection_spectral_angle in model_by_name[name].features
        @test :target_neighbor_spectral_angle ∉ model_by_name[name].features
        @test :target_protection_pair_decoy ∉ model_by_name[name].features
    end
    @test :complementary_pair_count ∉ model_by_name["SimpleLightGBM"].features
    @test :complementary_pair_count ∉ model_by_name["AdvancedLightGBM"].features
    @test :complementary_pair_count ∉ model_by_name["ProbitRegression"].features

end

@testset "Complex PSM complementary features" begin
    matches = [
        Pioneer.FragmentMatch{Float32}(1.0f0, 10.0f0, 100.0f0, 99.999f0, 1, UInt8(3), UInt8(1), UInt8(0), UInt8(1), false, UInt32(1), UInt8(1), UInt32(1), UInt32(1), UInt8(1)),
        Pioneer.FragmentMatch{Float32}(1.0f0, 20.0f0, 200.0f0, 200.004f0, 2, UInt8(4), UInt8(1), UInt8(0), UInt8(2), false, UInt32(1), UInt8(1), UInt32(1), UInt32(1), UInt8(1)),
        Pioneer.FragmentMatch{Float32}(1.0f0, 30.0f0, 300.0f0, 299.991f0, 3, UInt8(2), UInt8(1), UInt8(0), UInt8(1), false, UInt32(1), UInt8(1), UInt32(1), UInt32(1), UInt8(1)),
        Pioneer.FragmentMatch{Float32}(1.0f0, 40.0f0, 400.0f0, 400.008f0, 4, UInt8(2), UInt8(1), UInt8(0), UInt8(2), false, UInt32(1), UInt8(1), UInt32(1), UInt32(1), UInt8(1))
    ]
    id_to_col = Pioneer.ArrayDict(UInt32, UInt16, 10)
    Pioneer.update!(id_to_col, UInt32(1), UInt16(1))
    sequence_lengths = UInt8[7]

    pair_fractions, cleavage_coverages = Pioneer.getComplementaryIonFeatures(
        matches,
        id_to_col,
        sequence_lengths,
        length(matches),
        1
    )

    @test isapprox(Float32(pair_fractions[1]), 1.0f0 / 6.0f0; atol=1e-3)
    @test isapprox(Float32(cleavage_coverages[1]), 3.0f0 / 6.0f0; atol=1e-3)

    unscored = [Pioneer.ComplexUnscoredPSM{Float32}()]
    Pioneer.ScoreFragmentMatches!(
        unscored,
        id_to_col,
        matches,
        length(matches),
        Pioneer.MassErrorModel(0.0f0, (10.0f0, 10.0f0)),
        5
    )

    scored = Vector{Pioneer.ComplexScoredPSM{Float32, Float16}}(undef, 1)
    spectral_scores = [
        Pioneer.SpectralScoresComplex{Float16}(
            Float16(1), Float16(1), Float16(1), Float16(1), Float16(1),
            Float16(1), Float16(1), Float16(1), Float16(0)
        )
    ]

    last_val = Pioneer.Score!(
        scored,
        unscored,
        spectral_scores,
        Float32[1.0],
        id_to_col,
        1,
        0.5,
        0,
        1,
        1000.0f0,
        42,
        matches,
        length(matches),
        sequence_lengths;
        min_y_count=1,
        min_frag_count=1,
        max_best_rank=1,
        min_topn=1,
        target_protection_spectral_angle=Float16[0.25],
        target_protection_pair_decoy=Bool[true],
        target_neighbor_spectral_angle=Float16[0.75]
    )

    @test last_val == 1
    @test isapprox(Float32(scored[1].complementary_pair_fraction), 1.0f0 / 6.0f0; atol=1e-3)
    @test isapprox(Float32(scored[1].cleavage_coverage), 3.0f0 / 6.0f0; atol=1e-3)
    @test isapprox(Float32(scored[1].target_protection_spectral_angle), 0.25f0; atol=1e-3)
    @test scored[1].target_protection_pair_decoy
    @test isapprox(Float32(scored[1].target_neighbor_spectral_angle), 0.75f0; atol=1e-3)
end

@testset "Target protection spectral angle feature" begin
    H = Pioneer.SparseArray(Int64(10))
    H.n_vals = 5
    H.m = 3
    H.n = 4
    H.rowval[1:5] .= [1, 1, 2, 3, 3]
    H.nzval[1:5] .= Float32[1.0, 1.0, sqrt(3.0f0), 1.0, 1.0]
    H.colptr[1:5] .= [1, 2, 4, 5, 6]

    id_to_col = Pioneer.ArrayDict(UInt32, UInt16, 10)
    Pioneer.update!(id_to_col, UInt32(1), UInt16(1))
    Pioneer.update!(id_to_col, UInt32(2), UInt16(2))
    Pioneer.update!(id_to_col, UInt32(3), UInt16(3))
    Pioneer.update!(id_to_col, UInt32(4), UInt16(4))

    neighbor_angles = Pioneer.getMaxNeighborSpectralAngles(
        H,
        Bool[false, true, false, true],
        UInt32[1, 1, 2, 3],
        id_to_col
    )
    spectral_angles = neighbor_angles.target_protection_spectral_angle
    pair_decoys = neighbor_angles.target_protection_pair_decoy
    target_neighbor_angles = neighbor_angles.target_neighbor_spectral_angle

    expected_half_cosine_angle = 1.0f0 - (2.0f0 / Float32(pi)) * acos(0.5f0)
    @test isapprox(Float32(spectral_angles[1]), expected_half_cosine_angle; atol=1e-3)
    @test isapprox(Float32(spectral_angles[2]), 0.0f0; atol=1e-3)
    @test isapprox(Float32(spectral_angles[3]), 1.0f0; atol=1e-3)
    @test isapprox(Float32(spectral_angles[4]), 0.0f0; atol=1e-3)
    @test pair_decoys == Bool[true, false, false, false]
    @test isapprox(Float32(target_neighbor_angles[1]), 0.0f0; atol=1e-3)
    @test isapprox(Float32(target_neighbor_angles[2]), expected_half_cosine_angle; atol=1e-3)
    @test isapprox(Float32(target_neighbor_angles[3]), 0.0f0; atol=1e-3)
    @test isapprox(Float32(target_neighbor_angles[4]), 1.0f0; atol=1e-3)
end

@testset "Hypergeometric match score" begin
    expected = -log10(8 / 120)

    @test Pioneer.hypergeom_tail_neglog10(0, 2, 10, 3) == Float16(0)
    @test isapprox(Float32(Pioneer.hypergeom_tail_neglog10(2, 2, 10, 3)), expected; atol=2e-3)
    @test Pioneer.hypergeom_tail_neglog10(2, 2, 10, 3) > Pioneer.hypergeom_tail_neglog10(1, 2, 10, 3)

    score = Pioneer.get_hypergeom_match_score(
        2,
        2,
        3,
        100.0f0,
        100.05f0,
        100.0f0,
        Pioneer.MassErrorModel(0.0f0, (25.0f0, 25.0f0))
    )
    @test isapprox(Float32(score), expected; atol=2e-3)
end
