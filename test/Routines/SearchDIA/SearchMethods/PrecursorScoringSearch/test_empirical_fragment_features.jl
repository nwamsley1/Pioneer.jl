@testset "empirical fragment reference features use top-2 leave-one-out references" begin
    precursor_idx = UInt32[10, 10, 10, 10, 10, 10, 11]
    row_ids = UInt64[1, 2, 3, 4, 5, 6, 7]
    scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.7]
    peps = Float32[0.0, 0.75, 1.0, 1.0, 1.0, 0.0, 0.02]
    frag_cols = (
        Float32[100, 0, 0, 0, 0, 0, 25],
        Float32[0, 100, 100, 100, 100, 100, 25],
        Float32[0, 0, 0, 0, 0, 0, 0],
        Float32[0, 0, 0, 0, 0, 0, 0],
        Float32[0, 0, 0, 0, 0, 0, 0],
        Float32[0, 0, 0, 0, 0, 0, 0],
        Float32[0, 0, 0, 0, 0, 0, 0],
        Float32[0, 0, 0, 0, 0, 0, 0],
    )
    refs = Pioneer._empirical_fragment_topk_refs(
        precursor_idx,
        row_ids,
        scores,
        peps,
        frag_cols,
    )
    hellinger, ref_pep = Pioneer._empirical_fragment_reference_features(
        precursor_idx,
        row_ids,
        frag_cols,
        refs,
    )

    @test hellinger[1] ≈ 1.0f0
    @test ref_pep[1] ≈ 0.75f0
    @test hellinger[2] ≈ 1.0f0
    @test ref_pep[2] ≈ 0.0f0

    @test hellinger[6] ≈ 1.0f0
    @test ref_pep[6] ≈ 0.0f0

    @test hellinger[7] ≈ 1.0f0
    @test ref_pep[7] ≈ 1.0f0
end

@testset "empirical fragment features are cross-run only" begin
    @test :empirical_frag_best_hellinger in Pioneer.ADVANCED_FEATURE_SET
    @test !(:empirical_frag_ref_pep in Pioneer.ADVANCED_FEATURE_SET)
    @test all(feature -> !(feature in Pioneer.PRESCORE_FEATURES), Pioneer.EMPIRICAL_FRAGMENT_FEATURES)
end

@testset "empirical fragment references can use Huber shadow overrides" begin
    precursor_idx = UInt32[20, 20, 20]
    row_ids = UInt64[1, 2, 3]
    scores = Float32[0.9, 0.8, 0.7]
    peps = Float32[0.01, 0.02, 0.03]
    frag_cols = (
        Float32[100, 0, 0],
        Float32[0, 100, 100],
        Float32[0, 0, 0],
        Float32[0, 0, 0],
        Float32[0, 0, 0],
        Float32[0, 0, 0],
        Float32[0, 0, 0],
        Float32[0, 0, 0],
    )
    refs = Pioneer._empirical_fragment_topk_refs(
        precursor_idx,
        row_ids,
        scores,
        peps,
        frag_cols,
    )
    override = Pioneer._empirical_fragment_sqrt_tuple((
        0.0f0, 100.0f0, 0.0f0, 0.0f0,
        0.0f0, 0.0f0, 0.0f0, 0.0f0,
    ))

    stats = Pioneer._empirical_fragment_apply_reference_overrides!(
        refs,
        Dict(UInt64(1) => override),
    )
    hellinger, ref_pep = Pioneer._empirical_fragment_reference_features(
        precursor_idx,
        row_ids,
        frag_cols,
        refs,
    )

    @test stats.reference_rows == 2
    @test stats.replaced_rows == 1
    @test stats.fallback_rows == 1
    @test isapprox(hellinger[2], 0.0f0; atol = 1.0f-6)
    @test ref_pep[2] ≈ 0.01f0
end

@testset "empirical fragment features annotate fold files" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, DataFrame(
            precursor_idx = UInt32[20, 21],
            lgbm_prob = Float32[0.9, 0.7],
            main_pep = Float32[0.01, 0.02],
            frag1_int = Float32[100, 25],
            frag2_int = Float32[0, 25],
            frag3_int = Float32[0, 0],
            frag4_int = Float32[0, 0],
            frag5_int = Float32[0, 0],
            frag6_int = Float32[0, 0],
            frag7_int = Float32[0, 0],
            frag8_int = Float32[0, 0],
            shadow_frag1_int = Float32[0, 25],
            shadow_frag2_int = Float32[100, 25],
            shadow_frag3_int = Float32[0, 0],
            shadow_frag4_int = Float32[0, 0],
            shadow_frag5_int = Float32[0, 0],
            shadow_frag6_int = Float32[0, 0],
            shadow_frag7_int = Float32[0, 0],
            shadow_frag8_int = Float32[0, 0],
        ))
        Arrow.write(f2, DataFrame(
            precursor_idx = UInt32[20],
            lgbm_prob = Float32[0.8],
            main_pep = Float32[0.03],
            frag1_int = Float32[0],
            frag2_int = Float32[100],
            frag3_int = Float32[0],
            frag4_int = Float32[0],
            frag5_int = Float32[0],
            frag6_int = Float32[0],
            frag7_int = Float32[0],
            frag8_int = Float32[0],
            shadow_frag1_int = Float32[0],
            shadow_frag2_int = Float32[100],
            shadow_frag3_int = Float32[0],
            shadow_frag4_int = Float32[0],
            shadow_frag5_int = Float32[0],
            shadow_frag6_int = Float32[0],
            shadow_frag7_int = Float32[0],
            shadow_frag8_int = Float32[0],
        ))

        Pioneer.add_empirical_fragment_features_to_fold_files!([f1, f2])

        d1 = DataFrame(Arrow.Table(f1))
        d2 = DataFrame(Arrow.Table(f2))

        @test isapprox(d1.empirical_frag_best_hellinger[1], 0.0f0; atol = 1.0f-6)
        @test d1.empirical_frag_ref_pep[1] ≈ 0.03f0
        @test d1.empirical_frag_best_hellinger[2] ≈ 1.0f0
        @test d1.empirical_frag_ref_pep[2] ≈ 1.0f0
        @test isapprox(d2.empirical_frag_best_hellinger[1], 0.0f0; atol = 1.0f-6)
        @test d2.empirical_frag_ref_pep[1] ≈ 0.01f0
        @test !hasproperty(d1, :empirical_frag_corr_n)
        @test !hasproperty(d2, :empirical_frag_corr_n)
    end
end
