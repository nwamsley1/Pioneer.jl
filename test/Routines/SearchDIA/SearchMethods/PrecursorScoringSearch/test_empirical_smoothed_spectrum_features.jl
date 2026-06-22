import Pioneer
using Arrow
using DataFrames
using Test

function _empirical_smoothed_test_table(; precursor_idx, lgbm_prob, frag1, frag2)
    n = length(precursor_idx)
    return DataFrame(
        precursor_idx = UInt32.(precursor_idx),
        lgbm_prob = Float32.(lgbm_prob),
        frag1_smoothed_intensity = Float32.(frag1),
        frag2_smoothed_intensity = Float32.(frag2),
        frag3_smoothed_intensity = zeros(Float32, n),
        frag4_smoothed_intensity = zeros(Float32, n),
        frag5_smoothed_intensity = zeros(Float32, n),
        frag6_smoothed_intensity = zeros(Float32, n),
        frag7_smoothed_intensity = zeros(Float32, n),
        frag8_smoothed_intensity = zeros(Float32, n),
    )
end

@testset "empirical smoothed spectrum references use top-2 leave-one-out" begin
    precursor_idx = UInt32[10, 10, 10, 11]
    row_ids = UInt64[1, 2, 3, 4]
    scores = Float32[0.9, 0.8, 0.7, 0.6]
    frag_cols = (
        Float32[100, 0, 100, 50],
        Float32[0, 100, 0, 50],
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
        zeros(Float32, 4),
    )

    refs = Pioneer._empirical_smoothed_spectrum_topk_refs(
        precursor_idx,
        row_ids,
        scores,
        frag_cols,
    )
    hellinger = Pioneer._empirical_smoothed_spectrum_reference_features(
        precursor_idx,
        row_ids,
        frag_cols,
        refs,
    )

    @test hellinger[1] ≈ 1.0f0
    @test hellinger[2] ≈ 1.0f0
    @test hellinger[3] ≈ 0.0f0
    @test hellinger[4] ≈ 1.0f0
end

@testset "empirical smoothed spectrum feature annotates fold files" begin
    mktempdir() do dir
        f1 = joinpath(dir, "run1_fold0.arrow")
        f2 = joinpath(dir, "run2_fold0.arrow")

        Arrow.write(f1, _empirical_smoothed_test_table(
            precursor_idx = [20, 21],
            lgbm_prob = [0.9, 0.7],
            frag1 = [100, 25],
            frag2 = [0, 25],
        ))
        Arrow.write(f2, _empirical_smoothed_test_table(
            precursor_idx = [20],
            lgbm_prob = [0.8],
            frag1 = [100],
            frag2 = [0],
        ))

        Pioneer.add_empirical_smoothed_spectrum_features_to_fold_files!([f1, f2])

        d1 = DataFrame(Arrow.Table(f1))
        d2 = DataFrame(Arrow.Table(f2))

        @test isapprox(d1.empirical_smoothed_frag_hellinger[1], 0.0f0; atol = 1.0f-6)
        @test d1.empirical_smoothed_frag_hellinger[2] ≈ 1.0f0
        @test isapprox(d2.empirical_smoothed_frag_hellinger[1], 0.0f0; atol = 1.0f-6)
        @test :empirical_smoothed_frag_hellinger in Pioneer.ADVANCED_FEATURE_SET
        @test !hasproperty(d1, :empirical_frag_ref_pep)
        @test !hasproperty(d2, :empirical_frag_ref_pep)
    end
end
