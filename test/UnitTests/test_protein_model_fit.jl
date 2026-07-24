# Tests for src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/model_fit.jl.
#
# Targets the pure-numerics core that everything else in protein scoring
# depends on:
#   - fit_probit_model(X, y)         — argument validation + linear-separable recovery
#   - calculate_probit_scores(X, β)  — Φ(β·x) per row, in [0,1]
#   - end-to-end: fit then predict on linearly separable synthetic data,
#     training accuracy ≥ 0.95
#
# Skipped here: fit_probit_model_semisupervised (multi-iteration logic,
# intricate fixture), assign_protein_group_cv_folds!, apply_probit_scores_multifold!
# (DataFrame plumbing).

using DataFrames
using Random
using LinearAlgebra
using Distributions: Normal, cdf

using Pioneer: build_protein_semisupervised_training_set
using Pioneer: fit_probit_model, calculate_probit_scores

@testset "ProteinScoring/model_fit pure-numerics" begin

    @testset "fit_probit_model rejects empty feature matrix" begin
        X = Matrix{Float64}(undef, 5, 0)
        y = [true, false, true, false, true]
        @test_throws ArgumentError fit_probit_model(X, y)
    end

    @testset "fit_probit_model recovers sign on linearly separable data" begin
        # Single feature: positives have feature > 0, negatives < 0, well-separated.
        rng = MersenneTwister(13)
        n = 500
        X = zeros(Float64, n, 1)
        y = falses(n)
        @inbounds for i in 1:n
            is_pos = iseven(i)
            y[i] = is_pos
            X[i, 1] = is_pos ? 2.0 + 0.3 * randn(rng) : -2.0 + 0.3 * randn(rng)
        end
        β = fit_probit_model(X, y)
        # β is [intercept; feature_1]
        @test length(β) == 2
        @test β[2] > 0   # positive class has higher feature value → positive coeff
        @test all(isfinite, β)
    end

    @testset "calculate_probit_scores returns probabilities in [0,1]" begin
        X = randn(MersenneTwister(0), 100, 3)
        β = [0.1, 0.5, -0.3, 0.2]   # 1 intercept + 3 feature coeffs
        probs = calculate_probit_scores(X, β)
        @test length(probs) == 100
        @test all(0.0 .≤ probs .≤ 1.0)
    end

    @testset "calculate_probit_scores matches Φ(β·x) on hand-picked rows" begin
        # Construct a small example, compute Φ(β·x) directly, compare.
        β = [0.5, 1.0, -2.0]   # intercept + 2 features
        X = [
            0.0  0.0;
            1.0  1.0;
           -1.0 -1.0;
            2.0  0.0;
        ]
        probs = calculate_probit_scores(X, β)

        d = Normal(0, 1)
        for i in 1:size(X, 1)
            xb = β[1] + β[2] * X[i, 1] + β[3] * X[i, 2]
            @test probs[i] ≈ cdf(d, xb) atol=1f-6
        end
    end

    @testset "fit + predict round-trip: training accuracy ≥ 0.95 on separable data" begin
        rng = MersenneTwister(2026)
        n = 600
        # 2 features: positives are clustered around (1, 1), negatives around (-1, -1)
        X = zeros(Float64, n, 2)
        y = falses(n)
        for i in 1:n
            is_pos = rand(rng) > 0.5
            y[i] = is_pos
            μ = is_pos ? 1.0 : -1.0
            X[i, 1] = μ + 0.5 * randn(rng)
            X[i, 2] = μ + 0.5 * randn(rng)
        end

        β = fit_probit_model(X, y)
        probs = calculate_probit_scores(X, β)

        # Threshold at 0.5
        preds = probs .> 0.5
        accuracy = sum(preds .== y) / n
        @test accuracy ≥ 0.95
        # Both feature coefficients should be positive (positives have larger features)
        @test β[2] > 0
        @test β[3] > 0
    end

    @testset "calculate_probit_scores: extreme β values saturate near 0/1" begin
        # All-positive features × large positive β → probs near 1
        X = ones(Float64, 5, 1) .* 5.0
        β_pos = [0.0, 10.0]
        @test all(calculate_probit_scores(X, β_pos) .> 0.99)

        # Same features × large negative β → probs near 0
        β_neg = [0.0, -10.0]
        @test all(calculate_probit_scores(X, β_neg) .< 0.01)
    end

    @testset "ambiguous support protects both negative-mining paths" begin
        scores = Float32[0.9, 0.8, 0.2, 0.1, 0.05]
        targets = Bool[true, false, true, true, false]
        ambiguous_scores = Float32[0.0, 0.0, 2.0, 0.0, 0.0]
        neutral_shape = zeros(Float32, 5)
        nonsingletons = fill(Int64(2), 5)

        all_scores_weak = build_protein_semisupervised_training_set(
            scores,
            targets,
            neutral_shape,
            nonsingletons,
            ambiguous_scores;
            mined_negative_pep_threshold = 0.0f0,
        )
        @test all_scores_weak.mined_negative_mask ==
            Bool[true, false, false, true, false]
        @test all_scores_weak.positive_mask[3]

        low_shape_singletons = build_protein_semisupervised_training_set(
            scores,
            targets,
            Float32[-0.3, 0.0, -0.3, 0.0, 0.0],
            ones(Int64, 5),
            ambiguous_scores;
            mined_negative_pep_threshold = 1.1f0,
        )
        @test low_shape_singletons.mined_negative_mask ==
            Bool[true, false, false, false, false]
    end

end
