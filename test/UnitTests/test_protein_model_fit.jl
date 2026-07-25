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
using Statistics: mean

using Pioneer: build_protein_semisupervised_training_set
using Pioneer: build_lightgbm_classifier
using Pioneer: fit_probit_model, calculate_probit_scores
using Pioneer: fit_protein_lightgbm_semisupervised, lightgbm_predict
using Pioneer: log_probit_feature_importance
using Pioneer: RUN_LEVEL_PROTEIN_LGBM_HP
using Pioneer: _protein_lightgbm_monotone_constraints
import Pioneer: DEBUG_CONSOLE_LEVEL

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

    @testset "protein LightGBM learns a discriminating score" begin
        n_per_class = 500
        signal = vcat(
            fill(2.0f0, n_per_class),
            fill(-2.0f0, n_per_class),
        )
        feature_data = DataFrame(
            pg_score = signal,
            precursor_consensus_prefix_shape = signal,
            shared_precursor_consensus_prefix_shape = signal,
        )
        targets = vcat(
            trues(n_per_class),
            falses(n_per_class),
        )
        initial_scores = vcat(
            fill(2.0f0, n_per_class),
            fill(0.1f0, n_per_class),
        )
        model = fit_protein_lightgbm_semisupervised(
            feature_data,
            targets,
            initial_scores,
            zeros(Float32, 2 * n_per_class),
            fill(Int64(2), 2 * n_per_class);
            n_iterations = 1,
        )
        scores = lightgbm_predict(
            model,
            feature_data;
            output_type = Float32,
        )
        @test model.booster.monotone_constraints == Int[1, 1, 1]
        @test mean(scores[targets]) > mean(scores[.!targets])
    end

    @testset "protein LightGBM retains flexible parameters and monotonic directions" begin
        @test RUN_LEVEL_PROTEIN_LGBM_HP.num_iterations == 100
        @test RUN_LEVEL_PROTEIN_LGBM_HP.max_depth == 4
        @test RUN_LEVEL_PROTEIN_LGBM_HP.num_leaves == 15
        @test RUN_LEVEL_PROTEIN_LGBM_HP.min_data_in_leaf == 10
        @test RUN_LEVEL_PROTEIN_LGBM_HP.min_gain_to_split == 0.0
        @test RUN_LEVEL_PROTEIN_LGBM_HP.max_bin == 255
        @test RUN_LEVEL_PROTEIN_LGBM_HP.lambda_l1 == 1.0
        @test RUN_LEVEL_PROTEIN_LGBM_HP.lambda_l2 == 1.0
        @test RUN_LEVEL_PROTEIN_LGBM_HP.monotone_constraints_method ==
            "intermediate"

        feature_names = Symbol[
            :pg_score,
            :ambiguous_pg_score,
            :shared_peptide_coverage_logit,
            :shared_coverage_log_ratio,
            :peptide_coverage_logit,
            :any_common_peps,
            :coverage_log_ratio,
            :precursor_consensus_prefix_shape,
            :shared_precursor_consensus_prefix_shape,
            :mbr_recovered_peptides,
            :mbr_only_protein,
        ]
        constraints =
            _protein_lightgbm_monotone_constraints(feature_names)
        @test constraints == Int[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1]

        classifier = build_lightgbm_classifier(
            ;
            RUN_LEVEL_PROTEIN_LGBM_HP...,
            monotone_constraints = constraints,
        )
        @test classifier.monotone_constraints == constraints
        @test classifier.monotone_constraints_method == "intermediate"
    end

    @testset "protein LightGBM stops on target-ID gain" begin
        n_per_class = 500
        feature_data = DataFrame(
            signal = vcat(
                fill(2.0f0, n_per_class),
                fill(-2.0f0, n_per_class),
            ),
        )
        targets = vcat(trues(n_per_class), falses(n_per_class))
        initial_scores = vcat(
            fill(2.0f0, n_per_class),
            fill(0.1f0, n_per_class),
        )
        old_debug_level = DEBUG_CONSOLE_LEVEL[]

        try
            DEBUG_CONSOLE_LEVEL[] = 1
            output = mktemp() do _, io
                redirect_stdout(io) do
                    fit_protein_lightgbm_semisupervised(
                        feature_data,
                        targets,
                        initial_scores,
                        zeros(Float32, 2 * n_per_class),
                        fill(Int64(2), 2 * n_per_class);
                        n_iterations = 5,
                        context = "unit_test",
                    )
                end
                flush(io)
                seekstart(io)
                read(io, String)
            end

            @test occursin(
                "Protein LightGBM monotone constraints context=unit_test: " *
                "method=intermediate",
                output,
            )
            @test occursin(
                "Protein LightGBM semi-supervised iter 1 context=unit_test",
                output,
            )
            @test occursin(
                "current training: target positives=500, " *
                "mined target negatives=0, decoys=500, rows=1000",
                output,
            )
            @test occursin(
                "IDs at q≤0.01: targets=500, decoys=5",
                output,
            )
            @test occursin(
                "next training: target positives=500, " *
                "mined target negatives=0, decoys=500, dropped targets=0",
                output,
            )
            @test occursin(
                "Protein LightGBM semi-supervised continuing " *
                "context=unit_test after iter 1: " *
                "established target-ID baseline=500",
                output,
            )
            @test occursin(
                "Protein LightGBM semi-supervised stopping " *
                "context=unit_test after iter 2: " *
                "q≤0.01 target IDs=500 did not improve by 1.0% over 500; " *
                "using iter 2 with target IDs=500",
                output,
            )
        finally
            DEBUG_CONSOLE_LEVEL[] = old_debug_level
        end
    end

    @testset "protein probit feature importance uses level-1 debug logging" begin
        X = Float64[
            0.0  0.0;
            1.0  2.0;
            2.0  4.0;
        ]
        β = Float64[0.25, 0.5, -1.0]
        feature_names = Symbol[:feature_a, :feature_b]
        old_debug_level = DEBUG_CONSOLE_LEVEL[]

        try
            function capture_feature_importance_log(debug_level)
                DEBUG_CONSOLE_LEVEL[] = debug_level
                return mktemp() do _, io
                    redirect_stdout(io) do
                        log_probit_feature_importance(
                            feature_names,
                            β,
                            X;
                            context = "unit_test"
                        )
                    end
                    flush(io)
                    seekstart(io)
                    return read(io, String)
                end
            end

            @test isempty(capture_feature_importance_log(0))

            output = capture_feature_importance_log(1)
            @test occursin(
                "Protein probit model coefficients context=unit_test",
                output
            )
            @test occursin(
                "Protein probit feature importance context=unit_test rank=1 " *
                "feature=feature_b coefficient=-1.0",
                output
            )
            @test occursin(
                "Protein probit feature importance context=unit_test rank=2 " *
                "feature=feature_a coefficient=0.5",
                output
            )
        finally
            DEBUG_CONSOLE_LEVEL[] = old_debug_level
        end
    end

    @testset "negative mining uses non-MBR unique peptide count" begin
        scores = Float32[0.9, 0.8, 0.2, 0.1, 0.05, 0.3]
        targets = Bool[true, false, true, true, false, true]
        neutral_shape = zeros(Float32, 6)
        nonsingletons = fill(Int64(2), 6)

        all_scores_weak = build_protein_semisupervised_training_set(
            scores,
            targets,
            neutral_shape,
            nonsingletons;
            mined_negative_pep_threshold = 0.0f0,
        )
        @test all_scores_weak.mined_negative_mask ==
            Bool[true, false, true, true, false, true]

        low_shape_singletons = build_protein_semisupervised_training_set(
            scores,
            targets,
            Float32[-0.3, 0.0, -0.3, 0.0, 0.0, 0.0],
            ones(Int64, 6);
            mined_negative_pep_threshold = 1.1f0,
        )
        @test low_shape_singletons.mined_negative_mask ==
            Bool[true, false, true, false, false, false]

        low_shape_zero_non_mbr = build_protein_semisupervised_training_set(
            scores,
            targets,
            Float32[-0.3, 0.0, -0.3, 0.0, 0.0, 0.0],
            zeros(Int64, 6);
            mined_negative_pep_threshold = 1.1f0,
        )
        @test low_shape_zero_non_mbr.mined_negative_mask ==
            Bool[true, false, true, false, false, false]

        low_shape_non_mbr_nonsingletons =
            build_protein_semisupervised_training_set(
            scores,
            targets,
            Float32[-0.3, 0.0, -0.3, 0.0, 0.0, 0.0],
            nonsingletons;
            mined_negative_pep_threshold = 1.1f0,
        )
        @test !any(low_shape_non_mbr_nonsingletons.mined_negative_mask)
    end

end
