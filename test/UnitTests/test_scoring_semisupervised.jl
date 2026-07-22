using Test
using Pioneer

using Pioneer: _scoring_better_iteration_state, _scoring_semisupervised_train_mask
using Pioneer: _scoring_semisupervised_metrics_and_mask
using Pioneer: _scoring_target_gain_sufficient
using Pioneer: _select_global_scores
using Pioneer: get_qvalues!

@testset "ScoringSearch semi-supervised helpers" begin
    @testset "post-first-iteration training keeps all decoys and q-value passing targets" begin
        mask = _scoring_semisupervised_train_mask(
            Bool[true, false, true, false, true, true, true],
            Float32[0.005, 0.8, 0.02, 0.5, 0.01, 0.03, 0.031],
        )

        @test mask == Bool[true, true, true, true, true, true, false]
    end

    @testset "stopping requires configured target gain" begin
        @test _scoring_target_gain_sufficient(100, 101)
        @test !_scoring_target_gain_sufficient(100, 100)
        @test !_scoring_target_gain_sufficient(100, 100; min_fraction = 0.01f0)
        @test _scoring_target_gain_sufficient(0, 1)
        @test !_scoring_target_gain_sufficient(0, 0)
    end

    @testset "stopping keeps the iteration with more q-value passing targets" begin
        previous = (target_q01 = 7025, iter = 2, token = :previous)
        current = (target_q01 = 6935, iter = 3, token = :current)

        chosen = _scoring_better_iteration_state(previous, current)

        @test chosen.iter == 2
        @test chosen.token == :previous
        @test _scoring_better_iteration_state(current, previous).iter == 2
    end

    @testset "global score selection compares LightGBM with the empirical score" begin
        targets = Bool[true, false, false, true, true]
        better_scores = Float32[0.9, 0.2, 0.1, 0.8, 0.7]
        worse_scores = Float32[0.9, 0.8, 0.7, 0.6, 0.5]

        empirical = _select_global_scores(
            worse_scores,
            better_scores,
            targets;
            scoring_name = "Test",
            q_threshold = 0.5f0,
        )
        @test empirical.source == :empirical
        @test empirical.scores == better_scores
        @test empirical.empirical_target_count == 3
        @test empirical.model_target_count == 1

        model = _select_global_scores(
            better_scores,
            worse_scores,
            targets;
            scoring_name = "Test",
            q_threshold = 0.5f0,
        )
        @test model.source == :lightgbm
        @test model.scores == better_scores

        tied = _select_global_scores(
            better_scores,
            copy(better_scores),
            targets;
            scoring_name = "Test",
            q_threshold = 0.5f0,
        )
        @test tied.source == :empirical

        scaled = _select_global_scores(
            worse_scores,
            better_scores,
            targets;
            scoring_name = "Test",
            q_threshold = 0.5f0,
            fdr_scale_factor = 0.5f0,
        )
        @test scaled.model_target_count == 3
        @test scaled.source == :empirical
    end

    @testset "fused metrics use separate training and stopping q-value thresholds" begin
        scores = Float32[0.91, 0.88, 0.70, 0.60, 0.55, 0.20]
        targets = Bool[true, false, true, true, false, false]
        qvals = Vector{Float32}(undef, length(scores))
        get_qvalues!(scores, targets, qvals; doSort = true, fdr_scale_factor = 1.0f0)

        metrics = _scoring_semisupervised_metrics_and_mask(
            scores,
            targets;
            train_q_threshold = 0.75f0,
            stop_q_threshold = 0.50f0,
        )
        train_q_pass = qvals .<= 0.75f0
        stop_q_pass = qvals .<= 0.50f0

        @test metrics.training_mask == BitVector((.!targets) .| train_q_pass)
        @test metrics.target_q01 == count(stop_q_pass .& targets)
        @test metrics.decoy_q01 == count(stop_q_pass .& .!targets)
    end
end
