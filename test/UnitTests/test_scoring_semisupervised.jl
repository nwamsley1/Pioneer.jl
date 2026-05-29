using Test
using Pioneer

using Pioneer: _scoring_better_iteration_state, _scoring_semisupervised_train_mask
using Pioneer: _scoring_target_gain_sufficient

@testset "ScoringSearch semi-supervised helpers" begin
    @testset "post-first-iteration training keeps all decoys and q-value passing targets" begin
        mask = _scoring_semisupervised_train_mask(
            Bool[true, false, true, false, true],
            Float32[0.005, 0.8, 0.02, 0.5, 0.01],
        )

        @test mask == Bool[true, true, false, true, true]
    end

    @testset "patience-zero stopping requires at least one percent target gain" begin
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
end
