using Test
using DataFrames
using Pioneer

using Pioneer: MainSearchIrtCorrectionModel, MainSearchIrtRefinement, _compute_phase2_columns!
using Pioneer: _passing_precursor_targets, refine_mainsearch_irt_predictions!
using Pioneer: reapply_psm_classifier_and_select_best!, train_lgbm_and_select_best
using Pioneer: train_psm_classifier_with_fallback

struct MainSearchMockPrecursors <: Pioneer.LibraryPrecursors
    sequence::Vector{String}
    structural_mods::Vector{Union{Missing, String}}
    mz::Vector{Float32}
    irt::Vector{Float32}
end

@testset "MainSearch classifier fallback" begin
    @testset "low-data model selection accepts non-LightGBM candidates" begin
        psms = DataFrame(
            target = Bool[
                true, false, false, true,
                true, false, false, true,
                true, false, false, true,
            ],
            cv_fold = UInt8[
                1, 0, 1, 0,
                1, 0, 1, 0,
                1, 0, 1, 0,
            ],
            discriminant = Float32[
                0.9, 0.1, 0.2, 0.8,
                0.85, 0.15, 0.25, 0.75,
                0.88, 0.12, 0.22, 0.78,
            ],
        )

        tiny_lgbm_hp = (
            num_iterations = 3,
            learning_rate = 0.2,
            max_depth = 2,
            num_leaves = 3,
            min_data_in_leaf = 1,
            feature_fraction = 1.0,
            bagging_fraction = 1.0,
            bagging_freq = 0,
            is_unbalance = false,
            max_bin = 16,
            lambda_l1 = 0.0,
            lambda_l2 = 0.0,
        )

        scores, infold_scores, _last_classifier, info = train_psm_classifier_with_fallback(
            psms;
            features = [:discriminant],
            lgbm_hp = tiny_lgbm_hp,
        )

        @test length(scores) == nrow(psms)
        @test infold_scores === nothing
        @test info.low_data
        @test haskey(info.candidate_oof, "probit")
        @test all(isfinite, scores)
    end
end

Pioneer.getSequence(p::MainSearchMockPrecursors) = p.sequence
Pioneer.getStructuralMods(p::MainSearchMockPrecursors) = p.structural_mods
Pioneer.getMz(p::MainSearchMockPrecursors) = p.mz
Pioneer.getIrt(p::MainSearchMockPrecursors) = p.irt

@testset "MainSearch iRT refinement" begin
    @testset "passing targets use prescore q-value criteria" begin
        precursor_ids, irt_pred_inputs, irt_corrections = _passing_precursor_targets(
            UInt32[11, 22, 33],
            Bool[true, false, true],
            Float32[0.99, 0.98, 0.97],
            Float32[10.0, 20.0, 30.0],
            Float32[12.0, 18.0, 40.0],
            0.01f0,
        )

        @test precursor_ids == UInt32[11]
        @test irt_pred_inputs == Float32[10.0]
        @test irt_corrections == Float32[2.0]
    end

    @testset "refinement updates predictions and errors before classifier reapply" begin
        precursors = MainSearchMockPrecursors(
            ["AAAA", "CCCC", "AAAA", "CCCC", "DDDD", "DDDD"],
            fill(missing, 6),
            fill(500.0f0, 6),
            fill(10.0f0, 6),
        )
        psms = DataFrame(
            target = Bool[true, true, false, true, true, false],
            precursor_idx = UInt32[1, 2, 5, 3, 4, 6],
            irt_pred = fill(10.0f0, 6),
            irt_obs = Float32[4, 8, 6, 4, 8, 6],
            irt_error = fill(0.0f0, 6),
        )
        psms[!, :irt_error] = abs.(psms.irt_obs .- psms.irt_pred)
        scores = Float32[0.99, 0.98, 0.01, 0.99, 0.98, 0.01]
        best_psms = copy(psms)

        result = refine_mainsearch_irt_predictions!(
            psms,
            best_psms,
            scores,
            MainSearchIrtRefinement(precursors; q_value_threshold = 0.01f0, min_precursors = 2),
        )

        @test result.refined
        @test Set(result.training_target_precursors) == Set(UInt32[1, 2, 3, 4])
        @test psms.irt_pred[1] ≈ 4.0f0 atol = 1f-4
        @test psms.irt_pred[2] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_pred[4] ≈ 4.0f0 atol = 1f-4
        @test psms.irt_pred[5] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_error[1] ≈ 0.0f0 atol = 1f-4
        @test psms.irt_error[2] ≈ 0.0f0 atol = 1f-4
    end

    @testset "refinement trains out-of-fold models when cv_fold is present" begin
        precursors = MainSearchMockPrecursors(
            fill("AAAA", 6),
            fill(missing, 6),
            fill(500.0f0, 6),
            fill(10.0f0, 6),
        )
        psms = DataFrame(
            target = Bool[true, true, true, true, false, false],
            precursor_idx = UInt32[1, 2, 3, 4, 5, 6],
            cv_fold = UInt8[0, 0, 1, 1, 0, 1],
            irt_pred = fill(10.0f0, 6),
            irt_obs = Float32[8, 8, 13, 13, 8, 13],
            irt_error = fill(0.0f0, 6),
        )
        best_psms = copy(psms)
        scores = Float32[0.99, 0.98, 0.99, 0.98, 0.01, 0.01]

        result = refine_mainsearch_irt_predictions!(
            psms,
            best_psms,
            scores,
            MainSearchIrtRefinement(precursors; q_value_threshold = 0.01f0, min_precursors = 2),
        )

        @test result.refined
        @test result.model isa Dict{UInt8, MainSearchIrtCorrectionModel}
        @test sort(collect(keys(result.model))) == UInt8[0, 1]
        @test Set(result.training_target_precursors) == Set(UInt32[1, 2, 3, 4])
        @test psms.irt_pred[1] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[2] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[3] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_pred[4] ≈ 8.0f0 atol = 1f-4
        @test best_psms.irt_pred[1] ≈ 13.0f0 atol = 1f-4
        @test best_psms.irt_pred[3] ≈ 8.0f0 atol = 1f-4
    end

    @testset "phase two iRT difference uses refined prediction column" begin
        irt_diff_col = Vector{Float32}(undef, 2)
        prec_mz_col = Vector{Float32}(undef, 2)
        pair_id_col = Vector{UInt32}(undef, 2)
        entrap_col = Vector{UInt8}(undef, 2)

        _compute_phase2_columns!(
            UInt32[1, 2],
            Float32[12.0, 12.0],
            Float32[11.0, 15.0],
            Float32[100.0, 100.0],
            Float32[401.0, 402.0],
            UInt32[101, 102],
            UInt8[0, 1],
            irt_diff_col,
            prec_mz_col,
            pair_id_col,
            entrap_col,
        )

        @test irt_diff_col == Float32[1.0, 3.0]
        @test prec_mz_col == Float32[401.0, 402.0]
        @test pair_id_col == UInt32[101, 102]
        @test entrap_col == UInt8[0, 1]
    end

end
