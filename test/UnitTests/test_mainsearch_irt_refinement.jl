using Test
using DataFrames
using Pioneer

using Pioneer: MainSearchIrtCorrectionModel, MainSearchIrtRefinement, _compute_phase2_columns!
using Pioneer: _passing_precursor_targets, refine_mainsearch_irt_predictions!
using Pioneer: reapply_psm_classifier_and_select_best!, train_lgbm_and_select_best

struct MainSearchMockPrecursors <: Pioneer.LibraryPrecursors
    sequence::Vector{String}
    structural_mods::Vector{Union{Missing, String}}
    mz::Vector{Float32}
    irt::Vector{Float32}
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

    @testset "per-run classifier can reapply after iRT feature updates" begin
        n = 40
        psms = DataFrame(
            precursor_idx = UInt32.(repeat(1:20, inner = 2)),
            scan_idx = UInt32.(repeat([10, 20], 20)),
            rt = Float32.(repeat([1, 2], 20)),
            cv_fold = UInt8[isodd(i) ? 0 : 1 for i in 1:n],
            target = Bool[isodd(cld(i, 2)) for i in 1:n],
            x = Float32[isodd(cld(i, 2)) ? 5 + (i % 2) : -5 - (i % 2) for i in 1:n],
        )

        _, _, _, predictor = train_lgbm_and_select_best(
            psms;
            features = [:x],
        )
        first_pass_scores = copy(psms.lgbm_score)
        psms.x .= .-psms.x

        best_psms, scores, _ = reapply_psm_classifier_and_select_best!(
            psms,
            predictor,
        )

        @test nrow(best_psms) == 20
        @test length(scores) == 20
        @test psms.lgbm_score != first_pass_scores
        @test all(isfinite, scores)
    end

    @testset "MBR nearest-neighbor pairing uses refined iRT column" begin
        precursors = MainSearchMockPrecursors(
            ["AAAA", "BBBB", "CCCC", "DDDD"],
            fill(missing, 4),
            fill(500.0f0, 4),
            Float32[0.0, 50.0, 1.0, 100.0],
        )
        psms = DataFrame(
            precursor_idx = UInt32[1, 2, 3, 4],
            target = Bool[true, true, false, false],
            cv_fold = UInt8[0, 0, 0, 0],
            irt_pred = Float32[0.0, 50.0, 100.0, 1.0],
        )

        Pioneer.regenerate_counterfactual_partners_irt_nn!(psms, precursors)

        @test psms.counterfactual_partner_pid[1] == UInt32(4)
    end
end
