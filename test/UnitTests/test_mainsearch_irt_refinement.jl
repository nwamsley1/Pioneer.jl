using Test
using DataFrames
using Pioneer

using Pioneer: MainSearchIrtRefinement, _compute_phase2_columns!
using Pioneer: _passing_precursor_targets, refine_mainsearch_irt_predictions!
using Pioneer: _select_irt_refinement_psms, select_best_per_precursor
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
            num_enzymatic_termini = fill(UInt8(2), 12),
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
            features = [:discriminant, :num_enzymatic_termini],
            lgbm_hp = tiny_lgbm_hp,
        )

        @test length(scores) == nrow(psms)
        @test infold_scores === nothing
        @test info.low_data
        @test haskey(info.candidate_oof, "probit")
        @test all(isfinite, scores)
        @test :num_enzymatic_termini ∉ info.available_features

        psms.num_enzymatic_termini .= repeat(UInt8[1, 2], 6)
        _, _, _, varying_info = train_psm_classifier_with_fallback(
            psms;
            features = [:discriminant, :num_enzymatic_termini],
            lgbm_hp = tiny_lgbm_hp,
        )
        @test :num_enzymatic_termini in varying_info.available_features
    end

    @testset "OOM feature selection drops only constant enzymatic termini" begin
        mktempdir() do temp_dir
            first_path = joinpath(temp_dir, "first.arrow")
            second_path = joinpath(temp_dir, "second.arrow")
            Arrow.write(first_path, (
                precursor_idx = UInt32[1, 2],
                discriminant = Float32[0.1, 0.2],
                num_enzymatic_termini = UInt8[2, 2],
            ))
            Arrow.write(second_path, (
                precursor_idx = UInt32[3, 4],
                discriminant = Float32[0.3, 0.4],
                num_enzymatic_termini = UInt8[2, 2],
            ))
            requested = [:discriminant, :num_enzymatic_termini]
            @test Pioneer._resolve_available_features(
                [first_path, second_path], requested
            ) == [:discriminant]

            Arrow.write(second_path, (
                precursor_idx = UInt32[3, 4],
                discriminant = Float32[0.3, 0.4],
                num_enzymatic_termini = UInt8[1, 2],
            ))
            @test Pioneer._resolve_available_features(
                [first_path, second_path], requested
            ) == requested
        end
    end
end

Pioneer.getSequence(p::MainSearchMockPrecursors) = p.sequence
Pioneer.getStructuralMods(p::MainSearchMockPrecursors) = p.structural_mods
Pioneer.getMz(p::MainSearchMockPrecursors) = p.mz
Pioneer.getIrt(p::MainSearchMockPrecursors) = p.irt

@testset "MainSearch iRT refinement" begin
    @testset "lightweight selection preserves exact winners" begin
        sorted_psms = DataFrame(
            precursor_idx = UInt32[1, 1, 1, 2, 2, 3, 3],
            target = Bool[true, true, true, false, false, true, true],
            cv_fold = UInt8[0, 0, 0, 1, 1, 0, 0],
            lgbm_score = Float32[0.9, 0.8, 0.7, 0.5, 0.5, 0.6, 0.6],
            weight = Float32[1, 10, 100, 4, 8, 7, 7],
            irt_pred = Float32[10, 10, 10, 20, 20, 30, 30],
            irt_obs = Float32[101, 102, 103, 201, 202, 301, 302],
            rt = Float32[1, 2, 3, 4, 5, 6, 7],
            scan_idx = UInt32[1, 2, 3, 4, 5, 6, 7],
        )
        center_mzs = Union{Missing, Float32}[500, 500, 500, 500, 500, 500, 500]
        isolation_widths = Union{Missing, Float32}[4, 4, 4, 4, 4, 4, 4]

        full = select_best_per_precursor(
            sorted_psms;
            center_mzs = center_mzs,
            isolation_widths = isolation_widths,
        )
        lightweight = _select_irt_refinement_psms(sorted_psms)

        @test lightweight.precursor_idx == full.precursor_idx
        @test lightweight.irt_obs == full.irt_obs
        @test lightweight.lgbm_score == full.lgbm_score
        @test lightweight.irt_obs == Float32[102, 202, 301]
        @test propertynames(lightweight) == [
            :precursor_idx,
            :target,
            :irt_pred,
            :irt_obs,
            :cv_fold,
            :lgbm_score,
        ]
        @test !hasproperty(lightweight, :irt_fwhm)
        @test !hasproperty(lightweight, :smoothness)
    end

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

    @testset "refinement trains out-of-fold models and leaves its fit table unchanged" begin
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
        refinement_psms = copy(psms)
        refinement_psms[!, :lgbm_score] = Float32[0.99, 0.98, 0.99, 0.98, 0.01, 0.01]
        refinement_before = copy(refinement_psms)

        result = refine_mainsearch_irt_predictions!(
            psms,
            refinement_psms,
            MainSearchIrtRefinement(precursors; q_value_threshold = 0.01f0, min_precursors = 2),
        )

        @test result.refined
        @test Set(result.training_target_precursors) == Set(UInt32[1, 2, 3, 4])
        @test psms.irt_pred[1] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[2] ≈ 13.0f0 atol = 1f-4
        @test psms.irt_pred[3] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_pred[4] ≈ 8.0f0 atol = 1f-4
        @test psms.irt_error ≈ fill(5.0f0, 6) atol = 1f-4
        @test refinement_psms == refinement_before
    end

    @testset "insufficient fold training data leaves predictions unchanged" begin
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
            lgbm_score = Float32[0.99, 0.98, 0.99, 0.98, 0.01, 0.01],
            weight = ones(Float32, 6),
            scan_idx = UInt32[1, 2, 3, 4, 5, 6],
            rt = Float32[1, 2, 3, 4, 5, 6],
            irt_pred = fill(10.0f0, 6),
            irt_obs = Float32[8, 8, 13, 13, 8, 13],
            irt_error = Float32[2, 2, 3, 3, 2, 3],
        )
        refinement_psms = _select_irt_refinement_psms(psms)
        psms_before = copy(psms)

        result = refine_mainsearch_irt_predictions!(
            psms,
            refinement_psms,
            MainSearchIrtRefinement(
                precursors;
                q_value_threshold = 0.01f0,
                min_precursors = 3,
            ),
        )

        @test !result.refined
        @test psms == psms_before
        best_psms = select_best_per_precursor(
            psms;
            center_mzs = fill(500.0f0, 6),
            isolation_widths = fill(4.0f0, 6),
        )
        @test best_psms.lgbm_score == psms.lgbm_score
        @test best_psms.irt_fwhm == zeros(Float32, 6)
        @test all(
            column -> hasproperty(best_psms, column),
            (:irt_fwhm, :n_above_hm, :rt_fwhm, :best_rt, :smoothness),
        )
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
