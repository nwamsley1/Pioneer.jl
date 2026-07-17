# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

@testset "Global protein scoring" begin
    @testset "score-distribution feature table" begin
        key1 = ("P1", true, UInt8(0))
        key2 = ("P2", false, UInt8(1))
        inputs = Pioneer.GlobalProteinInputs(
            Dict(
                key1 => Float32[0.9, 0.8],
                key2 => Float32[0.7],
            ),
            Dict(key1 => UInt8(0), key2 => UInt8(1)),
        )
        table = Pioneer._build_global_protein_feature_table(inputs, 4)

        @test table.protein_name == ["P1", "P2"]
        @test table.target == Bool[true, false]
        @test table.entrap_id == UInt8[0, 1]
        @test table.cv_fold == UInt8[0, 1]
        @test Symbol.(names(table)[5:end]) == Pioneer.GLOBAL_PROTEIN_SCORE_FEATURES

        @test table.empirical_global_score[1] ≈ logodds(Float32[0.9, 0.8], 2)
        @test table.top1_pg_score[1] == 0.9f0
        @test table.top2_pg_score[1] == 0.8f0
        @test table.top3_pg_score[1] == 0.0f0
        @test table.top2_logodds_score[1] ≈ logodds(Float32[0.9, 0.8], 2)
        @test table.top3_logodds_score[1] ≈ logodds(Float32[0.9, 0.8], 3)
        @test table.mean_pg_score[1] ≈ 0.85f0
        @test table.median_pg_score[1] ≈ 0.85f0
        @test table.std_pg_score[1] ≈ sqrt(0.005f0)
        @test table.min_pg_score[1] == 0.8f0
        @test table.top1_top2_gap[1] ≈ 0.1f0
        @test table.top2_top3_gap[1] == 0.0f0
        @test table.n_runs_observed[1] == 2.0f0
        @test table.observed_run_fraction[1] == 0.5f0
        @test table.n_score_gt_0_5[1] == 2.0f0
        @test table.n_score_gt_0_9[1] == 0.0f0
        @test table.n_score_gt_0_99[1] == 0.0f0

        @test table.std_pg_score[2] == 0.0f0
        @test table.top1_top2_gap[2] == 0.0f0
        @test table.observed_run_fraction[2] == 0.25f0
    end

    @testset "run-level input collection" begin
        mktempdir() do directory
            run1_path = joinpath(directory, "run1.arrow")
            run2_path = joinpath(directory, "run2.arrow")
            Arrow.write(run1_path, DataFrame(
                protein_name = ["P1", "P2"],
                target = Bool[true, false],
                entrap_id = UInt8[0, 0],
                pg_score = Float32[0.9, 0.4],
            ))
            Arrow.write(run2_path, DataFrame(
                protein_name = ["P1", "P3"],
                target = Bool[true, false],
                entrap_id = UInt8[0, 1],
                pg_score = Float32[0.8, 0.3],
            ))
            refs = Pioneer.ProteinGroupFileReference[
                Pioneer.ProteinGroupFileReference(run1_path),
                Pioneer.ProteinGroupFileReference(run2_path),
            ]
            protein_to_cv_fold = Dictionary{
                String,
                @NamedTuple{best_score::Float32, cv_fold::UInt8},
            }()
            insert!(protein_to_cv_fold, "P1", (best_score = 0.9f0, cv_fold = UInt8(0)))
            insert!(protein_to_cv_fold, "P2", (best_score = 0.4f0, cv_fold = UInt8(1)))
            insert!(protein_to_cv_fold, "P3", (best_score = 0.3f0, cv_fold = UInt8(1)))

            inputs = Pioneer._collect_global_protein_inputs(
                refs,
                protein_to_cv_fold,
                3,
            )

            @test inputs.scores[("P1", true, UInt8(0))] == Float32[0.9, 0.8]
            @test inputs.scores[("P2", false, UInt8(0))] == Float32[0.4]
            @test inputs.scores[("P3", false, UInt8(1))] == Float32[0.3]
            @test inputs.folds == Dict(
                ("P1", true, UInt8(0)) => UInt8(0),
                ("P2", false, UInt8(0)) => UInt8(1),
                ("P3", false, UInt8(1)) => UInt8(1),
            )
        end
    end

    @testset "one fold uses empirical global protein scores" begin
        mktempdir() do directory
            run_scores = (
                Float32[0.9, 0.6],
                Float32[0.8, 0.5],
                Float32[0.7, 0.4],
                Float32[0.6, 0.3],
            )
            refs = Pioneer.ProteinGroupFileReference[]
            for (run_idx, scores) in enumerate(run_scores)
                path = joinpath(directory, "run$(run_idx).arrow")
                Arrow.write(path, DataFrame(
                    protein_name = ["P1", "P2"],
                    target = Bool[true, false],
                    entrap_id = UInt8[0, 0],
                    pg_score = scores,
                ))
                push!(refs, Pioneer.ProteinGroupFileReference(path))
            end
            protein_to_cv_fold = Dictionary{
                String,
                @NamedTuple{best_score::Float32, cv_fold::UInt8},
            }()
            insert!(protein_to_cv_fold, "P1", (best_score = 0.9f0, cv_fold = UInt8(0)))
            insert!(protein_to_cv_fold, "P2", (best_score = 0.6f0, cv_fold = UInt8(0)))

            scores, protein_scores = Pioneer.build_global_protein_score_dicts(
                refs,
                protein_to_cv_fold,
                2,
                4,
            )

            expected_p1 = logodds(Float32[0.9, 0.8, 0.7, 0.6], 2)
            expected_p2 = logodds(Float32[0.6, 0.5, 0.4, 0.3], 2)
            @test scores == Dict(
                ("P1", true, UInt8(0)) => expected_p1,
                ("P2", false, UInt8(0)) => expected_p2,
            )
            @test protein_scores == Dict(
                ProteinKey("P1", true, UInt8(0)) => expected_p1,
                ProteinKey("P2", false, UInt8(0)) => expected_p2,
            )
        end
    end

    @testset "protein features support out-of-fold LightGBM scoring" begin
        scores = Dict{Tuple{String, Bool, UInt8}, Vector{Float32}}()
        folds = Dict{Tuple{String, Bool, UInt8}, UInt8}()
        for group_idx in 1:40
            fold = UInt8(group_idx % 2)
            target_key = ("target_$(group_idx)", true, UInt8(0))
            decoy_key = ("decoy_$(group_idx)", false, UInt8(0))
            scores[target_key] = Float32[0.95 - 0.001 * group_idx, 0.9]
            scores[decoy_key] = Float32[0.05 + 0.001 * group_idx, 0.1]
            folds[target_key] = fold
            folds[decoy_key] = fold
        end
        table = Pioneer._build_global_protein_feature_table(
            Pioneer.GlobalProteinInputs(scores, folds),
            2,
        )
        test_hyperparameters = (
            num_iterations = 20,
            learning_rate = 0.2,
            max_depth = 2,
            num_leaves = 4,
            min_data_in_leaf = 1,
            min_gain_to_split = 0.0,
            feature_fraction = 1.0,
            bagging_fraction = 1.0,
            bagging_freq = 0,
            is_unbalance = false,
            max_bin = 31,
            lambda_l1 = 0.0,
            lambda_l2 = 0.0,
        )

        scored = Pioneer._score_global_features_oof(
            table,
            Pioneer.GLOBAL_PROTEIN_SCORE_FEATURES,
            UInt8[0, 1];
            scoring_name = "Global protein",
            lgbm_hp = test_hyperparameters,
            min_training_class_count = 1,
            max_train = 1_000,
            max_iterations = 3,
        )

        @test length(scored.scores) == nrow(table)
        @test length(scored.models) == 2
        @test scored.iter in 1:3
        @test mean(scored.scores[table.target]) > mean(scored.scores[.!table.target])
    end
end
