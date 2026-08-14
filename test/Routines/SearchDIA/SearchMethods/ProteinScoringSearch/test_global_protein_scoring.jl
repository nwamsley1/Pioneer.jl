# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

@testset "Global protein scoring" begin
    @testset "monotonic constraints follow feature names" begin
        feature_names = Symbol[
            :std_pg_score,
            :global_peptide_coverage,
            :empirical_global_score,
            :missing_run_similarity_mass_approx,
        ]
        constraints = Pioneer._global_lightgbm_monotone_constraints(
            feature_names,
            Pioneer.GLOBAL_PROTEIN_MONOTONE_INCREASING_FEATURES,
        )

        @test constraints == Int[0, 1, 1, 0]
        @test Pioneer.GLOBAL_PROTEIN_LGBM_HP.monotone_constraints_method ==
              "intermediate"
    end

    @testset "run-level protein model requires sufficient data in every fold" begin
        table = DataFrame(
            target = vcat(
                fill(true, 10),
                fill(false, 10),
                fill(true, 9),
                fill(false, 10),
            ),
            cv_fold = vcat(
                fill(UInt8(0), 20),
                fill(UInt8(1), 19),
            ),
        )

        @test !Pioneer._protein_model_training_data_sufficient(
            table,
            UInt8[0, 1],
        )
        push!(table, (target = true, cv_fold = UInt8(1)))
        @test Pioneer._protein_model_training_data_sufficient(
            table,
            UInt8[0, 1],
        )
    end

    @testset "unavailable model scoring uses initial probabilities for every fold" begin
        mktempdir() do directory
            path = joinpath(directory, "protein_groups.arrow")
            Arrow.write(path, DataFrame(
                protein_name = ["P0", "P1"],
                target = Bool[true, false],
                pg_score = Float32[0.5, 1.5],
                feature = Float32[1.0, -1.0],
            ))
            refs = Pioneer.ProteinGroupFileReference[
                Pioneer.ProteinGroupFileReference(path),
            ]
            protein_to_cv_fold = Dictionary{
                String,
                @NamedTuple{best_score::Float32, cv_fold::UInt8},
            }()
            insert!(protein_to_cv_fold, "P0", (best_score = 0.5f0, cv_fold = UInt8(0)))
            insert!(protein_to_cv_fold, "P1", (best_score = 1.5f0, cv_fold = UInt8(1)))

            Pioneer.apply_protein_scores_multifold!(
                refs,
                protein_to_cv_fold,
                Dict{UInt8, Pioneer.LightGBMModel}(),
                Symbol[:feature];
                use_model_scores = false,
            )

            table = DataFrame(Arrow.Table(path); copycols = true)
            sort!(table, :protein_name)
            @test table.old_pg_score == Float32[0.5, 1.5]
            @test table.pg_score ≈ Pioneer._initial_protein_probabilities(
                Float32[0.5, 1.5],
            )
        end
    end

    @testset "score-distribution feature table" begin
        key1 = ("P1", true, UInt8(0))
        key2 = ("P2", false, UInt8(1))
        inputs = Pioneer.GlobalProteinInputs(
            Dict(
                key1 => Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(1), 0.9f0),
                    Pioneer.GlobalProteinRunScore(UInt32(2), 0.8f0),
                ],
                key2 => Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(1), 0.7f0),
                ],
            ),
            Dict(
                key1 => Set(["PEPA", "PEPB", "PEPC"]),
                key2 => Set(["PEPX"]),
            ),
            Dict(key1 => 2, key2 => 1),
            Dict(key1 => 6, key2 => 4),
            Dict(key1 => UInt8(0), key2 => UInt8(1)),
        )
        run_similarity = Pioneer.build_run_similarity(Dict(
            UInt32(1) => UInt32[100, 101],
            UInt32(2) => UInt32[100, 101],
            UInt32(3) => UInt32[102],
            UInt32(4) => UInt32[],
        ))
        table = Pioneer._build_global_protein_feature_table(
            inputs,
            4,
            0.85f0,
            run_similarity = run_similarity,
        )

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
        @test table.n_unique_peptides_observed[1] == 3.0f0
        @test table.max_n_peptides_observed[1] == 2.0f0
        @test table.global_peptide_coverage[1] == 0.5f0
        @test table.max_peptide_coverage[1] ≈ 1.0f0 / 3.0f0
        @test table.global_all_peptide_coverage[1] == 0.5f0
        @test table.max_all_peptide_coverage[1] ≈ 1.0f0 / 3.0f0
        @test table.n_passing_runs[1] == 1.0f0
        @test !hasproperty(table, :n_runs_observed)
        @test table.n_score_gt_0_5[1] == 2.0f0
        @test table.n_score_gt_0_9[1] == 0.0f0
        @test table.n_score_gt_0_99[1] == 0.0f0
        @test table.observed_run_centrality_mean[1] ≈ 1.0f0 / 3.0f0
        @test table.observed_run_centrality_max[1] ≈ 1.0f0 / 3.0f0
        @test table.missing_run_similarity_mass_approx[1] ≈ 1.0f0

        @test table.std_pg_score[2] == 0.0f0
        @test table.top1_top2_gap[2] == 0.0f0
        @test table.n_unique_peptides_observed[2] == 1.0f0
        @test table.max_n_peptides_observed[2] == 1.0f0
        @test table.global_peptide_coverage[2] == 0.25f0
        @test table.max_peptide_coverage[2] == 0.25f0
        @test table.global_all_peptide_coverage[2] == 0.25f0
        @test table.max_all_peptide_coverage[2] == 0.25f0
        @test table.n_passing_runs[2] == 0.0f0
        @test table.observed_run_centrality_mean[2] == 0.0f0
        @test table.observed_run_centrality_max[2] == 0.0f0
        @test table.missing_run_similarity_mass_approx[2] == 0.0f0
    end

    @testset "global coverage separates common and all peptides" begin
        key = ("P1", true, UInt8(0))
        inputs = Pioneer.GlobalProteinInputs(
            Dict(
                key => Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(1), 0.9f0),
                    Pioneer.GlobalProteinRunScore(UInt32(2), 0.8f0),
                ],
            ),
            Dict(key => Set(["COMMON", "SEMI1", "SEMI2"])),
            Dict(key => Set(["COMMON"])),
            Dict(key => 2),
            Dict(key => 1),
            Dict(key => 6),
            Dict(key => 2),
            Dict(key => UInt8(0)),
        )
        table = Pioneer._build_global_protein_feature_table(inputs, 2, 0.5f0)

        @test table.global_peptide_coverage[1] == 0.5f0
        @test table.max_peptide_coverage[1] == 0.5f0
        @test table.global_all_peptide_coverage[1] == 0.5f0
        @test table.max_all_peptide_coverage[1] ≈ 1.0f0 / 3.0f0
        @test :global_all_peptide_coverage in
            Pioneer.GLOBAL_PROTEIN_MONOTONE_INCREASING_FEATURES
        @test :max_all_peptide_coverage in
            Pioneer.GLOBAL_PROTEIN_MONOTONE_INCREASING_FEATURES
    end

    @testset "run-level input collection" begin
        mktempdir() do directory
            run1_path = joinpath(directory, "run1.arrow")
            run2_path = joinpath(directory, "run2.arrow")
            Arrow.write(run1_path, DataFrame(
                protein_name = ["P1", "P2"],
                target = Bool[true, false],
                entrap_id = UInt8[0, 0],
                file_idx = UInt32[2, 2],
                pg_score = Float32[0.9, 0.4],
                n_peptides = Int[2, 1],
                n_common_peptides = Int[1, 0],
                n_possible_unique_peptides = Int[4, 2],
                n_possible_common_unique_peptides = Int[2, 1],
                peptide_list = ["PEPA;PEPB", "PEPX"],
                common_peptide_list = ["PEPA", ""],
            ))
            Arrow.write(run2_path, DataFrame(
                protein_name = ["P1", "P3"],
                target = Bool[true, false],
                entrap_id = UInt8[0, 1],
                file_idx = UInt32[4, 4],
                pg_score = Float32[0.8, 0.3],
                n_peptides = Int[2, 1],
                n_common_peptides = Int[1, 1],
                n_possible_unique_peptides = Int[4, 3],
                n_possible_common_unique_peptides = Int[2, 2],
                peptide_list = ["PEPB;PEPC", "PEPY"],
                common_peptide_list = ["PEPC", "PEPY"],
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

            @test inputs.run_scores[("P1", true, UInt8(0))] ==
                Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(2), 0.9f0),
                    Pioneer.GlobalProteinRunScore(UInt32(4), 0.8f0),
                ]
            @test inputs.run_scores[("P2", false, UInt8(0))] ==
                Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(2), 0.4f0),
                ]
            @test inputs.run_scores[("P3", false, UInt8(1))] ==
                Pioneer.GlobalProteinRunScore[
                    Pioneer.GlobalProteinRunScore(UInt32(4), 0.3f0),
                ]
            @test inputs.observed_peptides[("P1", true, UInt8(0))] ==
                  Set(["PEPA", "PEPB", "PEPC"])
            @test inputs.observed_peptides[("P2", false, UInt8(0))] == Set(["PEPX"])
            @test inputs.observed_peptides[("P3", false, UInt8(1))] == Set(["PEPY"])
            @test inputs.observed_common_peptides[("P1", true, UInt8(0))] ==
                Set(["PEPA", "PEPC"])
            @test isempty(inputs.observed_common_peptides[("P2", false, UInt8(0))])
            @test inputs.observed_common_peptides[("P3", false, UInt8(1))] ==
                Set(["PEPY"])
            @test inputs.max_n_peptides == Dict(
                ("P1", true, UInt8(0)) => 2,
                ("P2", false, UInt8(0)) => 1,
                ("P3", false, UInt8(1)) => 1,
            )
            @test inputs.n_possible_unique_peptides == Dict(
                ("P1", true, UInt8(0)) => 4,
                ("P2", false, UInt8(0)) => 2,
                ("P3", false, UInt8(1)) => 3,
            )
            @test inputs.max_n_common_peptides == Dict(
                ("P1", true, UInt8(0)) => 1,
                ("P2", false, UInt8(0)) => 0,
                ("P3", false, UInt8(1)) => 1,
            )
            @test inputs.n_possible_common_unique_peptides == Dict(
                ("P1", true, UInt8(0)) => 2,
                ("P2", false, UInt8(0)) => 1,
                ("P3", false, UInt8(1)) => 2,
            )
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
                    file_idx = fill(UInt32(run_idx), 2),
                    pg_score = scores,
                    n_peptides = Int[2, 1],
                    n_possible_unique_peptides = Int[4, 3],
                    peptide_list = ["PEPA;PEPB", "PEPX"],
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
                true,
                0.5f0,
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

    @testset "small training folds use empirical global protein scores" begin
        mktempdir() do directory
            run_scores = (
                Float32[0.9, 0.4, 0.8, 0.3],
                Float32[0.8, 0.3, 0.7, 0.2],
            )
            protein_names = ["T0", "D0", "T1", "D1"]
            targets = Bool[true, false, true, false]
            refs = Pioneer.ProteinGroupFileReference[]
            for (run_idx, scores) in enumerate(run_scores)
                path = joinpath(directory, "run$(run_idx).arrow")
                Arrow.write(path, DataFrame(
                    protein_name = protein_names,
                    target = targets,
                    entrap_id = zeros(UInt8, 4),
                    file_idx = fill(UInt32(run_idx), 4),
                    pg_score = scores,
                    n_peptides = fill(2, 4),
                    n_possible_unique_peptides = fill(4, 4),
                    peptide_list = fill("PEPA;PEPB", 4),
                ))
                push!(refs, Pioneer.ProteinGroupFileReference(path))
            end

            protein_to_cv_fold = Dictionary{
                String,
                @NamedTuple{best_score::Float32, cv_fold::UInt8},
            }()
            for (protein_name, score, fold) in zip(
                protein_names,
                run_scores[1],
                UInt8[0, 0, 1, 1],
            )
                insert!(
                    protein_to_cv_fold,
                    protein_name,
                    (best_score = score, cv_fold = fold),
                )
            end

            scores, _ = Pioneer.build_global_protein_score_dicts(
                refs,
                protein_to_cv_fold,
                4,
                2,
                true,
                0.5f0,
            )

            @test scores == Dict(
                ("T0", true, UInt8(0)) => logodds(Float32[0.9, 0.8], 1),
                ("D0", false, UInt8(0)) => logodds(Float32[0.4, 0.3], 1),
                ("T1", true, UInt8(0)) => logodds(Float32[0.8, 0.7], 1),
                ("D1", false, UInt8(0)) => logodds(Float32[0.3, 0.2], 1),
            )
        end
    end

    @testset "global protein training requires one hundred of each class" begin
        function training_inputs(n_per_class_per_fold)
            run_scores = Dict{
                Tuple{String, Bool, UInt8},
                Vector{Pioneer.GlobalProteinRunScore},
            }()
            folds = Dict{Tuple{String, Bool, UInt8}, UInt8}()
            for fold in UInt8[0, 1]
                for class_idx in 1:n_per_class_per_fold
                    target_key = ("T$(fold)_$(class_idx)", true, UInt8(0))
                    decoy_key = ("D$(fold)_$(class_idx)", false, UInt8(0))
                    run_scores[target_key] = Pioneer.GlobalProteinRunScore[
                        Pioneer.GlobalProteinRunScore(UInt32(1), 0.9f0),
                    ]
                    run_scores[decoy_key] = Pioneer.GlobalProteinRunScore[
                        Pioneer.GlobalProteinRunScore(UInt32(1), 0.1f0),
                    ]
                    folds[target_key] = fold
                    folds[decoy_key] = fold
                end
            end
            return Pioneer.GlobalProteinInputs(
                run_scores,
                Dict{Tuple{String, Bool, UInt8}, Set{String}}(),
                Dict{Tuple{String, Bool, UInt8}, Int}(),
                Dict{Tuple{String, Bool, UInt8}, Int}(),
                folds,
            )
        end

        @test !Pioneer._global_protein_training_data_sufficient(
            training_inputs(99),
            UInt8[0, 1],
        )
        @test Pioneer._global_protein_training_data_sufficient(
            training_inputs(100),
            UInt8[0, 1],
        )
        @test !Pioneer._global_protein_model_eligible(
            training_inputs(100),
            UInt8[0, 1],
            false,
        )
        @test Pioneer._global_protein_model_eligible(
            training_inputs(100),
            UInt8[0, 1],
            true,
        )
    end

    @testset "protein features support out-of-fold LightGBM scoring" begin
        run_scores = Dict{
            Tuple{String, Bool, UInt8},
            Vector{Pioneer.GlobalProteinRunScore},
        }()
        observed_peptides = Dict{Tuple{String, Bool, UInt8}, Set{String}}()
        max_n_peptides = Dict{Tuple{String, Bool, UInt8}, Int}()
        n_possible_unique_peptides = Dict{Tuple{String, Bool, UInt8}, Int}()
        folds = Dict{Tuple{String, Bool, UInt8}, UInt8}()
        for group_idx in 1:40
            fold = UInt8(group_idx % 2)
            target_key = ("target_$(group_idx)", true, UInt8(0))
            decoy_key = ("decoy_$(group_idx)", false, UInt8(0))
            run_scores[target_key] = Pioneer.GlobalProteinRunScore[
                Pioneer.GlobalProteinRunScore(
                    UInt32(1),
                    Float32(0.95 - 0.001 * group_idx),
                ),
                Pioneer.GlobalProteinRunScore(UInt32(2), 0.9f0),
            ]
            run_scores[decoy_key] = Pioneer.GlobalProteinRunScore[
                Pioneer.GlobalProteinRunScore(
                    UInt32(1),
                    Float32(0.05 + 0.001 * group_idx),
                ),
                Pioneer.GlobalProteinRunScore(UInt32(2), 0.1f0),
            ]
            observed_peptides[target_key] = Set(["PEPA", "PEPB"])
            observed_peptides[decoy_key] = Set(["PEPX", "PEPY"])
            max_n_peptides[target_key] = 2
            max_n_peptides[decoy_key] = 2
            n_possible_unique_peptides[target_key] = 4
            n_possible_unique_peptides[decoy_key] = 4
            folds[target_key] = fold
            folds[decoy_key] = fold
        end
        table = Pioneer._build_global_protein_feature_table(
            Pioneer.GlobalProteinInputs(
                run_scores,
                observed_peptides,
                max_n_peptides,
                n_possible_unique_peptides,
                folds,
            ),
            2,
            0.5f0,
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
