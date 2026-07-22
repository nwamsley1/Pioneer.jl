# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pioneer.jl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

@testset "Global precursor scoring" begin
    @testset "log-odds score summary" begin
        @test logodds([0.8], 1) == 0.8
        @test 0.0 < logodds([0.0, 0.0], 2) < 1.0
        @test 0.0 < logodds([1.0, 1.0], 2) < 1.0

        probabilities = [0.9, 0.8, 0.7, 0.6]
        @test logodds(probabilities, 2) ≈ Pioneer._logodds_from_sorted(
            sort(probabilities; rev = true),
            2,
        )
    end

    @testset "score-distribution feature table" begin
        inputs = Pioneer.GlobalPrecursorInputs(
            Dict(
                UInt32(1) => Float32[0.9, 0.8, 0.6],
                UInt32(2) => Float32[0.7],
            ),
            Dict(UInt32(1) => true, UInt32(2) => false),
            Dict(UInt32(1) => UInt8(0), UInt32(2) => UInt8(1)),
        )
        table = Pioneer._build_global_precursor_feature_table(inputs, 4)

        @test table.precursor_idx == UInt32[1, 2]
        @test table.target == Bool[true, false]
        @test table.cv_fold == UInt8[0, 1]
        @test Symbol.(names(table)[4:end]) == Pioneer.GLOBAL_PRECURSOR_SCORE_FEATURES

        @test table.empirical_global_score[1] ≈ logodds(Float32[0.9, 0.8, 0.6], 2)
        @test table.top1_prec_prob[1] == 0.9f0
        @test table.top2_prec_prob[1] == 0.8f0
        @test table.top3_prec_prob[1] == 0.6f0
        @test table.top2_logodds_score[1] ≈ logodds(Float32[0.9, 0.8, 0.6], 2)
        @test table.top3_logodds_score[1] ≈ logodds(Float32[0.9, 0.8, 0.6], 3)
        @test table.mean_prec_prob[1] ≈ mean(Float32[0.9, 0.8, 0.6])
        @test table.median_prec_prob[1] == 0.8f0
        @test table.std_prec_prob[1] ≈ std(Float32[0.9, 0.8, 0.6])
        @test table.min_prec_prob[1] == 0.6f0
        @test table.top1_top2_gap[1] ≈ 0.1f0
        @test table.top2_top3_gap[1] ≈ 0.2f0
        @test table.n_runs_observed[1] == 3.0f0
        @test table.n_prob_gt_0_5[1] == 3.0f0
        @test table.n_prob_gt_0_9[1] == 0.0f0
        @test table.n_prob_gt_0_99[1] == 0.0f0

        @test table.n_runs_observed[2] == 1.0f0
        @test table.top3_prec_prob[2] == 0.0f0
        @test table.median_prec_prob[2] == 0.7f0
        @test table.std_prec_prob[2] == 0.0f0
        @test table.min_prec_prob[2] == 0.7f0
        @test table.top1_top2_gap[2] == 0.0f0
        @test table.top2_top3_gap[2] == 0.0f0
        @test table.n_prob_gt_0_5[2] == 1.0f0
    end

    @testset "run-level input collection" begin
        mktempdir() do directory
            run1_path = joinpath(directory, "run1.arrow")
            run2_path = joinpath(directory, "run2.arrow")
            Arrow.write(run1_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                prec_prob = Float32[0.9, 0.7],
                target = Bool[true, false],
                cv_fold = UInt8[0, 1],
            ))
            Arrow.write(run2_path, DataFrame(
                precursor_idx = UInt32[1, 3],
                prec_prob = Float32[0.8, 0.6],
                target = Bool[true, false],
                cv_fold = UInt8[0, 1],
            ))
            refs = PSMFileReference[
                PSMFileReference(run1_path),
                PSMFileReference(run2_path),
            ]

            inputs = Pioneer._collect_global_precursor_inputs(refs, 3)

            @test inputs.probabilities[UInt32(1)] == Float32[0.9, 0.8]
            @test inputs.probabilities[UInt32(2)] == Float32[0.7]
            @test inputs.probabilities[UInt32(3)] == Float32[0.6]
            @test inputs.targets == Dict(
                UInt32(1) => true,
                UInt32(2) => false,
                UInt32(3) => false,
            )
            @test inputs.folds == Dict(
                UInt32(1) => UInt8(0),
                UInt32(2) => UInt8(1),
                UInt32(3) => UInt8(1),
            )
        end
    end

    @testset "single run uses individual precursor probabilities" begin
        mktempdir() do directory
            run_path = joinpath(directory, "single_run.arrow")
            Arrow.write(run_path, DataFrame(
                precursor_idx = UInt32[1, 2, 3],
                prec_prob = Float32[0.02, 0.5, 0.98],
                target = Bool[false, true, true],
            ))
            refs = PSMFileReference[PSMFileReference(run_path)]

            scores, targets = Pioneer.build_global_precursor_score_dicts(
                refs,
                3,
                1,
                1.0f0,
            )

            @test scores == Dict(
                UInt32(1) => 0.02f0,
                UInt32(2) => 0.5f0,
                UInt32(3) => 0.98f0,
            )
            @test targets == Dict(
                UInt32(1) => false,
                UInt32(2) => true,
                UInt32(3) => true,
            )
        end
    end

    @testset "two-fold out-of-fold LightGBM scoring" begin
        precursor_idx = UInt32[]
        target = Bool[]
        cv_fold = UInt8[]
        direct_evidence = Float32[]
        for group_idx in 1:40
            fold = UInt8(group_idx % 2)
            push!(precursor_idx, UInt32(2 * group_idx - 1), UInt32(2 * group_idx))
            push!(target, true, false)
            push!(cv_fold, fold, fold)
            push!(direct_evidence, 2.0f0 + 0.01f0 * group_idx)
            push!(direct_evidence, -2.0f0 - 0.01f0 * group_idx)
        end
        table = DataFrame(
            precursor_idx = precursor_idx,
            target = target,
            cv_fold = cv_fold,
            direct_evidence = direct_evidence,
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
            Symbol[:direct_evidence],
            UInt8[0, 1];
            scoring_name = "Global precursor",
            lgbm_hp = test_hyperparameters,
            min_training_class_count = 1,
            max_train = 1_000,
        )

        @test length(scored.scores) == nrow(table)
        @test length(scored.models) == 2
        @test scored.iter == 1
        @test mean(scored.scores[target]) > mean(scored.scores[.!target])
    end
end
