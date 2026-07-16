# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

@testset "Scoring Interface Tests" begin

    @testset "Logodds Function" begin
        # Test single probability
        @test logodds([0.8], 1) == 0.8

        # Test multiple probabilities
        probs = [0.9, 0.8, 0.7, 0.6]
        result = logodds(probs, 2)  # Take top 2
        @test result isa Float64
        @test 0.0 <= result <= 1.0

        # Test edge cases
        @test 0.0 < logodds([0.0, 0.0], 2) < 1.0  # Should handle zeros
        @test 0.0 < logodds([1.0, 1.0], 2) < 1.0  # Should handle ones

        # Test empty case
        @test logodds(Float64[], 1) == 0.0f0
    end

    @testset "Global precursor features summarize run-level scores" begin
        mktempdir() do dir
            run1_path = joinpath(dir, "run1.arrow")
            run2_path = joinpath(dir, "run2.arrow")
            Arrow.write(run1_path, DataFrame(
                precursor_idx = UInt32[1, 2],
                ms_file_idx = UInt32[1, 1],
                prec_prob = Float32[0.9, 0.7],
                weight = Float32[100.0, 25.0],
                target = Bool[true, false],
                cv_fold = UInt8[0, 1],
            ))
            Arrow.write(run2_path, DataFrame(
                precursor_idx = UInt32[1, 3],
                ms_file_idx = UInt32[2, 2],
                prec_prob = Float32[0.8, 0.6],
                weight = Float32[80.0, 20.0],
                target = Bool[true, false],
                cv_fold = UInt8[0, 1],
            ))
            refs = Pioneer.PSMFileReference[
                Pioneer.PSMFileReference(run1_path),
                Pioneer.PSMFileReference(run2_path),
            ]
            score_similarity = Pioneer.build_mbr_run_similarity_from_score_floor(
                refs,
                0.65f0,
            )
            @test score_similarity.total_weight_by_file[UInt32(1)] > 0.0f0
            @test score_similarity.total_weight_by_file[UInt32(2)] == 0.0f0
            @test score_similarity.cluster_by_precursor == UInt32[1, 1]
            @test score_similarity.cluster_n_runs_observed == UInt16[2, 1]

            inputs = Pioneer._collect_global_precursor_inputs(refs, 3)
            run_similarity = Pioneer._build_mbr_run_similarity_from_passed(Dict(
                UInt32(1) => BitSet((100, 101)),
                UInt32(2) => BitSet((100, 101)),
                UInt32(3) => BitSet((102,)),
            ))
            built = Pioneer._build_global_precursor_feature_table(
                inputs;
                sqrt_n_runs = 2,
                n_runs_total = 4,
                run_similarity = run_similarity,
            )

            row1 = findfirst(==(UInt32(1)), built.table.precursor_idx)
            row2 = findfirst(==(UInt32(2)), built.table.precursor_idx)
            @test built.features == Pioneer.GLOBAL_PRECURSOR_SCORE_FEATURES
            @test built.table.empirical_global_score[row1] ≈
                  Pioneer.logodds(Float32[0.9, 0.8], 2)
            @test built.table.top1_prec_prob[row1] == 0.9f0
            @test built.table.top2_prec_prob[row1] == 0.8f0
            @test built.table.top3_prec_prob[row1] == 0.0f0
            @test built.table.top2_logodds_score[row1] ≈
                  Pioneer.logodds(Float32[0.9, 0.8], 2)
            @test built.table.top3_logodds_score[row1] ≈
                  Pioneer.logodds(Float32[0.9, 0.8], 3)
            @test built.table.mean_prec_prob[row1] ≈ 0.85f0
            @test built.table.median_prec_prob[row1] ≈ 0.85f0
            @test built.table.std_prec_prob[row1] ≈ sqrt(0.005f0)
            @test built.table.min_prec_prob[row1] == 0.8f0
            @test built.table.top1_top2_gap[row1] ≈ 0.1f0
            @test built.table.top2_top3_gap[row1] == 0.0f0
            @test built.table.n_runs_observed[row1] == 2.0f0
            @test built.table.n_runs_passing_local_q[row1] == 2.0f0
            @test !hasproperty(built.table, :observed_run_fraction)
            @test built.table.n_prob_gt_0_5[row1] == 2.0f0
            @test built.table.n_prob_gt_0_9[row1] == 0.0f0
            @test built.table.n_prob_gt_0_99[row1] == 0.0f0
            @test built.table.observed_run_centrality_mean[row1] ≈ 1.0f0 / 3.0f0
            @test built.table.observed_run_centrality_max[row1] ≈ 1.0f0 / 3.0f0
            @test built.table.missing_run_similarity_mass_approx[row1] ≈ 2.0f0 / 3.0f0
            @test built.table.missing_run_similarity_mass_approx[row2] ≈ 1.0f0
            @test built.table.std_prec_prob[row2] == 0.0f0
            @test built.table.top1_top2_gap[row2] == 0.0f0
            @test !(:fitted_manhattan_distance_max in built.features)

            passing_only = Pioneer._build_global_precursor_feature_table(
                inputs;
                sqrt_n_runs = 2,
                n_runs_total = 4,
                run_similarity = run_similarity,
                run_score_floor = 0.85f0,
            ).table
            @test passing_only.n_runs_passing_local_q[row1] == 1.0f0
            @test passing_only.missing_run_similarity_mass_approx[row1] ≈ 1.0f0
            @test passing_only.n_runs_passing_local_q[row2] == 0.0f0
            @test passing_only.missing_run_similarity_mass_approx[row2] == 0.0f0
        end
    end

    @testset "Global precursor LightGBM produces out-of-fold scores" begin
        n_groups = 40
        precursor_idx = UInt32[]
        target = Bool[]
        cv_fold = UInt8[]
        direct_evidence = Float32[]
        for group_idx in 1:n_groups
            fold = UInt8(group_idx % 2)
            push!(precursor_idx, UInt32(2 * group_idx - 1), UInt32(2 * group_idx))
            push!(target, true, false)
            push!(cv_fold, fold, fold)
            push!(direct_evidence, 2.0f0 + 0.01f0 * group_idx, -2.0f0 - 0.01f0 * group_idx)
        end
        table = DataFrame(
            precursor_idx = precursor_idx,
            target = target,
            cv_fold = cv_fold,
            direct_evidence = direct_evidence,
        )
        test_hp = (
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
        scored = Pioneer._score_global_precursor_features_oof(
            table,
            Symbol[:direct_evidence];
            lgbm_hp = test_hp,
            min_training_class_count = 1,
            max_train = 1_000,
            max_iterations = 3,
        )

        @test scored !== nothing
        @test length(scored.models) == 2
        @test scored.iter in 1:3
        @test all(model -> Pioneer.importance(model) !== nothing, scored.models)
        @test sum(scored.scores[target]) / count(target) >
              sum(scored.scores[.!target]) / count(.!target)
        old_debug_level = Pioneer.DEBUG_CONSOLE_LEVEL[]
        debug_text = mktemp() do _, io
            try
                Pioneer.DEBUG_CONSOLE_LEVEL[] = 1
                redirect_stdout(io) do
                    Pioneer._log_global_precursor_feature_importances(
                        scored.models,
                        Symbol[:direct_evidence],
                    )
                end
                flush(io)
                seekstart(io)
                read(io, String)
            finally
                Pioneer.DEBUG_CONSOLE_LEVEL[] = old_debug_level
            end
        end
        @test occursin(
            "Global precursor LightGBM feature gains (summed across folds)",
            debug_text,
        )

        mktempdir() do dir
            path = joinpath(dir, "global_model.arrow")
            Arrow.write(path, DataFrame(
                precursor_idx = precursor_idx,
                prec_prob = Float32[is_target ? 0.95f0 : 0.05f0 for is_target in target],
                target = target,
                cv_fold = cv_fold,
            ))
            refs = Pioneer.PSMFileReference[Pioneer.PSMFileReference(path)]
            global_scores, global_targets = Pioneer.build_precursor_global_prob_dicts(
                refs,
                1,
                length(precursor_idx);
                n_runs_total = 1,
                lgbm_hp = test_hp,
                min_training_class_count = 1,
                max_train = 1_000,
                max_iterations = 3,
            )
            target_scores = Float32[global_scores[pid] for (pid, is_target) in global_targets if is_target]
            decoy_scores = Float32[global_scores[pid] for (pid, is_target) in global_targets if !is_target]
            @test length(global_scores) == length(precursor_idx)
            @test sum(target_scores) / length(target_scores) >
                  sum(decoy_scores) / length(decoy_scores)
        end
    end

    @testset "Global precursor scorer retains empirical legacy fallback" begin
        mktempdir() do dir
            path = joinpath(dir, "legacy.arrow")
            Arrow.write(path, DataFrame(
                precursor_idx = UInt32[1, 2],
                prec_prob = Float32[0.9, 0.2],
                target = Bool[true, false],
                fitted_manhattan_distance = Float32[1, 5],
            ))
            refs = Pioneer.PSMFileReference[Pioneer.PSMFileReference(path)]
            scores, targets = Pioneer.build_precursor_global_prob_dicts(refs, 1, 2)
            @test scores[UInt32(1)] ≈ Pioneer.logodds(Float32[0.9], 1)
            @test scores[UInt32(2)] ≈ Pioneer.logodds(Float32[0.2], 1)
            @test targets == Dict{UInt32, Bool}(UInt32(1) => true, UInt32(2) => false)
        end
    end
end
