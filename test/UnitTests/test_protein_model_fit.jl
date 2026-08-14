# Tests for src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/model_fit.jl.
#
# Covers run-level protein LightGBM fitting, monotonic constraints,
# semi-supervised stopping, and protein-specific negative mining.

using DataFrames
using Statistics: mean

using Pioneer: build_protein_semisupervised_training_set
using Pioneer: build_lightgbm_classifier
using Pioneer: fit_protein_lightgbm_semisupervised, lightgbm_predict
using Pioneer: RUN_LEVEL_PROTEIN_LGBM_HP
using Pioneer: _protein_lightgbm_monotone_constraints
import Pioneer: DEBUG_CONSOLE_LEVEL

@testset "ProteinScoring/model_fit" begin

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
        @test model.booster.monotone_constraints == Int[0, 1, 1]
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
            :all_peptide_coverage_logit,
            :any_common_peps,
            :coverage_log_ratio,
            :precursor_consensus_prefix_shape,
            :shared_precursor_consensus_prefix_shape,
            :single_non_mbr_peptide,
            :single_non_mbr_prefix_shape,
            :mbr_recovered_peptides,
            :mbr_only_protein,
        ]
        constraints =
            _protein_lightgbm_monotone_constraints(feature_names)
        @test constraints == Int[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, -1]

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

    @testset "negative mining uses non-MBR support" begin
        scores = Float32[0.9, 0.8, 0.2, 0.1, 0.05, 0.3]
        targets = Bool[true, false, true, true, false, true]
        neutral_shape = zeros(Float32, 6)
        low_shapes = Float32[-0.3, 0.0, -0.3, 0.0, 0.0, 0.0]

        all_scores_weak = build_protein_semisupervised_training_set(
            scores,
            targets,
            neutral_shape,
            fill(Int64(2), 6);
            mined_negative_pep_threshold = 0.0f0,
        )
        @test all_scores_weak.mined_negative_mask ==
            Bool[true, false, true, true, false, true]

        low_shape_one_non_mbr = build_protein_semisupervised_training_set(
            scores,
            targets,
            low_shapes,
            ones(Int64, 6);
            mined_negative_pep_threshold = 1.1f0,
        )
        @test low_shape_one_non_mbr.mined_negative_mask ==
            Bool[true, false, true, false, false, false]

        low_shape_zero_non_mbr =
            build_protein_semisupervised_training_set(
                scores,
                targets,
                low_shapes,
                zeros(Int64, 6);
                mined_negative_pep_threshold = 1.1f0,
            )
        @test low_shape_zero_non_mbr.mined_negative_mask ==
            Bool[true, false, true, false, false, false]

        low_shape_too_many_non_mbr = build_protein_semisupervised_training_set(
            scores,
            targets,
            low_shapes,
            fill(Int64(2), 6);
            mined_negative_pep_threshold = 1.1f0,
        )
        @test !any(low_shape_too_many_non_mbr.mined_negative_mask)
    end

end
