using Test
using DataFrames
using Statistics
using Arrow
using Distributions
using Pioneer

@testset "Protein Group Weight Coverage Features" begin
    empty_precursor_consensus = (
        relative_weight = Dict{Tuple{String, Bool, UInt8, UInt32}, Float32}(),
        mean_relative_weight = Dict{Tuple{String, Bool, UInt8}, Float32}(),
        profiled_precursor_count = Dict{Tuple{String, Bool, UInt8}, Int32}()
    )

    @testset "Quant Necessary Columns Include Weight" begin
        cols_nombr = Pioneer.get_quant_necessary_columns(false)
        cols_mbr = Pioneer.get_quant_necessary_columns(true)

        @test :weight in cols_nombr
        @test :weight in cols_mbr
    end

    @testset "Group PSMs By Protein Uses Max Weight Within Peptide" begin
        psms = DataFrame(
            inferred_protein_group = ["P1", "P1", "P1", "P2"],
            target = Bool[true, true, true, false],
            entrap_id = UInt8[1, 1, 1, 1],
            precursor_idx = UInt32[11, 12, 13, 21],
            sequence = ["PEP1", "PEP1", "PEP2", "PEP3"],
            use_for_protein_quant = Bool[true, true, true, true],
            missed_cleavage = Int64[0, 1, 0, 0],
            Mox = Int64[0, 0, 0, 0],
            prec_prob = Float32[0.8, 0.9, 0.6, 0.7],
            weight = Float32[100.0, 80.0, 50.0, 25.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = empty_precursor_consensus)

        @test hasproperty(grouped, :top_pep_weight)

        p1 = grouped[(grouped.protein_name .== "P1") .& grouped.target, :]
        @test nrow(p1) == 1
        @test p1.n_peptides[1] == 2
        @test p1.top_pep_weight[1] == 100.0f0
    end

    @testset "Grouped Protein Catalog Uses Union Peptide Set" begin
        protein_catalog = Dict(
            (protein_name = "A", target = true, entrap_id = UInt8(1)) => Set(["p1", "p2"]),
            (protein_name = "B", target = true, entrap_id = UInt8(1)) => Set(["p2", "p3"])
        )
        df = DataFrame(
            protein_name = ["A;B"],
            target = Bool[true],
            entrap_id = UInt8[1],
            n_peptides = Int64[2]
        )

        (_, op) = Pioneer.add_protein_features(protein_catalog)
        out = op(copy(df))

        @test out.n_possible_peptides[1] == 3
        @test out.peptide_coverage[1] ≈ (2f0 / 3f0)
    end

    @testset "Weight Calibration Uses Max Peptide Weights" begin
        psms = DataFrame(
            use_for_protein_quant = Bool[true, true, true, true, true, false],
            sequence = ["P1A", "P1A", "P1B", "P2A", "P2B", "P3A"],
            inferred_protein_group = Union{Missing, String}["P1", "P1", "P1", "P2", "P2", missing],
            target = Bool[true, true, true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1, 1, 1],
            weight = Float32[60.0f0, 40.0f0, 40.0f0, 80.0f0, 10.0f0, 500.0f0]
        )

        model = Pioneer.estimate_weight_detection_model(psms)

        expected_threshold = Float32(quantile(log.([60.0, 40.0, 80.0, 10.0]), 0.05))
        @test model.log_threshold ≈ expected_threshold
        @test length(model.rank_scale_profile) == model.profiled_rank_count
        @test all(x -> x > 0.0f0, model.rank_scale_profile)
        @test model.profiled_rank_count == 1
    end

    @testset "Coverage Match From Top Behaves as Expected" begin
        calibration = (
            log_threshold = Float32(log(10.0)),
            rank_drop_profile = Float32[0.0f0, -0.1f0, -0.8f0, -1.5f0],
            rank_scale_profile = Float32[1.0f0, 0.25f0, 0.35f0, 0.45f0],
            profiled_rank_count = 4
        )
        pg_df = DataFrame(
            protein_name = ["P_big", "P_small", "P_balanced", "P_invalid", "P_short"],
            target = Bool[true, true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1, 1],
            n_peptides = Int64[1, 1, 4, 1, 1],
            pg_score = Float32[2.0, 1.5, 0.5, 3.0, 1.0],
            n_possible_peptides = Int64[20, 2, 4, 12, 1],
            top_pep_weight = Float32[100.0, 100.0, 100.0, 0.0, 100.0]
        )

        (_, op) = Pioneer.add_weight_observation_features(calibration)
        out = op(copy(pg_df))
        (_, interaction_op) = Pioneer.add_pg_score_interaction_features()
        out = interaction_op(out)

        for col in (
            :expected_additional_from_top,
            :coverage_match_from_top,
            :pg_score_x_coverage_match_from_top
        )
            @test hasproperty(out, col)
        end

        @test out.expected_additional_from_top[1] > out.expected_additional_from_top[2]
        @test out.coverage_match_from_top[1] < out.coverage_match_from_top[2]
        @test out.coverage_match_from_top[3] ≈ 1.0f0

        # Invalid weight -> neutral values
        @test out.expected_additional_from_top[4] == 0.0f0
        @test out.coverage_match_from_top[4] == 1.0f0

        # N <= 1 -> neutral values
        @test out.expected_additional_from_top[5] == 0.0f0
        @test out.coverage_match_from_top[5] == 1.0f0

        expected_cov_match_1 = out.pg_score[1] * out.coverage_match_from_top[1]
        expected_cov_match_3 = out.pg_score[3] * out.coverage_match_from_top[3]

        @test out.pg_score_x_coverage_match_from_top[1] ≈ expected_cov_match_1
        @test out.pg_score_x_coverage_match_from_top[3] ≈ expected_cov_match_3
    end

    @testset "pg_score Coverage Match Interaction Preserves Good to Bad Ordering" begin
        df = DataFrame(
            pg_score = Float32[100.0, 100.0, 10.0, 10.0],
            coverage_match_from_top = Float32[0.9, 0.1, 0.9, 0.1]
        )

        (_, interaction_op) = Pioneer.add_pg_score_interaction_features()
        out = interaction_op(copy(df))

        @test out.pg_score_x_coverage_match_from_top[1] > out.pg_score_x_coverage_match_from_top[2]
        @test out.pg_score_x_coverage_match_from_top[2] > out.pg_score_x_coverage_match_from_top[3]
        @test out.pg_score_x_coverage_match_from_top[3] > out.pg_score_x_coverage_match_from_top[4]
    end

    @testset "Protein Semi-Supervised Training Set Mines Negatives by PEP" begin
        scores = Float32[
            0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90,
            0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.80,
            0.79, 0.78, 0.77, 0.76, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15,
            0.14, 0.13, 0.12, 0.11, 0.10, 0.09
        ]
        targets = Bool[
            true, true, true, true, true, true, true, true, true, true,
            true, true, true, true, true, true, true, true, false, false,
            false, false, true, true, false, false, false, false, false, false,
            false, false, false, false, false, false
        ]

        ss = Pioneer.build_protein_semisupervised_training_set(
            scores,
            targets;
            q_value_threshold = 0.01f0,
            min_pep_neg_threshold = 0.90f0
        )

        @test any(ss.positive_mask)
        @test any(ss.mined_negative_mask)
        @test all(ss.keep_mask[.!targets])
        @test sum(ss.labels) < sum(ss.positive_mask) + sum(ss.mined_negative_mask)
        @test sum(ss.labels) == sum(ss.positive_mask[ss.keep_mask])
        @test any(.!ss.labels)
    end

    @testset "Consensus Relative Weight Builder Ignores MBR Candidates and Uses Ranked pg_score Weighting" begin
        temp_dir = mktempdir()

        try
            file1 = DataFrame(
                inferred_protein_group = ["P", "P", "P"],
                target = Bool[true, true, true],
                entrap_id = UInt8[1, 1, 1],
                precursor_idx = UInt32[1, 2, 3],
                sequence = ["PEP1", "PEP2", "PEP3"],
                use_for_protein_quant = Bool[true, true, true],
                MBR_candidate = Bool[false, false, true],
                missed_cleavage = Int64[0, 0, 0],
                Mox = Int64[0, 0, 0],
                prec_prob = Float32[0.99, 0.98, 0.10],
                MBR_boosted_prec_prob = Float32[0.99, 0.98, 0.995],
                weight = Float32[100.0, 50.0, 1000.0]
            )

            file2 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[1, 2],
                sequence = ["PEP1", "PEP2"],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.10, 0.20],
                MBR_boosted_prec_prob = Float32[0.10, 0.20],
                weight = Float32[60.0, 120.0]
            )

            path1 = joinpath(temp_dir, "psms_001.arrow")
            path2 = joinpath(temp_dir, "psms_002.arrow")
            Arrow.write(path1, file1)
            Arrow.write(path2, file2)

            consensus = Pioneer.build_precursor_consensus([
                Pioneer.PSMFileReference(path1),
                Pioneer.PSMFileReference(path2)
            ])

            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(1))] >
                  consensus.relative_weight[("P", true, UInt8(1), UInt32(2))]
            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(1))] > 0.9f0
            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(2))] > 0.4f0
            @test !haskey(consensus.relative_weight, ("P", true, UInt8(1), UInt32(3)))
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Builder Keeps Only Top Five Protein Runs" begin
        temp_dir = mktempdir()

        try
            refs = Pioneer.PSMFileReference[]
            for run_idx in 1:6
                df = if run_idx <= 5
                    DataFrame(
                        inferred_protein_group = ["P"],
                        target = Bool[true],
                        entrap_id = UInt8[1],
                        precursor_idx = UInt32[1],
                        sequence = ["PEP1"],
                        use_for_protein_quant = Bool[true],
                        MBR_candidate = Bool[false],
                        missed_cleavage = Int64[0],
                        Mox = Int64[0],
                        prec_prob = Float32[0.95f0 - 0.01f0 * run_idx],
                        weight = Float32[100.0]
                    )
                else
                    DataFrame(
                        inferred_protein_group = ["P"],
                        target = Bool[true],
                        entrap_id = UInt8[1],
                        precursor_idx = UInt32[2],
                        sequence = ["PEP2"],
                        use_for_protein_quant = Bool[true],
                        MBR_candidate = Bool[false],
                        missed_cleavage = Int64[0],
                        Mox = Int64[0],
                        prec_prob = Float32[0.05],
                        weight = Float32[100.0]
                    )
                end

                path = joinpath(temp_dir, "psms_$(lpad(run_idx, 3, '0')).arrow")
                Arrow.write(path, df)
                push!(refs, Pioneer.PSMFileReference(path))
            end

            consensus = Pioneer.build_precursor_consensus(refs)

            @test consensus.profiled_precursor_count[("P", true, UInt8(1))] == 1
            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(1))] == 1.0f0
            @test !haskey(consensus.relative_weight, ("P", true, UInt8(1), UInt32(2)))
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Builder Applies Exponential Decay Across Selected Runs" begin
        temp_dir = mktempdir()

        try
            run1 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[2, 1],
                sequence = ["PEP2", "PEP1"],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.90, 0.90],
                weight = Float32[100.0, 80.0]
            )
            run2 = DataFrame(
                inferred_protein_group = ["P", "P"],
                target = Bool[true, true],
                entrap_id = UInt8[1, 1],
                precursor_idx = UInt32[1, 2],
                sequence = ["PEP1", "PEP2"],
                use_for_protein_quant = Bool[true, true],
                MBR_candidate = Bool[false, false],
                missed_cleavage = Int64[0, 0],
                Mox = Int64[0, 0],
                prec_prob = Float32[0.90, 0.90],
                weight = Float32[100.0, 80.0]
            )

            path1 = joinpath(temp_dir, "psms_001.arrow")
            path2 = joinpath(temp_dir, "psms_002.arrow")
            Arrow.write(path1, run1)
            Arrow.write(path2, run2)

            consensus = Pioneer.build_precursor_consensus([
                Pioneer.PSMFileReference(path1),
                Pioneer.PSMFileReference(path2)
            ])

            @test consensus.relative_weight[("P", true, UInt8(1), UInt32(2))] >
                  consensus.relative_weight[("P", true, UInt8(1), UInt32(1))]
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Relative Weight Support Uses Observed MBR Precursors" begin
        consensus = (
            relative_weight = Dict(
                ("P", true, UInt8(1), UInt32(11)) => 1.0f0,
                ("P", true, UInt8(1), UInt32(12)) => 0.2f0
            ),
            mean_relative_weight = Dict(
                ("P", true, UInt8(1)) => 0.4f0
            ),
            profiled_precursor_count = Dict(
                ("P", true, UInt8(1)) => Int32(5)
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P", "P"],
            target = Bool[true, true],
            entrap_id = UInt8[1, 1],
            precursor_idx = UInt32[11, 12],
            sequence = ["PEP1", "PEP2"],
            use_for_protein_quant = Bool[true, true],
            MBR_candidate = Bool[false, true],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            prec_prob = Float32[0.8, 0.7],
            MBR_boosted_prec_prob = Float32[0.8, 0.85],
            weight = Float32[100.0, 90.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)
        @test grouped.precursor_consensus_support[1] ≈ 1.2f0
        @test grouped.precursor_consensus_enrichment[1] ≈ (1.2f0 / (0.4f0 * 2.0f0)) atol = 1e-5

        (_, interaction_op) = Pioneer.add_precursor_consensus_interaction_features()
        out = interaction_op(copy(grouped))
        @test out.pg_score_x_precursor_consensus_support[1] ==
              out.pg_score[1] * out.precursor_consensus_support[1]
        @test out.pg_score_x_precursor_consensus_enrichment[1] ==
              out.pg_score[1] * out.precursor_consensus_enrichment[1]
    end

    @testset "Consensus Relative Weight Enrichment Rewards Top Support Among Many Possibilities" begin
        consensus = (
            relative_weight = Dict(
                ("P_many", true, UInt8(1), UInt32(101)) => 1.0f0,
                ("P_one", true, UInt8(1), UInt32(201)) => 1.0f0
            ),
            mean_relative_weight = Dict(
                ("P_many", true, UInt8(1)) => 0.2f0,
                ("P_one", true, UInt8(1)) => 1.0f0
            ),
            profiled_precursor_count = Dict(
                ("P_many", true, UInt8(1)) => Int32(10),
                ("P_one", true, UInt8(1)) => Int32(1)
            )
        )

        psms = DataFrame(
            inferred_protein_group = ["P_many", "P_one"],
            target = Bool[true, true],
            entrap_id = UInt8[1, 1],
            precursor_idx = UInt32[101, 201],
            sequence = ["PEP_MANY", "PEP_ONE"],
            use_for_protein_quant = Bool[true, true],
            MBR_candidate = Bool[false, false],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            prec_prob = Float32[0.8, 0.8],
            weight = Float32[100.0, 100.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms; precursor_consensus = consensus)

        many = grouped[grouped.protein_name .== "P_many", :]
        one = grouped[grouped.protein_name .== "P_one", :]

        @test many.precursor_consensus_support[1] == 1.0f0
        @test one.precursor_consensus_support[1] == 1.0f0
        @test many.precursor_consensus_enrichment[1] > one.precursor_consensus_enrichment[1]
    end

    @testset "Optional Probit Feature Columns Drop Cleanly" begin
        feature_names = Pioneer.protein_probit_feature_names()
        df = DataFrame(pg_score = Float32[0.1, 0.2, 0.3])

        Pioneer.remove_zero_variance_columns!(feature_names, df)

        @test isempty(feature_names)
    end

    @testset "Protein Probit Feature Names Include Consensus Support" begin
        @test Pioneer.protein_probit_feature_names() == [
            :peptide_coverage,
            :any_common_peps,
            :pg_score_x_coverage_match_from_top,
            :pg_score_x_precursor_consensus_support
        ]

        @test Pioneer.protein_probit_feature_names(include_n_possible_peptides = true) == [
            :peptide_coverage,
            :n_possible_peptides,
            :any_common_peps,
            :pg_score_x_coverage_match_from_top,
            :pg_score_x_precursor_consensus_support
        ]
    end

end
