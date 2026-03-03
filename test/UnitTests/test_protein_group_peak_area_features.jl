using Test
using DataFrames
using Statistics
using Arrow
using Pioneer

@testset "Protein Group Weight Coverage Features" begin
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
            sequence = ["PEP1", "PEP1", "PEP2", "PEP3"],
            use_for_protein_quant = Bool[true, true, true, true],
            missed_cleavage = Int64[0, 1, 0, 0],
            Mox = Int64[0, 0, 0, 0],
            prec_prob = Float32[0.8, 0.9, 0.6, 0.7],
            weight = Float32[100.0, 80.0, 50.0, 25.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms)

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
            weight = Float32[60.0f0, 40.0f0, 40.0f0, 80.0f0, 10.0f0, 500.0f0]
        )

        model = Pioneer.estimate_weight_detection_model(psms)

        expected_threshold = Float32(quantile(log.([60.0, 40.0, 80.0, 10.0]), 0.05))
        @test model.n_unique_peptides == 4
        @test model.log_threshold ≈ expected_threshold
        @test 0.25f0 <= model.sigma_log <= 2.5f0
        @test model.used_fallback == false
    end

    @testset "Coverage Surprise Features Behave as Expected" begin
        calibration = (log_threshold = Float32(log(10.0)), sigma_log = 1.0f0)
        pg_df = DataFrame(
            protein_name = ["P_big", "P_small", "P_balanced", "P_invalid", "P_short"],
            target = Bool[true, true, true, true, true],
            entrap_id = UInt8[1, 1, 1, 1, 1],
            n_peptides = Int64[1, 1, 6, 1, 1],
            pg_score = Float32[2.0, 1.5, 0.5, 3.0, 1.0],
            n_possible_peptides = Int64[20, 2, 11, 12, 1],
            top_pep_weight = Float32[100.0, 100.0, 10.0, 0.0, 100.0]
        )

        (_, op) = Pioneer.add_weight_observation_features(calibration)
        out = op(copy(pg_df))
        (_, interaction_op) = Pioneer.add_pg_score_interaction_features()
        out = interaction_op(out)

        for col in (
            :coverage_miss_pval,
            :coverage_miss_surprisal,
            :coverage_deficit_z,
            :top_weight_vs_threshold_z,
            :pg_score_x_coverage_miss_surprisal,
            :pg_score_x_coverage_deficit_z,
            :pg_score_x_top_weight_vs_threshold_z
        )
            @test hasproperty(out, col)
        end

        @test out.coverage_miss_pval[1] < out.coverage_miss_pval[2]
        @test out.coverage_miss_surprisal[1] > out.coverage_miss_surprisal[2]
        @test abs(out.coverage_deficit_z[3]) < 0.25f0
        @test out.top_weight_vs_threshold_z[1] > 0.0f0
        @test abs(out.top_weight_vs_threshold_z[3]) < 1e-5

        # Invalid weight -> neutral values
        @test out.coverage_miss_pval[4] == 1.0f0
        @test out.coverage_miss_surprisal[4] == 0.0f0
        @test out.coverage_deficit_z[4] == 0.0f0
        @test out.top_weight_vs_threshold_z[4] == 0.0f0

        # N <= 1 -> neutral values
        @test out.coverage_miss_pval[5] == 1.0f0
        @test out.coverage_miss_surprisal[5] == 0.0f0
        @test out.coverage_deficit_z[5] == 0.0f0
        @test out.top_weight_vs_threshold_z[5] == 0.0f0

        @test out.pg_score_x_coverage_miss_surprisal[1] == out.pg_score[1] * out.coverage_miss_surprisal[1]
        @test out.pg_score_x_coverage_deficit_z[1] == out.pg_score[1] * out.coverage_deficit_z[1]
        @test out.pg_score_x_top_weight_vs_threshold_z[1] == out.pg_score[1] * out.top_weight_vs_threshold_z[1]
        @test out.pg_score_x_coverage_miss_surprisal[3] == out.pg_score[3] * out.coverage_miss_surprisal[3]
        @test out.pg_score_x_coverage_deficit_z[3] == out.pg_score[3] * out.coverage_deficit_z[3]
        @test out.pg_score_x_top_weight_vs_threshold_z[3] == out.pg_score[3] * out.top_weight_vs_threshold_z[3]
    end

    @testset "Consensus Rank Builder Ignores MBR Candidates and Uses pg_score Weighting" begin
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

            consensus = Pioneer.build_precursor_rank_consensus([
                Pioneer.PSMFileReference(path1),
                Pioneer.PSMFileReference(path2)
            ])

            @test consensus.precursor_rank[("P", true, UInt8(1), UInt32(1))] == 1
            @test consensus.precursor_rank[("P", true, UInt8(1), UInt32(2))] == 2
            @test !haskey(consensus.precursor_rank, ("P", true, UInt8(1), UInt32(3)))
        finally
            rm(temp_dir, recursive = true, force = true)
        end
    end

    @testset "Consensus Rank Support Uses Observed MBR Precursors" begin
        consensus = (
            precursor_rank = Dict(
                ("P", true, UInt8(1), UInt32(11)) => Int32(1),
                ("P", true, UInt8(1), UInt32(12)) => Int32(5)
            ),
            protein_rank_count = Dict(
                ("P", true, UInt8(1)) => Int32(5)
            ),
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
        @test grouped.consensus_precursor_rank_support[1] ≈ 1.2f0
        @test grouped.consensus_precursor_rank_enrichment[1] ≈ (1.2f0 / ((2.2833333f0 / 5.0f0) * 2.0f0)) atol = 1e-5

        (_, interaction_op) = Pioneer.add_consensus_precursor_rank_interaction_feature()
        out = interaction_op(copy(grouped))
        @test out.pg_score_x_consensus_precursor_rank_support[1] ==
              out.pg_score[1] * out.consensus_precursor_rank_support[1]
        @test out.pg_score_x_consensus_precursor_rank_enrichment[1] ==
              out.pg_score[1] * out.consensus_precursor_rank_enrichment[1]
    end

    @testset "Consensus Rank Enrichment Rewards Top Rank Among Many Possibilities" begin
        consensus = (
            precursor_rank = Dict(
                ("P_many", true, UInt8(1), UInt32(101)) => Int32(1),
                ("P_one", true, UInt8(1), UInt32(201)) => Int32(1)
            ),
            protein_rank_count = Dict(
                ("P_many", true, UInt8(1)) => Int32(10),
                ("P_one", true, UInt8(1)) => Int32(1)
            ),
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

        @test many.consensus_precursor_rank_support[1] == 1.0f0
        @test one.consensus_precursor_rank_support[1] == 1.0f0
        @test many.consensus_precursor_rank_enrichment[1] > one.consensus_precursor_rank_enrichment[1]
    end

    @testset "Optional Probit Feature Columns Drop Cleanly" begin
        feature_names = Pioneer.protein_probit_feature_names()
        df = DataFrame(pg_score = Float32[0.1, 0.2, 0.3])

        Pioneer.remove_zero_variance_columns!(feature_names, df)

        @test feature_names == [:pg_score]
    end

    @testset "Protein Probit Feature Names Include Consensus Support" begin
        @test Pioneer.protein_probit_feature_names() == [
            :pg_score,
            :peptide_coverage,
            :any_common_peps,
            :consensus_precursor_rank_support,
            :consensus_precursor_rank_enrichment,
            :coverage_miss_surprisal,
            :coverage_deficit_z,
            :top_weight_vs_threshold_z,
            :pg_score_x_coverage_miss_surprisal,
            :pg_score_x_coverage_deficit_z,
            :pg_score_x_top_weight_vs_threshold_z,
            :pg_score_x_consensus_precursor_rank_support,
            :pg_score_x_consensus_precursor_rank_enrichment
        ]

        @test Pioneer.protein_probit_feature_names(include_n_possible_peptides = true) == [
            :pg_score,
            :peptide_coverage,
            :n_possible_peptides,
            :any_common_peps,
            :consensus_precursor_rank_support,
            :consensus_precursor_rank_enrichment,
            :coverage_miss_surprisal,
            :coverage_deficit_z,
            :top_weight_vs_threshold_z,
            :pg_score_x_coverage_miss_surprisal,
            :pg_score_x_coverage_deficit_z,
            :pg_score_x_top_weight_vs_threshold_z,
            :pg_score_x_consensus_precursor_rank_support,
            :pg_score_x_consensus_precursor_rank_enrichment
        ]
    end

end
