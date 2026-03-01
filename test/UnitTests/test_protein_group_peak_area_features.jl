using Test
using DataFrames
using Statistics
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
            weight = Union{Missing, Float32}[60.0f0, 40.0f0, 40.0f0, 80.0f0, 10.0f0, 500.0f0]
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
            peptide_coverage = Float32[0.05, 0.5, 0.55, 0.08, 1.0],
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
            :pg_score_x_peptide_coverage,
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

        @test out.pg_score_x_peptide_coverage[1] == out.pg_score[1] * out.peptide_coverage[1]
        @test out.pg_score_x_coverage_miss_surprisal[1] == out.pg_score[1] * out.coverage_miss_surprisal[1]
        @test out.pg_score_x_coverage_deficit_z[1] == out.pg_score[1] * out.coverage_deficit_z[1]
        @test out.pg_score_x_top_weight_vs_threshold_z[1] == out.pg_score[1] * out.top_weight_vs_threshold_z[1]
        @test out.pg_score_x_peptide_coverage[3] == out.pg_score[3] * out.peptide_coverage[3]
        @test out.pg_score_x_coverage_miss_surprisal[3] == out.pg_score[3] * out.coverage_miss_surprisal[3]
        @test out.pg_score_x_coverage_deficit_z[3] == out.pg_score[3] * out.coverage_deficit_z[3]
        @test out.pg_score_x_top_weight_vs_threshold_z[3] == out.pg_score[3] * out.top_weight_vs_threshold_z[3]
    end

    @testset "Optional Probit Feature Columns Drop Cleanly" begin
        feature_names = [
            :pg_score,
            :coverage_miss_surprisal,
            :coverage_deficit_z,
            :top_weight_vs_threshold_z,
            :pg_score_x_peptide_coverage,
            :pg_score_x_coverage_miss_surprisal,
            :pg_score_x_coverage_deficit_z,
            :pg_score_x_top_weight_vs_threshold_z
        ]
        df = DataFrame(pg_score = Float32[0.1, 0.2, 0.3])

        Pioneer.remove_zero_variance_columns!(feature_names, df)

        @test feature_names == [:pg_score]
    end

    @testset "Protein Feature QC Plots Are Written" begin
        df = DataFrame(
            target = Bool[true, true, false, false],
            pg_score = Float32[0.99, 0.85, 0.20, 0.05],
            top_pep_weight = Float32[100.0, 20.0, 15.0, 5.0],
            coverage_miss_pval = Float32[1e-3, 0.2, 0.6, 1.0],
            coverage_miss_surprisal = Float32[3.0, 0.7, 0.2, 0.0],
            coverage_deficit_z = Float32[-4.5, -1.0, 0.5, 0.0],
            top_weight_vs_threshold_z = Float32[3.0, 1.0, 0.5, 0.0]
        )

        qc_dir = mktempdir()
        try
            paths = Pioneer.generate_protein_feature_qc_plots(
                df,
                qc_dir;
                prefix = "protein_weight_feature_qc_test"
            )

            @test isfile(paths.combined_pdf)
            @test !isempty(paths.individual_pdfs)
            @test all(isfile, paths.individual_pdfs)
        finally
            rm(qc_dir, recursive = true, force = true)
        end
    end
end
