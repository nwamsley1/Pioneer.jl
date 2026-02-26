using Test
using DataFrames
using Pioneer

@testset "Protein Group Peak-Area Coverage Features" begin
    @testset "Quant Necessary Columns Include Peak Area" begin
        cols_nombr = Pioneer.get_quant_necessary_columns(false)
        cols_mbr = Pioneer.get_quant_necessary_columns(true)

        @test :peak_area in cols_nombr
        @test :peak_area in cols_mbr
    end

    @testset "Group PSMs By Protein Adds Peak-Area Summary" begin
        psms = DataFrame(
            inferred_protein_group = ["P1", "P1", "P1", "P2"],
            target = Bool[true, true, true, false],
            entrap_id = UInt8[1, 1, 1, 1],
            sequence = ["PEP1", "PEP1", "PEP2", "PEP3"],
            use_for_protein_quant = Bool[true, true, true, true],
            missed_cleavage = Int64[0, 1, 0, 0],
            Mox = Int64[0, 0, 0, 0],
            prec_prob = Float32[0.8, 0.9, 0.6, 0.7],
            peak_area = Float32[100.0, 80.0, 50.0, 25.0]
        )

        grouped = Pioneer.group_psms_by_protein(psms)

        @test hasproperty(grouped, :top_pep_peak_area)
        @test hasproperty(grouped, :has_valid_peak_area)

        p1 = grouped[(grouped.protein_name .== "P1") .& grouped.target, :]
        @test nrow(p1) == 1
        @test p1.n_peptides[1] == 2
        @test p1.top_pep_peak_area[1] == 100.0f0
        @test p1.has_valid_peak_area[1] == true
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

    @testset "Peak-Area Calibration Estimates Threshold and Sigma" begin
        psms = DataFrame(
            use_for_protein_quant = Bool[true, true, true, true, false],
            sequence = ["P1A", "P1B", "P2A", "P2B", "P3A"],
            inferred_protein_group = Union{Missing, String}["P1", "P1", "P2", "P2", missing],
            peak_area = Union{Missing, Float32}[100.0f0, 40.0f0, 80.0f0, 10.0f0, 500.0f0]
        )

        model = Pioneer.estimate_peak_area_detection_model(psms)

        @test model.n_unique_peptides == 4
        @test model.log_threshold > 0.0f0
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
            n_possible_peptides = Int64[20, 2, 11, 12, 1],
            top_pep_peak_area = Float32[100.0, 100.0, 10.0, 0.0, 100.0],
            has_valid_peak_area = Bool[true, true, true, false, true]
        )

        (_, op) = Pioneer.add_peak_area_observation_features(calibration)
        out = op(copy(pg_df))

        for col in (:coverage_miss_pval, :coverage_miss_surprisal, :coverage_deficit_z, :top_area_vs_threshold_z)
            @test hasproperty(out, col)
        end

        @test out.coverage_miss_pval[1] < out.coverage_miss_pval[2]
        @test out.coverage_miss_surprisal[1] > out.coverage_miss_surprisal[2]
        @test abs(out.coverage_deficit_z[3]) < 0.25f0
        @test out.top_area_vs_threshold_z[1] > 0.0f0
        @test abs(out.top_area_vs_threshold_z[3]) < 1e-5

        # Invalid peak area -> neutral values
        @test out.coverage_miss_pval[4] == 1.0f0
        @test out.coverage_miss_surprisal[4] == 0.0f0
        @test out.coverage_deficit_z[4] == 0.0f0
        @test out.top_area_vs_threshold_z[4] == 0.0f0

        # N <= 1 -> neutral values
        @test out.coverage_miss_pval[5] == 1.0f0
        @test out.coverage_miss_surprisal[5] == 0.0f0
        @test out.coverage_deficit_z[5] == 0.0f0
        @test out.top_area_vs_threshold_z[5] == 0.0f0
    end

    @testset "Optional Probit Feature Columns Drop Cleanly" begin
        feature_names = [:pg_score, :coverage_miss_surprisal, :coverage_deficit_z, :top_area_vs_threshold_z]
        df = DataFrame(pg_score = Float32[0.1, 0.2, 0.3])

        Pioneer.remove_zero_variance_columns!(feature_names, df)

        @test feature_names == [:pg_score]
    end
end
