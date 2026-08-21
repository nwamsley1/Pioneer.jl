using DataFrames
using Test

import Pioneer
using Pioneer: group_psms_by_protein, run_level_protein_feature_names

function _empty_precursor_consensus_for_mbr_feature_tests()
    return (
        selected_run_votes = Dict{Tuple{String, Bool, UInt8}, Vector{Any}}(),
        consensus_target_run_count = Dict{Tuple{String, Bool, UInt8}, Int32}(),
        cached_consensus_weight_sums = Dict{Tuple{String, Bool, UInt8}, Dict{UInt32, Float64}}(),
        cached_protein_total_vote = Dict{Tuple{String, Bool, UInt8}, Float64}(),
        shape_confidence_scale = 1.0f0
    )
end

@testset "ProteinScoring MBR feature rollup" begin
    @testset "Counterfactual protein controls select one MBR-only target precursor" begin
        psms = DataFrame(
            inferred_protein_group = [
                "P_MBR", "P_MBR", "P_MIXED", "P_MIXED", "DECOY_P",
            ],
            target = Bool[true, true, true, true, false],
            entrap_id = zeros(UInt8, 5),
            use_for_protein_quant = trues(5),
            qval = fill(0.001f0, 5),
            global_qval = fill(0.001f0, 5),
            mbr_recovered = Bool[true, true, true, false, true],
            precursor_idx = UInt32[10, 11, 20, 21, 30],
            prec_prob = fill(0.9f0, 5),
            mbr_counterfactual_decoy_prec_prob =
                Float32[0.2, 0.4, 0.8, NaN, 0.7],
        )

        shadows = Pioneer._mbr_counterfactual_shadow_psms(psms)

        @test nrow(shadows) == 1
        @test shadows.inferred_protein_group == ["P_MBR"]
        @test shadows.precursor_idx == UInt32[11]
        @test shadows.prec_prob == Float32[0.4]
    end

    @testset "Counterfactual protein controls are training-only negatives" begin
        actual = DataFrame(
            protein_name = ["P_SINGLE", "P_MULTI", "DECOY_P"],
            target = Bool[true, true, false],
            mbr_only_protein = Bool[true, true, true],
            mbr_recovered_peptides = Int64[1, 2, 1],
            pg_score = Float32[0.8, 0.9, 0.2],
        )
        shadows = DataFrame(
            protein_name = ["P_SINGLE"],
            target = Bool[true],
            mbr_only_protein = Bool[true],
            mbr_recovered_peptides = Int64[1],
            pg_score = Float32[0.7],
        )

        training = Pioneer.prepare_run_level_protein_training_rows(
            actual,
            shadows,
        )

        @test nrow(training) == 4
        @test training.protein_training_label ==
            Bool[true, true, false, false]
        @test training.protein_training_eligible ==
            Bool[false, true, true, true]
        @test training.protein_shadow_negative ==
            Bool[false, false, false, true]
        @test training.target == Bool[true, true, false, true]
    end

    @testset "MBR recovered support becomes run-level protein features" begin
        psms = DataFrame(
            inferred_protein_group = ["P_MBR", "P_MBR"],
            species = ["YEAST", "YEAST"],
            target = Bool[true, true],
            entrap_id = UInt8[0, 0],
            use_for_protein_quant = Bool[true, true],
            qval = Float32[0.0, 0.0],
            global_qval = Float32[0.001, 0.001],
            precursor_idx = UInt32[10, 11],
            prec_prob = Float32[0.12, 0.15],
            peak_area = Float32[1000.0, 800.0],
            base_pep_id = UInt32[101, 102],
            sequence = ["PEPTIDEA", "PEPTIDEB"],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            mbr_recovered = Bool[true, true],
        )

        protein_groups = group_psms_by_protein(
            psms;
            precursor_consensus = _empty_precursor_consensus_for_mbr_feature_tests(),
            q_value_threshold = 0.01f0
        )

        @test nrow(protein_groups) == 1
        @test protein_groups.mbr_recovered_peptides[1] == 2
        @test protein_groups.n_non_mbr_peptides[1] == 0
        @test !protein_groups.single_non_mbr_peptide[1]
        @test protein_groups.single_non_mbr_prefix_shape[1] == 0.0f0
        @test protein_groups.mbr_only_protein[1]
        @test protein_groups.any_common_peps[1]
    end

    @testset "Mixed support tracks retained MBR evidence without changing pg_score input" begin
        psms = DataFrame(
            inferred_protein_group = ["P_MIXED", "P_MIXED"],
            species = ["YEAST", "YEAST"],
            target = Bool[true, true],
            entrap_id = UInt8[0, 0],
            use_for_protein_quant = Bool[true, true],
            qval = Float32[0.0, 0.001],
            global_qval = Float32[0.001, 0.001],
            precursor_idx = UInt32[20, 21],
            prec_prob = Float32[0.20, 0.80],
            peak_area = Float32[500.0, 1500.0],
            base_pep_id = UInt32[201, 202],
            sequence = ["PEPTIDEC", "PEPTIDED"],
            missed_cleavage = Int64[0, 0],
            Mox = Int64[0, 0],
            mbr_recovered = Bool[true, false],
        )

        protein_groups = group_psms_by_protein(
            psms;
            precursor_consensus = _empty_precursor_consensus_for_mbr_feature_tests(),
            q_value_threshold = 0.01f0
        )

        @test nrow(protein_groups) == 1
        @test protein_groups.mbr_recovered_peptides[1] == 1
        @test protein_groups.n_non_mbr_peptides[1] == 1
        @test protein_groups.single_non_mbr_peptide[1]
        @test protein_groups.single_non_mbr_prefix_shape[1] ==
            protein_groups.precursor_consensus_prefix_shape[1]
        @test !protein_groups.mbr_only_protein[1]
        @test protein_groups.pg_score[1] > 0.0f0

        non_mbr_groups = group_psms_by_protein(
            select(psms, Not(:mbr_recovered));
            precursor_consensus =
                _empty_precursor_consensus_for_mbr_feature_tests(),
            q_value_threshold = 0.01f0,
        )
        @test non_mbr_groups.n_non_mbr_peptides[1] == 2
        @test !non_mbr_groups.single_non_mbr_peptide[1]
        @test non_mbr_groups.single_non_mbr_prefix_shape[1] == 0.0f0
        @test non_mbr_groups.mbr_recovered_peptides[1] == 0
    end

    @testset "Non-MBR peptide count deduplicates sequences independently" begin
        psms = DataFrame(
            inferred_protein_group = fill("P_OVERLAP", 3),
            species = fill("YEAST", 3),
            target = trues(3),
            entrap_id = zeros(UInt8, 3),
            use_for_protein_quant = trues(3),
            qval = fill(0.001f0, 3),
            global_qval = fill(0.001f0, 3),
            precursor_idx = UInt32[30, 31, 32],
            prec_prob = Float32[0.20, 0.30, 0.40],
            peak_area = Float32[500.0, 750.0, 1000.0],
            base_pep_id = UInt32[301, 301, 302],
            sequence = ["PEPTIDEE", "PEPTIDEE", "PEPTIDEF"],
            missed_cleavage = zeros(Int64, 3),
            Mox = zeros(Int64, 3),
            mbr_recovered = Bool[true, false, true],
        )

        protein_groups = group_psms_by_protein(
            psms;
            precursor_consensus =
                _empty_precursor_consensus_for_mbr_feature_tests(),
            q_value_threshold = 0.01f0,
        )

        @test protein_groups.n_peptides[1] == 2
        @test protein_groups.mbr_recovered_peptides[1] == 2
        @test protein_groups.n_non_mbr_peptides[1] == 1
        @test protein_groups.single_non_mbr_peptide[1]
        @test protein_groups.single_non_mbr_prefix_shape[1] ==
            protein_groups.precursor_consensus_prefix_shape[1]
        @test !protein_groups.mbr_only_protein[1]
    end

    @testset "Run-level model includes the lightweight MBR features" begin
        features = run_level_protein_feature_names()
        @test :single_non_mbr_peptide in features
        @test :single_non_mbr_prefix_shape in features
        @test :mbr_recovered_peptides in features
        @test :mbr_only_protein in features
        @test !in(:mbr_recovered_precursors, features)
        @test !in(:mbr_recovered_precursor_fraction, features)
        @test !in(:mbr_best_transfer_confidence, features)
        @test !in(:mbr_sum_transfer_evidence, features)
        @test !in(:mbr_best_pair_prob, features)
    end

    @testset "common peptide evidence uses specificity, cleavage, and variable modifications" begin
        psms = DataFrame(
            inferred_protein_group =
                ["P_SEMI", "P_FULL", "P_MISSED", "P_VARIABLE"],
            species = fill("YEAST", 4),
            target = trues(4),
            entrap_id = zeros(UInt8, 4),
            use_for_protein_quant = trues(4),
            qval = zeros(Float32, 4),
            global_qval = fill(0.001f0, 4),
            precursor_idx = UInt32[40, 41, 42, 43],
            prec_prob = fill(0.9f0, 4),
            peak_area = fill(1000.0f0, 4),
            base_pep_id = UInt32[401, 402, 403, 404],
            sequence = ["SEMIPEP", "FULLPEP", "MISSEDPEP", "PHOSPEP"],
            missed_cleavage = Int64[0, 0, 1, 0],
            num_enzymatic_termini = UInt8[1, 2, 2, 2],
            num_variable_modifications = UInt8[0, 0, 0, 1],
            Mox = zeros(Int64, 4),
            mbr_recovered = falses(4),
        )

        protein_groups = group_psms_by_protein(
            psms;
            precursor_consensus =
                _empty_precursor_consensus_for_mbr_feature_tests(),
            q_value_threshold = 0.01f0,
        )
        common_by_protein = Dict(
            String(row.protein_name) => Bool(row.any_common_peps)
            for row in eachrow(protein_groups)
        )

        @test !common_by_protein["P_SEMI"]
        @test common_by_protein["P_FULL"]
        @test !common_by_protein["P_MISSED"]
        @test !common_by_protein["P_VARIABLE"]

        full_row = only(eachrow(protein_groups[protein_groups.protein_name .== "P_FULL", :]))
        @test full_row.n_common_peptides == 1
        @test full_row.common_peptide_list == "FULLPEP"
        for protein_name in ("P_SEMI", "P_MISSED", "P_VARIABLE")
            row = only(eachrow(
                protein_groups[protein_groups.protein_name .== protein_name, :]
            ))
            @test row.n_common_peptides == 0
            @test isempty(row.common_peptide_list)
        end
    end

    @testset "PSM models include enzymatic specificity" begin
        @test :num_enzymatic_termini in Pioneer.PRESCORE_FEATURES
        @test :num_enzymatic_termini in Pioneer.ADVANCED_FEATURE_SET
    end

    @testset "Run-level model separates common and all-peptide coverage" begin
        features = run_level_protein_feature_names()
        @test :peptide_coverage_logit in features
        @test :all_peptide_coverage_logit in features
        @test :all_peptide_coverage_logit in
            Pioneer.PROTEIN_MONOTONE_INCREASING_FEATURES
    end
end
