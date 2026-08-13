using DataFrames
using Test

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
end
