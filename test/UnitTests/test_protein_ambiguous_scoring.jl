using DataFrames
using Dictionaries
using Test

using Pioneer: AMBIGUOUS_PROTEIN_SCORE_PSEUDOCOUNT,
    InferenceResult,
    PeptideKey,
    ProteinKey,
    ProteinPeptideOpportunityCounts,
    _register_protein_ambiguities!,
    add_ambiguous_pg_score!,
    add_protein_features,
    add_protein_ambiguity_id,
    count_protein_peptide_opportunities,
    protein_probit_feature_names

function _ambiguous_scoring_psms(; qval::Float32 = 0.001f0)
    return DataFrame(
        protein_ambiguity_id = UInt32[1],
        target = Bool[true],
        entrap_id = UInt8[0],
        precursor_idx = UInt32[101],
        prec_prob = Float32[0.8],
        peak_area = Float32[1000.0],
        base_pep_id = UInt32[11],
        sequence = ["SHARED"],
        qval = Float32[qval],
        global_qval = Float32[0.001]
    )
end

@testset "Ambiguous peptide protein scoring" begin
    protein_a = ProteinKey("A", true, UInt8(0))
    protein_b = ProteinKey("B", true, UInt8(0))
    candidates = Dict{UInt32, Vector{ProteinKey}}(1 => [protein_a, protein_b])

    @testset "raw unique scores determine conserved allocation" begin
        protein_groups = DataFrame(
            protein_name = ["A", "B"],
            target = Bool[true, true],
            entrap_id = UInt8[0, 0],
            pg_score = Float32[2.0, 8.0]
        )
        original_pg_score = copy(protein_groups.pg_score)

        add_ambiguous_pg_score!(
            protein_groups,
            _ambiguous_scoring_psms(),
            candidates
        )

        peptide_score = Float32(-log((1.0f0 - 0.8f0) + 0.001f0))
        support_a = 2.0f0 + AMBIGUOUS_PROTEIN_SCORE_PSEUDOCOUNT
        support_b = 8.0f0 + AMBIGUOUS_PROTEIN_SCORE_PSEUDOCOUNT
        total_support = support_a + support_b

        @test protein_groups.pg_score == original_pg_score
        @test protein_groups.ambiguous_pg_score[1] ≈
            peptide_score * support_a / total_support rtol = 1f-6
        @test protein_groups.ambiguous_pg_score[2] ≈
            peptide_score * support_b / total_support rtol = 1f-6
        @test sum(protein_groups.ambiguous_pg_score) ≈ peptide_score rtol = 1f-6
    end

    @testset "only eligible groups receive ambiguous evidence" begin
        protein_groups = DataFrame(
            protein_name = ["A"],
            target = Bool[true],
            entrap_id = UInt8[0],
            pg_score = Float32[2.0]
        )

        add_ambiguous_pg_score!(
            protein_groups,
            _ambiguous_scoring_psms(),
            candidates
        )

        peptide_score = Float32(-log((1.0f0 - 0.8f0) + 0.001f0))
        @test protein_groups.ambiguous_pg_score == Float32[peptide_score]
    end

    @testset "precursors are deduplicated before allocation" begin
        protein_groups = DataFrame(
            protein_name = ["A", "B"],
            target = Bool[true, true],
            entrap_id = UInt8[0, 0],
            pg_score = Float32[2.0, 8.0],
            n_peptides = Int64[1, 1]
        )
        psms = vcat(_ambiguous_scoring_psms(), _ambiguous_scoring_psms())
        psms.prec_prob = Float32[0.2, 0.8]

        add_ambiguous_pg_score!(protein_groups, psms, candidates)
        protein_peptide_opportunities = Dict(
            protein_a => ProteinPeptideOpportunityCounts(2, 2),
            protein_b => ProteinPeptideOpportunityCounts(4, 1)
        )
        add_protein_features(protein_peptide_opportunities).second(protein_groups)

        expected_score = Float32(-log((1.0f0 - 0.8f0) + 0.001f0))
        @test sum(protein_groups.ambiguous_pg_score) ≈ expected_score rtol = 1f-6
        @test protein_groups.n_possible_unique_peptides == Int64[2, 4]
        @test protein_groups.n_possible_shared_peptides == Int64[2, 1]
        @test protein_groups.peptide_coverage == Float32[0.5, 0.25]
        @test protein_groups.shared_peptide_coverage == Float32[0.5, 1.0]
        @test !hasproperty(protein_groups, :_ambiguous_peptide_count)
    end

    @testset "theoretical opportunities are classified against final groups" begin
        group_ab = ProteinKey("A;B", true, UInt8(0))
        group_c = ProteinKey("C", true, UInt8(0))
        decoy_d = ProteinKey("D", false, UInt8(0))
        entrap_e = ProteinKey("E", true, UInt8(1))
        final_groups = Set([group_ab, group_c, decoy_d, entrap_e])

        opportunities = count_protein_peptide_opportunities(
            [
                "A",
                "A;B",
                "C",
                "C",
                "B;C",
                "A",
                "C",
                "X",
                "D",
                "E",
            ],
            [
                "A_ONLY",
                "INTRA_GROUP",
                "C_ONLY",
                "C_ONLY",
                "BETWEEN_GROUPS",
                "SPLIT_SHARED",
                "SPLIT_SHARED",
                "IGNORED",
                "DECOY_ONLY",
                "ENTRAP_ONLY",
            ],
            Bool[
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                false,
                true,
                false,
            ],
            UInt8[0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            final_groups
        )

        @test opportunities[group_ab] == ProteinPeptideOpportunityCounts(2, 2)
        @test opportunities[group_c] == ProteinPeptideOpportunityCounts(1, 2)
        @test opportunities[decoy_d] == ProteinPeptideOpportunityCounts(1, 0)
        @test opportunities[entrap_e] == ProteinPeptideOpportunityCounts(1, 0)
    end

    @testset "zero shared opportunity denominator gives zero coverage" begin
        protein_groups = DataFrame(
            protein_name = ["A"],
            target = Bool[true],
            entrap_id = UInt8[0],
            n_peptides = Int64[1],
            _ambiguous_peptide_count = Int64[0]
        )
        add_protein_features(Dict(
            protein_a => ProteinPeptideOpportunityCounts(1, 0)
        )).second(protein_groups)

        @test protein_groups.peptide_coverage == Float32[1.0]
        @test protein_groups.shared_peptide_coverage == Float32[0.0]
    end

    @testset "allocation cannot cross target/decoy populations" begin
        decoy_b = ProteinKey("B", false, UInt8(0))
        mixed_candidates = Dict{UInt32, Vector{ProteinKey}}(
            1 => [protein_a, decoy_b]
        )
        protein_groups = DataFrame(
            protein_name = ["A", "B"],
            target = Bool[true, false],
            entrap_id = UInt8[0, 0],
            pg_score = Float32[2.0, 8.0]
        )

        add_ambiguous_pg_score!(
            protein_groups,
            _ambiguous_scoring_psms(),
            mixed_candidates
        )

        @test protein_groups.ambiguous_pg_score[1] > 0.0f0
        @test protein_groups.ambiguous_pg_score[2] == 0.0f0
    end

    @testset "confidence thresholds and absent ambiguity preserve zero contribution" begin
        protein_groups = DataFrame(
            protein_name = ["A", "B"],
            target = Bool[true, true],
            entrap_id = UInt8[0, 0],
            pg_score = Float32[2.0, 8.0]
        )

        add_ambiguous_pg_score!(
            protein_groups,
            _ambiguous_scoring_psms(qval = 0.02f0),
            candidates
        )
        @test protein_groups.ambiguous_pg_score == zeros(Float32, 2)
        @test protein_groups._ambiguous_peptide_count == zeros(Int64, 2)

        psms_without_ambiguity = select(
            _ambiguous_scoring_psms(),
            Not(:protein_ambiguity_id)
        )
        add_ambiguous_pg_score!(protein_groups, psms_without_ambiguity, candidates)
        @test protein_groups.ambiguous_pg_score == zeros(Float32, 2)
        @test protein_groups._ambiguous_peptide_count == zeros(Int64, 2)
    end

    @testset "in-memory ambiguity registration" begin
        peptide = PeptideKey("SHARED", true, UInt8(0))
        ambiguous = Dictionary{PeptideKey, Vector{ProteinKey}}()
        insert!(ambiguous, peptide, [protein_b, protein_a])
        result = InferenceResult(
            Dictionary{PeptideKey, ProteinKey}(),
            ambiguous
        )
        candidates_by_id = Dict{UInt32, Vector{ProteinKey}}()
        peptide_to_id, final_id = _register_protein_ambiguities!(
            candidates_by_id,
            result,
            zero(UInt32)
        )

        @test final_id == UInt32(1)
        @test peptide_to_id[peptide] == UInt32(1)
        @test candidates_by_id[UInt32(1)] == [protein_a, protein_b]

        annotation = DataFrame(
            sequence = ["SHARED", "UNASSIGNED"],
            is_decoy = Bool[false, false],
            entrap_id = UInt8[0, 0]
        )
        add_protein_ambiguity_id(peptide_to_id).second(annotation)
        @test annotation.protein_ambiguity_id == UInt32[1, 0]

    end

    @testset "probit uses only the requested ambiguity features" begin
        features = protein_probit_feature_names()
        @test :ambiguous_pg_score in features
        @test :shared_peptide_coverage in features
        @test !(:ambiguous_peptide_coverage in features)
        @test !(:augmented_pg_score in features)
        @test !(:effective_ambiguous_peptides in features)
        @test !(:ambiguous_score_fraction in features)
    end
end
