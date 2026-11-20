# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

using Test
using Pioneer

@testset "Protein grouping from targets with mirrored decoys" begin
    # Accessions 1 = A, 2 = B. Accession sets: 1 = {A}, 2 = {A, B}, 3 = {B}.
    accession_names = ["A", "B"]
    registry = Pioneer.ProteinGroupRegistry(length(accession_names))
    set_members = [UInt32[1], UInt32[1, 2], UInt32[2]]
    set_group = UInt32[Pioneer.intern_group!(registry, m) for m in set_members]
    # precursor:     1     2     3     4     5     6
    # kind:          T     T     D     D     D     T(entrap)
    sequence_id  = UInt32[10,   11,   20,   21,   22,   30]
    accession_id = UInt32[2,    1,    2,    1,    3,    1]
    is_decoy     = Bool[false, false, true, true, true, false]
    entrap_id    = UInt8[0,    0,    0,    0,    0,    1]
    base_pep_id  = UInt32[100,  101,  100,  555,  999,  100]
    inputs = Pioneer.IntegerInferenceInputs(sequence_id, accession_id, set_group, is_decoy, entrap_id,
                                            registry, accession_names)

    keys = Set(Pioneer._inference_key(inputs, p) for p in 1:6)
    target_keys = filter(k -> !k[3], keys)
    target_result = Pioneer._infer_integer_proteins(inputs, target_keys)
    mirrored = Pioneer._mirror_decoys(inputs, base_pep_id, target_result, target_keys, UInt32[3, 4, 5])

    group_of(seq, target, ent) = mirrored.peptide_to_protein[Pioneer.PeptideKey(seq, target, ent)]
    # Decoy 3 follows its source target (same base peptide) into that target's group.
    @test group_of(UInt32(20), false, 0x00).name == group_of(UInt32(10), true, 0x00).name
    @test !group_of(UInt32(20), false, 0x00).is_target
    # Decoy 4's base peptide has no target, so it takes the target group holding accession set {A}.
    @test group_of(UInt32(21), false, 0x00).name == group_of(UInt32(11), true, 0x00).name
    # Decoy 5 ({B}) has neither: it gets its own accession set's group.
    @test group_of(UInt32(22), false, 0x00).name == set_group[3]
    # Entrapment targets keep their own group (entrap_id 1), even with a shared base peptide.
    @test group_of(UInt32(30), true, 0x01).entrap_id == 0x01
    # Targets are untouched.
    for (seq, protein) in pairs(target_result.peptide_to_protein)
        @test mirrored.peptide_to_protein[seq] == protein
    end

    @testset "mode parsing" begin
        withenv("PIONEER_PROTEIN_GROUPING" => nothing) do
            @test Pioneer.protein_grouping_mode() == :default
        end
        withenv("PIONEER_PROTEIN_GROUPING" => "TARGET") do
            @test Pioneer.protein_grouping_mode() == :target
        end
        withenv("PIONEER_PROTEIN_GROUPING" => "bogus") do
            @test_throws ArgumentError Pioneer.protein_grouping_mode()
        end
    end
end
