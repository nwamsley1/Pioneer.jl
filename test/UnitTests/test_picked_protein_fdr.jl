# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

using Test
using Pioneer

@testset "Picked protein-group FDR" begin
    # pg_id => name. Target group "A;B" (id 1) has a decoy grouped only as "A" (id 2);
    # "C" (id 3) has a same-name decoy that beats it; "D" (id 4) is a target with no decoy;
    # "E" (id 5) is a decoy with no target. Entrapment group 1 reuses accession A.
    names = ["A;B", "A", "C", "D", "E"]
    scores = Dict{Tuple{UInt32, Bool, UInt8}, Float32}(
        (UInt32(1), true, 0x00) => 0.9f0,   # target A;B
        (UInt32(2), false, 0x00) => 0.8f0,  # decoy A (target A;B is present)
        (UInt32(3), true, 0x00) => 0.3f0,   # target C
        (UInt32(3), false, 0x00) => 0.5f0,  # decoy C beats its target
        (UInt32(4), true, 0x00) => 0.7f0,   # target D, no decoy
        (UInt32(5), false, 0x00) => 0.6f0,  # decoy E, no target
        (UInt32(2), true, 0x01) => 0.4f0,   # entrapment target A_trap
    )

    @testset "mode parsing" begin
        withenv("PIONEER_PROTEIN_PICK" => nothing) do
            @test Pioneer.protein_pick_mode() == :off
        end
        withenv("PIONEER_PROTEIN_PICK" => "Accession") do
            @test Pioneer.protein_pick_mode() == :accession
        end
        withenv("PIONEER_PROTEIN_PICK" => "bogus") do
            @test_throws ArgumentError Pioneer.protein_pick_mode()
        end
    end

    @testset "name pairing leaves the A;B / A mismatch unpaired" begin
        q = Pioneer.build_protein_global_qval_dict_picked(scores, :name; group_names = names)
        @test length(q) == length(scores)
        @test q[(UInt32(3), true, 0x00)] == 1f0           # target C loses to its decoy
        @test q[(UInt32(2), false, 0x00)] < 1f0           # decoy A has no same-name partner: kept
        @test q[(UInt32(1), true, 0x00)] < 1f0
    end

    @testset "accession pairing removes the decoy of a present protein" begin
        q = Pioneer.build_protein_global_qval_dict_picked(scores, :accession; group_names = names)
        @test q[(UInt32(2), false, 0x00)] == 1f0          # decoy A loses to target A;B (shares A)
        @test q[(UInt32(3), true, 0x00)] == 1f0           # target C loses to decoy C
        @test q[(UInt32(1), true, 0x00)] < 1f0
        @test q[(UInt32(4), true, 0x00)] < 1f0            # uncontested target kept
        @test q[(UInt32(5), false, 0x00)] < 1f0           # uncontested decoy kept
        # The entrapment group competes only within entrapment group 1 (no decoy there): kept.
        @test q[(UInt32(2), true, 0x01)] < 1f0
        # Fewer decoys survive than with name pairing, so the top target's q can only drop.
        qn = Pioneer.build_protein_global_qval_dict_picked(scores, :accession; group_names = names)
        @test q[(UInt32(1), true, 0x00)] <= qn[(UInt32(1), true, 0x00)]
    end

    @testset "ties go to the decoy" begin
        tie = Dict{Tuple{UInt32, Bool, UInt8}, Float32}(
            (UInt32(1), true, 0x00) => 0.5f0, (UInt32(1), false, 0x00) => 0.5f0)
        for mode in (:name, :accession)
            q = Pioneer.build_protein_global_qval_dict_picked(tie, mode; group_names = ["A"])
            @test q[(UInt32(1), true, 0x00)] == 1f0
        end
    end
end
