# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

using Test
using Dictionaries

# Include only the protein inference functions
# (types are already loaded by importScripts.jl in runtests.jl)
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "utils", "proteinInference.jl"))

@testset "Protein Inference Typed Implementation" begin
    
    @testset "Input Validation" begin
        # Test mismatched vector lengths
        proteins = [ProteinKey("A", true, UInt8(1))]
        peptides = [PeptideKey("pep1", true, UInt8(1)), PeptideKey("pep2", true, UInt8(1))]
        
        @test_throws ArgumentError infer_proteins(proteins, peptides)
        
        # Test empty inputs
        empty_proteins = ProteinKey[]
        empty_peptides = PeptideKey[]
        result = infer_proteins(empty_proteins, empty_peptides)

        @test length(result.peptide_to_protein) == 0
    end

    @testset "Case A: Distinct Proteins" begin
        # Protein A has peptides 1, 2
        # Protein B has peptides 3, 4
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # Check peptide-to-protein mappings
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A"
        @test result.peptide_to_protein[peptides[3]].name == "B"
        @test result.peptide_to_protein[peptides[4]].name == "B"

        # All peptides are unique and present (usable for quantification)
        @test length(result.peptide_to_protein) == 4

        # Check target status and entrapment IDs are preserved
        for (pep_key, prot_key) in pairs(result.peptide_to_protein)
            @test prot_key.is_target == pep_key.is_target
            @test prot_key.entrap_id == pep_key.entrap_id
        end
    end
    
    @testset "Case B: Differentiable Proteins" begin
        # Protein A has peptides 1
        # Proteins A and B share peptides 2, 3
        # Protein B has peptide 4
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]

        result = infer_proteins(proteins, peptides)

        # Only unique peptides should be in results (shared peptides excluded)
        @test length(result.peptide_to_protein) == 2

        # Check unique peptide assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[4]].name == "B"

        # Shared peptides should not be in results (not usable for quantification)
        @test !haskey(result.peptide_to_protein, peptides[2])  # Shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[3])  # Shared, excluded
    end
    
    @testset "Case C: Indistinguishable Proteins" begin
        # Proteins A and B have identical peptide sets
        proteins = [
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # All peptides should be assigned to the combined protein group
        @test length(result.peptide_to_protein) == 4
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A;B"
        end
    end
    
    @testset "Case D: Subset Proteins" begin
        # Protein A can explain all peptides (subset case)
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # All peptides should be assigned to A (minimal explanation)
        @test length(result.peptide_to_protein) == 4
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A"
        end
    end
    
    @testset "Case E: Subsumable Proteins" begin
        # Need both A and C to explain all peptides
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("B;C", true, UInt8(1)),
            ProteinKey("C", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # Check that minimal set A and C is chosen
        # All peptides get assigned to necessary proteins with maximum parsimony
        @test length(result.peptide_to_protein) == 4
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A"  # A can explain pep2
        @test result.peptide_to_protein[peptides[3]].name == "C"  # C can explain pep3
        @test result.peptide_to_protein[peptides[4]].name == "C"

        # All peptides are unique to their assigned proteins (usable for quantification)
    end
    
    @testset "Case F: Shared Peptides Only" begin
        # No protein has unique peptides
        proteins = [
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B;C", true, UInt8(1)),
            ProteinKey("A;B;C", true, UInt8(1)),
            ProteinKey("A;C", true, UInt8(1))
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # All peptides should be assigned to the combined group
        @test length(result.peptide_to_protein) == 4
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A;B;C"
        end
    end
    
    @testset "Entrapment Group Handling" begin
        # Test that different entrapment groups are handled separately
        proteins = [
            ProteinKey("A", true, UInt8(1)),   # Entrapment group 1
            ProteinKey("A", true, UInt8(2))    # Entrapment group 2
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep1", true, UInt8(2))  # Same sequence, different entrapment
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # Each peptide should be assigned to protein in matching entrapment group
        @test length(result.peptide_to_protein) == 2
        @test result.peptide_to_protein[peptides[1]].entrap_id == UInt8(1)
        @test result.peptide_to_protein[peptides[2]].entrap_id == UInt8(2)
    end
    
    @testset "Target/Decoy Handling" begin
        # Test that target and decoy proteins are handled separately
        proteins = [
            ProteinKey("A", true, UInt8(1)),    # Target
            ProteinKey("B", false, UInt8(1))    # Decoy
        ]
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),   # Target peptide
            PeptideKey("pep2", false, UInt8(1))   # Decoy peptide
        ]
        
        result = infer_proteins(proteins, peptides)
        
        # Check target/decoy status is preserved
        @test length(result.peptide_to_protein) == 2
        @test result.peptide_to_protein[peptides[1]].is_target == true
        @test result.peptide_to_protein[peptides[2]].is_target == false
    end
    
    @testset "Complex Case H: Multiple Scenarios" begin
        # Combination of all previous test cases
        proteins = [
            # Case A: Distinct proteins (A, B)
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1)),
            
            # Case B: Differentiable proteins (C, D)
            ProteinKey("C", true, UInt8(1)),
            ProteinKey("C;D", true, UInt8(1)),
            ProteinKey("C;D", true, UInt8(1)),
            ProteinKey("D", true, UInt8(1)),
            
            # Case C: Indistinguishable proteins (E, F)
            ProteinKey("E;F", true, UInt8(1)),
            ProteinKey("E;F", true, UInt8(1)),
            ProteinKey("E;F", true, UInt8(1)),
            ProteinKey("E;F", true, UInt8(1))
        ]
        
        peptides = [
            PeptideKey("pep1", true, UInt8(1)),
            PeptideKey("pep2", true, UInt8(1)),
            PeptideKey("pep3", true, UInt8(1)),
            PeptideKey("pep4", true, UInt8(1)),
            PeptideKey("pep5", true, UInt8(1)),
            PeptideKey("pep6", true, UInt8(1)),
            PeptideKey("pep7", true, UInt8(1)),
            PeptideKey("pep8", true, UInt8(1)),
            PeptideKey("pep9", true, UInt8(1)),
            PeptideKey("pep10", true, UInt8(1)),
            PeptideKey("pep11", true, UInt8(1)),
            PeptideKey("pep12", true, UInt8(1))
        ]
        
        result = infer_proteins(proteins, peptides)

        # Should have 10 peptides (12 total - 2 shared excluded: pep6, pep7)
        @test length(result.peptide_to_protein) == 10

        # Case A results (distinct proteins - all unique)
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A"
        @test result.peptide_to_protein[peptides[3]].name == "B"
        @test result.peptide_to_protein[peptides[4]].name == "B"

        # Case B-like results (differentiable proteins C, D)
        @test result.peptide_to_protein[peptides[5]].name == "C"
        @test result.peptide_to_protein[peptides[8]].name == "D"

        # Shared peptides should not be in results (not usable for quantification)
        @test !haskey(result.peptide_to_protein, peptides[6])  # Shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[7])  # Shared, excluded

        # Case C results (indistinguishable proteins)
        @test result.peptide_to_protein[peptides[9]].name == "E;F"
        @test result.peptide_to_protein[peptides[10]].name == "E;F"
        @test result.peptide_to_protein[peptides[11]].name == "E;F"
        @test result.peptide_to_protein[peptides[12]].name == "E;F"
    end

    @testset "Case I: Complex Component with Merge-First" begin
        # Test case demonstrating merge-first fix for indistinguishable proteins
        # This validates that B;C and E;F are properly merged during greedy set cover
        proteins = [
            ProteinKey("A", true, UInt8(1)),          # pep1
            ProteinKey("A;B", true, UInt8(1)),        # pep2
            ProteinKey("B;C", true, UInt8(1)),        # pep3
            ProteinKey("C;D", true, UInt8(1)),        # pep4
            ProteinKey("D", true, UInt8(1)),          # pep5
            ProteinKey("D;E", true, UInt8(1)),        # pep6
            ProteinKey("F;E", true, UInt8(1)),        # pep7
            ProteinKey("A;F", true, UInt8(1))         # pep8
        ]

        peptides = [
            PeptideKey("1", true, UInt8(1)),
            PeptideKey("2", true, UInt8(1)),
            PeptideKey("3", true, UInt8(1)),
            PeptideKey("4", true, UInt8(1)),
            PeptideKey("5", true, UInt8(1)),
            PeptideKey("6", true, UInt8(1)),
            PeptideKey("7", true, UInt8(1)),
            PeptideKey("8", true, UInt8(1))
        ]

        result = infer_proteins(proteins, peptides)

        # Only unique peptides should be in results (shared peptides excluded)
        @test length(result.peptide_to_protein) == 4

        # Unique peptide assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"      # pep1 unique to A
        @test result.peptide_to_protein[peptides[3]].name == "B;C"    # pep3 unique to B;C merged group
        @test result.peptide_to_protein[peptides[5]].name == "D"      # pep5 unique to D
        @test result.peptide_to_protein[peptides[7]].name == "E;F"    # pep7 unique to E;F merged group

        # Shared peptides are not in results (not usable for quantification)
        @test !haskey(result.peptide_to_protein, peptides[2])  # pep2 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[4])  # pep4 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[6])  # pep6 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[8])  # pep8 shared, excluded
    end

    @testset "Case J: Merge-First with Original Bug Case" begin
        # This is the original failing test case that exposed the merge-first bug
        # Protein-Peptide mapping:
        #   pep1 → {A}
        #   pep2 → {A, B, C}
        #   pep3 → {B, C}
        #   pep4 → {B, C, D}
        #   pep5 → {D}
        # Expected: A selected (unique pep1), D selected (unique pep5), B and C merged to B;C
        # Only unique peptides (pep1, pep3, pep5) should be in results
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A;B;C", true, UInt8(1)),
            ProteinKey("B;C", true, UInt8(1)),
            ProteinKey("B;C;D", true, UInt8(1)),
            ProteinKey("D", true, UInt8(1)),
        ]

        peptides = [
            PeptideKey("1", true, UInt8(1)),
            PeptideKey("2", true, UInt8(1)),
            PeptideKey("3", true, UInt8(1)),
            PeptideKey("4", true, UInt8(1)),
            PeptideKey("5", true, UInt8(1)),
        ]

        result = infer_proteins(proteins, peptides)

        # Only unique peptides should be in results (shared peptides excluded)
        @test length(result.peptide_to_protein) == 3

        # Unique peptide assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"      # pep1 unique to A
        @test result.peptide_to_protein[peptides[3]].name == "B;C"    # pep3 unique to B;C merged group
        @test result.peptide_to_protein[peptides[5]].name == "D"      # pep5 unique to D

        # Shared peptides are not in results (not usable for quantification)
        @test !haskey(result.peptide_to_protein, peptides[2])  # pep2 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[4])  # pep4 shared, excluded
    end

    @testset "Case K: Shared Peptides with Duplicate Protein Groups" begin
        # Test case with duplicate protein groups in input and shared peptides
        # Protein-Peptide mapping:
        #   pep1 → {A}
        #   pep2 → {A, B}
        #   pep3 → {A, C}
        #   pep4 → {B}
        #   pep5 → {A, B}
        #   pep6 → {A, B}
        #   pep7 → {A, B}
        # Expected: A selected (unique pep1), B selected (unique pep4)
        # Shared peptides: pep2, pep5, pep6, pep7 → A;B
        proteins = [
            ProteinKey("A", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;C", true, UInt8(1)),
            ProteinKey("B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1)),
            ProteinKey("A;B", true, UInt8(1))
        ]

        peptides = [
            PeptideKey("1", true, UInt8(1)),
            PeptideKey("2", true, UInt8(1)),
            PeptideKey("3", true, UInt8(1)),
            PeptideKey("4", true, UInt8(1)),
            PeptideKey("5", true, UInt8(1)),
            PeptideKey("6", true, UInt8(1)),
            PeptideKey("7", true, UInt8(1)),
        ]

        result = infer_proteins(proteins, peptides)

        # Only unique peptides should be in results (shared peptides excluded)
        @test length(result.peptide_to_protein) == 3

        # Unique peptide assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"      # pep1 unique to A
        @test result.peptide_to_protein[peptides[3]].name == "A"      # pep3 unique to A (C not selected)
        @test result.peptide_to_protein[peptides[4]].name == "B"      # pep4 unique to B

        # Shared peptides are not in results (not usable for quantification)
        @test !haskey(result.peptide_to_protein, peptides[2])  # pep2 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[5])  # pep5 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[6])  # pep6 shared, excluded
        @test !haskey(result.peptide_to_protein, peptides[7])  # pep7 shared, excluded
    end

    @testset "Real World Complex Example" begin 
        
        # Data from failing_component_2025-10-13T12-00-28.902.jl
        # Simplified protein/peptide names for readability:
        # Original → Simple:
        #   O60814→A, P33778→B, P02293→O, Q6DRA6→Q, Q93079→C, Q99880→D, P58876→E
        #   P23527→L, Q99879→F, Q5QNW6→G, P62807→H, P02294→P, P06899→M, Q99877→I
        #   Q96A08→S, Q16778→N, Q6DN03→R, Q8N257→J, P57053→K
        #   ESYSVYVYK→p1, QVHPDTGISSK→p2, KESYSVYVYK→p3, AMGIMNSFVNDIFER→p4, EIQTAVR→p5
        #   ESYSIYVYK→p6, KESYSIYVYK→p7, ETYSSYIYK→p8, SMSILNSFVNDIFER→p9
        #   ESYSIYIYK→p10, AMSIMNSFVTDIFER→p11

        #=
        Proteins: A-S (19 proteins organized into 6 groups)
        - Group 1: A,C,D,E,F,G,H,I,K → {p1,p2,p3,p4,p5}
        - Group 2: B,L,M,N → {p6,p2,p7,p4,p5}
        - Group 3: O,P (histones) → {p8,p9,p5}
        - Group 4: Q,R → {p6,p7}
        - Group 5: S → {p2,p10,p11,p5}
        - Group 6: J → {p6,p2,p7,p4} (no p5)
        =#

        proteins_realworld = [
            ProteinKey("A", true, UInt8(0)),  # A -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("A", true, UInt8(0)),
            ProteinKey("A", true, UInt8(0)),
            ProteinKey("A", true, UInt8(0)),
            ProteinKey("A", true, UInt8(0)),
            ProteinKey("B", true, UInt8(0)),  # B -> {p6,p2,p7,p4,p5} (Group 2)
            ProteinKey("B", true, UInt8(0)),
            ProteinKey("B", true, UInt8(0)),
            ProteinKey("B", true, UInt8(0)),
            ProteinKey("B", true, UInt8(0)),
            ProteinKey("O", true, UInt8(0)),  # O -> {p8,p9,p5} (Histone, Group 3)
            ProteinKey("O", true, UInt8(0)),
            ProteinKey("O", true, UInt8(0)),
            ProteinKey("Q", true, UInt8(0)),  # Q -> {p6,p7} (Group 4)
            ProteinKey("Q", true, UInt8(0)),
            ProteinKey("C", true, UInt8(0)),  # C -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("C", true, UInt8(0)),
            ProteinKey("C", true, UInt8(0)),
            ProteinKey("C", true, UInt8(0)),
            ProteinKey("C", true, UInt8(0)),
            ProteinKey("D", true, UInt8(0)),  # D -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("D", true, UInt8(0)),
            ProteinKey("D", true, UInt8(0)),
            ProteinKey("D", true, UInt8(0)),
            ProteinKey("D", true, UInt8(0)),
            ProteinKey("E", true, UInt8(0)),  # E -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("E", true, UInt8(0)),
            ProteinKey("E", true, UInt8(0)),
            ProteinKey("E", true, UInt8(0)),
            ProteinKey("E", true, UInt8(0)),
            ProteinKey("L", true, UInt8(0)),  # L -> {p6,p2,p7,p4,p5} (Group 2)
            ProteinKey("L", true, UInt8(0)),
            ProteinKey("L", true, UInt8(0)),
            ProteinKey("L", true, UInt8(0)),
            ProteinKey("L", true, UInt8(0)),
            ProteinKey("F", true, UInt8(0)),  # F -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("F", true, UInt8(0)),
            ProteinKey("F", true, UInt8(0)),
            ProteinKey("F", true, UInt8(0)),
            ProteinKey("F", true, UInt8(0)),
            ProteinKey("G", true, UInt8(0)),  # G -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("G", true, UInt8(0)),
            ProteinKey("G", true, UInt8(0)),
            ProteinKey("G", true, UInt8(0)),
            ProteinKey("G", true, UInt8(0)),
            ProteinKey("H", true, UInt8(0)),  # H -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("H", true, UInt8(0)),
            ProteinKey("H", true, UInt8(0)),
            ProteinKey("H", true, UInt8(0)),
            ProteinKey("H", true, UInt8(0)),
            ProteinKey("P", true, UInt8(0)),  # P -> {p8,p9,p5} (Histone, Group 3)
            ProteinKey("P", true, UInt8(0)),
            ProteinKey("P", true, UInt8(0)),
            ProteinKey("M", true, UInt8(0)),  # M -> {p6,p2,p7,p4,p5} (Group 2)
            ProteinKey("M", true, UInt8(0)),
            ProteinKey("M", true, UInt8(0)),
            ProteinKey("M", true, UInt8(0)),
            ProteinKey("M", true, UInt8(0)),
            ProteinKey("I", true, UInt8(0)),  # I -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("I", true, UInt8(0)),
            ProteinKey("I", true, UInt8(0)),
            ProteinKey("I", true, UInt8(0)),
            ProteinKey("I", true, UInt8(0)),
            ProteinKey("S", true, UInt8(0)),  # S -> {p2,p10,p11,p5} (Group 5)
            ProteinKey("S", true, UInt8(0)),
            ProteinKey("S", true, UInt8(0)),
            ProteinKey("S", true, UInt8(0)),
            ProteinKey("N", true, UInt8(0)),  # N -> {p6,p2,p7,p4,p5} (Group 2)
            ProteinKey("N", true, UInt8(0)),
            ProteinKey("N", true, UInt8(0)),
            ProteinKey("N", true, UInt8(0)),
            ProteinKey("N", true, UInt8(0)),
            ProteinKey("R", true, UInt8(0)),  # R -> {p6,p7} (Group 4)
            ProteinKey("R", true, UInt8(0)),
            ProteinKey("J", true, UInt8(0)),  # J -> {p6,p2,p7,p4} (Group 6, no p5!)
            ProteinKey("J", true, UInt8(0)),
            ProteinKey("J", true, UInt8(0)),
            ProteinKey("J", true, UInt8(0)),
            ProteinKey("K", true, UInt8(0)),  # K -> {p1,p2,p3,p4,p5} (Group 1)
            ProteinKey("K", true, UInt8(0)),
            ProteinKey("K", true, UInt8(0)),
            ProteinKey("K", true, UInt8(0)),
            ProteinKey("K", true, UInt8(0)),
        ]

        peptides_realworld = [
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR (highly shared!)
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p8", true, UInt8(0)),  # ETYSSYIYK
            PeptideKey("p9", true, UInt8(0)),  # SMSILNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p8", true, UInt8(0)),  # ETYSSYIYK
            PeptideKey("p9", true, UInt8(0)),  # SMSILNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p10", true, UInt8(0)),  # ESYSIYIYK
            PeptideKey("p11", true, UInt8(0)),  # AMSIMNSFVTDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
            PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
            PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
            PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
            PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
        ]

        result = infer_proteins(proteins_realworld, peptides_realworld)
            # Only unique peptides should be in results (shared peptides excluded)
        @test length(result.peptide_to_protein) == 8

        # Unique peptide assignments
        @test result.peptide_to_protein[PeptideKey("p1", true, UInt8(0))].name == "A;C;D;E;F;G;H;I;K"    
        @test !haskey(result.peptide_to_protein,PeptideKey("p2", true, UInt8(0)))  # pep2 shared, excluded  
        @test result.peptide_to_protein[PeptideKey("p3", true, UInt8(0))].name == "A;C;D;E;F;G;H;I;K"     
        @test !haskey(result.peptide_to_protein, PeptideKey("p4", true, UInt8(0)))
        @test !haskey(result.peptide_to_protein, PeptideKey("p5", true, UInt8(0)))
        @test result.peptide_to_protein[PeptideKey("p6", true, UInt8(0))].name == "B;J;L;M;N;Q;R"    
        @test result.peptide_to_protein[PeptideKey("p7", true, UInt8(0))].name == "B;J;L;M;N;Q;R"
        @test result.peptide_to_protein[PeptideKey("p8", true, UInt8(0))].name == "O;P"
        @test result.peptide_to_protein[PeptideKey("p9", true, UInt8(0))].name == "O;P"
        @test result.peptide_to_protein[PeptideKey("p10", true, UInt8(0))].name == "S"
        @test result.peptide_to_protein[PeptideKey("p11", true, UInt8(0))].name == "S"
    end
    @testset "InferenceResult Structure" begin
        # Test that the result structure is correctly formed
        proteins = [ProteinKey("A", true, UInt8(1))]
        peptides = [PeptideKey("pep1", true, UInt8(1))]

        result = infer_proteins(proteins, peptides)

        @test isa(result, InferenceResult)
        @test isa(result.peptide_to_protein, Dictionary{PeptideKey, ProteinKey})
        @test length(result.peptide_to_protein) == 1

        # Check that all unique peptides are present
        for peptide in peptides
            @test haskey(result.peptide_to_protein, peptide)
        end
    end

end