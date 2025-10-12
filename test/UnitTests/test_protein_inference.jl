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
        @test length(result.use_for_quant) == 0
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
        
        # All peptides should be used for quantification
        @test all(values(result.use_for_quant))
        
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
        
        # Check assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A;B"
        @test result.peptide_to_protein[peptides[3]].name == "A;B"
        @test result.peptide_to_protein[peptides[4]].name == "B"
        
        # Check quantification flags
        @test result.use_for_quant[peptides[1]] == true   # Unique to A
        @test result.use_for_quant[peptides[2]] == false  # Shared peptide
        @test result.use_for_quant[peptides[3]] == false  # Shared peptide
        @test result.use_for_quant[peptides[4]] == true   # Unique to B
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
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A;B"
            @test result.use_for_quant[peptide] == true
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
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A"
            @test result.use_for_quant[peptide] == true
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
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A"  # A can explain pep2
        @test result.peptide_to_protein[peptides[3]].name == "C"  # C can explain pep3  
        @test result.peptide_to_protein[peptides[4]].name == "C"
        
        # All peptides should be used for quantification in this case
        # since the necessary proteins (A and C) can uniquely explain all peptides
        @test result.use_for_quant[peptides[1]] == true
        @test result.use_for_quant[peptides[2]] == true  
        @test result.use_for_quant[peptides[3]] == true  
        @test result.use_for_quant[peptides[4]] == true
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
        for peptide in peptides
            @test result.peptide_to_protein[peptide].name == "A;B;C"
            @test result.use_for_quant[peptide] == true
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
        @test result.peptide_to_protein[peptides[1]].entrap_id == UInt8(1)
        @test result.peptide_to_protein[peptides[2]].entrap_id == UInt8(2)
        @test result.use_for_quant[peptides[1]] == true
        @test result.use_for_quant[peptides[2]] == true
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
        @test result.peptide_to_protein[peptides[1]].is_target == true
        @test result.peptide_to_protein[peptides[2]].is_target == false
        @test result.use_for_quant[peptides[1]] == true
        @test result.use_for_quant[peptides[2]] == true
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
        
        # Case A results
        @test result.peptide_to_protein[peptides[1]].name == "A"
        @test result.peptide_to_protein[peptides[2]].name == "A"
        @test result.peptide_to_protein[peptides[3]].name == "B"
        @test result.peptide_to_protein[peptides[4]].name == "B"
        @test all(result.use_for_quant[pep] for pep in peptides[1:4])
        
        # Case B-like results (differentiable proteins C, D - C has unique pep5, D has unique pep8)
        @test result.peptide_to_protein[peptides[5]].name == "C"
        @test result.peptide_to_protein[peptides[6]].name == "C;D"  # Shared peptide
        @test result.peptide_to_protein[peptides[7]].name == "C;D"  # Shared peptide
        @test result.peptide_to_protein[peptides[8]].name == "D"
        @test result.use_for_quant[peptides[5]] == true
        @test result.use_for_quant[peptides[6]] == false  # Shared, not used for quant
        @test result.use_for_quant[peptides[7]] == false  # Shared, not used for quant
        @test result.use_for_quant[peptides[8]] == true
        
        # Case C results
        @test result.peptide_to_protein[peptides[9]].name == "E;F"
        @test result.peptide_to_protein[peptides[10]].name == "E;F"
        @test result.peptide_to_protein[peptides[11]].name == "E;F"
        @test result.peptide_to_protein[peptides[12]].name == "E;F"
        @test all(result.use_for_quant[pep] for pep in peptides[9:12])
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

        # Peptide assignments
        @test result.peptide_to_protein[peptides[1]].name == "A"      # pep1 unique to A
        @test result.peptide_to_protein[peptides[2]].name == "A"      # pep2 uniquely explained by A in selected set
        @test result.peptide_to_protein[peptides[3]].name == "B;C"    # pep3 - B and C indistinguishable, merged
        @test result.peptide_to_protein[peptides[4]].name == "C;D"    # pep4 shared between B;C and D
        @test result.peptide_to_protein[peptides[5]].name == "D"      # pep5 unique to D
        @test result.peptide_to_protein[peptides[6]].name == "D"      # pep6 uniquely explained by D in selected set
        @test result.peptide_to_protein[peptides[7]].name == "E;F"    # pep7 - E and F indistinguishable, merged
        @test result.peptide_to_protein[peptides[8]].name == "A;F"    # pep8 shared between A and E;F

        # Quantification flags
        @test result.use_for_quant[peptides[1]] == true   # pep1 unique to A
        @test result.use_for_quant[peptides[2]] == true   # pep2 uniquely explained by A in selected set
        @test result.use_for_quant[peptides[3]] == true   # pep3 unique to B;C merged group
        @test result.use_for_quant[peptides[4]] == false  # pep4 shared between B;C and D
        @test result.use_for_quant[peptides[5]] == true   # pep5 unique to D
        @test result.use_for_quant[peptides[6]] == true   # pep6 uniquely explained by D in selected set
        @test result.use_for_quant[peptides[7]] == true   # pep7 unique to E;F merged group
        @test result.use_for_quant[peptides[8]] == false  # pep8 shared between A and E;F
    end

    @testset "InferenceResult Structure" begin
        # Test that the result structure is correctly formed
        proteins = [ProteinKey("A", true, UInt8(1))]
        peptides = [PeptideKey("pep1", true, UInt8(1))]

        result = infer_proteins(proteins, peptides)

        @test isa(result, InferenceResult)
        @test isa(result.peptide_to_protein, Dictionary{PeptideKey, ProteinKey})
        @test isa(result.use_for_quant, Dictionary{PeptideKey, Bool})
        @test length(result.peptide_to_protein) == 1
        @test length(result.use_for_quant) == 1

        # Check that all peptides have both mappings
        for peptide in peptides
            @test haskey(result.peptide_to_protein, peptide)
            @test haskey(result.use_for_quant, peptide)
        end
    end

end