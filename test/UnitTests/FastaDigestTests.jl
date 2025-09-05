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

@testset "FASTA Module Tests" begin
    #==========================================================================
    Tests for fasta_digest.jl
    ==========================================================================#
    @testset "digest_sequence" begin
        # Test basic tryptic digest with no missed cleavages
        sequence = "MKVGPKAFRVLTEDEMAKR"
        peptides, starts = digest_sequence(sequence, r"[KR][^P]", 20, 5, 1)
        @test Set(peptides) == Set(["MKVGPK", "VGPKAFR", "VLTEDEMAK", "VLTEDEMAKR", "AFRVLTEDEMAK"])
        @test Set(starts) == Set([1, 3, 7, 10])
        
        peptides, starts = digest_sequence(sequence, r"[KR][^P]", 20, 5, 0)
        @test Set(peptides) == Set(["VLTEDEMAK"])
        @test Set(starts) == Set([10])

        sequence = "MAKRTGKR"
        peptides, starts = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test Set(peptides) == Set(["MAK", "R", "TGK", "R"])
        @test Set(starts) == Set([1, 4, 5, 8])
        
        # Test with missed cleavages
        peptides, starts = digest_sequence(sequence, r"[KR]", 10, 1, 1)
        @test Set(peptides) == Set(["MAK", "MAKR", "R", "RTGK", "TGK", "TGKR", "R"])
        @test Set(starts) == Set([1, 4, 5, 8])
        
        # Test with length constraints
        peptides, starts = digest_sequence(sequence, r"[KR]", 5, 3, 0)
        @test Set(peptides) == Set(["MAK", "TGK"])
        @test Set(starts) == Set([1, 5])
        
        # Test with regex that has overlap=true
        sequence = "MAKRRMAK"
        peptides, starts = digest_sequence(sequence, r"R", 10, 1, 0)
        @test Set(peptides) == Set(["MAKR","R","MAK"])
        @test Set(starts) == Set([1, 5, 6])
        
        # Test with no matches
        sequence = "ACDEFGHIJ"
        peptides, starts = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test Set(peptides) == Set(["ACDEFGHIJ"])
        @test Set(starts) == Set([1])
        
        # Test with empty sequence
        sequence = ""
        peptides, starts = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test isempty(peptides)
        @test isempty(starts)
        
        # Test returned type is Vector{String}, not Vector{SubString}
        sequence = "MAKRTGKR"
        peptides, starts = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test eltype(peptides) === String
    end
    
    
    @testset "digest_fasta" begin
        # Create test FASTA entries
        fasta_entries = [
            FastaEntry(
                "P12345",
                "Test protein 1",
                "",
                "",
                "human",
                "test",
                "MAKRTGKRPEPT",
                UInt32(1),
                missing,
                missing,
                UInt8(0),
                UInt32(0),  # base_seq_id
                UInt32(0),  # base_pep_id
                UInt32(0),  # base_prec_id
                UInt8(0),   # entrapment_pair_id
                false,      # is_decoy
            ),
            FastaEntry(
                "P67890",
                "Test protein 2",
                "",
                "",
                "human",
                "test",
                "XBZMAKHUOPR", # Contains unusual AAs that should be filtered
                UInt32(1),
                missing,
                missing,
                UInt8(0),
                UInt32(0),  # base_seq_id
                UInt32(0),  # base_pep_id
                UInt32(0),  # base_prec_id
                UInt8(0),   # entrapment_pair_id
                false,      # is_decoy
            )
        ]
        
        # Test basic digestion
        peptides = digest_fasta(fasta_entries, "human", regex=r"[KR]", max_length=10, min_length=2, missed_cleavages=0)
        
        # Should only digest the first protein (second has unusual AAs)
        @test length(peptides) == 3
        
        # Check peptide sequences
        peptide_seqs = [get_sequence(p) for p in peptides]
        @test sort(peptide_seqs) == sort(["MAK", "PEPT", "TGK"])
        
        # Check properties are correctly set
        first_peptide = peptides[1]
        @test get_id(first_peptide) == "P12345"
        @test get_description(first_peptide) == "" # Description is empty to save memory
        @test get_proteome(first_peptide) == "human"
        @test ismissing(get_structural_mods(first_peptide))
        @test ismissing(get_isotopic_mods(first_peptide))
        @test get_charge(first_peptide) == 0
        @test get_base_pep_id(first_peptide) == 1
        @test get_base_prec_id(first_peptide) == 1
        @test get_entrapment_pair_id(first_peptide) == 0
        @test is_decoy(first_peptide) == false
        
        # Test base_pep_id and base_prec_id increment
        @test get_base_pep_id(peptides[2]) == 2
        @test get_base_prec_id(peptides[2]) == 2
        
        # Test with missed cleavages
        peptides = digest_fasta(fasta_entries, "human", regex=r"[KR]", max_length=10, min_length=2, missed_cleavages=1)
        
        # Should have original 4 peptides plus missed cleavage peptides
        @test length(peptides) == 7
        
        # Check some combined peptides exist
        peptide_seqs = [get_sequence(p) for p in peptides]
        @test Set(["MAK", "MAKR", "TGK", "RTGK", "TGKR", "PEPT", "RPEPT"]) == Set(peptide_seqs)
    end

    #==========================================================================
    Tests for fasta_parser.jl
    ==========================================================================#
    @testset "parse_fasta" begin
        # Create temporary FASTA files for testing
        temp_fasta = tempname() * ".fasta"
        temp_fasta_gz = tempname() * ".fasta.gz"
        
        # Write test data
        fasta_content = """
        >sp|P12345|TEST_HUMAN Test protein 1
        MAKRTGKRPEPT
        >P67890 Test protein 2
        ACDEFGHIJK
        """
        
        # Write uncompressed file
        open(temp_fasta, "w") do io
            write(io, fasta_content)
        end
        
        # Write compressed file
        open(temp_fasta_gz, "w") do io
            gz_stream = GzipCompressorStream(io)
            write(gz_stream, fasta_content)
            close(gz_stream)  # Important to close the stream to flush any buffered data
        end
        
        # Test uncompressed file parsing
        entries = parse_fasta(temp_fasta, "human")
        
        # Basic checks
        @test length(entries) == 2
        @test get_sequence(entries[1]) == "MAKRTGKRPEPT"
        @test get_sequence(entries[2]) == "ACDEFGHIJK"

        # Without regex the entire header should be kept
        @test get_description(entries[1]) == "sp|P12345|TEST_HUMAN Test protein 1"
        @test get_description(entries[2]) == "P67890 Test protein 2"
        
        # Without regex everything up to the first space should be kept
        @test get_id(entries[1]) == "sp|P12345|TEST_HUMAN"
        @test get_id(entries[2]) == "P67890"
        
        # Check proteome is set correctly
        @test all(e -> get_proteome(e) == "human", entries)
        
        # Check default values
        @test all(e -> ismissing(get_structural_mods(e)), entries)
        @test all(e -> ismissing(get_isotopic_mods(e)), entries)
        @test all(e -> get_charge(e) == 0, entries)
        @test all(e -> get_base_pep_id(e) == 0, entries)
        @test all(e -> get_base_prec_id(e) == 0, entries)
        @test all(e -> get_entrapment_pair_id(e) == 0, entries)
        @test all(e -> !is_decoy(e), entries)
        
        # Test compressed file parsing
        entries_gz = parse_fasta(temp_fasta_gz, "mouse")
        
        # Should have same data but different proteome
        @test length(entries_gz) == 2
        @test get_sequence(entries_gz[1]) == "MAKRTGKRPEPT"
        @test get_sequence(entries_gz[2]) == "ACDEFGHIJK"
        @test all(e -> get_proteome(e) == "mouse", entries_gz)
        
        # Clean up
        rm(temp_fasta)
        rm(temp_fasta_gz)
        
        # Test invalid file extension
        @test_throws ErrorException parse_fasta(tempname() * ".txt", "human")
    end
    
    #==========================================================================
    Tests for fasta_utils.jl
    ==========================================================================#
    @testset "PeptideSequenceSet" begin
        # Test empty constructor
        pss = PeptideSequenceSet()
        @test isempty(getSeqSet(pss))
        
        # Test adding a sequence
        charge = one(UInt8)
        push!(pss, "PEPTIDE", charge)
        @test length(getSeqSet(pss)) == 1
        @test ("PEPTIDE", charge) in pss
        
        # Test I/L equivalence
        push!(pss, "PEPTLDE", charge)
        @test length(getSeqSet(pss)) == 1  # Should still be 1 as I and L are equivalent
        @test ("PEPTIDE",charge) in pss
        @test ("PEPTLDE",charge) in pss
        @test ("PEPTLDE",UInt8(2)) ∉ pss
        # Test adding multiple sequences
        push!(pss, "ANOTHER", charge)
        @test length(getSeqSet(pss)) == 2
        @test ("ANOTHER", charge) in pss
        
        # Test constructor from FastaEntry vector
        entries = [
            FastaEntry("P1", "", "", "", "human", "test", "PEPTIDE", UInt32(1), missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt32(0), UInt8(0), false),
            FastaEntry("P2", "", "", "", "human", "test", "ANOTHER", UInt32(1), missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt32(0), UInt8(0), false),
            FastaEntry("P3", "", "", "", "human", "test", "PEPTLDE", UInt32(1), missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt32(0), UInt8(0), false)
        ]
        
        pss_from_entries = PeptideSequenceSet(entries)
        @test length(getSeqSet(pss_from_entries)) == 2  # PEPTIDE and PEPTLDE are considered the same
        @test ("PEPTIDE", zero(UInt8)) in pss_from_entries
        @test ("ANOTHER", zero(UInt8)) in pss_from_entries
    end
    
    @testset "add_entrapment_sequences" begin
        # Create test entries
        entries = [
            FastaEntry("P1", "", "", "", "human", "test", "PEPTIDEK", UInt32(1), missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "", "", "", "human", "test", "MAKEPROTEIN", UInt32(1), missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false)
        ]
        
        # Test with single entrapment sequence per target
        result = add_entrapment_sequences(entries, UInt8(1))
        
        @test length(result) == 4  # 2 original + 2 entrapment
        
        # Identify entrapment sequences
        entrapment = filter(e -> get_entrapment_pair_id(e) > 0, result)
        @test length(entrapment) == 2
        
        # Check properties
        for e in entrapment
            @test get_entrapment_pair_id(e) == 1
            @test !is_decoy(e)
            @test endswith(get_sequence(e), get_sequence(entries[get_base_pep_id(e)])[end])
        end
        
        # Test with multiple entrapment sequences
        result = add_entrapment_sequences(entries, UInt8(3))
        
        @test length(result) == 8  # 2 original + 6 entrapment
        
        # Count by entrapment group
        entrapment_groups = [filter(e -> get_entrapment_pair_id(e) == i, result) for i in 1:3]
        @test all(length.(entrapment_groups) .== 2)
    end
    
    @testset "add_decoy_sequences" begin
        # Create test entries
        entries = [
            FastaEntry("P1", "", "", "", "human", "test", "PEPTIDEK", UInt32(1), missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "", "", "", "human", "test", "MAKEPROTEIN", UInt32(1), missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false)
   ]
        
        # Test basic reversal
        result = add_decoy_sequences(entries)
        
        @test length(result) == 4  # 2 original + 2 decoy
        
        # Identify decoys
        decoys = filter(is_decoy, result)
        @test length(decoys) == 2
        
        # Check decoy sequences are reversed except last AA
        sort!(decoys, by = x-> x.base_prec_id)
        sort!(entries, by = x->  x.base_prec_id)
        for i in 1:2
            
            target_seq = get_sequence(entries[i])
            decoy = decoys[i]
            decoy_seq = get_sequence(decoy)
            @test decoy_seq[end] == target_seq[end]  # Last AA preserved
            @test all(decoy_seq .== target_seq) == false
            @test Set(decoy_seq[1:end-1]) == Set(target_seq[1:end-1])  # Rest is reversed
            
            # Check metadata preserved
            @test get_base_pep_id(decoy) == get_base_pep_id(entries[i])
            @test get_base_prec_id(decoy) == get_base_prec_id(entries[i])
            @test get_entrapment_pair_id(decoy) == get_entrapment_pair_id(entries[i])
        end
    end
    
    @testset "combine_shared_peptides" begin
        # Create test entries with shared sequences
        entries = [
            FastaEntry("P1", "desc1", "", "", "human", "human", "PEPTIDE", UInt32(1), missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "desc2", "", "", "human", "human", "PEPTIDE", UInt32(1), missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false),
            FastaEntry("P3", "desc3", "", "", "mouse", "mouse", "UNIQUE", UInt32(1), missing, missing, UInt8(0), UInt32(3), UInt32(3), UInt8(0), false),
            FastaEntry("P4", "desc4", "", "", "human", "human", "PEPTLDE", UInt32(1), missing, missing, UInt8(0), UInt32(4), UInt32(4), UInt8(0), false)
        ]
        
        # Test combination
        result = combine_shared_peptides(entries)
        
        # Should have 2 unique peptides (PEPTIDE/PEPTLDE are considered identical, UNIQUE is separate)
        @test length(result) == 2
        
        # Find combined entry
        combined = findfirst(e -> occursin(';', get_id(e)), result)
        @test !isnothing(combined)
        
        # Check combined properties
        combined_entry = result[combined]
        @test Set(split(get_id(combined_entry), ';')) == Set(["P1", "P2", "P4"])
        @test any(desc -> occursin(desc, get_description(combined_entry)), ["desc1", "desc2", "desc4"])
        
        # Find unique entry
        unique_entry = result[result .!= [combined_entry]][1]
        @test get_id(unique_entry) == "P3"
        @test get_description(unique_entry) == "desc3"
        @test get_proteome(unique_entry) == "mouse"
    end
    
end