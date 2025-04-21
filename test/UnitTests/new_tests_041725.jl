using Test
using DataFrames
using Combinatorics

@testset "chronologer_prep.jl" begin
    
    #==========================================================================
    Tests for adjustNCE
    ==========================================================================#

    @testset "adjustNCE" begin
        charge_facs = Float64[1.0, 0.9, 0.85, 0.8, 0.75]
        
        # Test adjusting NCE down for higher charge state
        @test adjustNCE(30.0, 2, 3, charge_facs) ≈ 30.0 * (charge_facs[2] / charge_facs[3])
        
        # Test adjusting NCE up for lower charge state
        @test adjustNCE(30.0, 3, 2, charge_facs) ≈ 30.0 * (charge_facs[3] / charge_facs[2])
        
        # Test no adjustment when default charge equals peptide charge
        @test adjustNCE(30.0, 2, 2, charge_facs) ≈ 30.0
        
        # Test edge case with highest charge state
        @test adjustNCE(30.0, 3, 5, charge_facs) ≈ 30.0 * (charge_facs[3] / charge_facs[5])
    end
    #==========================================================================
    Tests for getFixedMods!
    ==========================================================================#
    @testset "getFixedMods!" begin
        # Mock regex matches for testing
        mod_regex = r"C"
        sequence = "PEPTCIDE"
        mod_matches = eachmatch(mod_regex, sequence)
        mod_name = "Carbamidomethyl"
        
        # Test adding fixed mods
        fixed_mods = Vector{PeptideMod}()
        getFixedMods!(fixed_mods, mod_matches, mod_name)
        
        # Should have one fixed mod at position 5 (0-based indexing in regex)
        @test length(fixed_mods) == 1
        @test fixed_mods[1].position == 5
        @test fixed_mods[1].aa == 'C'
        @test fixed_mods[1].mod_name == "Carbamidomethyl"
        
        # Test multiple matches
        mod_regex = r"P"
        sequence = "PEPTIDE"
        mod_matches = eachmatch(mod_regex, sequence)
        mod_name = "Oxidation"
        
        fixed_mods = Vector{PeptideMod}()
        getFixedMods!(fixed_mods, mod_matches, mod_name)
        
        # Should have two fixed mods at positions 1 and 3
        @test length(fixed_mods) == 2
        @test fixed_mods[1].position == 1
        @test fixed_mods[1].aa == 'P'
        @test fixed_mods[1].mod_name == "Oxidation"
        @test fixed_mods[2].position == 3
        @test fixed_mods[2].aa == 'P'
        @test fixed_mods[2].mod_name == "Oxidation"
        
        # Test with no matches
        mod_regex = r"X"
        sequence = "PEPTIDE"
        mod_matches = eachmatch(mod_regex, sequence)
        mod_name = "Oxidation"
        
        fixed_mods = Vector{PeptideMod}()
        getFixedMods!(fixed_mods, mod_matches, mod_name)
        
        # Should have no fixed mods
        @test isempty(fixed_mods)
    end
    
    #==========================================================================
    Tests for matchVarMods
    ==========================================================================#
    @testset "matchVarMods" begin
        # Test with single modification pattern
        sequence = "PEPTCMSK"
        var_mods = [(p=r"M", r="Oxidation")]
        matches = matchVarMods(sequence, var_mods)
        
        @test length(matches) == 1
        @test matches[1].regex_match.match == "M"
        @test matches[1].regex_match.offset == 6
        @test matches[1].name == "Oxidation"
        
        # Test with multiple modification patterns
        var_mods = [(p=r"M", r="Oxidation"), (p=r"[ST]", r="Phospho"), (p=r"C", r="Carbamidomethyl")]
        matches = matchVarMods(sequence, var_mods)
        
        @test length(matches) == 4
        # Sort matches by offset for consistent testing
        sort!(matches, by = x -> x.regex_match.offset)
        
        @test matches[1].regex_match.match == "T"
        @test matches[1].regex_match.offset == 4
        @test matches[1].name == "Phospho"

        @test matches[2].regex_match.match == "C"
        @test matches[2].regex_match.offset == 5
        @test matches[2].name == "Carbamidomethyl"
        
        @test matches[3].regex_match.match == "M"
        @test matches[3].regex_match.offset == 6
        @test matches[3].name == "Oxidation"
        
        @test matches[4].regex_match.match == "S"
        @test matches[4].regex_match.offset == 7
        @test matches[4].name == "Phospho"

        # Test with no matches
        var_mods = [(p=r"X", r="Unknown")]
        matches = matchVarMods(sequence, var_mods)
        @test isempty(matches)
        
        # Test with regex that matches multiple instances
        var_mods = [(p=r"P", r="Proline-mod")]
        matches = matchVarMods(sequence, var_mods)
        @test length(matches) == 2
    end

    #==========================================================================
    Tests for countVarModCombinations
    ==========================================================================#
    @testset "countVarModCombinations" begin
        # Create mock var_mod_matches for testing
        sequence = "PEPTCMSK"
        var_mods = [(p=r"M", r="Oxidation"), (p=r"[ST]", r="Phospho"), (p=r"C", r="Carbamidomethyl")]
        var_mod_matches = matchVarMods(sequence, var_mods)
        
        # Test with max_var_mods = number of matches
        max_var_mods = 3
        n_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        # binomial(4, 3) + binomial(4, 2) + binomial(4, 1) + 1 = 15
        @test n_combinations == 15
        
        # Test with max_var_mods < number of matches
        max_var_mods = 2
        n_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        # binomial(4, 2) + binomial(4, 1) + 1 = 11
        @test n_combinations == 11
        
        # Test with max_var_mods > number of matches
        max_var_mods = 5
        n_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        # binomial(4, 4) + binomial(4, 3) + binomial(4, 2) + binomial(4, 1) + 1 = 16
        @test n_combinations == 16
        
        # Test with no matches
        var_mod_matches = Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}()
        n_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        # Should be 1 (no mods)
        @test n_combinations == 1
    end
    
    #==========================================================================
    Tests for fillVarModStrings!
    ==========================================================================#
    #==========================================================================
    Tests for fillVarModStrings!
    ==========================================================================#
    @testset "fillVarModStrings!" begin
        # Setup for testing
        sequence = "PEPTCMSK"
        var_mods = [(p=r"M", r="Oxidation"), (p=r"S", r="Phospho")]
        var_mod_matches = matchVarMods(sequence, var_mods)
        fixed_mods = [PeptideMod(UInt8(5), 'C', "Carbamidomethyl")]
        max_var_mods = 2
        
        # Calculate number of combinations
        n_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        # Should be 1 (no var mods) + 2 (choose 1) + 1 (choose 2) = 4
        @test n_combinations == 4
        
        # Initialize all_mods to hold all combinations
        all_mods = Vector{Vector{PeptideMod}}(undef, n_combinations)
        
        # Fill all_mods with modification combinations
        fillVarModStrings!(all_mods, var_mod_matches, fixed_mods, max_var_mods)
        
        # Check all combinations are generated
        @test length(all_mods) == 4
        
        # Check that each combination includes fixed mods
        @test all(mod_list -> any(mod -> mod.position == 5 && mod.aa == 'C' && mod.mod_name == "Carbamidomethyl", mod_list), all_mods)
        
        # Verify the combinations (should be fixed mods + various var mods)
        # Find combination with only fixed mods
        fixed_only = findfirst(mod_list -> length(mod_list) == 1, all_mods)
        @test !isnothing(fixed_only)
        @test length(all_mods[fixed_only]) == 1
        @test all_mods[fixed_only][1].position == 5
        @test all_mods[fixed_only][1].aa == 'C'
        @test all_mods[fixed_only][1].mod_name == "Carbamidomethyl"
        
        # Find combination with fixed mods + M oxidation
        m_oxidation = findfirst(mod_list -> 
            any(mod -> mod.position == 6 && mod.aa == 'M' && mod.mod_name == "Oxidation", mod_list) && 
            length(mod_list) == 2, 
            all_mods)
        @test !isnothing(m_oxidation)
        
        # Find combination with fixed mods + S phosphorylation
        s_phospho = findfirst(mod_list -> 
            any(mod -> mod.position == 7 && mod.aa == 'S' && mod.mod_name == "Phospho", mod_list) && 
            length(mod_list) == 2, 
            all_mods)
        @test !isnothing(s_phospho)
        
        # Find combination with all mods
        all_mod_combo = findfirst(mod_list -> length(mod_list) == 3, all_mods)
        @test !isnothing(all_mod_combo)
        @test any(mod -> mod.position == 6 && mod.aa == 'M' && mod.mod_name == "Oxidation", all_mods[all_mod_combo])
        @test any(mod -> mod.position == 7 && mod.aa == 'S' && mod.mod_name == "Phospho", all_mods[all_mod_combo])
        
        # Test with empty fixed mods
        empty_fixed_mods = Vector{PeptideMod}()
        all_mods = Vector{Vector{PeptideMod}}(undef, n_combinations)
        
        fillVarModStrings!(all_mods, var_mod_matches, empty_fixed_mods, max_var_mods)
        
        # Check all combinations are generated
        @test length(all_mods) == 4
        
        # Find combination with no mods
        no_mods = findfirst(mod_list -> isempty(mod_list), all_mods)
        @test !isnothing(no_mods)
        
        # Find combination with only M oxidation
        m_only = findfirst(mod_list -> 
            length(mod_list) == 1 && 
            mod_list[1].position == 6 && 
            mod_list[1].aa == 'M' && 
            mod_list[1].mod_name == "Oxidation", 
            all_mods)
        @test !isnothing(m_only)
        
        # Find combination with only S phosphorylation
        s_only = findfirst(mod_list -> 
            length(mod_list) == 1 && 
            mod_list[1].position == 7 && 
            mod_list[1].aa == 'S' && 
            mod_list[1].mod_name == "Phospho", 
            all_mods)
        @test !isnothing(s_only)
        
        # Find combination with both mods
        both_mods = findfirst(mod_list -> length(mod_list) == 2, all_mods)
        @test !isnothing(both_mods)
        @test any(mod -> mod.position == 6 && mod.aa == 'M' && mod.mod_name == "Oxidation", all_mods[both_mods])
        @test any(mod -> mod.position == 7 && mod.aa == 'S' && mod.mod_name == "Phospho", all_mods[both_mods])
    end
    
    #==========================================================================
    Tests for add_pair_indices!
    ==========================================================================#
    @testset "add_pair_indices!" begin
        # Test with pairs (two rows with same pair_id)
        df = DataFrame(pair_id = UInt32[1, 2, 2, 3, 4, 4, 5])
        add_pair_indices!(df)
        
        @test hasproperty(df, :partner_precursor_idx)
        @test ismissing(df[1, :partner_precursor_idx])
        @test df[2, :partner_precursor_idx] == 3
        @test df[3, :partner_precursor_idx] == 2
        @test ismissing(df[4, :partner_precursor_idx])
        @test df[5, :partner_precursor_idx] == 6
        @test df[6, :partner_precursor_idx] == 5
        @test ismissing(df[7, :partner_precursor_idx])
        
        # Test with no pairs (all unique pair_ids)
        df = DataFrame(pair_id = UInt32[1, 2, 3, 4, 5])
        add_pair_indices!(df)
        
        @test hasproperty(df, :partner_precursor_idx)
        @test all(ismissing.(df[!, :partner_precursor_idx]))
        
        # Test with only pairs (all rows paired)
        df = DataFrame(pair_id = UInt32[1, 1, 2, 2, 3, 3])
        add_pair_indices!(df)
        
        @test hasproperty(df, :partner_precursor_idx)
        @test df[1, :partner_precursor_idx] == 2
        @test df[2, :partner_precursor_idx] == 1
        @test df[3, :partner_precursor_idx] == 4
        @test df[4, :partner_precursor_idx] == 3
        @test df[5, :partner_precursor_idx] == 6
        @test df[6, :partner_precursor_idx] == 5
    end
end

@testset "FASTA Module Tests" begin
    #==========================================================================
    Tests for fasta_digest.jl
    ==========================================================================#
    @testset "digest_sequence" begin
        # Test basic tryptic digest with no missed cleavages
        sequence = "MKVGPKAFRVLTEDEMAKR"
        peptides = digest_sequence(sequence, r"[KR][^P]", 20, 5, 1)
        @test Set(peptides) == Set(["MKVGPK", "VGPKAFR", "VLTEDEMAK", "VLTEDEMAKR", "AFRVLTEDEMAK"])
        
        peptides = digest_sequence(sequence, r"[KR][^P]", 20, 5, 0)
        @test Set(peptides) == Set(["VLTEDEMAK"])

        sequence = "MAKRTGKR"
        peptides = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test Set(peptides) == Set(["MAK", "R", "TGK", "R"])
        
        # Test with missed cleavages
        peptides = digest_sequence(sequence, r"[KR]", 10, 1, 1)
        @test Set(peptides) == Set(["MAK", "MAKR", "R", "RTGK", "TGK", "TGKR", "R"])
        
        # Test with length constraints
        peptides = digest_sequence(sequence, r"[KR]", 5, 3, 0)
        @test Set(peptides) == Set(["MAK", "TGK"])
        
        # Test with regex that has overlap=true
        sequence = "MAKRRMAK"
        peptides = digest_sequence(sequence, r"R", 10, 1, 0)
        @test Set(peptides) == Set(["MAKR","R","MAK"])
        
        # Test with no matches
        sequence = "ACDEFGHIJ"
        peptides = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test Set(peptides) == Set(["ACDEFGHIJ"])
        
        # Test with empty sequence
        sequence = ""
        peptides = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test isempty(peptides)
        
        # Test returned type is Vector{String}, not Vector{SubString}
        sequence = "MAKRTGKR"
        peptides = digest_sequence(sequence, r"[KR]", 10, 1, 0)
        @test eltype(peptides) === String
    end
    
    
    @testset "digest_fasta" begin
        # Create test FASTA entries
        fasta_entries = [
            FastaEntry(
                "P12345", 
                "Test protein 1", 
                "test", 
                "MAKRTGKRPEPT", 
                missing, 
                missing, 
                UInt8(0), 
                UInt32(0), 
                UInt32(0), 
                UInt8(0), 
                false
            ),
            FastaEntry(
                "P67890", 
                "Test protein 2", 
                "test", 
                "XBZMAKHUOPR", # Contains unusual AAs that should be filtered
                missing, 
                missing, 
                UInt8(0), 
                UInt32(0), 
                UInt32(0), 
                UInt8(0), 
                false
            )
        ]
        
        # Test basic digestion
        peptides = digest_fasta(fasta_entries, "human", regex=r"[KR]", max_length=10, min_length=2, missed_cleavages=0)
        
        # Should only digest the first protein (second has unusual AAs)
        @test length(peptides) == 4
        
        # Check peptide sequences
        peptide_seqs = [get_sequence(p) for p in peptides]
        @test sort(peptide_seqs) == sort(["MA", "PEPT", "R", "TG"])
        
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
        @test get_entrapment_group_id(first_peptide) == 0
        @test is_decoy(first_peptide) == false
        
        # Test base_pep_id and base_prec_id increment
        @test get_base_pep_id(peptides[2]) == 2
        @test get_base_prec_id(peptides[2]) == 2
        
        # Test with missed cleavages
        peptides = digest_fasta(fasta_entries, "human", regex=r"[KR]", max_length=10, min_length=2, missed_cleavages=1)
        
        # Should have original 4 peptides plus missed cleavage peptides
        @test length(peptides) > 4
        
        # Check some combined peptides exist
        peptide_seqs = [get_sequence(p) for p in peptides]
        @test "MAR" in peptide_seqs
        @test "RTGK" in peptide_seqs
    end
    #=
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
        GZip.open(temp_fasta_gz, "w") do io
            write(io, fasta_content)
        end
        
        # Test uncompressed file parsing
        entries = parse_fasta(temp_fasta, "human")
        
        # Basic checks
        @test length(entries) == 2
        @test get_sequence(entries[1]) == "MAKRTGKRPEPT"
        @test get_sequence(entries[2]) == "ACDEFGHIJK"
        
        # Check ID parsing (should extract P12345 from UniProt format)
        @test get_id(entries[1]) == "P12345"
        @test get_id(entries[2]) == "P67890"
        
        # Check proteome is set correctly
        @test all(e -> get_proteome(e) == "human", entries)
        
        # Check default values
        @test all(e -> ismissing(get_structural_mods(e)), entries)
        @test all(e -> ismissing(get_isotopic_mods(e)), entries)
        @test all(e -> get_charge(e) == 0, entries)
        @test all(e -> get_base_pep_id(e) == 0, entries)
        @test all(e -> get_base_prec_id(e) == 0, entries)
        @test all(e -> get_entrapment_group_id(e) == 0, entries)
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
        push!(pss, "PEPTIDE")
        @test length(getSeqSet(pss)) == 1
        @test "PEPTIDE" in pss
        
        # Test I/L equivalence
        push!(pss, "PEPTLDE")
        @test length(getSeqSet(pss)) == 1  # Should still be 1 as I and L are equivalent
        @test "PEPTIDE" in pss
        @test "PEPTLDE" in pss
        
        # Test adding multiple sequences
        push!(pss, "ANOTHER")
        @test length(getSeqSet(pss)) == 2
        @test "ANOTHER" in pss
        
        # Test constructor from FastaEntry vector
        entries = [
            FastaEntry("P1", "", "test", "PEPTIDE", missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt8(0), false),
            FastaEntry("P2", "", "test", "ANOTHER", missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt8(0), false),
            FastaEntry("P3", "", "test", "PEPTLDE", missing, missing, UInt8(0), UInt32(0), UInt32(0), UInt8(0), false)
        ]
        
        pss_from_entries = PeptideSequenceSet(entries)
        @test length(getSeqSet(pss_from_entries)) == 2  # PEPTIDE and PEPTLDE are considered the same
        @test "PEPTIDE" in pss_from_entries
        @test "ANOTHER" in pss_from_entries
    end
    
    @testset "shuffle_fast" begin
        # Test basic shuffling preserves C-terminal
        s = "PEPTIDEK"
        shuffled = shuffle_fast(s)
        @test shuffled[end] == 'K'  # Last AA preserved
        @test length(shuffled) == length(s)  # Length preserved
        @test Set(shuffled) == Set(s)  # Same character composition
        
        # Ensure some shuffling actually occurred
        @test shuffled != s
        
        # Test with short strings
        s = "AB"
        shuffled = shuffle_fast(s)
        @test shuffled == s  # Can't shuffle a 2-letter string with last preserved
        
        # Test with single character
        s = "A"
        @test shuffle_fast(s) == "A"
        
        # Test with empty string
        s = ""
        @test shuffle_fast(s) == ""
    end
    
    @testset "add_entrapment_sequences" begin
        # Create test entries
        entries = [
            FastaEntry("P1", "", "test", "PEPTIDEK", missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "", "test", "MAKEPROTEIN", missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false)
        ]
        
        # Test with single entrapment sequence per target
        result = add_entrapment_sequences(entries, UInt8(1))
        
        @test length(result) == 4  # 2 original + 2 entrapment
        
        # Identify entrapment sequences
        entrapment = filter(e -> get_entrapment_group_id(e) > 0, result)
        @test length(entrapment) == 2
        
        # Check properties
        for e in entrapment
            @test get_entrapment_group_id(e) == 1
            @test !is_decoy(e)
            @test endswith(get_sequence(e), get_sequence(entries[get_base_pep_id(e)])[end])
        end
        
        # Test with multiple entrapment sequences
        result = add_entrapment_sequences(entries, UInt8(3))
        
        @test length(result) == 8  # 2 original + 6 entrapment
        
        # Count by entrapment group
        entrapment_groups = [filter(e -> get_entrapment_group_id(e) == i, result) for i in 1:3]
        @test all(length.(entrapment_groups) .== 2)
    end
    
    @testset "add_reverse_decoys" begin
        # Create test entries
        entries = [
            FastaEntry("P1", "", "test", "PEPTIDEK", missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "", "test", "MAKEPROTEIN", missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false)
        ]
        
        # Test basic reversal
        result = add_reverse_decoys(entries)
        
        @test length(result) == 4  # 2 original + 2 decoy
        
        # Identify decoys
        decoys = filter(is_decoy, result)
        @test length(decoys) == 2
        
        # Check decoy sequences are reversed except last AA
        for i in 1:2
            target_seq = get_sequence(entries[i])
            decoy = decoys[i]
            decoy_seq = get_sequence(decoy)
            
            @test decoy_seq[end] == target_seq[end]  # Last AA preserved
            @test decoy_seq[1:end-1] == reverse(target_seq[1:end-1])  # Rest is reversed
            
            # Check metadata preserved
            @test get_base_pep_id(decoy) == get_base_pep_id(entries[i])
            @test get_base_prec_id(decoy) == get_base_prec_id(entries[i])
            @test get_entrapment_group_id(decoy) == get_entrapment_group_id(entries[i])
        end
        
        # Test with modifications (would need to mock PeptideMod)
        # This is more complex and would require additional setup
    end
    
    @testset "combine_shared_peptides" begin
        # Create test entries with shared sequences
        entries = [
            FastaEntry("P1", "desc1", "human", "PEPTIDE", missing, missing, UInt8(0), UInt32(1), UInt32(1), UInt8(0), false),
            FastaEntry("P2", "desc2", "human", "PEPTIDE", missing, missing, UInt8(0), UInt32(2), UInt32(2), UInt8(0), false),
            FastaEntry("P3", "desc3", "mouse", "UNIQUE", missing, missing, UInt8(0), UInt32(3), UInt32(3), UInt8(0), false),
            FastaEntry("P4", "desc4", "human", "PEPTLDE", missing, missing, UInt8(0), UInt32(4), UInt32(4), UInt8(0), false)
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
    =#
end