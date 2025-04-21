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

