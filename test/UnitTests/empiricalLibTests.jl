
@testset "BasicEmpiricalLibrary Tests" begin
    # Test constructor and column renaming
    println("dir ", @__DIR__)
    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    df = getDF(lib)
    
    @test names(df)[1] == "precursor_idx"  # Test precursor_idx is first column
    @test "prec_mz" in names(df)  # Test PrecursorMz was renamed
    @test "irt" in names(df)  # Test Tr_recalibrated was renamed
    
    # Test precursor_idx assignment
    @test length(unique(df.precursor_idx)) == 3  # Should have 3 unique precursors
    @test all(df[df.modified_sequence .== "PEPT(+80)IDE", :precursor_idx] .== 
              first(df[df.modified_sequence .== "PEPT(+80)IDE", :precursor_idx]))

    # Test nested sorting with different rt_bin_tolerances
    nestedLibrarySort!(lib, rt_bin_tol = 1.5)
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [500.5, 500.5, 500.5, 600.7, 600.7, 300.2],
        ModifiedPeptide = ["PEPT(+80)IDE", "PEPT(+80)IDE", "PEPT(+80)IDE", "SAMPL(+16)E", "SAMPL(+16)E", "TEST"],
        PeptideSequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE", "SAMPLE", "SAMPLE", "TEST"],
        PrecursorCharge = [2, 2, 2, 3, 3, 2],
        Tr_recalibrated = [10.5, 10.5, 10.5, 9.2, 9.2, 20.1],
        LibraryIntensity = [1000.0, 500.0, 100.0, 800.0, 200.0, 600.0],
        ProductMz = [200.1, 300.2, 400.3, 250.1, 350.2, 450.3],
        FragmentType = ["y", "b", "y", "y", "b", "y"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [2, 3, 4, 2, 3, 4]
    )
    @test all(df.prec_mz .≈ test_df_a.PrecursorMz)
    @test all(df.modified_sequence .== test_df_a.ModifiedPeptide)
    @test all(df.prec_charge .≈ test_df_a.PrecursorCharge)
    @test all(df.irt .≈ test_df_a.Tr_recalibrated)
    @test all(df.library_intensity .≈ test_df_a.LibraryIntensity)
    @test all(df.frag_mz .≈ test_df_a.ProductMz)
    @test all(df.frag_type .== test_df_a.FragmentType)
    @test all(df.frag_charge .≈ test_df_a.FragmentCharge)
    @test all(df.frag_series_number .≈ test_df_a.FragmentSeriesNumber)

    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    nestedLibrarySort!(lib, rt_bin_tol = 1.0)
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [600.7, 600.7, 500.5, 500.5, 500.5, 300.2],
        ModifiedPeptide = ["SAMPL(+16)E", "SAMPL(+16)E", "PEPT(+80)IDE", "PEPT(+80)IDE", "PEPT(+80)IDE", "TEST"],
        PeptideSequence = ["SAMPLE", "SAMPLE", "PEPTIDE", "PEPTIDE", "PEPTIDE", "TEST"],
        PrecursorCharge = [3, 3, 2, 2, 2, 2],
        Tr_recalibrated = [9.2, 9.2, 10.5, 10.5, 10.5, 20.1],
        LibraryIntensity = [800.0, 200.0, 1000.0, 500.0, 100.0, 600.0],
        ProductMz = [250.1, 350.2, 200.1, 300.2, 400.3, 450.3],
        FragmentType = ["y", "b", "y", "b", "y", "y"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [2, 3, 2, 3, 4, 4]
     )

    @test all(df.prec_mz .≈ test_df_a.PrecursorMz)
    @test all(df.modified_sequence .== test_df_a.ModifiedPeptide)
    @test all(df.prec_charge .≈ test_df_a.PrecursorCharge)
    @test all(df.irt .≈ test_df_a.Tr_recalibrated)
    @test all(df.library_intensity .≈ test_df_a.LibraryIntensity)
    @test all(df.frag_mz .≈ test_df_a.ProductMz)
    @test all(df.frag_type .== test_df_a.FragmentType)
    @test all(df.frag_charge .≈ test_df_a.FragmentCharge)
    @test all(df.frag_series_number .≈ test_df_a.FragmentSeriesNumber)

    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    nestedLibrarySort!(lib, rt_bin_tol = typemax(Float64))
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [300.2, 500.5, 500.5, 500.5, 600.7, 600.7],
        ModifiedPeptide = ["TEST", "PEPT(+80)IDE", "PEPT(+80)IDE", "PEPT(+80)IDE", "SAMPL(+16)E", "SAMPL(+16)E"],
        PeptideSequence = ["TEST", "PEPTIDE", "PEPTIDE", "PEPTIDE", "SAMPLE", "SAMPLE"],
        PrecursorCharge = [2, 2, 2, 2, 3, 3],
        Tr_recalibrated = [20.1, 10.5, 10.5, 10.5, 9.2, 9.2],
        LibraryIntensity = [600.0, 1000.0, 500.0, 100.0, 800.0, 200.0],
        ProductMz = [450.3, 200.1, 300.2, 400.3,  250.1, 350.2],
        FragmentType = ["y", "y", "b", "y", "y", "b"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [4, 2, 3, 4, 2, 3]
     )

    @test all(df.prec_mz .≈ test_df_a.PrecursorMz)
    @test all(df.modified_sequence .== test_df_a.ModifiedPeptide)
    @test all(df.prec_charge .≈ test_df_a.PrecursorCharge)
    @test all(df.irt .≈ test_df_a.Tr_recalibrated)
    @test all(df.library_intensity .≈ test_df_a.LibraryIntensity) #
    @test all(df.frag_mz .≈ test_df_a.ProductMz) #
    @test all(df.frag_type .== test_df_a.FragmentType) #
    @test all(df.frag_charge .≈ test_df_a.FragmentCharge)
    @test all(df.frag_series_number .≈ test_df_a.FragmentSeriesNumber) #
end

@testset "ParseSpecLib Tests" begin
    @testset "matchVarMods" begin
        sequence = "PEPMTIDME"
        var_mods = [(p=r"M", r="Unimod:35")]
        matches = matchVarMods(sequence, var_mods)
        @test length(matches) == 2  # Two M residues
        @test matches[1].regex_match.offset == 4  # First M at position 4
        @test matches[2].regex_match.offset == 8  # Second M at position 8
        @test matches[1].name == "Unimod:35"
        @test matches[2].name == "Unimod:35"
        
        # Multiple modification patterns
        var_mods = [(p=r"M", r="Unimod:35"), (p=r"[ST]", r="Unimod:21")]
        matches = matchVarMods(sequence, var_mods)
        @test length(matches) == 3  # Two M and one T residue
        # Find the T modification
        t_mod = findfirst(m -> m.regex_match.match == "T", matches)
        @test !isnothing(t_mod)
        @test matches[t_mod].name == "Unimod:21"
    end
    
    @testset "countVarModCombinations" begin
        sequence = "PEPMTIDME"
        var_mods = [(p=r"M", r="Unimod:35")]
        matches = matchVarMods(sequence, var_mods)
        
        # 2 modifications (2 M residues) with max 2 allowed
        n_combinations = countVarModCombinations(matches, 2)
        @test n_combinations == 4  # unmodified + 2 single mods + 1 double mod
        
        # 2 modifications with max 1 allowed
        n_combinations = countVarModCombinations(matches, 1)
        @test n_combinations == 3  # unmodified + 2 single mods
    end
    
    @testset "fillVarModStrings!" begin
        sequence = "PEPMTIDME"
        var_mods = [(p=r"M", r="Unimod:35")]
        matches = matchVarMods(sequence, var_mods)
        n_combinations = countVarModCombinations(matches, 2)
        var_mod_strings = Vector{String}(undef, n_combinations)
        fillVarModStrings!(var_mod_strings, matches, "", 2)
        
        # Sort for reliable comparison
        sort!(var_mod_strings)
        
        @test "" in var_mod_strings  # Unmodified version
        @test "(4,M,Unimod:35)" in var_mod_strings  # Single mod on first M
        @test "(8,M,Unimod:35)" in var_mod_strings  # Single mod on second M
        @test "(4,M,Unimod:35)(8,M,Unimod:35)" in var_mod_strings  # Both Ms modified
        
        # Test with fixed modifications
        var_mod_strings = Vector{String}(undef, n_combinations)
        fillVarModStrings!(var_mod_strings, matches, "(1,P,Unimod:4)", 2)
        
        # Sort for reliable comparison
        sort!(var_mod_strings)
        
        @test "(1,P,Unimod:4)" in var_mod_strings  # Only fixed mod
        @test "(1,P,Unimod:4)(4,M,Unimod:35)" in var_mod_strings  # Fixed + first M
        @test "(1,P,Unimod:4)(8,M,Unimod:35)" in var_mod_strings  # Fixed + second M
        @test "(1,P,Unimod:4)(4,M,Unimod:35)(8,M,Unimod:35)" in var_mod_strings  # Fixed + both Ms
    end
    
    @testset "parseEmpiricalLibraryMods" begin
        # Test case 1: Standard modification
        sequence = "I(tag6)LSISADI(Unimod:35)ETIGEILK(tag6)"
        result = parseEmpiricalLibraryMods(sequence)
        @test result == "(1,I,tag6)(8,I,Unimod:35)(16,K,tag6)"
        
        # Test case 2: N-terminal modification
        sequence = "n(tag6)PEPTIDE(tag6)"
        result = parseEmpiricalLibraryMods(sequence)
        @test result == "(1,n,tag6)(7,E,tag6)"
        
        # Test case 3: N and C-terminal modifications
        sequence = "n(tag6)PEPTIDEc(tag6)"
        result = parseEmpiricalLibraryMods(sequence)
        @test result == "(1,n,tag6)(7,c,tag6)"
    end
    
    @testset "getMassOffset" begin
        mod_masses = [229.163f0, 0.0f0, 0.0f0, 15.995f0, 229.163f0]
        
        # Test single position
        @test getMassOffset(1, 1, mod_masses) ≈ 229.163f0
        
        # Test range
        @test getMassOffset(1, 3, mod_masses) ≈ 229.163f0
        @test getMassOffset(3, 5, mod_masses) ≈ 15.995f0 + 229.163f0
        @test getMassOffset(1, 5, mod_masses) ≈ 229.163f0 + 15.995f0 + 229.163f0
    end
    
    @testset "addIsoMods" begin
        # Test case 1: Basic indexing
        structural_mods = "(1,I,tag6)(5,L,tag5)(16,K,tag6)"
        isotopic_mods = ""
        result = addIsoMods(isotopic_mods, structural_mods, "tag6", (channel="d0", mass=0.0f0))
        @test result == "(1, d0)(3, d0)"
        
        # Test case 2: With existing isotopic mods
        structural_mods = "(1,I,tag6)(5,L,tag5)(16,K,tag6)"
        isotopic_mods = "(2, d4)"  # For the second mod (tag5)
        result = addIsoMods(isotopic_mods, structural_mods, "tag6", (channel="d0", mass=0.0f0))
        @test result == "(1, d0)(2, d4)(3, d0)"
        
        # Test case 3: Missing target modification
        structural_mods = "(1,I,tag7)(5,L,tag5)"  # No tag6
        isotopic_mods = "(2, d4)"
        result = addIsoMods(isotopic_mods, structural_mods, "tag6", (channel="d0", mass=0.0f0))
        @test result == "(2, d4)"  # Unchanged
        
        # Test case 4: Empty structural mods
        structural_mods = ""
        isotopic_mods = ""
        result = addIsoMods(isotopic_mods, structural_mods, "tag6", (channel="d0", mass=0.0f0))
        @test result == ""  # Empty
    end
    
    @testset "getIsoModMasses!" begin
        iso_mods_dict = Dict(
            "tag6" => Dict("d0" => 0.0f0, "d4" => 4.0f0),
            "tag5" => Dict("d0" => 0.0f0, "d4" => 4.0f0)
        )
        
        # Test 1: Basic indexing
        structural_mods = "(1,I,tag6)(5,L,tag5)(7,K,tag6)"
        isotopic_mods = "(3, d4)"  # Third mod gets d4
        iso_mod_masses = zeros(Float32, 255)
        getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
        @test iso_mod_masses[7] ≈ 4.0f0  # Position 7 (K) has mass 4.0
        @test iso_mod_masses[1] ≈ 0.0f0  # Position 1 (I) has mass 0.0
        
        # Test 2: Multiple mods
        structural_mods = "(1,n,tag6)(5,L,tag5)(16,c,tag6)"
        isotopic_mods = "(1, d0)(2, d4)"  # First mod gets d0, second gets d4
        iso_mod_masses = zeros(Float32, 255)
        getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
        @test iso_mod_masses[1] ≈ 0.0f0  # Position 1 (n) has mass 0.0
        @test iso_mod_masses[5] ≈ 4.0f0  # Position 5 (L) has mass 4.0
    end
    
    @testset "get_aa_masses!" begin
        # Test for a simple tripeptide
        sequence = "PAK"
        aa_masses = zeros(Float32, 3)
        get_aa_masses!(aa_masses, sequence)
        @test aa_masses ≈ [97.05276f0, 71.03711f0, 128.09496f0]
        
        # Test for a sequence with all standard amino acids
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        aa_masses = zeros(Float32, 20)
        get_aa_masses!(aa_masses, sequence)
        @test aa_masses[1] ≈ 71.03711f0  # A
        @test aa_masses[8] ≈ 113.08406f0  # I
        @test aa_masses[20] ≈ 163.06333f0  # Y
    end
    
    @testset "reverseSequence" begin
        # Test 1: Basic reversal with terminal modifications
        sequence = "PEPTIDE"
        mods = "(4,T,Phospho)(7,E,Acetyl)"
        rev_seq, rev_mods = reverseSequence(sequence, mods)
        @test rev_seq == "DITPEPE"
        @test rev_mods == "(3,T,Phospho)(7,E,Acetyl)"

        # Test 2: Reversal with N-terminal modification
        sequence = "PEPTIDE"
        mods = "(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)"
        rev_seq, rev_mods = reverseSequence(sequence, mods)
        @test rev_seq == "DITPEPE"
        @test rev_mods == "(3,T,Phospho)(6,P,mymod)(7,E,Acetyl)"

        # Test 3: Reversal with terminal specifiers
        sequence = "PEPTIDE"
        mods = "(1,n,mymod-nterm)(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
        rev_seq, rev_mods = reverseSequence(sequence, mods)
        @test rev_seq == "DITPEPE"
        @test rev_mods == "(1,n,mymod-nterm)(3,T,Phospho)(6,P,mymod)(7,E,Acetyl)(7,c,mymod-cterm)"
    end
    
    @testset "shuffleSequence" begin
        using Random
        # Ensure deterministic results for testing
        Random.seed!(1844)
        
        # Test with terminal modifications
        sequence = "PEPTIDE"
        mods = "(1,n,mymod-nterm)(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
        shuffled_seq, shuffled_mods = shuffleSequence(sequence, mods)
        
        # With this seed, sequence should be "DPTPIEE"
        @test shuffled_seq == "DPTPIEE"
        @test shuffled_mods == "(1,n,mymod-nterm)(2,P,mymod)(3,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
        
        # Reset seed for consistency
        Random.seed!(1844)
        
        # Test with no modifications
        sequence = "ABCDEF"
        mods = ""
        shuffled_seq, shuffled_mods = shuffleSequence(sequence, mods)
        @test length(shuffled_seq) == length(sequence)
        @test last(shuffled_seq) == last(sequence)  # Last character should be preserved
        @test shuffled_mods == ""  # No mods
    end


    @testset "get_aa_masses!" begin
        # Test for a simple tripeptide
        sequence = "PAK"
        aa_masses = zeros(Float32, 3)
        get_aa_masses!(aa_masses, sequence)
        @test aa_masses ≈ [97.05276f0, 71.03711f0, 128.09496f0]
        
        # Test for a sequence with all standard amino acids
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        aa_masses = zeros(Float32, 20)
        get_aa_masses!(aa_masses, sequence)
        @test aa_masses[1] ≈ 71.03711f0  # A
        @test aa_masses[8] ≈ 113.08406f0  # I
        @test aa_masses[20] ≈ 163.06333f0  # Y
    end
    
    @testset "get_structural_mod_masses!" begin
        sequence = "PEPTIDE"
        structural_mods = "(1,P,Phospho)(7,E,Acetyl)"
        mod_masses = zeros(Float32, length(sequence))
        mod_to_mass = Dict("Phospho" => 79.966331f0, "Acetyl" => 42.010565f0)
        
        get_structural_mod_masses!(mod_masses, structural_mods, mod_to_mass)
        
        @test mod_masses[1] ≈ 79.966331f0  # Phospho on P
        @test mod_masses[7] ≈ 42.010565f0  # Acetyl on E
        @test all(iszero.(mod_masses[2:6]))  # No mods on other positions
        
        # Test with empty mods
        mod_masses = zeros(Float32, length(sequence))
        get_structural_mod_masses!(mod_masses, "", mod_to_mass)
        @test all(iszero.(mod_masses))  # All zeros
        
        # Test with unknown mod
        mod_masses = zeros(Float32, length(sequence))
        structural_mods = "(3,P,UnknownMod)(7,E,Acetyl)"
        get_structural_mod_masses!(mod_masses, structural_mods, mod_to_mass)
        @test mod_masses[3] ≈ 0.0f0  # Unknown mod
        @test mod_masses[7] ≈ 42.010565f0  # Acetyl on E
    end
    
    @testset "get_sulfur_counts!" begin
        sequence = "ACMPSTWY"  # C and M contain sulfur
        structural_mods = "(1,A,SulfurMod)(5,S,NonSulfurMod)"
        mods_to_sulfur_diff = Dict{String, Int8}(
            "SulfurMod" => Int8(1),  # Adds 1 sulfur
            "NonSulfurMod" => Int8(0)  # Adds 0 sulfur
        )
        
        sulfur_counts = zeros(Int8, length(sequence))
        get_sulfur_counts!(sulfur_counts, sequence, structural_mods, mods_to_sulfur_diff)
        
        @test sulfur_counts[1] == 1  # A (0) + SulfurMod (1)
        @test sulfur_counts[2] == 1  # C (1) + no mod
        @test sulfur_counts[3] == 1  # M (1) + no mod
        @test sulfur_counts[4] == 0  # P (0) + no mod
        @test sulfur_counts[5] == 0  # S (0) + NonSulfurMod (0)
        
        # Test with empty modifications
        sulfur_counts = zeros(Int8, length(sequence))
        get_sulfur_counts!(sulfur_counts, sequence, "", mods_to_sulfur_diff)
        @test sulfur_counts[2] == 1  # C (1)
        @test sulfur_counts[3] == 1  # M (1)
        @test all(iszero.(sulfur_counts[[1,4,5,6,7,8]]))  # No sulfur in other AAs
    end

    @testset "get_fragment_mz" begin
        # Mock the constants that would be in Pioneer.jl
        H2O = 18.010565f0
        PROTON = 1.007276f0
        
        # Setup for test
        aa_masses = Float32[71.03711, 97.05276, 128.09496]  # A, P, K
        structural_mod_masses = Float32[79.966331, 0.0, 42.010565]  # Phospho on A, nothing on P, Acetyl on K
        iso_mod_masses = Float32[0.0, 4.0, 0.0]  # Nothing on A, deuterium on P, nothing on K
        
        # Test b-ion calculation
        b_frag_mz = get_fragment_mz(1, 2, 'b', UInt8(1), aa_masses, structural_mod_masses, iso_mod_masses)
        # Expected: (A + Phospho + P + deuterium + proton) / charge
        expected_b_mass = (71.03711 + 79.966331 + 97.05276 + 4.0 + PROTON) / 1
        @test b_frag_mz ≈ expected_b_mass
        
        # Test y-ion calculation
        y_frag_mz = get_fragment_mz(2, 3, 'y', UInt8(2), aa_masses, structural_mod_masses, iso_mod_masses)
        # Expected: (P + deuterium + K + Acetyl + water + 2 protons) / 2
        expected_y_mass = (97.05276 + 4.0 + 128.09496 + 42.010565 + H2O + 2*PROTON) / 2
        @test y_frag_mz ≈ expected_y_mass
    end
    
    @testset "get_precursor_mz" begin
        # Mock the constants 
        H2O = 18.010565f0
        PROTON = 1.007276f0
        
        # Setup test data
        aa_masses = Float32[71.03711, 97.05276, 128.09496]  # A, P, K
        structural_mod_masses = Float32[79.966331, 0.0, 42.010565]  # Phospho on A, nothing on P, Acetyl on K
        iso_mod_masses = Float32[0.0, 4.0, 0.0]  # Nothing on A, deuterium on P, nothing on K
        
        # Test precursor m/z calculation with charge 1
        prec_mz = get_precursor_mz(3, UInt8(1), aa_masses, structural_mod_masses, iso_mod_masses)
        # Expected: (A + Phospho + P + deuterium + K + Acetyl + water + proton) / 1
        expected_prec_mass = (71.03711 + 79.966331 + 97.05276 + 4.0 + 128.09496 + 42.010565 + H2O + PROTON) / 1
        @test prec_mz ≈ expected_prec_mass
        
        # Test with charge 2
        prec_mz = get_precursor_mz(3, UInt8(2), aa_masses, structural_mod_masses, iso_mod_masses)
        # Expected: (A + Phospho + P + deuterium + K + Acetyl + water + 2*proton) / 2
        expected_prec_mass = (71.03711 + 79.966331 + 97.05276 + 4.0 + 128.09496 + 42.010565 + H2O + 2*PROTON) / 2
        @test prec_mz ≈ expected_prec_mass
    end
    
    @testset "adjust_masses!" begin
        # Mock DataFrame for testing
        using DataFrames
        
        # Mock the constants
        H2O = 18.010565f0
        PROTON = 1.007276f0
        
        # Mock fragment indices function
        global function get_fragment_indices(base_type::Char, frag_series_number::Integer, seq_length::Integer)
            if base_type == 'b'
                return 1, frag_series_number
            else # 'y' type
                return seq_length - frag_series_number + 1, seq_length
            end
        end
        
        # Create test data
        df = DataFrame(
            precursor_idx = [1, 1, 2, 2],
            sequence = ["APK", "APK", "LMR", "LMR"],
            structural_mods = ["(1,A,tag6)", "(1,A,tag6)", "(2,M,tag5)", "(2,M,tag5)"],
            isotopic_mods = ["(1, d0)", "(1, d0)", "(1, d4)", "(1, d4)"],
            prec_mz = [150.0f0, 150.0f0, 200.0f0, 200.0f0],
            frag_mz = [80.0f0, 120.0f0, 100.0f0, 150.0f0],
            prec_charge = [UInt8(1), UInt8(1), UInt8(2), UInt8(2)],
            frag_charge = [UInt8(1), UInt8(1), UInt8(1), UInt8(1)],
            frag_type = ["b2", "y2", "b2", "y2"],
            frag_series_number = [UInt8(2), UInt8(2), UInt8(2), UInt8(2)]
        )
        
        # Save original values to compare after adjustment
        original_prec_mz = copy(df.prec_mz)
        original_frag_mz = copy(df.frag_mz)
        
        # Run the function
        adjust_masses!(df)
        
        # Check that values were modified
        @test any(df.prec_mz .!= original_prec_mz)
        @test any(df.frag_mz .!= original_frag_mz)
        
        # Check that all rows for the same precursor have consistent adjustments
        @test df.prec_mz[1] == df.prec_mz[2]  # Same precursor should have same adjusted mass
        @test df.prec_mz[3] == df.prec_mz[4]  # Same precursor should have same adjusted mass
        
        # Different fragment types should have different adjustments
        @test df.frag_mz[1] != df.frag_mz[2]  # b2 vs y2 for same precursor
    end
    
    @testset "calculate_mz_and_sulfur_count!" begin
        # Mock DataFrame for testing
        using DataFrames
        
        # Mock the constants
        H2O = 18.010565f0
        PROTON = 1.007276f0
        
        # Mock fragment indices function again if needed
        if !@isdefined(get_fragment_indices)
            global function get_fragment_indices(base_type::Char, frag_series_number::Integer, seq_length::Integer)
                if base_type == 'b'
                    return 1, frag_series_number
                else # 'y' type
                    return seq_length - frag_series_number + 1, seq_length
                end
            end
        end
        
        # Setup dictionaries
        structural_mod_to_mass = Dict{String, Float32}(
            "Phospho" => 79.966331f0,
            "Acetyl" => 42.010565f0,
            "tag6" => 229.163f0
        )
        
        iso_mods_dict = Dict{String, Dict{String, Float32}}(
            "tag6" => Dict("d0" => 0.0f0, "d4" => 4.0f0, "d8" => 8.0f0),
            "tag5" => Dict("d0" => 0.0f0, "d4" => 4.0f0, "d8" => 8.0f0)
        )
        
        mods_to_sulfur_diff = Dict{String, Int8}(
            "Phospho" => Int8(0),
            "Acetyl" => Int8(0),
            "tag6" => Int8(1)  # Adds 1 sulfur
        )
        
        # Create test data with sequences containing C and M (sulfur-containing amino acids)
        df = DataFrame(
            sequence = ["ACMPK", "APKMC"],
            structural_mods = ["(1,A,tag6)", "(3,K,tag6)"],
            isotopic_mods = ["(1, d0)", "(1, d4)"],
            prec_mz = [300.0f0, 300.0f0],
            frag_mz = [120.0f0, 150.0f0],
            prec_charge = [UInt8(2), UInt8(2)],
            frag_charge = [UInt8(1), UInt8(1)],
            frag_type = ["b3", "y2"],
            frag_series_number = [UInt8(3), UInt8(2)],
            prec_sulfur_count = [Int8(0), Int8(0)],
            frag_sulfur_count = [Int8(0), Int8(0)]
        )
        
        # Run the function
        calculate_mz_and_sulfur_count!(df, structural_mod_to_mass, iso_mods_dict, mods_to_sulfur_diff)
        
        # Check sulfur counts
        # First sequence: ACMPK with tag6 on A
        # Expected precursor sulfur: 1 (from tag6) + 1 (C) + 1 (M) = 3
        @test df.prec_sulfur_count[1] == 3
        
        # b3 fragment of ACMPK covers ACM
        # Expected fragment sulfur: 1 (from tag6) + 1 (C) + 1 (M) = 3
        @test df.frag_sulfur_count[1] == 3
        
        # Second sequence: APKMC with tag6 on K
        # Expected precursor sulfur: 1 (from tag6) + 1 (M) + 1 (C) = 3
        @test df.prec_sulfur_count[2] == 3
        
        # y2 fragment of APKMC covers MC
        # Expected fragment sulfur: 1 (M) + 1 (C) = 2 (tag6 is on K which isn't in y2)
        @test df.frag_sulfur_count[2] == 2
        
        # Check that m/z values were recalculated
        @test df.prec_mz[1] != 300.0f0
        @test df.frag_mz[1] != 120.0f0
    end
end