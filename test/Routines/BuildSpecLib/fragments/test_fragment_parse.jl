
@testset "Fragment Parse Tests" begin
    
    @testset "parseInternalIon" begin
        @testset "Basic parsing" begin
            @test parseInternalIon("\"Int/PDTAENLPAmY-CO/7\"") == "-CO\""
            @test parseInternalIon("\"Int/PDTAENLPAmY-CO^2+i/7\"") == "-CO^2+i\""
            @test parseInternalIon("\"Int/PEPTIDE-H2O-NH3/4\"") == "-H2O-NH3\""
        end
        
        @testset "No modifications" begin
            @test parseInternalIon("\"Int/SEQUENCE/3\"") == ""
            @test parseInternalIon("\"Int/PEPTIDE/7\"") == ""
        end
        
        @testset "Complex patterns" begin
            @test parseInternalIon("\"Int/ABCDEFG+NH3/2\"") == "+NH3\""
            @test parseInternalIon("\"Int/SEQUENCE-2H2O^3+2i/5\"") == "-2H2O^3+2i\""
            @test parseInternalIon("\"Int/PEP+CO-H2O^2/1\"") == "+CO-H2O^2\""
        end
    end
    
    @testset "getIonAnnotationSet" begin
        @testset "Regular and internal ions" begin
            frag_names = [
                "b9-H2O^2+3i",
                "Int/PDTAENLPAmY-CO/7",
                "Int/PEPTIDE-CO/5",  # Should merge to same as above
                "y10+NH3"
            ]
            
            ion_set = getIonAnnotationSet(frag_names)
            
            @test length(ion_set) == 3
            @test "b9-H2O^2+3i" in ion_set
            @test "Int/-CO\"" in ion_set
            @test "y10+NH3" in ion_set
        end
    end
    
    @testset "getIonAnnotationDict" begin
        test_set = Set(["b9^2", "y8-H2O", "Int/-CO^2+i", "a3", "c5^3"])
        
        annotation_dict = getIonAnnotationDict(test_set)
        
        @test length(annotation_dict) == 5
        @test all(v isa UInt32 for v in values(annotation_dict))
        
        # Check ordering (alphabetical)
        sorted_annotations = sort(collect(test_set))
        for (i, ann) in enumerate(sorted_annotations)
            @test annotation_dict[ann] == UInt32(i)
        end
    end
    
    @testset "countSulfurLoss" begin
        @testset "Single sulfur losses" begin
            @test countSulfurLoss("-CH3SOH") == -1
            @test countSulfurLoss("+CH3SOH") == 1
            @test countSulfurLoss("-SO2") == -1
            @test countSulfurLoss("+SO3") == 1
        end
        
        @testset "Multiple sulfur losses" begin
            @test countSulfurLoss("-2CH3SOH") == -2
            @test countSulfurLoss("+3SO2") == 3
            @test countSulfurLoss("-S2O3") == -2
        end
        
        @testset "No sulfur" begin
            @test countSulfurLoss("-H2O") == 0
            @test countSulfurLoss("+NH3") == 0
            @test countSulfurLoss("-CO") == 0
        end
        
        @testset "Complex formulas" begin
            @test countSulfurLoss("-CH3S2OH") == -2
            @test countSulfurLoss("+2S3O4") == 6  # 2 * 3 sulfurs
        end
    end
    
    @testset "count_sulfurs!" begin
        seq_idx_to_sulfur = zeros(UInt8, 20)
        
        @testset "Basic sequence counting" begin
            sequence = "ACDEFMHIKLMNPQRSTVWY"
            mods_iterator = Base.RegexMatchIterator(r"(?<=\().*?(?=\))", "", false)  # No mods
            mods_to_sulfur = Dict{String, Int8}()
            
            sulfur_count = count_sulfurs!(seq_idx_to_sulfur, sequence, mods_iterator, mods_to_sulfur)
            
            @test sulfur_count == 3  # C at position 2, M at position 6, M at position 12
            @test seq_idx_to_sulfur[2] == 1  # C
            @test seq_idx_to_sulfur[6] == 1  # M
            @test seq_idx_to_sulfur[11] == 1  # M
            @test sum(seq_idx_to_sulfur) == 3
        end
        
        @testset "With modifications" begin
            fill!(seq_idx_to_sulfur, 0)
            sequence = "PEPTIDE"
            # Mock modification at position 3 that adds 2 sulfurs
            mods_str = "1(3,P,TestMod)"
            mods_iterator = Base.RegexMatchIterator(r"(?<=\().*?(?=\))", mods_str, false)
            mods_to_sulfur = Dict("TestMod" => Int8(2))
            
            sulfur_count = count_sulfurs!(seq_idx_to_sulfur, sequence, mods_iterator, mods_to_sulfur)
            
            @test sulfur_count == 2
            @test seq_idx_to_sulfur[3] == 2  # Modified position
        end
        
        @testset "Multiple cysteines and methionines" begin
            fill!(seq_idx_to_sulfur, 0)
            sequence = "MCMCMCMC"
            mods_iterator = Base.RegexMatchIterator(r"(?<=\().*?(?=\))", "", false)
            mods_to_sulfur = Dict{String, Int8}()
            
            sulfur_count = count_sulfurs!(seq_idx_to_sulfur, sequence, mods_iterator, mods_to_sulfur)
            
            @test sulfur_count == 8  # All positions have sulfur
            @test all(seq_idx_to_sulfur[1:8] .== 1)
        end
    end
    
    @testset "fill_isotope_mods!" begin
        seq_idx_to_iso_mod = zeros(Float32, 20)
        
        @testset "Single isotope modification" begin
            # Mock modification string: position 5 with mass 2.0
            mods_str = "1(5,K,Heavy)"
            mods_iterator = Base.RegexMatchIterator(r"(?<=\().*?(?=\))", mods_str, false)
            iso_mod_to_mass = Dict("Heavy" => Float32(2.0))
            
            fill_isotope_mods!(seq_idx_to_iso_mod, mods_iterator, iso_mod_to_mass)
            
            @test seq_idx_to_iso_mod[5] == 2.0
            @test sum(seq_idx_to_iso_mod) == 2.0
        end
        
        @testset "Multiple isotope modifications" begin
            fill!(seq_idx_to_iso_mod, 0)
            mods_str = "3(2,R,Heavy1)(7,K,Heavy2)(10,R,Heavy1)"
            mods_iterator = Base.RegexMatchIterator(r"(?<=\().*?(?=\))", mods_str, false)
            iso_mod_to_mass = Dict("Heavy1" => Float32(1.5), "Heavy2" => Float32(3.0))
            
            fill_isotope_mods!(seq_idx_to_iso_mod, mods_iterator, iso_mod_to_mass)
            
            @test seq_idx_to_iso_mod[2] == 1.5
            @test seq_idx_to_iso_mod[7] == 3.0
            @test seq_idx_to_iso_mod[10] == 1.5
            @test sum(seq_idx_to_iso_mod) ≈ 6.0
        end
    end
    
    @testset "get_fragment_indices" begin
        sequence_length = UInt8(20)
        
        @testset "N-terminal ions" begin
            @test get_fragment_indices('b', UInt8(5), sequence_length) == (UInt8(1), UInt8(5))
            @test get_fragment_indices('a', UInt8(10), sequence_length) == (UInt8(1), UInt8(10))
            @test get_fragment_indices('c', UInt8(3), sequence_length) == (UInt8(1), UInt8(3))
        end
        
        @testset "C-terminal ions" begin
            @test get_fragment_indices('y', UInt8(5), sequence_length) == (UInt8(16), UInt8(20))
            @test get_fragment_indices('x', UInt8(10), sequence_length) == (UInt8(11), UInt8(20))
            @test get_fragment_indices('z', UInt8(3), sequence_length) == (UInt8(18), UInt8(20))
        end
        
        @testset "Precursor ion" begin
            @test get_fragment_indices('p', UInt8(0), sequence_length) == (UInt8(1), UInt8(20))
            @test get_fragment_indices('p', UInt8(99), sequence_length) == (UInt8(1), UInt8(20))
        end
    end
    
    @testset "apply_isotope_mod" begin
        seq_idx_to_iso_mod = Float32[0, 0, 2.0, 0, 1.5, 0, 0, 3.0, 0, 0]
        
        @testset "Single charge" begin
            # Positions 3-8, charge 1
            mass_shift = apply_isotope_mod(UInt8(1), seq_idx_to_iso_mod, UInt8(3), UInt8(8))
            @test mass_shift ≈ 6.5  # 2.0 + 1.5 + 3.0
        end
        
        @testset "Multiple charges" begin
            # Same range, charge 2
            mass_shift = apply_isotope_mod(UInt8(2), seq_idx_to_iso_mod, UInt8(3), UInt8(8))
            @test mass_shift ≈ 3.25  # 6.5 / 2
            
            # Charge 3
            mass_shift = apply_isotope_mod(UInt8(3), seq_idx_to_iso_mod, UInt8(3), UInt8(8))
            @test mass_shift ≈ 6.5 / 3
        end
        
        @testset "No modifications in range" begin
            mass_shift = apply_isotope_mod(UInt8(1), seq_idx_to_iso_mod, UInt8(1), UInt8(2))
            @test mass_shift == 0.0
        end
    end
    
    @testset "getNumeric" begin
        @testset "Basic extraction" begin
            @test getNumeric("b11", UInt8(0)) == UInt8(11)
            @test getNumeric("y123", UInt8(0)) == UInt8(123)
            @test getNumeric("x5y", UInt8(0)) == UInt8(5)
        end
        
        @testset "Multiple numbers" begin
            @test getNumeric("x101x102", UInt8(0)) == UInt8(101)  # First match
            @test getNumeric("123abc456", UInt8(0)) == UInt8(123)
        end
        
        @testset "No numeric content" begin
            @test getNumeric("xyz", UInt8(0)) == UInt8(0)
            @test getNumeric("abc", UInt8(1)) == UInt8(1)
            @test getNumeric("", UInt8(99)) == UInt8(99)
        end
        
        @testset "Different numeric types" begin
            @test getNumeric("42", Float32(0)) == Float32(42)
            @test getNumeric("3.14", Int32(0)) == Int32(3)  # Parse as int
        end
    end
    
    @testset "get_immonium_sulfur_dict" begin
        temp_dir = mktempdir()
        immonium_path = joinpath(temp_dir, "immonium.txt")
        
        # Create test immonium file
        open(immonium_path, "w") do f
            println(f, "IA\tC2H3N")
            println(f, "IC\tC2H3NS")  # Contains 1 sulfur
            println(f, "IM\tC4H7NS")  # Contains 1 sulfur
            println(f, "IW\tC9H6N2")
            println(f, "ICCAM\tC5H8N2OS2")  # Contains 2 sulfurs
        end
        
        try
            sulfur_dict = get_immonium_sulfur_dict(immonium_path)
            
            @test sulfur_dict["IA"] == 0
            @test sulfur_dict["IC"] == 1
            @test sulfur_dict["IM"] == 1
            @test sulfur_dict["IW"] == 0
            @test sulfur_dict["ICCAM"] == 2
        finally
            rm(temp_dir, recursive=true)
        end
    end
    
    @testset "Integration test: parse_koina_fragments" begin
        temp_dir = mktempdir()
        
        # Create mock precursor table
        precursor_data = DataFrame(
            sequence = ["PEPTIDE", "SEQUENCE", "FRAGMENT"],
            mods = ["", "1(3,Q,Oxidation)", ""],
            isotope_mods = ["", "1(5,E,Heavy)", ""],
            charge = [2, 3, 2]
        )
        
        # Create mock fragment table
        fragment_data = DataFrame(
            precursor_idx = [1, 1, 1, 2, 2, 3, 3, 3],
            annotation = ["b3", "y4^2", "b5-H2O", "y7^3+i", "b2", "a3", "y5-NH3", "p-H2O"],
            mz = Float32[300.1, 400.2, 500.3, 600.4, 200.1, 350.2, 450.3, 550.4],
            intensities = Float32[0.8, 0.9, 0.7, 0.6, 0.95, 0.5, 0.85, 0.4]
        )
        
        # Save as Arrow files
        precursor_path = joinpath(temp_dir, "precursors.arrow")
        fragment_path = joinpath(temp_dir, "fragments.arrow")
        Arrow.write(precursor_path, precursor_data)
        Arrow.write(fragment_path, fragment_data)
        
        # Create test ion annotation set
        ion_set = Set(["b3", "y4^2", "b5-H2O", "y7^3+i", "b2", "a3", "y5-NH3", "p-H2O"])
        frag_name_to_idx = Dict(ann => UInt16(i) for (i, ann) in enumerate(sort(collect(ion_set))))
        
        # Mock parameters
        precursor_batch_size = 2
        immonium_data_path = joinpath(temp_dir, "immonium.txt")
        
        # Create empty immonium file
        open(immonium_data_path, "w") do f
            println(f, "IA\tC2H3N")
        end
        
        mods_to_sulfur_diff = Dict{String, Int8}("Oxidation" => 0)
        iso_mod_to_mass = Dict{String, Float32}("Heavy" => 1.0)
        
        try
            # Test parsing
            precursor_table = Arrow.Table(precursor_path)
            fragment_table = Arrow.Table(fragment_path)
            
            result = parse_koina_fragments(
                precursor_table,
                fragment_table,
                UniSpecFragAnnotation(""),  # Dummy annotation type
                ion_set,
                frag_name_to_idx,
                precursor_batch_size,
                immonium_data_path,
                temp_dir,
                mods_to_sulfur_diff,
                iso_mod_to_mass,
                InstrumentSpecificModel("unispec")
            )
            
            # Verify output files were created
            @test isfile(joinpath(temp_dir, "fragments_table.arrow"))
            @test isfile(joinpath(temp_dir, "prec_to_frag.arrow"))
            @test isfile(joinpath(temp_dir, "frag_name_to_idx.jld2"))
            @test isfile(joinpath(temp_dir, "ion_annotations.jld2"))
            
            # Verify result dictionary
            @test length(result) == length(ion_set)
            @test all(v isa PioneerFragAnnotation for v in values(result))
        finally
            rm(temp_dir, recursive=true)
        end
    end
end