
@testset "Fragment Annotation Tests" begin
    
    @testset "get_ion_annotation_set" begin
        @testset "Regular ions" begin
            frag_names = ["b9-H2O^2+3i", "y17-H2O-NH3^4+i", "y16-CH3SOH^4", "b11^2"]
            ion_set = get_ion_annotation_set(frag_names)
            
            @test length(ion_set) == 4
            @test "b9-H2O^2+3i" in ion_set
            @test "y17-H2O-NH3^4+i" in ion_set
            @test "y16-CH3SOH^4" in ion_set
            @test "b11^2" in ion_set
        end
        
        @testset "Internal ions" begin
            # Internal ions should be processed to remove sequence
            frag_names = [
                "\"Int/PDTAENLPAmY-CO/7\"",
                "\"Int/PDTAENLPAmY-CO^2+i/7\"",
                "\"Int/PEPTIDE-H2O^2+i/4\"",
                "\"Int/SEQUENCE/3\""  # No modifications
            ]
            ion_set = get_ion_annotation_set(frag_names)
            
            @test length(ion_set) == 4
            @test "Int/-CO\"" in ion_set
            @test "Int/-CO^2+i\"" in ion_set
            @test "Int/-H2O^2+i\"" in ion_set
            @test "Int/\"" in ion_set  # Empty when no modifications
        end
        
        @testset "Mixed ion types" begin
            frag_names = [
                "b9-H2O^2+3i",
                "\"Int/PDTAENLPAmY-CO/7\"",
                "IEA",  # Immonium
                "p-2CH3SOH-H2O"  # Precursor
            ]
            ion_set = get_ion_annotation_set(frag_names)
            
            @test length(ion_set) == 4
            @test "b9-H2O^2+3i" in ion_set
            @test "Int/-CO\"" in ion_set
            @test "IEA" in ion_set
            @test "p-2CH3SOH-H2O" in ion_set
        end
        
        @testset "Duplicate handling" begin
            # Same annotations should not create duplicates
            frag_names = [
                "b9-H2O^2+3i",
                "b9-H2O^2+3i",
                "\"Int/PEPTIDE1-CO/7\"",
                "\"Int/PEPTIDE2-CO/5\""  # Different sequence, same modification
            ]
            ion_set = get_ion_annotation_set(frag_names)
            
            @test length(ion_set) == 2
            @test "b9-H2O^2+3i" in ion_set
            @test "Int/-CO\"" in ion_set
        end
    end
    
    @testset "parse_internal_ion" begin
        @testset "Basic internal ions" begin
            @test parse_internal_ion("\"Int/PDTAENLPAmY-CO/7\"") == "-CO\""
            @test parse_internal_ion("\"Int/PDTAENLPAmY-CO^2+i/7\"") == "-CO^2+i\""
            @test parse_internal_ion("\"Int/PEPTIDE-H2O-NH3/4\"") == "-H2O-NH3\""
        end
        
        @testset "Charge and isotope states" begin
            @test parse_internal_ion("\"Int/SEQUENCE^2/3\"") == "^2\""
            @test parse_internal_ion("\"Int/SEQUENCE+i/3\"") == "+i\""
            @test parse_internal_ion("\"Int/SEQUENCE^3+2i/3\"") == "^3+2i\""
        end
        
        @testset "No modifications" begin
            @test parse_internal_ion("\"Int/SEQUENCE/3\"") == "\""
            @test parse_internal_ion("\"Int/PEPTIDE/7\"") == "\""
        end
        
        @testset "Complex modifications" begin
            @test parse_internal_ion("\"Int/PEPTIDE-2H2O+NH3^2+3i/5\"") == "-2H2O+NH3^2+3i\""
            @test parse_internal_ion("\"Int/SEQUENCE+CH3-CO^3/2\"") == "+CH3-CO^3\""
        end
    end
    
    @testset "create_ion_annotation_index" begin
        @testset "Basic indexing" begin
            annotations = Set(["b9^2", "y8-H2O", "Int/-CO^2+i", "IEA"])
            index = create_ion_annotation_index(annotations)
            
            @test length(index) == 4
            @test all(v isa UInt16 for v in values(index))
            @test all(k in annotations for k in keys(index))
            
            # Check that indices are unique
            @test length(unique(values(index))) == 4
            
            # Check that indices start from 1
            @test minimum(values(index)) == 1
            @test maximum(values(index)) == 4
        end
        
        @testset "Deterministic ordering" begin
            # Same input should produce same indices
            annotations = Set(["z5", "a3", "b2", "y7"])
            index1 = create_ion_annotation_index(annotations)
            index2 = create_ion_annotation_index(annotations)
            
            @test index1 == index2
            
            # Sorted order should be maintained
            sorted_annotations = sort(collect(annotations))
            for (i, ann) in enumerate(sorted_annotations)
                @test index1[ann] == i
            end
        end
        
        @testset "Empty set" begin
            annotations = Set{String}()
            index = create_ion_annotation_index(annotations)
            
            @test isempty(index)
        end
        
        @testset "Large annotation set" begin
            # Test with many annotations
            annotations = Set{String}()
            for ion_type in ["b", "y"]
                for i in 1:20
                    for charge in 1:3
                        push!(annotations, "$(ion_type)$(i)^$(charge)")
                    end
                end
            end
            
            index = create_ion_annotation_index(annotations)
            
            @test length(index) == length(annotations)
            @test all(v <= typemax(UInt16) for v in values(index))
        end
    end
    
    @testset "parse_fragment_annotation with UniSpecFragAnnotation" begin
        # Create mock dictionaries
        immonium_to_sulfur = Dict{String, Int8}(
            "IEA" => 0,
            "ICM" => 1,  # Cysteine/Methionine immonium
            "ICCAM" => 2  # Modified cysteine
        )
        
        @testset "Basic b/y ions" begin
            # Mock annotation iterator for "b11-2H2O^3+2i"
            annotation = UniSpecFragAnnotation("b11-2H2O^3+2i")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.base_type == 'b'
            @test result.frag_index == 11
            @test result.charge == 3
            @test result.isotope == 2
            @test result.internal == false
            @test result.immonium == false
            @test result.is_gain_loss == true
            @test result.sulfur_diff == 0  # H2O has no sulfur
        end
        
        @testset "Sulfur-containing losses" begin
            # Test with sulfur-containing neutral loss
            annotation = UniSpecFragAnnotation("y16-CH3SOH^2")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.base_type == 'y'
            @test result.frag_index == 16
            @test result.charge == 2
            @test result.isotope == 0
            @test result.sulfur_diff == -1  # Lost one sulfur
        end
        
        @testset "Multiple sulfur losses" begin
            annotation = UniSpecFragAnnotation("b9-2CH3SOH+i")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.sulfur_diff == -2  # Lost two sulfurs
            @test result.isotope == 1
        end
        
        @testset "Immonium ions" begin
            annotation = UniSpecFragAnnotation("ICM")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.base_type == '_'
            @test result.frag_index == 1  # Default
            @test result.immonium == true
            @test result.sulfur_diff == 1  # From dictionary
        end
        
        @testset "Internal ions" begin
            annotation = UniSpecFragAnnotation("Int/PEPTIDE-CO^2+i/7")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.internal == true
            @test result.charge == 2
            @test result.isotope == 1
            @test result.is_gain_loss == true
        end
        
        @testset "Precursor ions" begin
            annotation = UniSpecFragAnnotation("p-H2O-NH3^3")
            result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
            
            @test result.base_type == 'p'
            @test result.frag_index == 0  # Precursor has no index
            @test result.charge == 3
            @test result.is_gain_loss == true
        end
        
        @testset "Other ion types" begin
            for (ion_type, expected_char) in [("a5", 'a'), ("c12^2", 'c'), ("x3+i", 'x'), ("z7-NH3", 'z')]
                annotation = UniSpecFragAnnotation(ion_type)
                result = parse_fragment_annotation(annotation; immonium_to_sulfur_count=immonium_to_sulfur)
                @test result.base_type == expected_char
            end
        end
    end
    
    @testset "getAnnotationToID" begin
        # Create test dictionary
        test_dict = Dict{String, PioneerFragAnnotation}(
            "\"p-2CH3SOH-H2O\"" => PioneerFragAnnotation('p', 0x00, 0x01, 0x00, false, false, true, -2),
            "\"IEA\"" => PioneerFragAnnotation('_', 0x00, 0x01, 0x00, false, true, false, 0),
            "\"y1-NH3\"" => PioneerFragAnnotation('y', 0x01, 0x01, 0x00, false, false, true, 0),
            "\"Int/-CO^2+i\"" => PioneerFragAnnotation('_', 0x00, 0x02, 0x01, true, false, true, 0)
        )
        
        id_to_annotation, annotation_to_id = getAnnotationToID(test_dict)
        
        @testset "Bidirectional mapping" begin
            @test length(id_to_annotation) == 4
            @test length(annotation_to_id) == 4
            
            # Check that all annotations are mapped
            for (ann, frag_data) in test_dict
                @test haskey(annotation_to_id, ann)
                id = annotation_to_id[ann]
                @test id_to_annotation[id] == frag_data
            end
        end
        
        @testset "ID range" begin
            # IDs should be 1-indexed UInt16
            ids = collect(values(annotation_to_id))
            @test minimum(ids) == 1
            @test maximum(ids) == 4
            @test all(id isa UInt16 for id in ids)
        end
    end
    
    @testset "get_altimeter_ion_dict" begin
        # Create temporary test file
        temp_dir = mktempdir()
        ion_table_path = joinpath(temp_dir, "ion_table.tsv")
        
        # Write test data
        open(ion_table_path, "w") do f
            println(f, "b1\t100.0\t1\tinfo1")
            println(f, "y1\t200.0\t1\tinfo2")
            println(f, "b2^2\t150.0\t2\tinfo3")
            println(f, "y3-H2O\t300.0\t1\tinfo4")
        end
        
        try
            ion_dict = get_altimeter_ion_dict(ion_table_path)
            
            @test length(ion_dict) == 4
            @test ion_dict[Int32(0)] == "b1"
            @test ion_dict[Int32(1)] == "y1"
            @test ion_dict[Int32(2)] == "b2^2"
            @test ion_dict[Int32(3)] == "y3-H2O"
            
            # Test that indices are 0-based
            @test haskey(ion_dict, Int32(0))
            @test !haskey(ion_dict, Int32(4))
        finally
            rm(temp_dir, recursive=true)
        end
    end
end