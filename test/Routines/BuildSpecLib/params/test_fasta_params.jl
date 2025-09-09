# Test for FASTA input enhancement in GetBuildLibParams
using Test
using Pioneer
using JSON
using DataStructures: OrderedDict
using CodecZlib

@testset "GetBuildLibParams FASTA Input Enhancement" begin
    
    # Setup test data directory
    test_data_dir = joinpath(dirname(@__DIR__), "..", "..", "..", "data", "test_fasta_params")
    
    # Clean up existing parameter JSON files if directory exists
    if isdir(test_data_dir)
        for item in readdir(test_data_dir, join=true)
            if endswith(item, ".json") || isdir(item)
                rm(item, force=true, recursive=true)
            end
        end
    else
        mkpath(test_data_dir)
    end
    
    # Setup FASTA test files directory (permanent, not temp)
    fasta_test_dir = joinpath(dirname(@__DIR__), "..", "..", "..", "data", "test_fasta_files")
    
    # Create directory if it doesn't exist, but don't delete existing files
    # This preserves files like minimal_protein.fasta used by other tests
    if !isdir(fasta_test_dir)
        mkpath(fasta_test_dir)
    end
    
    # Create test directory structure
    dir1 = joinpath(fasta_test_dir, "fastas_dir1")
    dir2 = joinpath(fasta_test_dir, "fastas_dir2")
    mkpath(dir1)
    mkpath(dir2)
    
    # Create sample FASTA content
    fasta1_content = ">sp|P12345|PROT1 Protein 1 OS=Homo sapiens GN=GENE1\nMKLLSSIEQACDICRLKKLKCSKEKPKCAKCLKNNWECRYSPKTKRSPLTRAHLTEVESRLERL\n"
    fasta2_content = ">sp|P67890|PROT2 Protein 2 OS=Homo sapiens GN=GENE2\nMSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFL\n"
    fasta3_content = ">sp|Q11111|PROT3 Protein 3 OS=Mus musculus GN=GENE3\nMASSHLLLVLLCLGLTWGLRASPPGGSSGSGGAPLAGSPRSLPSSPTTYLSLAPLNPKVAPGA\n"
    custom_fasta_content = ">CustomID_001 Custom protein description gene=CUSTOM1\nMEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQNSY\n"
    
    # Create minimal_protein.fasta if it doesn't exist (used by BuildSpecLib tests)
    minimal_protein_path = joinpath(fasta_test_dir, "minimal_protein.fasta")
    if !isfile(minimal_protein_path)
        open(minimal_protein_path, "w") do f
            write(f, ">TEST_PROTEIN Test protein for missed cleavage testing\nMARALYKDERPK")
        end
    end
    
    # Write FASTA files only if they don't exist
    proteins1_path = joinpath(dir1, "proteins1.fasta")
    if !isfile(proteins1_path)
        open(proteins1_path, "w") do f
            write(f, fasta1_content)
        end
    end
    
    proteins2_path = joinpath(dir1, "proteins2.fasta.gz")
    if !isfile(proteins2_path)
        # Actually compress the file since it has .gz extension
        open(GzipCompressorStream, proteins2_path, "w") do f
            write(f, fasta2_content)
        end
    end
    
    proteins3_path = joinpath(dir2, "proteins3.fasta")
    if !isfile(proteins3_path)
        open(proteins3_path, "w") do f
            write(f, fasta3_content)
        end
    end
    
    single_file = joinpath(fasta_test_dir, "single_protein.fasta")
    if !isfile(single_file)
        open(single_file, "w") do f
            write(f, custom_fasta_content)
        end
    end
    
    # Setup output paths
    output_dir = joinpath(test_data_dir, "output")
    mkpath(output_dir)
    lib_name = joinpath(output_dir, "test_lib")
    
    @testset "Single Directory Input (backward compatibility)" begin
        params_file = joinpath(test_data_dir, "single_dir_params.json")
        result = GetBuildLibParams(output_dir, lib_name, dir1; params_path=params_file)
        
        @test result == params_file
        @test isfile(params_file)
        
        # Verify JSON content
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 2
        @test length(config["fasta_names"]) == 2
        @test "PROTEINS1" in config["fasta_names"]
        @test "PROTEINS2" in config["fasta_names"]
        
        # Check regex arrays are properly expanded
        @test length(config["fasta_header_regex_accessions"]) == 2
        @test length(config["fasta_header_regex_genes"]) == 2
        @test length(config["fasta_header_regex_proteins"]) == 2
        @test length(config["fasta_header_regex_organisms"]) == 2
        
        # All should have the same default regex
        @test config["fasta_header_regex_accessions"][1] == config["fasta_header_regex_accessions"][2]
    end
    
    @testset "Single File Input" begin
        params_file = joinpath(test_data_dir, "single_file_params.json")
        result = GetBuildLibParams(output_dir, lib_name, single_file; params_path=params_file)
        
        @test result == params_file
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 1
        @test config["fasta_paths"][1] == single_file
        @test length(config["fasta_names"]) == 1
        @test config["fasta_names"][1] == "SINGLE_PROTEIN"
        
        # Check regex arrays
        @test length(config["fasta_header_regex_accessions"]) == 1
    end
    
    @testset "Multiple Directory Input" begin
        params_file = joinpath(test_data_dir, "multi_dir_params.json")
        result = GetBuildLibParams(output_dir, lib_name, [dir1, dir2]; params_path=params_file)
        
        @test result == params_file
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 3  # 2 from dir1, 1 from dir2
        @test length(config["fasta_names"]) == 3
        
        # Check all regex arrays have correct length
        @test all(length(config[k]) == 3 for k in ["fasta_header_regex_accessions", 
                                                     "fasta_header_regex_genes",
                                                     "fasta_header_regex_proteins",
                                                     "fasta_header_regex_organisms"])
    end
    
    @testset "Mixed Directory and File Input" begin
        params_file = joinpath(test_data_dir, "mixed_input_params.json")
        result = GetBuildLibParams(output_dir, lib_name, [dir1, single_file]; params_path=params_file)
        
        @test result == params_file
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 3  # 2 from dir1, 1 single file
        @test single_file in config["fasta_paths"]
        @test "SINGLE_PROTEIN" in config["fasta_names"]
    end
    
    @testset "Custom Regex Codes - Single Set" begin
        custom_regex = Dict(
            "accessions" => "^>(\\S+)",
            "genes" => "GN=(\\S+)",
            "proteins" => "\\s+(.+?)\\s+OS=",
            "organisms" => "OS=(.+?)\\s+GN="
        )
        
        params_file = joinpath(test_data_dir, "custom_regex_single_params.json")
        result = GetBuildLibParams(output_dir, lib_name, [dir1, dir2];
                                          params_path=params_file,
                                          regex_codes=custom_regex)
        
        @test result == params_file
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        
        # All files should have the same custom regex
        num_files = length(config["fasta_paths"])
        @test all(config["fasta_header_regex_accessions"][i] == custom_regex["accessions"] 
                 for i in 1:num_files)
        @test all(config["fasta_header_regex_genes"][i] == custom_regex["genes"] 
                 for i in 1:num_files)
    end
    
    @testset "Custom Regex Codes - Positional Mapping" begin
        regex_set1 = Dict(
            "accessions" => "^\\w+\\|(\\w+)\\|",
            "genes" => " GN=(\\S+)",
            "proteins" => "^\\w+\\|\\w+\\|\\S+ (.+?) OS=",
            "organisms" => " OS=(.+?)( GN=|\$)"
        )
        regex_set2 = Dict(
            "accessions" => "^>(\\S+)",
            "genes" => "gene=(\\S+)",
            "proteins" => "^>\\S+ (.+?)( gene=|\$)",
            "organisms" => "organism=(\\S+)"
        )
        
        params_file = joinpath(test_data_dir, "custom_regex_positional_params.json")
        result = GetBuildLibParams(output_dir, lib_name, [dir1, single_file];
                                          params_path=params_file,
                                          regex_codes=[regex_set1, regex_set2])
        
        @test result == params_file
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        
        # First 2 files (from dir1) should have regex_set1
        @test config["fasta_header_regex_accessions"][1] == regex_set1["accessions"]
        @test config["fasta_header_regex_accessions"][2] == regex_set1["accessions"]
        @test config["fasta_header_regex_genes"][1] == regex_set1["genes"]
        
        # Last file (single_file) should have regex_set2
        @test config["fasta_header_regex_accessions"][3] == regex_set2["accessions"]
        @test config["fasta_header_regex_genes"][3] == regex_set2["genes"]
    end
    
    @testset "Error Cases" begin
        # Non-existent path
        @test_throws ErrorException GetBuildLibParams(output_dir, lib_name, "/nonexistent/path")
        
        # Non-FASTA file
        non_fasta = joinpath(fasta_test_dir, "not_fasta.txt")
        if !isfile(non_fasta)
            open(non_fasta, "w") do f
                write(f, "This is not a FASTA file")
            end
        end
        @test_throws ErrorException GetBuildLibParams(output_dir, lib_name, non_fasta)
        
        # Mismatched regex codes (more regex sets than inputs)
        # Now issues a warning instead of error
        result = GetBuildLibParams(
            output_dir, lib_name, [dir1, dir2];
            regex_codes = [Dict(), Dict(), Dict()]
        )
        @test isfile(result)  # Should still create params file despite warning
        
        # Mismatched regex codes (more than 1 but not matching inputs)
        # Now issues a warning instead of error
        result2 = GetBuildLibParams(
            output_dir, lib_name, [dir1, dir2, single_file];
            regex_codes = [Dict(), Dict()]
        )
        @test isfile(result2)  # Should still create params file despite warning
    end
    
    @testset "Empty Directory Warning" begin
        empty_dir = joinpath(fasta_test_dir, "empty_dir")
        if !isdir(empty_dir)
            mkpath(empty_dir)
        end
        
        # Should not error but should warn and produce empty arrays
        params_file = joinpath(test_data_dir, "empty_dir_params.json")
        
        # This now issues a warning instead of error for empty directory
        result = GetBuildLibParams(
            output_dir, lib_name, empty_dir;
            params_path=params_file
        )
        @test isfile(params_file)  # Should still create params file
        
        # Verify the params file has empty fasta arrays
        config = JSON.parsefile(params_file)
        @test isempty(config["fasta_paths"])
        @test isempty(config["fasta_names"])
    end
    
    # Don't cleanup - leave test files for inspection
end

println("✓ GetBuildLibParams FASTA enhancement tests completed")