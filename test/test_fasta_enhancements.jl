#!/usr/bin/env julia

# Test script for FASTA input enhancements and Koina warning reduction
using Pioneer
using Test
using JSON
using DataStructures: OrderedDict

# Test the new GetBuildLibParams function with various input types
@testset "FASTA Input Enhancement Tests" begin
    
    # Create temporary test directory structure
    test_dir = mktempdir()
    
    # Create test FASTA files
    dir1 = joinpath(test_dir, "dir1")
    dir2 = joinpath(test_dir, "dir2")
    mkpath(dir1)
    mkpath(dir2)
    
    # Create sample FASTA files
    fasta1_content = ">sp|P12345|PROT1 Protein 1 OS=Homo sapiens GN=GENE1\nMKLLSSIEQACDICRLKKLKCSKEKPKCAKCLKNNWECRYSPKTKRSPLTRAHLTEVESRLERL\n"
    fasta2_content = ">sp|P67890|PROT2 Protein 2 OS=Homo sapiens GN=GENE2\nMSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFL\n"
    fasta3_content = ">sp|Q11111|PROT3 Protein 3 OS=Mus musculus GN=GENE3\nMASSHLLLVLLCLGLTWGLRASPPGGSSGSGGAPLAGSPRSLPSSPTTYLSLAPLNPKVAPGA\n"
    
    open(joinpath(dir1, "proteins1.fasta"), "w") do f
        write(f, fasta1_content)
    end
    open(joinpath(dir1, "proteins2.fasta.gz"), "w") do f
        # For testing, just create a file with .fasta.gz extension (not actually compressed)
        write(f, fasta2_content)
    end
    open(joinpath(dir2, "proteins3.fasta"), "w") do f
        write(f, fasta3_content)
    end
    single_file = joinpath(test_dir, "single.fasta")
    open(single_file, "w") do f
        write(f, fasta1_content)
    end
    
    output_dir = joinpath(test_dir, "output")
    mkpath(output_dir)
    lib_name = joinpath(output_dir, "test_lib")
    
    @testset "Single Directory Input (backward compatibility)" begin
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, dir1)
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 2
        @test length(config["fasta_names"]) == 2
        @test all(length(config[k]) == 2 for k in ["fasta_header_regex_accessions", 
                                                     "fasta_header_regex_genes",
                                                     "fasta_header_regex_proteins",
                                                     "fasta_header_regex_organisms"])
        rm(params_file)
    end
    
    @testset "Single File Input" begin
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, single_file)
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 1
        @test config["fasta_paths"][1] == single_file
        @test length(config["fasta_names"]) == 1
        @test all(length(config[k]) == 1 for k in ["fasta_header_regex_accessions", 
                                                     "fasta_header_regex_genes",
                                                     "fasta_header_regex_proteins",
                                                     "fasta_header_regex_organisms"])
        rm(params_file)
    end
    
    @testset "Multiple Directory Input" begin
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, [dir1, dir2])
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 3  # 2 from dir1, 1 from dir2
        @test length(config["fasta_names"]) == 3
        @test all(length(config[k]) == 3 for k in ["fasta_header_regex_accessions", 
                                                     "fasta_header_regex_genes",
                                                     "fasta_header_regex_proteins",
                                                     "fasta_header_regex_organisms"])
        rm(params_file)
    end
    
    @testset "Mixed Directory and File Input" begin
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, [dir1, single_file])
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        @test length(config["fasta_paths"]) == 3  # 2 from dir1, 1 single file
        @test single_file in config["fasta_paths"]
        rm(params_file)
    end
    
    @testset "Custom Regex Codes - Single Set" begin
        custom_regex = Dict(
            "accessions" => "^>(\\S+)",
            "genes" => "GN=(\\S+)",
            "proteins" => "\\s+(.+?)\\s+OS=",
            "organisms" => "OS=(.+?)\\s+GN="
        )
        
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, [dir1, dir2];
                                               regex_codes = custom_regex)
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        # All files should have the same custom regex
        @test all(config["fasta_header_regex_accessions"][i] == custom_regex["accessions"] 
                 for i in 1:length(config["fasta_paths"]))
        rm(params_file)
    end
    
    @testset "Custom Regex Codes - Positional Mapping" begin
        regex_set1 = Dict(
            "accessions" => "^>(\\S+)_set1",
            "genes" => "GN=(\\S+)_set1",
            "proteins" => "\\s+(.+?)\\s+OS=_set1",
            "organisms" => "OS=(.+?)\\s+GN=_set1"
        )
        regex_set2 = Dict(
            "accessions" => "^>(\\S+)_set2",
            "genes" => "GN=(\\S+)_set2",
            "proteins" => "\\s+(.+?)\\s+OS=_set2",
            "organisms" => "OS=(.+?)\\s+GN=_set2"
        )
        
        params_file = Pioneer.GetBuildLibParams(output_dir, lib_name, [dir1, dir2];
                                               regex_codes = [regex_set1, regex_set2])
        @test isfile(params_file)
        
        config = JSON.parsefile(params_file, dicttype=OrderedDict)
        # First 2 files (from dir1) should have regex_set1
        @test config["fasta_header_regex_accessions"][1] == regex_set1["accessions"]
        @test config["fasta_header_regex_accessions"][2] == regex_set1["accessions"]
        # Last file (from dir2) should have regex_set2
        @test config["fasta_header_regex_accessions"][3] == regex_set2["accessions"]
        rm(params_file)
    end
    
    @testset "Error Cases" begin
        # Non-existent path
        @test_throws ErrorException Pioneer.GetBuildLibParams(output_dir, lib_name, "/nonexistent/path")
        
        # Non-FASTA file
        non_fasta = joinpath(test_dir, "not_fasta.txt")
        open(non_fasta, "w") do f
            write(f, "This is not a FASTA file")
        end
        @test_throws ErrorException Pioneer.GetBuildLibParams(output_dir, lib_name, non_fasta)
        
        # Mismatched regex codes
        @test_throws ErrorException Pioneer.GetBuildLibParams(output_dir, lib_name, [dir1, dir2];
                                                             regex_codes = [Dict(), Dict(), Dict()])
    end
    
    # Clean up
    rm(test_dir, recursive=true)
end

println("All FASTA enhancement tests passed!")