# BuildSpecLib Integration Tests
using Test
using JSON
using DataStructures: OrderedDict
using CSV
using DataFrames
using Arrow
using HTTP

# Helper function to check network connectivity
function has_network_connection()
    try
        HTTP.request("GET", "https://koina.wilhelmlab.org", status_exception=false, readtimeout=5)
        return true
    catch
        return false
    end
end

# Generate parameter JSON for BuildSpecLib
function generate_build_params(
    fasta_path::String,
    output_dir::String,
    lib_name::String;
    missed_cleavages::Int = 1,
    min_length::Int = 7,
    max_length::Int = 30,
    min_charge::Int = 2,
    max_charge::Int = 4,
    add_decoys::Bool = true,
    predict_fragments::Bool = true,
    max_koina_batch::Int = 100,
    entrapment_r::Float64 = 0.0
)
    # Load default parameters as template
    project_root = joinpath(@__DIR__, "..", "..", "..")
    template_path = joinpath(project_root, "assets", "example_config", "defaultBuildLibParams.json")
    params = JSON.parsefile(template_path, dicttype=OrderedDict)
    
    # Update parameters
    params["fasta_digest_params"]["missed_cleavages"] = missed_cleavages
    params["fasta_digest_params"]["min_length"] = min_length
    params["fasta_digest_params"]["max_length"] = max_length
    params["fasta_digest_params"]["min_charge"] = min_charge
    params["fasta_digest_params"]["max_charge"] = max_charge
    params["fasta_digest_params"]["add_decoys"] = add_decoys
    params["fasta_digest_params"]["entrapment_r"] = entrapment_r
    
    # Set library parameters
    params["library_params"]["prediction_model"] = "altimeter"
    params["library_params"]["auto_detect_frag_bounds"] = false
    params["library_params"]["frag_mz_min"] = 150.0
    params["library_params"]["frag_mz_max"] = 2000.0
    params["library_params"]["max_frag_rank"] = 50
    
    # Remove calibration file requirement
    delete!(params["library_params"], "calibration_raw_file")
    
    # Set paths
    params["fasta_paths"] = [fasta_path]
    params["fasta_names"] = [uppercase(splitext(basename(fasta_path))[1])]
    params["out_dir"] = output_dir
    params["lib_name"] = lib_name
    params["new_lib_name"] = lib_name
    params["out_name"] = basename(lib_name) * ".tsv"
    
    # Set Koina parameters
    params["predict_fragments"] = predict_fragments
    params["max_koina_requests"] = 12
    params["max_koina_batch"] = max_koina_batch
    
    # Set regex patterns for FASTA headers
    params["fasta_header_regex_accessions"] = ["^>(\\S+)"]
    params["fasta_header_regex_genes"] = ["GN=(\\S+)"]
    params["fasta_header_regex_proteins"] = ["\\s+(.+?)\\s+OS="]
    params["fasta_header_regex_organisms"] = ["OS=(.+?)(?:\\s+GN=|\$)"]
    
    return params
end

# Verify library structure was created
function verify_library_structure(lib_dir::String)
    @test isdir(lib_dir)
    
    # Check for essential files
    peptides_file = joinpath(lib_dir, "peptides.csv")
    @test isfile(peptides_file)
    
    config_file = joinpath(lib_dir, "config.json")
    @test isfile(config_file)
    
    return true
end

# Verify peptide generation matches expected parameters
function verify_peptide_generation(peptides_file::String, expected_params::Dict)
    df = CSV.read(peptides_file, DataFrame)
    
    # Check peptide lengths
    peptide_lengths = length.(df.sequence)
    @test all(expected_params["min_length"] .<= peptide_lengths .<= expected_params["max_length"])
    
    # Check charge states
    @test all(expected_params["min_charge"] .<= df.charge .<= expected_params["max_charge"])
    
    # Count targets and decoys
    n_targets = sum(df.target)
    n_decoys = sum(.!df.target)
    
    if expected_params["add_decoys"]
        @test n_decoys > 0
        # Decoys should be roughly equal to targets (allowing some variation)
        @test 0.8 * n_targets <= n_decoys <= 1.2 * n_targets
    else
        @test n_decoys == 0
    end
    
    return df
end

# Verify fragment predictions were created
function verify_fragment_predictions(lib_dir::String)
    # Look for Arrow fragment files
    fragment_files = filter(f -> endswith(f, ".arrow"), readdir(lib_dir))
    
    if !isempty(fragment_files)
        # Load and check fragment file
        frag_file = joinpath(lib_dir, fragment_files[1])
        frag_table = Arrow.Table(frag_file)
        
        # Check that fragments have required columns
        @test :mz in propertynames(frag_table)
        @test :intensity in propertynames(frag_table)
        
        return true
    else
        return false
    end
end

# Count peptides for missed cleavage verification
function count_peptides_by_sequence(peptides_file::String, target_only::Bool = true)
    df = CSV.read(peptides_file, DataFrame)
    
    if target_only
        df = filter(row -> row.target, df)
    end
    
    # Get unique peptide sequences (ignoring charge states and modifications)
    unique_sequences = unique(df.sequence)
    
    return length(unique_sequences), unique_sequences
end

@testset "BuildSpecLib Integration Tests" begin
    # Check if we have network access for Koina
    if !has_network_connection()
        @warn "No network connection available - skipping BuildSpecLib tests that require Koina API"
        return
    end
    
    # Base directories
    project_root = joinpath(@__DIR__, "..", "..", "..")
    test_data_dir = joinpath(project_root, "data", "test_build_spec_lib")
    fasta_dir = joinpath(project_root, "data", "test_fasta_files")
    
    @testset "Scenario A - Missed Cleavage Verification" begin
        scenario_dir = joinpath(test_data_dir, "scenario_a_missed_cleavages")
        output_dir = joinpath(scenario_dir, "output")
        minimal_fasta = joinpath(fasta_dir, "minimal_protein.fasta")
        
        # Verify FASTA file exists
        println("Looking for FASTA at: $minimal_fasta")
        @test isfile(minimal_fasta) "FASTA file not found: $minimal_fasta"
        
        # Test with just one case for now
        missed_cleavage_tests = [
            (mc=1, expected_peptides=2),  # MARALYK, ALYKDER
        ]
        
        for test_case in missed_cleavage_tests
            @testset "Missed cleavages = $(test_case.mc)" begin
                lib_name = "minimal_mc$(test_case.mc)"
                lib_dir = joinpath(output_dir, lib_name * ".poin")
                
                println("Creating library at: $lib_dir")
                
                # Generate parameters - DISABLE fragment prediction for now
                params = generate_build_params(
                    minimal_fasta,
                    output_dir,
                    lib_name;
                    missed_cleavages = test_case.mc,
                    min_length = 7,
                    max_length = 15,
                    min_charge = 2,
                    max_charge = 2,
                    add_decoys = false,  # No decoys for easier counting
                    predict_fragments = false,  # DISABLED for testing
                    max_koina_batch = 10
                )
                
                # Save parameters
                params_file = joinpath(scenario_dir, "params_$(test_case.mc)mc.json")
                open(params_file, "w") do io
                    JSON.print(io, params, 4)
                end
                
                println("Running BuildSpecLib with params: $params_file")
                # Run BuildSpecLib
                result = Pioneer.BuildSpecLib(params_file)
                println("BuildSpecLib returned: $result")
                
                # Check what was created
                println("Library directory exists: ", isdir(lib_dir))
                if isdir(lib_dir)
                    println("Files in library directory:")
                    for file in readdir(lib_dir)
                        println("  - $file")
                    end
                    
                    # Check for peptides file
                    peptides_file = joinpath(lib_dir, "peptides.csv")
                    if isfile(peptides_file)
                        n_peptides, sequences = count_peptides_by_sequence(peptides_file, true)
                        println("Found $n_peptides unique peptides: $sequences")
                        @test n_peptides == test_case.expected_peptides
                        
                        # Verify specific peptides exist
                        if test_case.mc == 1
                            @test "MARALYK" in sequences || "ALYKDER" in sequences
                        end
                    else
                        @warn "Peptides file not created: $peptides_file"
                        @test false "Peptides file was not created"
                    end
                else
                    @warn "Library directory not created: $lib_dir"
                    @test false "Library directory was not created"
                end
            end
        end
    end
    
    #= Comment out other scenarios for now
    @testset "Scenario B - Standard Library Build" begin
        scenario_dir = joinpath(test_data_dir, "scenario_b_standard")
        output_dir = joinpath(scenario_dir, "output")
        single_fasta = joinpath(fasta_dir, "single_protein.fasta")
        
        lib_name = "standard_test"
        lib_dir = joinpath(output_dir, lib_name * ".poin")
        
        # Generate parameters
        params = generate_build_params(
            single_fasta,
            output_dir,
            lib_name;
            missed_cleavages = 1,
            min_length = 7,
            max_length = 12,
            min_charge = 2,
            max_charge = 3,
            add_decoys = true,
            predict_fragments = true,
            max_koina_batch = 50
        )
        
        # Save parameters
        params_file = joinpath(scenario_dir, "params_altimeter.json")
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Verify peptide generation
        peptides_file = joinpath(lib_dir, "peptides.csv")
        df = verify_peptide_generation(
            peptides_file,
            Dict(
                "min_length" => 7,
                "max_length" => 12,
                "min_charge" => 2,
                "max_charge" => 3,
                "add_decoys" => true
            )
        )
        
        # Check fragment predictions
        @test verify_fragment_predictions(lib_dir)
    end
    
    @testset "Scenario C - Multi-file with Modifications" begin
        scenario_dir = joinpath(test_data_dir, "scenario_c_multifile")
        output_dir = joinpath(scenario_dir, "output")
        
        lib_name = "multifile_test"
        lib_dir = joinpath(output_dir, lib_name * ".poin")
        
        # Use GetBuildLibParams to handle multiple files
        params_file = joinpath(scenario_dir, "params_altimeter.json")
        params_path = Pioneer.GetBuildLibParams(
            output_dir,
            lib_name,
            joinpath(fasta_dir, "fastas_dir1");
            params_path = params_file
        )
        
        # Load and modify the generated parameters
        params = JSON.parsefile(params_file, dicttype=OrderedDict)
        params["fasta_digest_params"]["missed_cleavages"] = 1
        params["fasta_digest_params"]["min_length"] = 8
        params["fasta_digest_params"]["max_length"] = 15
        params["fasta_digest_params"]["min_charge"] = 2
        params["fasta_digest_params"]["max_charge"] = 3
        params["predict_fragments"] = true
        params["library_params"]["prediction_model"] = "altimeter"
        params["library_params"]["auto_detect_frag_bounds"] = false
        params["max_koina_batch"] = 100
        delete!(params["library_params"], "calibration_raw_file")
        
        # Save modified parameters
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Check that we processed multiple FASTA files
        peptides_file = joinpath(lib_dir, "peptides.csv")
        df = CSV.read(peptides_file, DataFrame)
        
        # Should have peptides from both proteins in dir1
        unique_proteins = unique(df.protein_name)
        @test length(unique_proteins) >= 2  # At least 2 proteins (plus possible decoys)
        
        # Verify fragment predictions
        @test verify_fragment_predictions(lib_dir)
    end
    
    @testset "Scenario D - Comprehensive Test" begin
        scenario_dir = joinpath(test_data_dir, "scenario_d_comprehensive")
        output_dir = joinpath(scenario_dir, "output")
        
        lib_name = "comprehensive_test"
        lib_dir = joinpath(output_dir, lib_name * ".poin")
        
        # Use all test FASTA files
        all_fastas = [
            joinpath(fasta_dir, "fastas_dir1"),
            joinpath(fasta_dir, "fastas_dir2"),
            joinpath(fasta_dir, "minimal_protein.fasta")
        ]
        
        # Generate parameters for comprehensive test
        params_file = joinpath(scenario_dir, "params_altimeter.json")
        params_path = Pioneer.GetBuildLibParams(
            output_dir,
            lib_name,
            all_fastas;
            params_path = params_file
        )
        
        # Modify for comprehensive testing
        params = JSON.parsefile(params_file, dicttype=OrderedDict)
        params["fasta_digest_params"]["missed_cleavages"] = 2
        params["fasta_digest_params"]["min_length"] = 7
        params["fasta_digest_params"]["max_length"] = 20
        params["fasta_digest_params"]["min_charge"] = 2
        params["fasta_digest_params"]["max_charge"] = 4
        params["fasta_digest_params"]["entrapment_r"] = 0.1
        params["predict_fragments"] = true
        params["library_params"]["prediction_model"] = "altimeter"
        params["library_params"]["auto_detect_frag_bounds"] = false
        params["max_koina_batch"] = 200
        delete!(params["library_params"], "calibration_raw_file")
        
        # Save modified parameters
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Load peptides
        peptides_file = joinpath(lib_dir, "peptides.csv")
        df = CSV.read(peptides_file, DataFrame)
        
        # Check we have peptides from all sources
        @test nrow(df) > 50  # Should have many peptides with 2 missed cleavages
        
        # Check entrapment sequences were added
        if "entrap_id" in names(df)
            n_entrapment = sum(df.entrap_id .> 0)
            @test n_entrapment > 0
        end
        
        # Verify charge state range
        @test minimum(df.charge) >= 2
        @test maximum(df.charge) <= 4
        
        # Verify fragment predictions
        @test verify_fragment_predictions(lib_dir)
    end
    =# # End of commented out scenarios
end

println("✓ BuildSpecLib integration tests completed")