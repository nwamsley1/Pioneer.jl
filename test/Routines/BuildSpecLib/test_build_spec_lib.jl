# BuildSpecLib Integration Tests

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
# Note: Pioneer uses cleavage regex [KR][^_|$] by default, which allows cleavage
# after K/R even when followed by P (unlike standard trypsin [KR][^P]).
# This is intentional to provide more flexible digestion patterns.
# The regex prevents cleavage only when K/R is followed by underscore, pipe, or dollar sign.
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
    entrapment_r::Float64 = 0.0,
    prec_mz_min::Float64 = 390.0,
    prec_mz_max::Float64 = 1010.0
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
    params["library_params"]["prec_mz_min"] = prec_mz_min
    params["library_params"]["prec_mz_max"] = prec_mz_max
    params["library_params"]["max_frag_rank"] = 50

    # Set calibration_raw_file to empty for testing (will trigger warning but that's expected)
    params["calibration_raw_file"] = ""

    # Set library_path properly (replaces placeholder from template)
    params["library_path"] = joinpath(output_dir, lib_name)

    # Set paths
    params["fasta_paths"] = [fasta_path]
    params["fasta_names"] = [uppercase(splitext(basename(fasta_path))[1])]
    
    # Set Koina parameters
    params["predict_fragments"] = predict_fragments
    params["max_koina_requests"] = 12
    params["max_koina_batch"] = max_koina_batch
    params["include_contaminants"] = false  # Don't add contaminants for testing
    
    # Set regex patterns for FASTA headers
    params["fasta_header_regex_accessions"] = ["^>(\\S+)"]
    params["fasta_header_regex_genes"] = ["GN=(\\S+)"]
    params["fasta_header_regex_proteins"] = ["\\s+(.+?)\\s+OS="]
    params["fasta_header_regex_organisms"] = ["OS=(.+?)(?:\\s+GN=|\$)"]
    
    # Clear variable modifications for simpler testing
    params["variable_mods"] = Dict(
        "pattern" => [],
        "mass" => [],
        "name" => []
    )
    
    return params
end

# Generate parameters with variable modifications
function generate_build_params_with_var_mods(
    fasta_path::String,
    output_dir::String,
    lib_name::String;
    max_var_mods::Int = 1,
    include_oxidation::Bool = true,
    include_phospho::Bool = false,
    length_to_frag_count_multiple::Float64 = 2.0,
    missed_cleavages::Int = 1,
    min_length::Int = 7,
    max_length::Int = 30,
    min_charge::Int = 2,
    max_charge::Int = 4,
    add_decoys::Bool = true,
    predict_fragments::Bool = true,
    max_koina_batch::Int = 100,
    prec_mz_min::Float64 = 390.0,
    prec_mz_max::Float64 = 1010.0
)
    # Start with base parameters
    params = generate_build_params(
        fasta_path, output_dir, lib_name;
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
        min_charge=min_charge,
        max_charge=max_charge,
        add_decoys=add_decoys,
        predict_fragments=predict_fragments,
        max_koina_batch=max_koina_batch,
        prec_mz_min=prec_mz_min,
        prec_mz_max=prec_mz_max
    )
    
    # Set max variable modifications
    params["fasta_digest_params"]["max_var_mods"] = max_var_mods
    
    # Set fragment count multiplier
    params["library_params"]["length_to_frag_count_multiple"] = length_to_frag_count_multiple
    
    # Configure variable modifications
    var_mods = Dict(
        "pattern" => String[],
        "mass" => Float64[],
        "name" => String[]
    )
    
    if include_oxidation && max_var_mods > 0
        push!(var_mods["pattern"], "M")
        push!(var_mods["mass"], 15.994915)
        push!(var_mods["name"], "Unimod:35")
    end
    
    if include_phospho && max_var_mods > 0
        push!(var_mods["pattern"], "C")
        push!(var_mods["mass"], 57.021464)
        push!(var_mods["name"], "Unimod:4")
    end
    
    params["variable_mods"] = var_mods
    
    return params
end

# Calculate expected number of modification combinations
function calculate_expected_combinations(
    sequence::String,
    mod_patterns::Vector{String},
    max_var_mods::Int
)
    if max_var_mods == 0 || isempty(mod_patterns)
        return 1  # Only unmodified
    end
    
    # Count occurrences of each modifiable residue
    mod_counts = Dict{String, Int}()
    for pattern in mod_patterns
        mod_counts[pattern] = count(pattern, sequence)
    end
    
    total_combinations = 1  # Start with unmodified
    
    # Calculate combinations for each number of modifications
    for n_mods in 1:min(max_var_mods, sum(values(mod_counts)))
        if length(mod_patterns) == 1
            # Single modification type
            pattern = mod_patterns[1]
            n_sites = mod_counts[pattern]
            if n_mods <= n_sites
                total_combinations += binomial(n_sites, n_mods)
            end
        else
            # Multiple modification types - more complex calculation
            # For simplicity, we'll calculate some common cases
            if n_mods == 1
                # Single mod from any type
                for pattern in mod_patterns
                    total_combinations += mod_counts[pattern]
                end
            elseif n_mods == 2 && length(mod_patterns) == 2
                # Two mods - could be same type or different
                p1, p2 = mod_patterns
                n1, n2 = mod_counts[p1], mod_counts[p2]
                
                # Two of same type
                if n1 >= 2
                    total_combinations += binomial(n1, 2)
                end
                if n2 >= 2
                    total_combinations += binomial(n2, 2)
                end
                
                # One of each type
                total_combinations += n1 * n2
            end
            # For higher complexity, would need recursive calculation
        end
    end
    
    return total_combinations
end

# Verify fragment counts in library
function verify_fragment_counts(
    lib_dir::String,
    expected_multiplier::Float64,
    max_frag_rank::Int = 255
)
    # Load precursors to get peptide lengths
    precursors_file = joinpath(lib_dir, "precursors_table.arrow")
    precursors = Arrow.Table(precursors_file)
    
    # Load fragments from JLD2 files
    fragments_file = joinpath(lib_dir, "detailed_fragments.jld2")
    indices_file = joinpath(lib_dir, "precursor_to_fragment_indices.jld2")
    
    if !isfile(fragments_file) || !isfile(indices_file)
        @warn "Fragment files not found in $lib_dir"
        return false
    end
    
    # Load the fragment data
    fragments = load(fragments_file)["data"]
    pid_to_fid = load(indices_file)["pid_to_fid"]
    
    all_valid = true
    
    # Verify fragment counts for each precursor
    for prec_idx in 1:length(precursors.sequence)
        # Get fragment range for this precursor
        # pid_to_fid is 1-indexed array where element i gives start index for precursor i
        start_idx = pid_to_fid[prec_idx]
        end_idx = pid_to_fid[prec_idx + 1] - 1
        frag_count = end_idx - start_idx + 1
        
        # Check expected fragment count
        seq_length = length(precursors.sequence[prec_idx])
        expected_max = min(max_frag_rank, round(Int, seq_length * expected_multiplier) + 1)
        
        if frag_count > expected_max
            @warn "Precursor $prec_idx has $frag_count fragments, expected max $expected_max"
            all_valid = false
        end
    end
    
    return all_valid
end

# Count unique precursors with modifications
function count_precursors_with_mods(peptides_file::String)
    df = CSV.read(peptides_file, DataFrame)
    
    # Filter to targets only
    df = filter(row -> row.target, df)
    
    # Count by modification status
    unmodified = 0
    modified = 0
    
    for row in eachrow(df)
        if haskey(row, :structural_mods) && !ismissing(row.structural_mods) && row.structural_mods != ""
            modified += 1
        else
            unmodified += 1
        end
    end
    
    return (total=nrow(df), unmodified=unmodified, modified=modified)
end

# Verify library structure was created
function verify_library_structure(lib_dir::String)
    @test isdir(lib_dir)
    
    # Check for essential files
    precursors_file = joinpath(lib_dir, "precursors_table.arrow")
    @test isfile(precursors_file)
    
    config_file = joinpath(lib_dir, "config.json")
    @test isfile(config_file)
    
    proteins_file = joinpath(lib_dir, "proteins_table.arrow")
    @test isfile(proteins_file)
    
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
    
    for charge in range(expected_params["min_charge"], expected_params["max_charge"])
        @test any(df.charge .== charge)
    end
    # Count targets and decoys
    n_targets = sum(df.target)
    n_decoys = sum(.!df.target)
    
    if expected_params["add_decoys"]
        @test n_decoys > 0
        # Decoys should be roughly equal to targets (allowing some variation)
        @test n_targets == n_decoys
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
        @test isfile(minimal_fasta)
        
        # Update expected peptides for new sequence MCMKALYKMSRPKMCER
        # With [KR][^_|$] regex, this should cleave after K and R (including RP)
        # Cleaving at positions: K4, K8, R11(before P), K14, R17
        # Note: Need to set prec_mz_min=0 and prec_mz_max=5000 to avoid filtering
        # small peptides by their m/z values
        
        missed_cleavage_tests = [
            (mc=1, expected_peptides=9),  # 5 base + 4 with 1 mc
            #(mc=2, expected_peptides=12), # Above + 3 with 2 mc
            #(mc=3, expected_peptides=14), # Above + 2 with 3 mc
        ]
        
        for test_case in missed_cleavage_tests
            @testset "Missed cleavages = $(test_case.mc)" begin
                lib_name = "minimal_mc$(test_case.mc)"
                lib_dir = joinpath(output_dir, lib_name * ".poin")
                
                println("Creating library at: $lib_dir")
                
                # Generate parameters with wide m/z range to avoid filtering
                params = generate_build_params(
                    minimal_fasta,
                    output_dir,
                    lib_name;
                    missed_cleavages = test_case.mc,
                    min_length = 1,
                    max_length = 15,
                    min_charge = 2,
                    max_charge = 2,
                    add_decoys = false,  # No decoys for easier counting
                    predict_fragments = true,  # Required for BuildSpecLib
                    max_koina_batch = 10,
                    prec_mz_min = 0.0,  # Set to 0 to avoid filtering small peptides
                    prec_mz_max = 5000.0  # Set high to avoid filtering any peptides
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
                    
                    # Check precursors table
                    precursors_file = joinpath(lib_dir, "precursors_table.arrow")
                    if isfile(precursors_file)
                        # Read Arrow table and extract sequences
                        precursors = Arrow.Table(precursors_file)
                        sequences = unique(precursors[:sequence])
                        println("Found $(length(sequences)) unique peptide sequences")
                        println("First 10 sequences: $(sequences[1:min(10, length(sequences))])")
                        
                        # For minimal protein MCMKALYKMSRPKMCER:
                        # Using cleavage regex [KR][^_|$] (allows cleavage after K/R even when followed by P):
                        # Cleavage sites: K4, K8, R11, K14, R17
                        # With 0 missed cleavages (5 peptides): MCMK, ALYK, MSR, PK, MCER
                        # With 1 missed cleavage (4 additional): MCMKALYK, ALYKMSR, MSRPK, PKMCER
                        # Total with mc=1: 9 peptides
                        if test_case.mc == 1
                            # Check we got the expected peptides with 1 missed cleavage
                            @test length(sequences) == test_case.expected_peptides
                            # Check for some expected peptides (0-MC)
                            @test "MCMK" in sequences || "MCmK" in sequences  # May have lowercase m for oxidation
                            @test "ALYK" in sequences
                            @test "MSR" in sequences || "mSR" in sequences
                            # Check for some 1-MC peptides
                            @test "MCMKALYK" in sequences || "MCmKALYK" in sequences || "MCMKALYK" in sequences
                            @test "PKMCER" in sequences || "PKmCER" in sequences || "PKMCER" in sequences
                        end
                    else
                        @warn "Precursors file not created: $precursors_file"
                        @test false
                    end
                else
                    @warn "Library directory not created: $lib_dir"
                    @test false
                end
            end
        end
    end
    
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
        
        # Check precursors table
        precursors_file = joinpath(lib_dir, "precursors_table.arrow")
        if isfile(precursors_file)
            precursors = Arrow.Table(precursors_file)
            sequences = unique(precursors[:sequence])
            println("Scenario B: Found $(length(sequences)) unique peptide sequences")
            # With decoys enabled, we should have both targets and decoys
            @test length(sequences) > 0
        else
            @warn "Precursors file not created: $precursors_file"
            @test false
        end
    end
    
    @testset "Scenario C - Multi-file with Modifications" begin
        scenario_dir = joinpath(test_data_dir, "scenario_c_multifile")
        output_dir = joinpath(scenario_dir, "output")
        
        lib_name = "multifile_test"
        lib_dir = joinpath(output_dir, lib_name * ".poin")
        
        # Use GetBuildLibParams to handle multiple files
        # Use full template for complex scenarios
        params_file = joinpath(scenario_dir, "params_altimeter.json")
        params_path = Pioneer.GetBuildLibParams(
            output_dir,
            lib_name,
            joinpath(fasta_dir, "fastas_dir1");
            params_path = params_file,
            simplified = false  # Use full template for complex test
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
        params["include_contaminants"] = false  # Don't add contaminants for testing
        delete!(params["library_params"], "calibration_raw_file")
        # Clear variable modifications for simpler testing
        params["variable_mods"] = Dict("pattern" => [], "mass" => [], "name" => [])
        
        # Save modified parameters
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Check precursors table for multiple FASTA files
        precursors_file = joinpath(lib_dir, "precursors_table.arrow")
        if isfile(precursors_file)
            precursors = Arrow.Table(precursors_file)
            sequences = unique(precursors[:sequence])
            proteome_ids = unique(precursors[:proteome_identifiers])
            println("Scenario C: Found $(length(sequences)) unique sequences from $(length(proteome_ids)) proteomes")
            # Should have peptides from both files in fastas_dir1
            @test length(proteome_ids) >= 2  # At least 2 proteomes
        else
            @warn "Precursors file not created: $precursors_file"
            @test false
        end
    end
    
    @testset "Variable Modification Tests" begin
        scenario_base_dir = joinpath(test_data_dir, "scenario_var_mods")
        minimal_fasta = joinpath(fasta_dir, "minimal_protein.fasta")
        
        # Test different max_var_mods values
        var_mod_tests = [
            (max_var_mods=0, include_oxidation=true, include_phospho=true, desc="No mods allowed"),
            (max_var_mods=1, include_oxidation=true, include_phospho=false, desc="Max 1 mod - M oxidation only"),
            (max_var_mods=2, include_oxidation=true, include_phospho=true, desc="Max 2 mods - M ox and C carbamidomethyl"),
            (max_var_mods=3, include_oxidation=true, include_phospho=true, desc="Max 3 mods"),
            (max_var_mods=4, include_oxidation=true, include_phospho=true, desc="Max 4 mods"),
        ]
        
        for test_case in var_mod_tests
            @testset "$(test_case.desc)" begin
                output_dir = joinpath(scenario_base_dir, "max_$(test_case.max_var_mods)_mods", "output")
                lib_name = "var_mods_$(test_case.max_var_mods)"
                lib_dir = joinpath(output_dir, lib_name * ".poin")
                
                # Generate parameters with variable modifications
                params = generate_build_params_with_var_mods(
                    minimal_fasta,
                    output_dir,
                    lib_name;
                    max_var_mods = test_case.max_var_mods,
                    include_oxidation = test_case.include_oxidation,
                    include_phospho = test_case.include_phospho,
                    missed_cleavages = 1,
                    min_length = 1,  # Lower to catch all peptides
                    max_length = 20,
                    add_decoys = false,  # Simplify for testing
                    predict_fragments = true,  # Required for library to build successfully
                    max_koina_batch = 50,
                    prec_mz_min = 0.0,  # Wide range to avoid filtering
                    prec_mz_max = 5000.0  # Wide range to avoid filtering
                )
                
                # Save parameters
                params_file = joinpath(dirname(output_dir), "params_var_$(test_case.max_var_mods).json")
                mkpath(dirname(params_file))
                open(params_file, "w") do io
                    JSON.print(io, params, 4)
                end
                
                # Run BuildSpecLib
                @test Pioneer.BuildSpecLib(params_file) === nothing
                
                # Verify precursor generation with modifications
                precursors_file = joinpath(lib_dir, "precursors_table.arrow")
                @info "precursors_file: $precursors_file"
                if isfile(precursors_file)
                    precursors = Arrow.Table(precursors_file)
                    
                    # Count unique sequences and modifications
                    unique_seqs = unique(precursors[:sequence])
                    println("Max var mods = $(test_case.max_var_mods): Found $(length(unique_seqs)) unique sequences")
                    
                    # Verify modification combinations
                    # For sequence MCMKALYKMSRPKMCER with 4 M's and 3 C's
                    if test_case.max_var_mods == 0
                        # Should have no modifications
                        @test !any(contains.(String.(precursors[:structural_mods]), "["))
                    elseif test_case.max_var_mods == 1 && test_case.include_oxidation && !test_case.include_phospho
                        # Should have unmodified + 4 positions for M oxidation
                        # Each peptide can have 0 or 1 oxidation
                        # Count precursors with oxidation marker (Unimod:35)
                        has_ox = sum(contains.(String.(precursors[:structural_mods]), "Unimod:35"))
                        println("  Found $has_ox precursors with oxidation")
                        @test has_ox > 0
                    elseif test_case.max_var_mods >= 2 && test_case.include_oxidation && test_case.include_phospho
                        # Should have combinations of M oxidation and C carbamidomethylation
                        has_ox = sum(contains.(String.(precursors[:structural_mods]), "Unimod:35"))
                        has_carbamid = sum(contains.(String.(precursors[:structural_mods]), "Unimod:4"))
                        println("  Found $has_ox with oxidation, $has_carbamid with carbamidomethylation")
                        @test has_ox > 0 || has_carbamid > 0
                    end
                else
                    @warn "Precursors file not created: $precursors_file"
                    @test false
                end
            end
        end
    end
    
    
    @testset "Fragment Count Rule Tests" begin
        scenario_base_dir = joinpath(test_data_dir, "scenario_frag_count")
        minimal_fasta = joinpath(fasta_dir, "minimal_protein.fasta")
        
        # Test different fragment count multipliers
        frag_multiplier_tests = [
            (multiplier=1.0, desc="1x peptide length"),
            (multiplier=2.0, desc="2x peptide length (default)"),
            (multiplier=3.0, desc="3x peptide length"),
            (multiplier=4.0, desc="4x peptide length"),
        ]
        
        for test_case in frag_multiplier_tests
            @testset "$(test_case.desc)" begin
                output_dir = joinpath(scenario_base_dir, "mult_$(Int(test_case.multiplier))", "output")
                lib_name = "frag_mult_$(Int(test_case.multiplier))"
                lib_dir = joinpath(output_dir, lib_name * ".poin")
                
                # Generate parameters with specific fragment multiplier
                params = generate_build_params_with_var_mods(
                    minimal_fasta,
                    output_dir,
                    lib_name;
                    length_to_frag_count_multiple = test_case.multiplier,
                    max_var_mods = 0,  # No mods to simplify fragment testing
                    missed_cleavages = 1,
                    min_length = 7,
                    max_length = 20,
                    add_decoys = false,
                    predict_fragments = true,
                    max_koina_batch = 50,
                    prec_mz_min = 0.0,  # Wide range to avoid filtering
                    prec_mz_max = 5000.0  # Wide range to avoid filtering
                )
                
                # Set max_frag_rank to a reasonable value
                params["library_params"]["max_frag_rank"] = 100
                
                # Save parameters
                params_file = joinpath(dirname(output_dir), "params_frag_$(Int(test_case.multiplier)).json")
                mkpath(dirname(params_file))
                open(params_file, "w") do io
                    JSON.print(io, params, 4)
                end
                
                # Run BuildSpecLib
                @test Pioneer.BuildSpecLib(params_file) === nothing
                
                # Verify fragment counts
                if verify_library_structure(lib_dir)
                    is_valid = verify_fragment_counts(lib_dir, test_case.multiplier, 100)
                    @test is_valid
                    
                    # Additional verification - check specific precursor
                    precursors_file = joinpath(lib_dir, "precursors_table.arrow")
                    if isfile(precursors_file)
                        precursors = Arrow.Table(precursors_file)
                        if length(precursors) > 0
                            first_seq = precursors[:sequence][1]
                            seq_length = length(first_seq)
                            expected_frags = min(100, round(Int, seq_length * test_case.multiplier) + 1)
                            println("  First precursor '$(first_seq)' (length $seq_length)")
                            println("  Expected max fragments: $expected_frags with multiplier $(test_case.multiplier)")
                        end
                    end
                else
                    @warn "Library structure verification failed for multiplier $(test_case.multiplier)"
                    @test false
                end
            end
        end
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
            params_path = params_file,
            simplified = false  # Use full template for comprehensive test
        )
        
        # Modify for comprehensive testing
        params = JSON.parsefile(params_file, dicttype=OrderedDict)
        params["fasta_digest_params"]["missed_cleavages"] = 2
        params["fasta_digest_params"]["min_length"] = 7
        params["fasta_digest_params"]["max_length"] = 20
        params["fasta_digest_params"]["min_charge"] = 2
        params["fasta_digest_params"]["max_charge"] = 4
        params["fasta_digest_params"]["entrapment_r"] = 1  # Integer ratio for entrapment
        params["predict_fragments"] = true
        params["library_params"]["prediction_model"] = "altimeter"
        params["library_params"]["auto_detect_frag_bounds"] = false
        params["max_koina_batch"] = 200
        params["include_contaminants"] = false  # Don't add contaminants for testing
        delete!(params["library_params"], "calibration_raw_file")
        # Clear variable modifications for simpler testing
        params["variable_mods"] = Dict("pattern" => [], "mass" => [], "name" => [])
        
        # Save modified parameters
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Check comprehensive precursors table
        precursors_file = joinpath(lib_dir, "precursors_table.arrow")
        if isfile(precursors_file)
            precursors = Arrow.Table(precursors_file)
            sequences = unique(precursors[:sequence])
            proteome_ids = unique(precursors[:proteome_identifiers])
            entrapment_groups = unique(precursors[:entrapment_group_id])
            println("Scenario D: Found $(length(sequences)) sequences from $(length(proteome_ids)) proteomes")
            println("Entrapment groups: $entrapment_groups")
            
            # Should have peptides from all sources (3 FASTA sources)
            @test length(proteome_ids) >= 3
            
            # Check entrapment sequences were added (group_id > 0)
            n_entrapment = sum(precursors[:entrapment_group_id] .> 0)
            @test n_entrapment > 0  # Should have entrapment sequences
            
            # Verify charge state range
            charges = unique(precursors[:prec_charge])
            @test minimum(charges) >= 2
            @test maximum(charges) <= 4
        else
            @warn "Precursors file not created: $precursors_file"
            @test false
        end
    end
    
    @testset "Prosit Model Integration Test" begin
        scenario_dir = joinpath(test_data_dir, "scenario_prosit")
        output_dir = joinpath(scenario_dir, "output")
        minimal_fasta = joinpath(fasta_dir, "minimal_protein.fasta")
        
        # Verify FASTA file exists
        @test isfile(minimal_fasta)
        
        lib_name = "prosit_test"
        lib_dir = joinpath(output_dir, lib_name * ".poin")
        
        # Generate parameters using Prosit model instead of altimeter
        params = generate_build_params(
            minimal_fasta,
            output_dir,
            lib_name;
            missed_cleavages = 1,
            min_length = 7,
            max_length = 15,
            min_charge = 2,
            max_charge = 3,
            add_decoys = true,
            predict_fragments = true,
            max_koina_batch = 50,
            prec_mz_min = 390.0,
            prec_mz_max = 1010.0
        )
        
        # Override to use Prosit instead of altimeter
        params["library_params"]["prediction_model"] = "prosit_2020_hcd"
        
        # Save parameters
        params_file = joinpath(scenario_dir, "params_prosit.json")
        mkpath(dirname(params_file))
        open(params_file, "w") do io
            JSON.print(io, params, 4)
        end
        
        # Run BuildSpecLib with Prosit
        @test Pioneer.BuildSpecLib(params_file) === nothing
        
        # Verify outputs
        @test verify_library_structure(lib_dir)
        
        # Check precursors table
        precursors_file = joinpath(lib_dir, "precursors_table.arrow")
        if isfile(precursors_file)
            precursors = Arrow.Table(precursors_file)
            sequences = unique(precursors[:sequence])
            println("Prosit Integration: Found $(length(sequences)) unique peptide sequences")
            
            # Should have peptides from minimal protein
            @test length(sequences) > 0
            
            # Verify charge states
            charges = unique(precursors[:prec_charge])
            @test minimum(charges) >= 2
            @test maximum(charges) <= 3
            
            # Should have both targets and decoys
            n_targets = sum(.!precursors[:is_decoy])  # Targets are non-decoys
            n_decoys = sum(precursors[:is_decoy])     # Decoys are marked as is_decoy=true
            @test n_targets > 0
            @test n_decoys > 0
        else
            @warn "Precursors file not created: $precursors_file"
            @test false
        end
        
        println("✓ Prosit integration test completed successfully")
    end
    
end

println("✓ BuildSpecLib integration tests completed")