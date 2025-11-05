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

# Entry point for PackageCompiler
function main_GetSearchParams(argv=ARGS)::Cint
    
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "library_path"
            help = "Path to spectral library (.poin)"
            arg_type = String
        "ms_data_path"
            help = "Directory containing MS data"
            arg_type = String
        "results_path"
            help = "Output directory for search results"
            arg_type = String
        "--params-path"
            help = "Output path for generated parameters file"
            arg_type = String
            default = joinpath(pwd(), "search_parameters.json")
        "--full"
            help = "Generate full parameter template with all advanced options"
            action = :store_true
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    
    # Determine template type (simplified is default)
    simplified = !parsed_args[:full]
    
    params_path = parsed_args[:params_path]
    try
       GetSearchParams(parsed_args[:library_path],
                       parsed_args[:ms_data_path],
                       parsed_args[:results_path];
                       params_path=params_path,
                       simplified=simplified)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end


# Entry point for PackageCompiler
function main_GetBuildLibParams(argv=ARGS)::Cint
    
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "out_dir"
            help = "Output directory for library"
            arg_type = String
        "fasta_path"
            help = "Directory containing FASTA files or path to a specific FASTA file"
            arg_type = String
        "--params-path"
            help = "Output path for generated parameters file"
            arg_type = String
            default = joinpath(pwd(), "buildspeclib_params.json")
        "--full"
            help = "Generate full parameter template with all advanced options"
            action = :store_true
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)

    # Determine template type (simplified is default)
    simplified = !parsed_args[:full]

    params_path = parsed_args[:params_path]
    try
        # Extract library name and parent directory from output path
        out_dir = dirname(parsed_args[:out_dir])
        lib_name = basename(parsed_args[:out_dir])

        GetBuildLibParams(out_dir,
                          lib_name,
                          parsed_args[:fasta_path];
                          params_path=params_path,
                          simplified=simplified)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

# Entry point for PackageCompiler
function main_GetParseSpecLibParams(argv=ARGS)::Cint
    
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "input_lib_path"
            help = "Input empirical library TSV"
            arg_type = String
        "output_lib_path"
            help = "Output path for processed library"
            arg_type = String
        "--params-path"
            help = "Output path for generated parameters file"
            arg_type = String
            default = joinpath(pwd(), "parsespeclib_params.json")
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    
    params_path = parsed_args[:params_path]
    try
        GetParseSpecLibParams(parsed_args[:input_lib_path], 
                              parsed_args[:output_lib_path];
                              params_path = params_path)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end


"""
    GetSearchParams(lib_path::String, ms_data_path::String, results_path::String; 
                   params_path::Union{String, Missing} = missing,
                   simplified::Bool = true)

Creates a search parameter configuration file with user-specified paths.

The function loads default parameters from either the simplified or full JSON template
(from assets/example_config/) and creates a customized parameter file with the user's
file paths. All other parameters retain their default values and can be modified later.

Arguments:
- lib_path: Path to the spectral library file (.poin)
- ms_data_path: Path to the MS data directory  
- results_path: Path where search results will be stored
- params_path: Output path for the parameter file. Can be a directory (creates search_parameters.json) 
  or full file path. Defaults to "search_parameters.json" in current directory.
- simplified: If true (default), uses simplified template with essential parameters only. 
  If false, uses full template with all advanced options.

Returns:
- String: Path to the newly created search parameters file

Templates used:
- Simplified: `defaultSearchParamsSimplified.json` (basic parameters)
- Full: `defaultSearchParams.json` (all advanced parameters)

Example:
```julia
# Create simplified parameter file
output_path = GetSearchParams(
    "/path/to/speclib.poin",
    "/path/to/ms/data/dir", 
    "/path/to/results/dir"
)

# Create full parameter file with custom output location
output_path = GetSearchParams(
    "/path/to/speclib.poin",
    "/path/to/ms/data/dir",
    "/path/to/results/dir";
    params_path = "/custom/path/my_params.json",
    simplified = false
)
```
"""
function GetSearchParams(lib_path::String, ms_data_path::String, results_path::String; 
                        params_path::Union{String, Missing} = missing,
                        simplified::Bool = true)
    # Clean up any old file handlers in case the program crashed
    GC.gc()
    
    if ismissing(params_path)
        output_path = joinpath(pwd(), "search_parameters.json")
    else
        params_path = expanduser(params_path)
        name, ext = splitext(params_path)
        if isempty(ext)
            mkpath(params_path)
            output_path = joinpath(params_path, "search_parameters.json")
        else
            output_path = params_path
        end
    end
    
    # Choose template based on simplified flag
    template_name = simplified ? "defaultSearchParamsSimplified.json" : "defaultSearchParams.json"
    
    # Read the JSON template and convert to OrderedDict
    config_text = read(asset_path("example_config", template_name), String)
    config = JSON.parse(config_text, dicttype=OrderedDict)

        
    # Update paths in the configuration
    if !isdir(results_path)
        mkdir(results_path)
    end

    config["paths"] = Dict(
        "library" => lib_path,
        "ms_data" => ms_data_path,
        "results" => results_path
    )
    

    # Write the modified configuration to search_parameters.json in current directory
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end

"""
    GetBuildLibParams(out_dir::String, lib_name::String, fasta_inputs; 
                     params_path::Union{String, Missing} = missing,
                     regex_codes::Union{Missing, Dict, Vector} = missing,
                     simplified::Bool = true)

Creates a library building parameter configuration file with user-specified paths and FASTA files.

The function loads default parameters from either the simplified or full JSON template 
(from assets/example_config/) and creates a customized parameter file with the user's
paths and automatically discovered FASTA files. All other parameters retain their
default values and can be modified later.

Arguments:
- out_dir: Output directory path where the library will be built
- lib_name: Name for the spectral library (used for directory and file naming)
- fasta_inputs: FASTA file specification. Can be:
  - A single directory path (String) - searches for .fasta/.fasta.gz files
  - A single FASTA file path (String) 
  - An array of directories and/or FASTA file paths
- params_path: Output path for the parameter file. Can be a directory (creates buildspeclib_params.json)
  or full file path. Defaults to "buildspeclib_params.json" in current directory.
- regex_codes: Optional FASTA header regex patterns for protein annotation extraction. Can be:
  - A single Dict with keys: "accessions", "genes", "proteins", "organisms" (applied to all FASTA files)
  - A Vector of Dicts for positional mapping to fasta_inputs
  - If missing, uses default patterns from the template
- simplified: If true (default), uses simplified template with essential parameters only.
  If false, uses full template with all advanced library building options.

Returns:
- String: Path to the newly created library building parameters file

Templates used:
- Simplified: `defaultBuildLibParamsSimplified.json` (basic parameters)
- Full: `defaultBuildLibParams.json` (all advanced parameters)

The function automatically:
- Discovers FASTA files in specified directories
- Generates appropriate library names from FASTA filenames
- Expands regex patterns to match the number of FASTA files found
- Validates that all specified paths exist and are accessible

Example:
```julia
# Create simplified parameter file with directory of FASTA files
output_path = GetBuildLibParams(
    "/path/to/output", 
    "my_library",
    "/path/to/fasta/directory"
)

# Create full parameter file with specific FASTA files and custom regex
output_path = GetBuildLibParams(
    "/path/to/output",
    "my_library", 
    ["/path/to/human.fasta", "/path/to/yeast.fasta"];
    params_path = "/custom/path/build_params.json",
    regex_codes = Dict("accessions" => "^sp\\|(\\w+)\\|", "genes" => " GN=(\\S+)"),
    simplified = false
)
```
"""
function GetBuildLibParams(out_dir::String, lib_name::String, fasta_inputs; 
                         params_path::Union{String, Missing} = missing,
                         regex_codes::Union{Missing, Dict, Vector} = missing,
                         simplified::Bool = true)
    # Clean up any old file handlers in case the program crashed
    GC.gc()

    if ismissing(params_path)
        output_path = joinpath(pwd(), "buildspeclib_params.json")
    else
        params_path = expanduser(params_path)
        name, ext = splitext(params_path)
        if isempty(ext)
            mkpath(params_path)
            output_path = joinpath(params_path, "buildspeclib_params.json")
        else
            output_path = params_path
        end
    end

    # Choose template based on simplified flag
    template_name = simplified ? "defaultBuildLibParamsSimplified.json" : "defaultBuildLibParams.json"
    
    # Parse JSON
    config_text = read(asset_path("example_config", template_name), String)
    config = JSON.parse(config_text, dicttype=OrderedDict)

    # Process fasta_inputs to get list of FASTA files
    fasta_files = String[]
    fasta_names = String[]
    input_to_files_map = OrderedDict{String, Vector{String}}()  # Maps each input to its FASTA files
    
    # Normalize inputs to always be a vector
    normalized_inputs = if isa(fasta_inputs, String)
        [fasta_inputs]
    elseif isa(fasta_inputs, Vector)
        fasta_inputs
    else
        error("fasta_inputs must be a String or Vector of Strings")
    end
    
    # Process each input (directory or file)
    for input in normalized_inputs
        input_files = String[]
        
        if isdir(input)
            # It's a directory - scan for FASTA files
            for file in readdir(input, join=true)
                if endswith(lowercase(file), ".fasta") || endswith(lowercase(file), ".fasta.gz")
                    push!(input_files, file)
                    push!(fasta_files, file)
                    base_name = uppercase(splitext(basename(splitext(file)[1]))[1])
                    push!(fasta_names, base_name)
                end
            end
            if isempty(input_files)
                @user_warn "No FASTA files found in directory: $input"
            end
        elseif isfile(input)
            # It's a file - validate it's a FASTA file
            if endswith(lowercase(input), ".fasta") || endswith(lowercase(input), ".fasta.gz")
                push!(input_files, input)
                push!(fasta_files, input)
                base_name = uppercase(splitext(basename(splitext(input)[1]))[1])
                push!(fasta_names, base_name)
            else
                error("File is not a FASTA file (.fasta or .fasta.gz): $input")
            end
        else
            error("Path does not exist or is not accessible: $input")
        end
        
        input_to_files_map[input] = input_files
    end
    
    # Handle regex codes expansion/mapping
    if !ismissing(regex_codes)
        # Get default regex fields from template
        default_accessions = config["fasta_header_regex_accessions"][1]
        default_genes = config["fasta_header_regex_genes"][1]
        default_proteins = config["fasta_header_regex_proteins"][1]
        default_organisms = config["fasta_header_regex_organisms"][1]
        
        # Process regex_codes based on type
        if isa(regex_codes, Dict)
            # Single regex set - apply to all FASTA files
            accessions = get(regex_codes, "accessions", default_accessions)
            genes = get(regex_codes, "genes", default_genes)
            proteins = get(regex_codes, "proteins", default_proteins)
            organisms = get(regex_codes, "organisms", default_organisms)
            
            config["fasta_header_regex_accessions"] = fill(accessions, length(fasta_files))
            config["fasta_header_regex_genes"] = fill(genes, length(fasta_files))
            config["fasta_header_regex_proteins"] = fill(proteins, length(fasta_files))
            config["fasta_header_regex_organisms"] = fill(organisms, length(fasta_files))
            
        elseif isa(regex_codes, Vector)
            # Multiple regex sets - apply based on mapping logic
            num_inputs = length(normalized_inputs)
            num_regex_sets = length(regex_codes)
            
            if num_regex_sets == 1
                # Single regex set for all inputs
                accessions = get(regex_codes[1], "accessions", default_accessions)
                genes = get(regex_codes[1], "genes", default_genes)
                proteins = get(regex_codes[1], "proteins", default_proteins)
                organisms = get(regex_codes[1], "organisms", default_organisms)
                
                config["fasta_header_regex_accessions"] = fill(accessions, length(fasta_files))
                config["fasta_header_regex_genes"] = fill(genes, length(fasta_files))
                config["fasta_header_regex_proteins"] = fill(proteins, length(fasta_files))
                config["fasta_header_regex_organisms"] = fill(organisms, length(fasta_files))
                
            elseif num_regex_sets == num_inputs
                # Positional mapping: each regex set maps to corresponding input
                expanded_accessions = String[]
                expanded_genes = String[]
                expanded_proteins = String[]
                expanded_organisms = String[]
                
                for (i, input) in enumerate(normalized_inputs)
                    regex_set = regex_codes[i]
                    accessions = get(regex_set, "accessions", default_accessions)
                    genes = get(regex_set, "genes", default_genes)
                    proteins = get(regex_set, "proteins", default_proteins)
                    organisms = get(regex_set, "organisms", default_organisms)
                    
                    # Apply this regex set to all files from this input
                    num_files = length(input_to_files_map[input])
                    append!(expanded_accessions, fill(accessions, num_files))
                    append!(expanded_genes, fill(genes, num_files))
                    append!(expanded_proteins, fill(proteins, num_files))
                    append!(expanded_organisms, fill(organisms, num_files))
                end
                
                config["fasta_header_regex_accessions"] = expanded_accessions
                config["fasta_header_regex_genes"] = expanded_genes
                config["fasta_header_regex_proteins"] = expanded_proteins
                config["fasta_header_regex_organisms"] = expanded_organisms
                
            else
                @user_warn "Number of regex code sets ($num_regex_sets) must be either 1 or match the number of inputs ($num_inputs)"
            end
        else
            @user_warn "regex_codes must be either a Dict or Vector of Dicts"
        end
    else
        # Use default regex from template, expanded to match number of files
        if haskey(config, "fasta_header_regex_accessions") && length(config["fasta_header_regex_accessions"]) == 1
            config["fasta_header_regex_accessions"] = fill(config["fasta_header_regex_accessions"][1], length(fasta_files))
        end
        if haskey(config, "fasta_header_regex_genes") && length(config["fasta_header_regex_genes"]) == 1
            config["fasta_header_regex_genes"] = fill(config["fasta_header_regex_genes"][1], length(fasta_files))
        end
        if haskey(config, "fasta_header_regex_proteins") && length(config["fasta_header_regex_proteins"]) == 1
            config["fasta_header_regex_proteins"] = fill(config["fasta_header_regex_proteins"][1], length(fasta_files))
        end
        if haskey(config, "fasta_header_regex_organisms") && length(config["fasta_header_regex_organisms"]) == 1
            config["fasta_header_regex_organisms"] = fill(config["fasta_header_regex_organisms"][1], length(fasta_files))
        end
    end
    
    # Update values while maintaining structure
    config["fasta_paths"] = fasta_files
    config["fasta_names"] = fasta_names

    # Set library_path from out_dir and lib_name (replaces placeholder)
    config["library_path"] = joinpath(out_dir, lib_name)

    # Set calibration_raw_file to empty string (optional parameter)
    config["calibration_raw_file"] = ""

    # Write output using the same formatting as template
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end


"""
    GetParseSpecLibParams(input_lib_path::String, output_lib_path::String)

Create a new parameter file for `ParseSpecLib` with updated input and output paths.
Uses a default template from `assets/example_config/defaultParseSpecLibParams.json`.

Arguments:
- input_lib_path: Path to the empirical library TSV file
- output_lib_path: Path where the processed library will be written
- params_path: Optional path to folder or JSON file where the parameters will be
  saved. Defaults to `parsespeclib_params.json` in the working directory.

Returns:
- String: Path to the newly created parameters file
"""
function GetParseSpecLibParams(input_lib_path::String, output_lib_path::String; params_path::Union{String, Missing}=missing)
    GC.gc()

    if ismissing(params_path)
        output_path = joinpath(pwd(), "parsespeclib_params.json")
    else
        params_path = expanduser(params_path)
        name, ext = splitext(params_path)
        if isempty(ext)
            mkpath(params_path)
            output_path = joinpath(params_path, "parsespeclib_params.json")
        else
            output_path = params_path
        end
    end

    config_text = read(asset_path("example_config", "defaultParseEmpiricalLibParams.json"), String)
    config = JSON.parse(config_text, dicttype=OrderedDict)

    config["library_params"]["input_lib_path"] = input_lib_path
    config["library_params"]["output_lib_path"] = output_lib_path

    open(output_path, "w") do io
        JSON.print(io, config, 4)
    end

    return output_path
end