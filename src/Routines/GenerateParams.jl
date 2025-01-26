
"""
    getSearchParams(template_path::String, lib_path::String, ms_data_path::String, results_path::String)

Creates a new search parameter file based on a template, with updated file paths.

The function reads a template JSON configuration file and creates a new 'search_parameters.json' 
in the current working directory with updated paths while preserving all other settings.

Arguments:
- lib_path: Path to the library file (.poin)
- ms_data_path: Path to the MS data directory
- results_path: Path where results will be stored
- params_path: Path to folder or .json file in which to write the template parameters file. Defaults to joinpath(pwd(), "./search_parameters.json")

Returns:
- String: Path to the newly created search parameters file

Example:
```julia
output_path = getSearchParams(
    "/path/to/speclib.poin",
    "/path/to/ms/data/dir",
    "/path/to/output/dir"
)
```
"""
function GetSearchParams(lib_path::String, ms_data_path::String, results_path::String; params_path::Union{String, Missing} = missing)
    if ismissing(params_path)
        output_path = joinpath(pwd(), "search_parameters.json")
    elseif isdir(params_path)
        output_path = joinpath(params_path, "search_parameters.json")
    elseif isfile(params_path)
        output_path = params_path
    else
        throw("User supplied `params_path` is not a path to a file or a directory")
    end
    # Read the JSON template
    config = JSON.parsefile(joinpath(@__DIR__, "../../data/example_config/defaultSearchParams.json"))
    
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
    @info "Writing default parameters .json to: $output_path"
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end

"""
    GetBuildLibParams(out_dir::String, lib_name::String, fasta_dir::String)

Creates a new library build parameter file with updated paths and automatically discovered FASTA files.
Uses a default template from data/example_config/defaultBuildLibParams.json.

Arguments:
- out_dir: Output directory path
- lib_name: Library name path
- fasta_dir: Directory to search for FASTA files
- params_path: Path to folder or .json file in which to write the template parameters file. Defaults to joinpath(pwd(), "./buildspeclib_params.json")

Returns:
- String: Path to the newly created parameters file
"""
function GetBuildLibParams(out_dir::String, lib_name::String, fasta_dir::String; params_path::Union{String, Missing} = missing)
    if ismissing(params_path)
        output_path = joinpath(pwd(), "buildspeclib_params.json")
    elseif isdir(params_path)
        output_path = joinpath(params_path, "buildspeclib_params.json")
    elseif isfile(params_path)
        output_path = params_path
    else
        throw("User supplied `params_path` is not a path to a file or a directory")
    end
    # Read the template file
    template_path = joinpath(@__DIR__, "../../data/example_config/defaultBuildLibParams.json")
    template_text = read(template_path, String)
    
    # Parse JSON
    config = JSON.parse(template_text)
    
    # Find all FASTA files in the specified directory
    fasta_files = String[]
    fasta_names = String[]
    
    for file in readdir(fasta_dir, join=true)
        if endswith(lowercase(file), ".fasta") || endswith(lowercase(file), ".fasta.gz")
            push!(fasta_files, file)
            base_name = uppercase(splitext(basename(splitext(file)[1]))[1])
            push!(fasta_names, base_name)
        end
    end
    
    # Update values while maintaining structure
    config["out_dir"] = out_dir
    config["lib_name"] = lib_name
    config["new_lib_name"] = lib_name
    config["fasta_paths"] = fasta_files
    config["fasta_names"] = fasta_names
    config["out_name"] = basename(lib_name) * ".tsv"
    
    # Write output using the same formatting as template
    output_path = joinpath(pwd(), "build_parameters.json")
    open(output_path, "w") do io
        # Extract indentation from template
        indent_match = match(r"\n(\s+)\"", template_text)
        indent = indent_match === nothing ? "    " : indent_match[1]
        
        # Write with matching format
        JSON.print(io, config, length(indent))
    end
    
    return output_path
end