
"""
    getSearchParams(template_path::String, lib_path::String, ms_data_path::String, results_path::String)

Creates a new search parameter file based on a template, with updated file paths.

The function reads a template JSON configuration file and creates a new 'search_parameters.json' 
in the current working directory with updated paths while preserving all other settings.

Arguments:
- template_path: Path to the template JSON configuration file
- lib_path: Path to the library file (.poin)
- ms_data_path: Path to the MS data directory
- results_path: Path where results will be stored

Returns:
- String: Path to the newly created search parameters file

Example:
```julia
output_file = getBuildLibParams(
    "/path/to/template.json",
    "/path/to/output/dir",
    "/path/to/library/name",
    "/path/to/fasta/directory"
)
```

Note: The function preserves all configuration settings from the template 
(fragment settings, search settings, etc.) while only modifying the paths section.
"""
function getSearchParams(lib_path::String, ms_data_path::String, results_path::String)
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
    output_path = joinpath(pwd(), "search_parameters.json")
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end

"""
    getBuildLibParams(template_path::String, out_dir::String, lib_name::String, fasta_dir::String)

Creates a new library build parameter file based on a template, with updated paths and FASTA files.

Arguments:
- template_path: Path to the template JSON file
- out_dir: New output directory path
- lib_name: New library name path
- fasta_dir: Directory to search for FASTA files

Returns:
- Path to the newly created parameters file

Example:
```julia
output_file = getSearchParams(
    "/path/to/template.json",
    "/path/to/library.poin",
    "/path/to/ms_data",
    "/path/to/results"
)
```
"""
function getBuildLibParams(template_path::String, out_dir::String, lib_name::String, fasta_dir::String)
    # Read the JSON template
    config = JSON.parsefile(template_path)
    
    # Find all FASTA files in the specified directory
    fasta_files = String[]
    fasta_names = String[]
    
    for file in readdir(fasta_dir, join=true)
        if endswith(lowercase(file), ".fasta") || endswith(lowercase(file), ".fasta.gz")
            push!(fasta_files, file)
            # Get base name without extension and directory
            base_name = uppercase(splitext(basename(splitext(file)[1]))[1])
            push!(fasta_names, base_name)
        end
    end
    
    # Update paths and FASTA information in the configuration
    config["out_dir"] = out_dir
    config["lib_name"] = lib_name
    config["new_lib_name"] = lib_name  # Update both lib_name and new_lib_name
    config["fasta_paths"] = fasta_files
    config["fasta_names"] = fasta_names
    
    # Set output name based on the library name
    config["out_name"] = basename(lib_name) * ".tsv"
    
    # Write the modified configuration to build_parameters.json in current directory
    output_path = joinpath(pwd(), "build_parameters.json")
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end

