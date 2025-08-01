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
    settings = ArgParseSettings()
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
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    params_path = parsed_args[:"params-path"]
    try
       GetSearchParams(parsed_args[:library_path],
                       parsed_args[:ms_data_path],
                       parsed_args[:results_path];
                       params_path=params_path)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end


# Entry point for PackageCompiler
function main_GetBuildLibParams(argv=ARGS)::Cint
    settings = ArgParseSettings()
    @add_arg_table! settings begin
        "out_dir"
            help = "Output directory for library"
            arg_type = String
        "lib_name"
            help = "Name of the library"
            arg_type = String
        "fasta_dir"
            help = "Directory containing FASTA files"
            arg_type = String
        "--params-path"
            help = "Output path for generated parameters file"
            arg_type = String
            default = joinpath(pwd(), "buildspeclib_params.json")
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)
    params_path = parsed_args[:"params-path"]
    try
        GetBuildLibParams(parsed_args[:out_dir],
                          parsed_args[:lib_name],
                          parsed_args[:fasta_dir];
                          params_path=params_path)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

# Entry point for PackageCompiler
function main_GetParseSpecLibParams(argv=ARGS)::Cint
    settings = ArgParseSettings()
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
    params_path = parsed_args[:"params-path"]
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
    
    # Read the JSON template and convert to OrderedDict
    config_text = read(asset_path("example_config", "defaultSearchParams.json"), String)
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
    @info "Writing default parameters .json to: $output_path"
    open(output_path, "w") do io
        JSON.print(io, config, 4)  # indent with 4 spaces for readability
    end
    
    return output_path
end

"""
    GetBuildLibParams(out_dir::String, lib_name::String, fasta_dir::String)

Creates a new library build parameter file with updated paths and automatically discovered FASTA files.
Uses a default template from assets/example_config/defaultBuildLibParams.json.

Arguments:
- out_dir: Output directory path
- lib_name: Library name path
- fasta_dir: Directory to search for FASTA files
- params_path: Path to folder or .json file in which to write the template parameters file. Defaults to joinpath(pwd(), "./buildspeclib_params.json")

Returns:
- String: Path to the newly created parameters file
"""
function GetBuildLibParams(out_dir::String, lib_name::String, fasta_dir::String; params_path::Union{String, Missing} = missing)
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

    # Parse JSON
    config_text = read(asset_path("example_config", "defaultBuildLibParams.json"), String)
    config = JSON.parse(config_text, dicttype=OrderedDict)

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
    @info "Writing default parameters .json to: $output_path"
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

    @info "Writing default parameters .json to: $output_path"
    open(output_path, "w") do io
        JSON.print(io, config, 4)
    end

    return output_path
end