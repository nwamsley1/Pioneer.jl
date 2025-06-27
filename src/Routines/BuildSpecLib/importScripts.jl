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

function importScriptsSpecLib(files_loaded::Set{String} = Set{String}())
    package_root = dirname(dirname(dirname(@__DIR__)))
    
    # Safe import system with shared files_loaded set
    total_attempts = 0
    successful_includes = 0
    conflicts_skipped = String[]
    
    function get_julia_files(dir::String)
        julia_files = String[]
        for (root, _, files) in walkdir(dir)
            for file in files
                if endswith(file, ".jl")
                    push!(julia_files, joinpath(root, file))
                end
            end
        end
        return julia_files
    end
    
    function safe_include!(files_loaded::Set{String}, file_path::String)::Bool
        total_attempts += 1
        
        # Check if file already loaded
        if file_path in files_loaded
            push!(conflicts_skipped, file_path)
            return false
        end
        
        # Check if file exists
        if !isfile(file_path)
            @warn "File not found: $file_path"
            return false
        end
        
        # Include the file
        try
            include(file_path)
            push!(files_loaded, file_path)
            successful_includes += 1
            return true
        catch e
            @warn "Failed to include $file_path: $e"
            return false
        end
    end
    
    function include_files!(files_loaded::Set{String}, file_dir::String, file_names::Vector{String})
        file_paths = [joinpath(file_dir, fname) for fname in file_names]
        successful_count = 0
        for fpath in file_paths
            if safe_include!(files_loaded, fpath)
                successful_count += 1
            end
        end
        return successful_count
    end
    
    function safe_include_directory!(files_loaded::Set{String}, dir::String)
        files = get_julia_files(dir)
        successful_count = 0
        for jfile in files
            if safe_include!(files_loaded, jfile)
                successful_count += 1
            end
        end
        return successful_count
    end
    
    root_path = joinpath(package_root, "src", "Routines", "BuildSpecLib")
    
    # Include KoinaStructs directory
    safe_include_directory!(files_loaded, joinpath(package_root, "src", "structs", "KoinaStructs"))
    
    # Include FileOperations (only if not already loaded)
    fileops_path = joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl")
    safe_include!(files_loaded, fileops_path)
    
    # FASTA processing
    safe_include!(files_loaded, joinpath(root_path, "structs", "mods.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fasta", "fasta_parser.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fasta", "fasta_digest.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fasta", "fasta_utils.jl"))
    
    # Fragment handling
    safe_include!(files_loaded, joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fragments", "fragment_parse.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fragments", "fragment_index.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fragments", "fragment_annotation.jl"))
    safe_include!(files_loaded, joinpath(root_path, "fragments", "fragment_predict.jl"))
    
    # Koina integration
    safe_include!(files_loaded, joinpath(root_path, "koina", "koina_api.jl"))
    safe_include!(files_loaded, joinpath(root_path, "koina", "koina_batch_prep.jl"))
    safe_include!(files_loaded, joinpath(root_path, "koina", "koina_batch_parse.jl"))
    
    # Utilities
    safe_include!(files_loaded, joinpath(root_path, "utils", "io.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "estimate_collision_ev.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "math.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "get_mz.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "parse_isotope_mods.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "check_params.jl"))
    
    # Structs 
    safe_include!(files_loaded, joinpath(root_path, "structs", "EmpiricalLibrary.jl"))
    safe_include!(files_loaded, joinpath(root_path, "utils", "parse_mods.jl"))
    
    # Library building
    safe_include!(files_loaded, joinpath(root_path, "build", "build_poin_lib.jl"))
    
    # Chronologer Methods 
    safe_include!(files_loaded, joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    safe_include!(files_loaded, joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    safe_include!(files_loaded, joinpath(root_path, "chronologer", "chronologer_parse.jl"))
    
    # Profiling
    safe_include!(files_loaded, joinpath(package_root, "src", "utils", "profile.jl"))

    return files_loaded
end