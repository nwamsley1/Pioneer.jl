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

#!/usr/bin/env julia

# Script to add copyright headers to all .jl files in the Pioneer.jl repository

const COPYRIGHT_HEADER = """
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

"""

# Directories to skip
const SKIP_DIRS = [".git", "deps", "build", ".julia", "docs/build"]

# Check if a file already has a copyright header
function has_copyright_header(filepath::String)
    try
        content = read(filepath, String)
        # Check if the file starts with copyright information
        return occursin(r"^#\s*Copyright.*Nathan\s*Wamsley"i, content) ||
               occursin(r"^#\s*This file is part of Pioneer\.jl"i, content)
    catch
        println("Warning: Could not read file $filepath")
        return true  # Skip files we can't read
    end
end

# Add copyright header to a file
function add_header_to_file(filepath::String; dry_run::Bool=false)
    if has_copyright_header(filepath)
        return false
    end
    
    try
        original_content = read(filepath, String)
        new_content = COPYRIGHT_HEADER * original_content
        
        if !dry_run
            write(filepath, new_content)
        end
        
        return true
    catch e
        println("Error processing file $filepath: $e")
        return false
    end
end

# Find all .jl files in the repository
function find_julia_files(root_dir::String)
    julia_files = String[]
    
    for (root, dirs, files) in walkdir(root_dir)
        # Remove directories we want to skip from the traversal
        for skip_dir in SKIP_DIRS
            if skip_dir in dirs
                deleteat!(dirs, findfirst(==(skip_dir), dirs))
            end
        end
        
        # Skip if we're in a directory that contains any of the skip patterns
        should_skip = false
        for skip_pattern in SKIP_DIRS
            if occursin(skip_pattern, root)
                should_skip = true
                break
            end
        end
        
        if should_skip
            continue
        end
        
        # Add .jl files to our list
        for file in files
            if endswith(file, ".jl")
                push!(julia_files, joinpath(root, file))
            end
        end
    end
    
    return julia_files
end

# Main function
function main(; dry_run::Bool=false)
    # Get the repository root (assuming script is run from repo root)
    repo_root = pwd()
    
    println("Searching for Julia files in: $repo_root")
    println("Dry run mode: $dry_run")
    println()
    
    # Find all Julia files
    julia_files = find_julia_files(repo_root)
    println("Found $(length(julia_files)) Julia files")
    
    # Process each file
    modified_count = 0
    skipped_count = 0
    
    for filepath in julia_files
        rel_path = relpath(filepath, repo_root)
        
        if add_header_to_file(filepath; dry_run=dry_run)
            println("✓ $(dry_run ? "Would add" : "Added") header to: $rel_path")
            modified_count += 1
        else
            skipped_count += 1
        end
    end
    
    # Print summary
    println()
    println("Summary:")
    println("  Files that $(dry_run ? "would be" : "were") modified: $modified_count")
    println("  Files skipped (already have header): $skipped_count")
    println("  Total files processed: $(length(julia_files))")
    
    if dry_run && modified_count > 0
        println()
        println("This was a dry run. To actually modify the files, run:")
        println("  julia add_copyright_headers.jl --no-dry-run")
    end
end

# Parse command line arguments
if abspath(PROGRAM_FILE) == @__FILE__
    dry_run = true
    
    if "--no-dry-run" in ARGS
        dry_run = false
    elseif "--help" in ARGS || "-h" in ARGS
        println("Usage: julia add_copyright_headers.jl [--no-dry-run]")
        println()
        println("Options:")
        println("  --no-dry-run    Actually modify the files (default is dry run)")
        println("  --help, -h      Show this help message")
        exit(0)
    end
    
    main(dry_run=dry_run)
end