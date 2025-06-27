#!/usr/bin/env julia

# Script to remove coverage-related files from git tracking

using Base: shell_escape

# File extensions to remove
const EXTENSIONS_TO_REMOVE = [".cov", ".gcov", ".html"]

# Function to find all files with specified extensions
function find_files_to_remove(root_dir::String)
    files_to_remove = String[]
    
    for (root, dirs, files) in walkdir(root_dir)
        # Skip .git directory
        if ".git" in dirs
            deleteat!(dirs, findfirst(==(".git"), dirs))
        end
        
        for file in files
            for ext in EXTENSIONS_TO_REMOVE
                if endswith(file, ext)
                    push!(files_to_remove, joinpath(root, file))
                    break
                end
            end
        end
    end
    
    return files_to_remove
end

# Function to check if a file is tracked by git
function is_git_tracked(filepath::String)
    try
        # Use git ls-files to check if file is tracked
        cmd = `git ls-files --error-unmatch $(filepath)`
        run(pipeline(cmd, stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end

# Main function
function main(; dry_run::Bool=false)
    repo_root = pwd()
    
    println("Searching for coverage files in: $repo_root")
    println("Looking for extensions: ", join(EXTENSIONS_TO_REMOVE, ", "))
    println("Dry run mode: $dry_run")
    println()
    
    # Find all files to remove
    files_to_remove = find_files_to_remove(repo_root)
    
    if isempty(files_to_remove)
        println("No files found with extensions: ", join(EXTENSIONS_TO_REMOVE, ", "))
        return
    end
    
    println("Found $(length(files_to_remove)) files to process")
    println()
    
    # Separate tracked and untracked files
    tracked_files = String[]
    untracked_files = String[]
    
    for filepath in files_to_remove
        if is_git_tracked(filepath)
            push!(tracked_files, filepath)
        else
            push!(untracked_files, filepath)
        end
    end
    
    # Process tracked files with git rm
    if !isempty(tracked_files)
        println("Files tracked by git (will use 'git rm'):")
        for file in tracked_files
            rel_path = relpath(file, repo_root)
            println("  - $rel_path")
        end
        
        if !dry_run
            println()
            print("Removing tracked files from git... ")
            flush(stdout)
            
            # Use git rm for tracked files
            for file in tracked_files
                try
                    run(`git rm -f $(file)`)
                catch e
                    println("\nError removing $file: $e")
                end
            end
            println("done!")
        end
    end
    
    # Process untracked files with regular rm
    if !isempty(untracked_files)
        println()
        println("Files not tracked by git (will use regular 'rm'):")
        for file in untracked_files
            rel_path = relpath(file, repo_root)
            println("  - $rel_path")
        end
        
        if !dry_run
            println()
            print("Removing untracked files... ")
            flush(stdout)
            
            for file in untracked_files
                try
                    rm(file, force=true)
                catch e
                    println("\nError removing $file: $e")
                end
            end
            println("done!")
        end
    end
    
    # Add patterns to .gitignore
    gitignore_path = joinpath(repo_root, ".gitignore")
    gitignore_patterns = ["*.cov", "*.gcov", "*.html"]
    
    println()
    println("Checking .gitignore...")
    
    existing_patterns = String[]
    if isfile(gitignore_path)
        existing_patterns = readlines(gitignore_path)
    end
    
    patterns_to_add = String[]
    for pattern in gitignore_patterns
        if !(pattern in existing_patterns)
            push!(patterns_to_add, pattern)
        end
    end
    
    if !isempty(patterns_to_add)
        println("Patterns to add to .gitignore:")
        for pattern in patterns_to_add
            println("  - $pattern")
        end
        
        if !dry_run
            open(gitignore_path, "a") do io
                if !isempty(existing_patterns) && !isempty(strip(existing_patterns[end]))
                    println(io)  # Add blank line if needed
                end
                println(io, "\n# Coverage files")
                for pattern in patterns_to_add
                    println(io, pattern)
                end
            end
            println("Updated .gitignore")
        end
    else
        println("All patterns already in .gitignore")
    end
    
    # Print summary
    println()
    println("Summary:")
    println("  Tracked files $(dry_run ? "to be" : "") removed: $(length(tracked_files))")
    println("  Untracked files $(dry_run ? "to be" : "") removed: $(length(untracked_files))")
    println("  Total files processed: $(length(files_to_remove))")
    
    if dry_run && !isempty(files_to_remove)
        println()
        println("This was a dry run. To actually remove the files, run:")
        println("  julia remove_coverage_files.jl --no-dry-run")
    elseif !dry_run && !isempty(tracked_files)
        println()
        println("Don't forget to commit the changes:")
        println("  git commit -m \"Remove coverage files from tracking\"")
    end
end

# Parse command line arguments
if abspath(PROGRAM_FILE) == @__FILE__
    dry_run = true
    
    if "--no-dry-run" in ARGS
        dry_run = false
    elseif "--help" in ARGS || "-h" in ARGS
        println("Usage: julia remove_coverage_files.jl [--no-dry-run]")
        println()
        println("Removes .cov, .gcov, and .html files from the repository.")
        println("Uses 'git rm' for tracked files and regular 'rm' for untracked files.")
        println()
        println("Options:")
        println("  --no-dry-run    Actually remove the files (default is dry run)")
        println("  --help, -h      Show this help message")
        exit(0)
    end
    
    main(dry_run=dry_run)
end