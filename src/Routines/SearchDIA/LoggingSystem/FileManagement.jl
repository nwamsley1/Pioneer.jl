# FileManagement.jl
# File management utilities for the logging system

using Dates

export ensure_log_directory, archive_old_logs, cleanup_logs
export get_log_paths, estimate_log_size

"""
    ensure_log_directory(output_dir::String) -> String

Ensure the log directory exists and return the path.
"""
function ensure_log_directory(output_dir::String)
    log_dir = joinpath(output_dir, "logs")
    mkpath(log_dir)
    return log_dir
end

"""
    get_log_paths(output_dir::String) -> NamedTuple

Get standard log file paths for the given output directory.
"""
function get_log_paths(output_dir::String)
    return (
        simple = joinpath(output_dir, "pioneer_search_log.txt"),
        full = joinpath(output_dir, "pioneer_debug.log"),
        warnings = joinpath(output_dir, "warnings.txt")
    )
end

"""
    archive_old_logs(output_dir::String; days_old::Int = 7)

Archive log files older than specified days.
"""
function archive_old_logs(output_dir::String; days_old::Int = 7)
    archive_dir = joinpath(output_dir, "archived_logs")
    
    # Find old log files
    cutoff_date = now() - Day(days_old)
    
    for file in readdir(output_dir, join=true)
        if endswith(file, ".log") || endswith(file, ".txt")
            if stat(file).mtime < datetime2unix(cutoff_date)
                # Create archive directory if needed
                if !isdir(archive_dir)
                    mkdir(archive_dir)
                end
                
                # Move to archive with timestamp
                timestamp = Dates.format(unix2datetime(stat(file).mtime), "yyyymmdd")
                archived_name = "$(timestamp)_$(basename(file))"
                mv(file, joinpath(archive_dir, archived_name), force=true)
            end
        end
    end
end

"""
    cleanup_logs(output_dir::String; keep_last::Int = 5)

Clean up old log files, keeping only the most recent ones.
"""
function cleanup_logs(output_dir::String; keep_last::Int = 5)
    # Group log files by type
    log_groups = Dict{String, Vector{String}}()
    
    for file in readdir(output_dir, join=true)
        if endswith(file, ".log") || endswith(file, ".txt")
            # Extract base name without timestamp
            base = replace(basename(file), r"_\d{8}_\d{6}" => "")
            if !haskey(log_groups, base)
                log_groups[base] = String[]
            end
            push!(log_groups[base], file)
        end
    end
    
    # Keep only the most recent files for each type
    for (base, files) in log_groups
        if length(files) > keep_last
            # Sort by modification time
            sorted_files = sort(files, by=f->stat(f).mtime, rev=true)
            
            # Remove old files
            for file in sorted_files[keep_last+1:end]
                rm(file, force=true)
            end
        end
    end
end

"""
    estimate_log_size(level::Symbol, num_operations::Int) -> Float64

Estimate the log file size in MB based on verbosity level and operations.
"""
function estimate_log_size(level::Symbol, num_operations::Int)
    # Rough estimates bytes per operation
    bytes_per_op = Dict(
        :silent => 10,
        :minimal => 50,
        :normal => 200,
        :verbose => 500,
        :debug => 2000,
        :trace => 5000
    )
    
    bytes = get(bytes_per_op, level, 200) * num_operations
    return bytes / (1024 * 1024)  # Convert to MB
end

"""
    create_log_summary(output_dir::String) -> String

Create a summary of all log files in the directory.
"""
function create_log_summary(output_dir::String)
    summary = IOBuffer()
    
    println(summary, "Log Files Summary")
    println(summary, "=" ^ 50)
    println(summary, "Directory: $(output_dir)")
    println(summary, "Generated: $(now())")
    println(summary)
    
    total_size = 0
    file_count = 0
    
    for file in readdir(output_dir, join=true)
        if endswith(file, ".log") || endswith(file, ".txt")
            size = filesize(file)
            total_size += size
            file_count += 1
            
            # Get file info
            mtime = unix2datetime(stat(file).mtime)
            size_mb = size / (1024 * 1024)
            
            println(summary, "$(basename(file)):")
            println(summary, "  Size: $(round(size_mb, digits=2)) MB")
            println(summary, "  Modified: $(mtime)")
            
            # Count lines for text files
            if endswith(file, ".txt") || endswith(file, ".log")
                try
                    line_count = countlines(file)
                    println(summary, "  Lines: $(line_count)")
                catch
                    # Skip if file can't be read
                end
            end
            println(summary)
        end
    end
    
    println(summary, "-" ^ 50)
    println(summary, "Total files: $(file_count)")
    println(summary, "Total size: $(round(total_size / (1024 * 1024), digits=2)) MB")
    
    return String(take!(summary))
end

"""
    compress_log_file(filepath::String)

Compress a log file using gzip.
"""
function compress_log_file(filepath::String)
    if !isfile(filepath)
        return
    end
    
    using CodecZlib
    
    compressed_path = filepath * ".gz"
    
    open(filepath, "r") do input
        open(compressed_path, "w") do output
            stream = GzipCompressorStream(output)
            write(stream, read(input))
            close(stream)
        end
    end
    
    # Remove original file after successful compression
    if isfile(compressed_path)
        rm(filepath)
    end
end