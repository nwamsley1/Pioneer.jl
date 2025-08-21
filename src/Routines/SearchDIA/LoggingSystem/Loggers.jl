# Loggers.jl
# Factory functions for creating the three logger types

export create_console_logger, create_simplified_logger, create_full_logger

# Log level mappings
const VERBOSITY_LEVELS = Dict(
    :silent => Logging.Error + 1,  # Only critical errors
    :minimal => Logging.Warn,       # Warnings and errors
    :normal => Logging.Info,        # Standard user messages
    :verbose => Logging.Info - 100, # Include @user_info and similar
    :debug => Logging.Debug         # All debug messages
)

"""
    create_console_logger(config::LoggingConfig) -> AbstractLogger

Create a console logger with progress bar awareness and runtime-adjustable verbosity.
"""
function create_console_logger(config::LoggingConfig)
    level = get(VERBOSITY_LEVELS, config.console_level, Logging.Info)
    
    # Create base console logger
    console_logger = ConsoleLogger(
        stdout,
        level;
        meta_formatter = console_meta_formatter,
        show_limited = true,
        right_justify = 0
    )
    
    # Add progress-aware formatting
    if config.enable_progress
        console_logger = TransformerLogger(console_logger) do log
            # Check if progress bar is active (we'll set this in context)
            if get(log.kwargs, :progress_active, false)
                # Suppress or queue message when progress bar is active
                return nothing
            end
            return log
        end
    end
    
    # Add custom formatting for user messages
    console_logger = TransformerLogger(console_logger) do log
        # Format based on custom log level
        if haskey(log.kwargs, :is_user_message) && log.kwargs[:is_user_message]
            # Clean user message formatting
            formatted_msg = format_user_message(log.level, log.message, log.kwargs)
            return merge(log, (message = formatted_msg,))
        end
        return log
    end
    
    return MinLevelLogger(console_logger, level)
end

"""
    create_simplified_logger(filepath::String, config::LoggingConfig) -> AbstractLogger

Create a simplified file logger similar to current dual_println output.
"""
function create_simplified_logger(filepath::String, config::LoggingConfig)
    # Open file in append mode
    io = open(filepath, "a")
    
    # Create base file logger
    file_logger = SimpleLogger(
        io,
        Logging.Info;  # Only user-facing messages
        format = simple_format
    )
    
    # Add timestamp formatting
    file_logger = TransformerLogger(file_logger) do log
        # Add timestamp
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        formatted_msg = "[$(timestamp)] $(log.message)"
        return merge(log, (message = formatted_msg,))
    end
    
    # Filter to only user messages and warnings
    file_logger = EarlyFilteredLogger(file_logger) do log
        return log.level >= Logging.Info || 
               (haskey(log.kwargs, :is_user_message) && log.kwargs[:is_user_message])
    end
    
    return file_logger
end

"""
    create_full_logger(filepath::String, config::LoggingConfig) -> AbstractLogger

Create a comprehensive debug file logger with all metadata.
"""
function create_full_logger(filepath::String, config::LoggingConfig)
    # Check for rotation needs
    if should_rotate(filepath, config.rotation_size_mb)
        rotate_log_file(filepath)
    end
    
    # Open file in append mode
    io = open(filepath, "a")
    
    # Create comprehensive file logger
    file_logger = FormatLogger(io) do io, args
        # Full format with all metadata
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        level_str = string(args.level)
        module_str = string(args._module)
        file_line = "$(basename(string(args.file))):$(args.line)"
        
        # Format main message
        println(io, "[$(timestamp)] [$(level_str)] [$(module_str)] [$(file_line)]")
        println(io, "  Message: ", args.message)
        
        # Add all kwargs
        if !isempty(args.kwargs)
            println(io, "  Context:")
            for (k, v) in args.kwargs
                if k != :is_user_message && k != :progress_active
                    println(io, "    $(k): $(v)")
                end
            end
        end
        
        # Add separator for readability
        println(io, "  " * "-"^60)
    end
    
    return MinLevelLogger(file_logger, Logging.Debug - 1000)  # Capture everything
end

# Helper functions

function console_meta_formatter(level, _module, group, id, file, line)
    color = Logging.default_logcolor(level)
    prefix = level == Logging.Warn ? "WARNING" : 
             level == Logging.Error ? "ERROR" : 
             level == Logging.Debug ? "DEBUG" : ""
    
    if !isempty(prefix)
        return color, prefix, ""
    else
        return :normal, "", ""
    end
end

function simple_format(io, args)
    # Clean, simple format for the simplified log
    println(io, args.message)
end

function format_user_message(level, message, kwargs)
    # Format user-facing messages cleanly
    if level == Logging.Warn
        return "⚠️  $(message)"
    elseif level == Logging.Error
        return "❌ $(message)"
    elseif haskey(kwargs, :progress_step)
        return "▶ $(message)"
    else
        return message
    end
end

function should_rotate(filepath::String, max_size_mb::Float64)
    if !isfile(filepath)
        return false
    end
    
    file_size_mb = filesize(filepath) / (1024 * 1024)
    return file_size_mb >= max_size_mb
end

function rotate_log_file(filepath::String)
    if !isfile(filepath)
        return
    end
    
    # Generate rotation name with timestamp
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    dir = dirname(filepath)
    base = basename(filepath)
    name, ext = splitext(base)
    
    rotated_path = joinpath(dir, "$(name)_$(timestamp)$(ext)")
    
    # Move current file
    mv(filepath, rotated_path, force=true)
    
    # Keep only last 5 rotated files
    pattern = Regex("^$(name)_\\d{8}_\\d{6}$(ext)\$")
    rotated_files = filter(f -> occursin(pattern, f), readdir(dir))
    
    if length(rotated_files) > 5
        # Sort by timestamp and remove oldest
        sort!(rotated_files)
        for f in rotated_files[1:end-5]
            rm(joinpath(dir, f), force=true)
        end
    end
end