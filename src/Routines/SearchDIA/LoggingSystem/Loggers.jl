# Loggers.jl
# Factory functions for creating the three logger types

export create_console_logger, create_simplified_logger, create_full_logger

# Log level mappings (using integer values directly)
const VERBOSITY_LEVELS = Dict(
    :silent => LogLevel(2001),    # Error + 1 (only critical)
    :minimal => LogLevel(1000),   # Warn level
    :normal => LogLevel(-100),    # User info level to include user messages
    :verbose => LogLevel(-500),   # More detailed info
    :debug => LogLevel(-1000)     # Debug level
)

"""
    create_console_logger(config::LoggingConfig) -> AbstractLogger

Create a console logger with progress bar awareness and runtime-adjustable verbosity.
"""
function create_console_logger(config::LoggingConfig)
    level = get(VERBOSITY_LEVELS, config.console_level, Logging.Info)
    
    # Debug: Print the actual level being used
    println("DEBUG: Console logger level set to: $(level) (value: $(Int(level))) for config level: $(config.console_level)")
    
    # Create base console logger with custom formatting for colored prefixes
    console_logger = ConsoleLogger(
        stdout,
        level;  # This is the minimum level for ConsoleLogger
        meta_formatter = custom_meta_formatter,  # Use our custom formatter for colored prefixes
        show_limited = false,
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
    
    # Add custom formatting for user messages and clean kwargs
    console_logger = TransformerLogger(console_logger) do log
        # Check if this is a no-prefix message - handle it specially
        if get(log.kwargs, :no_prefix, false)
            # For no-prefix messages, create a simple logger that just prints
            # We bypass the normal ConsoleLogger formatting
            println(stdout, log.message)
            # Return a modified log that won't be processed further
            # Set level to below any threshold so it gets filtered out
            return merge(log, (level = LogLevel(-10000),))
        end
        
        # For normal messages, clean kwargs
        clean_kwargs = Dict{Symbol,Any}()
        for (k, v) in log.kwargs
            if k != :is_user_message && k != :_group && k != :progress_step && k != :progress_active && k != :no_prefix
                clean_kwargs[k] = v
            end
        end
        
        return merge(log, (kwargs = clean_kwargs,))
    end
    
    # ConsoleLogger already filters by level, so we don't need MinLevelLogger wrapper
    return console_logger
end

"""
    create_simplified_logger(filepath::String, config::LoggingConfig) -> AbstractLogger

Create a simplified file logger similar to current dual_println output.
"""
function create_simplified_logger(filepath::String, config::LoggingConfig)
    # Open file in append mode
    io = open(filepath, "a")
    
    # Create base file logger using FormatLogger
    file_logger = FormatLogger(io) do io, args
        # Simple format - just output the message
        println(io, args.message)
        flush(io)  # Ensure immediate write to disk
    end
    
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
    
    # Wrap with MinLevelLogger to respect minimum level
    return MinLevelLogger(file_logger, Logging.Info)
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
        flush(io)  # Ensure immediate write to disk
    end
    
    return MinLevelLogger(file_logger, LogLevel(-2000))  # Capture everything (Debug - 1000)
end

# Helper functions

"""
    custom_meta_formatter(level, _module, group, id, file, line)

Provides colored prefixes for different log levels, including custom Pioneer levels.
"""
function custom_meta_formatter(level, _module, group, id, file, line)
    level_value = Logging.LogLevel(level)
    
    # Handle standard Julia log levels
    if level_value == Logging.Error
        return (:error, "Error", "")
    elseif level_value == Logging.Warn
        return (:yellow, "Warning", "")
    elseif level_value == Logging.Info
        return (:cyan, "Info", "")  # Use cyan for standard Info too
    # Handle custom Pioneer log levels (using numeric values)
    elseif level_value == LogLevel(-100)  # USER_INFO_LEVEL
        return (:cyan, "Info", "")
    elseif level_value == LogLevel(900)  # USER_WARN_LEVEL
        return (:yellow, "Warning", "")
    elseif level_value == LogLevel(-1000) || level_value == LogLevel(-1100) || level_value == LogLevel(-1200)  # DEBUG levels
        return (:blue, "Debug", "")
    elseif level_value == LogLevel(-2000)  # TRACE_LEVEL
        return (:light_black, "Trace", "")
    else
        return (:normal, "", "")
    end
end

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