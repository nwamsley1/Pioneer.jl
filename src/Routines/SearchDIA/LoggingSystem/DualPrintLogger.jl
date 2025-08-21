# DualPrintLogger.jl
# Minimal logger that exactly mimics the original dual_println behavior
# No locks, no flush, just simple println to avoid Arrow.Table conflicts

export DualPrintLogger, create_dual_print_logger, close_dual_print_logger

# Minimal logger that writes directly to file like dual_println
# Uses lazy file opening - opens and closes file for each write to avoid Arrow.Table conflicts
mutable struct DualPrintLogger <: AbstractLogger
    filepath::String  # Store filepath instead of open handle
    min_level::LogLevel
    formatter::Function
end

"""
    create_dual_print_logger(filepath, min_level; formatter)

Create a minimal logger that mimics dual_println behavior.
Opens and closes file for each write to avoid conflicts with Arrow.Table async operations.
"""
function create_dual_print_logger(
    filepath::String, 
    min_level::LogLevel;
    formatter::Function = simple_formatter
)
    # Don't open file here - will open/close for each write
    return DualPrintLogger(filepath, min_level, formatter)
end

# Simple formatter - just timestamp and message
function simple_formatter(level, message, timestamp)
    ts = Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS")
    return "[$(ts)] $(message)"
end

# Full formatter with metadata
function full_formatter(level, message, timestamp, _module, file, line, kwargs)
    ts = Dates.format(timestamp, "yyyy-mm-dd HH:MM:SS.sss")
    level_str = string(level)
    module_str = string(_module)
    file_line = "$(basename(string(file))):$(line)"
    
    output = IOBuffer()
    println(output, "[$(ts)] [$(level_str)] [$(module_str)] [$(file_line)]")
    println(output, "  Message: ", message)
    
    if !isempty(kwargs)
        println(output, "  Context:")
        for (k, v) in kwargs
            if k != :is_user_message && k != :progress_active && k != :no_prefix
                println(output, "    $(k): $(v)")
            end
        end
    end
    
    println(output, "  " * "-"^60)
    return String(take!(output))
end

# Implement AbstractLogger interface
function Logging.shouldlog(logger::DualPrintLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.min_enabled_level(logger::DualPrintLogger)
    return logger.min_level
end

function Logging.handle_message(
    logger::DualPrintLogger, level, message, _module, group, id, file, line; kwargs...
)
    if level >= logger.min_level
        # Format based on logger type
        formatted = if logger.formatter === simple_formatter
            logger.formatter(level, string(message), now())
        else
            # Full formatter needs more args
            logger.formatter(level, string(message), now(), _module, file, line, kwargs)
        end
        
        # Open, write, close - exactly like dual_println did
        # This avoids keeping file handles open during Arrow.Table operations
        try
            open(logger.filepath, "a") do io
                println(io, formatted)
            end
        catch e
            # If file write fails, don't crash - just skip this message
            # This matches dual_println's resilient behavior
        end
    end
end

# Close the logger
function close_dual_print_logger(logger::DualPrintLogger)
    # Nothing to close since we open/close for each write
    # This function exists for API compatibility
    return nothing
end