# SimpleFileLogger.jl
# Simple thread-safe file logger that writes directly to file (like dual_println)

export SimpleFileLogger, create_simple_file_logger, close_simple_logger

# Simple file logger with thread-safe direct writing
mutable struct SimpleFileLogger <: AbstractLogger
    filepath::String
    io::Union{IOStream, Nothing}
    lock::ReentrantLock
    min_level::LogLevel
    formatter::Function
    last_flush_time::Ref{Float64}
    flush_interval::Float64  # seconds between flushes
end

"""
    create_simple_file_logger(filepath, min_level; formatter, flush_interval)

Create a simple thread-safe file logger that writes directly to file.
Similar to the original dual_println approach but integrated with logging system.
"""
function create_simple_file_logger(
    filepath::String, 
    min_level::LogLevel;
    formatter::Function = default_formatter,
    flush_interval::Float64 = 1.0  # Flush every second by default
)
    # Open file for appending
    io = open(filepath, "a")
    
    return SimpleFileLogger(
        filepath,
        io,
        ReentrantLock(),
        min_level,
        formatter,
        Ref(time()),  # Initialize last flush time
        flush_interval
    )
end

# Default formatter for messages
function default_formatter(msg::NamedTuple)
    timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS")
    return "[$(timestamp)] [$(msg.level)] $(msg.message)"
end

# Implement AbstractLogger interface
function Logging.shouldlog(logger::SimpleFileLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.min_enabled_level(logger::SimpleFileLogger)
    return logger.min_level
end

function Logging.handle_message(
    logger::SimpleFileLogger, level, message, _module, group, id, file, line; kwargs...
)
    # Only process if we should log this level
    if level >= logger.min_level && logger.io !== nothing
        # Create message data
        msg = (
            level = level,
            message = string(message),
            _module = _module,
            _group = group,
            id = id,
            file = file,
            line = line,
            kwargs = kwargs,
            timestamp = now()
        )
        
        # Format the message
        formatted = logger.formatter(msg)
        
        # Thread-safe write WITHOUT immediate flush (prevents Arrow.Table deadlock)
        lock(logger.lock) do
            try
                println(logger.io, formatted)
                
                # Only flush periodically to avoid blocking I/O
                current_time = time()
                if current_time - logger.last_flush_time[] > logger.flush_interval
                    flush(logger.io)
                    logger.last_flush_time[] = current_time
                end
            catch e
                # If write fails, print to stderr but don't crash
                println(stderr, "Failed to write to log file: ", e)
            end
        end
    end
end

# Properly close the logger
function close_simple_logger(logger::SimpleFileLogger)
    lock(logger.lock) do
        if logger.io !== nothing
            try
                flush(logger.io)
                close(logger.io)
            catch
                # Ignore errors on close
            end
            logger.io = nothing
        end
    end
end

# Formatter for simplified log output
function simple_formatter(msg::NamedTuple)
    timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS")
    return "[$(timestamp)] $(msg.message)"
end

# Formatter for full debug output
function full_formatter(msg::NamedTuple)
    timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS.sss")
    level_str = string(msg.level)
    module_str = string(msg._module)
    file_line = "$(basename(string(msg.file))):$(msg.line)"
    
    output = IOBuffer()
    println(output, "[$(timestamp)] [$(level_str)] [$(module_str)] [$(file_line)]")
    println(output, "  Message: ", msg.message)
    
    # Add all kwargs
    if !isempty(msg.kwargs)
        println(output, "  Context:")
        for (k, v) in msg.kwargs
            if k != :is_user_message && k != :progress_active && k != :no_prefix
                println(output, "    $(k): $(v)")
            end
        end
    end
    
    println(output, "  " * "-"^60)
    return String(take!(output))
end