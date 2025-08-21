# AsyncFileLogger.jl
# Asynchronous file logger implementation to prevent I/O blocking

export AsyncFileLogger, create_async_file_logger, close_async_logger, LogMessage

const DEFAULT_CHANNEL_SIZE = 10000
const FLUSH_INTERVAL = 1.0  # seconds
const FLUSH_MESSAGE_COUNT = 100

# Message struct to pass through channel
struct LogMessage
    level::LogLevel
    message::String
    _module::Module
    _group::Symbol
    id::Symbol
    file::String
    line::Int
    kwargs::Dict{Symbol,Any}
    timestamp::DateTime
end

# Async file logger implementation
mutable struct AsyncFileLogger <: AbstractLogger
    filepath::String
    io::Union{IOStream, Nothing}
    message_channel::Channel{LogMessage}
    writer_task::Union{Task, Nothing}
    min_level::LogLevel
    is_active::Ref{Bool}
    formatter::Function
    last_flush::Ref{Float64}
    message_count::Ref{Int}
    channel_size::Int
end

"""
    create_async_file_logger(filepath, min_level; kwargs...)

Create an asynchronous file logger that writes to a background task.
"""
function create_async_file_logger(
    filepath::String, 
    min_level::LogLevel;
    formatter::Function = default_formatter,
    channel_size::Int = DEFAULT_CHANNEL_SIZE
)
    # Create channel for messages
    message_channel = Channel{LogMessage}(channel_size)
    
    # Create logger instance (don't open file yet)
    logger = AsyncFileLogger(
        filepath,
        nothing,  # io - will be opened in writer task
        message_channel,
        nothing,  # task - will be started separately
        min_level,
        Ref(true),
        formatter,
        Ref(time()),
        Ref(0),
        channel_size
    )
    
    # Start the async writer task
    logger.writer_task = start_writer_task(logger)
    
    return logger
end

# Default formatter for messages
function default_formatter(msg::LogMessage)
    timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS")
    return "[$(timestamp)] [$(msg.level)] $(msg.message)"
end

# Start the background writer task
function start_writer_task(logger::AsyncFileLogger)
    return @async begin
        try
            # Open file in the writer task (avoids handle conflicts)
            logger.io = open(logger.filepath, "a")
            
            while logger.is_active[]
                try
                    # Wait for message
                    msg = take!(logger.message_channel)
                    
                    # Write formatted message
                    write_message(logger, msg)
                    
                    # Check if we should flush
                    maybe_flush(logger)
                    
                catch e
                    if e isa InvalidStateException && !logger.is_active[]
                        break  # Channel closed, expected during shutdown
                    else
                        # Log error but continue (don't crash the writer)
                        println(stderr, "Async logger error: ", e)
                    end
                end
            end
            
        catch e
            println(stderr, "Fatal async logger error: ", e)
        finally
            # Drain remaining messages
            drain_channel(logger)
            
            # Close file
            if logger.io !== nothing
                flush(logger.io)
                close(logger.io)
                logger.io = nothing
            end
        end
    end
end

# Write a message to the file
function write_message(logger::AsyncFileLogger, msg::LogMessage)
    if logger.io !== nothing
        formatted = logger.formatter(msg)
        println(logger.io, formatted)
        logger.message_count[] += 1
    end
end

# Check if we should flush the file
function maybe_flush(logger::AsyncFileLogger)
    current_time = time()
    time_since_flush = current_time - logger.last_flush[]
    
    if time_since_flush >= FLUSH_INTERVAL || logger.message_count[] >= FLUSH_MESSAGE_COUNT
        if logger.io !== nothing
            flush(logger.io)
        end
        logger.last_flush[] = current_time
        logger.message_count[] = 0
    end
end

# Drain any remaining messages in the channel
function drain_channel(logger::AsyncFileLogger)
    while isready(logger.message_channel)
        try
            msg = take!(logger.message_channel)
            write_message(logger, msg)
        catch
            break
        end
    end
    if logger.io !== nothing
        flush(logger.io)
    end
end

# Implement AbstractLogger interface
function Logging.shouldlog(logger::AsyncFileLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.min_enabled_level(logger::AsyncFileLogger)
    return logger.min_level
end

function Logging.handle_message(
    logger::AsyncFileLogger, level, message, _module, group, id, file, line; kwargs...
)
    # Create LogMessage struct
    msg = LogMessage(
        level, string(message), _module, group, id, file, line,
        Dict{Symbol,Any}(kwargs), now()
    )
    
    # Try to put message in channel (non-blocking)
    if !put_nowait(logger.message_channel, msg)
        # Channel full - try brief blocking, then fallback
        if !put_with_timeout(logger.message_channel, msg, 0.01)  # 10ms timeout
            # Last resort: write to stderr
            println(stderr, "[DROPPED LOG] ", logger.formatter(msg))
        end
    end
end

# Non-blocking put to channel
function put_nowait(channel::Channel{T}, item::T) where T
    try
        if length(channel.data) < channel.sz_max
            put!(channel, item)
            return true
        end
    catch
        # Channel might be closed or full
    end
    return false
end

# Put with timeout
function put_with_timeout(channel::Channel{T}, item::T, timeout::Float64) where T
    start_time = time()
    while time() - start_time < timeout
        if put_nowait(channel, item)
            return true
        end
        yield()  # Give other tasks a chance
    end
    return false
end

# Properly close the async logger
function close_async_logger(logger::AsyncFileLogger)
    # Signal shutdown
    logger.is_active[] = false
    
    # Close channel to wake up writer
    try
        close(logger.message_channel)
    catch
        # Channel might already be closed
    end
    
    # Wait for writer to finish (with timeout)
    if logger.writer_task !== nothing
        try
            wait_with_timeout(logger.writer_task, 5.0)
        catch e
            @warn "Async logger shutdown timeout" error=e
        end
    end
end

# Wait for task with timeout
function wait_with_timeout(task::Task, timeout::Float64)
    start_time = time()
    while !istaskdone(task) && time() - start_time < timeout
        sleep(0.01)
    end
    if !istaskdone(task)
        @warn "Task did not complete within timeout"
    end
end

# Get statistics about the logger
function get_logger_stats(logger::AsyncFileLogger)
    return (
        channel_size = logger.channel_size,
        messages_queued = isopen(logger.message_channel) ? length(logger.message_channel.data) : 0,
        is_active = logger.is_active[],
        task_state = logger.writer_task !== nothing ? logger.writer_task.state : :not_started,
        file_open = logger.io !== nothing
    )
end