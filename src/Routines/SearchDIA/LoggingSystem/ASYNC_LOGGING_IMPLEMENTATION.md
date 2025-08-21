# Asynchronous File Logging Implementation Plan

## Overview

Implement non-blocking, asynchronous file logging to prevent I/O conflicts and hanging issues. This solution addresses the root cause of file I/O blocking rather than working around it.

## Core Architecture

### 1. Message Queue System

```julia
mutable struct AsyncFileLogger <: AbstractLogger
    filepath::String
    io::IOStream
    message_channel::Channel{LogMessage}
    writer_task::Task
    min_level::LogLevel
    is_active::Ref{Bool}
    formatter::Function
end

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
```

### 2. Async Writer Task

The writer task runs in the background, consuming messages from the channel:

```julia
function start_async_writer(logger::AsyncFileLogger)
    return @async begin
        try
            while logger.is_active[]
                # Wait for message with timeout to allow graceful shutdown
                msg = take!(logger.message_channel)
                
                # Format and write message
                formatted = logger.formatter(msg)
                println(logger.io, formatted)
                
                # Periodic flush (every N messages or time interval)
                if should_flush(logger)
                    flush(logger.io)
                end
            end
        catch e
            if !(e isa InvalidStateException)  # Channel closed
                @error "Async logger writer failed" error=e
            end
        finally
            # Ensure remaining messages are written
            while isready(logger.message_channel)
                msg = take!(logger.message_channel)
                println(logger.io, logger.formatter(msg))
            end
            flush(logger.io)
            close(logger.io)
        end
    end
end
```

## Implementation Details

### Phase 1: Core Async Logger

**File: `AsyncFileLogger.jl`**

```julia
# AsyncFileLogger.jl
export AsyncFileLogger, create_async_file_logger

const DEFAULT_CHANNEL_SIZE = 10000
const FLUSH_INTERVAL = 1.0  # seconds
const FLUSH_MESSAGE_COUNT = 100

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

function start_writer_task(logger::AsyncFileLogger)
    return @async begin
        try
            # Open file in the writer task (avoids handle conflicts)
            logger.io = open(logger.filepath, "a")
            
            while logger.is_active[]
                try
                    # Wait for message with timeout
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

function write_message(logger::AsyncFileLogger, msg::LogMessage)
    formatted = logger.formatter(msg)
    println(logger.io, formatted)
    logger.message_count[] += 1
end

function maybe_flush(logger::AsyncFileLogger)
    current_time = time()
    time_since_flush = current_time - logger.last_flush[]
    
    if time_since_flush >= FLUSH_INTERVAL || logger.message_count[] >= FLUSH_MESSAGE_COUNT
        flush(logger.io)
        logger.last_flush[] = current_time
        logger.message_count[] = 0
    end
end

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
        # Channel full - we have options:
        # 1. Drop message (fastest, may lose logs)
        # 2. Block briefly with timeout
        # 3. Write directly to stderr as fallback
        
        # Option 2: Brief blocking with timeout
        try
            put_with_timeout(logger.message_channel, msg, 0.01)  # 10ms timeout
        catch
            # Option 3: Fallback to stderr
            println(stderr, "[DROPPED] ", logger.formatter(msg))
        end
    end
end

function put_nowait(channel::Channel{T}, item::T) where T
    if isready(channel) && length(channel.data) < channel.sz_max
        put!(channel, item)
        return true
    end
    return false
end

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

function close_async_logger(logger::AsyncFileLogger)
    # Signal shutdown
    logger.is_active[] = false
    
    # Close channel to wake up writer
    close(logger.message_channel)
    
    # Wait for writer to finish (with timeout)
    try
        wait_with_timeout(logger.writer_task, 5.0)
    catch
        @warn "Async logger shutdown timeout"
    end
end

function wait_with_timeout(task::Task, timeout::Float64)
    start_time = time()
    while !istaskdone(task) && time() - start_time < timeout
        sleep(0.01)
    end
    if !istaskdone(task)
        throw(TimeoutError("Task did not complete within timeout"))
    end
end
```

### Phase 2: Integration with Existing Loggers

**Update `Loggers.jl`:**

```julia
function create_simplified_logger(filepath::String, config::LoggingConfig)
    # Use async logger instead of FormatLogger
    formatter = (msg) -> begin
        timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS")
        "[$(timestamp)] $(msg.message)"
    end
    
    logger = create_async_file_logger(
        filepath,
        Logging.Info;
        formatter = formatter,
        channel_size = 10000
    )
    
    # Add filtering for user messages
    return EarlyFilteredLogger(logger) do log
        return log.level >= Logging.Info || 
               (haskey(log.kwargs, :is_user_message) && log.kwargs[:is_user_message])
    end
end

function create_full_logger(filepath::String, config::LoggingConfig)
    formatter = (msg) -> begin
        timestamp = Dates.format(msg.timestamp, "yyyy-mm-dd HH:MM:SS.sss")
        level_str = string(msg.level)
        module_str = string(msg._module)
        file_line = "$(basename(string(msg.file))):$(msg.line)"
        
        output = IOBuffer()
        println(output, "[$(timestamp)] [$(level_str)] [$(module_str)] [$(file_line)]")
        println(output, "  Message: ", msg.message)
        
        if !isempty(msg.kwargs)
            println(output, "  Context:")
            for (k, v) in msg.kwargs
                if k != :is_user_message && k != :progress_active
                    println(output, "    $(k): $(v)")
                end
            end
        end
        
        println(output, "  " * "-"^60)
        String(take!(output))
    end
    
    return create_async_file_logger(
        filepath,
        LogLevel(-2000);
        formatter = formatter,
        channel_size = 50000  # Larger for debug logger
    )
end
```

### Phase 3: Graceful Shutdown

**Update `LoggingSystem.jl`:**

```julia
struct LoggerState
    warning_logger::WarningCapturingLogger
    console_logger::Ref{AbstractLogger}
    simplified_logger::Union{Nothing, AsyncFileLogger}  # Now AsyncFileLogger
    full_logger::Union{Nothing, AsyncFileLogger}        # Now AsyncFileLogger
    config::LoggingConfig
    start_time::DateTime
end

function finalize_logging()
    state = LOGGER_STATE[]
    if state === nothing
        return
    end
    
    # ... existing warning handling ...
    
    # Close async file loggers with proper shutdown
    if state.simplified_logger !== nothing && isa(state.simplified_logger, AsyncFileLogger)
        close_async_logger(state.simplified_logger)
    end
    if state.full_logger !== nothing && isa(state.full_logger, AsyncFileLogger)
        close_async_logger(state.full_logger)
    end
    
    # Reset global state
    LOGGER_STATE[] = nothing
    global_logger(ConsoleLogger())
end
```

## Benefits of This Approach

### 1. Non-Blocking I/O
- Main thread never blocks on file writes
- Library loading can proceed without I/O conflicts
- Better performance for high-frequency logging

### 2. Graceful Degradation
- If channel fills up, messages can be dropped or redirected
- System continues running even if logging has issues
- Configurable overflow behavior

### 3. Better Resource Management
- Single file handle per logger (no contention)
- Controlled flushing reduces I/O overhead
- Proper cleanup on shutdown

### 4. Future-Proof
- Solves current hanging issue
- Prevents future I/O conflicts
- Scalable for high-throughput scenarios

## Testing Strategy

### 1. Unit Tests

```julia
# test/test_async_logger.jl
@testset "AsyncFileLogger" begin
    # Test basic functionality
    @test begin
        logger = create_async_file_logger("test.log", Logging.Info)
        Logging.handle_message(logger, Logging.Info, "Test message", @__MODULE__, :test, :id, "file.jl", 1)
        close_async_logger(logger)
        
        # Check file was written
        content = read("test.log", String)
        occursin("Test message", content)
    end
    
    # Test high-throughput
    @test begin
        logger = create_async_file_logger("test2.log", Logging.Debug)
        for i in 1:10000
            Logging.handle_message(logger, Logging.Info, "Message $i", @__MODULE__, :test, :id, "file.jl", i)
        end
        close_async_logger(logger)
        
        # All messages should be written
        content = read("test2.log", String)
        count(occursin("Message", content)) == 10000
    end
    
    # Test overflow handling
    @test begin
        logger = create_async_file_logger("test3.log", Logging.Info; channel_size=10)
        # Flood with messages
        for i in 1:100
            Logging.handle_message(logger, Logging.Info, "Flood $i", @__MODULE__, :test, :id, "file.jl", i)
        end
        close_async_logger(logger)
        true  # Should not crash
    end
end
```

### 2. Integration Test

```julia
# Test with actual SearchDIA
@test begin
    # Should not hang at library loading
    SearchDIA("test_params.json")
    
    # Check log files exist and have content
    isfile("pioneer_search_log.txt") && filesize("pioneer_search_log.txt") > 0
end
```

### 3. Performance Benchmark

```julia
using BenchmarkTools

# Compare sync vs async logging
sync_time = @benchmark begin
    logger = create_sync_logger("sync.log")
    for i in 1:1000
        log_message(logger, "Message $i")
    end
end

async_time = @benchmark begin
    logger = create_async_file_logger("async.log")
    for i in 1:1000
        log_message(logger, "Message $i")
    end
    close_async_logger(logger)
end

# Async should be significantly faster
@test median(async_time).time < median(sync_time).time * 0.5
```

## Rollout Plan

### Step 1: Implement AsyncFileLogger
- Create new file `AsyncFileLogger.jl`
- Add unit tests
- Verify basic functionality

### Step 2: Integration
- Update `Loggers.jl` to use async loggers
- Update `LoggingSystem.jl` for proper shutdown
- Test with small dataset

### Step 3: Production Testing
- Test with full dataset
- Monitor memory usage
- Check for message loss

### Step 4: Optimization
- Tune channel sizes
- Adjust flush intervals
- Add monitoring/metrics

## Monitoring and Diagnostics

```julia
function get_logger_stats(logger::AsyncFileLogger)
    return (
        channel_size = logger.channel_size,
        messages_queued = length(logger.message_channel.data),
        is_active = logger.is_active[],
        task_state = logger.writer_task.state,
        file_open = logger.io !== nothing
    )
end

# Add to debug output
@debug_l1 "Logger stats" stats=get_logger_stats(logger)
```

## Configuration Options

Add to `LoggingConfig`:

```julia
struct LoggingConfig
    # ... existing fields ...
    
    # Async logger settings
    async_channel_size::Int           # Default: 10000
    async_flush_interval::Float64     # Default: 1.0 seconds
    async_flush_count::Int            # Default: 100 messages
    async_overflow_policy::Symbol     # :drop, :block, :stderr
end
```

## Conclusion

This asynchronous approach provides a robust, permanent solution that:
1. Eliminates I/O blocking issues
2. Prevents the library loading hang
3. Scales well for future needs
4. Maintains all logging functionality
5. Provides better performance

Unlike the delayed initialization approach, this solves the root cause and prevents any future I/O conflicts throughout the entire pipeline.