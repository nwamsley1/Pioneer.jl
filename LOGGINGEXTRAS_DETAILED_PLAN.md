# Pioneer.jl LoggingExtras-Based Logging System

## Overview

This plan leverages LoggingExtras.jl to create a sophisticated, composable logging system for Pioneer.jl with multiple output streams, custom macros, and flexible console verbosity controls.

## Architecture Design

### Core Components

```julia
using LoggingExtras, Logging, ProgressBars
using Dates, Printf

# Custom log levels for Pioneer
const PIONEER_USER = LogLevel(-100)     # User-facing messages (replace dual_println)
const PIONEER_DEBUG1 = LogLevel(100)    # Basic debugging  
const PIONEER_DEBUG2 = LogLevel(200)    # Detailed debugging
const PIONEER_TRACE = LogLevel(300)     # Function tracing
const PIONEER_PROGRESS = LogLevel(-200) # Progress indicators
```

### LoggingExtras-Based Architecture

```julia
struct PioneerLoggingSystem
    console_logger::AbstractLogger
    simplified_logger::AbstractLogger  
    full_logger::AbstractLogger
    composed_logger::AbstractLogger
    console_verbosity::Ref{LogLevel}    # Runtime adjustable console verbosity
    show_progress::Ref{Bool}
end
```

## Implementation Options

### Option A: Symmetric Multi-Logger (Recommended)

Uses TeeLogger to distribute messages to three parallel streams with independent filtering:

```julia
function create_pioneer_logging_system(
    log_dir::String;
    console_level::LogLevel = PIONEER_USER,
    simplified_level::LogLevel = PIONEER_USER,
    full_level::LogLevel = PIONEER_TRACE,
    console_verbosity_override::Union{Nothing, LogLevel} = nothing
)
    # Console logger with runtime verbosity control
    console_base = ConsoleLogger(stdout, console_verbosity_override === nothing ? console_level : console_verbosity_override)
    console_logger = TransformerLogger(console_base) do log_args
        # Add console formatting and progress bar coordination
        format_for_console(log_args)
    end
    
    # Simplified log (equivalent to dual_println output)
    simplified_path = joinpath(log_dir, "pioneer_simple.log")
    simplified_base = FileLogger(simplified_path)
    simplified_logger = MinLevelLogger(
        TransformerLogger(simplified_base) do log_args
            format_simplified_log(log_args)
        end,
        simplified_level
    )
    
    # Full debug log (everything including @info, debug levels)
    full_path = joinpath(log_dir, "pioneer_full.log") 
    full_base = DatetimeRotatingFileLogger(full_path, dateformat"yyyy-mm-dd_HH")
    full_logger = MinLevelLogger(
        TransformerLogger(full_base) do log_args
            format_full_log(log_args)
        end,
        full_level
    )
    
    # Compose all loggers with TeeLogger
    composed_logger = TeeLogger(
        console_logger,
        simplified_logger, 
        full_logger
    )
    
    return PioneerLoggingSystem(
        console_logger,
        simplified_logger,
        full_logger,
        composed_logger,
        Ref(console_level),
        Ref(true)
    )
end
```

### Option B: Hierarchical Filter Chain

Uses nested filters to create a waterfall of increasingly detailed logging:

```julia
function create_hierarchical_logging(log_dir::String)
    # Start with most verbose logger (full debug)
    full_logger = DatetimeRotatingFileLogger(
        joinpath(log_dir, "pioneer_full.log"),
        dateformat"yyyy-mm-dd_HH"
    )
    
    # Filter for simplified log (only user messages)
    simplified_logger = TeeLogger(
        full_logger,
        EarlyFilteredLogger(
            FileLogger(joinpath(log_dir, "pioneer_simple.log"))
        ) do log_args
            log_args.level <= PIONEER_USER
        end
    )
    
    # Filter for console (user messages + warnings, with verbosity override)
    console_logger = TeeLogger(
        simplified_logger,
        ActiveFilteredLogger(
            ConsoleLogger(stdout)
        ) do log_args
            should_show_on_console(log_args)
        end
    )
    
    return console_logger
end
```

### Option C: Dynamic Router Logger

Custom logger that routes messages based on runtime configuration:

```julia
struct DynamicRoutingLogger <: AbstractLogger
    loggers::Dict{Symbol, AbstractLogger}
    routing_rules::Dict{LogLevel, Vector{Symbol}}
    console_verbosity::Ref{LogLevel}
    config::Ref{Dict{String, Any}}
end

function handle_message(logger::DynamicRoutingLogger, level, message, _module, group, id, file, line; kwargs...)
    # Route message based on level and current configuration
    destinations = get_destinations(logger, level)
    
    for dest in destinations
        if dest == :console && level > logger.console_verbosity[]
            continue  # Skip console if above verbosity threshold
        end
        
        target_logger = logger.loggers[dest]
        Logging.handle_message(target_logger, level, message, _module, group, id, file, line; kwargs...)
    end
end
```

## Custom Macros Design

### Core Pioneer Macros

```julia
# Replace dual_println with structured logging
macro user_info(msg, args...)
    quote
        @logmsg PIONEER_USER $(esc(msg)) $(map(esc, args)...)
    end
end

macro user_info_fmt(fmt, args...)
    quote
        @logmsg PIONEER_USER Printf.@sprintf($(esc(fmt)), $(map(esc, args)...))
    end
end

# Debug levels for replacing @info
macro debug_l1(msg, args...)
    quote
        @logmsg PIONEER_DEBUG1 $(esc(msg)) $(map(esc, args)...)
    end
end

macro debug_l2(msg, args...)
    quote
        @logmsg PIONEER_DEBUG2 $(esc(msg)) $(map(esc, args)...)
    end
end

# Function tracing
macro trace(func_name, args...)
    quote
        @logmsg PIONEER_TRACE "→ $($(string(func_name)))" $(map(esc, args)...)
    end
end

macro trace_exit(func_name, result=nothing)
    quote
        if $(esc(result)) !== nothing
            @logmsg PIONEER_TRACE "← $($(string(func_name))): $($(esc(result)))"
        else
            @logmsg PIONEER_TRACE "← $($(string(func_name)))"
        end
    end
end
```

### Progress Integration Macros

```julia
# Progress-aware logging that coordinates with ProgressBars
macro progress_info(msg, args...)
    quote
        if current_progress_active()
            # Log without interfering with progress bar
            @logmsg PIONEER_PROGRESS $(esc(msg)) $(map(esc, args)...)
        else
            @user_info $(esc(msg)) $(map(esc, args)...)
        end
    end
end

# Smart iteration that chooses between progress bar and logging
macro smart_iterate(iterable, description)
    quote
        items = $(esc(iterable))
        desc = $(esc(description))
        
        if should_show_progress(length(items)) && get_show_progress()
            ProgressBar(items; description=desc)
        else
            @user_info "Processing $desc ($(length(items)) items)"
            items
        end
    end
end
```

## Console Verbosity Control

### Runtime Verbosity Adjustment

```julia
# Global verbosity control
const CONSOLE_VERBOSITY = Ref{LogLevel}(PIONEER_USER)
const SHOW_PROGRESS = Ref{Bool}(true)

function set_console_verbosity!(level::LogLevel)
    CONSOLE_VERBOSITY[] = level
    update_console_logger_level!(level)
    @user_info "Console verbosity set to $(level)"
end

function increase_console_verbosity!()
    current = CONSOLE_VERBOSITY[]
    new_level = if current == PIONEER_USER
        PIONEER_DEBUG1
    elseif current == PIONEER_DEBUG1  
        PIONEER_DEBUG2
    elseif current == PIONEER_DEBUG2
        PIONEER_TRACE
    else
        current
    end
    set_console_verbosity!(new_level)
end

function decrease_console_verbosity!()
    current = CONSOLE_VERBOSITY[]
    new_level = if current == PIONEER_TRACE
        PIONEER_DEBUG2
    elseif current == PIONEER_DEBUG2
        PIONEER_DEBUG1
    elseif current == PIONEER_DEBUG1
        PIONEER_USER
    else
        current
    end
    set_console_verbosity!(new_level)
end
```

### Command Line Interface

```julia
# Command line arguments for SearchDIA/BuildSpecLib
struct LoggingArgs
    console_verbosity::String      # "user", "debug1", "debug2", "trace"
    log_dir::String               # Directory for log files
    show_progress::Bool           # Show progress bars
    simplified_log::Bool          # Enable simplified log file
    full_log::Bool               # Enable full debug log file
    console_override::Bool        # Override console level temporarily
end

function parse_logging_args(parsed_args)
    LoggingArgs(
        get(parsed_args, :console_verbosity, "user"),
        get(parsed_args, :log_dir, "logs"),
        get(parsed_args, :show_progress, true),
        get(parsed_args, :simplified_log, true),
        get(parsed_args, :full_log, true),
        get(parsed_args, :console_override, false)
    )
end

# Usage:
# julia SearchDIA.jl params.json --console-verbosity=debug1 --show-progress=false
# julia BuildSpecLib.jl params.json --console-verbosity=trace --log-dir=debug_logs
```

### Interactive Verbosity Control

```julia
# Interactive commands during execution
function setup_interactive_logging()
    @async begin
        while true
            input = readline()
            if input == "v+"
                increase_console_verbosity!()
            elseif input == "v-"
                decrease_console_verbosity!()
            elseif input == "p"
                toggle_progress_display!()
            elseif startswith(input, "vset:")
                level_str = input[6:end]
                try
                    level = parse_verbosity_level(level_str)
                    set_console_verbosity!(level)
                catch
                    @warn "Invalid verbosity level: $level_str"
                end
            end
        end
    end
end
```

## Formatters and Transformers

### Console Formatter

```julia
function format_for_console(log_args)
    timestamp = Dates.format(now(), "HH:MM:SS")
    level_str = format_level_for_console(log_args.level)
    
    # Different formatting based on level
    formatted_message = if log_args.level == PIONEER_USER
        "$(log_args.message)"  # Clean user messages
    elseif log_args.level == PIONEER_PROGRESS
        "[PROGRESS] $(log_args.message)"
    elseif log_args.level >= PIONEER_DEBUG1
        "[$timestamp] $level_str $(log_args.message)"  # Detailed debug format
    else
        "[$timestamp] $level_str $(log_args.message)"
    end
    
    # Return modified log_args
    merge(log_args, (message = formatted_message,))
end
```

### Simplified Log Formatter

```julia
function format_simplified_log(log_args)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Only format user-level messages (dual_println equivalent)
    if log_args.level <= PIONEER_USER
        formatted = "[$timestamp] $(log_args.message)"
        return merge(log_args, (message = formatted,))
    else
        # Filter out non-user messages
        return nothing  # This will skip the message
    end
end
```

### Full Debug Log Formatter

```julia
function format_full_log(log_args)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
    level_str = string(log_args.level)
    module_str = string(log_args._module)
    location = "$(log_args.file):$(log_args.line)"
    
    # Rich formatting with all context
    formatted = "[$timestamp] [$level_str] [$module_str] $(log_args.message)"
    
    # Add kwargs if present
    if !isempty(log_args.kwargs)
        kwargs_str = join(["$k=$v" for (k,v) in log_args.kwargs], ", ")
        formatted *= " | $kwargs_str"
    end
    
    # Add location for debug levels
    if log_args.level >= PIONEER_DEBUG1
        formatted *= " @ $location"
    end
    
    return merge(log_args, (message = formatted,))
end
```

## Progress Bar Integration

### Progress-Aware Console Logger

```julia
struct ProgressAwareConsoleLogger <: AbstractLogger
    base_logger::ConsoleLogger
    active_progress::Ref{Bool}
    min_level::Ref{LogLevel}
end

function handle_message(logger::ProgressAwareConsoleLogger, level, message, args...)
    # Coordinate with active progress bars
    if logger.active_progress[] && level < Logging.Warn
        # Defer non-warning messages while progress is active
        defer_message(level, message, args...)
    else
        # Use base logger
        Logging.handle_message(logger.base_logger, level, message, args...)
    end
end

# Progress bar wrapper
function with_progress_logging(f, iterable, description; show_progress=true)
    if show_progress && should_show_progress(length(iterable))
        set_progress_active!(true)
        try
            pbar = ProgressBar(iterable; description=description)
            return f(pbar)
        finally
            set_progress_active!(false)
            flush_deferred_messages()
        end
    else
        @user_info "Starting: $description ($(length(iterable)) items)"
        result = f(iterable)
        @user_info "Completed: $description"
        return result
    end
end
```

## Configuration System

### JSON Configuration

```json
{
    "logging": {
        "console_verbosity": "user",
        "simplified_log": {
            "enabled": true,
            "level": "user",
            "file": "pioneer_simple.log"
        },
        "full_log": {
            "enabled": true, 
            "level": "trace",
            "file": "pioneer_full.log",
            "rotation": "daily"
        },
        "progress": {
            "show_bars": true,
            "min_items_for_bar": 10
        },
        "formatting": {
            "timestamps": true,
            "location_info": true,
            "color": true
        }
    }
}
```

### Environment Variable Overrides

```julia
function load_logging_config()
    # Start with defaults
    config = default_logging_config()
    
    # Override with environment variables
    if haskey(ENV, "PIONEER_CONSOLE_VERBOSITY")
        config.console_verbosity = parse_verbosity_level(ENV["PIONEER_CONSOLE_VERBOSITY"])
    end
    
    if haskey(ENV, "PIONEER_LOG_DIR")
        config.log_dir = ENV["PIONEER_LOG_DIR"]
    end
    
    if haskey(ENV, "PIONEER_SHOW_PROGRESS")
        config.show_progress = parse(Bool, ENV["PIONEER_SHOW_PROGRESS"])
    end
    
    return config
end

# Usage:
# PIONEER_CONSOLE_VERBOSITY=debug1 julia SearchDIA.jl params.json
# PIONEER_LOG_DIR=/tmp/pioneer_logs julia BuildSpecLib.jl params.json
```

## Migration Strategy

### Phase 1: Infrastructure Setup (Week 1-2)

```julia
# Create the logging module
module PioneerLogging
    using LoggingExtras, Logging, Dates, Printf
    
    # Export custom macros and functions
    export @user_info, @debug_l1, @debug_l2, @trace
    export setup_pioneer_logging!, set_console_verbosity!
    export with_progress_logging, @smart_iterate
    
    # Implementation here...
end
```

### Phase 2: Backward Compatibility (Week 3)

```julia
# Provide transition functions
function dual_println(args...)
    @user_info join(string.(args), " ")
end

function dual_print(args...)
    @user_info join(string.(args), "")
end

# Mark as deprecated
Base.@deprecate dual_println(args...) @user_info(join(string.(args), " "))
Base.@deprecate dual_print(args...) @user_info(join(string.(args), ""))
```

### Phase 3: Gradual Migration (Week 4-6)

1. Replace `dual_println` calls with `@user_info`
2. Replace existing `@info` calls with appropriate `@debug_l1`/`@debug_l2`
3. Add `@trace` calls to key functions
4. Update progress bar usage with `@smart_iterate`

### Phase 4: Advanced Features (Week 7-8)

1. Implement interactive verbosity control
2. Add log file rotation and management
3. Performance optimization
4. Comprehensive testing

## Usage Examples

### Basic Setup

```julia
# In SearchDIA.jl or BuildSpecLib.jl
using PioneerLogging

function main_SearchDIA(argv=ARGS)
    # Parse command line arguments including logging options
    parsed_args = parse_args(argv, settings; as_symbols = true)
    logging_args = parse_logging_args(parsed_args)
    
    # Setup logging system
    setup_pioneer_logging!(
        logging_args.log_dir;
        console_verbosity = parse_verbosity_level(logging_args.console_verbosity),
        show_progress = logging_args.show_progress
    )
    
    @user_info "Starting SearchDIA analysis"
    @debug_l1 "Command line arguments" args=argv
    
    try
        SearchDIA(parsed_args[:params_path])
        @user_info "SearchDIA completed successfully"
    catch e
        @error "SearchDIA failed" exception=e
        return 1
    end
    return 0
end
```

### In Search Methods

```julia
function process_files!(results, params, context, files)
    @trace process_files!
    @debug_l1 "Processing $(length(files)) files"
    
    for (i, file) in @smart_iterate(files, "Processing MS files")
        @debug_l2 "Processing file" file=file index=i
        
        try
            process_single_file!(results, params, context, file)
            @debug_l1 "File processed successfully" file=file
        catch e
            @warn "Failed to process file" file=file error=e
            continue
        end
    end
    
    @trace_exit process_files! length(results)
end
```

### Console Verbosity Demo

```bash
# Normal operation (user messages only)
$ julia SearchDIA.jl params.json
14:30:00 Starting SearchDIA analysis
14:30:01 Processing MS files (5 files)
[████████████████████████████████████████] 100%
14:32:15 SearchDIA completed successfully

# Debug level 1
$ julia SearchDIA.jl params.json --console-verbosity=debug1  
[14:30:00] USER Starting SearchDIA analysis
[14:30:00] DEBUG1 Command line arguments: ["params.json", "--console-verbosity=debug1"]
[14:30:01] DEBUG1 Processing 5 files
[14:30:01] DEBUG1 File processed successfully: file1.mzML
[14:30:30] DEBUG1 File processed successfully: file2.mzML
...
[14:32:15] USER SearchDIA completed successfully

# Full trace mode
$ julia SearchDIA.jl params.json --console-verbosity=trace
[14:30:00] USER Starting SearchDIA analysis
[14:30:00] DEBUG1 Command line arguments: ["params.json", "--console-verbosity=trace"]
[14:30:00] TRACE → process_files!
[14:30:01] DEBUG1 Processing 5 files  
[14:30:01] DEBUG2 Processing file: file1.mzML, index: 1
[14:30:01] TRACE → process_single_file!
[14:30:01] DEBUG2 Fragment match details: {...}
[14:30:01] TRACE ← process_single_file!
...
[14:32:15] TRACE ← process_files!: 1250
[14:32:15] USER SearchDIA completed successfully
```

## Performance Considerations

### Lazy Message Evaluation

```julia
# Use Julia's lazy evaluation for expensive debug messages
@debug_l2 "Expensive computation result" result=expensive_computation() _group=:performance

# Or conditional evaluation
if should_log(PIONEER_DEBUG2)
    result = expensive_computation()
    @debug_l2 "Expensive computation result" result=result
end
```

### Asynchronous Logging Option

```julia
struct AsyncPioneerLogger <: AbstractLogger
    sync_logger::AbstractLogger
    async_buffer::Channel{LogMessage}
    worker_task::Task
end

function create_async_logger(base_logger::AbstractLogger; buffer_size=1000)
    buffer = Channel{LogMessage}(buffer_size)
    
    worker = @async begin
        while true
            msg = take!(buffer)
            handle_message(base_logger, msg.level, msg.message, msg.args...)
        end
    end
    
    return AsyncPioneerLogger(base_logger, buffer, worker)
end
```

## Testing Strategy

### Unit Tests

```julia
@testset "PioneerLogging Tests" begin
    @testset "Custom Macros" begin
        # Test macro expansion and level assignment
        io = IOBuffer()
        logger = ConsoleLogger(io, PIONEER_DEBUG2)
        
        with_logger(logger) do
            @user_info "test user message"
            @debug_l1 "test debug message"
            @debug_l2 "test detailed debug"
        end
        
        output = String(take!(io))
        @test contains(output, "test user message")
        @test contains(output, "test debug message") 
        @test contains(output, "test detailed debug")
    end
    
    @testset "Console Verbosity Control" begin
        # Test runtime verbosity changes
    end
    
    @testset "Log File Output" begin
        # Test simplified vs full log content
    end
end
```

### Integration Tests

```julia
@testset "SearchDIA Logging Integration" begin
    mktempdir() do log_dir
        setup_pioneer_logging!(log_dir; console_verbosity=PIONEER_DEBUG1)
        
        # Run a small SearchDIA test
        SearchDIA("data/ecoli_test/ecoli_test_params.json")
        
        # Verify log files were created and contain expected content
        simplified_log = joinpath(log_dir, "pioneer_simple.log")
        full_log = joinpath(log_dir, "pioneer_full.log")
        
        @test isfile(simplified_log)
        @test isfile(full_log)
        
        # Check content
        simplified_content = read(simplified_log, String)
        full_content = read(full_log, String)
        
        @test contains(simplified_content, "Starting SearchDIA")
        @test length(split(full_content, "\n")) > length(split(simplified_content, "\n"))
    end
end
```

## Comparison of Options

### Option A: Symmetric Multi-Logger
**Pros:**
- Clear separation of concerns
- Easy to understand and debug
- Flexible configuration per destination
- Good performance with parallel writes

**Cons:**  
- Slightly more complex setup
- Three separate logger configurations to manage

### Option B: Hierarchical Filter Chain
**Pros:**
- Single logger configuration
- Messages flow through natural hierarchy
- Simpler conceptually

**Cons:**
- Less flexible - changing one filter affects downstream
- Harder to customize formatting per destination

### Option C: Dynamic Router Logger
**Pros:**
- Most flexible - can change routing at runtime
- Single configuration point
- Supports complex routing rules

**Cons:**
- More complex implementation
- Potential performance overhead
- Harder to debug routing issues

## Recommendation

**Option A (Symmetric Multi-Logger)** is recommended because:

1. **LoggingExtras Strengths**: Leverages TeeLogger and TransformerLogger effectively
2. **Clear Architecture**: Each destination has independent configuration
3. **Runtime Flexibility**: Console verbosity can be adjusted without affecting file logs
4. **Performance**: Parallel writes to different destinations
5. **Maintainability**: Easy to understand and modify individual components
6. **Testing**: Each logger component can be tested independently

The implementation provides all requested features:
- Multiple output streams (console, simplified log, full debug log)
- Custom macros for different message types
- Runtime console verbosity control
- Progress bar integration
- Backward compatibility with existing dual_println

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"id": "1", "content": "Research LoggingExtras.jl capabilities and API", "status": "completed"}, {"id": "2", "content": "Design LoggingExtras-based architecture with custom macros", "status": "completed"}, {"id": "3", "content": "Create detailed implementation plan with console verbosity options", "status": "completed"}]