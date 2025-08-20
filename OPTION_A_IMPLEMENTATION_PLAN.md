# Pioneer.jl Option A Logging Implementation Plan
## Symmetric Multi-Logger with LoggingExtras.jl

### Table of Contents
1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Module Structure](#module-structure)
4. [Core Implementation](#core-implementation)
5. [Custom Macros](#custom-macros)
6. [Configuration System](#configuration-system)
7. [Progress Bar Integration](#progress-bar-integration)
8. [Warnings Tracking System](#warnings-tracking-system)
9. [File Management](#file-management)
10. [Migration Strategy](#migration-strategy)
11. [Testing Plan](#testing-plan)
12. [Performance Considerations](#performance-considerations)
13. [Usage Examples](#usage-examples)

---

## Overview

This plan implements a comprehensive logging system using LoggingExtras.jl's TeeLogger to create three independent, parallel logging streams with integrated warnings tracking:

1. **Console Logger**: User messages, progress bars, warnings (with runtime verbosity control)
2. **Simplified Logger**: Clean log file equivalent to current dual_println output
3. **Full Debug Logger**: Complete debugging output with all levels and metadata
4. **Warnings Tracker**: Automatic warning collection and end-of-run reporting

### Key Design Principles
- **Independence**: Each logger operates independently with its own configuration
- **Flexibility**: Console verbosity can be changed between SearchDIA runs
- **Performance**: Parallel processing, no cascading filters
- **Maintainability**: Clear separation of concerns using LoggingExtras patterns
- **Backward Compatibility**: Smooth migration from existing dual_println system
- **Comprehensive Warning Management**: Automatic collection, categorization, and reporting of warnings

---

## Architecture

### High-Level Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   Log Message   │───▶│WarningCapturing │───▶│  Three Parallel │
│   (@user_info,  │    │    Logger       │    │   Destinations  │
│   @debug_l1,    │    │  (TeeLogger +   │    │   + Warning     │
│   @debug_l2,    │    │   Warning       │    │   Tracking)     │
│   @trace, @warn,│    │   Tracker)      │    │                 │
│   etc.)         │    │                  │    │                 │
└─────────────────┘    └──────────────────┘    └─────────────────┘
                                │
                ┌───────────────┼───────────────┐
                │               │               │
                ▼               ▼               ▼
        ┌───────────────┐ ┌─────────────┐ ┌─────────────┐
        │ Console       │ │ Simplified  │ │ Full Debug  │
        │ Logger        │ │ Logger      │ │ Logger      │
        │               │ │             │ │             │
        │ - User msgs   │ │ - User msgs │ │ - All levels│
        │ - Warnings    │ │ - Timestmp  │ │ - Metadata  │
        │ - Progress    │ │ - Clean fmt │ │ - Location  │
        │ - Runtime     │ │             │ │ - Rotation  │
        │   verbosity   │ │             │ │             │
        └───────────────┘ └─────────────┘ └─────────────┘
                │               │               │
                ▼               ▼               ▼
        ┌───────────────┐ ┌─────────────┐ ┌─────────────┐
        │    stdout     │ │ simple.log  │ │ full.log    │
        └───────────────┘ └─────────────┘ └─────────────┘
                                               
                        ┌─────────────────────┐
                        │  Warning System     │
                        │  ─────────────────  │
                        │ • Auto Collection   │
                        │ • Categorization    │
                        │ • End-of-run Report │
                        │ • Console Summary   │
                        └─────────────────────┘
                                 │
                                 ▼
                        ┌─────────────────────┐
                        │   warnings.txt      │
                        └─────────────────────┘
```

### Component Breakdown

```julia
# Core logging system structure
PioneerLoggingSystem
├── WarningCapturingLogger (wraps TeeLogger + warning tracking)
│   ├── WarningTracker (automatic collection and categorization)
│   └── TeeLogger (distributes to all destinations)
│       ├── ConsoleLogger (with TransformerLogger for formatting)
│       │   ├── MinLevelLogger (runtime adjustable level)
│       │   ├── ProgressAwareFormatter 
│       │   └── ColorFormatter
│       ├── SimplifiedLogger (FileLogger with TransformerLogger)
│       │   ├── MinLevelLogger (PIONEER_USER level)
│       │   ├── SimpleFormatter (timestamp + message)
│       │   └── SafeFileWriter
│       └── FullDebugLogger (DatetimeRotatingFileLogger with TransformerLogger)
│           ├── MinLevelLogger (PIONEER_TRACE level)
│           ├── DebugFormatter (full metadata)
│           └── RotationManager
└── WarningReportGenerator (end-of-run summary and file generation)
```

---

## Module Structure

### File Organization

```
src/utils/PioneerLogging/
├── PioneerLogging.jl          # Main module file
├── core/
│   ├── LogLevels.jl          # Custom log level definitions
│   ├── LoggingSystem.jl      # Main system struct and setup
│   └── Configuration.jl      # Configuration management
├── loggers/
│   ├── ConsoleLogger.jl      # Console-specific logger
│   ├── SimplifiedLogger.jl   # Simplified file logger
│   └── FullDebugLogger.jl    # Full debug file logger
├── formatters/
│   ├── ConsoleFormatter.jl   # Console message formatting
│   ├── SimpleFormatter.jl    # Simplified log formatting
│   └── DebugFormatter.jl     # Full debug formatting
├── macros/
│   ├── UserMacros.jl         # @user_info, @warn_with_context, etc.
│   ├── DebugMacros.jl        # @debug_l1, @debug_l2
│   └── TraceMacros.jl        # @trace, @trace_exit
├── progress/
│   ├── ProgressIntegration.jl # Progress bar coordination
│   └── SmartIteration.jl     # @smart_iterate implementation
├── warnings/
│   ├── WarningTracker.jl     # Core warning tracking infrastructure
│   ├── WarningCapture.jl     # Warning capturing logger
│   ├── WarningReports.jl     # Report generation and console summaries
│   └── WarningCategories.jl  # Warning categorization system
├── utils/
│   ├── FileManagement.jl     # Log file rotation, cleanup
│   ├── LevelUtils.jl         # Level parsing and conversion
│   └── SessionManagement.jl  # Multiple session support
└── testing/
    ├── TestUtils.jl          # Testing utilities
    ├── MockLoggers.jl        # Mock loggers for testing
    └── WarningTestUtils.jl   # Warning system testing utilities
```

### Module Dependencies

```julia
# Main module imports
using LoggingExtras, Logging
using Dates, Printf
using ProgressBars
using UUIDs  # For session ID generation
```

---

## Core Implementation

### 1. Custom Log Levels

```julia
# src/utils/PioneerLogging/core/LogLevels.jl

"""
Custom log levels for Pioneer.jl

Log level hierarchy (lower values have higher priority):
- PIONEER_PROGRESS (-200): Progress indicators  
- PIONEER_USER (-100): User-facing messages (replaces dual_println)
- Info (0): Standard Julia info level
- Warn (1000): Standard Julia warning level
- PIONEER_DEBUG1 (1100): Basic debugging information
- PIONEER_DEBUG2 (1200): Detailed debugging information  
- PIONEER_TRACE (1300): Function tracing and detailed execution flow
"""

const PIONEER_PROGRESS = LogLevel(-200)
const PIONEER_USER = LogLevel(-100)
const PIONEER_DEBUG1 = LogLevel(1100)
const PIONEER_DEBUG2 = LogLevel(1200)
const PIONEER_TRACE = LogLevel(1300)

# Level name mapping for display
const LEVEL_NAMES = Dict(
    PIONEER_PROGRESS => "PROGRESS",
    PIONEER_USER => "USER",
    Logging.Info => "INFO",
    Logging.Warn => "WARN",
    Logging.Error => "ERROR",
    PIONEER_DEBUG1 => "DEBUG1",
    PIONEER_DEBUG2 => "DEBUG2",
    PIONEER_TRACE => "TRACE"
)

function level_name(level::LogLevel)::String
    return get(LEVEL_NAMES, level, string(level))
end

function parse_level(level_str::String)::LogLevel
    level_map = Dict(
        "progress" => PIONEER_PROGRESS,
        "user" => PIONEER_USER,
        "info" => Logging.Info,
        "warn" => Logging.Warn,
        "error" => Logging.Error,
        "debug1" => PIONEER_DEBUG1,
        "debug2" => PIONEER_DEBUG2,
        "trace" => PIONEER_TRACE
    )
    
    normalized = lowercase(strip(level_str))
    haskey(level_map, normalized) || error("Unknown log level: $level_str")
    return level_map[normalized]
end
```

### 2. Configuration System

```julia
# src/utils/PioneerLogging/core/Configuration.jl

"""
Configuration for Pioneer logging system
"""
Base.@kwdef struct PioneerLoggingConfig
    # Log directory
    log_dir::String = "logs"
    
    # Console settings
    console_level::LogLevel = PIONEER_USER
    console_color::Bool = true
    show_progress::Bool = true
    
    # Simplified log settings  
    simplified_enabled::Bool = true
    simplified_level::LogLevel = PIONEER_USER
    simplified_filename::String = "pioneer_simple.log"
    
    # Full debug log settings
    full_enabled::Bool = true
    full_level::LogLevel = PIONEER_TRACE
    full_filename::String = "pioneer_full.log"
    full_rotation::Bool = true
    full_rotation_size::Int = 50_000_000  # 50MB
    full_max_files::Int = 5
    
    # Progress settings
    progress_min_items::Int = 10
    progress_update_interval::Float64 = 0.1
    
    # Warning tracking settings
    track_warnings::Bool = true
    max_warnings::Int = 10000
    warnings_filename::String = "warnings.txt"
    console_summary_threshold::Int = 5  # Show summary if more than 5 warnings
    critical_warning_threshold::Int = 1  # Alert if any critical warnings
    ignored_warning_patterns::Vector{String} = String[]  # Warning patterns to ignore
    
    # Performance settings
    async_logging::Bool = false
    buffer_size::Int = 1000
end

function load_config_from_json(path::String)::PioneerLoggingConfig
    if !isfile(path)
        @warn "Config file not found: $path, using defaults"
        return PioneerLoggingConfig()
    end
    
    try
        json_data = JSON.parsefile(path)
        logging_config = get(json_data, "logging", Dict())
        
        return PioneerLoggingConfig(
            log_dir = get(logging_config, "log_dir", "logs"),
            console_level = parse_level(get(logging_config, "console_level", "user")),
            console_color = get(logging_config, "console_color", true),
            show_progress = get(logging_config, "show_progress", true),
            simplified_enabled = get(logging_config, "simplified_enabled", true),
            simplified_level = parse_level(get(logging_config, "simplified_level", "user")),
            simplified_filename = get(logging_config, "simplified_filename", "pioneer_simple.log"),
            full_enabled = get(logging_config, "full_enabled", true),
            full_level = parse_level(get(logging_config, "full_level", "trace")),
            full_filename = get(logging_config, "full_filename", "pioneer_full.log"),
            # ... other fields
        )
    catch e
        @warn "Error loading config from $path: $e, using defaults"
        return PioneerLoggingConfig()
    end
end

function load_config_from_env()::PioneerLoggingConfig
    config = PioneerLoggingConfig()
    
    # Override with environment variables
    if haskey(ENV, "PIONEER_LOG_DIR")
        config = @set config.log_dir = ENV["PIONEER_LOG_DIR"]
    end
    
    if haskey(ENV, "PIONEER_CONSOLE_LEVEL")
        config = @set config.console_level = parse_level(ENV["PIONEER_CONSOLE_LEVEL"])
    end
    
    if haskey(ENV, "PIONEER_SHOW_PROGRESS")
        config = @set config.show_progress = parse(Bool, ENV["PIONEER_SHOW_PROGRESS"])
    end
    
    return config
end
```

### 3. Main Logging System

```julia
# src/utils/PioneerLogging/core/LoggingSystem.jl

"""
Main Pioneer logging system using TeeLogger for symmetric multi-destination logging
"""
mutable struct PioneerLoggingSystem
    config::PioneerLoggingConfig
    console_logger::AbstractLogger
    simplified_logger::Union{AbstractLogger, Nothing}
    full_logger::Union{AbstractLogger, Nothing}
    tee_logger::AbstractLogger
    
    # Runtime state
    console_level::Ref{LogLevel}
    progress_active::Ref{Bool}
    session_id::String
    created_at::DateTime
    
    # File handles (for proper cleanup)
    simplified_file::Union{IO, Nothing}
    full_file::Union{IO, Nothing}
end

function PioneerLoggingSystem(config::PioneerLoggingConfig)
    # Ensure log directory exists
    mkpath(config.log_dir)
    
    # Generate unique session ID
    session_id = string(uuid4())[1:8]
    created_at = now()
    
    # Create console logger
    console_logger = create_console_logger(config, session_id)
    
    # Create simplified logger (if enabled)
    simplified_logger, simplified_file = if config.simplified_enabled
        create_simplified_logger(config, session_id)
    else
        nothing, nothing
    end
    
    # Create full debug logger (if enabled)
    full_logger, full_file = if config.full_enabled
        create_full_logger(config, session_id)
    else
        nothing, nothing
    end
    
    # Create TeeLogger to combine all active loggers
    active_loggers = filter(!isnothing, [console_logger, simplified_logger, full_logger])
    tee_logger = TeeLogger(active_loggers...)
    
    return PioneerLoggingSystem(
        config,
        console_logger,
        simplified_logger,
        full_logger,
        tee_logger,
        Ref(config.console_level),
        Ref(false),
        session_id,
        created_at,
        simplified_file,
        full_file
    )
end

# Global system instance
const GLOBAL_LOGGING_SYSTEM = Ref{Union{PioneerLoggingSystem, Nothing}}(nothing)

function setup_pioneer_logging!(
    log_dir::String = "logs";
    console_level::LogLevel = PIONEER_USER,
    simplified_enabled::Bool = true,
    full_enabled::Bool = true,
    show_progress::Bool = true,
    config_file::Union{String, Nothing} = nothing
)
    # Load configuration
    config = if config_file !== nothing
        load_config_from_json(config_file)
    else
        load_config_from_env()
    end
    
    # Override with function arguments
    config = @set config.log_dir = log_dir
    config = @set config.console_level = console_level
    config = @set config.simplified_enabled = simplified_enabled
    config = @set config.full_enabled = full_enabled
    config = @set config.show_progress = show_progress
    
    # Clean up previous system if it exists
    cleanup_logging_system!()
    
    # Create new system
    system = PioneerLoggingSystem(config)
    GLOBAL_LOGGING_SYSTEM[] = system
    
    # Set as global logger
    global_logger(system.tee_logger)
    
    @info "Pioneer logging system initialized" session_id=system.session_id log_dir=config.log_dir
    
    return system
end

function cleanup_logging_system!()
    system = GLOBAL_LOGGING_SYSTEM[]
    if system !== nothing
        # Close file handles
        if system.simplified_file !== nothing
            close(system.simplified_file)
        end
        if system.full_file !== nothing
            close(system.full_file)
        end
        
        @debug "Cleaned up logging system" session_id=system.session_id
        GLOBAL_LOGGING_SYSTEM[] = nothing
    end
end

function get_logging_system()::Union{PioneerLoggingSystem, Nothing}
    return GLOBAL_LOGGING_SYSTEM[]
end
```

### 4. Individual Logger Implementations

```julia
# src/utils/PioneerLogging/loggers/ConsoleLogger.jl

function create_console_logger(config::PioneerLoggingConfig, session_id::String)
    base_console = ConsoleLogger(stdout, config.console_level)
    
    # Add transformer for console formatting
    transformer = TransformerLogger(base_console) do log_args
        format_console_message(log_args, config)
    end
    
    # Add level filtering with runtime adjustability
    level_filter = RuntimeLevelLogger(transformer, config.console_level)
    
    return level_filter
end

"""
Logger that allows runtime level adjustment
"""
mutable struct RuntimeLevelLogger <: AbstractLogger
    base_logger::AbstractLogger
    min_level::Ref{LogLevel}
end

function RuntimeLevelLogger(base_logger::AbstractLogger, initial_level::LogLevel)
    return RuntimeLevelLogger(base_logger, Ref(initial_level))
end

function Logging.min_enabled_level(logger::RuntimeLevelLogger)
    return logger.min_level[]
end

function Logging.shouldlog(logger::RuntimeLevelLogger, level, _module, group, id)
    return level >= logger.min_level[] && 
           Logging.shouldlog(logger.base_logger, level, _module, group, id)
end

function Logging.handle_message(logger::RuntimeLevelLogger, level, message, _module, group, id, file, line; kwargs...)
    return Logging.handle_message(logger.base_logger, level, message, _module, group, id, file, line; kwargs...)
end

function set_console_level!(level::LogLevel)
    system = get_logging_system()
    if system !== nothing
        system.console_level[] = level
        # Update the runtime level logger
        if system.console_logger isa RuntimeLevelLogger
            system.console_logger.min_level[] = level
        end
        @user_info "Console verbosity set to $(level_name(level))"
    else
        @warn "No Pioneer logging system active"
    end
end
```

```julia
# src/utils/PioneerLogging/loggers/SimplifiedLogger.jl

function create_simplified_logger(config::PioneerLoggingConfig, session_id::String)
    # Create unique filename with session ID
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(timestamp)_$(session_id)_$(config.simplified_filename)"
    filepath = joinpath(config.log_dir, filename)
    
    # Open file for writing
    file_handle = open(filepath, "w")
    
    # Create file logger
    base_logger = FileLogger(file_handle)
    
    # Add transformer for simplified formatting
    transformer = TransformerLogger(base_logger) do log_args
        format_simplified_message(log_args, config)
    end
    
    # Add level filter
    level_filter = MinLevelLogger(transformer, config.simplified_level)
    
    return level_filter, file_handle
end
```

```julia
# src/utils/PioneerLogging/loggers/FullDebugLogger.jl

function create_full_logger(config::PioneerLoggingConfig, session_id::String)
    # Create rotating file logger
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "$(timestamp)_$(session_id)_$(config.full_filename)"
    filepath = joinpath(config.log_dir, filename)
    
    if config.full_rotation
        # Use DatetimeRotatingFileLogger for automatic rotation
        base_logger = DatetimeRotatingFileLogger(
            filepath,
            dateformat"yyyy-mm-dd_HH";
            rotation_callback = cleanup_old_logs
        )
        file_handle = nothing  # Managed by DatetimeRotatingFileLogger
    else
        # Use simple file logger
        file_handle = open(filepath, "w")
        base_logger = FileLogger(file_handle)
    end
    
    # Add transformer for full debug formatting
    transformer = TransformerLogger(base_logger) do log_args
        format_debug_message(log_args, config)
    end
    
    # Add level filter
    level_filter = MinLevelLogger(transformer, config.full_level)
    
    return level_filter, file_handle
end

function cleanup_old_logs(rotated_file_path::String)
    # This callback is called when a log file is rotated
    # Implement cleanup of old files based on config.full_max_files
    log_dir = dirname(rotated_file_path)
    log_pattern = basename(rotated_file_path)
    
    # Find all matching log files
    all_logs = filter(readdir(log_dir)) do filename
        startswith(filename, split(log_pattern, '_')[1])  # Match date prefix
    end
    
    # Sort by modification time
    sort!(all_logs, by = f -> stat(joinpath(log_dir, f)).mtime, rev = true)
    
    # Remove old files if we exceed max_files
    system = get_logging_system()
    max_files = system !== nothing ? system.config.full_max_files : 5
    
    if length(all_logs) > max_files
        for old_file in all_logs[max_files+1:end]
            old_path = joinpath(log_dir, old_file)
            try
                rm(old_path)
                @debug "Removed old log file" file=old_path
            catch e
                @warn "Failed to remove old log file" file=old_path error=e
            end
        end
    end
end
```

---

## Custom Macros

### 1. User-Facing Macros

```julia
# src/utils/PioneerLogging/macros/UserMacros.jl

"""
User-facing logging macros that replace dual_println functionality
"""

macro user_info(msg, args...)
    quote
        @logmsg PIONEER_USER $(esc(msg)) $(map(esc, args)...)
    end
end

macro user_info_fmt(fmt, args...)
    quote
        formatted_msg = Printf.@sprintf($(esc(fmt)), $(map(esc, args)...))
        @logmsg PIONEER_USER formatted_msg
    end
end

# Convenience macros for common patterns
macro user_info_with_time(msg, args...)
    quote
        timestamp = Dates.format(now(), "HH:MM:SS")
        @logmsg PIONEER_USER "[$timestamp] $($(esc(msg)))" $(map(esc, args)...)
    end
end

macro user_separator(char='-', length=80)
    quote
        separator = repeat($(esc(char)), $(esc(length)))
        @logmsg PIONEER_USER separator
    end
end

# Backward compatibility - marked as deprecated
macro dual_println(args...)
    quote
        Base.depwarn("@dual_println is deprecated, use @user_info instead", :dual_println)
        message = join(string.($(map(esc, args))), " ")
        @logmsg PIONEER_USER message
    end
end
```

### 2. Debug Macros

```julia
# src/utils/PioneerLogging/macros/DebugMacros.jl

"""
Debug logging macros for different verbosity levels
"""

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

# Conditional debug logging to avoid expensive computations
macro debug_l1_if(condition, msg, args...)
    quote
        if $(esc(condition))
            @logmsg PIONEER_DEBUG1 $(esc(msg)) $(map(esc, args)...)
        end
    end
end

macro debug_l2_if(condition, msg, args...)
    quote
        if $(esc(condition))
            @logmsg PIONEER_DEBUG2 $(esc(msg)) $(map(esc, args)...)
        end
    end
end

# Debug with timing
macro debug_l1_timed(msg, expr)
    quote
        start_time = time()
        result = $(esc(expr))
        elapsed = time() - start_time
        @logmsg PIONEER_DEBUG1 "$($(esc(msg))) ($(round(elapsed, digits=3))s)"
        result
    end
end

# Debug with memory allocation tracking
macro debug_l2_alloc(msg, expr)
    quote
        gc_before = Base.gc_num()
        result = $(esc(expr))
        gc_after = Base.gc_num()
        alloc_diff = gc_after.allocd - gc_before.allocd
        @logmsg PIONEER_DEBUG2 "$($(esc(msg))) (allocated $(Base.format_bytes(alloc_diff)))"
        result
    end
end
```

### 3. Trace Macros

```julia
# src/utils/PioneerLogging/macros/TraceMacros.jl

"""
Function tracing macros for detailed execution flow
"""

macro trace(func_name, args...)
    quote
        func_str = $(string(func_name))
        if !isempty($(map(esc, args)))
            args_str = join([string(k) * "=" * string(v) for (k,v) in $(map(esc, args))], ", ")
            @logmsg PIONEER_TRACE "→ $func_str($args_str)"
        else
            @logmsg PIONEER_TRACE "→ $func_str"
        end
    end
end

macro trace_exit(func_name, result=nothing)
    quote
        func_str = $(string(func_name))
        if $(esc(result)) !== nothing
            @logmsg PIONEER_TRACE "← $func_str: $($(esc(result)))"
        else
            @logmsg PIONEER_TRACE "← $func_str"
        end
    end
end

# Trace a function call with automatic entry/exit
macro trace_function(func_call)
    quote
        func_name = $(string(func_call.args[1]))
        @trace $func_name
        try
            result = $(esc(func_call))
            @trace_exit $func_name result
            result
        catch e
            @logmsg PIONEER_TRACE "← $func_name (exception: $e)"
            rethrow(e)
        end
    end
end

# Trace with stack depth indication
let trace_depth = Ref(0)
    global function get_trace_depth()
        return trace_depth[]
    end
    
    global function increment_trace_depth!()
        trace_depth[] += 1
    end
    
    global function decrement_trace_depth!()
        trace_depth[] = max(0, trace_depth[] - 1)
    end
end

macro trace_depth(func_name, args...)
    quote
        increment_trace_depth!()
        depth_str = "  " ^ get_trace_depth()
        func_str = $(string(func_name))
        
        if !isempty($(map(esc, args)))
            args_str = join([string(k) * "=" * string(v) for (k,v) in $(map(esc, args))], ", ")
            @logmsg PIONEER_TRACE "$(depth_str)→ $func_str($args_str)"
        else
            @logmsg PIONEER_TRACE "$(depth_str)→ $func_str"
        end
        
        try
            result = yield()  # This would need to be part of a larger macro pattern
            depth_str = "  " ^ get_trace_depth()
            @logmsg PIONEER_TRACE "$(depth_str)← $func_str"
            result
        finally
            decrement_trace_depth!()
        end
    end
end
```

---

## Progress Bar Integration

### 1. Progress-Aware Logging

```julia
# src/utils/PioneerLogging/progress/ProgressIntegration.jl

"""
Coordinate logging with ProgressBars to avoid output conflicts
"""

function is_progress_active()::Bool
    system = get_logging_system()
    return system !== nothing && system.progress_active[]
end

function set_progress_active!(active::Bool)
    system = get_logging_system()
    if system !== nothing
        system.progress_active[] = active
    end
end

function should_show_progress(item_count::Int)::Bool
    system = get_logging_system()
    if system === nothing
        return item_count >= 10  # Default threshold
    end
    
    return system.config.show_progress && 
           item_count >= system.config.progress_min_items
end

# Deferred message queue for messages that arrive during progress display
const DEFERRED_MESSAGES = Ref{Vector{Any}}(Any[])

function defer_message(level, message, module_, group, id, file, line; kwargs...)
    if is_progress_active() && level < Logging.Warn
        # Store non-warning messages for later display
        push!(DEFERRED_MESSAGES[], (level, message, module_, group, id, file, line, kwargs))
        return
    end
    
    # Let warnings and errors through immediately
    @logmsg level message _module=module_ _group=group _id=id _file=file _line=line kwargs...
end

function flush_deferred_messages!()
    for (level, message, module_, group, id, file, line, kwargs) in DEFERRED_MESSAGES[]
        @logmsg level message _module=module_ _group=group _id=id _file=file _line=line kwargs...
    end
    empty!(DEFERRED_MESSAGES[])
end

function with_progress_logging(f, iterable, description::String; show_progress::Bool=true)
    items = collect(iterable)  # Ensure we can get length
    
    if show_progress && should_show_progress(length(items))
        @user_info "Starting: $description ($(length(items)) items)"
        
        set_progress_active!(true)
        try
            pbar = ProgressBar(items; description=description)
            result = f(pbar)
            set_progress_active!(false)
            flush_deferred_messages!()
            
            @user_info "Completed: $description"
            return result
        catch e
            set_progress_active!(false)
            flush_deferred_messages!()
            @error "Failed during: $description" exception=e
            rethrow(e)
        end
    else
        @user_info "Processing: $description ($(length(items)) items)"
        result = f(items)
        @user_info "Completed: $description"
        return result
    end
end
```

### 2. Smart Iteration Macro

```julia
# src/utils/PioneerLogging/progress/SmartIteration.jl

"""
Smart iteration that automatically chooses between progress bars and simple logging
"""

macro smart_iterate(iterable, description, args...)
    # Parse optional arguments
    show_progress = true
    for arg in args
        if arg isa Expr && arg.head == :(=) && arg.args[1] == :show_progress
            show_progress = arg.args[2]
        end
    end
    
    quote
        items = $(esc(iterable))
        desc = $(esc(description))
        show_prog = $(esc(show_progress))
        
        if show_prog && should_show_progress(length(items))
            set_progress_active!(true)
            pbar = ProgressBar(items; description=desc)
            # Return the progress bar iterator
            pbar
        else
            @user_info "Processing $desc ($(length(items)) items)"
            # Return the original iterable
            items
        end
    end
end

# Usage pattern:
# for item in @smart_iterate(data, "Processing files")
#     # work with item
# end
# if progress was shown, need to clean up:
# if is_progress_active()
#     set_progress_active!(false)
#     flush_deferred_messages!()
# end

# Better version with automatic cleanup:
macro smart_iterate_block(iterable, description, block)
    quote
        items = $(esc(iterable))
        desc = $(esc(description))
        
        if should_show_progress(length(items))
            @user_info "Starting: $desc ($(length(items)) items)"
            set_progress_active!(true)
            try
                for item in ProgressBar(items; description=desc)
                    $(esc(block))(item)
                end
                set_progress_active!(false)
                flush_deferred_messages!()
                @user_info "Completed: $desc"
            catch e
                set_progress_active!(false)
                flush_deferred_messages!()
                @error "Failed during: $desc" exception=e
                rethrow(e)
            end
        else
            @user_info "Processing: $desc ($(length(items)) items)"
            for item in items
                $(esc(block))(item)
            end
            @user_info "Completed: $desc"
        end
    end
end

# Usage:
# @smart_iterate_block files "Processing MS files" item -> begin
#     process_file(item)
# end
```

---

## Warnings Tracking System

### Overview

The warnings tracking system provides comprehensive warning management for Pioneer.jl by automatically collecting, categorizing, and reporting warnings that occur during SearchDIA and BuildSpecLib execution. At the end of each run, users receive intelligent notifications based on warning severity and count, with detailed reports saved to dedicated warning files.

### Key Features

1. **Automatic Warning Collection**: Captures all warnings during execution without user intervention
2. **Intelligent Categorization**: Groups warnings by type (Data Quality, Parameter Issues, Performance, etc.) with severity levels
3. **End-of-Run Summary**: Smart console notifications based on warning count and severity
4. **Detailed Warning Reports**: Complete reports saved to dedicated warnings file with timestamps and context
5. **Context-Aware Warnings**: Enhanced macros that automatically capture processing context

### 1. Warning Tracking Infrastructure

```julia
# src/utils/PioneerLogging/warnings/WarningTracker.jl

"""
Structure to hold individual warning information
"""
struct WarningRecord
    message::String
    source_module::Module
    source_file::String
    source_line::Int
    timestamp::DateTime
    context::Dict{Symbol, Any}  # Additional context like file being processed
    stack_trace::Vector{Base.StackTraces.StackFrame}
end

"""
Warning category for grouping similar warnings
"""
struct WarningCategory
    name::String
    pattern::Regex
    severity::Symbol  # :low, :medium, :high, :critical
    description::String
end

"""
Central warning tracking system
"""
mutable struct WarningTracker
    warnings::Vector{WarningRecord}
    categories::Vector{WarningCategory}
    warning_counts::Dict{String, Int}  # Count by warning message
    category_counts::Dict{String, Int}  # Count by category
    session_start::DateTime
    max_warnings::Int  # Limit to prevent memory issues
    enabled::Bool
end

function WarningTracker(; max_warnings::Int = 10000)
    return WarningTracker(
        WarningRecord[],
        get_default_warning_categories(),
        Dict{String, Int}(),
        Dict{String, Int}(),
        now(),
        max_warnings,
        true
    )
end

# Global warning tracker instance
const GLOBAL_WARNING_TRACKER = Ref{Union{WarningTracker, Nothing}}(nothing)

function get_warning_tracker()::Union{WarningTracker, Nothing}
    return GLOBAL_WARNING_TRACKER[]
end

function initialize_warning_tracker!(; max_warnings::Int = 10000)
    GLOBAL_WARNING_TRACKER[] = WarningTracker(; max_warnings = max_warnings)
end

function reset_warning_tracker!()
    tracker = get_warning_tracker()
    if tracker !== nothing
        empty!(tracker.warnings)
        empty!(tracker.warning_counts)
        empty!(tracker.category_counts)
        tracker.session_start = now()
    end
end
```

### 2. Warning Categories

```julia
# src/utils/PioneerLogging/warnings/WarningCategories.jl

"""
Define standard warning categories for Pioneer.jl
"""
function get_default_warning_categories()::Vector{WarningCategory}
    return [
        # Data Quality Warnings
        WarningCategory(
            "Low Quality Spectra",
            r"low quality|poor quality|insufficient peaks",
            :medium,
            "Spectra with insufficient data quality for reliable identification"
        ),
        
        WarningCategory(
            "Missing Data",
            r"missing|not found|file not found|empty",
            :high,
            "Missing required data files or empty datasets"
        ),
        
        # Parameter Issues
        WarningCategory(
            "Parameter Convergence",
            r"convergence|failed to converge|optimization",
            :high,
            "Parameter optimization did not converge properly"
        ),
        
        WarningCategory(
            "Tolerance Issues",
            r"tolerance|mass error|retention time",
            :medium,
            "Mass or retention time tolerance issues"
        ),
        
        # Performance Warnings
        WarningCategory(
            "Memory Issues",
            r"memory|allocation|GC|garbage|out of memory",
            :high,
            "Memory allocation or garbage collection issues"
        ),
        
        WarningCategory(
            "Performance",
            r"slow|timeout|exceeded.*time|performance",
            :medium,
            "Performance-related warnings"
        ),
        
        # Library and Model Issues
        WarningCategory(
            "Library Issues",
            r"library|spectral.*library|koina|model",
            :high,
            "Issues with spectral library or prediction models"
        ),
        
        # File I/O Issues
        WarningCategory(
            "File I/O",
            r"file.*error|write.*error|read.*error|permission",
            :high,
            "File input/output errors"
        ),
        
        # Statistical Issues
        WarningCategory(
            "Statistical",
            r"FDR|false discovery|p-value|statistical",
            :medium,
            "Statistical analysis warnings"
        ),
        
        # Configuration Issues
        WarningCategory(
            "Configuration",
            r"config|parameter|setting|invalid.*value",
            :high,
            "Configuration or parameter setting issues"
        ),
        
        # General catch-all
        WarningCategory(
            "General",
            r".*",  # Matches everything as fallback
            :low,
            "General warnings not categorized elsewhere"
        )
    ]
end

function categorize_warning(message::String, categories::Vector{WarningCategory})::WarningCategory
    for category in categories
        if occursin(category.pattern, lowercase(message))
            return category
        end
    end
    # Should never reach here due to "General" catch-all, but just in case
    return categories[end]
end
```

### 3. Warning Capture Integration

```julia
# src/utils/PioneerLogging/warnings/WarningCapture.jl

"""
Custom logger that captures warnings for tracking while passing them through
"""
struct WarningCapturingLogger <: AbstractLogger
    base_logger::AbstractLogger
    tracker::WarningTracker
end

function WarningCapturingLogger(base_logger::AbstractLogger)
    tracker = get_warning_tracker()
    if tracker === nothing
        throw(ArgumentError("Warning tracker not initialized. Call initialize_warning_tracker!() first."))
    end
    return WarningCapturingLogger(base_logger, tracker)
end

function Logging.min_enabled_level(logger::WarningCapturingLogger)
    return Logging.min_enabled_level(logger.base_logger)
end

function Logging.shouldlog(logger::WarningCapturingLogger, level, _module, group, id)
    return Logging.shouldlog(logger.base_logger, level, _module, group, id)
end

function Logging.handle_message(
    logger::WarningCapturingLogger, 
    level, 
    message, 
    _module, 
    group, 
    id, 
    file, 
    line; 
    kwargs...
)
    # Capture warnings for tracking
    if level >= Logging.Warn && logger.tracker.enabled
        capture_warning!(logger.tracker, message, _module, file, line; kwargs...)
    end
    
    # Pass through to base logger
    return Logging.handle_message(
        logger.base_logger, level, message, _module, group, id, file, line; kwargs...
    )
end

function capture_warning!(
    tracker::WarningTracker, 
    message::String, 
    source_module::Module, 
    source_file::String, 
    source_line::Int;
    kwargs...
)
    # Prevent memory issues by limiting stored warnings
    if length(tracker.warnings) >= tracker.max_warnings
        # Remove oldest 10% of warnings
        num_to_remove = max(1, tracker.max_warnings ÷ 10)
        deleteat!(tracker.warnings, 1:num_to_remove)
    end
    
    # Create warning record
    warning = WarningRecord(
        message,
        source_module,
        source_file,
        source_line,
        now(),
        Dict{Symbol, Any}(kwargs...),
        Base.stacktrace()
    )
    
    push!(tracker.warnings, warning)
    
    # Update counts
    tracker.warning_counts[message] = get(tracker.warning_counts, message, 0) + 1
    
    # Categorize and count by category
    category = categorize_warning(message, tracker.categories)
    tracker.category_counts[category.name] = get(tracker.category_counts, category.name, 0) + 1
end
```

### 4. Enhanced Macros with Warning Context

```julia
# Addition to src/utils/PioneerLogging/macros/UserMacros.jl

"""
Enhanced warning macro that automatically adds context
"""
macro warn_with_context(msg, context...)
    quote
        # Add current processing context if available
        ctx_dict = Dict{Symbol, Any}()
        
        # Add provided context
        for (i, ctx) in enumerate($(map(esc, context)))
            if ctx isa Pair
                ctx_dict[ctx.first] = ctx.second
            else
                ctx_dict[Symbol("context_$i")] = ctx
            end
        end
        
        # Add automatic context from current execution
        if haskey(task_local_storage(), :current_file)
            ctx_dict[:processing_file] = task_local_storage(:current_file)
        end
        
        if haskey(task_local_storage(), :current_method)
            ctx_dict[:search_method] = task_local_storage(:current_method)
        end
        
        @logmsg Logging.Warn $(esc(msg)) ctx_dict...
    end
end

# Convenience macros for common warning scenarios
macro warn_file_processing(file_name, msg)
    quote
        @warn_with_context $(esc(msg)) processing_file=$(esc(file_name))
    end
end

macro warn_parameter_tuning(msg, params...)
    quote
        @warn_with_context $(esc(msg)) context=:parameter_tuning $(map(esc, params)...)
    end
end

macro warn_data_quality(msg, quality_metrics...)
    quote
        @warn_with_context $(esc(msg)) context=:data_quality $(map(esc, quality_metrics)...)
    end
end
```

### 5. Warning Report Generation

```julia
# src/utils/PioneerLogging/warnings/WarningReports.jl

"""
Generate comprehensive warning reports
"""
function generate_warning_report()
    tracker = get_warning_tracker()
    if tracker === nothing || isempty(tracker.warnings)
        @user_info "✓ No warnings detected during execution."
        return
    end
    
    system = get_logging_system()
    if system === nothing
        @warn "Cannot generate warning report: no logging system active"
        return
    end
    
    # Create warnings file
    warnings_file = joinpath(system.config.log_dir, "warnings_$(system.session_id).txt")
    
    try
        open(warnings_file, "w") do io
            write_warning_report(io, tracker)
        end
        
        # Generate console summary
        generate_console_warning_summary(tracker, warnings_file)
        
    catch e
        @error "Failed to generate warning report" error=e file=warnings_file
    end
end

function generate_console_warning_summary(tracker::WarningTracker, warnings_file::String)
    total_warnings = length(tracker.warnings)
    
    if total_warnings == 0
        @user_info "✓ No warnings detected during execution."
        return
    end
    
    # Determine severity level
    critical_count = get(tracker.category_counts, "Critical", 0)
    high_count = sum(get(tracker.category_counts, cat, 0) for cat in keys(tracker.category_counts) 
                    if any(c -> c.name == cat && c.severity == :high, tracker.categories))
    
    # Generate appropriate message based on warning count and severity
    if total_warnings <= 2 && critical_count == 0 && high_count == 0
        @user_info "⚠ $total_warnings warning(s) detected. See details in: $warnings_file"
    elseif total_warnings <= 5 && critical_count == 0
        @warn "⚠ $total_warnings warnings detected during execution."
        print_brief_warning_summary(tracker)
        @user_info "📋 Complete warning report saved to: $warnings_file"
    else
        @warn "⚠ SIGNIFICANT WARNINGS DETECTED: $total_warnings total warnings"
        if critical_count > 0
            @warn "🚨 $critical_count CRITICAL warnings detected!"
        end
        if high_count > 0
            @warn "❗ $high_count HIGH severity warnings detected!"
        end
        
        print_brief_warning_summary(tracker)
        @warn "📋 IMPORTANT: Review the complete warning report at: $warnings_file"
        @user_info "Consider investigating warnings before proceeding with results."
    end
end

function print_brief_warning_summary(tracker::WarningTracker)
    @user_info "Warning Summary:"
    
    # Show top 3 most frequent warnings
    sorted_warnings = sort(collect(tracker.warning_counts), by=x->x[2], rev=true)
    for (i, (message, count)) in enumerate(sorted_warnings)
        if i > 3  # Show top 3
            break
        end
        
        # Truncate long messages
        display_message = length(message) > 60 ? message[1:57] * "..." : message
        @user_info "  • [$count×] $display_message"
    end
    
    if length(sorted_warnings) > 3
        remaining = length(sorted_warnings) - 3
        @user_info "  • ... and $remaining more warning types"
    end
end
```

### 6. Integration with PioneerLoggingSystem

```julia
# Update to src/utils/PioneerLogging/core/LoggingSystem.jl

function setup_pioneer_logging!(
    log_dir::String = "logs";
    console_level::LogLevel = PIONEER_USER,
    simplified_enabled::Bool = true,
    full_enabled::Bool = true,
    show_progress::Bool = true,
    track_warnings::Bool = true,  # New parameter
    max_warnings::Int = 10000,    # New parameter
    config_file::Union{String, Nothing} = nothing
)
    # Initialize warning tracking if enabled
    if track_warnings
        initialize_warning_tracker!(; max_warnings = max_warnings)
    end
    
    # ... existing configuration loading code ...
    
    # Create new system
    system = PioneerLoggingSystem(config)
    
    # Wrap the tee logger with warning capturing if tracking is enabled
    final_logger = if track_warnings
        WarningCapturingLogger(system.tee_logger)
    else
        system.tee_logger
    end
    
    GLOBAL_LOGGING_SYSTEM[] = system
    
    # Set as global logger
    global_logger(final_logger)
    
    @info "Pioneer logging system initialized" session_id=system.session_id log_dir=config.log_dir track_warnings=track_warnings
    
    return system
end

# Update cleanup function
function cleanup_logging_system!()
    system = GLOBAL_LOGGING_SYSTEM[]
    if system !== nothing
        # Generate warning report before cleanup
        generate_warning_report()
        
        # Close file handles
        if system.simplified_file !== nothing
            close(system.simplified_file)
        end
        if system.full_file !== nothing
            close(system.full_file)
        end
        
        @debug "Cleaned up logging system" session_id=system.session_id
        GLOBAL_LOGGING_SYSTEM[] = nothing
    end
    
    # Reset warning tracker
    GLOBAL_WARNING_TRACKER[] = nothing
end
```

### 7. Usage in Search Methods

```julia
# Example integration in ParameterTuningSearch.jl
function process_file!(results, params, context, file_idx, spectra)
    @trace process_file! file_idx=file_idx
    
    # Set processing context for warnings
    task_local_storage(:current_file, "file_$file_idx")
    task_local_storage(:current_method, "ParameterTuningSearch")
    
    try
        # Existing processing logic...
        
        if length(good_psms) < minimum_required_psms
            @warn_parameter_tuning "Insufficient PSMs for reliable parameter tuning" file_idx=file_idx psm_count=length(good_psms) required=minimum_required_psms
        end
        
        if convergence_failed
            @warn_parameter_tuning "Parameter optimization failed to converge" file_idx=file_idx iterations=max_iterations final_error=final_error
        end
        
    catch e
        @warn_file_processing "file_$file_idx" "Failed to process file during parameter tuning: $e"
        rethrow(e)
    finally
        # Clean up context
        delete!(task_local_storage(), :current_file)
        delete!(task_local_storage(), :current_method)
    end
end
```

### 8. Expected Console Output Examples

```bash
# Case 1: No warnings
$ julia SearchDIA.jl params.json
...
16:32:15 SearchDIA completed successfully
16:32:15 ✓ No warnings detected during execution.

# Case 2: Few minor warnings  
$ julia SearchDIA.jl params.json
...
16:32:15 SearchDIA completed successfully
16:32:15 ⚠ 3 warning(s) detected. See details in: logs/warnings_a1b2c3d4.txt

# Case 3: Multiple warnings
$ julia SearchDIA.jl params.json
...
16:32:15 SearchDIA completed successfully
16:32:15 ⚠ 12 warnings detected during execution.
16:32:15 Warning Summary:
16:32:15   • [4×] Low quality spectra detected in file processing
16:32:15   • [3×] Parameter optimization failed to converge for some files
16:32:15   • [2×] Insufficient PSMs for reliable parameter tuning
16:32:15   • ... and 3 more warning types
16:32:15 📋 Complete warning report saved to: logs/warnings_a1b2c3d4.txt

# Case 4: Critical warnings
$ julia SearchDIA.jl params.json
...
16:32:15 SearchDIA completed successfully
16:32:15 ⚠ SIGNIFICANT WARNINGS DETECTED: 25 total warnings
16:32:15 🚨 2 CRITICAL warnings detected!
16:32:15 ❗ 8 HIGH severity warnings detected!
16:32:15 Warning Summary:
16:32:15   • [8×] Missing spectral library entries for precursors
16:32:15   • [6×] Memory allocation exceeded expected limits
16:32:15   • [4×] File I/O errors during data writing
16:32:15   • ... and 5 more warning types
16:32:15 ⚠ 📋 IMPORTANT: Review the complete warning report at: logs/warnings_a1b2c3d4.txt
16:32:15 Consider investigating warnings before proceeding with results.
```

---

## File Management

### 1. Log File Utilities

```julia
# src/utils/PioneerLogging/utils/FileManagement.jl

"""
Utilities for managing log files, rotation, and cleanup
"""

function ensure_log_directory(log_dir::String)
    if !isdir(log_dir)
        mkpath(log_dir)
        @debug "Created log directory" path=log_dir
    end
end

function get_log_file_size(filepath::String)::Int
    try
        return stat(filepath).size
    catch
        return 0
    end
end

function rotate_log_file(filepath::String, max_size::Int, max_files::Int)
    if !isfile(filepath) || get_log_file_size(filepath) < max_size
        return filepath
    end
    
    # Generate rotated filename
    base_name = splitext(filepath)[1]
    extension = splitext(filepath)[2]
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    rotated_path = "$(base_name)_$(timestamp)$(extension)"
    
    # Move current file to rotated name
    mv(filepath, rotated_path)
    
    # Clean up old files
    cleanup_old_log_files(dirname(filepath), basename(base_name), max_files)
    
    @debug "Rotated log file" old=filepath new=rotated_path
    return filepath  # Return original path for new log
end

function cleanup_old_log_files(log_dir::String, base_name::String, max_files::Int)
    pattern = Regex("^$(base_name)_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}-\\d{2}")
    
    old_files = filter(readdir(log_dir)) do filename
        occursin(pattern, filename)
    end
    
    if length(old_files) <= max_files
        return
    end
    
    # Sort by modification time (newest first)
    sort!(old_files, by = f -> stat(joinpath(log_dir, f)).mtime, rev = true)
    
    # Remove oldest files
    for old_file in old_files[max_files+1:end]
        old_path = joinpath(log_dir, old_file)
        try
            rm(old_path)
            @debug "Removed old log file" file=old_path
        catch e
            @warn "Failed to remove old log file" file=old_path error=e
        end
    end
end

function compress_old_logs(log_dir::String, days_old::Int=7)
    cutoff_time = now() - Day(days_old)
    
    for filename in readdir(log_dir)
        filepath = joinpath(log_dir, filename)
        if isfile(filepath) && 
           endswith(filename, ".log") && 
           stat(filepath).mtime < datetime2unix(cutoff_time)
            
            compressed_path = filepath * ".gz"
            if !isfile(compressed_path)
                try
                    # Use gzip to compress old log
                    run(`gzip $filepath`)
                    @debug "Compressed old log file" original=filepath compressed=compressed_path
                catch e
                    @warn "Failed to compress log file" file=filepath error=e
                end
            end
        end
    end
end
```

### 2. Session Management

```julia
# src/utils/PioneerLogging/utils/SessionManagement.jl

"""
Manage multiple SearchDIA/BuildSpecLib sessions with different logging configs
"""

struct LoggingSession
    id::String
    config::PioneerLoggingConfig
    system::PioneerLoggingSystem
    started_at::DateTime
end

const LOGGING_SESSIONS = Ref{Vector{LoggingSession}}(LoggingSession[])

function start_logging_session(
    log_dir::String = "logs";
    console_level::LogLevel = PIONEER_USER,
    session_name::Union{String, Nothing} = nothing,
    kwargs...
)::String
    # Generate session ID
    session_id = session_name !== nothing ? session_name : string(uuid4())[1:8]
    
    # Create unique log directory for this session
    session_log_dir = joinpath(log_dir, "session_$(session_id)")
    
    # Setup logging for this session
    system = setup_pioneer_logging!(session_log_dir; console_level=console_level, kwargs...)
    
    # Create session record
    session = LoggingSession(session_id, system.config, system, now())
    push!(LOGGING_SESSIONS[], session)
    
    @user_info "Started logging session" session_id=session_id log_dir=session_log_dir console_level=level_name(console_level)
    
    return session_id
end

function end_logging_session(session_id::String)
    sessions = LOGGING_SESSIONS[]
    session_idx = findfirst(s -> s.id == session_id, sessions)
    
    if session_idx === nothing
        @warn "Session not found" session_id=session_id
        return
    end
    
    session = sessions[session_idx]
    duration = now() - session.started_at
    
    @user_info "Ending logging session" session_id=session_id duration=duration
    
    # Cleanup the system
    cleanup_logging_system!()
    
    # Remove from active sessions
    deleteat!(sessions, session_idx)
end

function list_active_sessions()::Vector{String}
    return [s.id for s in LOGGING_SESSIONS[]]
end

function switch_to_session(session_id::String)
    sessions = LOGGING_SESSIONS[]
    session = findfirst(s -> s.id == session_id, sessions)
    
    if session === nothing
        @warn "Session not found" session_id=session_id
        return false
    end
    
    # Set the session's system as global logger
    global_logger(session.system.tee_logger)
    GLOBAL_LOGGING_SYSTEM[] = session.system
    
    @user_info "Switched to logging session" session_id=session_id
    return true
end

# Convenience functions for multiple SearchDIA runs
function run_with_logging(f, log_dir::String="logs"; console_level::LogLevel=PIONEER_USER, kwargs...)
    session_id = start_logging_session(log_dir; console_level=console_level, kwargs...)
    try
        return f()
    finally
        end_logging_session(session_id)
    end
end

# Usage:
# run_with_logging("logs/run1"; console_level=PIONEER_DEBUG1) do
#     SearchDIA("params1.json")
# end
#
# run_with_logging("logs/run2"; console_level=PIONEER_USER) do
#     SearchDIA("params2.json")  
# end
```

---

## Converting Existing Logging Patterns

### Overview

Pioneer.jl currently uses a mix of logging patterns that need systematic conversion to the new unified system:

1. **dual_println/dual_print**: Local functions that write to both console and simplified log
2. **@info macros**: Julia's built-in info logging (currently goes to stdout) 
3. **@warn/@error macros**: Julia's warning and error logging
4. **println**: Direct console output (rare, mostly in utility functions)
5. **Progress bars**: ProgressBars.jl for long-running operations

### Conversion Rules and Mapping

#### Rule 1: dual_println → @user_info (Simplified Log)
**Current Pattern:**
```julia
# In SearchDIA.jl and BuildSpecLib.jl
open(log_path, "w") do log_file
    function dual_println(args...)
        println(stdout, args...)
        println(log_file, args...)
    end
    
    dual_println("Starting SearchDIA")
    dual_println("Output directory: ", params.paths[:results])
    dual_println("Processing ", length(files), " files")
end
```

**New Pattern:**
```julia
# Setup logging once at start
setup_pioneer_logging!(log_dir; console_level=PIONEER_USER)

# Convert all dual_println calls
@user_info "Starting SearchDIA"
@user_info "Output directory: $(params.paths[:results])"
@user_info "Processing $(length(files)) files"
```

**Additional Requirements:**
- Add Pioneer version to simplified log header
- Include session timestamp and configuration summary

#### Rule 2: Major Process @info → @user_info (Console + Simplified Log)
**Current Pattern:**
```julia
@info "Executing Parameter Tuning..."
@info "Executing NCE Tuning..."
@info "Executing First Pass Search..."
@info "Executing Scoring..."
```

**New Pattern:**
```julia
@user_info "Starting Parameter Tuning Search"
@user_info "Starting NCE Tuning Search" 
@user_info "Starting First Pass Search"
@user_info "Starting Scoring Search"
```

**Identification Criteria for Major Processes:**
- Top-level search method initiation (9 main search methods)
- Library loading/building phases
- Major pipeline transitions
- Final results generation

#### Rule 3: Detailed @info → @debug_l2 (Full Log Only)
**Current Pattern:**
```julia
@info "Loading Parameters..."
@info "Initializing Search Context..."
@info "Step 1 completed in $(round(step1_time, digits=2)) seconds"
@info "Found $(n_matches) spectral matches"
@info "Processing file $file_idx of $total_files"
```

**New Pattern:**
```julia
@debug_l2 "Loading Parameters..."
@debug_l2 "Initializing Search Context..."
@debug_l2 "Step 1 completed in $(round(step1_time, digits=2)) seconds"
@debug_l2 "Found $(n_matches) spectral matches"
@debug_l2 "Processing file $file_idx of $total_files"
```

**Identification Criteria for Detailed Info:**
- Step-by-step progress within major processes
- Detailed timing information
- File-by-file processing status
- Technical metrics and counts
- Intermediate calculations

#### Rule 4: @warn/@error → Enhanced with Context
**Current Pattern:**
```julia
@warn "Low quality spectra detected"
@error "Failed to process file: $filename"
```

**New Pattern:**
```julia
@warn_data_quality "Low quality spectra detected" file_idx=file_idx spectrum_count=count
@warn_file_processing filename "Failed to process file: $error_msg"
```

#### Rule 5: Progress Bars → Smart Integration
**Current Pattern:**
```julia
pbar = ProgressBar(files; description="Processing files")
for file in pbar
    # process file
end
```

**New Pattern:**
```julia
for file in @smart_iterate(files, "Processing MS files")
    # process file
end
# Automatic cleanup handled by macro
```

### Detailed Conversion Mapping

#### SearchDIA.jl Main Function
**Before:**
```julia
function SearchDIA(params_path::String)
    # Setup dual_println
    open(log_path, "w") do log_file
        function dual_println(args...)
            println(stdout, args...)
            println(log_file, args...)
        end
        
        dual_println("\n", repeat("=", 90))
        dual_println("Sarting SearchDIA")
        dual_println(repeat("=", 90))
        dual_println("\nStarting search at: ", Dates.now())
        dual_println("Output directory: ", params.paths[:results])
        
        @info "Loading Parameters..."
        @info "Loading Spectral Library..."
        @info "Initializing Search Context..."
        
        for (name, search) in searches
            @info "Executing $name..."
        end
        
        print_performance_report(timings, MS_TABLE_PATHS, SEARCH_CONTEXT, dual_println)
    end
end
```

**After:**
```julia
function SearchDIA(params_path::String)
    # Parse parameters to get log directory
    params = parse_pioneer_parameters(params_path)
    log_dir = params.paths[:results]
    
    # Setup integrated logging system
    setup_pioneer_logging!(log_dir; 
                          console_level=PIONEER_USER, 
                          track_warnings=true)
    
    # Add version info to simplified log
    write_simplified_log_header(params)
    
    @user_info "Starting SearchDIA"
    @user_info "Starting search at: $(Dates.now())"
    @user_info "Output directory: $(params.paths[:results])"
    
    @debug_l2 "Loading Parameters..."
    @debug_l2 "Loading Spectral Library..."  
    @debug_l2 "Initializing Search Context..."
    
    for (name, search) in searches
        @user_info "Starting $name"
        @debug_l1 "Executing $name..."
    end
    
    print_performance_report(timings, MS_TABLE_PATHS, SEARCH_CONTEXT)
    
    # Automatic warning report generation in cleanup
    cleanup_logging_system!()
end
```

#### Search Methods Conversion
**Before (ScoringSearch.jl):**
```julia
@info "Step 1 completed in $(round(step1_time, digits=2)) seconds"
@info "Step 2 completed in $(round(step2_time, digits=2)) seconds"
@info "Found $(length(psms)) PSMs passing FDR"
```

**After (ScoringSearch.jl):**
```julia
@debug_l2 "Step 1 completed in $(round(step1_time, digits=2)) seconds"
@debug_l2 "Step 2 completed in $(round(step2_time, digits=2)) seconds"  
@debug_l1 "Found $(length(psms)) PSMs passing FDR"
```

#### Warning Integration in Search Methods
**Before:**
```julia
if length(psms) < minimum_required
    @warn "Insufficient PSMs for reliable analysis"
end
```

**After:**
```julia
# Set processing context
task_local_storage(:current_method, "ScoringSearch")
task_local_storage(:current_file, filename)

if length(psms) < minimum_required
    @warn_parameter_tuning "Insufficient PSMs for reliable analysis" psm_count=length(psms) required=minimum_required
end

# Clean up context when done
delete!(task_local_storage(), :current_method)
delete!(task_local_storage(), :current_file)
```

#### Performance Report Conversion
**Before:**
```julia
function print_performance_report(timings, ms_table_paths, search_context, println_func)
    println_func("\n" * repeat("=", 90))
    println_func("DIA Search Performance Report") 
    println_func(repeat("=", 90))
    # ... detailed reporting
end
```

**After:**
```julia
function print_performance_report(timings, ms_table_paths, search_context)
    @user_info "DIA Search Performance Report"
    @user_info repeat("=", 50)
    
    # Summary info goes to simplified log
    @user_info "Total runtime: $(round(total_time, digits=2)) seconds"
    @user_info "Peak memory: $(round(peak_memory/1024^3, digits=2)) GB"
    
    # Detailed breakdown goes to debug log
    @debug_l1 "Detailed Step Analysis:"
    for step in sorted_steps
        timing = timings[step]
        @debug_l1 "$(rpad(step, 30)): $(round(timing.time, digits=2))s"
    end
end
```

### Simplified Log Header Format
**New requirement**: Include Pioneer version and configuration at top of simplified log

```julia
function write_simplified_log_header(params::PioneerParameters)
    # This gets written to simplified log only
    @user_info repeat("=", 80)
    @user_info "PIONEER.JL v$(get_pioneer_version()) - DIA SEARCH LOG"
    @user_info repeat("=", 80)
    @user_info "Session started: $(Dates.now())"
    @user_info "Configuration file: $(params.config_path)"
    @user_info "Julia version: $(VERSION)"
    @user_info "Threads: $(Threads.nthreads())"
    @user_info "Platform: $(Sys.MACHINE)"
    @user_info repeat("=", 80)
    @user_info ""  # Blank line for readability
end

function get_pioneer_version()::String
    try
        # Read from Project.toml or VERSION file
        project_toml = joinpath(@__DIR__, "../../../Project.toml")
        if isfile(project_toml)
            content = read(project_toml, String)
            m = match(r"version\s*=\s*\"([^\"]+)\"", content)
            return m !== nothing ? m[1] : "unknown"
        end
        return "unknown"
    catch
        return "unknown"
    end
end
```

### File-by-File Conversion Plan

#### Phase 1 Files (Core Entry Points):
1. **src/Routines/SearchDIA.jl** - Main SearchDIA function
   - Convert all dual_println → @user_info  
   - Convert major process @info → @user_info
   - Convert detailed @info → @debug_l2
   - Add simplified log header
   - Integrate warnings system

2. **src/Routines/BuildSpecLib.jl** - Library building
   - Same pattern as SearchDIA.jl
   - Focus on major phases vs detailed steps

#### Phase 2 Files (Search Methods):
3. **SearchMethods/ParameterTuningSearch/**
   - Major: "Starting Parameter Tuning Search" → @user_info
   - Details: All current @info → @debug_l2
   - Warnings: Add context-aware warning macros

4. **SearchMethods/FirstPassSearch/**
   - Major: "Starting First Pass Search" → @user_info
   - Details: File processing info → @debug_l1
   - Progress: Convert to @smart_iterate

5. **SearchMethods/ScoringSearch/**
   - Major: "Starting Scoring Search" → @user_info
   - Steps: "Step N completed" → @debug_l2
   - Important: PSM counts → @debug_l1
   - ML: Model training progress → @debug_l1

6. **SearchMethods/MaxLFQSearch/**
   - Major: "Starting MaxLFQ Search" → @user_info
   - Quantification: Results summary → @debug_l1

#### Phase 3 Files (Utilities):
7. **CommonSearchUtils/** - Peak matching, RT alignment
   - Most current @info → @debug_l2
   - Warnings: Add appropriate context

8. **utils/** - General utilities
   - Keep existing @info for now
   - Add warnings context where appropriate

### Migration Testing Strategy

#### Validate Conversion Success:
1. **Console Output Test**: Should only show user messages, warnings, and progress
2. **Simplified Log Test**: Should contain clean timeline of major events
3. **Full Log Test**: Should contain all debugging information
4. **Warning Report Test**: Should capture and categorize warnings appropriately

#### Test Script:
```julia
function test_logging_conversion()
    # Run SearchDIA with each verbosity level
    test_cases = [
        ("user", PIONEER_USER),
        ("debug1", PIONEER_DEBUG1), 
        ("debug2", PIONEER_DEBUG2),
        ("trace", PIONEER_TRACE)
    ]
    
    for (name, level) in test_cases
        log_dir = "test_logs_$name"
        setup_pioneer_logging!(log_dir; console_level=level)
        
        # Capture console output
        console_output = capture_stdout() do
            SearchDIA("test_params.json")
        end
        
        # Analyze outputs
        simplified_log = read(joinpath(log_dir, "pioneer_simple.log"), String)
        full_log = read(joinpath(log_dir, "pioneer_full.log"), String) 
        
        # Validate expectations for each level
        validate_console_output(console_output, level)
        validate_simplified_log(simplified_log)
        validate_full_log(full_log, level)
        
        cleanup_logging_system!()
    end
end
```

### Backward Compatibility During Migration

#### Deprecation Strategy:
```julia
# Provide transitional dual_println for gradual migration
function dual_println(args...)
    Base.depwarn("dual_println is deprecated, use @user_info instead", :dual_println)
    @user_info join(string.(args), " ")
end

function dual_print(args...)
    Base.depwarn("dual_print is deprecated, use @user_info instead", :dual_print)  
    @user_info join(string.(args), "")
end
```

This comprehensive conversion plan ensures that:
1. All existing output patterns are systematically converted
2. The simplified log contains a clean timeline with version info
3. Console output is appropriate for each verbosity level
4. Warnings are properly tracked and contextualized
5. Major process milestones are clearly visible to users
6. Detailed debugging information is available when needed

---

## Migration Strategy

### Phase 1: Infrastructure Setup (Week 1-2)

#### Week 1: Core Infrastructure
1. **Create Module Structure**
   ```bash
   mkdir -p src/utils/PioneerLogging/{core,loggers,formatters,macros,progress,utils,testing}
   ```

2. **Implement Core Components**
   - [ ] Create `LogLevels.jl` with custom levels
   - [ ] Implement `Configuration.jl` with config loading
   - [ ] Build `LoggingSystem.jl` with TeeLogger setup
   - [ ] Add basic error handling and validation

3. **Basic Testing**
   - [ ] Unit tests for configuration loading
   - [ ] Tests for log level parsing
   - [ ] Validation of TeeLogger setup

#### Week 2: Logger Implementation
1. **Individual Loggers**
   - [ ] Implement `ConsoleLogger.jl` with runtime level control
   - [ ] Build `SimplifiedLogger.jl` with file output
   - [ ] Create `FullDebugLogger.jl` with rotation

2. **Integration Testing**
   - [ ] Test message routing to all three destinations
   - [ ] Verify independent level filtering
   - [ ] Test file creation and writing

### Phase 2: Formatters and Macros (Week 3-4)

#### Week 3: Message Formatting
1. **Formatter Implementation**
   - [ ] Create `ConsoleFormatter.jl` for clean console output
   - [ ] Build `SimpleFormatter.jl` for dual_println equivalent
   - [ ] Implement `DebugFormatter.jl` with full metadata

2. **Formatter Testing**
   - [ ] Test message formatting consistency
   - [ ] Verify timestamp and metadata inclusion
   - [ ] Test color output (if enabled)

#### Week 4: Custom Macros
1. **Macro Implementation**
   - [ ] Create `@user_info` and related user macros
   - [ ] Implement `@debug_l1`, `@debug_l2` debug macros
   - [ ] Build `@trace`, `@trace_exit` tracing macros

2. **Macro Testing**
   - [ ] Test macro expansion and level assignment
   - [ ] Verify argument handling and performance
   - [ ] Test conditional debug macros

### Phase 3: Progress Integration (Week 5-6)

#### Week 5: Progress Bar Coordination
1. **Progress Integration**
   - [ ] Implement progress-aware logging
   - [ ] Create message deferral system
   - [ ] Build `with_progress_logging` function

2. **Smart Iteration**
   - [ ] Implement `@smart_iterate` macro
   - [ ] Create automatic progress/logging selection
   - [ ] Add progress bar cleanup handling

#### Week 6: Testing and Refinement
1. **Integration Testing**
   - [ ] Test progress bar + logging coordination
   - [ ] Verify message deferral and flushing
   - [ ] Test smart iteration in real scenarios

### Phase 4: File Management and Advanced Features (Week 7-8)

#### Week 7: File Management and Warnings System
1. **File Utilities**
   - [ ] Implement log rotation and cleanup
   - [ ] Create session management system
   - [ ] Add compression for old logs

2. **Warnings System**
   - [ ] Implement warning tracking infrastructure
   - [ ] Create warning categorization system
   - [ ] Build warning capture logger
   - [ ] Implement warning report generation

3. **Session Support**
   - [ ] Enable multiple SearchDIA runs with different configs
   - [ ] Implement session switching
   - [ ] Add session cleanup and management

#### Week 8: Performance and Polish
1. **Performance Optimization**
   - [ ] Profile logging performance impact
   - [ ] Implement lazy evaluation for expensive debug messages
   - [ ] Optimize formatter performance

2. **Documentation and Examples**
   - [ ] Create comprehensive usage examples
   - [ ] Document all macros and functions
   - [ ] Write migration guide

### Phase 5: Migration and Deployment (Week 9-10)

#### Week 9: Code Migration
1. **Backward Compatibility**
   - [ ] Create `dual_println` compatibility functions
   - [ ] Add deprecation warnings
   - [ ] Implement gradual migration helpers

2. **SearchDIA Integration**
   - [ ] Update main SearchDIA.jl to use new logging
   - [ ] Migrate ParameterTuningSearch to new macros
   - [ ] Update other search methods
   - [ ] Add warning tracking calls throughout codebase

#### Week 10: Testing and Validation
1. **End-to-End Testing**
   - [ ] Run full SearchDIA pipeline with new logging
   - [ ] Verify all log files are created correctly
   - [ ] Test different verbosity levels
   - [ ] Validate warning tracking and reporting

2. **Performance Validation**
   - [ ] Benchmark before/after performance
   - [ ] Ensure no significant performance regression
   - [ ] Validate memory usage
   - [ ] Test warning system performance impact

---

## Testing Plan

### Unit Tests

```julia
# test/UnitTests/test_pioneer_logging.jl

using Test, PioneerLogging, LoggingExtras, Logging

@testset "PioneerLogging Tests" begin
    
    @testset "Log Levels" begin
        @test PIONEER_USER < Logging.Info
        @test PIONEER_DEBUG1 > Logging.Info  
        @test PIONEER_DEBUG2 > PIONEER_DEBUG1
        @test PIONEER_TRACE > PIONEER_DEBUG2
        
        @test parse_level("user") == PIONEER_USER
        @test parse_level("debug1") == PIONEER_DEBUG1
        @test_throws Exception parse_level("invalid")
        
        @test level_name(PIONEER_USER) == "USER"
        @test level_name(PIONEER_DEBUG1) == "DEBUG1"
    end
    
    @testset "Configuration" begin
        # Test default configuration
        config = PioneerLoggingConfig()
        @test config.console_level == PIONEER_USER
        @test config.log_dir == "logs"
        @test config.show_progress == true
        
        # Test JSON loading
        mktempdir() do tmpdir
            config_file = joinpath(tmpdir, "test_config.json")
            write(config_file, """
            {
                "logging": {
                    "console_level": "debug1",
                    "log_dir": "/tmp/test_logs",
                    "show_progress": false
                }
            }
            """)
            
            config = load_config_from_json(config_file)
            @test config.console_level == PIONEER_DEBUG1
            @test config.log_dir == "/tmp/test_logs"
            @test config.show_progress == false
        end
    end
    
    @testset "Logging System Setup" begin
        mktempdir() do tmpdir
            # Test system creation
            system = setup_pioneer_logging!(tmpdir; console_level=PIONEER_DEBUG1)
            @test system isa PioneerLoggingSystem
            @test system.config.log_dir == tmpdir
            @test system.console_level[] == PIONEER_DEBUG1
            
            # Test log directory creation
            @test isdir(tmpdir)
            
            # Test global logger is set
            @test global_logger() == system.tee_logger
            
            cleanup_logging_system!()
        end
    end
    
    @testset "Runtime Console Level Control" begin
        mktempdir() do tmpdir
            setup_pioneer_logging!(tmpdir; console_level=PIONEER_USER)
            
            # Test level changes
            set_console_level!(PIONEER_DEBUG1)
            system = get_logging_system()
            @test system.console_level[] == PIONEER_DEBUG1
            
            set_console_level!(PIONEER_TRACE)
            @test system.console_level[] == PIONEER_TRACE
            
            cleanup_logging_system!()
        end
    end
    
    @testset "Custom Macros" begin
        # Capture log output
        io = IOBuffer()
        logger = ConsoleLogger(io, PIONEER_TRACE)
        
        with_logger(logger) do
            @user_info "test user message"
            @debug_l1 "test debug message"
            @debug_l2 "test detailed debug"
            @trace "test_function"
        end
        
        output = String(take!(io))
        @test contains(output, "test user message")
        @test contains(output, "test debug message")
        @test contains(output, "test detailed debug")
        @test contains(output, "test_function")
    end
    
    @testset "Message Filtering" begin
        mktempdir() do tmpdir
            setup_pioneer_logging!(tmpdir; 
                console_level=PIONEER_USER,
                simplified_enabled=true,
                full_enabled=true
            )
            
            # Test messages at different levels
            @user_info "user message"
            @debug_l1 "debug message"
            @debug_l2 "detailed debug"
            
            # Allow time for file writing
            sleep(0.1)
            
            # Check simplified log (should only have user messages)
            simplified_files = filter(f -> contains(f, "simple"), readdir(tmpdir))
            @test length(simplified_files) == 1
            
            simplified_content = read(joinpath(tmpdir, simplified_files[1]), String)
            @test contains(simplified_content, "user message")
            @test !contains(simplified_content, "debug message")
            
            # Check full log (should have all messages)
            full_files = filter(f -> contains(f, "full"), readdir(tmpdir))
            @test length(full_files) == 1
            
            full_content = read(joinpath(tmpdir, full_files[1]), String)
            @test contains(full_content, "user message")
            @test contains(full_content, "debug message")
            @test contains(full_content, "detailed debug")
            
            cleanup_logging_system!()
        end
    end
    
    @testset "Progress Integration" begin
        mktempdir() do tmpdir
            setup_pioneer_logging!(tmpdir; show_progress=true)
            
            # Test progress detection
            @test should_show_progress(100) == true
            @test should_show_progress(5) == false
            
            # Test progress state management
            @test is_progress_active() == false
            set_progress_active!(true)
            @test is_progress_active() == true
            set_progress_active!(false)
            @test is_progress_active() == false
            
            cleanup_logging_system!()
        end
    end
    
    @testset "Session Management" begin
        mktempdir() do tmpdir
            # Test session creation
            session1 = start_logging_session(tmpdir; console_level=PIONEER_USER, session_name="test1")
            @test session1 == "test1"
            @test "test1" in list_active_sessions()
            
            # Test session switching
            session2 = start_logging_session(tmpdir; console_level=PIONEER_DEBUG1, session_name="test2")
            @test length(list_active_sessions()) == 2
            
            @test switch_to_session("test1") == true
            system = get_logging_system()
            @test system.config.console_level == PIONEER_USER
            
            @test switch_to_session("test2") == true
            system = get_logging_system()
            @test system.config.console_level == PIONEER_DEBUG1
            
            # Test session cleanup
            end_logging_session("test1")
            @test !("test1" in list_active_sessions())
            
            end_logging_session("test2")
            @test length(list_active_sessions()) == 0
        end
    end
    
    @testset "Warning Tracking" begin
        mktempdir() do tmpdir
            setup_pioneer_logging!(tmpdir; track_warnings=true, max_warnings=100)
            
            # Test warning capture
            @warn "Test warning message"
            @warn "Parameter convergence failed"
            @warn "Low quality spectra detected"
            
            tracker = get_warning_tracker()
            @test tracker !== nothing
            @test length(tracker.warnings) == 3
            @test haskey(tracker.warning_counts, "Test warning message")
            @test tracker.warning_counts["Test warning message"] == 1
            
            # Test warning categorization
            @test haskey(tracker.category_counts, "Parameter Convergence")
            @test haskey(tracker.category_counts, "Low Quality Spectra") 
            
            # Test warning report generation
            generate_warning_report()
            
            # Check that warnings file was created
            warnings_files = filter(f -> startswith(f, "warnings_"), readdir(tmpdir))
            @test length(warnings_files) == 1
            
            warnings_content = read(joinpath(tmpdir, warnings_files[1]), String)
            @test contains(warnings_content, "PIONEER.JL WARNING REPORT")
            @test contains(warnings_content, "Test warning message")
            @test contains(warnings_content, "Parameter convergence failed")
            
            cleanup_logging_system!()
        end
    end
    
    @testset "Warning Context Macros" begin
        mktempdir() do tmpdir
            setup_pioneer_logging!(tmpdir; track_warnings=true)
            
            # Test context-aware warning macros
            task_local_storage(:current_file, "test_file.mzML")
            task_local_storage(:current_method, "TestMethod")
            
            @warn_with_context "Test warning with context" test_param=42
            @warn_parameter_tuning "Optimization failed" iterations=100 error=0.05
            @warn_data_quality "Poor spectrum quality" peak_count=5 snr=1.2
            
            # Clean up context
            delete!(task_local_storage(), :current_file)
            delete!(task_local_storage(), :current_method)
            
            tracker = get_warning_tracker()
            @test length(tracker.warnings) == 3
            
            # Verify context was captured
            first_warning = tracker.warnings[1]
            @test haskey(first_warning.context, :processing_file)
            @test first_warning.context[:processing_file] == "test_file.mzML"
            @test haskey(first_warning.context, :search_method)
            @test first_warning.context[:search_method] == "TestMethod"
            
            cleanup_logging_system!()
        end
    end
end
```

### Integration Tests

```julia
# test/IntegrationTests/test_searchdia_logging.jl

@testset "SearchDIA Logging Integration" begin
    mktempdir() do tmpdir
        # Test with different verbosity levels
        for (level_name, level) in [("user", PIONEER_USER), ("debug1", PIONEER_DEBUG1), ("trace", PIONEER_TRACE)]
            @testset "Console Level: $level_name" begin
                log_dir = joinpath(tmpdir, "test_$level_name")
                
                run_with_logging(log_dir; console_level=level) do
                    # Run a small test case
                    SearchDIA("data/ecoli_test/ecoli_test_params.json")
                end
                
                # Verify log files exist
                @test isdir(log_dir)
                log_files = readdir(log_dir)
                
                simplified_files = filter(f -> contains(f, "simple"), log_files)
                full_files = filter(f -> contains(f, "full"), log_files)
                
                @test length(simplified_files) >= 1
                @test length(full_files) >= 1
                
                # Verify content differences
                simplified_content = read(joinpath(log_dir, simplified_files[1]), String)
                full_content = read(joinpath(log_dir, full_files[1]), String)
                
                # Full log should have more content
                @test length(split(full_content, "\n")) >= length(split(simplified_content, "\n"))
                
                # Both should contain user messages
                @test contains(simplified_content, "SearchDIA")
                @test contains(full_content, "SearchDIA")
            end
        end
    end
end

@testset "BuildSpecLib Logging Integration" begin
    mktempdir() do tmpdir
        log_dir = joinpath(tmpdir, "buildlib_test")
        
        run_with_logging(log_dir; console_level=PIONEER_DEBUG1) do
            # Run a small BuildSpecLib test
            BuildSpecLib("data/example_config/defaultBuildLibParams.json")
        end
        
        # Verify log files and content
        @test isdir(log_dir)
        log_files = readdir(log_dir)
        
        @test any(f -> contains(f, "simple"), log_files)
        @test any(f -> contains(f, "full"), log_files)
    end
end
```

---

## Performance Considerations

### 1. Lazy Evaluation

```julia
# Expensive debug computations should be evaluated lazily
macro debug_l2_lazy(msg, expensive_expr)
    quote
        if Logging.min_enabled_level(current_logger()) <= PIONEER_DEBUG2
            result = $(esc(expensive_expr))
            @logmsg PIONEER_DEBUG2 $(esc(msg)) computed_value=result
        end
    end
end

# Usage:
# @debug_l2_lazy "Complex analysis result" compute_expensive_statistics(data)
```

### 2. Conditional Logging

```julia
# Check if logging level is enabled before expensive operations
function should_log_debug1()::Bool
    return Logging.min_enabled_level(current_logger()) <= PIONEER_DEBUG1
end

function should_log_debug2()::Bool
    return Logging.min_enabled_level(current_logger()) <= PIONEER_DEBUG2
end

# Usage in performance-critical code:
if should_log_debug2()
    detailed_stats = compute_detailed_statistics(data)
    @debug_l2 "Detailed statistics" stats=detailed_stats
end
```

### 3. Asynchronous Logging (Optional)

```julia
# For high-throughput scenarios, implement async logging
struct AsyncLoggingBuffer
    messages::Channel{LogMessage}
    worker_task::Task
    base_logger::AbstractLogger
end

function AsyncLoggingBuffer(base_logger::AbstractLogger; buffer_size=1000)
    messages = Channel{LogMessage}(buffer_size)
    
    worker = @async begin
        while isopen(messages)
            try
                msg = take!(messages)
                handle_message(base_logger, msg.level, msg.message, 
                             msg._module, msg.group, msg.id, msg.file, msg.line; msg.kwargs...)
            catch e
                if e isa InvalidStateException && !isopen(messages)
                    break  # Channel was closed
                end
                # Log the error somehow without recursion
                println(stderr, "Async logging error: $e")
            end
        end
    end
    
    return AsyncLoggingBuffer(messages, worker, base_logger)
end
```

### 4. Memory Management

```julia
# Periodic cleanup of deferred messages
function cleanup_deferred_messages!()
    deferred = DEFERRED_MESSAGES[]
    if length(deferred) > 1000  # Prevent unbounded growth
        # Keep only the most recent 500 messages
        DEFERRED_MESSAGES[] = deferred[end-499:end]
        @warn "Deferred message buffer was full, dropped $(length(deferred) - 500) old messages"
    end
end

# Call this periodically in long-running operations
function periodic_logging_maintenance()
    cleanup_deferred_messages!()
    
    # Force garbage collection of closed file handles
    GC.gc()
end
```

---

## Usage Examples

### 1. SearchDIA Integration

```julia
# In SearchDIA.jl main function
function main_SearchDIA(argv=ARGS)::Cint
    # Parse arguments including logging options
    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "params_path"
            help = "Path to search parameters JSON file"
            arg_type = String
        "--console-verbosity"
            help = "Console verbosity level (user, debug1, debug2, trace)"
            arg_type = String
            default = "user"
        "--log-dir"
            help = "Directory for log files"
            arg_type = String
            default = "logs"
        "--no-progress"
            help = "Disable progress bars"
            action = :store_true
    end
    
    parsed_args = parse_args(argv, settings; as_symbols = true)
    
    # Setup logging
    console_level = parse_level(parsed_args[:console_verbosity])
    log_dir = parsed_args[:log_dir]
    show_progress = !parsed_args[:no_progress]
    
    setup_pioneer_logging!(log_dir; 
                          console_level=console_level, 
                          show_progress=show_progress,
                          track_warnings=true)
    
    @user_info "Starting SearchDIA analysis"
    @debug_l1 "Command line arguments" args=argv
    @debug_l1 "Console verbosity level" level=level_name(console_level)
    
    try
        SearchDIA(parsed_args[:params_path])
        @user_info "SearchDIA completed successfully"
        
        # Warning report is automatically generated in cleanup_logging_system!()
        
        return 0
    catch e
        @error "SearchDIA failed" exception=(e, catch_backtrace())
        return 1
    finally
        cleanup_logging_system!()  # This triggers warning report generation
    end
end
```

### 2. Search Method Implementation

```julia
# In ParameterTuningSearch.jl
function process_file!(results, params, context, file_idx, spectra)
    @trace process_file! file_idx=file_idx n_spectra=length(spectra)
    @debug_l1 "Processing file for parameter tuning" file_idx=file_idx
    
    # Set processing context for warnings
    task_local_storage(:current_file, "file_$file_idx")
    task_local_storage(:current_method, "ParameterTuningSearch")
    
    try
        # Process with progress tracking
        psms = Vector{SimpleScoredPSM}()
        
        for (scan_idx, spectrum) in @smart_iterate(spectra, "Scanning spectra") 
            @debug_l2 "Processing spectrum" scan_idx=scan_idx mz_count=length(spectrum.mz)
            
            # Core processing logic here
            matches = find_matches(spectrum, context.lib)
            
            if !isempty(matches)
                @debug_l2 "Found matches" scan_idx=scan_idx n_matches=length(matches)
                append!(psms, matches)
            elseif length(spectrum.mz) < minimum_peaks
                @warn_data_quality "Spectrum has insufficient peaks for reliable identification" scan_idx=scan_idx peak_count=length(spectrum.mz) minimum=minimum_peaks
            end
        end
        
        # Clean up progress state if active
        if is_progress_active()
            set_progress_active!(false)
            flush_deferred_messages!()
        end
        
        # Check for potential issues and warn appropriately
        if length(psms) < minimum_required_psms
            @warn_parameter_tuning "Insufficient PSMs for reliable parameter tuning" file_idx=file_idx psm_count=length(psms) required=minimum_required_psms
        end
        
        @debug_l1 "File processing completed" file_idx=file_idx total_psms=length(psms)
        results.psms[file_idx] = psms
        
        @trace_exit process_file! length(psms)
        
    catch e
        if is_progress_active()
            set_progress_active!(false)
            flush_deferred_messages!()
        end
        
        @warn_file_processing "file_$file_idx" "Failed to process file during parameter tuning: $e"
        @error "Failed to process file" file_idx=file_idx error=e
        rethrow(e)
    finally
        # Clean up context
        delete!(task_local_storage(), :current_file)
        delete!(task_local_storage(), :current_method)
    end
end
```

### 3. Multiple Session Usage

```julia
# Running multiple SearchDIA analyses with different verbosity
function run_comparative_analysis()
    # Run 1: Silent mode for production
    @user_info "Starting production run"
    run_with_logging("logs/production"; console_level=PIONEER_USER) do
        SearchDIA("configs/production_params.json")
    end
    
    # Run 2: Debug mode for troubleshooting
    @user_info "Starting debug run"
    run_with_logging("logs/debug"; console_level=PIONEER_DEBUG1) do
        SearchDIA("configs/debug_params.json")
    end
    
    # Run 3: Full trace for detailed analysis
    @user_info "Starting trace run"
    run_with_logging("logs/trace"; console_level=PIONEER_TRACE, show_progress=false) do
        SearchDIA("configs/trace_params.json")
    end
    
    @user_info "All analyses completed"
end

# Interactive session with level changes
function interactive_development()
    setup_pioneer_logging!("logs/dev"; console_level=PIONEER_USER)
    
    # First run with user-level output
    @user_info "Running initial analysis"
    SearchDIA("test_params.json")
    
    # Switch to debug for detailed investigation
    set_console_level!(PIONEER_DEBUG1)
    @user_info "Running with debug output"
    SearchDIA("debug_params.json")
    
    # Switch to trace for deep dive
    set_console_level!(PIONEER_TRACE)
    @user_info "Running with full trace"
    SearchDIA("trace_params.json")
    
    cleanup_logging_system!()
end
```

### 4. Console Output Examples

```bash
# Normal user mode (no warnings)
$ julia SearchDIA.jl params.json
16:30:00 Starting SearchDIA analysis
16:30:01 Loading spectral library: /path/to/lib.arrow
16:30:05 Processing MS files (3 files)
[████████████████████████████████████████] 100%
16:30:45 Parameter tuning completed
16:31:20 First pass search completed
16:32:10 Scoring search completed  
16:32:15 SearchDIA completed successfully
16:32:15 ✓ No warnings detected during execution.

# Normal user mode (with some warnings)
$ julia SearchDIA.jl params.json
16:30:00 Starting SearchDIA analysis
16:30:01 Loading spectral library: /path/to/lib.arrow
16:30:05 Processing MS files (3 files)
[████████████████████████████████████████] 100%
16:30:45 Parameter tuning completed
16:31:20 First pass search completed
16:32:10 Scoring search completed  
16:32:15 SearchDIA completed successfully
16:32:15 ⚠ 8 warnings detected during execution.
16:32:15 Warning Summary:
16:32:15   • [3×] Low quality spectra detected in file processing
16:32:15   • [2×] Parameter optimization failed to converge for some files
16:32:15   • [2×] Insufficient PSMs for reliable parameter tuning
16:32:15   • ... and 1 more warning types
16:32:15 📋 Complete warning report saved to: logs/warnings_a1b2c3d4.txt

# Debug level 1
$ julia SearchDIA.jl params.json --console-verbosity=debug1
[16:30:00] USER Starting SearchDIA analysis
[16:30:00] DEBUG1 Command line arguments: ["params.json", "--console-verbosity=debug1"]
[16:30:01] USER Loading spectral library: /path/to/lib.arrow
[16:30:01] DEBUG1 Library loaded successfully: 50000 precursors, 500000 fragments
[16:30:05] USER Processing MS files (3 files)
[16:30:05] DEBUG1 Processing file 1/3: file1.arrow
[16:30:15] DEBUG1 File 1 processed: 1250 PSMs found
[16:30:15] DEBUG1 Processing file 2/3: file2.arrow
[16:30:25] DEBUG1 File 2 processed: 1180 PSMs found
[16:30:25] DEBUG1 Processing file 3/3: file3.arrow
[16:30:35] DEBUG1 File 3 processed: 1320 PSMs found
[16:30:45] USER Parameter tuning completed
[16:30:45] DEBUG1 Optimized parameters: mass_tol=5.2ppm, rt_tol=0.8min
...

# Full trace mode
$ julia SearchDIA.jl params.json --console-verbosity=trace --no-progress
[16:30:00] USER Starting SearchDIA analysis
[16:30:00] DEBUG1 Command line arguments: ["params.json", "--console-verbosity=trace", "--no-progress"]
[16:30:00] TRACE → SearchDIA(params_path="params.json")
[16:30:00] TRACE → load_search_parameters
[16:30:00] DEBUG2 Parsing JSON parameters from params.json
[16:30:00] TRACE ← load_search_parameters: PioneerParameters(...)
[16:30:01] USER Loading spectral library: /path/to/lib.arrow
[16:30:01] TRACE → load_spectral_library
[16:30:01] DEBUG2 Reading Arrow file: /path/to/lib.arrow
[16:30:01] DEBUG2 Loaded library tables: precursors(50000 rows), fragments(500000 rows)
[16:30:01] TRACE ← load_spectral_library: SpectralLibrary(...)
...
```

<function_calls>
<invoke name="TodoWrite">
<parameter name="todos">[{"id": "1", "content": "Complete detailed Option A implementation plan", "status": "completed"}, {"id": "2", "content": "Design complete module structure and API", "status": "completed"}, {"id": "3", "content": "Write comprehensive implementation guide", "status": "completed"}, {"id": "4", "content": "Create migration strategy and testing plan", "status": "completed"}]