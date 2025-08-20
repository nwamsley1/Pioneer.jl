# Pioneer.jl Logging Refactor Plan

## Current State Analysis

### Existing Logging Mechanisms
1. **dual_println/dual_print**: Local functions that write to both console and log file
2. **@info macros**: Julia's built-in info logging (currently goes to stdout)
3. **@warn macros**: Julia's built-in warning logging (currently goes to stderr) 
4. **@error macros**: Julia's built-in error logging (currently goes to stderr)
5. **println**: Direct console output
6. **ProgressBars**: Visual progress indicators for console

### Current Configuration (src/Pioneer.jl:52-59)
```julia
stdout_logger = ConsoleLogger(stdout)
stderr_logger = MinLevelLogger(ConsoleLogger(stderr), Logging.Warn)
split_logger = TeeLogger(stdout_logger, stderr_logger)
global_logger(split_logger)
```

### Problems with Current System
- No centralized logging control
- `dual_println` functions are locally defined in each main routine
- No hierarchical log levels beyond warn/error split
- Limited debugging capabilities
- No structured logging for complex data
- No log rotation or file management

## Proposed Solutions

### Option 1: Julia Logging System Extension (Recommended)

#### Architecture Overview
Create a centralized logging system that extends Julia's built-in `Logging` module with custom loggers and formatters.

```julia
# New logging levels
const LOG_PROGRESS = LogLevel(-500)    # For progress indicators
const LOG_USER = LogLevel(0)           # User-facing messages  
const LOG_DEBUG_L1 = LogLevel(100)     # Basic debugging
const LOG_DEBUG_L2 = LogLevel(200)     # Detailed debugging
const LOG_TRACE = LogLevel(300)        # Function tracing
```

#### Multi-Destination Logging
```julia
struct PioneerLoggingConfig
    console_level::LogLevel        # What goes to console
    simplified_log_level::LogLevel # What goes to simplified log
    full_log_level::LogLevel       # What goes to full debug log
    show_progress::Bool            # Show progress bars on console
    show_warnings::Bool            # Show warnings on console
end

# Three output destinations:
# 1. Console: progress bars, user messages, warnings only
# 2. Simplified log: dual_println equivalent messages  
# 3. Full debug log: everything including @info, @debug levels
```

#### Implementation Strategy
1. **PioneerLogger**: Custom logger that routes messages based on level
2. **Structured macros**: `@user_info`, `@debug_l1`, `@debug_l2`, `@trace`
3. **Progress integration**: Seamless ProgressBars with logging
4. **File management**: Automatic log rotation and cleanup

#### Usage Examples
```julia
@user_info "Starting SearchDIA analysis"              # → Console + simplified log
@debug_l1 "Processing file $(filename)"               # → Full log only  
@debug_l2 "Fragment match details: $(match_data)"     # → Full log only
@warn "Low quality spectra detected"                  # → Console + both logs
@trace "Entering function matchPeaks"                 # → Full log only
```

### Option 2: Structured Logging with JSON Output

#### Architecture Overview
Implement structured logging that outputs JSON-formatted log entries for machine parsing while maintaining human-readable console output.

```julia
struct StructuredLogger
    console_formatter::Function    # Human-readable formatting
    file_formatter::Function      # JSON formatting
    destinations::Vector{IO}      # Multiple output streams
end
```

#### Benefits
- Machine-parseable logs for analysis
- Rich metadata capture (timestamps, thread IDs, function names)
- Integration with log aggregation systems
- Structured search and filtering capabilities

#### Implementation
```julia
@log_structured :user_info "Analysis started" file_count=length(files) timestamp=now()
# Console: "2025-01-20 14:30:00 [INFO] Analysis started (5 files)"
# File: {"level":"INFO","message":"Analysis started","file_count":5,"timestamp":"2025-01-20T14:30:00"}
```

### Option 3: Custom Logging Framework

#### Architecture Overview
Build a Pioneer-specific logging framework from scratch with complete control over formatting, routing, and behavior.

```julia
abstract type PioneerLogLevel end
struct UserLevel <: PioneerLogLevel end
struct ProgressLevel <: PioneerLogLevel end  
struct DebugLevel <: PioneerLogLevel 
    depth::Int  # 1, 2, 3 for increasing detail
end
struct TraceLevel <: PioneerLogLevel end

struct PioneerLogger
    outputs::Dict{Type{<:PioneerLogLevel}, Vector{IO}}
    formatters::Dict{Type{<:PioneerLogLevel}, Function}
    filters::Dict{Type{<:PioneerLogLevel}, Function}
end
```

#### Benefits
- Complete control over behavior
- No dependency on Julia's logging system
- Custom filtering and formatting
- Performance optimization for high-throughput logging

#### Drawbacks
- More complex to implement and maintain
- Less integration with Julia ecosystem tools
- Potential performance overhead

### Option 4: Hierarchical Logger Composition

#### Architecture Overview
Create a composition of multiple specialized loggers that can be configured independently.

```julia
struct PioneerLoggingSystem
    console_logger::ConsoleLogger
    simplified_logger::FileLogger
    full_logger::FileLogger
    progress_manager::ProgressManager
    router::LogRouter
end

function route_message(system::PioneerLoggingSystem, level, message, metadata...)
    # Intelligent routing based on level, content, and configuration
end
```

#### Benefits
- Modular design allows independent configuration of each logger
- Easy to extend with new output destinations
- Clear separation of concerns
- Testable components

## Detailed Implementation Recommendation: Option 1

### Phase 1: Core Infrastructure

#### 1.1 Create Logging Module
```julia
# src/utils/Logging/PioneerLogging.jl
module PioneerLogging
    using Logging, LoggingExtras, Dates
    
    # Custom log levels
    const LOG_PROGRESS = LogLevel(-500)
    const LOG_USER = LogLevel(0)
    const LOG_DEBUG_L1 = LogLevel(100)
    const LOG_DEBUG_L2 = LogLevel(200)
    const LOG_TRACE = LogLevel(300)
    
    # Configuration
    struct LoggingConfig
        console_level::LogLevel
        simplified_log_level::LogLevel
        full_log_level::LogLevel
        show_progress::Bool
        show_warnings::Bool
        log_dir::String
        max_log_size::Int
        max_log_files::Int
    end
    
    # Main logger implementation
    struct PioneerLogger <: AbstractLogger
        config::LoggingConfig
        console_io::IO
        simplified_io::IO
        full_io::IO
        progress_enabled::Bool
    end
end
```

#### 1.2 Logger Implementation
```julia
function handle_message(logger::PioneerLogger, level, message, _module, group, id, file, line; kwargs...)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Route to appropriate destinations based on level
    if level <= logger.config.console_level
        write_console_message(logger, level, message, timestamp, kwargs...)
    end
    
    if level <= logger.config.simplified_log_level
        write_simplified_log(logger, level, message, timestamp, kwargs...)
    end
    
    if level <= logger.config.full_log_level
        write_full_log(logger, level, message, _module, file, line, timestamp, kwargs...)
    end
end
```

#### 1.3 Custom Macros
```julia
macro user_info(msg, args...)
    quote
        @logmsg LOG_USER $(esc(msg)) $(map(esc, args)...)
    end
end

macro debug_l1(msg, args...)
    quote
        @logmsg LOG_DEBUG_L1 $(esc(msg)) $(map(esc, args)...)
    end
end

macro debug_l2(msg, args...)
    quote
        @logmsg LOG_DEBUG_L2 $(esc(msg)) $(map(esc, args)...)
    end
end

macro trace(func_name)
    quote
        @logmsg LOG_TRACE "Entering function: $($(string(func_name)))"
    end
end
```

### Phase 2: Progress Bar Integration

#### 2.1 Progress-Aware Logger
```julia
struct ProgressAwareLogger <: AbstractLogger
    base_logger::PioneerLogger
    active_progress_bars::Set{String}
    progress_channel::Channel{ProgressUpdate}
end

# Coordinate with ProgressBars.jl to avoid output conflicts
function with_progress_logging(f, logger, description)
    # Temporarily adjust console output behavior
    # Execute function with progress bar
    # Restore normal logging
end
```

#### 2.2 Smart Progress Integration
```julia
# Automatically detect when to show progress vs log messages
function smart_iterate(iterable, description; show_progress=true)
    if show_progress && should_show_progress(length(iterable))
        return ProgressBar(iterable, description=description)
    else
        @user_info "Processing $description ($(length(iterable)) items)"
        return iterable
    end
end
```

### Phase 3: Migration Strategy

#### 3.1 Backward Compatibility Layer
```julia
# Provide dual_println equivalent during migration
function dual_println(args...)
    @user_info join(args, " ")
end

function dual_print(args...)
    @user_info join(args, "")
end
```

#### 3.2 Gradual Migration
1. **Week 1-2**: Implement core logging infrastructure
2. **Week 3-4**: Add progress bar integration
3. **Week 5-6**: Replace dual_println calls with @user_info
4. **Week 7-8**: Add @debug_l1/@debug_l2 calls to replace existing @info
5. **Week 9-10**: Add @trace calls for detailed debugging
6. **Week 11-12**: Testing and refinement

### Phase 4: Advanced Features

#### 4.1 Configuration Management
```julia
# Default configuration
const DEFAULT_CONFIG = LoggingConfig(
    console_level = LOG_USER,           # Only user messages and warnings on console
    simplified_log_level = LOG_USER,    # Only user messages in simplified log
    full_log_level = LOG_TRACE,         # Everything in full debug log
    show_progress = true,               # Show progress bars
    show_warnings = true,               # Show warnings on console
    log_dir = "logs",                   # Log directory
    max_log_size = 50_000_000,         # 50MB max per log file
    max_log_files = 5                   # Keep 5 log files
)

# Configuration from JSON/environment
function load_logging_config(config_path::String)::LoggingConfig
    # Load from file or use defaults
end
```

#### 4.2 Log File Management
```julia
function setup_log_files(config::LoggingConfig)
    ensure_log_directory(config.log_dir)
    
    simplified_path = rotate_log_file(
        joinpath(config.log_dir, "pioneer_simple.log"),
        config.max_log_size,
        config.max_log_files
    )
    
    full_path = rotate_log_file(
        joinpath(config.log_dir, "pioneer_full.log"),
        config.max_log_size,
        config.max_log_files
    )
    
    return open(simplified_path, "a"), open(full_path, "a")
end
```

#### 4.3 Performance Considerations
```julia
# Async logging for high-throughput scenarios
struct AsyncPioneerLogger <: AbstractLogger
    base_logger::PioneerLogger
    message_queue::Channel{LogMessage}
    worker_task::Task
end

function start_async_logging(config::LoggingConfig)
    # Start background task for log writing
    # Use thread-safe queue for message passing
end
```

## Alternative Approaches Considered

### Approach A: Extend LoggingExtras.jl
- **Pros**: Leverages existing robust library
- **Cons**: May not provide enough customization for Pioneer's needs

### Approach B: Use Memento.jl
- **Pros**: Mature logging library with hierarchical loggers
- **Cons**: Additional dependency, learning curve

### Approach C: Simple File-Based Approach
- **Pros**: Minimal complexity
- **Cons**: Limited features, no structured logging

## Testing Strategy

### Unit Tests
```julia
@testset "PioneerLogging Tests" begin
    @testset "Message Routing" begin
        # Test that messages go to correct destinations
    end
    
    @testset "Level Filtering" begin
        # Test that only appropriate levels are logged
    end
    
    @testset "Progress Integration" begin
        # Test progress bar coordination
    end
    
    @testset "File Management" begin
        # Test log rotation and cleanup
    end
end
```

### Integration Tests
```julia
@testset "SearchDIA Logging Integration" begin
    config = LoggingConfig(...)
    with_test_logging(config) do
        SearchDIA("test_params.json")
        # Verify correct messages in each log file
    end
end
```

## Migration Checklist

- [ ] Implement core PioneerLogging module
- [ ] Add configuration management
- [ ] Implement progress bar integration
- [ ] Create backward compatibility layer
- [ ] Replace dual_println calls
- [ ] Add debug-level logging throughout codebase
- [ ] Add trace logging for key functions
- [ ] Implement log file rotation
- [ ] Add comprehensive tests
- [ ] Update documentation
- [ ] Performance optimization
- [ ] User acceptance testing

## Configuration Examples

### Development Configuration
```json
{
    "logging": {
        "console_level": "USER",
        "simplified_log_level": "USER", 
        "full_log_level": "DEBUG_L2",
        "show_progress": true,
        "show_warnings": true,
        "log_dir": "logs",
        "max_log_size": 10000000,
        "max_log_files": 3
    }
}
```

### Production Configuration
```json
{
    "logging": {
        "console_level": "USER",
        "simplified_log_level": "USER",
        "full_log_level": "DEBUG_L1", 
        "show_progress": true,
        "show_warnings": true,
        "log_dir": "/var/log/pioneer",
        "max_log_size": 100000000,
        "max_log_files": 10
    }
}
```

### Debug Configuration
```json
{
    "logging": {
        "console_level": "DEBUG_L1",
        "simplified_log_level": "DEBUG_L1",
        "full_log_level": "TRACE",
        "show_progress": false,
        "show_warnings": true,
        "log_dir": "debug_logs",
        "max_log_size": 50000000,
        "max_log_files": 20
    }
}
```

## Expected Benefits

1. **Improved Debugging**: Multi-level debug logging enables better troubleshooting
2. **Better User Experience**: Clean console output with progress indicators
3. **Operational Monitoring**: Structured logs enable monitoring and alerting
4. **Development Efficiency**: Trace logging helps understand execution flow
5. **Configuration Flexibility**: Different logging levels for different environments
6. **Performance**: Async logging reduces impact on main computation
7. **Maintenance**: Automatic log rotation prevents disk space issues

## Implementation Timeline

- **Phase 1** (2 weeks): Core infrastructure and basic functionality
- **Phase 2** (2 weeks): Progress bar integration and testing  
- **Phase 3** (3 weeks): Migration of existing code
- **Phase 4** (2 weeks): Advanced features and optimization
- **Phase 5** (1 week): Documentation and user testing

**Total Estimated Time**: 10 weeks

This plan provides a robust, extensible logging system that addresses all the identified requirements while maintaining backward compatibility and providing a clear migration path.