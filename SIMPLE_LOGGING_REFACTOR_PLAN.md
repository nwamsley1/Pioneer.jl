# Simple Logging System Refactoring Plan for Pioneer.jl

## Executive Summary

LoggingExtras.jl components (TeeLogger, TransformerLogger, etc.) appear to conflict with Arrow.Table's internal @sync/@async operations, causing the application to hang. With 150+ Arrow.Table calls across 42 files, we need a simple, robust logging solution that works everywhere.

The original `dual_println` pattern works because it's just simple, synchronous I/O without any complex Julia logging infrastructure. We'll create a macro-based system that mimics this simplicity while providing the features we need.

## Core Design Principles

1. **No LoggingExtras dependency** - Avoid TeeLogger, TransformerLogger, etc.
2. **Simple synchronous I/O** - Just like dual_println
3. **Macro-based interface** - For consistency and ease of use
4. **Thread-safe but simple** - Use basic locks only when necessary
5. **Global state management** - Simple module-level variables

## Proposed Architecture

### 1. Global State Module

```julia
module SimpleLogging

# Global state - simple mutable refs
const LOG_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const LOG_PATH = Ref{String}("")
const CONSOLE_LEVEL = Ref{Symbol}(:normal)
const FILE_LEVEL = Ref{Symbol}(:debug)
const WARNINGS = Ref{Vector{String}}(String[])

# Lock for thread-safe file writes (only if needed)
const FILE_LOCK = ReentrantLock()

end
```

### 2. Initialization Function

```julia
function init_logging(output_dir::String; 
                     console_level::Symbol = :normal,
                     file_level::Symbol = :debug)
    # Create output directory
    mkpath(output_dir)
    
    # Set log file path
    LOG_PATH[] = joinpath(output_dir, "pioneer_search_log.txt")
    
    # Open log file
    LOG_FILE[] = open(LOG_PATH[], "w")
    
    # Set levels
    CONSOLE_LEVEL[] = console_level
    FILE_LEVEL[] = file_level
    
    # Write header
    write_both("="^90)
    write_both("Pioneer Search Log - $(Dates.now())")
    write_both("="^90)
end
```

### 3. Core Writing Functions (Like dual_println)

```julia
# Direct write to both console and file - no Julia logging infrastructure
function write_both(msg::String)
    println(stdout, msg)
    if LOG_FILE[] !== nothing
        println(LOG_FILE[], msg)
        # Note: No flush! Let OS handle buffering
    end
end

function write_console(msg::String)
    println(stdout, msg)
end

function write_file(msg::String)
    if LOG_FILE[] !== nothing
        println(LOG_FILE[], msg)
    end
end
```

### 4. Level-Aware Macros

```julia
# Define log levels as simple integers
const LOG_LEVELS = Dict(
    :error => 40,
    :warn => 30,
    :info => 20,
    :debug => 10,
    :trace => 0
)

const LEVEL_COLORS = Dict(
    :error => :red,
    :warn => :yellow,
    :info => :cyan,
    :debug => :light_black,
    :trace => :light_black
)

const LEVEL_PREFIXES = Dict(
    :error => "ERROR",
    :warn => "Warning",
    :info => "Info",
    :debug => "Debug",
    :trace => "Trace"
)

# Check if message should be shown
function should_log(level::Symbol, target::Symbol)
    return LOG_LEVELS[level] >= LOG_LEVELS[target]
end

# Format message with color and prefix
function format_message(level::Symbol, msg::String; 
                        with_timestamp::Bool = true,
                        with_prefix::Bool = true)
    parts = String[]
    
    if with_timestamp
        push!(parts, "[$(Dates.format(now(), "HH:MM:SS"))]")
    end
    
    if with_prefix && haskey(LEVEL_PREFIXES, level)
        prefix = LEVEL_PREFIXES[level]
        color = get(LEVEL_COLORS, level, :normal)
        
        # Console gets color, file gets plain text
        console_prefix = "[ $(prefix)"
        push!(parts, console_prefix)
    end
    
    push!(parts, msg)
    return join(parts, " ")
end
```

### 5. User-Facing Macros

```julia
# Main logging macro
macro log(level, msg)
    quote
        local lvl = $(esc(level))
        local message = string($(esc(msg)))
        
        # Console output
        if should_log(lvl, CONSOLE_LEVEL[])
            formatted = format_message(lvl, message)
            if lvl in [:error, :warn]
                # Use color for warnings/errors
                printstyled(stdout, formatted, "\n"; 
                           color = LEVEL_COLORS[lvl])
            else
                println(stdout, formatted)
            end
        end
        
        # File output (always unformatted)
        if LOG_FILE[] !== nothing && should_log(lvl, FILE_LEVEL[])
            println(LOG_FILE[], "[$(Dates.now())] [$(lvl)] $(message)")
        end
    end
end

# Convenience macros
macro user_info(msg)
    :(@log :info $(esc(msg)))
end

macro user_warn(msg, category=nothing)
    quote
        local message = string($(esc(msg)))
        local cat = $(esc(category))
        
        # Add to warnings collection
        push!(WARNINGS[], message)
        
        # Format with category if provided
        if cat !== nothing
            message = "$(message) [category: $(cat)]"
        end
        
        @log :warn message
    end
end

macro user_error(msg)
    :(@log :error $(esc(msg)))
end

macro debug_l1(msg)
    :(@log :debug $(esc(msg)))
end

# Special macro for decorative output (no prefix)
macro user_print(msg)
    quote
        local message = string($(esc(msg)))
        println(stdout, message)
        if LOG_FILE[] !== nothing
            println(LOG_FILE[], message)
        end
    end
end
```

### 6. Cleanup Function

```julia
function close_logging()
    if LOG_FILE[] !== nothing
        # Write warnings summary
        if !isempty(WARNINGS[])
            write_both("\n" * "="^60)
            write_both("WARNINGS SUMMARY")
            write_both("Total warnings: $(length(WARNINGS[]))")
            write_both("="^60)
        end
        
        # Close file
        close(LOG_FILE[])
        LOG_FILE[] = nothing
    end
end
```

## Implementation Strategy

### Phase 1: Create New Simple Logging Module

1. **Create SimpleLogging.jl** in src/Routines/SearchDIA/
2. **No dependencies** on LoggingExtras or complex Julia logging
3. **Pure Julia I/O** operations only

### Phase 2: Replace Existing System

1. **Remove LoggingSystem folder** entirely
2. **Update importScripts.jl** to load SimpleLogging.jl
3. **Update SearchDIA.jl** to use new init/close functions
4. **Remove all LoggingExtras imports** from Pioneer.jl

### Phase 3: Update All Files

1. **Keep existing macro calls** - they'll work with new implementation
2. **Test with Arrow.Table operations** throughout codebase
3. **Verify no hanging** at any point

## Key Differences from Current System

| Current System | New Simple System |
|---------------|------------------|
| LoggingExtras (TeeLogger, etc.) | No external dependencies |
| Complex logger hierarchy | Single global file handle |
| Async operations possible | Pure synchronous I/O |
| Julia logging infrastructure | Direct println calls |
| Multiple logger objects | Simple module variables |
| Complex initialization | Simple open/close |

## Example Usage

```julia
# In SearchDIA.jl
using SimpleLogging

function SearchDIA(params_path::String)
    # Initialize logging
    SimpleLogging.init_logging(params.paths[:results])
    
    try
        @user_info "Starting SearchDIA"
        @user_info "Loading Spectral Library..."
        
        # Arrow.Table operations work fine
        SPEC_LIB = loadSpectralLibrary(SPEC_LIB_DIR, params)
        
        @debug_l1 "Library loaded successfully"
        
    finally
        SimpleLogging.close_logging()
    end
end
```

## Testing Plan

1. **Unit Tests**
   - Test each macro independently
   - Verify level filtering works
   - Check file/console output separation

2. **Integration Tests**
   - Run with ecoli_test dataset
   - Verify no hanging with Arrow.Table
   - Check log file contents

3. **Stress Tests**
   - Run with large datasets
   - Multiple threads writing simultaneously
   - Many Arrow.Table operations

## Migration Checklist

- [ ] Create SimpleLogging.jl module
- [ ] Implement core functions (init, write, close)
- [ ] Implement logging macros
- [ ] Remove LoggingSystem folder
- [ ] Update importScripts.jl
- [ ] Update Pioneer.jl (remove LoggingExtras)
- [ ] Update SearchDIA.jl initialization
- [ ] Test with ecoli dataset
- [ ] Test with full dataset
- [ ] Update documentation

## Risk Mitigation

1. **Backup current system** - Keep LoggingSystem folder initially
2. **Gradual testing** - Test on small datasets first
3. **Fallback option** - Can revert to println-only if issues
4. **Minimal changes** - Macros stay the same, only backend changes

## Conclusion

This simple logging system:
- **Eliminates Arrow.Table conflicts** by avoiding LoggingExtras
- **Maintains all functionality** users expect
- **Simplifies the codebase** significantly
- **Matches the proven dual_println pattern**
- **Works everywhere** in the 42 files with Arrow operations

The key insight is that the original dual_println worked not because it was opening/closing files, but because it was using simple, direct I/O without any complex logging infrastructure. By returning to this simplicity while adding macro convenience, we get the best of both worlds.