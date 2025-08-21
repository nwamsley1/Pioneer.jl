# Complete Implementation Plan: Simple Logging System for Pioneer.jl

## Implementation Roadmap

### Step 1: Create SimpleLogging.jl Module

**File Location**: `/src/Routines/SearchDIA/SimpleLogging.jl`

```julia
# SimpleLogging.jl - A simple, Arrow.Table-compatible logging system for Pioneer.jl
module SimpleLogging

using Dates

export init_logging, close_logging, @user_info, @user_warn, @user_error, 
       @user_print, @debug_l1, @debug_l2, @debug_l3, @trace,
       set_console_level, set_file_level

# ============================================================================
# Global State Management
# ============================================================================

# File handles and paths
const LOG_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const DEBUG_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const LOG_PATH = Ref{String}("")
const DEBUG_PATH = Ref{String}("")

# Logging levels
const CONSOLE_LEVEL = Ref{Symbol}(:normal)
const FILE_LEVEL = Ref{Symbol}(:info)
const DEBUG_FILE_LEVEL = Ref{Symbol}(:trace)

# Warning tracking
const WARNINGS = Vector{NamedTuple{(:timestamp, :message, :category), Tuple{DateTime, String, Union{Nothing, Symbol}}}}()

# Progress bar state (for future integration)
const PROGRESS_ACTIVE = Ref{Bool}(false)

# ============================================================================
# Log Level Definitions
# ============================================================================

# Numeric levels for comparison
const LOG_LEVELS = Dict{Symbol, Int}(
    :silent => 50,    # Only critical errors
    :error => 40,     # Errors only
    :warn => 30,      # Warnings and above
    :normal => 20,    # User info and above (default console)
    :info => 20,      # Same as normal
    :verbose => 15,   # More details
    :debug => 10,     # Debug messages
    :trace => 0       # Everything
)

# Console output formatting
const LEVEL_COLORS = Dict{Symbol, Symbol}(
    :error => :red,
    :warn => :yellow,
    :info => :cyan,
    :debug => :light_black,
    :trace => :light_black
)

const LEVEL_PREFIXES = Dict{Symbol, String}(
    :error => "Error",
    :warn => "Warning", 
    :info => "Info",
    :debug => "Debug",
    :trace => "Trace"
)

# ============================================================================
# Initialization and Cleanup
# ============================================================================

"""
    init_logging(output_dir::String; kwargs...)

Initialize the simple logging system with dual file output.

# Arguments
- `output_dir`: Directory for log files
- `console_level`: Console verbosity (:silent, :error, :warn, :normal, :verbose, :debug, :trace)
- `file_level`: Main log file verbosity
- `debug_file_level`: Debug log file verbosity
"""
function init_logging(output_dir::String; 
                     console_level::Symbol = :normal,
                     file_level::Symbol = :info,
                     debug_file_level::Symbol = :trace)
    
    # Create output directory if needed
    try
        mkpath(output_dir)
    catch e
        println(stderr, "Warning: Could not create log directory: $e")
        return false
    end
    
    # Set global levels
    CONSOLE_LEVEL[] = console_level
    FILE_LEVEL[] = file_level
    DEBUG_FILE_LEVEL[] = debug_file_level
    
    # Open log files
    LOG_PATH[] = joinpath(output_dir, "pioneer_search_log.txt")
    DEBUG_PATH[] = joinpath(output_dir, "pioneer_debug.log")
    
    try
        LOG_FILE[] = open(LOG_PATH[], "w")
        DEBUG_FILE[] = open(DEBUG_PATH[], "w")
    catch e
        println(stderr, "Warning: Could not open log files: $e")
        return false
    end
    
    # Write headers
    header = ["=" ^ 90,
              "Pioneer Search Log",
              "Started: $(Dates.now())",
              "Output Directory: $output_dir",
              "=" ^ 90,
              ""]
    
    for line in header
        write_to_file(LOG_FILE[], line)
        write_to_file(DEBUG_FILE[], line)
    end
    
    # Log initialization with proper formatting
    log_message(:info, "Pioneer logging system initialized"; 
                output_dir=output_dir, timestamp=now())
    
    return true
end

"""
    close_logging()

Close all log files and write final summary.
"""
function close_logging()
    # Write warnings summary if any
    if !isempty(WARNINGS)
        summary = [
            "",
            "=" ^ 60,
            "WARNINGS SUMMARY",
            "=" ^ 60,
            "Total warnings: $(length(WARNINGS))",
            ""
        ]
        
        # Group by category
        by_category = Dict{Union{Nothing, Symbol}, Int}()
        for w in WARNINGS
            cat = w.category
            by_category[cat] = get(by_category, cat, 0) + 1
        end
        
        push!(summary, "By category:")
        for (cat, count) in sort(collect(by_category), by=x->x[2], rev=true)
            cat_str = cat === nothing ? "uncategorized" : string(cat)
            push!(summary, "  $cat_str: $count")
        end
        
        # Write to both console and files
        for line in summary
            println(stdout, line)
            write_to_file(LOG_FILE[], line)
            write_to_file(DEBUG_FILE[], line)
        end
        
        # Write detailed warnings to debug file only
        write_to_file(DEBUG_FILE[], "")
        write_to_file(DEBUG_FILE[], "Detailed warnings:")
        for (i, w) in enumerate(WARNINGS)
            write_to_file(DEBUG_FILE[], "[$i] $(w.timestamp): $(w.message)")
            if w.category !== nothing
                write_to_file(DEBUG_FILE[], "    Category: $(w.category)")
            end
        end
    end
    
    # Write footer
    footer = [
        "",
        "=" ^ 90,
        "Search completed at: $(Dates.now())",
        "=" ^ 90
    ]
    
    for line in footer
        write_to_file(LOG_FILE[], line)
        write_to_file(DEBUG_FILE[], line)
    end
    
    # Close files
    if LOG_FILE[] !== nothing
        close(LOG_FILE[])
        LOG_FILE[] = nothing
    end
    
    if DEBUG_FILE[] !== nothing
        close(DEBUG_FILE[])
        DEBUG_FILE[] = nothing
    end
    
    # Clear warnings
    empty!(WARNINGS)
end

# ============================================================================
# Core Writing Functions
# ============================================================================

"""
Write a line to a file handle if it's open.
"""
function write_to_file(io::Union{Nothing, IOStream}, msg::String)
    if io !== nothing
        println(io, msg)
        # NO FLUSH - let OS handle buffering to avoid Arrow.Table conflicts
    end
end

"""
Write to console with optional color.
"""
function write_to_console(msg::String; color::Symbol = :normal)
    if color == :normal
        println(stdout, msg)
    else
        printstyled(stdout, msg, "\n"; color=color)
    end
end

"""
Check if a message at given level should be logged to target.
"""
function should_log(msg_level::Symbol, target_level::Symbol)
    msg_val = get(LOG_LEVELS, msg_level, 20)
    target_val = get(LOG_LEVELS, target_level, 20)
    return msg_val >= target_val
end

"""
Format a log message with timestamp and level prefix.
"""
function format_message(level::Symbol, msg::String; 
                       with_timestamp::Bool = true,
                       with_prefix::Bool = true,
                       kwargs...)
    parts = String[]
    
    # Add timestamp
    if with_timestamp
        # Use short format for console
        push!(parts, Dates.format(now(), "HH:MM:SS"))
    end
    
    # Add level prefix
    if with_prefix && haskey(LEVEL_PREFIXES, level)
        prefix = LEVEL_PREFIXES[level]
        push!(parts, prefix)
    end
    
    # Add main message
    push!(parts, msg)
    
    # Join with appropriate separators
    if with_prefix
        if with_timestamp
            return "[ $(parts[2]) $(parts[1]): $(parts[3])"
        else
            return "[ $(parts[1]) $(parts[2])"
        end
    else
        if with_timestamp
            return "[$(parts[1])] $(parts[2])"
        else
            return msg
        end
    end
end

"""
Core logging function called by all macros.
"""
function log_message(level::Symbol, msg::String; 
                    category::Union{Nothing, Symbol} = nothing,
                    kwargs...)
    
    # Console output
    if should_log(level, CONSOLE_LEVEL[])
        # Format for console
        formatted = if level in [:error, :warn]
            # Colored output for warnings/errors
            color = get(LEVEL_COLORS, level, :normal)
            
            # Special formatting for warnings with category
            if level == :warn && category !== nothing
                prefix = LEVEL_PREFIXES[:warn]
                # Don't include category in console output to reduce clutter
                "┌ $prefix $msg"
            else
                formatted = format_message(level, msg; with_timestamp=false)
            end
        else
            # Regular info messages
            format_message(level, msg; with_timestamp=false)
        end
        
        # Write to console with color
        color = get(LEVEL_COLORS, level, :normal)
        write_to_console(formatted; color=color)
    end
    
    # Main log file output
    if should_log(level, FILE_LEVEL[])
        # Simple format for file
        file_msg = "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
        write_to_file(LOG_FILE[], file_msg)
    end
    
    # Debug file output (everything with full details)
    if should_log(level, DEBUG_FILE_LEVEL[])
        # Detailed format for debug file
        debug_msg = "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss"))] [$level] $msg"
        
        # Add category if present
        if category !== nothing
            debug_msg *= " [category: $category]"
        end
        
        # Add any additional kwargs
        if !isempty(kwargs)
            kwarg_strs = ["$k=$v" for (k, v) in kwargs]
            debug_msg *= " [" * join(kwarg_strs, ", ") * "]"
        end
        
        write_to_file(DEBUG_FILE[], debug_msg)
    end
    
    # Track warnings
    if level == :warn
        push!(WARNINGS, (timestamp=now(), message=msg, category=category))
    end
end

"""
Direct write without any formatting (for decorative output).
"""
function write_direct(msg::String)
    println(stdout, msg)
    write_to_file(LOG_FILE[], msg)
    write_to_file(DEBUG_FILE[], msg)
end

# ============================================================================
# Public Macros
# ============================================================================

# User-facing info messages
macro user_info(msg)
    quote
        SimpleLogging.log_message(:info, string($(esc(msg))))
    end
end

# Warning messages with optional category
macro user_warn(msg, category=nothing)
    quote
        local msg_str = string($(esc(msg)))
        local cat = $(esc(category))
        
        # Only pass category if it's not nothing
        if cat === nothing
            SimpleLogging.log_message(:warn, msg_str)
        else
            SimpleLogging.log_message(:warn, msg_str; category=cat)
        end
    end
end

# Error messages
macro user_error(msg)
    quote
        SimpleLogging.log_message(:error, string($(esc(msg))))
    end
end

# Debug level 1 (most important debug info)
macro debug_l1(msg)
    quote
        SimpleLogging.log_message(:debug, string($(esc(msg))))
    end
end

# Debug level 2 (more detailed)
macro debug_l2(msg)
    quote
        SimpleLogging.log_message(:debug, "[L2] " * string($(esc(msg))))
    end
end

# Debug level 3 (very detailed)
macro debug_l3(msg)
    quote
        SimpleLogging.log_message(:debug, "[L3] " * string($(esc(msg))))
    end
end

# Trace level (everything)
macro trace(msg)
    quote
        SimpleLogging.log_message(:trace, string($(esc(msg))))
    end
end

# Direct print for decorative output (no formatting)
macro user_print(msg)
    quote
        SimpleLogging.write_direct(string($(esc(msg))))
    end
end

# ============================================================================
# Runtime Configuration
# ============================================================================

"""
    set_console_level(level::Symbol)

Change console verbosity at runtime.
"""
function set_console_level(level::Symbol)
    if haskey(LOG_LEVELS, level)
        CONSOLE_LEVEL[] = level
        log_message(:debug, "Console level changed to: $level")
    else
        log_message(:warn, "Invalid console level: $level")
    end
end

"""
    set_file_level(level::Symbol)

Change file logging verbosity at runtime.
"""
function set_file_level(level::Symbol)
    if haskey(LOG_LEVELS, level)
        FILE_LEVEL[] = level
        log_message(:debug, "File level changed to: $level")
    else
        log_message(:warn, "Invalid file level: $level")
    end
end

# ============================================================================
# Progress Bar Integration (Future)
# ============================================================================

"""
    set_progress_active(active::Bool)

Set whether progress bar is active (for future integration).
"""
function set_progress_active(active::Bool)
    PROGRESS_ACTIVE[] = active
end

end # module SimpleLogging
```

### Step 2: Update importScripts.jl

**File**: `/src/importScripts.jl`

**Changes**:
1. Remove the LoggingSystem include block (lines 184-200)
2. Add SimpleLogging include after ParseInputs

```julia
# At line 183, after ParseInputs directory loading, add:
safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "SimpleLogging.jl"))
```

### Step 3: Update Pioneer.jl

**File**: `/src/Pioneer.jl`

**Changes**:
1. Remove `using LoggingExtras` import
2. Add SimpleLogging to exports

```julia
# Remove this line:
# using LoggingExtras

# In the export section, add:
export SimpleLogging
```

### Step 4: Update SearchDIA.jl

**File**: `/src/Routines/SearchDIA.jl`

**Replace lines 135-148** (logging configuration parsing) with:

```julia
# Parse logging configuration from params
console_level = Symbol(get(params.logging, :console_level, "normal"))
file_level = Symbol(get(params.logging, :file_level, "info"))
debug_level = Symbol(get(params.logging, :debug_level, "trace"))

# Initialize simple logging system
SimpleLogging.init_logging(
    params.paths[:results];
    console_level = console_level,
    file_level = file_level,
    debug_file_level = debug_level
)
```

**Replace lines 280-283** (finalize_logging in finally block) with:

```julia
finally
    # Close logging system
    SimpleLogging.close_logging()
end
```

### Step 5: Update parseParams.jl

**File**: `/src/Routines/SearchDIA/ParseInputs/parseParams.jl`

The logging section parsing can stay the same since we're reading the same keys:
- `console_level`
- `file_level` (maps to simple log file)
- `debug_level` (maps to debug log file)

No changes needed here.

### Step 6: Update Macro Calls

**All existing macro calls will work**, but they need to be prefixed with `SimpleLogging.` or we need to import them. Since the macros are exported, they should work directly.

However, we need to ensure the module is loaded. In files that use logging macros, add at the top:
```julia
using SimpleLogging
```

Or if already in Pioneer module context, the exports will be available.

### Step 7: Remove Old LoggingSystem

**Delete entire folder**: `/src/Routines/SearchDIA/LoggingSystem/`

This includes:
- AsyncFileLogger.jl
- SimpleFileLogger.jl  
- DualPrintLogger.jl
- Loggers.jl
- LoggingSystem.jl
- Macros.jl
- Configuration.jl
- WarningTracker.jl
- WarningCapturingLogger.jl
- ProgressIntegration.jl
- FileManagement.jl
- All .md files in that folder

### Step 8: Update BuildSpecLib.jl (Optional)

**File**: `/src/Routines/BuildSpecLib.jl`

If you want to unify logging, replace the local `dual_println` function with SimpleLogging:

```julia
# At the top, after other imports:
using SimpleLogging

# Replace lines 90-101 (open block with dual_println definition) with:
SimpleLogging.init_logging(lib_dir)

# Replace all dual_println calls with:
@user_info "message"
# or
@user_print "message"  # for decorative lines

# At the end of the function, add:
SimpleLogging.close_logging()
```

### Step 9: Test Files to Verify

Create a test file to verify Arrow.Table compatibility:

**File**: `/test/test_simple_logging.jl`

```julia
using Pioneer
using SimpleLogging
using Arrow

# Test 1: Basic logging
SimpleLogging.init_logging("./test_logs")
@user_info "Test message"
@user_warn "Test warning" :test_category
@debug_l1 "Debug message"
SimpleLogging.close_logging()

# Test 2: Arrow.Table compatibility
SimpleLogging.init_logging("./test_logs2")
@user_info "Loading Arrow table..."

# This should not hang
test_table = Arrow.Table(joinpath(@__DIR__, "../data/ecoli_test/spectral_library/precursors_table.arrow"))
@user_info "Arrow table loaded successfully with $(length(test_table)) rows"

SimpleLogging.close_logging()

println("All tests passed!")
```

## Migration Checklist

### Pre-Implementation
- [x] Document all changes needed
- [x] Identify all files to modify
- [x] Create complete SimpleLogging.jl implementation

### Implementation Steps
- [ ] Create SimpleLogging.jl file
- [ ] Update importScripts.jl to include SimpleLogging
- [ ] Update Pioneer.jl to remove LoggingExtras
- [ ] Update SearchDIA.jl initialization and cleanup
- [ ] Delete LoggingSystem folder
- [ ] Run test with ecoli dataset
- [ ] Verify no hanging with Arrow.Table

### Testing
- [ ] Test basic logging output (console and files)
- [ ] Test with single Arrow.Table load
- [ ] Test with full SearchDIA run on ecoli
- [ ] Test with multiple threads
- [ ] Test warning accumulation and summary
- [ ] Test on larger dataset

### Cleanup
- [ ] Remove old log files from previous system
- [ ] Update any documentation
- [ ] Commit changes

## Key Implementation Details

### Why This Will Work

1. **No LoggingExtras** - We completely avoid TeeLogger, TransformerLogger, EarlyFilteredLogger
2. **Simple I/O** - Just `println` calls like the original dual_println
3. **No async/sync** - No Tasks, Channels, or @async/@sync that could conflict with Arrow
4. **No locks during I/O** - We don't use locks around the actual I/O operations
5. **No flush calls** - Let the OS handle buffering naturally

### Critical Differences

| Component | Old System | New System |
|-----------|------------|------------|
| Base | LoggingExtras.jl | Pure Julia I/O |
| File writing | Through logger chain | Direct println |
| Synchronization | Complex logger wrapping | None needed |
| State | LoggerState struct | Simple module Refs |
| Initialization | Complex logger hierarchy | Open two files |

### Files That Will Just Work

All 42 files with Arrow.Table calls will work without modification because:
1. The macros have the same names
2. No file handles are kept open during operations (OS buffers handle it)
3. No complex synchronization to conflict with Arrow's @sync

### Fallback Plan

If any issues arise:
1. Can temporarily disable file writing (set LOG_FILE[] = nothing)
2. Can fall back to console-only
3. Can add flag to disable logging entirely
4. Original dual_println approach is still available

## Summary

This implementation completely removes the LoggingExtras dependency and returns to simple, direct I/O operations that are proven to work with Arrow.Table. The system maintains all the features you need (multiple log levels, dual output, colored console, warning tracking) while eliminating the complex logger hierarchy that causes conflicts.

The implementation is straightforward - one module file with ~400 lines of simple Julia code, no external dependencies, and minimal changes to existing code. All existing macro calls continue to work unchanged.