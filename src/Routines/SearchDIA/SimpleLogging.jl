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
                format_message(level, msg; with_timestamp=false)
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