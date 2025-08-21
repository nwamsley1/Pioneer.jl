# Macros.jl
# Custom logging macros for Pioneer

# Export all custom macros
export @user_info, @user_warn, @user_error
export @debug_l1, @debug_l2, @debug_l3
export @trace
export @progress_start, @progress_update, @progress_end

# Define custom log levels
const USER_INFO_LEVEL = LogLevel(Int(Logging.Info) - 100)
const DEBUG_L1_LEVEL = Logging.Debug
const DEBUG_L2_LEVEL = LogLevel(Int(Logging.Debug) - 100)
const DEBUG_L3_LEVEL = LogLevel(Int(Logging.Debug) - 200)
const TRACE_LEVEL = LogLevel(Int(Logging.Debug) - 1000)

"""
    @user_info msg [key=value...]

Log a user-facing informational message that appears in all outputs.
"""
macro user_info(msg, args...)
    quote
        @logmsg USER_INFO_LEVEL $(esc(msg)) $(esc.(args)...) is_user_message=true _group=:user_info
    end
end

"""
    @user_warn msg [key=value...]

Log a user-facing warning that appears in all outputs and is tracked.
"""
macro user_warn(msg, args...)
    # Extract category if provided
    category_expr = nothing
    other_args = []
    
    for arg in args
        if isa(arg, Expr) && arg.head == :(=) && arg.args[1] == :category
            category_expr = arg.args[2]
        else
            push!(other_args, arg)
        end
    end
    
    quote
        @warn $(esc(msg)) $(esc.(other_args)...) is_user_message=true category=$(esc(category_expr))
    end
end

"""
    @user_error msg [key=value...]

Log a user-facing error that appears in all outputs.
"""
macro user_error(msg, args...)
    quote
        @error $(esc(msg)) $(esc.(args)...) is_user_message=true
    end
end

"""
    @debug_l1 msg [key=value...]

Debug level 1 - High-level debugging information.
"""
macro debug_l1(msg, args...)
    quote
        @logmsg DEBUG_L1_LEVEL $(esc(msg)) $(esc.(args)...) _group=:debug_l1
    end
end

"""
    @debug_l2 msg [key=value...]

Debug level 2 - Detailed debugging information.
"""
macro debug_l2(msg, args...)
    quote
        @logmsg DEBUG_L2_LEVEL $(esc(msg)) $(esc.(args)...) _group=:debug_l2
    end
end

"""
    @debug_l3 msg [key=value...]

Debug level 3 - Very detailed debugging information.
"""
macro debug_l3(msg, args...)
    quote
        @logmsg DEBUG_L3_LEVEL $(esc(msg)) $(esc.(args)...) _group=:debug_l3
    end
end

"""
    @trace msg [key=value...]

Trace level - Extremely detailed trace information.
"""
macro trace(msg, args...)
    quote
        @logmsg TRACE_LEVEL $(esc(msg)) $(esc.(args)...) _group=:trace
    end
end

# Progress-related macros

"""
    @progress_start name total

Start a progress tracking operation.
"""
macro progress_start(name, total)
    quote
        @info "Starting: $($(esc(name)))" progress_step=true total=$(esc(total))
    end
end

"""
    @progress_update name current total

Update progress for an operation.
"""
macro progress_update(name, current, total)
    quote
        @logmsg USER_INFO_LEVEL "Progress: $($(esc(name)))" progress_active=true current=$(esc(current)) total=$(esc(total))
    end
end

"""
    @progress_end name

Complete a progress tracking operation.
"""
macro progress_end(name)
    quote
        @info "Completed: $($(esc(name)))" progress_step=true
    end
end

# Helper macro for timing operations

"""
    @timed_operation name expr

Execute an expression and log its execution time.
"""
macro timed_operation(name, expr)
    quote
        local start_time = time()
        local result = $(esc(expr))
        local elapsed = time() - start_time
        @debug_l1 "Operation completed" operation=$(esc(name)) elapsed_seconds=elapsed
        result
    end
end

# Conditional logging helpers

"""
    @log_if condition level msg [args...]

Log a message only if condition is true.
"""
macro log_if(condition, level, msg, args...)
    quote
        if $(esc(condition))
            @logmsg $(esc(level)) $(esc(msg)) $(esc.(args)...)
        end
    end
end

"""
    @warn_once msg [args...]

Log a warning only once per session (deduplicated by message).
"""
const WARNED_MESSAGES = Set{String}()

macro warn_once(msg, args...)
    quote
        local msg_str = string($(esc(msg)))
        if !(msg_str in WARNED_MESSAGES)
            push!(WARNED_MESSAGES, msg_str)
            @warn $(esc(msg)) $(esc.(args)...)
        end
    end
end