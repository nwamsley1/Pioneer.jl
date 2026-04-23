const ESSENTIAL_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const CONSOLE_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const DEBUG_FILE = Ref{Union{Nothing, IOStream}}(nothing)
const WARNINGS_FILE = Ref{Union{Nothing, IOStream}}(nothing)

const DEBUG_CONSOLE_LEVEL = Ref{Int}(0)
const MAX_LOG_MSG_BYTES = Ref{Int}(4096)

const ESSENTIAL_PATTERNS = [
    r"^Starting search at:",
    r"^Output directory:",
    r"^Loading Parameters",
    r"^Loading Spectral Library",
    r"^Initializing Search Context",
    r"^Executing .+\.\.\.",
]

function truncate_for_log(msg::String; max_bytes::Int=MAX_LOG_MSG_BYTES[])
    nb = ncodeunits(msg)
    if nb <= max_bytes
        return msg
    end

    n = 0
    i = firstindex(msg)
    stop = lastindex(msg)
    last_end = i
    while i <= stop
        j = nextind(msg, i)
        cu = j - i
        if n + cu > max_bytes
            break
        end
        n += cu
        last_end = j
        i = j
    end

    if last_end <= firstindex(msg)
        return string(msg[1:prevind(msg, lastindex(msg))], " ... [truncated]")
    end

    prefix = String(@view msg[1:prevind(msg, last_end)])
    omitted = nb - n
    return string(prefix, " ... [truncated ", omitted, " bytes]")
end

function is_essential_message(msg::String)
    for pattern in ESSENTIAL_PATTERNS
        if occursin(pattern, msg)
            return true
        end
    end
    return false
end

function user_info(msg::String)
    msg_trunc = truncate_for_log(msg)

    printstyled("[ ", bold=true, color=:cyan)
    printstyled("Info:", bold=true, color=:cyan)
    println(" ", msg_trunc)

    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    if is_essential_message(msg_trunc) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] $msg_trunc")
        flush(ESSENTIAL_FILE[])
    end

    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [INFO] $msg_trunc")
        flush(CONSOLE_FILE[])
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [info] $msg_trunc")
        flush(DEBUG_FILE[])
    end
end

function user_warn(msg::String, file::String="", line::String="", mod::String="")
    msg_trunc = truncate_for_log(msg)

    printstyled("┌ ", bold=true, color=:yellow)
    printstyled("Warning:", bold=true, color=:yellow)
    println(" ", msg_trunc)

    if !isempty(file) && !isempty(line) && DEBUG_CONSOLE_LEVEL[] > 0
        printstyled("└ ", color=:yellow)
        println("@ $mod $file:$line")
    end

    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    source_loc = isempty(file) || isempty(line) ? "" : " @ $mod $file:$line"

    if is_essential_message(msg_trunc) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] WARNING: $msg_trunc$source_loc")
        flush(ESSENTIAL_FILE[])
    end

    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [WARN] $msg_trunc")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                        └$source_loc")
        end
        flush(CONSOLE_FILE[])
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [warn] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end

    if WARNINGS_FILE[] !== nothing
        println(WARNINGS_FILE[], "[$timestamp] $msg_trunc$source_loc")
        flush(WARNINGS_FILE[])
    end
end

function user_error(msg::String, file::String="", line::String="", mod::String="")
    msg_trunc = truncate_for_log(msg)

    printstyled("┌ ", color=:red)
    printstyled("Error:", bold=true, color=:red)
    println(" ", msg_trunc)

    if !isempty(file) && !isempty(line) && DEBUG_CONSOLE_LEVEL[] > 0
        printstyled("└ ", color=:red)
        println("@ $mod $file:$line")
    end

    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    source_loc = isempty(file) || isempty(line) ? "" : " @ $mod $file:$line"

    if ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] ERROR: $msg_trunc$source_loc")
        flush(ESSENTIAL_FILE[])
    end

    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [ERROR] $msg_trunc")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                         └$source_loc")
        end
        flush(CONSOLE_FILE[])
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [error] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end
end

function user_print(msg::String)
    msg_trunc = truncate_for_log(msg)
    println(msg_trunc)

    if ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], msg_trunc)
        flush(ESSENTIAL_FILE[])
    end

    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], msg_trunc)
        flush(CONSOLE_FILE[])
    end

    if DEBUG_FILE[] !== nothing
        println(DEBUG_FILE[], msg_trunc)
        flush(DEBUG_FILE[])
    end
end

function debug_l1(msg::String, file::String="", line::String="", mod::String="")
    if DEBUG_CONSOLE_LEVEL[] < 1
        return nothing
    end

    msg_trunc = truncate_for_log(msg)
    printstyled("┌ ", bold=true, color=:blue)
    printstyled("Debug:", bold=true, color=:blue)
    println(" ", msg_trunc)

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG1] $msg_trunc")
        flush(DEBUG_FILE[])
    end

    return nothing
end

function debug_l2(msg::String, file::String="", line::String="", mod::String="")
    if DEBUG_CONSOLE_LEVEL[] < 2
        return nothing
    end

    msg_trunc = truncate_for_log(msg)
    printstyled("┌ ", bold=true, color=:blue)
    printstyled("Debug:", bold=true, color=:blue)
    println(" ", msg_trunc)

    if !isempty(file) && !isempty(line)
        printstyled("└ ", color=:blue)
        println("@ $mod $file:$line")
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        source_loc = isempty(file) || isempty(line) ? "" : " @ $mod $file:$line"
        println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG2] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end

    return nothing
end

function debug_l3(msg::String, file::String="", line::String="", mod::String="")
    if DEBUG_CONSOLE_LEVEL[] < 3
        return nothing
    end

    msg_trunc = truncate_for_log(msg)
    printstyled("┌ ", bold=true, color=:blue)
    printstyled("Debug:", bold=true, color=:blue)
    println(" ", msg_trunc)

    if !isempty(file) && !isempty(line)
        printstyled("└ ", color=:blue)
        println("@ $mod $file:$line")
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        source_loc = isempty(file) || isempty(line) ? "" : " @ $mod $file:$line"
        println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG3] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end

    return nothing
end

function trace_msg(msg::String, file::String="", line::String="", mod::String="")
    if DEBUG_CONSOLE_LEVEL[] < 4
        return nothing
    end

    msg_trunc = truncate_for_log(msg)
    printstyled("┌ ", bold=true, color=:blue)
    printstyled("Debug:", bold=true, color=:blue)
    println(" ", msg_trunc)

    if !isempty(file) && !isempty(line)
        printstyled("└ ", color=:blue)
        println("@ $mod $file:$line")
    end

    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        source_loc = isempty(file) || isempty(line) ? "" : " @ $mod $file:$line"
        println(DEBUG_FILE[], "[$debug_timestamp] [TRACE] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end

    return nothing
end

macro user_info(msg)
    runtime_module = @__MODULE__
    :($runtime_module.user_info(string($(esc(msg)))))
end

macro user_warn(msg, kwargs...)
    runtime_module = @__MODULE__
    quote
        $runtime_module.user_warn(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro user_error(msg)
    runtime_module = @__MODULE__
    quote
        $runtime_module.user_error(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro user_print(msg)
    runtime_module = @__MODULE__
    :($runtime_module.user_print(string($(esc(msg)))))
end

macro debug_l1(msg)
    runtime_module = @__MODULE__
    quote
        $runtime_module.debug_l1(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l2(msg)
    runtime_module = @__MODULE__
    quote
        $runtime_module.debug_l2(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l3(msg)
    runtime_module = @__MODULE__
    quote
        $runtime_module.debug_l3(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro trace(msg)
    runtime_module = @__MODULE__
    quote
        $runtime_module.trace_msg(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end
