# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

module Pioneer

using Arrow, ArrowTypes, ArgParse, Dates
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, Combinatorics, CodecZlib
using Serialization
using DataFrames, DataStructures, Dictionaries, Distributions
using EzXML
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML, Logging
using Measures
using NumericalIntegration
using Optim
using Plots, Polynomials, ProgressBars, Printf
using Tables
using StatsPlots, SentinelArrays
using Random
import RobustModels: rlm, TauEstimator, TukeyLoss
import StatsModels: @formula
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
ENV["LIGHTGBM_LOG_LIBRARY_DISCOVERY"] = "false"
using LightGBM
import MLJModelInterface: fit, predict
using KernelDensity
using FastGaussQuadrature
using LaTeXStrings, Printf
using Dates
using InlineStrings
using HTTP


# Simple console logger - detailed logging handled by custom logging system
global_logger(ConsoleLogger())


"""
Type alias for m/z to eV interpolation functions.
Uses GriddedInterpolation with linear interpolation and line extrapolation.
"""
const InterpolationTypeAlias = Interpolations.Extrapolation{
    Float32,  # Value type
    1,        # Dimension
    Interpolations.GriddedInterpolation{
        Float32,                            # Value type
        1,                                  # Dimension
        Vector{Float32},                    # Values
        Gridded{Linear{Throw{OnGrid}}},     # Method
        Tuple{Vector{Float32}}              # Grid type
    },
    Gridded{Linear{Throw{OnGrid}}},         # Method
    Line{Nothing}                           # Extrapolation
}

# ============================================================================
# LOGGING SYSTEM - Global state and functions
# ============================================================================

# Global logging state - four file handles
const ESSENTIAL_FILE = Ref{Union{Nothing, IOStream}}(nothing)  # Clean log (dual_println style)
const CONSOLE_FILE = Ref{Union{Nothing, IOStream}}(nothing)    # Mirror of console
const DEBUG_FILE = Ref{Union{Nothing, IOStream}}(nothing)      # Everything including debug
const WARNINGS_FILE = Ref{Union{Nothing, IOStream}}(nothing)   # All warnings

# Global debug level setting (0 = no debug on console, 1-3 = show debug levels 1-3)
const DEBUG_CONSOLE_LEVEL = Ref{Int}(0)
# Verbosity of pioneer_search_debug.log, independent of the console. Previously the file write lived
# inside the console-level check, so at the default console level of 0 the debug log received nothing
# the console log did not already have -- 6,365 vs 6,263 bytes on a 6-file run, i.e. a duplicate.
# Defaulting this to 1 gives a genuinely useful debug log while the terminal stays quiet.
const DEBUG_FILE_LEVEL = Ref{Int}(1)

# Max bytes of a single log message's content before truncation.
# This applies to the message text; suffix indicating truncation may exceed this cap.
const MAX_LOG_MSG_BYTES = Ref{Int}(4096)

"""
    get_pioneer_version()

Return the Pioneer version string used for user-facing logs/CLI output.
Lookup order:
1. `ENV["PIONEER_VERSION"]` (if set)
2. Bundled `VERSION` file near the executable/app root
3. `Project.toml` version (source/development fallback)
"""
function get_pioneer_version()
    env_version = strip(get(ENV, "PIONEER_VERSION", ""))
    if !isempty(env_version)
        return env_version
    end

    version_candidates = String[]

    # Candidate paths near the currently executing program.
    if !isempty(PROGRAM_FILE)
        exe = PROGRAM_FILE
        if !isabspath(exe)
            exe_full = Sys.which(exe)
            exe = exe_full === nothing ? exe : exe_full
        end
        if ispath(exe)
            exe_real = realpath(exe)
            exe_dir = dirname(exe_real)
            push!(version_candidates, joinpath(exe_dir, "VERSION"))
            push!(version_candidates, joinpath(exe_dir, "..", "VERSION"))
        end
    end

    # Source/development fallbacks.
    push!(version_candidates, joinpath(pkgdir(Pioneer), "VERSION"))
    push!(version_candidates, joinpath(pkgdir(Pioneer), "Project.toml"))

    for path in unique(version_candidates)
        if !isfile(path)
            continue
        end

        try
            if endswith(path, "Project.toml")
                content = read(path, String)
                version_match = match(r"^version[ \t]*=[ \t]*\"([^\"]+)\""m, content)
                if version_match !== nothing
                    return strip(version_match[1])
                end
            else
                v = strip(read(path, String))
                if !isempty(v)
                    return first(split(v, '\n'))
                end
            end
        catch
            # Ignore malformed or unreadable candidates and continue fallback chain.
        end
    end

    return "unknown"
end

"""
    truncate_for_log(msg::String; max_bytes::Int=MAX_LOG_MSG_BYTES[])

Safely truncate a message by bytes without splitting UTF-8 characters.
Adds a suffix indicating how many bytes were omitted.
Returns `msg` unchanged when within the byte limit.
"""
function truncate_for_log(msg::String; max_bytes::Int=MAX_LOG_MSG_BYTES[])
    nb = ncodeunits(msg)
    if nb <= max_bytes
        return msg
    end
    # Walk valid codepoint indices and accumulate bytes up to max_bytes
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
    # Ensure we have at least one full character; otherwise, return minimal safe truncation
    if last_end <= firstindex(msg)
        return string(msg[1:prevind(msg, lastindex(msg))], " … [truncated]")
    end
    prefix = String(@view msg[1:prevind(msg, last_end)])
    omitted = nb - n
    return string(prefix, " … [truncated ", omitted, " bytes]")
end

# List of message patterns that are "essential" (like dual_println)
# These match the messages that would have been output by dual_println
const ESSENTIAL_PATTERNS = [
    r"^Starting search at:",
    r"^Output directory:",
    r"^Loading Parameters",
    r"^Loading Spectral Library",
    r"^Initializing Search Context",
    r"^Executing .+\.\.\.",  # Matches all "Executing [Method]..." messages
    # Note: Performance report and decorative outputs handled by @user_print
]

function is_essential_message(msg::String)
    for pattern in ESSENTIAL_PATTERNS
        if occursin(pattern, msg)
            return true
        end
    end
    return false
end

# Core logging functions - these do the actual work
function user_info(msg::String)
    msg_trunc = truncate_for_log(msg)
    # Console output - match Julia's [ Info: format
    printstyled("[ ", bold=true, color=:cyan)
    printstyled("Info:", bold=true, color=:cyan)
    println(" ", msg_trunc)
    
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Essential file - only for key messages
    if is_essential_message(msg_trunc) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] $msg_trunc")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - all info messages
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [INFO] $msg_trunc")
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [info] $msg_trunc")
        flush(DEBUG_FILE[])
    end
end

function user_warn(msg::String, file::String="", line::String="", mod::String="")
    msg_trunc = truncate_for_log(msg)
    # Console output. Leading "\n" drops past any in-flight ProgressBars line
    # so the warning isn't clobbered by the next bar refresh; trailing newline
    # leaves the cursor on a fresh line for the bar to redraw on.
    print("\n")
    printstyled("┌ ", bold=true, color=:yellow)
    printstyled("Warning:", bold=true, color=:yellow)
    println(" ", msg_trunc)
    
    # Add source location line if available and debug level > 0
    if !isempty(file) && !isempty(line) && DEBUG_CONSOLE_LEVEL[] > 0
        printstyled("└ ", color=:yellow)
        println("@ $mod $file:$line")
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Build source location string for files
    source_loc = ""
    if !isempty(file) && !isempty(line)
        source_loc = " @ $mod $file:$line"
    end
    
    # Essential file - only critical warnings
    if is_essential_message(msg_trunc) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] WARNING: $msg_trunc$source_loc")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console format exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [WARN] $msg_trunc")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                        └$source_loc")  # Indented continuation
        end
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything with full details
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [warn] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end
    
    # Warnings file - all warnings with source for easy debugging
    if WARNINGS_FILE[] !== nothing
        println(WARNINGS_FILE[], "[$timestamp] $msg_trunc$source_loc")
        flush(WARNINGS_FILE[])
    end
end

function user_error(msg::String, file::String="", line::String="", mod::String="")
    msg_trunc = truncate_for_log(msg)
    # Leading "\n" drops past any in-flight ProgressBars line so the error
    # text isn't clobbered by the next bar refresh.
    print("\n")
    printstyled("┌ ", color=:red)
    printstyled("Error:", bold=true, color=:red)
    println(" ", msg_trunc)
    
    # Add source location line if available and debug level > 0
    if !isempty(file) && !isempty(line) && DEBUG_CONSOLE_LEVEL[] > 0
        printstyled("└ ", color=:red)
        println("@ $mod $file:$line")
    end
    
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Build source location string for files
    source_loc = ""
    if !isempty(file) && !isempty(line)
        source_loc = " @ $mod $file:$line"
    end
    
    # Essential file - errors are always essential
    if ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] ERROR: $msg_trunc$source_loc")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console format exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [ERROR] $msg_trunc")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                         └$source_loc")  # Indented continuation
        end
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything with full details
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [error] $msg_trunc$source_loc")
        flush(DEBUG_FILE[])
    end
end

function user_print(msg::String)
    msg_trunc = truncate_for_log(msg)
    # Direct output without formatting
    println(msg_trunc)
    
    # Essential file - ALWAYS gets @user_print messages (like dual_println)
    if ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], msg_trunc)
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], msg_trunc)
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything
    if DEBUG_FILE[] !== nothing
        println(DEBUG_FILE[], msg_trunc)
        flush(DEBUG_FILE[])
    end
end

# Debug logging functions - console output based on DEBUG_CONSOLE_LEVEL
function debug_l1(msg::String, file::String="", line::String="", mod::String="")
    to_console = DEBUG_CONSOLE_LEVEL[] >= 1
    to_file    = _debug_file_level()   >= 1 && DEBUG_FILE[] !== nothing
    (to_console || to_file) || return nothing
    msg_trunc = truncate_for_log(msg)
    if to_console
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg_trunc)
    end
    if to_file
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG1] $msg_trunc")
        # Flushed per message deliberately: the point of this log is post-mortem
        # diagnosis, so buffered lines must not be lost if the run dies.
        flush(DEBUG_FILE[])
    end
    return nothing
end

function debug_l2(msg::String, file::String="", line::String="", mod::String="")
    to_console = DEBUG_CONSOLE_LEVEL[] >= 2
    to_file    = _debug_file_level()   >= 2 && DEBUG_FILE[] !== nothing
    if to_console || to_file
        msg_trunc = truncate_for_log(msg)
        source_loc = (!isempty(file) && !isempty(line)) ? " @ $mod $file:$line" : ""
        if to_console
            printstyled("┌ ", bold=true, color=:blue)
            printstyled("Debug:", bold=true, color=:blue)
            println(" ", msg_trunc)
            isempty(source_loc) || (printstyled("└ ", color=:blue); println(source_loc))
        end
        if to_file
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG2] $msg_trunc$source_loc")
        # Flushed per message deliberately: the point of this log is post-mortem
        # diagnosis, so buffered lines must not be lost if the run dies.
        flush(DEBUG_FILE[])
        end
    end
    # If debug level < 2, no output at all
end

function debug_l3(msg::String, file::String="", line::String="", mod::String="")
    to_console = DEBUG_CONSOLE_LEVEL[] >= 3
    to_file    = _debug_file_level()   >= 3 && DEBUG_FILE[] !== nothing
    if to_console || to_file
        msg_trunc = truncate_for_log(msg)
        source_loc = (!isempty(file) && !isempty(line)) ? " @ $mod $file:$line" : ""
        if to_console
            printstyled("┌ ", bold=true, color=:blue)
            printstyled("Debug:", bold=true, color=:blue)
            println(" ", msg_trunc)
            isempty(source_loc) || (printstyled("└ ", color=:blue); println(source_loc))
        end
        if to_file
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG3] $msg_trunc$source_loc")
        # Flushed per message deliberately: the point of this log is post-mortem
        # diagnosis, so buffered lines must not be lost if the run dies.
        flush(DEBUG_FILE[])
        end
    end
    # If debug level < 3, no output at all
end

function trace_msg(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows (4+)
    if DEBUG_CONSOLE_LEVEL[] >= 4
        msg_trunc = truncate_for_log(msg)
        # Console output WITH line numbers
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg_trunc)
        
        if !isempty(file) && !isempty(line)
            printstyled("└ ", color=:blue)
            println("@ $mod $file:$line")
        end
        
        # File output WITH line numbers (only when debug level allows)
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            source_loc = ""
            if !isempty(file) && !isempty(line)
                source_loc = " @ $mod $file:$line"
            end
            println(DEBUG_FILE[], "[$debug_timestamp] [TRACE] $msg_trunc$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 4, no output at all
end

"""
    _debug_enabled(level) -> Bool

True when a debug message of `level` would reach either destination. Lets the macros skip building
their message string entirely when nothing is listening.
"""
@inline _debug_file_level() = max(DEBUG_FILE_LEVEL[], DEBUG_CONSOLE_LEVEL[])

@inline _debug_enabled(level::Int) =
    DEBUG_CONSOLE_LEVEL[] >= level || (_debug_file_level() >= level && DEBUG_FILE[] !== nothing)

# MACROS - defined once, used everywhere
# These expand at parse time to function calls
macro user_info(msg)
    :(Pioneer.user_info(string($(esc(msg)))))
end

macro user_warn(msg, kwargs...)
    # For now, just ignore extra arguments (like exception=e)
    return quote
        Pioneer.user_warn(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro user_error(msg)
    return quote
        Pioneer.user_error(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro user_print(msg)
    :(Pioneer.user_print(string($(esc(msg)))))
end

macro debug_l1(msg)
    # Guard BEFORE interpolating: the message used to be built at every call site regardless of
    # level, so a suppressed debug line still paid for its own string.
    return quote
        Pioneer._debug_enabled(1) && Pioneer.debug_l1(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l2(msg)
    # Guard BEFORE interpolating: the message used to be built at every call site regardless of
    # level, so a suppressed debug line still paid for its own string.
    return quote
        Pioneer._debug_enabled(2) && Pioneer.debug_l2(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l3(msg)
    # Guard BEFORE interpolating: the message used to be built at every call site regardless of
    # level, so a suppressed debug line still paid for its own string.
    return quote
        Pioneer._debug_enabled(3) && Pioneer.debug_l3(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro trace(msg)
    return quote
        Pioneer.trace_msg(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

# Export the macros for use throughout the codebase
export @user_info, @user_warn, @user_error, @user_print, @debug_l1, @debug_l2, @debug_l3, @trace

"""
    init_pioneer_logging(out_dir, banner_title;
                         debug_console_level=0,
                         max_message_bytes=4096,
                         essential_filename="pioneer_search_report.txt",
                         console_filename="pioneer_search_log.log",
                         debug_filename="pioneer_search_debug.log",
                         warnings_filename="pioneer_warnings.log")

Open the four Pioneer log files inside `out_dir`, write a uniform header
to each, and configure `DEBUG_CONSOLE_LEVEL[]` / `MAX_LOG_MSG_BYTES[]`.

Returns the absolute path of the warnings file (used by
`finalize_pioneer_logging` to count and report warnings).

Used by `SearchDIA` and `BuildSpecLib` so both routines emit the same
log artefacts: a clean essential log, a console mirror, a full debug
trace, and a warnings sidecar.
"""
function init_pioneer_logging(out_dir::AbstractString,
                              banner_title::AbstractString;
                              debug_console_level::Integer = 0,
                              debug_file_level::Integer = 1,
                              max_message_bytes::Integer = 4096,
                              essential_filename::AbstractString = "pioneer_search_report.txt",
                              console_filename::AbstractString = "pioneer_search_log.log",
                              debug_filename::AbstractString = "pioneer_search_debug.log",
                              warnings_filename::AbstractString = "pioneer_warnings.log")
    mkpath(out_dir)

    # Honor an env-var override for the message-byte cap.
    max_bytes = Int(max_message_bytes)
    if haskey(ENV, "PIONEER_MAX_LOG_MSG_BYTES")
        env_val = tryparse(Int, ENV["PIONEER_MAX_LOG_MSG_BYTES"])
        if env_val !== nothing
            max_bytes = env_val
        end
    end
    DEBUG_CONSOLE_LEVEL[] = Int(debug_console_level)
    DEBUG_FILE_LEVEL[]    = Int(debug_file_level)
    MAX_LOG_MSG_BYTES[]   = clamp(max_bytes, 1024, 1048576)

    essential_path = joinpath(out_dir, essential_filename)
    console_path   = joinpath(out_dir, console_filename)
    debug_path     = joinpath(out_dir, debug_filename)
    warnings_path  = joinpath(out_dir, warnings_filename)

    ESSENTIAL_FILE[] = open(essential_path, "w")
    CONSOLE_FILE[]   = open(console_path,   "w")
    DEBUG_FILE[]     = open(debug_path,     "w")
    WARNINGS_FILE[]  = open(warnings_path,  "w")

    pioneer_version = get_pioneer_version()
    banner = repeat("=", 90)
    essential_header = [
        banner,
        banner_title,
        "Version: $pioneer_version",
        "Started: $(Dates.now())",
        "Output Directory: $out_dir",
        banner,
        ""
    ]
    debug_header = [
        banner,
        "$banner_title (Full Trace)",
        "Version: $pioneer_version",
        "Started: $(Dates.now())",
        "Output Directory: $out_dir",
        "Julia Version: $(VERSION)",
        "Threads: $(Threads.nthreads())",
        banner,
        ""
    ]
    warnings_header = [
        banner,
        "$banner_title — Warnings",
        "Started: $(Dates.now())",
        banner,
        ""
    ]
    for line in essential_header
        println(ESSENTIAL_FILE[], line)
        println(CONSOLE_FILE[],   line)
    end
    for line in debug_header
        println(DEBUG_FILE[], line)
    end
    for line in warnings_header
        println(WARNINGS_FILE[], line)
    end

    @user_info "Pioneer logging system initialized"
    return abspath(warnings_path)
end

"""
    finalize_pioneer_logging(warnings_full_path; banner_title="Run completed")

Close the four log files opened by `init_pioneer_logging`, write a
uniform footer (with warning count) to each, and emit a single console
summary line if warnings were generated. Safe to call from a `finally`
block — handles partial init.
"""
function finalize_pioneer_logging(warnings_full_path::AbstractString = "";
                                  banner_title::AbstractString = "Run completed")
    warning_count = 0
    if WARNINGS_FILE[] !== nothing
        close(WARNINGS_FILE[])
        if !isempty(warnings_full_path) && isfile(warnings_full_path)
            # 5-line header (see init_pioneer_logging) is subtracted.
            warning_count = max(0, countlines(warnings_full_path) - 5)
        end
        WARNINGS_FILE[] = nothing
    end

    banner = repeat("=", 90)
    function _close_with_footer!(handle_ref::Ref{Union{Nothing, IOStream}}, include_warning_line::Bool)
        if handle_ref[] !== nothing
            footer = ["", banner]
            if include_warning_line && warning_count > 0
                push!(footer, "⚠️  $warning_count warnings were generated during the run")
            end
            push!(footer, "$banner_title at: $(Dates.now())")
            push!(footer, banner)
            for line in footer
                println(handle_ref[], line)
            end
            close(handle_ref[])
            handle_ref[] = nothing
        end
    end
    _close_with_footer!(ESSENTIAL_FILE, true)
    _close_with_footer!(CONSOLE_FILE,   true)
    _close_with_footer!(DEBUG_FILE,     false)

    if warning_count > 0 && !isempty(warnings_full_path)
        printstyled("┌ ", color=:yellow)
        printstyled("Warning:", bold=true, color=:yellow)
        println(" $warning_count warnings were generated - see $warnings_full_path")
    end
    return warning_count
end

#Set Seed 
Random.seed!(1776);

#Import Pioneer Files 
include("importScripts.jl")
files_loaded = importScripts()

#importScriptsSpecLib(files_loaded)
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))

# H2O, PROTON, and isotope-spacing constants are defined in get_mz.jl and available via importScripts()

# Spectral deconvolution solver defaults (Huber / OLS / Poisson MM coordinate descent)
const DECONV_MAX_ITER::Int64 = Int64(1000)
const DECONV_CONVERGENCE_TOL::Float32 = Float32(0.01)

# AA_to_mass is defined in get_mz.jl and available via importScripts()



const MODEL_CONFIGS = Dict(
    "altimeter" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = SplineCoefficientModel("altimeter"),
        instruments = Set([])
    ),
)


const KOINA_URLS = Dict(
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer",
    "altimeter" => "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_splines_index/infer",#"http://127.0.0.1:8000/v2/models/Altimeter_2024_splines_index/infer"
)

function __init__()
    # Don't initialize gr() immediately - let it be initialized when first used
    ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
    # Force GR's off-screen workstation so plot creation never opens a gksqt
    # window. Plots are rendered to in-memory bits and then combined into the
    # multi-page PDFs by save_multipage_pdf — opening Qt windows here both
    # slows down PDF writes substantially and disrupts the user's desktop
    # during long searches. Set both legacy and modern env var names.
    # "100" (nul) suppresses *all* output, so savefig writes 0-byte PNGs and
    # PDF assembly then fails parsing them. "png" is off-screen but still
    # renders to file.
    get!(ENV, "GKSwstype", "png")
    get!(ENV, "GKS_WSTYPE", "png")
end

export SearchDIA, BuildSpecLib, GetSearchParams, GetBuildLibParams, convertMzML,
       get_pioneer_version,
       @user_info, @user_warn, @user_error, @user_print, @debug_l1, @debug_l2, @debug_l3, @trace
end
