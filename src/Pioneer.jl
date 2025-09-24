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
#using Profile
#using PProf
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, Combinatorics, CodecZlib
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
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
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
    # Console output - match Julia's [ Info: format
    printstyled("[ ", bold=true, color=:cyan)
    printstyled("Info:", bold=true, color=:cyan)
    println(" ", msg)
    
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    
    # Essential file - only for key messages
    if is_essential_message(msg) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] $msg")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - all info messages
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [INFO] $msg")
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [info] $msg")
        flush(DEBUG_FILE[])
    end
end

function user_warn(msg::String, file::String="", line::String="", mod::String="")
    # Console output
    printstyled("┌ ", bold=true, color=:yellow)
    printstyled("Warning:", bold=true, color=:yellow)
    println(" ", msg)
    
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
    if is_essential_message(msg) && ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], "[$timestamp] WARNING: $msg$source_loc")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console format exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [WARN] $msg")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                        └$source_loc")  # Indented continuation
        end
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything with full details
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [warn] $msg$source_loc")
        flush(DEBUG_FILE[])
    end
    
    # Warnings file - all warnings with source for easy debugging
    if WARNINGS_FILE[] !== nothing
        println(WARNINGS_FILE[], "[$timestamp] $msg$source_loc")
        flush(WARNINGS_FILE[])
    end
end

function user_error(msg::String, file::String="", line::String="", mod::String="")
    # Console output
    printstyled("┌ ", color=:red)
    printstyled("Error:", bold=true, color=:red)
    println(" ", msg)
    
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
        println(ESSENTIAL_FILE[], "[$timestamp] ERROR: $msg$source_loc")
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console format exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], "[$timestamp] [ERROR] $msg")
        if !isempty(source_loc)
            println(CONSOLE_FILE[], "                         └$source_loc")  # Indented continuation
        end
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything with full details
    if DEBUG_FILE[] !== nothing
        debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
        println(DEBUG_FILE[], "[$debug_timestamp] [error] $msg$source_loc")
        flush(DEBUG_FILE[])
    end
end

function user_print(msg::String)
    # Direct output without formatting
    println(msg)
    
    # Essential file - ALWAYS gets @user_print messages (like dual_println)
    if ESSENTIAL_FILE[] !== nothing
        println(ESSENTIAL_FILE[], msg)
        flush(ESSENTIAL_FILE[])
    end
    
    # Console file - mirror console exactly
    if CONSOLE_FILE[] !== nothing
        println(CONSOLE_FILE[], msg)
        flush(CONSOLE_FILE[])
    end
    
    # Debug file - everything
    if DEBUG_FILE[] !== nothing
        println(DEBUG_FILE[], msg)
        flush(DEBUG_FILE[])
    end
end

# Debug logging functions - console output based on DEBUG_CONSOLE_LEVEL
function debug_l1(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 1
        # Console output WITHOUT line numbers
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        # NO LINE NUMBER OUTPUT for debug_l1
        
        # File output WITHOUT line numbers (only when debug level allows)
        if DEBUG_FILE[] !== nothing
            debug_timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG1] $msg")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 1, no output at all
end

function debug_l2(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 2
        # Console output WITH line numbers
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
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
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG2] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 2, no output at all
end

function debug_l3(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows
    if DEBUG_CONSOLE_LEVEL[] >= 3
        # Console output WITH line numbers
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
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
            println(DEBUG_FILE[], "[$debug_timestamp] [DEBUG3] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 3, no output at all
end

function trace_msg(msg::String, file::String="", line::String="", mod::String="")
    # Only process if debug level allows (4+)
    if DEBUG_CONSOLE_LEVEL[] >= 4
        # Console output WITH line numbers
        printstyled("┌ ", bold=true, color=:blue)
        printstyled("Debug:", bold=true, color=:blue)
        println(" ", msg)
        
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
            println(DEBUG_FILE[], "[$debug_timestamp] [TRACE] $msg$source_loc")
            flush(DEBUG_FILE[])
        end
    end
    # If debug level < 4, no output at all
end

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
    return quote
        Pioneer.debug_l1(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l2(msg)
    return quote
        Pioneer.debug_l2(
            string($(esc(msg))),
            $(string(__source__.file)),
            $(string(__source__.line)),
            $(string(__module__))
        )
    end
end

macro debug_l3(msg)
    return quote
        Pioneer.debug_l3(
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

#Set Seed 
Random.seed!(1776);

#Import Pioneer Files 
include("importScripts.jl")
files_loaded = importScripts()

#importScriptsSpecLib(files_loaded)
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")
const CHARGE_ADJUSTMENT_FACTORS = Float64[1, 0.9, 0.85, 0.8, 0.75]

# H2O, PROTON, NEUTRON constants are defined in get_mz.jl and available via importScripts()
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)

# AA_to_mass is defined in get_mz.jl and available via importScripts()



const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE","QEHFX","LUMOS","ELITE","VELOS","NONE"])
    ),
    "altimeter" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = SplineCoefficientModel("altimeter"),
        instruments = Set([])
    ),
    "prosit_2020_hcd" => (
        annotation_type = GenericFragAnnotation("y1+1"), 
        model_type = InstrumentAgnosticModel("prosit_2020_hcd"),
        instruments = Set([])
    ),
    "AlphaPeptDeep" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentSpecificModel("AlphaPeptDeep"),
        instruments = Set(["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"])
    )
)


const KOINA_URLS = Dict(
    "unispec" => "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
    "prosit_2020_hcd" => "https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
    "AlphaPeptDeep" => "https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer",
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer",
    "altimeter" => "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024_splines_index/infer",#"http://127.0.0.1:8000/v2/models/Altimeter_2024_splines_index/infer"
)

function __init__()
    # Don't initialize gr() immediately - let it be initialized when first used
    ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
end

export SearchDIA, BuildSpecLib, GetSearchParams, GetBuildLibParams, convertMzML, # ParseSpecLib, GetParseSpecLibParams, # COMMENTED OUT: ParseSpecLib has loading issues
       @user_info, @user_warn, @user_error, @user_print, @debug_l1, @debug_l2, @debug_l3, @trace
end
