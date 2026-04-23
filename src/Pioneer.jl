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

include(joinpath(@__DIR__, "shared", "version_utils.jl"))
include(joinpath(@__DIR__, "shared", "asset_utils.jl"))
include(joinpath(@__DIR__, "shared", "plot_specs.jl"))
include(joinpath(@__DIR__, "shared", "runtime_core.jl"))
include(joinpath(@__DIR__, "shared", "logging_utils.jl"))

#Set Seed 
Random.seed!(1776);

#Import Pioneer Files 
include("importScripts.jl")
files_loaded = importScripts()

#importScriptsSpecLib(files_loaded)
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")

# H2O, PROTON, NEUTRON constants are defined in get_mz.jl and available via importScripts()
include(joinpath(@__DIR__, "shared", "koina_model_constants.jl"))

function __init__()
    # Don't initialize gr() immediately - let it be initialized when first used
    ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
end

export SearchDIA, BuildSpecLib, GetSearchParams, GetBuildLibParams, convertMzML, # ParseSpecLib, GetParseSpecLibParams, # COMMENTED OUT: ParseSpecLib has loading issues
       @user_info, @user_warn, @user_error, @user_print, @debug_l1, @debug_l2, @debug_l3, @trace,
       get_pioneer_version
end
