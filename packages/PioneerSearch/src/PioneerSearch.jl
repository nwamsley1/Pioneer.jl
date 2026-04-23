module PioneerSearch

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
using LaTeXStrings
using InlineStrings
using HTTP
using PioneerCommon

const Pioneer = PioneerCommon
const SEARCH_APP_NAME = "SearchDIA"
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

global_logger(ConsoleLogger())
Random.seed!(1776)

include(joinpath(@__DIR__, "bootstrap.jl"))
include(joinpath(@__DIR__, "owned", "GenerateParams.jl"))
include(joinpath(@__DIR__, "owned", "SearchDIA.jl"))

export SEARCH_APP_NAME
export GetParseSpecLibParams
export GetSearchParams
export SearchDIA
export main_GetParseSpecLibParams
export main_GetSearchParams
export main_SearchDIA

end
