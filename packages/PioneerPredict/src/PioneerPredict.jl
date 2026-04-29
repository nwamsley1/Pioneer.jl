module PioneerPredict

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
using LinearAlgebra, LinearSolve, LightXML, Logging
using NumericalIntegration
using Optim
using Polynomials, ProgressBars, Printf
using Tables
using SentinelArrays
using Random
import RobustModels: rlm, TauEstimator, TukeyLoss
import StatsModels: @formula
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
ENV["LIGHTGBM_LOG_LIBRARY_DISCOVERY"] = "false"
using LightGBM
import MLJModelInterface: fit, predict
using KernelDensity
using FastGaussQuadrature
using InlineStrings
using HTTP
using PioneerCommon
using PioneerParams

const Pioneer = PioneerCommon
const PREDICT_APP_NAME = "BuildSpecLib"
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

global_logger(ConsoleLogger())
Random.seed!(1776)

function __init__()
    load_pioneer_simd!()
    return nothing
end

GetBuildLibParams(args...; kwargs...) = PioneerParams.GetBuildLibParams(args...; kwargs...)
main_GetBuildLibParams(argv=ARGS)::Cint = PioneerParams.main_GetBuildLibParams(argv)

include(joinpath(@__DIR__, "bootstrap.jl"))
include(joinpath(REPO_ROOT, "src", "shared", "koina_model_constants.jl"))
include(joinpath(@__DIR__, "routines", "BuildSpecLib.jl"))

export PREDICT_APP_NAME
export BuildSpecLib
export GetBuildLibParams
export main_BuildSpecLib
export main_GetBuildLibParams

end
