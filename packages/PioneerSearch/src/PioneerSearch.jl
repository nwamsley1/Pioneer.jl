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
using PioneerPlots

const Runtime = PioneerCommon
const Pioneer = @__MODULE__
const SEARCH_APP_NAME = "SearchDIA"
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const Plots = PioneerPlots.Plots

global_logger(ConsoleLogger())
Random.seed!(1776)

function __init__()
    load_pioneer_simd!()
    return nothing
end

plot(args...; kwargs...) = Plots.plot(args...; kwargs...)
plot!(args...; kwargs...) = Plots.plot!(args...; kwargs...)
scatter!(args...; kwargs...) = Plots.scatter!(args...; kwargs...)
histogram(args...; kwargs...) = Plots.histogram(args...; kwargs...)
savefig(args...; kwargs...) = Plots.savefig(args...; kwargs...)
annotate!(args...; kwargs...) = Plots.annotate!(args...; kwargs...)
vline!(args...; kwargs...) = Plots.vline!(args...; kwargs...)
hline!(args...; kwargs...) = Plots.hline!(args...; kwargs...)
xlims(args...; kwargs...) = Plots.xlims(args...; kwargs...)
ylims(args...; kwargs...) = Plots.ylims(args...; kwargs...)
xlims!(args...; kwargs...) = Plots.xlims!(args...; kwargs...)
ylims!(args...; kwargs...) = Plots.ylims!(args...; kwargs...)
text(args...; kwargs...) = Plots.text(args...; kwargs...)

GetParseSpecLibParams(args...; kwargs...) = PioneerParams.GetParseSpecLibParams(args...; kwargs...)
GetSearchParams(args...; kwargs...) = PioneerParams.GetSearchParams(args...; kwargs...)
main_GetParseSpecLibParams(argv=ARGS)::Cint = PioneerParams.main_GetParseSpecLibParams(argv)
main_GetSearchParams(argv=ARGS)::Cint = PioneerParams.main_GetSearchParams(argv)

include(joinpath(@__DIR__, "bootstrap.jl"))
include(joinpath(@__DIR__, "routines", "SearchDIA.jl"))

export SEARCH_APP_NAME
export GetParseSpecLibParams
export GetSearchParams
export SearchDIA
export main_GetParseSpecLibParams
export main_GetSearchParams
export main_SearchDIA

end
