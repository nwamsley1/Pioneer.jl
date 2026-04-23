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

const Runtime = PioneerCommon
const Pioneer = @__MODULE__
const SEARCH_APP_NAME = "SearchDIA"
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PIONEER_PLOTS_PKGID = Base.PkgId(Base.UUID("9adec8df-eb60-4bac-9b0a-4137287d3bbd"), "PioneerPlots")
const PIONEER_PLOTS_MODULE = Ref{Union{Nothing, Module}}(nothing)

global_logger(ConsoleLogger())
Random.seed!(1776)

function get_pioneer_plots_module()::Module
    plots_module = PIONEER_PLOTS_MODULE[]
    if plots_module === nothing
        Base.require(PIONEER_PLOTS_PKGID)
        plots_module = Base.root_module(PIONEER_PLOTS_PKGID)
        PIONEER_PLOTS_MODULE[] = plots_module
    end
    return plots_module
end

module LazyPlotsProxy
function plot(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.plot(args...; kwargs...)
end
function plot!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.plot!(args...; kwargs...)
end
function scatter!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.scatter!(args...; kwargs...)
end
function histogram(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.histogram(args...; kwargs...)
end
function savefig(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.savefig(args...; kwargs...)
end
function annotate!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.annotate!(args...; kwargs...)
end
function vline!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.vline!(args...; kwargs...)
end
function hline!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.hline!(args...; kwargs...)
end
function xlims(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.xlims(args...; kwargs...)
end
function ylims(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.ylims(args...; kwargs...)
end
function xlims!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.xlims!(args...; kwargs...)
end
function ylims!(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.ylims!(args...; kwargs...)
end
function text(args...; kwargs...)
    return parentmodule(@__MODULE__).get_pioneer_plots_module().Plots.text(args...; kwargs...)
end
end

const Plots = LazyPlotsProxy

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
plotRTAlign(args...; kwargs...) = get_pioneer_plots_module().plotRTAlign(args...; kwargs...)
qcPlots(args...; kwargs...) = get_pioneer_plots_module().qcPlots(args...; kwargs...)
save_multipage_pdf(args...; kwargs...) = get_pioneer_plots_module().save_multipage_pdf(args...; kwargs...)

include(joinpath(@__DIR__, "bootstrap.jl"))
include(joinpath(@__DIR__, "owned", "SearchDIA.jl"))

export SEARCH_APP_NAME
export GetParseSpecLibParams
export GetSearchParams
export SearchDIA
export main_GetParseSpecLibParams
export main_GetSearchParams
export main_SearchDIA

end
