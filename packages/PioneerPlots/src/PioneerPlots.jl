module PioneerPlots

using Arrow
using DataFrames
using Dictionaries
using GR
using Plots
using PioneerCommon
using SentinelArrays
using StatsPlots
using Tables

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(REPO_ROOT, "src", "structs", "MassSpecData.jl"))
include(joinpath(REPO_ROOT, "src", "structs", "MassErrorModel.jl"))
include(joinpath(REPO_ROOT, "src", "utils", "safeFileOps.jl"))

include(joinpath(@__DIR__, "owned", "pdf_utils.jl"))
include(joinpath(@__DIR__, "owned", "plot_rt_alignment.jl"))
include(joinpath(@__DIR__, "owned", "qc_plots.jl"))

export qcPlots
export plotRTAlign
export save_multipage_pdf

end
