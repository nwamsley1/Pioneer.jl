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

include(joinpath(@__DIR__, "routines", "pdf_utils.jl"))
include(joinpath(@__DIR__, "routines", "plot_rt_alignment.jl"))
include(joinpath(@__DIR__, "routines", "qc_plots.jl"))

export qcPlots
export plotRTAlign
export save_multipage_pdf

end
