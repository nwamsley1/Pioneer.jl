using Arrow, ArrowTypes, ArgParse
#using BSplineKit Don't need this imports anymore?
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, CategoricalArrays, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries #, Distributions 
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures
using NumericalIntegration
using Optim
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots, SentinelArrays
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
using FastGaussQuadrature
using LaTeXStrings, Printf
using SparseArrays
using Dates 

main_dir = joinpath(@__DIR__, "../src")
include(joinpath(dirname(@__DIR__), "src", "Routines","SearchDIA","importScripts.jl"))
importScripts()

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

include(joinpath(dirname(@__DIR__), "src", "Routines","BuildSpecLib","importScripts.jl"))
importScripts()
#include(joinpath(main_dir, "Routines","LibrarySearch","methods","loadSpectralLibrary.jl"))
#const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(dirname(@__DIR__), "src", "Routines","SearchDIA.jl"))
include(joinpath(dirname(@__DIR__), "src", "Routines","BuildSpecLib.jl"))
include(joinpath(dirname(@__DIR__), "src", "Routines","ParseSpecLib.jl"))
const CHARGE_ADJUSTMENT_FACTORS = Float64[1, 0.9, 0.85, 0.8, 0.75]

const H2O::Float64 = Float64(18.010565)
const PROTON::Float64 = Float64(1.0072764)
const NEUTRON::Float64 = Float64(1.00335)
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)


const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE","QEHFX","LUMOS","ELITE","VELOS","NONE"])
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
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
)

export SearchDIA, BuildSpecLib
#This is an important alias.
const Pioneer = Main
using Test
results_dir = joinpath(@__DIR__, "../data/ecoli_test/ecoli_test_results")
if isdir(results_dir)
    # Delete all files and subdirectories within the directory
    #for item in readdir(results_dir, join=true)
    #    rm(item, force=true, recursive=true)
    #end
end
@testset "Pioneer.jl" begin
    println("dir ", @__DIR__)
    @testset "process_test" begin 
        @test SearchDIA("./../data/ecoli_test/ecoli_test_params.json")===nothing
    end
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
    include("./UnitTests/empiricalLibTests.jl")
end
