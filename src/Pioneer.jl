module Pioneer
#__precompile__(false)
using Arrow, ArrowTypes, ArgParse
using Profile
using PProf
#using BSplineKit Don't need this imports anymore?
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, CategoricalArrays, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries#, Distributions 
using EzXML
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures
using NumericalIntegration
using Optim
using Plots, PrettyPrinting, Polynomials, PDFmerger, Profile, ProgressBars, Pkg, Printf
using Tables, Test
using StatsPlots, SentinelArrays
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
using XGBoost
using KernelDensity
using FastGaussQuadrature
using LaTeXStrings, Printf
using Dates
using InlineStrings
using HTTP
#Set Seed 
Random.seed!(1776);

#Import Pioneer Files 
include(joinpath(@__DIR__, "Routines","SearchDIA","importScripts.jl"))
files_loaded = importScripts()

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


include(joinpath(@__DIR__, "Routines","BuildSpecLib","importScripts.jl"))
importScriptsSpecLib(files_loaded)
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(@__DIR__, "Routines","SearchDIA.jl"))
include(joinpath(@__DIR__, "Routines","BuildSpecLib.jl"))
include(joinpath(@__DIR__, "Routines","ParseSpecLib.jl"))
include(joinpath(@__DIR__, "Routines","GenerateParams.jl"))
include(joinpath(@__DIR__, "Routines","mzmlConverter","convertMzML.jl"))
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
    "altimeter" => "http://127.0.0.1:8000/v2/models/Altimeter_2024_splines_index/infer",
)



export SearchDIA, BuildSpecLib, ParseSpecLib, GetSearchParams, GetBuildLibParams, convertMzML, safeRm
end