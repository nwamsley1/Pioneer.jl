module Pioneer
__precompile__(false)
using Arrow, ArrowTypes, ArgParse
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
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots, SentinelArrays
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
using FastGaussQuadrature
using LaTeXStrings, Printf
using Dates
using InlineStrings
using HTTP
#Set Seed 
Random.seed!(1776);

#Inport Pioneer Files 
include(joinpath(@__DIR__, "Routines","SearchDIA","importScripts.jl"))
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


include(joinpath(@__DIR__, "Routines","BuildSpecLib","importScripts.jl"))
importScripts()
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(@__DIR__, "Routines","SearchDIA.jl"))
include(joinpath(@__DIR__, "Routines","BuildSpecLib.jl"))
include(joinpath(@__DIR__, "Routines","ParseSpecLib.jl"))
include(joinpath(@__DIR__, "Routines","GenerateParams.jl"))
include(joinpath(@__DIR__, "Routines","mzmlConverter","convertMzML.jl"))
const CHARGE_ADJUSTMENT_FACTORS = Float64[1, 0.9, 0.85, 0.8, 0.75]

const H2O::Float64 = Float64(18.010565)
const PROTON::Float64 = Float64(1.0072764)
const NEUTRON::Float64 = Float64(1.00335)
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)

const AA_to_mass::Dict{Char, Float64} = Dict{Char, Float64}(
        'A' => 71.03711,
        'R' => 156.10111,
        'N' => 114.04293,
        'D' => 115.02694,
        'C' => 103.00919,
        'E' => 129.04259,
        'Q' => 128.05858,
        'G' => 57.02146,
        'H' => 137.05891,
        'I' => 113.08406,
        'L' => 113.08406,
        'K' => 128.09496,
        'M' => 131.04049,
        'F' => 147.06841,
        'P' => 97.05276,
        'S' => 87.03203,
        'T' => 101.04768,
        'W' => 186.07931,
        'Y' => 163.06333,
        'V' => 99.06841,
        'U' => 150.95363,
        'O' => 237.14773
        )



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
    "altimeter" => "http://127.0.0.1:8000/v2/models/Altimeter_2024_splines/infer",
)



export SearchDIA, BuildSpecLib, ParseSpecLib, GetSearchParams, GetBuildLibParams, convertMzML
end