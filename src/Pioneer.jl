module Pioneer
__precompile__(false)
using Arrow, ArrowTypes, ArgParse
#using BSplineKit Don't need this imports anymore?
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, CategoricalArrays, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries, Distributions 
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures
using NumericalIntegration
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
#create_app("../Pioneer","../Pioneer_Compiled", force = true)
#Inport Pioneer Files 
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(@__DIR__, "Routines","SearchDIA.jl"))
include(joinpath(@__DIR__, "Routines","ThreeProteomeAnalysis.jl"))
const charge_facs = Float64[1, 0.9, 0.85, 0.8, 0.75]
const prediction_model_options =  Set(["unispec","prosit_2020_hcd","AlphaPeptDeep"])
const prediction_model_to_annotation_type = Dict(
    "unispec" => UniSpecFragAnnotation("y1^1"),
    "prosit_2020_hcd" => GenericFragAnnotation("y1+1"),
    "AlphaPeptDeep" => GenericFragAnnotation("y1+1")
)
const prediction_model_to_model_type = Dict(
    "unispec" => InstrumentSpecificModel("unispec"),
    "prosit_2020_hcd" => InstrumentAgnosticModel("prosit_2020_hcd"),
    "AlphaPeptDeep" => InstrumentSpecificModel("AlphaPeptDeep")
)
const H2O::Float64 = Float64(18.010565)
const PROTON::Float64 = Float64(1.0072764)
const NEUTRON::Float64 = Float64(1.00335)



export SearchDIA, ThreeProteomeAnalysis, BuildSpecLib
end
