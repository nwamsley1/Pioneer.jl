# Define abstract types and traits
#=
using Arrow, ArrowTypes, ArgParse
using LaTeXStrings
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
using Optim
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
using FastGaussQuadrature
include("structs/Ion.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("utils/ML/uniformBasisCubicSpline.jl")
 include("structs/RetentionTimeConversionModel.jl")
include("utils/quadTransmissionModeling/quadTransmissionModel.jl")
include("utils/isotopeSplines.jl")
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchTypes.jl")
importScripts()
include("Routines/LibrarySearch/SearchMethods/NceTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/SearchMethods.jl")
include("Routines/LibrarySearch/SearchMethods/ParameterTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/FirstPassSearch.jl")
include("structs/Ion.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("utils/ML/uniformBasisCubicSpline.jl")
 include("structs/RetentionTimeConversionModel.jl")
include("utils/quadTransmissionModeling/quadTransmissionModel.jl")
include("utils/isotopeSplines.jl")
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchTypes.jl")
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchMethods.jl")
include("Routines/LibrarySearch/SearchMethods/ParameterTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/FirstPassSearch.jl")

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
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)

=#