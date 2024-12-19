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
using Optim
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
using FastGaussQuadrature
using LaTeXStrings, Printf
using SparseArrays

main_dir = joinpath(@__DIR__, "../src")
include(joinpath(dirname(@__DIR__), "src", "Routines","SearchDIA","importScripts.jl"))
importScripts()
#include(joinpath(main_dir, "Routines","LibrarySearch","methods","loadSpectralLibrary.jl"))
#const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(dirname(@__DIR__), "src", "Routines","SearchDIA.jl"))
include(joinpath(dirname(@__DIR__), "src", "Routines","ThreeProteomeAnalysis.jl"))
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
export SearchDIA, ThreeProteomeAnalysis, BuildSpecLib
#This is an important alias.
const Pioneer = Main
using Test
results_dir = joinpath(@__DIR__, "../data/ecoli_test/ecoli_test_results")
if isdir(results_dir)
    # Delete all files and subdirectories within the directory
    for item in readdir(results_dir, join=true)
        rm(item, force=true, recursive=true)
    end
end
@testset "Pioneer.jl" begin
   @testset "process_test" begin 
        @test SearchDIA("./data/ecoli_test/ecoli_test_params.json")===nothing
    end
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
end
