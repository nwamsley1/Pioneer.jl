#=
using Arrow, ArgParse
#using BSplineKit Don't need this imports anymore?


using Base64
using Base.Order
using Base.Iterators: partition
using CSV, CategoricalArrays, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries, Distributions 
using ExpectationMaximization
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures
using NumericalIntegration
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars
using Tables
using StatsPlots
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
using XGBoost

using Pioneer
using Test


@testset "Pioneer.jl" begin
    # Write your tests here.
    @testset "process_test" begin 
        @test SearchDIA("data/ecoli_test/ecoli_test_params.json")==10
    end
    include(joinpath(dirname(@__DIR__), "src","Routines","LibrarySearch","methods","importScripts.jl"))
    importScripts()
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
end
=#

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
using StaticArrays, StatsBase, SpecialFunctions, Statistics, SparseArrays
using XGBoost


main_dir = joinpath(@__DIR__, "../src")
include(joinpath(main_dir, "Routines","LibrarySearch","methods","importScripts.jl"))
importScripts()
include(joinpath(main_dir, "Routines","LibrarySearch","methods","loadSpectralLibrary.jl"))
const methods_path = joinpath(main_dir, "Routines","LibrarySearch")       
include(joinpath(main_dir, "Routines","SearchDIA.jl"))
include(joinpath(main_dir, "Routines","ThreeProteomeAnalysis.jl"))
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

using Test
@testset "Pioneer.jl" begin
    # Write your tests here.
    @testset "process_test" begin 
        @test SearchDIA("data/ecoli_test/ecoli_test_params.json")==10
    end
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
end
