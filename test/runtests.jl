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

include(joinpath(dirname(@__DIR__), "src","Routines","LibrarySearch","methods","importScripts.jl"))
importScripts()
@testset "Pioneer.jl" begin
    # Write your tests here.
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
end
