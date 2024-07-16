module Pioneer
    using Arrow, ArgParse
    using BSplineKit, Base64
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
    using StaticArrays, StatsBase, SpecialFunctions, Statistics
    using XGBoost
    #create_app("../Pioneer","../Pioneer_Compiled", force = true)


    [include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in [
                                                                    "ChromObject.jl",
                                                                    "ArrayDict.jl",
                                                                    "Counter.jl",
                                                                    "Ion.jl",
                                                                    "LibraryIon.jl",
                                                                    "MatchIon.jl",
                                                                    "LibraryFragmentIndex.jl",
                                                                    "SparseArray.jl"]];

    #Utilities
    [include(joinpath("Utils", jl_file)) for jl_file in [

                                                                    "isotopes.jl",
                                                                    "globalConstants.jl",
                                                                    "uniformBasisCubicSpline.jl",
                                                                    "isotopeSplines.jl",
                                                                    "normalizeQuant.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "SpectralDeconvolution.jl",
                                                                    "percolatorSortOf.jl",
                                                                    "plotRTAlignment.jl",
                                                                    "probitRegression.jl",
                                                                    "partitionThreadTasks.jl",
                                                                    "LFQ.jl",
                                                                    "scoreProteinGroups.jl",
                                                                    "wittakerHendersonSmoothing.jl",
                                                                    "getBestTrace.jl"]];

    [include(joinpath("PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

    #Files needed for PRM routines
    [include(joinpath("Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl"]];
    include("Routines/LibrarySearch/MAIN.jl")     
    const methods_path = joinpath(pwd(), "src","Routines","LibrarySearch")                  
    function julia_main()::Cint
        pioneer_main()
        return 0
    end

end
