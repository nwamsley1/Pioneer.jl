function importScripts()
    package_root = dirname(dirname(dirname(@__DIR__)))
    #package_root = dirname(@__DIR__)
    function get_julia_files(dir::String)
        julia_files = String[]
        for (root, _, files) in walkdir(dir)
            for file in files
                if endswith(file, ".jl")
                    push!(julia_files, joinpath(root, file))
                end
            end
        end
        return julia_files
    end
    [include(joinpath(package_root, "src","utils", "quadTransmissionModeling", jl_file)) for jl_file in [
        "quadTransmissionModel.jl",
        "generalGaussModel.jl",
        "NoQuadModel.jl",
        "RazoQuadModel.jl",
        "SplineQuadModel.jl",
        "binIsotopeRatioData.jl",
        "SquareQuadModel.jl",
    
    ]];   

    [include(joinpath(package_root, "src", "structs", jl_file)) for jl_file in [
                                                                    "ChromObject.jl",
                                                                    "ArrayDict.jl",
                                                                    "Counter.jl",
                                                                    "Ion.jl",
                                                                    "LibraryIon.jl",
                                                                    "LibraryFragmentIndex.jl",                                                                    
                                                                    "IsotopeTraceType.jl",
                                                                    "MatchIon.jl",
                                                                    "SparseArray.jl",
                                                                    "fastaEntry.jl",
                                                                    "fragBoundModel.jl",
                                                                    "modelTypes.jl",
                                                                    "RetentionTimeIndex.jl",
                                                                    "MassErrorModel.jl",                                                                  
                                                                    "RetentionTimeConversionModel.jl"
                                                                    ]];
  
        #Utilities
        [include(joinpath(package_root, "src","utils", "ML", jl_file)) for jl_file in [
                "percolatorSortOf.jl",
                "piecewiseLinearFunction.jl",
                "probitRegression.jl",
                "spectralLinearRegression.jl",
                "uniformBasisCubicSpline.jl",
                "wittakerHendersonSmoothing.jl",
                "libraryBSpline.jl"
                ]]; 

    [include(joinpath(package_root, "src","utils", jl_file)) for jl_file in [
        "isotopes.jl",
        "isotopeSplines.jl",
        "maxLFQ.jl"
    ]];   

    [include(joinpath(package_root,"src","Routines","SearchDIA","PSMs", jl_file)) for jl_file in [
        "PSM.jl",
        "spectralDistanceMetrics.jl",
        "UnscoredPSMs.jl",
        "ScoredPSMs.jl"]];

        
    #Search Method 
    include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils"))]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "ParseInputs"))]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods"))]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "WriteOutputs"))]

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "LibrarySearch.jl"))
end