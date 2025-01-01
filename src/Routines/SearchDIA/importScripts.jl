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
    function include_files!(files_loded::Set{String}, file_dir::String, file_names::Vector{String})
        file_paths = [joinpath(file_dir, fname) for fname in file_names]
        [include(fpath) for fpath in file_paths if fpath ∉ files_loaded]
        push!(files_loaded, file_paths...)
        return nothing 
    end
    files_loaded = Set{String}()
    include_files!(
        files_loaded, 
        joinpath(package_root, "src","utils", "quadTransmissionModeling"),
        [
            "quadTransmissionModel.jl",
            "generalGaussModel.jl",
            "noQuadModel.jl",
            "RazoQuadModel.jl",
            "SplineQuadModel.jl",
            "binIsotopeRatioData.jl",
            "squareQuadModel.jl",
        ]
    )
    println("hey")
    include(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))

    include_files!(
        files_loaded, 
        joinpath(package_root, "src","structs"),
        [
            "MassSpecData.jl",
            "ChromObject.jl",
            "ArrayDict.jl",
            "Counter.jl",
            "Ion.jl",
            "LibraryIon.jl",
            "LibraryFragmentIndex.jl",                                                                    
            "IsotopeTraceType.jl",
            "MatchIon.jl",
            "SparseArray.jl",
            "fragBoundModel.jl",
            "RetentionTimeIndex.jl",
            "MassErrorModel.jl",                                                                  
            "RetentionTimeConversionModel.jl"
            ]
    )

    
    # Utilities/ML
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "utils", "ML"),
        [
            "percolatorSortOf.jl",
            "piecewiseLinearFunction.jl",
            "probitRegression.jl",
            "spectralLinearRegression.jl",
            "uniformBasisCubicSpline.jl",
            "wittakerHendersonSmoothing.jl",
            "libraryBSpline.jl"
        ]
    )


    # Utils
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "utils"),
        [
            "isotopes.jl",
            "isotopeSplines.jl",
            "maxLFQ.jl",
                        "writeArrow.jl"
        ]
    )

    # PSMs
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "Routines", "SearchDIA", "PSMs"),
        [
            "PSM.jl",
            "spectralDistanceMetrics.jl",
            "UnscoredPSMs.jl",
            "ScoredPSMs.jl"
        ]
    )

        
    #Search Method 
    include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))

    #[println(fpath) for fpath in collect(files_loaded)]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "ParseInputs")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "WriteOutputs")) if jfile ∉ files_loaded]

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "LibrarySearch.jl"))
end