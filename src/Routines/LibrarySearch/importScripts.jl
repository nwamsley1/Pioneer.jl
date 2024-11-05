function importScripts()
    package_root = dirname(dirname(dirname(@__DIR__)))
    
    [include(joinpath(package_root, "src","utils", "quadTransmissionModeling", jl_file)) for jl_file in [
        "quadTransmissionModel.jl",
        "generalGaussModel.jl",
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
                                                                    "MatchIon.jl",
                                                                    "LibraryFragmentIndex.jl",
                                                                    "SparseArray.jl",
                                                                    "fastaEntry.jl",
                                                                    "fragBoundModel.jl",
                                                                    "modelTypes.jl",
                                                                    "RetentionTimeIndex.jl",
                                                                    "MassErrorModel.jl"]];
    #Utilities
    [include(joinpath(package_root, "src","utils", "ML", jl_file)) for jl_file in [
        "percolatorSortOf.jl",
        "probitRegression.jl",
        "spectralLinearRegression.jl",
        "uniformBasisCubicSpline.jl",
        "wittakerHendersonSmoothing.jl"
    ]];   

    [include(joinpath(package_root, "src","utils", jl_file)) for jl_file in [
        "isotopes.jl",
        "isotopeSplines.jl",
        "maxLFQ.jl"
    ]];   



    [include(joinpath(package_root,"src","Routines","LibrarySearch","PSMs", jl_file)) for jl_file in [
        "PSM.jl",
        "spectralDistanceMetrics.jl",
        "UnscoredPSMs.jl",
        "ScoredPSMs.jl"]];

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","FirstSearch",jl_file)) for jl_file in [
            "firstSearch.jl",
            "getBestPSMs.jl",
            "scoreMainSearch.jl"]];


    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "normalizeQuant.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "queryFragmentIndex.jl"]];
    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","ParameterTuning",jl_file)) for jl_file in [
            "addPreSearchColumns.jl",
            "huberLossSearch.jl",
            "mapLibraryToEmpiricalRT.jl",
            "massErrorEstimation.jl",
            "parameterTuningSearch.jl",
            "quadTuningSearch.jl"]];

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","ParseInputs",jl_file)) for jl_file in [
            "getCVFolds.jl",
            "loadSpectralLibrary.jl",
            "parseFileNames.jl",
            "parseParams.jl"]];

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","ProteinGroups",jl_file)) for jl_file in [
            "proteinQuant.jl",
            "scoreProteinGroups.jl"]];

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","QuantitativeSearch",jl_file)) for jl_file in [
            "addSecondSearchColumns.jl",
            "integrateChroms.jl",
            "quantitativeSearch.jl",
            "samplePsmsForXgboost.jl",
            "secondQuant.jl"]];

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","SummarizeRuns",jl_file)) for jl_file in [
            "buildRTIndex.jl",
            "entrapmentAnalysis.jl",
            "getBestPrecursorsAccrossRuns.jl",
            "getBestTrace.jl",
            "getIrtErrs.jl",
            "getPSMsPassingQVal.jl",
            "makeOutputDirectories.jl",
            "mergePsmTables.jl",
            "summarizeToProtein.jl"]];    

    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch","WriteOutputs",jl_file)) for jl_file in [
            "plotRTAlignment.jl",
            "qcPlots.jl",
            "writeCSVTables.jl"]];                    
    #Files needed for PRM routines
    [include(joinpath(package_root,"src","Routines","LibrarySearch",jl_file)) for jl_file in [
                                                                                    "paramsChecks.jl",
                                                                                    "partitionThreadTasks.jl",
                                                                                    "scoreTraces.jl",
                                                                                    "searchRAW.jl"]];                                                                                      

                                                                                                                                
    [include(joinpath(package_root,"src","Routines","BuildSpecLib",jl_file)) for jl_file in [
        "PioneerLib.jl",
        "buildPioneerLib.jl",  
        "buildUniSpecInput.jl",
        "estimateCollisionEv.jl",
        "fragBounds.jl",
        "getIonAnnotations.jl",
        "getMZ.jl",
        "koinaRequests.jl",
        "paramsChecks.jl",
        "parseChronologerResults.jl",
        "parseFasta.jl",
        "parseIonAnnotations.jl",
        "parseIsotopeMods.jl",
        "parseKoinaFragments.jl",
        "prepareChronologerInput.jl"
        ]];

end