function importScripts()

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
    [include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in [

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

    [include(joinpath(pwd(), "src","PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

    #Files needed for PRM routines
    [include(joinpath(pwd(), "src", "Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "parseFileNames.jl",
                                                                                    "makeOutputDirectories.jl",
                                                                                    "parseParams.jl",
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl"]];

    #Files needed for PRM routines
    [include(joinpath(pwd(), "src", "Routines","LibrarySearch",jl_file)) for jl_file in [
                                                                                    "MAIN.jl",
                                                                                    "parameterTuningSearch.jl",
                                                                                    "firstSearch.jl",
                                                                                    "quantitativeSearch.jl",
                                                                                    "scoreTraces.jl",
                                                                                    "secondQuant.jl",
                                                                                    "proteinQuant.jl",
                                                                                    "qcPlots.jl"]]; 

end