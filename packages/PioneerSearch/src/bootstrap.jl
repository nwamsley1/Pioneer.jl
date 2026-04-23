function load_search_sources()
    files_loaded = Set{String}()

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

    function safe_include!(file_path::String)::Bool
        if file_path in files_loaded
            return false
        end

        if !isfile(file_path)
            @user_warn "File not found: $file_path"
            return false
        end

        try
            include(file_path)
            push!(files_loaded, file_path)
            return true
        catch e
            @user_warn "Failed to include $file_path: $e"
            return false
        end
    end

    function include_files!(file_dir::String, file_names::Vector{String})
        for fname in file_names
            safe_include!(joinpath(file_dir, fname))
        end
        return nothing
    end

    function safe_include_directory!(dir::String)
        for jfile in get_julia_files(dir)
            safe_include!(jfile)
        end
        return nothing
    end

    include_files!(
        joinpath(REPO_ROOT, "src", "utils", "quadTransmissionModeling"),
        [
            "quadTransmissionModel.jl",
            "generalGaussModel.jl",
            "noQuadModel.jl",
            "RazoQuadModel.jl",
            "SplineQuadModel.jl",
            "binIsotopeRatioData.jl",
            "SquareQuadModel.jl",
        ],
    )

    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "ParseInputs", "paramDefaults.jl"))
    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib", "utils", "buildParamDefaults.jl"))
    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "ParseInputs", "parseParams.jl"))
    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib", "structs", "mods.jl"))

    include_files!(
        joinpath(REPO_ROOT, "src", "structs"),
        [
            "MassSpecData.jl",
            "FilteredMassSpecData.jl",
            "ChromObject.jl",
            "ArrayDict.jl",
            "Counter.jl",
            "Ion.jl",
            "LibraryIon.jl",
            "LibraryProteins.jl",
            "LibraryFragmentIndex.jl",
            "IsotopeTraceType.jl",
            "MatchIon.jl",
            "SparseArray.jl",
            "FragBoundModel.jl",
            "RetentionTimeIndex.jl",
            "MassErrorModel.jl",
            "RetentionTimeConversionModel.jl",
            "protein_inference_types.jl",
        ],
    )

    include_files!(
        joinpath(REPO_ROOT, "src", "utils", "ML"),
        [
            "fdrUtilities.jl",
            "ftrUtilities.jl",
            "lightgbm_utils.jl",
            "probitRegression.jl",
            "piecewiseLinearFunction.jl",
            "spectralLinearRegression.jl",
            "uniformBasisCubicSpline.jl",
            "wittakerHendersonSmoothing.jl",
            "libraryBSpline.jl",
            "psm_container.jl",
            "arrow_psm_container.jl",
            "scoring_traits.jl",
            "scoring_config.jl",
            "pairing.jl",
            "model_training.jl",
            "training_selection.jl",
            "feature_selection.jl",
            "iteration_scheme.jl",
            "mbr_update.jl",
            "scoring_workspace.jl",
            "percolator_generic.jl",
            "percolatorSortOf.jl",
        ],
    )

    include_files!(
        joinpath(REPO_ROOT, "src", "utils"),
        [
            "serialization.jl",
            "isotopes.jl",
            "isotopeSplines.jl",
            "maxLFQ.jl",
            "writeArrow.jl",
            "safeFileOps.jl",
            "proteinInference.jl",
            "profile.jl",
        ],
    )

    safe_include!(joinpath(REPO_ROOT, "src", "utils", "FileOperations", "FileOperations.jl"))

    include_files!(
        joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "PSMs"),
        [
            "PSM.jl",
            "spectralDistanceMetrics.jl",
            "UnscoredPSMs.jl",
            "ScoredPSMs.jl",
        ],
    )

    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))
    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))

    safe_include_directory!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "CommonSearchUtils"))
    safe_include_directory!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "ParseInputs"))

    search_methods_dir = joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "SearchMethods")

    include_files!(
        joinpath(search_methods_dir, "ParameterTuningSearch"),
        ["types.jl", "ParameterTuningSearch.jl", "utils.jl"],
    )

    include_files!(
        joinpath(search_methods_dir, "PrecursorScoringSearch"),
        ["utils.jl", "model_config.jl", "score_psms.jl", "scoring_interface.jl", "PrecursorScoringSearch.jl"],
    )

    include_files!(
        joinpath(search_methods_dir, "ProteinInferenceSearch"),
        ["utils.jl", "ProteinInferenceSearch.jl"],
    )

    include_files!(
        joinpath(search_methods_dir, "ProteinScoringSearch"),
        ["utils.jl", "ProteinScoringSearch.jl"],
    )

    for (root, _, files) in walkdir(search_methods_dir)
        root_basename = basename(root)
        if root_basename in ("ParameterTuningSearch", "PrecursorScoringSearch", "ProteinInferenceSearch", "ProteinScoringSearch")
            continue
        end
        for file in files
            if endswith(file, ".jl")
                safe_include!(joinpath(root, file))
            end
        end
    end

    include_files!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "WriteOutputs"), ["writeCSVTables.jl"])
    safe_include!(joinpath(REPO_ROOT, "src", "Routines", "SearchDIA", "LibrarySearch.jl"))

    return nothing
end

load_search_sources()
