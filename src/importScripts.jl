# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

function importScripts()
    package_root = dirname(@__DIR__)
    
    # Simple safe import system without world age issues
    files_loaded = Set{String}()
    total_attempts = 0
    successful_includes = 0
    conflicts_skipped = String[]
    
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
        total_attempts += 1
        
        # Check if file already loaded
        if file_path in files_loaded
            return false
        end
        
        # Check if file exists
        if !isfile(file_path)
            @warn "File not found: $file_path"
            return false
        end
        
        # Include the file
        try
            include(file_path)
            push!(files_loaded, file_path)
            successful_includes += 1
            return true
        catch e
            @warn "Failed to include $file_path: $e"
            return false
        end
    end
    
    function include_files!(file_dir::String, file_names::Vector{String})
        file_paths = [joinpath(file_dir, fname) for fname in file_names]
        successful_count = 0
        for fpath in file_paths
            if safe_include!(fpath)
                successful_count += 1
            end
        end
        return successful_count
    end
    
    function safe_include_directory!(dir::String)
        files = get_julia_files(dir)
        successful_count = 0
        for jfile in files
            if safe_include!(jfile)
                successful_count += 1
            end
        end
        return successful_count
    end
    
    # Include files using the safe import system
    include_files!(
        joinpath(package_root, "src","utils", "quadTransmissionModeling"),
        [
            "quadTransmissionModel.jl",
            "generalGaussModel.jl",
            "noQuadModel.jl",
            "RazoQuadModel.jl",
            "SplineQuadModel.jl",
            "binIsotopeRatioData.jl",
            "SquareQuadModel.jl"
        ]
    )

    safe_include!(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))
    safe_include!(joinpath(package_root, "src", "Routines","BuildSpecLib", "structs", "mods.jl"))
    
    include_files!(
        joinpath(package_root, "src","structs"),
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
            "protein_inference_types.jl"
            ]
    )

    
    # Utilities/ML
    include_files!(
        joinpath(package_root, "src", "utils", "ML"),
        [
            "fdrUtilities.jl",
            "ftrUtilities.jl",
            "percolatorSortOf.jl",
            "piecewiseLinearFunction.jl",
            "probitRegression.jl",
            "spectralLinearRegression.jl",
            "uniformBasisCubicSpline.jl",
            "wittakerHendersonSmoothing.jl",
            "libraryBSpline.jl"
        ]
    )


    # Utils (must load writeArrow before FileOperations)
    include_files!(
        joinpath(package_root, "src", "utils"),
        [
            "isotopes.jl",
            "isotopeSplines.jl",
            "maxLFQ.jl",
            "writeArrow.jl",
            "safeFileOps.jl",
            "proteinInference.jl",
            "profile.jl",
            "pdfUtils.jl"
        ]
    )

    # Include new FileOperations module from utils (after writeArrow is loaded)
    safe_include!(joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl"))

    # PSMs
    include_files!(
        joinpath(package_root, "src", "Routines", "SearchDIA", "PSMs"),
        [
            "PSM.jl",
            "spectralDistanceMetrics.jl",
            "UnscoredPSMs.jl",
            "ScoredPSMs.jl"
        ]
    )

        
    #Search Method 
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))

    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))

    # Include remaining files using safe import for directories
    safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils"))
    safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "ParseInputs"))
    
    # LoggingSystem - Load BEFORE SearchMethods so macros are available
    include_files!(
        joinpath(package_root, "src", "Routines", "SearchDIA", "LoggingSystem"),
        [
            "WarningTracker.jl",
            "WarningCapturingLogger.jl",
            "Configuration.jl",
            "AsyncFileLogger.jl",  # Keep for future use
            "SimpleFileLogger.jl",  # Keep for comparison
            "DualPrintLogger.jl",  # New minimal logger like dual_println
            "Loggers.jl",
            "Macros.jl",
            "ProgressIntegration.jl",
            "FileManagement.jl",
            "LoggingSystem.jl"  # Main orchestration file should be last
        ]
    )
    
    # SearchMethods (excluding the old FileReferences.jl and FileOperations.jl files)
    search_methods_dir = joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods")
    
    # Load ParameterTuningSearch files in dependency order
    include_files!(
        joinpath(search_methods_dir, "ParameterTuningSearch"),
        [
            "types.jl",                    # All type definitions to avoid circular dependencies
            "diagnostics.jl",              # Diagnostic functions (types moved to types.jl)
            "cross_run_learning.jl",       # Cross-run learning functions (types moved to types.jl)
            "ParameterTuningSearch.jl",    # Main implementation (types moved to types.jl)
            "bias_detection.jl",           # Uses ParameterTuningSearchParameters from types.jl
            "boundary_sampling.jl",        # Uses MassErrorModel
            "utils.jl"                     # Uses all types - NOTE: MS2CHROM dependency temporarily commented out
        ]
    )
    
    # Include remaining SearchMethods files (excluding old FileReferences and FileOperations)
    # Skip ParameterTuningSearch since we loaded it above
    for (root, dirs, files) in walkdir(search_methods_dir)
        # Skip ParameterTuningSearch directory
        if occursin("ParameterTuningSearch", root)
            continue
        end
        for file in files
            if endswith(file, ".jl")
                safe_include!(joinpath(root, file))
            end
        end
    end
    
    safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "WriteOutputs"))

    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "LibrarySearch.jl"))




    # BuildSpecLib
    root_path = joinpath(package_root, "src", "Routines", "BuildSpecLib")
    
    # Include KoinaStructs directory
    safe_include_directory!(joinpath(package_root, "src", "structs", "KoinaStructs"))
    
    # Include FileOperations (only if not already loaded)
    fileops_path = joinpath(package_root, "src", "utils", "FileOperations", "FileOperations.jl")
    safe_include!(fileops_path)
    
    # FASTA processing
    safe_include!(joinpath(root_path, "structs", "mods.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_parser.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_digest.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_utils.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_protein_table.jl"))
    
    # Fragment handling
    safe_include!(joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_parse.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_index.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_annotation.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_predict.jl"))
    
    # Koina integration
    safe_include!(joinpath(root_path, "koina", "koina_api.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_prep.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_parse.jl"))
    
    # Utilities
    safe_include!(joinpath(root_path, "utils", "io.jl"))
    safe_include!(joinpath(root_path, "utils", "estimate_collision_ev.jl"))
    safe_include!(joinpath(root_path, "utils", "math.jl"))
    safe_include!(joinpath(root_path, "utils", "get_mz.jl"))
    safe_include!(joinpath(root_path, "utils", "parse_isotope_mods.jl"))
    safe_include!(joinpath(root_path, "utils", "check_params.jl"))
    
    # Structs 
    safe_include!(joinpath(root_path, "structs", "EmpiricalLibrary.jl"))
    safe_include!(joinpath(root_path, "utils", "parse_mods.jl"))
    
    # Library building
    safe_include!(joinpath(root_path, "build", "build_poin_lib.jl"))
    
    # Chronologer Methods
    safe_include!(joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_parse.jl"))
    
    # Profiling
    safe_include!(joinpath(package_root, "src", "utils", "profile.jl"))
    safe_include!(joinpath(package_root, "src", "utils", "pdfUtils.jl"))

    #importSpecLibScripts()

    return files_loaded
end

function importSpecLibScripts()
    include("types.jl")
    include("constants.jl")
    include("../../structs/Ion.jl")
    include("../../structs/LibraryIon.jl")
    include("../../structs/LibraryFragmentIndex.jl")
    # FASTA processing
    include("fasta/fasta_types.jl")
    include("fasta/fasta_parser.jl") 
    include("fasta/fasta_digest.jl")
    include("fasta/fasta_utils.jl")
    # Fragment handling
    include("fragments/fragment_types.jl")
    include("fragments/get_frag_bounds.jl")
    include("fragments/fragment_parse.jl")
    include("fragments/fragment_index.jl")
    include("fragments/fragment_annotation.jl")
    include("fragments/fragment_predict.jl")

    # Koina integration
    include("koina/koina_types.jl")
    include("koina/koina_api.jl")
    include("koina/koina_batch_prep.jl")
    include("koina/koina_batch_parse.jl")

    # Utilities
    include("utils/io.jl")
    include("utils/estimate_collision_ev.jl")
    include("utils/math.jl")
    include("utils/modifications.jl")
    include("utils/get_mz.jl")
    include("utils/parse_isotope_mods.jl")
    include("utils/check_params.jl")

    # Library building
    include("build/build_types.jl")
    include("build/build_lib.jl")
    include("build/build_index.jl")
    include("build/build_poin_lib.jl")


    include("chronologer/chronologer_types.jl")
    include("chronologer/chronologer_prep.jl")
    include("chronologer/chronologer_predict.jl")
    include("chronologer/chronologer_parse.jl")
end