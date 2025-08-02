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
    package_root = dirname(dirname(dirname(@__DIR__)))
    
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
            "profile.jl"
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
    

    return files_loaded
end