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
            @user_warn "File not found: $file_path"
            return false
        end
        
        # Include the file
        try
            include(file_path)
            push!(files_loaded, file_path)
            successful_includes += 1
            return true
        catch e
            @user_warn "Failed to include $file_path: $e"
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
            "binIsotopeRatioData.jl",
            "SquareQuadModel.jl"
        ]
    )

    # Load SearchDIA early to make asset_path available to other modules
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA.jl"))
    
    safe_include!(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "paramDefaults.jl"))
    safe_include!(joinpath(package_root, "src", "Routines","BuildSpecLib", "utils", "buildParamDefaults.jl"))
    safe_include!(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))
    safe_include!(joinpath(package_root, "src", "Routines","BuildSpecLib", "structs", "mods.jl"))
    
    # Logging is now handled directly in Pioneer.jl
    
    include_files!(
        joinpath(package_root, "src","structs"),
        [
            "MassSpecData/types.jl",
            "MassSpecData/getters.jl",
            "MassSpecData/FilteredMassSpecData.jl",
            "MassSpecData/IndexedMassSpecData.jl",
            "ChromObject.jl",
            "Counter.jl",
            "Ion.jl",
            "SpectralLibrary/fragment_types.jl",
            "SpectralLibrary/fragment_lookup.jl",
            "SpectralLibrary/precursors.jl",
            "SpectralLibrary/proteins.jl",
            "SpectralLibrary/PartitionedFragmentIndex/types.jl",
            "SpectralLibrary/SpectralLibrary.jl",
            "IsotopeTraceType.jl",
            "MatchIon.jl",
            "SparseArray.jl",
            "SparseArrayFused.jl",
            "PrecEstimation.jl",
            "PrecursorMap.jl",
            "FragBoundModel.jl",
            "RetentionTimeIndex.jl",
            "MassErrorModel.jl",
            "RetentionTimeConversionModel.jl",
            "IntensityMassErrorModel.jl",
            "protein_inference_types.jl"
            ]
    )

    
    # Sort utilities (needed by ML and FileOperations)
    safe_include!(joinpath(package_root, "src", "utils", "sortUtils.jl"))

    # Parallel utilities (used by threading patterns across SearchDIA)
    safe_include!(joinpath(package_root, "src", "utils", "parallelUtils.jl"))

    # Utilities/SpectralDeconvolution - OLS and Poisson MM solvers
    include_files!(
        joinpath(package_root, "src", "utils", "SpectralDeconvolution"),
        [
            "solveOLS.jl",          # Non-negative OLS coordinate descent
            "solvePoissonMM.jl",    # Poisson MLE coordinate descent + dispatch
            "solveHuber.jl",        # Robust Huber coordinate descent
        ]
    )

    # Utilities/Splines - B-spline fitting and evaluation
    include_files!(
        joinpath(package_root, "src", "utils", "Splines"),
        [
            "uniformBasisCubicSpline.jl",  # RT alignment, quant normalization
            "libraryBSpline.jl",           # NCE intensity evaluation (de Boor)
        ]
    )

    # Utilities/ML - general ML utilities
    include_files!(
        joinpath(package_root, "src", "utils", "ML"),
        [
            "fdrUtilities.jl",
            "ftrUtilities.jl",
            "probitRegression.jl",
            "piecewiseLinearFunction.jl",
            "wittakerHendersonSmoothing.jl",
        ]
    )

    # Utilities/ML/PSMScoring - trait-based PSM scoring pipeline (in dependency order)
    include_files!(
        joinpath(package_root, "src", "utils", "ML", "PSMScoring"),
        [
            "psm_container.jl",       # AbstractPSMContainer abstraction
            "arrow_psm_container.jl", # ArrowFilePSMContainer (file-backed OOM)
            "types.jl",               # Abstract types (PSMScoringModel, etc.)
            "config.jl",              # ScoringConfig struct
            "lightgbm_utils.jl",      # LightGBM API wrapper
            "model_training.jl",      # PSMScoringModel train/predict dispatch
            "training_selection.jl",  # TrainingDataStrategy implementations
            "feature_selection.jl",   # FeatureSelectionStrategy implementations
            "workspace.jl",           # AbstractScoringWorkspace, CV fold setup
            "scoring.jl",             # percolator_scoring! entry point
        ]
    )


    # Utils
    include_files!(
        joinpath(package_root, "src", "utils"),
        [
            "serialization.jl",
            "isotopes.jl",
            "isotopeSplines.jl",
            "maxLFQ.jl",
            "normalizeQuant.jl",
            "proteinInference.jl",
            "profile.jl",
            "pdfUtils.jl"
        ]
    )

    # FileOperations module (includes writeArrow, safeFileOps, Arrow I/O, streaming, pipeline)
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
    # fusedScan defines FusedScratch, which SearchTypes references as a field
    # on SimpleLibrarySearch — must load first.
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "fusedScan.jl"))
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))

    # Include remaining files using safe import for directories
    safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils"))
    safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "ParseInputs"))

    # Partitioned fragment index (build + search depend on CommonSearchUtils types)
    include_files!(
        joinpath(package_root, "src", "structs", "SpectralLibrary", "PartitionedFragmentIndex"),
        ["build.jl", "search.jl"]
    )
    
    # SearchMethods (excluding the old FileReferences.jl and FileOperations.jl files)
    search_methods_dir = joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods")
    
    # Load ParameterTuningSearch files in dependency order
    include_files!(
        joinpath(search_methods_dir, "ParameterTuningSearch"),
        [
            "constants.jl",                # Hardcoded tuning constants (shared by NCE/Quad tuning too)
            "types.jl",                    # All type definitions to avoid circular dependencies
            "ParameterTuningSearch.jl",    # Main implementation (types moved to types.jl)
            "utils.jl",                    # Uses all types - NOTE: MS2CHROM dependency temporarily commented out
            "fit_intensity_mass_error.jl", # IntensityMassErrorModel fitting pipeline
            "ms1_diagnostic.jl"            # MS1 ppm-residual diagnostic histogram
        ]
    )
    
    # Load PrecursorScoringSearch files in dependency order
    include_files!(
        joinpath(search_methods_dir, "PrecursorScoringSearch"),
        [
            "utils.jl",                        # get_qvalue_spline + other helpers
            "model_config.jl",                 # Model configuration
            "mbr_pairing.jl",                  # MBR Phase 1: 1:1 pair regeneration with cloning
            "mbr_streaming.jl",                # MBR Phase 2: donor dict + per-file MBR sidecars
            "mbr_ftr.jl",                      # MBR Phase 4: FTR controller on MBR-boosted score
            "pass1_oom.jl",                    # Out-of-memory Pass-1 training (stream + reservoir sample + per-file predict)
            "score_psms.jl",                   # PSM scoring functions
            "wide_window_features.jl",          # Raw same-window MS1/MS2 support features
            "scoring_interface.jl",            # Interface functions
            "build_rt_indices.jl",             # RT index construction for IntegrateChromatogramsSearch
            "PrecursorScoringSearch.jl"        # Main implementation - depends on utils.jl
        ]
    )

    # ProteinInferenceSearch (annotates passing PSMs with inferred protein groups)
    include_files!(
        joinpath(search_methods_dir, "ProteinInferenceSearch"),
        [
            "utils.jl",
            "ProteinInferenceSearch.jl"
        ]
    )

    # ProteinScoringSearch — utils.jl cascades the other 7 sibling files
    include_files!(
        joinpath(search_methods_dir, "ProteinScoringSearch"),
        [
            "utils.jl",
            "ProteinScoringSearch.jl"
        ]
    )
    
    # MainSearch (fragment index search + deconvolution + prescore scoring)
    include_files!(joinpath(search_methods_dir, "MainSearch"), [
        "prescore_aggregation.jl",       # _logodds_combine (no deps)
        "types.jl",                      # MainSearchParameters (no deps)
        "deconvolution.jl",              # deconvolve_spectra, deconvolve_scans! (thin wrapper)
        "features.jl",                   # prepare_psm_features!, add_features! (uses types)
        "irt_refinement.jl",             # predicted iRT refinement between LGBM passes
        "scoring.jl",                    # train_lgbm_and_select_best (uses features)
        "utils.jl",                      # recalibrate_rt!
        "MainSearch.jl"                  # struct + interface (uses everything above)
    ])

    # Unified scan loop (needs MainSearchParameters + CommonSearchUtils + PSM scoring)
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "process_scans.jl"))
    # Fused variant for MainSearch (loaded after process_scans.jl — uses its dispatch helpers)
    safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "process_scans_fused.jl"))

    # Chromatogram integration (explicit order so new files are precompile-tracked)
    include_files!(joinpath(search_methods_dir, "IntegrateChromatogramsSearch"), [
        "integrate_chrom.jl",
        "IntegrateChromatogramsSearch.jl",
        "utils.jl"
    ])

    # Include remaining SearchMethods files (excluding explicitly loaded directories)
    for (root, dirs, files) in walkdir(search_methods_dir)
        root_basename = basename(root)
        if root_basename == "ParameterTuningSearch" || root_basename == "PrecursorScoringSearch" ||
           root_basename == "HuberTuningSearch" ||
           root_basename == "ProteinInferenceSearch" || root_basename == "ProteinScoringSearch" ||
           occursin("MainSearch", root)
            continue
        end
        for file in files
            if endswith(file, ".jl")
                safe_include!(joinpath(root, file))
            end
        end
    end

    # Huber calibration depends on the chromatogram-integration fused RT-window
    # utilities, so load it after the remaining SearchMethods pass above.
    include_files!(
        joinpath(search_methods_dir, "HuberTuningSearch"),
        [
            "HuberTuningSearch.jl",
            "utils.jl"
        ]
    )
    
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
    safe_include!(joinpath(root_path, "fasta", "diann_decoys.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_protein_table.jl"))
    
    # Fragment handling
    safe_include!(joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    safe_include!(joinpath(root_path, "fragments", "site_determining.jl"))
    safe_include!(joinpath(root_path, "fragments", "isomer_decoys.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_parse.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_annotation.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_predict.jl"))

    # Koina integration
    safe_include!(joinpath(root_path, "koina", "koina_client.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_api.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_prep.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_parse.jl"))

    # Utilities
    safe_include!(joinpath(root_path, "utils", "io.jl"))
    safe_include!(joinpath(root_path, "utils", "math.jl"))
    safe_include!(joinpath(root_path, "utils", "get_mz.jl"))
    safe_include!(joinpath(root_path, "utils", "check_params.jl"))
    safe_include!(joinpath(root_path, "utils", "essential_mods.jl"))
    safe_include!(joinpath(root_path, "utils", "parse_mods.jl"))

    # Library building
    safe_include!(joinpath(root_path, "build", "build_poin_lib.jl"))

    # Chronologer Methods
    safe_include!(joinpath(root_path, "chronologer", "pair_decoys.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_parse.jl"))

    # Profiling
    safe_include!(joinpath(package_root, "src", "utils", "profile.jl"))
    safe_include!(joinpath(package_root, "src", "utils", "pdfUtils.jl"))

    # Main routines that use logging macros - load at the end
    # SearchDIA.jl already loaded early to provide asset_path function
    safe_include!(joinpath(package_root, "src", "Routines", "BuildSpecLib.jl"))
    safe_include!(joinpath(package_root, "src", "Routines", "GenerateParams.jl"))
    safe_include!(joinpath(package_root, "src", "Routines", "mzmlConverter", "convertMzML.jl"))

    return files_loaded
end
