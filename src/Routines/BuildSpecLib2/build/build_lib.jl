"""
Main function to build a spectral library
"""
function build_library(params_path::String)
    # Read and validate parameters
    params = validate_build_params(read(params_path, String))
    
    # Set up output directories
    lib_dir = setup_output_directories(params)
    
    # Get fragment m/z bounds
    frag_bounds, prec_mz_min, prec_mz_max = get_fragment_bounds(
        params.library_params["auto_detect_frag_bounds"],
        params.library_params["calibration_raw_file"],
        (Float32(params.library_params["frag_mz_min"]), 
         Float32(params.library_params["frag_mz_max"])),
        (Float32(params.library_params["prec_mz_min"]), 
         Float32(params.library_params["prec_mz_max"]))
    )
    
    # Process FASTA files and predict spectra if needed
    if params.predict_fragments
        process_fasta_files(params, lib_dir, frag_bounds, prec_mz_min, prec_mz_max)
    end
    
    # Build fragment index
    build_fragment_indices(lib_dir, 
                         params.library_params["frag_bin_tol_ppm"],
                         params.library_params["rt_bin_tol"])
                         
    return lib_dir
end

"""
Setup output directories for library building
"""
function setup_output_directories(params::LibraryBuildParams)
    # Create main output directory
    mkpath(params.out_dir)
    
    # Create library directory
    lib_dir = joinpath(params.out_dir, params.lib_name * ".poin")
    mkpath(lib_dir)
    
    # Create temporary directory for chronologer
    chronologer_dir = joinpath(lib_dir, "chronologer_temp")
    mkpath(chronologer_dir)
    
    return lib_dir
end

