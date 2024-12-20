# src/BuildSpecLib.jl

"""
Main function to build a spectral library from parameters
"""
function BuildSpecLib(params_path::String)
    # Read and validate parameters
    #=
    include("Routines/BuildSpecLib2/types.jl")
    include("Routines/BuildSpecLib2/build/build_types.jl")
    include("Routines/BuildSpecLib2/build/build_lib.jl")
    =#
    params_path = "/Users/n.t.wamsley/RIS_temp/koina_testing/config.json"
    params = checkParams(read(params_path, String));
    # Create output directories
    lib_out_dir = params["out_dir"]
    mkpath(lib_out_dir)

    # Library directory (.poin extension)
    lib_dir = joinpath(lib_out_dir, params["lib_name"] * ".poin")
    mkpath(lib_dir)

    # Temporary chronologer directory
    chronologer_dir = joinpath(lib_dir, "chronologer_temp")
    mkpath(chronologer_dir)
    chronologer_out_path = joinpath(chronologer_dir, "precursors_for_chronologer.arrow")

    # Format parameters into named tuple for easier access
    _params = (
        fasta_digest_params = Dict{String, Any}(k => v for (k, v) in params["fasta_digest_params"]),
        nce_params = Dict{String, Any}(k => v for (k, v) in params["nce_params"]),
        library_params = Dict{String, Any}(k => v for (k, v) in params["library_params"])
    )

    # Get fragment bounds
    frag_bounds, prec_mz_min, prec_mz_max = get_fragment_bounds(
        _params.library_params["auto_detect_frag_bounds"],
        _params.library_params["calibration_raw_file"],
        (Float32(_params.library_params["frag_mz_min"]), Float32(_params.library_params["frag_mz_max"])),
        (Float32(_params.library_params["prec_mz_min"]), Float32(_params.library_params["prec_mz_max"]))
    )

    if params["predict_fragments"]
        # Validate prediction model
        prediction_model = params["library_params"]["prediction_model"]
        if !(prediction_model ∈ PREDICTION_MODELS)
            error("Invalid prediction model: $prediction_model. Valid options: $(join(PREDICTION_MODELS, ", "))")
        end
        println("Using prediction model $prediction_model")

        # Get annotation type and model type
        frag_annotation_type = MODEL_CONFIGS[prediction_model].annotation_type
        koina_model_type = MODEL_CONFIGS[prediction_model].model_type

        # Get collision energy interpolator if needed
        mz_to_ev_interp = missing
        if occursin("unispec", prediction_model)
            try
                mz_to_ev_interp = get_mz_to_ev_interp(
                    _params.library_params["calibration_raw_file"],
                    lib_dir
                )
            catch
                @warn "Could not estimate mz to ev conversion. Using default NCE."
            end
        end

        # Set instrument type
        instrument_type = ismissing(mz_to_ev_interp) ? 
                         _params.library_params["instrument_type"] : 
                         "NONE"

        # Prepare input for retention time prediction
        println("Parsing FASTA files and preparing peptide sequences...")
        prepare_chronologer_input(params,
                                mz_to_ev_interp,
                                prec_mz_min,
                                prec_mz_max,
                                chronologer_out_path)

        # Predict retention times
        println("Predicting retention times...")
        predict_retention_times(chronologer_out_path)

        # Parse results and prepare for fragment prediction
        iso_mod_to_mass = Dict{String, Float32}()
        precursors_arrow_path = parse_chronologer_output(
            chronologer_out_path,
            lib_dir,
            Dict{String,Int8}(), # mods_to_sulfur_diff
            iso_mod_to_mass,
            params["isotope_mod_groups"],
            Float32(_params.library_params["rt_bin_tol"])
        )

        # Clean up temporary files
        rm(chronologer_out_path, force=true)
        dir, filename = splitdir(precursors_arrow_path)
        raw_fragments_arrow_path = joinpath(dir, "raw_fragments.arrow")
        rm(raw_fragments_arrow_path, force=true)

        # Predict fragments
        println("Predicting fragment ion intensities...")
        predict_fragments(
            precursors_arrow_path,
            raw_fragments_arrow_path,
            koina_model_type,
            instrument_type,
            params["max_koina_requests"],
            params["max_koina_batch"],
            prediction_model
        )

        # Load tables for processing
        precursors_table = Arrow.Table(precursors_arrow_path)
        fragments_table = Arrow.Table(raw_fragments_arrow_path)

        # Get ion annotations
        ion_annotation_set = get_ion_annotation_set(fragments_table[:annotation])
        frag_name_to_idx = Dict(ion => UInt16(i) for (i, ion) in enumerate(ion_annotation_set))

        # Parse and process fragments
        parse_koina_fragments(
            precursors_table,
            fragments_table,
            frag_annotation_type,
            ion_annotation_set,
            frag_name_to_idx,
            10000, # batch size
            joinpath(@__DIR__, "../../../data/immonium.txt"),
            lib_dir,
            Dict{String, Int8}(), # mods_to_sulfur_dict
            iso_mod_to_mass,
            koina_model_type
        )

        # Cleanup and rename tables
        fragments_table = nothing
        precursors_table = DataFrame(precursors_table)
        rename!(precursors_table, [
            :accession_number => :accession_numbers,
            :precursor_charge => :prec_charge,
            :decoy => :is_decoy,
            :mods => :structural_mods,
            :isotope_mods => :isotopic_mods
        ])

        # Convert types
        precursors_table[!, :missed_cleavages] = UInt8.(precursors_table[!, :missed_cleavages])
        precursors_table[!, :prec_charge] = UInt8.(precursors_table[!, :prec_charge])
        precursors_table[!, :mz] = Float32.(precursors_table[!, :mz])
        precursors_table[!, :irt] = Float32.(precursors_table[!, :irt])

        # Save processed precursor table
        Arrow.write(
            joinpath(lib_dir, "precursors_table.arrow"),
            precursors_table
        )

        GC.gc()
    end

    # Verify required files exist
    required_files = ["fragments_table.arrow", "prec_to_frag.arrow", "precursors_table.arrow"]
    if !all(isfile.(joinpath.(lib_dir, required_files)))
        error("Missing required files in $lib_dir. Try running with predict_fragments=true")
    end

    # Build final indices
    buildPionLib(
        lib_dir,
        UInt8(_params.library_params["y_start_index"]),
        UInt8(_params.library_params["y_start"]),
        UInt8(_params.library_params["b_start_index"]),
        UInt8(_params.library_params["b_start"]),
        _params.library_params["include_p_index"],
        _params.library_params["include_p"],
        _params.library_params["include_isotope"],
        _params.library_params["include_immonium"],
        _params.library_params["include_internal"],
        _params.library_params["include_neutral_diff"],
        UInt8(_params.library_params["max_frag_charge"]),
        UInt8(_params.library_params["max_frag_rank"]),
        Float32(_params.library_params["min_frag_intensity"]),
        UInt8.(_params.library_params["rank_to_score"]),
        frag_bounds,
        Float32(_params.library_params["frag_bin_tol_ppm"]),
        Float32(_params.library_params["rt_bin_tol"]),
    )

    return nothing
end