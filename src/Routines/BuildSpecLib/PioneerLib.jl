function BuildSpecLib(params_path::String)
    #Read library building parameters json
    # params_path = "/Users/n.t.wamsley/Projects/Pioneer.jl/data/example_config/fragIndex111914.json"
    params = checkParams(read(params_path, String));
    println("chronologer_dir ", joinpath(@__DIR__, "../../../chronologer/Predict_RT.py"))
    #Directory in which to place library
    lib_out_dir = params["out_dir"]
    if !isdir(lib_out_dir)
        mkpath(lib_out_dir)
    end
    #Library is a folder with ".pion" extention 
    lib_dir = joinpath(lib_out_dir, params["lib_name"]*".poin")
    if !isdir(lib_dir)
        mkpath(lib_dir)
    end
    #Temporary folder in the .pion library to hold the chronologer results 
    chronologer_dir = joinpath(lib_dir, "chronologer_temp")
    if !isdir(chronologer_dir)
        mkpath(chronologer_dir)
    end
    chronologer_out_path = joinpath(chronologer_dir, "precursors_for_chronologer.tsv")

    #Reformat key parameters
    _params = (
        fasta_digest_params = Dict{String, Any}(k => v for (k, v) in params["fasta_digest_params"]),
        nce_params = Dict{String, Any}(k => v for (k, v) in params["nce_params"]),
        library_params = Dict{String, Any}(k => v for (k, v) in params["library_params"]),
    )
    ###########
    #Parse precursor and fragment m/z ramges 
    frag_bounds, prec_mz_min, prec_mz_max = getFragBounds(
        _params[:library_params]["auto_detect_frag_bounds"],
        _params[:library_params]["calibration_raw_file"],
        (Float32(_params[:library_params]["frag_mz_min"]), Float32(_params[:library_params]["frag_mz_max"])),
        (Float32(_params[:library_params]["prec_mz_min"]), Float32(_params[:library_params]["prec_mz_max"]))
    )
    #Test if the library fragments have already been predicted.
    #If they have not been, digest the fasta and predict retention times 
    #and spectra for each precursor 
    if params["predict_fragments"]
        frag_annotation_type = missing
        prediction_model = params["library_params"]["prediction_model"]
        if !(prediction_model ∈ prediction_model_options)
            error("User supplied prediction model `$prediction_model` not in available options: 
            $prediciton_model_options")
        else
            println("Using prediction model $prediction_model")
        end
        #msp convetion is more complex. Also denotes charges by `^`
        #prosit convention is simple. `y1+1`
        frag_annotation_type = prediction_model_to_annotation_type[prediction_model]
        #Determines how collision energy is handled. 
        koina_model_type = prediction_model_to_model_type[prediction_model]
        #Digest fasta files to get a table of sequences with mods, charge states, 
        #colliion energies
        mz_to_ev_interp = missing
        try
            #Unispec can take raw ev values rather than the normalized collision energies. 
            #If there is an example raw file, this should be more accurate than using NCE
            if (occursin("unispec", params["library_params"]["prediction_model"]))
                mz_to_ev_interp = getMzToEvInterp(
                    _params[:library_params]["calibration_raw_file"],
                    lib_dir
                )
            end
        catch
            default_nce = params["nce_params"]["nce"]
            default_charge = params["nce_params"]["default_charge"]
            @warn "Could not estimate mz to ev conversion from raw file...
            Using default nce of $default_nce at default charge of $default_charge
            and assuming Thermo instrument"
        end
        instrument_type = missing
        if ismissing(mz_to_ev_interp)
            instrument_type = _params[:library_params]["instrument_type"]
        else
            instrument_type = "NONE"
        end
        println("Parsing fasta files and preparing peptide sequences...")
        #Digest fasta into precursors and build a dataframe with one row per precursor
        prepareChronologerInput(params,
                                mz_to_ev_interp,
                                prec_mz_min,
                                prec_mz_max,
                                chronologer_out_path)

        #Predict the retention times for each precursor in the `fasta_df` and 
        #write to a .csv file
        println("Predicting Retention Times with chronologer (Wilburn et al. 2023)...")
        chronologer_table = Arrow.Table(chronologer_out_path)
        #run(`python3.9 ../chronologer/Predict_RT.py $chronologer_out_path $chronologer_out_path`)
        chronologer_dir = joinpath(@__DIR__, "../../../chronologer/Predict_RT.py")
        run(`python3.9 $chronologer_dir $chronologer_out_path $chronologer_out_path`)

        #Parse isotope mod groups into this dictionary
        iso_mod_to_mass = Dict{String, Float32}()
        #Convert to an arrow table and write to the main library folder 
        precursors_arrow_path = parseChronologerOutput(chronologer_out_path,
                                                lib_dir,
                                                Dict{String,Int8}(), #mods_to_sulfur_diff
                                                iso_mod_to_mass,
                                                params["isotope_mod_groups"],
                                                Float32(_params[:library_params]["rt_bin_tol"]))
        ##########
        #Manage file paths 
        rm(chronologer_out_path, force = true)
        dir, filename = splitdir(precursors_arrow_path)
        name, _ = splitext(filename)
        raw_fragments_arrow_path = joinpath(dir, name*"_raw_fragments.arrow")
        rm(raw_fragments_arrow_path, force = true)
        ##########
        #Request fragment prediction from koina in batches 
        println("Predicting fragment ion intensities using koina server...")
        predictFragments(precursors_arrow_path,
                        raw_fragments_arrow_path,
                        koina_model_type,
                        instrument_type,
                        params["max_koina_requests"],
                        params["max_koina_batch"],
                        prediction_model)
        #=
        altimeter_dir = "/Volumes/d.goldfarb/Active/Backpack/libraries/astral/Altimeter101324_MixedSpecies_OlsenAstral_NoEntrapment_101324_bright-u10"
        predictFragments(
                        raw_fragments_arrow_path,
                        basedir)    


        frags_out_path = raw_fragments_arrow_path

        println("Getting json fpaths...")
            ordered_altimeter_json_paths = joinpath.(
                altimeter_dir,
                sort(
                    readdir(altimeter_dir), 
                    by = x->parse(Int64, String(split(split(x,'_')[end], '.')[1])) #Names in format of "Altimiter_####.json"
                    )
            )
            precursor_idx = one(UInt32)
            println("Reading Altimeter Outputs...")
            batch_dfs = DataFrame()
            batch_counter = 0

            batch_path = first(ordered_altimeter_json_paths)
            df, frags_per_prec, knot_vec = parseBatchToTable(JSON.parse(read(batch_path, String))["outputs"], SplineCoefficientModel("test"))
            filterEachPrecursor!(df, model_type, intensity_threshold = intensity_threshold)


            JSON.parse(read(batch_path, String))["outputs"][4]


            JSON.parse(read(batch_path, String))["outputs"][3]["data"][1]
            JSON.parse(read(batch_path, String))["outputs"][4]["data"][[1, 50001, 100001, 150001]]


            n = 1
            JSON.parse(read(batch_path, String))["outputs"][3]["data"][n]
            JSON.parse(read(batch_path, String))["outputs"][4]["data"][[n, 50000 + n, 100000 + n, 150000 + n]]
            n += 1

            flat_frags = JSON.parse(read(batch_path, String))["outputs"][4]["data"]

            i = 1
            test_coef = Vector{NTuple{n_coef, Float32}}(undef, 50000)
            test_idx = Vector{NTuple{n_coef, Int64}}(undef, 50000)
            @inbounds for i in range(1, 50000)
                test_coef[i] = ntuple(j -> flat_frags[i + (j - 1)*50000], 4)
                test_idx[i] = ntuple(j -> i + (j - 1)*50000, 4)
            end

                        test_coef = Vector{NTuple{n_coef, Float32}}(undef, 50000)
            test_idx = Vector{NTuple{n_coef, Int64}}(undef, 50000)
            @inbounds for i in range(1, 50000)
                test_coef[i] = ntuple(j -> flat_frags[(i-1)*4 + j], 4)
                test_idx[i] = ntuple(j -> (i-1)*4 + j, 4)
            end



                    test_coef = Vector{NTuple{n_coef, Float32}}(undef, 50000)
            test_idx = Vector{NTuple{n_coef, Int64}}(undef, 50000)
            @inbounds for i in range(1, 50000)
                test_coef[i] = ntuple(j -> flat_frags[i + (j-1)*50], 4)
                test_idx[i] = ntuple(j -> i + (j-1)*50, 4)
            end



            knots = Float32.(JSON.parse(read(batch_path, String))["outputs"][1]["data"]);
            tus = UniformSpline(
                 SVector(df[1,:coefficients]),
                 3,
                 6.0f0,
                 55.0f0,
                 7.0f0
            )
            
            
        =#            

        precursors_table = Arrow.Table(precursors_arrow_path)
        fragments_table = Arrow.Table(raw_fragments_arrow_path)
        ion_annotation_set =  getIonAnnotationSet(fragments_table[:annotation])
        frag_name_to_idx = Dict(zip(collect(ion_annotation_set), collect(range(one(UInt16), UInt16(length(ion_annotation_set))))))
        
        parseKoinaFragments(
            precursors_table,
            fragments_table,
            frag_annotation_type,
            ion_annotation_set,
            frag_name_to_idx,
            params["match_lib_build_batch"],
            joinpath(@__DIR__, "../data/immonium2.txt"),
            lib_dir,
            Dict{String, Int8}(), #mods to sulfur dict 
            iso_mod_to_mass,
        );

        fragments_table = nothing
        #fragments_table = Arrow.Table(joinpath(lib_dir,"fragment_table.arrow"));
        #prec_to_frag = Arrow.Table(joinpath(lib_dir, "prec_to_frag.arrow"));
        precursors_table = DataFrame(precursors_table)
        rename!(precursors_table, :accession_number=>:accession_numbers)
        rename!(precursors_table, :precursor_charge=>:prec_charge)
        rename!(precursors_table, :decoy=>:is_decoy)
        rename!(precursors_table, :mods=>:structural_mods)
        rename!(precursors_table, :isotope_mods=>:isotopic_mods)
        precursors_table[!,:missed_cleavages] = UInt8.(precursors_table[!,:missed_cleavages])
        Arrow.write(
            joinpath(lib_dir, "precursors_table.arrow"),
            precursors_table
        )
        precursors_table = nothing
        GC.gc()
    end

    if (!isfile(joinpath(lib_dir,"fragments_table.arrow")) | (
        !isfile(joinpath(lib_dir,"prec_to_frag.arrow"))) | (
        !isfile(joinpath(lib_dir,"precursors_table.arrow"))))

        lib_dir_files = readdir(lib_dir, join=true)
        error("Cannot build fragment index without precursors_table.arrow,
        fragments_table.arrow, and prec_to_frag.arrow. `lib_dir` $lib_dir
        only contains $lib_dir_files. Try running again with params['predict_fragments'] equal
        to true")
    end
    new_lib_dir = nothing
    if !params["predict_fragments"]

        println("Copy old tables")
        new_lib_dir = joinpath(lib_out_dir, params["new_lib_name"]*".poin")
        if !isdir(new_lib_dir)
            mkpath(new_lib_dir)
        end
        #Copy the previous predicted fragments to the new folder 
        if new_lib_dir != lib_dir
            #Copy predicted fragments and precursors to teh new directory
            println("Copy fragments table...")
            #rm(joinpath(new_lib_dir, "fragments_table.arrow"), force = true)
            old_fragments_table = Arrow.Table( joinpath(lib_dir, "fragments_table.arrow"))
            Arrow.write( joinpath(new_lib_dir, "fragments_table.arrow"),
            old_fragments_table)
            old_fragments_table = nothing

            old_prec_to_frag = Arrow.Table( joinpath(lib_dir, "prec_to_frag.arrow"))
            Arrow.write( joinpath(new_lib_dir, "prec_to_frag.arrow"),
            old_prec_to_frag)
            old_prec_to_frag = nothing


            old_precursors_table = Arrow.Table( joinpath(lib_dir, "precursors_table.arrow"))
            Arrow.write( joinpath(new_lib_dir, "precursors_table.arrow"),
            old_precursors_table)
            old_precursors_table = nothing
            GC.gc()
        end
    else
        new_lib_dir = lib_dir
    end
    rm(joinpath(new_lib_dir, "config.json"), force = true)
    cp(
        params_path,
        joinpath(new_lib_dir, "config.json")
    )
    println("build pioneer library")
    buildPionLib(
        new_lib_dir,
        UInt8(_params[:library_params]["y_start_index"]),
        UInt8(_params[:library_params]["y_start"]),
        UInt8(_params[:library_params]["b_start_index"]),
        UInt8(_params[:library_params]["b_start"]),
        _params[:library_params]["include_p_index"],
        _params[:library_params]["include_p"],
        _params[:library_params]["include_isotope"],
        _params[:library_params]["include_immonium"],
        _params[:library_params]["include_internal"],
        _params[:library_params]["include_neutral_diff"],
        UInt8(_params[:library_params]["max_frag_charge"]),
        UInt8(_params[:library_params]["max_frag_rank"]),
        Float32(_params[:library_params]["min_frag_intensity"]),
        UInt8.(_params[:library_params]["rank_to_score"]),
        frag_bounds,
        Float32(_params[:library_params]["frag_bin_tol_ppm"]),
        Float32(_params[:library_params]["rt_bin_tol"])
    )
    return nothing
end