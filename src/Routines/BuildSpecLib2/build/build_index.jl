"""
Build fragment indices for the spectral library
"""
function build_fragment_indices(lib_dir::String,
                              frag_bin_tol_ppm::Float32,
                              rt_bin_tol::Float32)
    
    # Load required tables
    fragments_table = Arrow.Table(joinpath(lib_dir, "fragments_table.arrow"))
    prec_to_frag = Arrow.Table(joinpath(lib_dir, "prec_to_frag.arrow"))
    precursors_table = Arrow.Table(joinpath(lib_dir, "precursors_table.arrow"))
    
    # Get simple fragments for indexing
    simple_frags = get_simple_fragments(
        fragments_table,
        precursors_table,
        prec_to_frag
    )
    
    # Build main fragment index
    build_fragment_index!(
        lib_dir,
        simple_frags,
        frag_bin_tol_ppm,
        rt_bin_tol,
        index_name = ""
    )
    
    # Build presearch fragment index (no RT binning)
    build_fragment_index!(
        lib_dir,
        simple_frags,
        frag_bin_tol_ppm,
        typemax(Float32),
        index_name = "presearch_"
    )
end

"""
Process FASTA files and predict spectra
"""
function process_fasta_files(params::LibraryBuildParams,
                           lib_dir::String,
                           frag_bounds::FragBoundModel,
                           prec_mz_min::Float32,
                           prec_mz_max::Float32)
    
    # Parse FASTA files
    fasta_entries = Vector{FastaEntry}()
    for (proteome_name, fasta_path) in zip(params.fasta_names, params.fasta_paths)
        append!(fasta_entries, 
            digest_fasta(
                parse_fasta(fasta_path, proteome_name),
                proteome_name,
                regex = Regex(params.fasta_params["cleavage_regex"]),
                max_length = params.fasta_params["max_length"],
                min_length = params.fasta_params["min_length"],
                missed_cleavages = params.fasta_params["missed_cleavages"]
            )
        )
    end
    
    # Add entrapment sequences if specified
    if params.fasta_params["entrapment_r"] > 0
        fasta_entries = add_entrapment_sequences(
            fasta_entries,
            UInt8(params.fasta_params["entrapment_r"])
        )
    end
    
    # Add decoys
    if params.fasta_params["add_decoys"]
        fasta_entries = add_reverse_decoys(fasta_entries)
    end
    
    # Combine shared peptides
    fasta_entries = combine_shared_peptides(fasta_entries)
    
    # Create precursor dataframe
    precursor_df = build_unispec_input(
        fasta_entries,
        fixed_mods = params.fixed_mods,
        var_mods = params.var_mods,
        mod_to_mass_dict = params.mod_to_mass_dict,
        max_var_mods = params.fasta_params["max_var_mods"],
        nce = params.nce_params["nce"],
        default_charge = params.nce_params["default_charge"],
        dynamic_nce = params.nce_params["dynamic_nce"],
        min_charge = params.fasta_params["min_charge"],
        max_charge = params.fasta_params["max_charge"]
    )
    
    # Filter by mass range
    filter!(row -> row.mz >= prec_mz_min && row.mz <= prec_mz_max, precursor_df)
    
    # Predict retention times
    predict_retention_times(precursor_df, lib_dir)
    
    # Predict fragments
    predict_fragments(precursor_df, lib_dir, params)
end