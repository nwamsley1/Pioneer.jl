InterpolationTypeAlias = Interpolations.Extrapolation{Float32, 1,
Interpolations.GriddedInterpolation{Float32, 1, Vector{Float32}, 
Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float32}}}, 
Gridded{Linear{Throw{OnGrid}}}, Line{Nothing}}


function prepareChronologerInput(
    params::Dict{String, Any},
    mz_to_ev_interp::Union{Missing, InterpolationTypeAlias},
    prec_mz_min::Float32,
    prec_mz_max::Float32,
    chronologer_out_path::String)

    _params = (
        fasta_digest_params = Dict{String, Any}(k => v for (k, v) in params["fasta_digest_params"]),
        nce_params = Dict{String, Any}(k => v for (k, v) in params["nce_params"]),
        library_params = Dict{String, Any}(k => v for (k, v) in params["library_params"]),
    )
    ###########
    #Parse Variable and fixed mods from json
    n_var_mods = length(params["variable_mods"]["pattern"])
    n_fixed_mods = length(params["fixed_mods"]["pattern"])
    var_mods = Vector{@NamedTuple{p::Regex, r::String}}()
    mod_to_mass_dict = Dict{String, String}()
    for n in range(1, n_var_mods)
        pattern = Regex(params["variable_mods"]["pattern"][n])
        mass = Float32(params["variable_mods"]["mass"][n])
        name = params["variable_mods"]["name"][n]
        push!(var_mods,
        (
            p=pattern, 
            r=name
        ))
        mod_to_mass_dict[name] = string(mass)
    end
    fixed_mods = Vector{@NamedTuple{p::Regex, r::String}}()
    for n in range(1, n_fixed_mods)
        pattern = Regex(params["fixed_mods"]["pattern"][n])
        mass = Float32(params["fixed_mods"]["mass"][n])
        name = params["fixed_mods"]["name"][n]
        push!(fixed_mods,
        (
            p=pattern, 
            r=name
        ))
        mod_to_mass_dict[name] = string(mass)
    end
    mod_to_mass_float = Dict{String, Float64}()
    for (key, value) in pairs(mod_to_mass_dict)
        mod_to_mass_float[key] = parse(Float64, value)
    end
    #Why the "../../" ?
    fasta_paths = params["fasta_paths"]
    fasta_names = params["fasta_names"]
    #Get each sequence in the fasta files as a `FastaEntry` and digest according
    #to the cleavage regex and number of missed cleavages 
    fasta_entries = Vector{FastaEntry}()
    for (proteome_name, fasta) in zip(fasta_names, fasta_paths)
    append!(fasta_entries, digestFasta(
            parseFasta(fasta, proteome_name),
            proteome_name,
            #add_decoys = _params[:fasta_digest_params]["add_decoys"],
            regex = Regex(_params[:fasta_digest_params]["cleavage_regex"]),
            max_length = _params[:fasta_digest_params]["max_length"],
            min_length = _params[:fasta_digest_params]["min_length"],
            missed_cleavages = _params[:fasta_digest_params]["missed_cleavages"])
            )
    end
    fasta_entries = addEntrapmentSequences(
        fasta_entries, 
        UInt8(_params[:fasta_digest_params]["entrapment_r"])
    )
    fasta_entries = addReverseDecoys(
        fasta_entries
    )
    #combine "duplicate/shared peptides". That is,
    #peptpdies shared by one or more proteins should be combined into a single entry
    fasta_entries = combineSharedPeptides(fasta_entries);
    #Create a dataframe with each peptide. Apply variable mods 
    fasta_df = buildUniSpecInput(fasta_entries, 
                                    fixed_mods = fixed_mods,
                                    var_mods = var_mods,
                                    mod_to_mass_dict = mod_to_mass_dict,
                                    max_var_mods = _params[:fasta_digest_params]["max_var_mods"],
                                    nce = _params[:nce_params]["nce"],
                                    default_charge = _params[:nce_params]["default_charge"],
                                    dynamic_nce = _params[:nce_params]["dynamic_nce"],
                                    min_charge = _params[:fasta_digest_params]["min_charge"],
                                    max_charge = _params[:fasta_digest_params]["max_charge"]
                                    );
    ##########
    #Add columns and filter 
    fasta_df[!,:mz] =  getMZs(fasta_df[!,:sequence], fasta_df[!,:mods], fasta_df[!,:precursor_charge], mod_to_mass_float);
    fasta_df[!,:length] = UInt8.(length.(fasta_df[!,:sequence]))
    fasta_df[!,:missed_cleavages] =  zeros(UInt8, size(fasta_df, 1))
    cleavage_regex = Regex(_params[:fasta_digest_params]["cleavage_regex"])
    for i in range(1, size(fasta_df, 1))
        fasta_df[i,:missed_cleavages] = count(cleavage_regex, fasta_df[i,:sequence])
    end
    filter!(x->(x.mz>=prec_mz_min)&(x.mz<=prec_mz_max),fasta_df)
    #For astral. Need to be able to detect this from the arrow file?
    if !ismissing(mz_to_ev_interp)
        fasta_df[!,:collision_energy] = mz_to_ev_interp.(fasta_df[!,:mz])
    elseif occursin("prosit", params["library_params"]["prediction_model"])
        for (i, precursor_charge) in enumerate(fasta_df[!,:precursor_charge])
            fasta_df[i,:collision_energy] = adjustNCE(
                params["nce_params"]["nce"],
                params["nce_params"]["default_charge"],
                precursor_charge,
                charge_facs
            )
        end
    else
        fasta_df[!,:collision_energy] .= Float32(params["nce_params"]["nce"])
    end
    ##########
    #Write precursors out for chronologer irt prediction 
    Arrow.write(chronologer_out_path, fasta_df)
    return nothing
end