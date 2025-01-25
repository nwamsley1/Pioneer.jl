function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer, charge_facs::Vector{T}) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
end

"""
Prepare input data for retention time prediction with chronologer.

Parameters:
- params: Dictionary containing all build parameters
- mz_to_ev_interp: Optional interpolation function for m/z to eV conversion
- prec_mz_min: Minimum precursor m/z to include
- prec_mz_max: Maximum precursor m/z to include 
- chronologer_out_path: Path to save the prepared data

Returns:
Writes prepared data to chronologer_out_path
"""
function prepare_chronologer_input(
    params::Dict{String, Any},
    mz_to_ev_interp::Union{Missing, InterpolationTypeAlias},
    prec_mz_min::Float32,
    prec_mz_max::Float32,
    chronologer_out_path::String)

    # Parse parameters into structured format
    _params = (
        fasta_digest_params = Dict{String, Any}(k => v for (k, v) in params["fasta_digest_params"]),
        nce_params = Dict{String, Any}(k => v for (k, v) in params["nce_params"]),
        library_params = Dict{String, Any}(k => v for (k, v) in params["library_params"])
    )

    # Parse modification configurations
    var_mods = Vector{@NamedTuple{p::Regex, r::String}}()
    mod_to_mass_dict = Dict{String, String}()

    # Process variable modifications
    for i in 1:length(params["variable_mods"]["pattern"])
        pattern = Regex(params["variable_mods"]["pattern"][i])
        mass = Float32(params["variable_mods"]["mass"][i])
        name = params["variable_mods"]["name"][i]
        push!(var_mods, (p=pattern, r=name))
        mod_to_mass_dict[name] = string(mass)
    end

    # Process fixed modifications
    fixed_mods = Vector{@NamedTuple{p::Regex, r::String}}()
    for i in 1:length(params["fixed_mods"]["pattern"])
        pattern = Regex(params["fixed_mods"]["pattern"][i])
        mass = Float32(params["fixed_mods"]["mass"][i])
        name = params["fixed_mods"]["name"][i]
        push!(fixed_mods, (p=pattern, r=name))
        mod_to_mass_dict[name] = string(mass)
    end

    # Convert mass dictionary to float values
    mod_to_mass_float = Dict(k => parse(Float64, v) for (k, v) in mod_to_mass_dict)

    # Process FASTA files
    fasta_entries = Vector{FastaEntry}()
    for (proteome_name, fasta) in zip(params["fasta_names"], params["fasta_paths"])
        append!(fasta_entries, 
            digest_fasta(
                parse_fasta(fasta, proteome_name),
                proteome_name,
                regex = Regex(_params.fasta_digest_params["cleavage_regex"]),
                max_length = _params.fasta_digest_params["max_length"],
                min_length = _params.fasta_digest_params["min_length"],
                missed_cleavages = _params.fasta_digest_params["missed_cleavages"]
            )
        )
    end

    # Add entrapment sequences if specified
    fasta_entries = add_entrapment_sequences(
        fasta_entries,
        UInt8(_params.fasta_digest_params["entrapment_r"])
    )

    # Combine shared peptides
    fasta_entries = combine_shared_peptides(fasta_entries)
    
    # Add decoy sequences
    fasta_entries = add_reverse_decoys(fasta_entries)

    # Build UniSpec input dataframe
    fasta_df = add_mods_and_filter(
        fasta_entries,
        fixed_mods = fixed_mods,
        var_mods = var_mods,
        mod_to_mass_dict = mod_to_mass_dict,
        max_var_mods = _params.fasta_digest_params["max_var_mods"],
        nce = _params.nce_params["nce"],
        default_charge = _params.nce_params["default_charge"],
        dynamic_nce = _params.nce_params["dynamic_nce"],
        min_charge = _params.fasta_digest_params["min_charge"],
        max_charge = _params.fasta_digest_params["max_charge"]
    )

    # Calculate precursor m/z values
    fasta_df[!, :mz] = getMZs(
        fasta_df[!, :sequence],
        fasta_df[!, :mods],
        fasta_df[!, :precursor_charge],
        mod_to_mass_float
    )

    # Add length and missed cleavages columns
    fasta_df[!, :length] = UInt8.(length.(fasta_df[!, :sequence]))
    fasta_df[!, :missed_cleavages] = zeros(UInt8, size(fasta_df, 1))
    cleavage_regex = Regex(_params.fasta_digest_params["cleavage_regex"])
    
    for i in 1:size(fasta_df, 1)
        fasta_df[i, :missed_cleavages] = count(cleavage_regex, fasta_df[i, :sequence])
    end

    # Filter by mass range
    filter!(x -> (x.mz >= prec_mz_min) & (x.mz <= prec_mz_max), fasta_df)

    # Handle collision energy
    if !ismissing(mz_to_ev_interp)
        fasta_df[!, :collision_energy] = mz_to_ev_interp.(fasta_df[!, :mz])
    elseif occursin("prosit", params["library_params"]["prediction_model"])
        # For Prosit models, apply charge-based NCE adjustment
        for (i, precursor_charge) in enumerate(fasta_df[!, :precursor_charge])
            fasta_df[i, :collision_energy] = adjustNCE(
                params["nce_params"]["nce"],
                params["nce_params"]["default_charge"],
                precursor_charge,
                CHARGE_ADJUSTMENT_FACTORS
            )
        end
    else
        fasta_df[!, :collision_energy] .= Float32(params["nce_params"]["nce"])
    end

    # Write output
    Arrow.write(chronologer_out_path, fasta_df)

    return nothing
end



function add_mods_and_filter(fasta_peptides::Vector{FastaEntry};
                           min_charge::Int = 2, max_charge::Int = 4,
                           fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(), 
                           var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = [(p=r"M", r="Unimod:4")],
                           mod_to_mass_dict::Dict{String, String} = Dict("Unimod:4" => "16.000"),
                           max_var_mods::Int = 2,
                           nce::Float64 = 30.0,
                           default_charge::Int = 3,
                           dynamic_nce::Bool = true
                           )

                        
    function matchVarMods(sequence::String, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
        var_mod_matches = Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}()
        for mod in var_mods
            for mod_match in eachmatch(mod[:p], sequence)
                push!(var_mod_matches, (regex_match=mod_match, name=mod[:r]))
            end
        end
        return var_mod_matches
    end

    function countVarModCombinations(var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                                    max_var_mods::Int)
        n_var_mods = length(var_mod_matches)
        n_var_mod_combinations = 0
        for N in 1:min(max_var_mods, n_var_mods)
            n_var_mod_combinations += binomial(n_var_mods, N)
        end
        n_var_mod_combinations += 1 #For the sequence with no variable modifications 
        return n_var_mod_combinations
    end
                  
    function buildFixedModsString!(mod_string::String,
                            mod_matches::Base.RegexMatchIterator,
                            mod_name::String)
        for mod_match in mod_matches
            index = UInt8(mod_match.offset)
            aa = mod_match.match
            mod_string *= "("*string(index)*","*aa*","*mod_name*")"
        end
        return mod_string
    end

    function fillVarModStrings!(
                                var_mod_strings::Vector{String},
                                var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                                fixed_mods_string::String,
                                max_var_mods::Int
                                )
        i = 1
        for N in 1:min(max_var_mods, length(var_mod_matches))
            for mods in combinations(var_mod_matches, N)
                mods_string = fixed_mods_string
                for mod in mods
                    index = UInt8(mod[:regex_match].offset)
                    aa = mod[:regex_match].match
                    mods_string *= "("*string(index)*","*aa*","*mod[:name]*")"
                end
                var_mod_strings[i] = mods_string
                i += 1
            end
        end
        #Version with 0 variable mods
        var_mod_strings[end] = fixed_mods_string
    end

    function getKoinaSeqs(
        sequences::Vector{String},
        mods::Vector{Union{Missing, String}},
        mod_name_to_mass::Dict{String, String}
        )   
        mod_parse_regex = r"\(([0-9]+?),.*?,(.*?)\)"
        mod_sequences = Vector{String}(undef, length(sequences))
        for i in range(1, length(sequences))
            seq = sequences[i]
            if !ismissing(mods[i])
                #parsed_mods = eachmatch(mod_parse_regex, mods[i])
                parsed_mods = [x for x in parseMods(mods[i])]
                sort!(parsed_mods, by = x->getModIndex(x.match))
                chronologer_seq = ""
                start_idx, stop_idx = 1, length(seq)
                for mod in parsed_mods
                    stop_idx = getModIndex(mod.match)#parse(UInt8, first(mod.captures))
                    chronologer_seq *= seq[start_idx:stop_idx]
                    #chronologer_seq *= "[+"*mod_name_to_mass[getModName(mod.match)]*"]"
                    chronologer_seq *= "["*uppercase(getModName(mod.match))*"]"
                    start_idx = stop_idx += one(UInt8)
                end
                chronologer_seq *= seq[start_idx:length(seq)]
                mod_sequences[i] = chronologer_seq
            else
                mod_sequences[i] = seq
            end
    
        end
        return mod_sequences
    end
    charge_facs = Float64[1, 0.9, 0.85, 0.8, 0.75]

    #Number of precursors to pre-allocate memory for
    prec_alloc_size = 0
    for peptide in fasta_peptides
        sequence = get_sequence(peptide)
        var_mod_matches = matchVarMods(sequence, var_mods)
        n_var_mod_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        prec_alloc_size += n_var_mod_combinations*(max_charge - min_charge + 1)
    end

    #Pre-allocate columns
    _upid = Vector{String}(undef, prec_alloc_size)
    _accession_number = Vector{String}(undef, prec_alloc_size)
    _sequence = Vector{String}(undef, prec_alloc_size)
    _mods = Vector{Union{String, Missing}}(undef, prec_alloc_size)
    _precursor_charge = Vector{UInt8}(undef, prec_alloc_size)
    _collision_energy = Vector{Float32}(undef, prec_alloc_size )
    _decoy = Vector{Bool}(undef, prec_alloc_size)  
    _entrapment_group_id = Vector{UInt8}(undef, prec_alloc_size)
    _base_pep_id = Vector{UInt32}(undef, prec_alloc_size)
    n = 0
    for peptide in fasta_peptides
        sequence = get_sequence(peptide)

        #Check for illegal amino acid characters
        if (occursin("[H", sequence)) | (occursin("U", sequence)) | (occursin("O", sequence)) |  (occursin("X", sequence)) | occursin("Z", sequence) | occursin("B", sequence)
            continue
        end

        #Apply fixed mods to the sequence 
        fixed_mods_string = ""
        for mod in fixed_mods
            fixed_mods_string = buildFixedModsString!(
                fixed_mods_string,
                eachmatch(mod[:p], sequence),
                mod[:r]
            )
        end

        #Get each instance of a variable mod
        var_mod_matches = matchVarMods(sequence, var_mods)
        #Count number of unique variable mod combinations 
        n_var_mod_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        var_mod_strings = Vector{String}(undef, n_var_mod_combinations)
        #Build modification strings for all combinations of variable mods 
        fillVarModStrings!(var_mod_strings,
                            var_mod_matches,
                            fixed_mods_string,
                            max_var_mods)

        #Get unique combinations of variable mods from 0-max_var_mods
        accession_id = get_id(peptide)
        decoy = is_decoy(peptide) 
        entrapment_group_id = get_entrapment_group_id(peptide)
        base_pep_id = get_base_pep_id(peptide)
        #if (decoy == false)
        for charge in range(min_charge, max_charge)
            for var_mod_string in var_mod_strings
                n += 1
                if var_mod_string == ""
                    var_mod_string = missing
                end
                NCE = nce
                if dynamic_nce
                    NCE = adjustNCE(NCE, default_charge, charge, charge_facs)
                end
                _upid[n] = get_proteome(peptide)
                _accession_number[n] = accession_id
                _sequence[n] = sequence
                _mods[n] = var_mod_string
                _precursor_charge[n] = charge
                _collision_energy[n] = NCE
                _decoy[n] = decoy
                _entrapment_group_id[n] = entrapment_group_id
                _base_pep_id[n] = base_pep_id
            end
        end

    end

    seq_df = DataFrame(
        (upid = _upid[1:n],
         accession_number = _accession_number[1:n],
         sequence = _sequence[1:n],
         mods = _mods[1:n],
         precursor_charge = _precursor_charge[1:n],
         collision_energy = _collision_energy[1:n],
         decoy = _decoy[1:n],
         entrapment_group_id = _entrapment_group_id[1:n],
         base_pep_id = _base_pep_id[1:n])
    )
    
    seq_df[!,:koina_sequence] = getKoinaSeqs(
                                                    seq_df[!,:sequence],
                                                    seq_df[!,:mods],
                                                    mod_to_mass_dict
                                                    )
    
    return seq_df
end

