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
    """
        reverseVarModKeys(var_mod_key::Vector{@NamedTuple{position::UInt8, aa::Char, mod_name::String}}, seq_length::Int)

    Reverses the positions of variable modifications in a peptide sequence according to specific rules:

    - Modifications at the last position (`seq_length`) remain in their original position
    - Modifications where amino acid (`aa`) is 'n' or 'c' (terminal modifications) remain in their original position
    - All other modifications have their positions reversed relative to the sequence length

    # Arguments
    - `var_mod_key::Vector{@NamedTuple{position::UInt8, aa::Char, mod_name::String}}`: Vector of modification specifications
    - `seq_length::Int`: Length of the peptide sequence

    # Returns
    - A new sorted vector of modification specifications with appropriate positions reversed

    # Example
    ```julia
    var_mod_key = [
        (position=UInt8(2), aa='S', mod_name="Phospho"),
        (position=UInt8(5), aa='K', mod_name="Acetyl"),
        (position=UInt8(1), aa='n', mod_name="N-term"),
        (position=UInt8(10), aa='Y', mod_name="Nitration")
    ]
    reversed = reverseVarModKeys(var_mod_key, 10)
    # Result: N-term mod stays at position 1, 
    # position 2 becomes position 9,
    # position 5 becomes position 6,
    # position 10 stays at position 10
    ```
    """
    function reverseVarModKeys(
        var_mod_key::Vector{@NamedTuple{position::UInt8, aa::Char, mod_name::String}},
        seq_length::Int
    )
        # Separate mods into those that need to stay in place and those to be reversed
        static_mods = filter(mod -> 
            mod.position == seq_length || lowercase(string(mod.aa)) in ["n", "c"], 
            var_mod_key)
        
        to_reverse_mods = filter(mod -> 
            mod.position != seq_length && !(lowercase(string(mod.aa)) in ["n", "c"]), 
            var_mod_key)
        
        # For mods to be reversed, calculate their new positions
        # New position = seq_length + 1 - old_position
        reversed_mods = map(to_reverse_mods) do mod
            new_position = UInt8(seq_length - mod.position)
            return (position=new_position, aa=mod.aa, mod_name=mod.mod_name)
        end
        
        # Combine the static and reversed mods, then sort by position
        result = vcat(static_mods, reversed_mods)
        sort!(result, by = mod -> mod.position)
        
        return result
    end

    function fillVarModStrings!(
                                var_mod_strings::Vector{String},
                                var_mod_keys::Vector{String},
                                var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                                fixed_mods_string::String,
                                max_var_mods::Int,
                                seq_length::Int,
                                seq_is_decoy::Bool
                                )
        i = 1
        for N in 1:min(max_var_mods, length(var_mod_matches))
            for mods in combinations(var_mod_matches, N)
                mods_string = fixed_mods_string
                #Need a unique key that can map between the target and reverse decoy verion of a sequence. 
                var_mod_tuples = Vector{@NamedTuple{position::UInt8, aa::Char, mod_name::String}}()
                for mod in mods
                    index = UInt8(mod[:regex_match].offset)
                    aa = mod[:regex_match].match
                    mods_string *= "("*string(index)*","*aa*","*mod[:name]*")"
                    push!(var_mod_tuples, (position=index, aa=aa[1], mod_name=mod[:name]))
                end
                if seq_is_decoy
                    var_mod_tuples = reverseVarModKeys(var_mod_tuples, seq_length)
                end
                var_mod_keys[i] = join(["("*string(mod.position)*","*mod.aa*","*mod.mod_name*")" for mod in var_mod_tuples])
                var_mod_strings[i] = mods_string
                i += 1
            end
        end
        #Version with 0 variable mods
        var_mod_strings[end] = fixed_mods_string
        var_mod_keys[end] = "" 
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
    _mod_keys= Vector{String}(undef, prec_alloc_size)
    _precursor_charge = Vector{UInt8}(undef, prec_alloc_size)
    _collision_energy = Vector{Float32}(undef, prec_alloc_size )
    _decoy = Vector{Bool}(undef, prec_alloc_size)  
    _entrapment_group_id = Vector{UInt8}(undef, prec_alloc_size)
    _base_pep_id = Vector{UInt32}(undef, prec_alloc_size)
    _base_prec_id = Vector{UInt32}(undef, prec_alloc_size)
    _pair_id = Vector{UInt32}(undef, prec_alloc_size)
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
        #Simple key that will be the same for the target and decoy version of the peptide. 
        var_mod_keys = Vector{String}(undef, n_var_mod_combinations)
        #Build modification strings for all combinations of variable mods 
        fillVarModStrings!(var_mod_strings,
                            var_mod_keys,
                            var_mod_matches,
                            fixed_mods_string,
                            max_var_mods,
                            length(sequence),
                            is_decoy(peptide) )

        #Get unique combinations of variable mods from 0-max_var_mods
        accession_id = get_id(peptide)
        decoy = is_decoy(peptide) 
        entrapment_group_id = get_entrapment_group_id(peptide)
        base_pep_id = get_base_pep_id(peptide)
        base_prec_id = get_base_prec_id(peptide)
        #if (decoy == false)
        for charge in range(min_charge, max_charge)
            for (var_mod_string, var_mod_key) in zip(var_mod_strings, var_mod_keys)
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
                _mod_keys[n] = var_mod_key
                _precursor_charge[n] = charge
                _collision_energy[n] = NCE
                _decoy[n] = decoy
                _entrapment_group_id[n] = entrapment_group_id
                _base_pep_id[n] = base_pep_id
                _base_prec_id[n] = base_prec_id
            end
        end

    end

    pair_id = UInt32(0)
    seq_to_pair_id = Dictionary{@NamedTuple{base_prec_id::UInt32, charge::UInt8, mod_key::String}, UInt32}()
    @info "Building pair IDs..."
    for i in ProgressBar(collect(range(1, n)))
        base_prec_id = _base_prec_id[i]
        charge = _precursor_charge[i]
        mod_key = _mod_keys[i]
        key = (base_prec_id=base_prec_id, charge=charge, mod_key=mod_key)
        if !haskey(seq_to_pair_id, key)
            pair_id += one(UInt32)
            insert!(seq_to_pair_id, key, pair_id)
        end
        _pair_id[i] = seq_to_pair_id[key]
    end

    
    seq_df = DataFrame(
        (upid = _upid[1:n],
         accession_number = _accession_number[1:n],
         sequence = _sequence[1:n],
         mods = _mods[1:n],
         mod_key = _mod_keys[1:n],
         precursor_charge = _precursor_charge[1:n],
         collision_energy = _collision_energy[1:n],
         decoy = _decoy[1:n],
         entrapment_group_id = _entrapment_group_id[1:n],
         base_pep_id = _base_pep_id[1:n],
         base_prec_id = _base_prec_id[1:n],
         pair_id = _pair_id[1:n])
    )

    seq_df[!,:koina_sequence] = getKoinaSeqs(
                                                    seq_df[!,:sequence],
                                                    seq_df[!,:mods],
                                                    mod_to_mass_dict
                                                    )
    
    return seq_df
end

"""
   add_pair_indexes!(df) -> DataFrame

Add a `precursor_pair_idx` column to a DataFrame that contains `pair_id` values.

For each row, the function looks for another row with the same `pair_id`:
- If found, `precursor_pair_idx` contains the row index of the paired entry
- If the `pair_id` is unique (no partner found), `precursor_pair_idx` is set to `missing`

# Arguments
- `df`: DataFrame containing a `pair_id` column where each value appears at most twice

# Returns
- The input DataFrame with an added `precursor_pair_idx` column of type `Union{Int, Missing}`

# Example
```julia
seq_df = DataFrame(
   pair_id = UInt32[1, 2, 2, 3]
)
add_pair_indices!(seq_df)
# Result: seq_df now contains a precursor_pair_idx column with values [missing, 3, 2, missing]
Note: This function modifies the input DataFrame in-place.
"""

function add_pair_indices!(df)
    n = nrow(df)
    
    # Create a dictionary mapping pair_ids to their row indices
    pair_id_to_rows = Dict{UInt32, Vector{Int}}()
    for i in 1:n
        pair_id = df.pair_id[i]
        if !haskey(pair_id_to_rows, pair_id)
            pair_id_to_rows[pair_id] = [i]
        else
            push!(pair_id_to_rows[pair_id], i)
        end
    end
    
    # Create the precursor_pair_idx column
    precursor_pair_idx = Vector{Union{Int, Missing}}(missing, n)
    
    for i in 1:n
        pair_id = df.pair_id[i]
        rows = pair_id_to_rows[pair_id]
        
        # If there are exactly 2 rows with this pair_id
        if length(rows) == 2
            # Set the index to point to the other row
            other_idx = (rows[1] == i) ? rows[2] : rows[1]
            precursor_pair_idx[i] = other_idx
        end
        # If only 1 row with this pair_id, leave as missing
    end
    
    # Add the column to the DataFrame
    df.precursor_pair_idx = precursor_pair_idx
    
    return nothing
end
