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
    # Combine shared peptides
    fasta_entries = combine_shared_peptides(fasta_entries)
    # Add entrapment sequences if specified
    fasta_entries = add_entrapment_sequences(
        fasta_entries,
        UInt8(_params.fasta_digest_params["entrapment_r"])
    )
    #Fasta entries with the fixed and variable mods added 
    fasta_entries = add_mods(
        fasta_entries, 
        fixed_mods, 
        var_mods,
        _params.fasta_digest_params["max_var_mods"])
    # Add decoy sequences
    fasta_entries = add_reverse_decoys(fasta_entries)

    fasta_entries = add_charge(
        fasta_entries,
        _params.fasta_digest_params["min_charge"],
        _params.fasta_digest_params["max_charge"]
    )
    # Build UniSpec input dataframe
    fasta_df = build_fasta_df(
        fasta_entries,
        mod_to_mass_dict = mod_to_mass_dict,
        nce = _params.nce_params["nce"],
        default_charge = _params.nce_params["default_charge"],
        dynamic_nce = _params.nce_params["dynamic_nce"]
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

function getFixedMods!(
                        fixed_mods::Vector{PeptideMod},
                        mod_matches::Base.RegexMatchIterator,
                        mod_name::String)
    for mod_match in mod_matches
        index = UInt8(mod_match.offset)
        aa = mod_match.match
        push!(fixed_mods, PeptideMod(index, first(aa), mod_name))
    end
    return nothing
end

function fillVarModStrings!(
                            all_mods::Vector{Vector{PeptideMod}},
                            var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                            fixed_mods::Vector{PeptideMod},
                            max_var_mods::Int
                            )
    i = 1
    for N in 1:min(max_var_mods, length(var_mod_matches))
        for mods in combinations(var_mod_matches, N)
            #Need a unique key that can map between the target and reverse decoy verion of a sequence. 
            peptide_all_mods = copy(fixed_mods)
            for mod in mods
                index = UInt8(mod[:regex_match].offset)
                aa = mod[:regex_match].match
                push!(peptide_all_mods, PeptideMod(index, first(aa), mod[:name]))
            end
            all_mods[i] = peptide_all_mods
            i += 1
        end
    end
    #Version with 0 variable mods
    all_mods[end] = fixed_mods
    #sort!(x->(x.position, x.mod_name, x.aa), all_mods)
    sort!(all_mods)
    return nothing
end

function add_charge(
    fasta_peptides::Vector{FastaEntry},
    min_charge::Int,
    max_charge::Int,
)
    min_charge, max_charge = UInt8(min_charge), UInt8(max_charge)
    fasta_peptides_wcharge = Vector{FastaEntry}()
    base_prec_id = one(UInt32)
    for fasta_peptide in fasta_peptides
        for charge in range(UInt8(min_charge), UInt8(max_charge))
            push!(fasta_peptides_wcharge,
                FastaEntry(
                    get_id(fasta_peptide),
                    get_description(fasta_peptide),
                    get_proteome(fasta_peptide),
                    get_sequence(fasta_peptide),
                    get_structural_mods(fasta_peptide),
                    get_isotopic_mods(fasta_peptide),
                    charge,
                    get_base_pep_id(fasta_peptide),
                    base_prec_id,
                    get_entrapment_group_id(fasta_peptide),
                    is_decoy(fasta_peptide)
                ))
            base_prec_id += one(UInt32)
        end
    end
    return fasta_peptides_wcharge
end
function add_mods(
    fasta_peptides::Vector{FastaEntry},
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
    var_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
    max_var_mods::Int)


    fasta_mods = Vector{FastaEntry}()
    for fasta_peptide in fasta_peptides
        sequence = get_sequence(fasta_peptide)

        #Apply fixed mods to the sequence 
        fixed_mods_vector = Vector{PeptideMod}()
        for mod in fixed_mod_names
            getFixedMods!(
                fixed_mods_vector,
                eachmatch(mod[:p], sequence),
                mod[:r]
            )
        end

        #Get each instance of a variable mod
        var_mod_matches = matchVarMods(sequence, var_mod_names)
        #Count number of unique variable mod combinations 
        n_var_mod_combinations = countVarModCombinations(var_mod_matches, max_var_mods)
        
        var_mods = Vector{Vector{PeptideMod}}(undef, n_var_mod_combinations)
        #Build modification strings for all combinations of variable mods 
        fillVarModStrings!(var_mods,
                            var_mod_matches,
                            fixed_mods_vector,
                            max_var_mods
                            )

        for var_mod in var_mods
            push!(fasta_mods,
                FastaEntry(
                    get_id(fasta_peptide),
                    get_description(fasta_peptide),
                    get_proteome(fasta_peptide),
                    get_sequence(fasta_peptide),
                    var_mod,
                    get_isotopic_mods(fasta_peptide),
                    get_charge(fasta_peptide),
                    get_base_pep_id(fasta_peptide),
                    zero(UInt32),
                    get_entrapment_group_id(fasta_peptide),
                    is_decoy(fasta_peptide),
                ))
        end
    end
    return fasta_mods 
end

function build_fasta_df(fasta_peptides::Vector{FastaEntry};
                           mod_to_mass_dict::Dict{String, String} = Dict("Unimod:4" => "16.000"),
                           nce::Float64 = 30.0,
                           default_charge::Int = 3,
                           dynamic_nce::Bool = true
                           )

  
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
    prec_alloc_size = length(fasta_peptides)

    #Pre-allocate columns
    _upid = Vector{String}(undef, prec_alloc_size)
    _accession_number = Vector{String}(undef, prec_alloc_size)
    _sequence = Vector{String}(undef, prec_alloc_size)
    _structural_mods = Vector{Union{String, Missing}}(undef, prec_alloc_size)
    _isotopic_mods = Vector{Union{String, Missing}}(undef, prec_alloc_size)
    _precursor_charge = Vector{UInt8}(undef, prec_alloc_size)
    _collision_energy = Vector{Float32}(undef, prec_alloc_size )
    _decoy = Vector{Bool}(undef, prec_alloc_size)  
    _entrapment_group_id = Vector{UInt8}(undef, prec_alloc_size)
    _base_pep_id = Vector{UInt32}(undef, prec_alloc_size)
    _base_prec_id = Vector{UInt32}(undef, prec_alloc_size)
    _pair_id = Vector{UInt32}(undef, prec_alloc_size)
    for (n, peptide) in enumerate(fasta_peptides)
        sequence = get_sequence(peptide)
        #Get unique combinations of variable mods from 0-max_var_mods
        accession_id = get_id(peptide)
        decoy = is_decoy(peptide) 
        entrapment_group_id = get_entrapment_group_id(peptide)
        NCE = nce
        if dynamic_nce
            NCE = adjustNCE(NCE, default_charge, get_charge(peptide), charge_facs)
        end
        _upid[n] = get_proteome(peptide)
        _accession_number[n] = accession_id
        _sequence[n] = sequence
        _structural_mods[n] = getModString(get_structural_mods(peptide))
        _isotopic_mods[n] = getModString(get_isotopic_mods(peptide))
        _precursor_charge[n] = get_charge(peptide)
        _collision_energy[n] = NCE
        _decoy[n] = decoy
        _entrapment_group_id[n] = entrapment_group_id
        _base_pep_id[n] = get_base_pep_id(peptide)
        _pair_id[n] = get_base_prec_id(peptide)
    end
    n = length(fasta_peptides)
    #=
    pair_id = UInt32(0)
    base_prec_id_to_pair_id = Dictionary{UInt32, UInt32}()
    @info "Building pair IDs..."
    for i in ProgressBar(collect(range(1, n)))
        base_prec_id = _base_prec_id[i]
        if !haskey(base_prec_id_to_pair_id , base_prec_id)
            pair_id += one(UInt32)
            insert!(base_prec_id_to_pair_id, base_prec_id, pair_id)
        end
        _pair_id[i] = base_prec_id_to_pair_id[base_prec_id]
    end
    =#
    
    seq_df = DataFrame(
        (upid = _upid[1:n],
         accession_number = _accession_number[1:n],
         sequence = _sequence[1:n],
         mods = _structural_mods[1:n],
         isotopic_mods = _isotopic_mods[1:n],
         precursor_charge = _precursor_charge[1:n],
         collision_energy = _collision_energy[1:n],
         decoy = _decoy[1:n],
         entrapment_group_id = _entrapment_group_id[1:n],
         base_pep_id = _base_pep_id[1:n],
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
    df[!,:partner_precursor_idx] = precursor_pair_idx
    
    return nothing
end
