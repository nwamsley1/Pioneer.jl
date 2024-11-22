using Combinatorics, Dictionaries

function combineSharedPeptides(peptides::Vector{FastaEntry})
    seq_to_fasta_entry = Dictionary{String, FastaEntry}()
    n = 0
    a = 0
    for peptide in peptides
        sequence = getSeq(peptide)
        if haskey(seq_to_fasta_entry, sequence)
            a += 1
            fasta_entry = seq_to_fasta_entry[sequence]
            accession = getID(peptide)*";"*getID(fasta_entry)
            seq_to_fasta_entry[sequence] = FastaEntry(accession, 
                                                        getDescription(fasta_entry), 
                                                        getProteome(fasta_entry),
                                                        getSeq(fasta_entry),
                                                        getBasePepId(fasta_entry),
                                                        getEntrapmentGroupId(fasta_entry), 
                                                        isDecoy(fasta_entry)
                                                        )
        else
            n += 1
            insert!(seq_to_fasta_entry, sequence, peptide)
        end
    end
    fasta_entries = Vector{FastaEntry}(undef, length(seq_to_fasta_entry))
    i = 1
    for (key, value) in pairs(seq_to_fasta_entry)
        fasta_entries[i] = value
        i += 1
    end

    return fasta_entries
end

function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer, charge_facs::Vector{T}) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
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

function getChronologerSeqs(
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

function buildUniSpecInput(fasta_peptides::Vector{FastaEntry};
                           min_charge::Int = 2, max_charge::Int = 4,
                           fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(), 
                           var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = [(p=r"M", r="Unimod:4")],
                           mod_to_mass_dict::Dict{String, String} = Dict("Unimod:4" => "16.000"),
                           max_var_mods::Int = 2,
                           nce::Float64 = 30.0,
                           default_charge::Int = 3,
                           dynamic_nce::Bool = true
                           )

    charge_facs = Float64[1, 0.9, 0.85, 0.8, 0.75]

    #Number of precursors to pre-allocate memory for
    prec_alloc_size = 0
    for peptide in fasta_peptides
        sequence = getSeq(peptide)
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
        sequence = getSeq(peptide)

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
        accession_id = getID(peptide)
        is_decoy = isDecoy(peptide) 
        entrapment_group_id = getEntrapmentGroupId(peptide)
        base_pep_id = getBasePepId(peptide)
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
                _upid[n] = getProteome(peptide)
                _accession_number[n] = accession_id
                _sequence[n] = sequence
                _mods[n] = var_mod_string
                _precursor_charge[n] = charge
                _collision_energy[n] = NCE
                _decoy[n] = is_decoy
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

    seq_df[!,:chronologer_sequence] = getChronologerSeqs(
                                                    seq_df[!,:sequence],
                                                    seq_df[!,:mods],
                                                    mod_to_mass_dict
                                                    )

    return seq_df
end

