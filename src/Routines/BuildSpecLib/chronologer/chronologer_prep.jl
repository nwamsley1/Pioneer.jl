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

"""
    adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer, charge_facs::Vector{T}) where {T<:AbstractFloat}

Adjust the normalized collision energy (NCE) based on precursor charge state. This is useful for thermo instruments when 
predicting a spectral library using Prosit. 

# Parameters
- `NCE::T`: Base normalized collision energy value
- `default_charge::Integer`: Default charge state for which the base NCE is calibrated
- `peptide_charge::Integer`: Actual charge state of the peptide
- `charge_facs::Vector{T}`: Vector of charge-specific adjustment factors

# Returns
- Adjusted NCE value scaled according to the charge state

# Notes
This function applies charge-specific adjustment factors to the base NCE value.
"""
function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer, charge_facs::Vector{T}) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
end

"""
    prepare_chronologer_input(
        params::Dict{String, Any},
        mz_to_ev_interp::Union{Missing, InterpolationTypeAlias},
        prec_mz_min::Float32,
        prec_mz_max::Float32,
        chronologer_out_path::String)

Prepare input data for retention time prediction with Chronologer.

# Parameters
- `params::Dict{String, Any}`: Dictionary containing all build parameters including:
  - `fasta_digest_params`: Settings for FASTA digestion
  - `nce_params`: Settings for collision energy calculation
  - `library_params`: General spectral library settings
  - `variable_mods`: Variable modification specifications
  - `fixed_mods`: Fixed modification specifications
  - `fasta_paths`: Paths to FASTA files
  - `fasta_names`: Names for the proteomes in FASTA files
- `mz_to_ev_interp::Union{Missing, InterpolationTypeAlias}`: Optional interpolation function for m/z to eV conversion
- `prec_mz_min::Float32`: Minimum precursor m/z to include
- `prec_mz_max::Float32`: Maximum precursor m/z to include 
- `chronologer_out_path::String`: Path to save the prepared data

# Returns
`nothing` - Writes prepared data to chronologer_out_path in Arrow format

# Notes
This function processes FASTA files according to the specified parameters, generating
peptides with proper modifications and charge states. It applies fixed and variable modifications
to the sequences, and uses a regular expression to determine the enzymatic cleavage site. 
It outputs a file suitable for input to the Chronologer retention time prediction algorithm.

Order of operations is: 
digest fasta -> 
add entrapment sequences -> 
add mods -> 
add charge states -> 
add decoys -> 
build dataframe for chronologer 
"""
function prepare_chronologer_input(
    params::Dict{String, Any},
    mz_to_ev_interp::Union{Missing, InterpolationTypeAlias},
    prec_mz_min::Float32,
    prec_mz_max::Float32,
    chronologer_out_path::String,
    proteins_out_path::String)

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

    # Build per-FASTA regex lists; empty strings map to `nothing`
    n_fastas = length(params["fasta_paths"])
    function build_regex_list(key)
        if haskey(params, key)
            [ isempty(r) ? nothing : Regex(r) for r in params[key] ]
        else
            fill(nothing, n_fastas)
        end
    end

    accession_rgxs = build_regex_list("fasta_header_regex_accessions")
    gene_rgxs = build_regex_list("fasta_header_regex_genes")
    protein_rgxs = build_regex_list("fasta_header_regex_proteins")
    organism_rgxs = build_regex_list("fasta_header_regex_organisms")

    # Process FASTA files
    fasta_entries = Vector{FastaEntry}()
    protein_entries = Vector{FastaEntry}()

    for (proteome_name, fasta, acc_rgx, gene_rgx, prot_rgx, org_rgx) in zip(
            params["fasta_names"],
            params["fasta_paths"],
            accession_rgxs,
            gene_rgxs,
            protein_rgxs,
            organism_rgxs,
        )
        parsed = parse_fasta(
            fasta,
            proteome_name;
            accession_regex = acc_rgx,
            gene_regex = gene_rgx,
            protein_regex = prot_rgx,
            organism_regex = org_rgx,
        )
        append!(protein_entries, parsed)
        append!(fasta_entries,
            digest_fasta(
                parsed,
                proteome_name,
                regex = Regex(_params.fasta_digest_params["cleavage_regex"]),
                max_length = _params.fasta_digest_params["max_length"],
                min_length = _params.fasta_digest_params["min_length"],
                missed_cleavages = _params.fasta_digest_params["missed_cleavages"]
            )
        )
    end

    protein_df = build_protein_df(protein_entries)
    Arrow.write(proteins_out_path, protein_df)

    # Step 1: Combine shared peptides (I/L equivalence)
    fasta_entries = combine_shared_peptides(fasta_entries)
    
    # Step 2: Assign base_seq_id for sequence tracking
    seq_entries_processed = assign_base_seq_ids!(fasta_entries)
    seq_count = length(unique([get_base_seq_id(e) for e in fasta_entries]))
    @user_info "Step 2 - Base Seq ID Assignment: $seq_entries_processed entries processed"
    @user_info "  Unique base_seq_ids: $seq_count"
    @user_info "  Purpose: Track base sequences through modification and entrapment variants"
    
    # Step 3: Add modifications (creates peptide variants before entrapment)
    fasta_entries = add_mods(
        fasta_entries, 
        fixed_mods, 
        var_mods,
        _params.fasta_digest_params["max_var_mods"])
    @user_info "Step 3 - Modifications added (creates peptide variants for entrapment shuffling)"

    # Step 4: Assign base_pep_id for peptide tracking (after modifications)
    pep_entries_processed = assign_base_pep_ids!(fasta_entries)
    pep_count = length(unique([get_base_pep_id(e) for e in fasta_entries]))
    @user_info "Step 4 - Base Pep ID Assignment: $pep_entries_processed entries processed"
    @user_info "  Unique base_pep_ids: $pep_count"
    @user_info "  Purpose: Track peptides (sequence + modifications) through entrapment variants"
    
    # Step 5: Add entrapment sequences (now handles modifications properly)
    fasta_entries = add_entrapment_sequences(
        fasta_entries,
        UInt8(_params.fasta_digest_params["entrapment_r"])
    )
    @user_info "Step 5 - Entrapment sequences added with shuffled modifications"
    
    # Step 6: Assign base_entrap_id values for entrapment grouping
    assign_base_entrap_ids!(fasta_entries)
    @user_info "Step 6 - base_entrap_id assigned for entrapment tracking"

    # Step 7: Add decoys (inherit base_seq_id and base_pep_id for pairing)
    if _params.fasta_digest_params["add_decoys"]
        fasta_entries = add_decoy_sequences(fasta_entries)
        @user_info "Step 7 - Decoy sequences added (inherit base_seq_id and base_pep_id)"
    end
        
    # Step 8: Add charges (creates precursor variants)
    fasta_entries = add_charge(
        fasta_entries,
        _params.fasta_digest_params["min_charge"],
        _params.fasta_digest_params["max_charge"]
    )
    @user_info "Step 8 - Charge states added (creates precursor variants)"

    # Step 9: Assign base_prec_id for precursor tracking
    prec_entries_processed = assign_base_prec_ids!(fasta_entries)
    prec_count = length(unique([get_base_prec_id(e) for e in fasta_entries]))
    @user_info "Step 9 - Base Prec ID Assignment: $prec_entries_processed entries processed"
    @user_info "  Unique base_prec_ids: $prec_count"
    @user_info "  Purpose: Track unique precursors for target-decoy pairing"
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
    println("🔍 PARTNER INDEX DEBUG:")
    nrow_before = nrow(fasta_df)
    println("   Before m/z filtering: $nrow_before precursors")
    filter!(x -> (x.mz >= prec_mz_min) & (x.mz <= prec_mz_max), fasta_df)
    nrow_after = nrow(fasta_df)
    percent_removed = round((1 - nrow_after/nrow_before) * 100, digits=1)
    println("   After m/z filtering: $nrow_after precursors ($percent_removed% removed)")
    
    # Apply charge-specific target-decoy pairing AFTER all filtering is complete
    # This ensures partner_precursor_idx values are valid row indices
    fasta_df = add_charge_specific_partner_columns!(fasta_df)
    
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

"""
    add_mods(
        fasta_peptides::Vector{FastaEntry},
        fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
        var_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
        max_var_mods::Int)

Generate peptide variants with all valid combinations of modifications.

# Parameters
- `fasta_peptides::Vector{FastaEntry}`: Input peptides
- `fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Fixed modification patterns and names
  - `:p`: Regex pattern to match modification sites
  - `:r`: Name of the modification to apply
- `var_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Variable modification patterns and names
- `max_var_mods::Int`: Maximum number of variable modifications allowed per peptide

# Returns
- `Vector{FastaEntry}`: New vector containing all peptides with all valid modification combinations

# Notes
This function combines fixed modifications (always applied) with all valid combinations
of variable modifications (subject to `max_var_mods` limit). The number of entries in
the output will typically be larger than the input as each peptide can generate
multiple modification variants.
"""
function add_mods(
    fasta_peptides::Vector{FastaEntry},
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
    var_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
    max_var_mods::Int)


    fasta_mods = Vector{FastaEntry}()
    # NOTE: base_pep_id will be preserved from original peptides (not made unique per modification)
    # base_prec_id will be incremented to make each modification variant unique for pairing
    current_base_prec_id = UInt32(0)
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
            current_base_prec_id += 1  # Unique base_prec_id for each modification variant
            # NOTE: Preserve original base_pep_id to maintain link with entrapment sequences
            push!(fasta_mods,
                FastaEntry(
                    get_id(fasta_peptide),
                    get_description(fasta_peptide),
                    get_gene(fasta_peptide),
                    get_protein(fasta_peptide),
                    get_organism(fasta_peptide),
                    get_proteome(fasta_peptide),
                    get_sequence(fasta_peptide),
                    get_start_idx(fasta_peptide),
                    var_mod,
                    get_isotopic_mods(fasta_peptide),
                    get_charge(fasta_peptide),
                    get_base_seq_id(fasta_peptide),   # preserve base_seq_id
                    get_base_entrap_id(fasta_peptide), # preserve base_entrap_id
                    get_base_pep_id(fasta_peptide),  # Preserve original base_pep_id
                    current_base_prec_id,  # Unique base_prec_id for pairing
                    get_entrapment_pair_id(fasta_peptide),
                    is_decoy(fasta_peptide),
                ))
        end
    end
    return fasta_mods 
end

"""
    add_charge(
        fasta_peptides::Vector{FastaEntry},
        min_charge::Int,
        max_charge::Int)

Create multiple precursor variants with different charge states for each peptide.

# Parameters
- `fasta_peptides::Vector{FastaEntry}`: Input peptides
- `min_charge::Int`: Minimum charge state to apply
- `max_charge::Int`: Maximum charge state to apply

# Returns
- `Vector{FastaEntry}`: New vector containing all peptides with all specified charge states

# Notes
This function expands the input peptide set by creating variants with each possible
charge state between `min_charge` and `max_charge`, inclusive. The returned vector
will have length = original_length × number_of_charge_states.
"""
function add_charge(
    fasta_peptides::Vector{FastaEntry},
    min_charge::Int,
    max_charge::Int,
)
    min_charge, max_charge = UInt8(min_charge), UInt8(max_charge)
    fasta_peptides_wcharge = Vector{FastaEntry}()
    # NOTE: Preserve base_prec_id from modification variants (don't increment per charge)
    for fasta_peptide in fasta_peptides
        for charge in range(UInt8(min_charge), UInt8(max_charge))
            push!(fasta_peptides_wcharge,
                FastaEntry(
                    get_id(fasta_peptide),
                    get_description(fasta_peptide),
                    get_gene(fasta_peptide),
                    get_protein(fasta_peptide),
                    get_organism(fasta_peptide),
                    get_proteome(fasta_peptide),
                    get_sequence(fasta_peptide),
                    get_start_idx(fasta_peptide),
                    get_structural_mods(fasta_peptide),
                    get_isotopic_mods(fasta_peptide),
                    charge,
                    get_base_seq_id(fasta_peptide),   # preserve base_seq_id
                    get_base_entrap_id(fasta_peptide), # preserve base_entrap_id
                    get_base_pep_id(fasta_peptide),   # preserve base_pep_id
                    get_base_prec_id(fasta_peptide),  # Preserve from modification variant
                    get_entrapment_pair_id(fasta_peptide),
                    is_decoy(fasta_peptide)
                ))
            # NOTE: No increment - all charge states share same base_prec_id
        end
    end
    return fasta_peptides_wcharge
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

"""
    build_fasta_df(
        fasta_peptides::Vector{FastaEntry};
        mod_to_mass_dict::Dict{String, String} = Dict("Unimod:4" => "16.000"),
        nce::Float64 = 30.0,
        default_charge::Int = 3,
        dynamic_nce::Bool = true)

Create a DataFrame containing peptide information from FASTA entries.

# Parameters
- `fasta_peptides::Vector{FastaEntry}`: Processed FASTA entries with modifications and charge states
- `mod_to_mass_dict::Dict{String, String}`: Mapping from modification names to mass shifts (as strings)
- `nce::Float64`: Base normalized collision energy value
- `default_charge::Int`: Default charge state for NCE calculation
- `dynamic_nce::Bool`: Whether to adjust NCE based on charge state

# Returns
- `DataFrame`: DataFrame containing peptide information including:
  - Precursor identification info (upid, accession_number)
  - Sequence information (sequence, mods, isotopic_mods)
  - Charge and energy (precursor_charge, collision_energy)
  - Group information (decoy, entrapment_group_id, base_pep_id, pair_id)
  - Koina-compatible sequence representation

# Notes
The returned DataFrame is formatted for compatibility
Koina-specific sequence representations that contain embedded modification
information.
"""
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
    _start_idx = Vector{UInt32}(undef, prec_alloc_size)
    _precursor_charge = Vector{UInt8}(undef, prec_alloc_size)
    _collision_energy = Vector{Float32}(undef, prec_alloc_size )
    _decoy = Vector{Bool}(undef, prec_alloc_size)  
    _entrapment_group_id = Vector{UInt8}(undef, prec_alloc_size)
    _base_seq_id = Vector{UInt32}(undef, prec_alloc_size)
    _base_entrap_id = Vector{UInt32}(undef, prec_alloc_size)
    _base_pep_id = Vector{UInt32}(undef, prec_alloc_size)
    _base_prec_id = Vector{UInt32}(undef, prec_alloc_size)
    _pair_id = Vector{UInt32}(undef, prec_alloc_size)
    for (n, peptide) in enumerate(fasta_peptides)
        sequence = get_sequence(peptide)
        #Get unique combinations of variable mods from 0-max_var_mods
        accession_id = get_id(peptide)
        decoy = is_decoy(peptide) 
        entrapment_group_id = get_entrapment_pair_id(peptide)
        NCE = nce
        if dynamic_nce
            NCE = adjustNCE(NCE, default_charge, get_charge(peptide), charge_facs)
        end
        _upid[n] = get_proteome(peptide)
        _accession_number[n] = accession_id
        _sequence[n] = sequence
        _start_idx[n] = get_start_idx(peptide)
        _structural_mods[n] = getModString(get_structural_mods(peptide))
        _isotopic_mods[n] = getModString(get_isotopic_mods(peptide))
        _precursor_charge[n] = get_charge(peptide)
        _collision_energy[n] = NCE
        _decoy[n] = decoy
        _entrapment_group_id[n] = entrapment_group_id
        _base_seq_id[n] = get_base_seq_id(peptide)
        _base_entrap_id[n] = get_base_entrap_id(peptide)
        _base_pep_id[n] = get_base_pep_id(peptide)
        _base_prec_id[n] = get_base_prec_id(peptide)
        _pair_id[n] = get_base_prec_id(peptide)  # Use base_prec_id for unique pairing per precursor variant
    end
    n = length(fasta_peptides)
    
    seq_df = DataFrame(
        (upid = _upid[1:n],
         accession_number = _accession_number[1:n],
         sequence = _sequence[1:n],
         start_idx = _start_idx[1:n],
         mods = _structural_mods[1:n],
         isotopic_mods = _isotopic_mods[1:n],
         precursor_charge = _precursor_charge[1:n],
         collision_energy = _collision_energy[1:n],
         decoy = _decoy[1:n],
         entrapment_group_id = _entrapment_group_id[1:n],
         base_seq_id = _base_seq_id[1:n],
         base_entrap_id = _base_entrap_id[1:n],
         base_pep_id = _base_pep_id[1:n],
         base_prec_id = _base_prec_id[1:n],
         pair_id = _pair_id[1:n])
    )

    seq_df[!,:koina_sequence] = getKoinaSeqs(
                                                    seq_df[!,:sequence],
                                                    seq_df[!,:mods],
                                                    mod_to_mass_dict
                                                    )
    
    # NOTE: Pairing moved to after filtering to ensure valid partner_precursor_idx values
    # seq_df = add_charge_specific_partner_columns!(seq_df)
    
    return seq_df
end

