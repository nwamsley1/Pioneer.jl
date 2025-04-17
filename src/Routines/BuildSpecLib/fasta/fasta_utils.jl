struct PeptideSequenceSet
    sequences::Set{String}
    # Constructor
    function PeptideSequenceSet()
        new(Set{String}())
    end
    function PeptideSequenceSet(fasta_entries::Vector{FastaEntry})
        pss = PeptideSequenceSet()
        # Add each sequence to the set
        for entry in fasta_entries
            push!(pss, get_sequence(entry))
        end
        return pss
    end
end

getSeqSet(s::PeptideSequenceSet) = s.sequences

import Base: push!
function push!(pss::PeptideSequenceSet, seq::AbstractString)
    # Replace 'I' with 'L' in the sequence and add to the set
    push!(pss.sequences, replace(seq, 'I' => 'L'))
    return pss
end

import Base: in
function in(seq::AbstractString, pss::PeptideSequenceSet)
    # Check if the modified sequence is in the set
    return replace(seq, 'I' => 'L') ∈ getSeqSet(pss)
end


"""
Fast sequence shuffling algorithm that preserves terminal amino acids.
"""
function shuffle_fast(s::String)
    ss = sizeof(s)
    l = length(s) - 1  # Preserve last amino acid

    # Create positions vector
    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end

    # Shuffle middle positions
    p = pointer(s)
    u = Vector{UInt8}(undef, ss)
    k = 1
    for i in randperm(l)
        u[k] = unsafe_load(p, v[i])
        k += 1
    end
    u[end] = unsafe_load(p, ss)  # Keep last amino acid
    
    return String(u)
end

"""
Add entrapment sequences to a set of target peptides.
Creates shuffled decoy sequences while preserving terminal amino acids.

Parameters:
- target_fasta_entries: Vector of target peptide entries
- entrapment_r: Number of entrapment sequences per target
- max_shuffle_attempts: Maximum attempts to generate unique shuffled sequence

Returns:
Vector{FastaEntry} with original and entrapment sequences
"""
function add_entrapment_sequences(
    target_fasta_entries::Vector{FastaEntry}, 
    entrapment_r::UInt8;
    max_shuffle_attempts::Int64 = 20
)::Vector{FastaEntry}
    
    # Pre-allocate output vector
    entrapment_fasta_entries = Vector{FastaEntry}(
        undef, 
        length(target_fasta_entries) * entrapment_r
    )
    
    # Track unique sequences
    sequences_set = PeptideSequenceSet(target_fasta_entries)#Set{String}()
    #sizehint!(sequences_set, length(entrapment_fasta_entries) + length(target_fasta_entries))
    #union!(sequences_set, target_sequences)
    
    n = 1
    maximum_prec_id = zero(UInt32)
    for tfe = target_fasta_entries
        if get_base_prec_id(tfe) > maximum_prec_id
            maximum_prec_id = get_base_prec_id(tfe)
        end
    end
    base_prec_id = maximum_prec_id + one(UInt32)
    for target_entry in target_fasta_entries
        for entrapment_group_id in 1:entrapment_r
            n_shuffle_attempts = 0
            
            while n_shuffle_attempts < max_shuffle_attempts
                new_sequence = shuffle_fast(get_sequence(target_entry))
                #Make sure the entrapment sequence is unique (I and L are equivalent)
                if new_sequence ∉ sequences_set
                    entrapment_fasta_entries[n] = FastaEntry(
                        get_id(target_entry),
                        get_description(target_entry),
                        get_proteome(target_entry),
                        new_sequence,
                        missing, #structural_mods 
                        missing, #istopic_mods 
                        get_charge(target_entry),
                        get_base_pep_id(target_entry),
                        base_prec_id,
                        entrapment_group_id,
                        false
                    )
                    base_prec_id += one(UInt32)
                    n += 1
                    push!(sequences_set, new_sequence)
                    break
                end
                n_shuffle_attempts += 1
            end
            
            if n_shuffle_attempts >= max_shuffle_attempts
                @warn "Max shuffle attempts exceeded for $(get_sequence(target_entry))"
            end
        end
    end
    
    return vcat(target_fasta_entries, entrapment_fasta_entries[1:n-1])
end

"""
Enhanced version of shuffle_fast that tracks position changes using a position vector.
"""
function shuffle_fast_with_positions(s::String, positions::Vector{UInt8})
    ss = sizeof(s)
    l = length(s) - 1  # Preserve last amino acid

    # Create indices vector
    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end

    # Generate random permutation
    perm = randperm(l)
    
    # Update positions vector based on the permutation
    # Store the original positions temporarily
    temp_positions = copy(positions[1:l])
    for (new_idx, old_idx) in enumerate(perm)
        positions[new_idx] = temp_positions[old_idx]
    end
    
    # Shuffle middle positions
    p = pointer(s)
    u = Vector{UInt8}(undef, ss)
    k = 1
    for i in perm
        u[k] = unsafe_load(p, v[i])
        k += 1
    end
    u[end] = unsafe_load(p, ss)  # Keep last amino acid
    
    return String(u)
end

function adjust_mod_positions(
    mods::Union{Missing, Vector{PeptideMod}}, 
    positions::Vector{UInt8},
    seq_length::UInt8
)::Union{Missing, Vector{PeptideMod}}
    # If no modifications or missing, return as is
    if ismissing(mods) || isempty(mods)
        return mods
    end
    
    # Create new vector for adjusted modifications
    adjusted_mods = Vector{PeptideMod}(undef, length(mods))
    
    # Create a reverse mapping for efficiency
    # This maps original positions to new positions
    reverse_mapping = Vector{UInt8}(undef, seq_length)
    for new_pos in 1:seq_length
        orig_pos = positions[new_pos]
        if orig_pos <= seq_length
            reverse_mapping[orig_pos] = UInt8(new_pos)
        end
    end
    
    for (i, mod) in enumerate(mods)
        position = mod.position
        aa = mod.aa
        mod_name = mod.mod_name
        
        # Special case for N-terminal and C-terminal modifications
        if aa == 'n'
            # N-terminal modifications stay at position 1
            adjusted_mods[i] = PeptideMod(UInt8(1), 'n', mod_name)
        elseif aa == 'c'
            # C-terminal modifications stay at the end
            adjusted_mods[i] = PeptideMod(seq_length, 'c', mod_name)
        else
            # For normal residue modifications, use the reverse mapping
            if position <= seq_length
                adjusted_mods[i] = PeptideMod(reverse_mapping[position], aa, mod_name)
            else
                # Edge case handling if position is somehow out of bounds
                adjusted_mods[i] = mod
            end
        end
    end
    sort!(adjusted_mods)
    return adjusted_mods
end


"""
    addReverseDecoys(target_fasta_entries::Vector{FastaEntry}; max_shuffle_attempts::Int64 = 20)::Vector{FastaEntry}

Creates decoy sequences for target peptides by reversing all but the last amino acid. If reversal
creates a duplicate sequence, falls back to shuffling.

Parameters:
- target_fasta_entries::Vector{FastaEntry}: Vector of target peptide entries to generate decoys for

Keyword Arguments:
- max_shuffle_attempts::Int64 = 20: Maximum number of attempts to generate unique shuffled sequence
  when reversal creates a duplicate

Returns:
- Vector{FastaEntry}: Sorted vector containing both original entries and their decoys. Decoy entries
  have is_decoy=true and maintain their original base_pep_id and entrapment_group_id.

Notes:
- Preserves C-terminal amino acid to maintain enzymatic cleavage properties
- Orders entries by sequence for efficient lookup
- Maintains original metadata (base_pep_id, entrapment_group_id) for tracking
- Warns if unable to generate unique decoy after max attempts
- Uses efficient pointer-based sequence manipulation

Example:
```julia
# Original entries
entries = [FastaEntry("P1", "", "human", "PEPTIDER", 1, 0, false),
          FastaEntry("P2", "", "human", "SAMPLEK", 2, 0, false)]

# Add decoys
with_decoys = addReverseDecoys(entries)
# Results in 4 entries: original 2 plus:
# "EPITPEDER" (reverse of "PEPTIDER" keeping R)
# "LPMASEK" (reverse of "SAMPLEK" keeping K)
```
"""
function add_reverse_decoys(target_fasta_entries::Vector{FastaEntry}; max_shuffle_attempts::Int64 = 20)
    # Pre-allocate space for decoy entries
    decoy_fasta_entries = Vector{FastaEntry}(undef, length(target_fasta_entries))
    
    # Set to track unique sequences
    sequences_set = PeptideSequenceSet(target_fasta_entries)
    
    # Initialize position tracking vector (max peptide length of 255 should be sufficient)
    positions = Vector{UInt8}(undef, 255)
    
    n = 1
    for target_entry in target_fasta_entries
        target_sequence = get_sequence(target_entry)
        seq_length = UInt8(length(target_sequence))
        
        # Reset positions vector to initial positions
        for i in 1:seq_length
            positions[i] = i
        end
        
        # For reversal, modify positions vector accordingly (except last position)
        for i in 1:(seq_length-1)
            positions[i] = seq_length - i
        end
        
        # Create reversed sequence (keeping last amino acid)
        decoy_sequence = reverse(target_sequence[1:(end-1)]) * target_sequence[end]
        
        n_shuffle_attempts = 0
        
        # If reversal creates a duplicate, try shuffling
        if decoy_sequence ∈ sequences_set
            # Reset positions vector before shuffling
            for i in 1:seq_length
                positions[i] = i
            end
            
            while n_shuffle_attempts < max_shuffle_attempts
                # Use enhanced shuffle function that updates positions vector
                decoy_sequence = shuffle_fast_with_positions(target_sequence, positions)
                
                if decoy_sequence ∉ sequences_set
                    break
                end
                n_shuffle_attempts += 1
                
                # Reset positions vector before next attempt
                for i in 1:seq_length
                    positions[i] = i
                end
            end
        end
        
        if n_shuffle_attempts >= max_shuffle_attempts
            @warn "Exceeded max shuffle attempts for $(get_sequence(target_entry))"
        else
            # Adjust modification positions based on sequence manipulation
            adjusted_structural_mods = adjust_mod_positions(
                get_structural_mods(target_entry),
                positions,
                seq_length
            )
            
            adjusted_isotopic_mods = adjust_mod_positions(
                get_isotopic_mods(target_entry),
                positions,
                seq_length
            )
            
            # Create decoy entry with adjusted modifications
            decoy_fasta_entries[n] = FastaEntry(
                get_id(target_entry),
                get_description(target_entry),
                get_proteome(target_entry),
                decoy_sequence,
                adjusted_structural_mods,
                adjusted_isotopic_mods,
                get_charge(target_entry),
                get_base_pep_id(target_entry),
                get_base_prec_id(target_entry),
                get_entrapment_group_id(target_entry),
                true  # This is a decoy sequence
            )
            
            n += 1
            push!(sequences_set, decoy_sequence)
        end
    end
    # Sort the peptides by sequence
    return sort(vcat(target_fasta_entries, decoy_fasta_entries[1:n-1]), by = x -> get_sequence(x))
end

"""
    combineSharedPeptides(peptides::Vector{FastaEntry})::Vector{FastaEntry}

Combines entries that share identical peptide sequences by concatenating their protein accessions.
Used to handle peptides that map to multiple proteins while maintaining unique sequences
in the library.

Parameters:
- peptides::Vector{FastaEntry}: Vector of peptide entries that may contain duplicates

Returns:
- Vector{FastaEntry}: Vector of unique peptide entries where shared peptides have concatenated
  accession numbers separated by semicolons. Other metadata from the first encountered instance 
  of each peptide is preserved.

Notes:
- Uses Dictionary for efficient sequence lookup
- Concatenates accession numbers with semicolon separators
- Preserves first encountered description, proteome, base_pep_id, etc.
- Memory efficient by pre-allocating final vector

Example:
```julia
# Original entries with shared sequence
entries = [
    FastaEntry("P1", "desc1", "human", "PEPTIDE", 1, 0, false),
    FastaEntry("P2", "desc2", "human", "PEPTIDE", 2, 0, false),
    FastaEntry("P3", "desc3", "human", "UNIQUE", 3, 0, false)
]

# Combine shared peptides
combined = combineSharedPeptides(entries)
# Results in 2 entries:
# 1. FastaEntry("P1;P2", "desc1", "human", "PEPTIDE", 1, 0, false)
# 2. FastaEntry("P3", "desc3", "human", "UNIQUE", 3, 0, false)
```
"""
function combine_shared_peptides(peptides::Vector{FastaEntry})
    seq_to_fasta_entry = Dictionary{String, FastaEntry}()
    n = 0
    a = 0
    for peptide in peptides
        sequence = get_sequence(peptide)
        sequence_il_equiv = replace(sequence, 'I' => 'L')
        if haskey(seq_to_fasta_entry, sequence_il_equiv)
            a += 1
            fasta_entry = seq_to_fasta_entry[sequence_il_equiv]
            accession = get_id(peptide)*";"*get_id(fasta_entry)
            proteome = get_proteome(peptide)*";"*get_proteome(fasta_entry)
            description = get_description(peptide)*";"*get_description(fasta_entry)
            seq_to_fasta_entry[sequence_il_equiv] = FastaEntry(accession, 
                                                        description, 
                                                        proteome,
                                                        get_sequence(fasta_entry),
                                                        get_structural_mods(fasta_entry),
                                                        get_isotopic_mods(fasta_entry),
                                                        get_charge(fasta_entry),
                                                        get_base_pep_id(fasta_entry),
                                                        get_base_prec_id(fasta_entry),
                                                        get_entrapment_group_id(fasta_entry), 
                                                        is_decoy(fasta_entry)
                                                        )
        else
            n += 1
            insert!(seq_to_fasta_entry, sequence_il_equiv, peptide)
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