# src/fasta/fasta_utils.jl
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
    
    # Get original sequences for uniqueness checking
    target_sequences = map(get_sequence, target_fasta_entries)
    
    # Pre-allocate output vector
    entrapment_fasta_entries = Vector{FastaEntry}(
        undef, 
        length(target_fasta_entries) * entrapment_r
    )
    
    # Track unique sequences
    sequences_set = Set{String}()
    sizehint!(sequences_set, length(entrapment_fasta_entries) + length(target_fasta_entries))
    union!(sequences_set, target_sequences)
    
    n = 1
    for target_entry in target_fasta_entries
        for entrapment_group_id in 1:entrapment_r
            n_shuffle_attempts = 0
            
            while n_shuffle_attempts < max_shuffle_attempts
                new_sequence = shuffle_fast(get_sequence(target_entry))
                
                if new_sequence ∉ sequences_set
                    entrapment_fasta_entries[n] = FastaEntry(
                        get_id(target_entry),
                        get_description(target_entry),
                        get_proteome(target_entry),
                        new_sequence,
                        get_base_pep_id(target_entry),
                        entrapment_group_id,
                        false
                    )
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
    println("add_reverse_decoys")
    #Get the sequences for the target entries 
    target_sequences = (map(x->get_sequence(x), target_fasta_entries))
    #Pre-allocate space for entrapment fasta entries 
    decoy_fasta_entries = Vector{FastaEntry}(undef, length(target_fasta_entries))
    #Set to keep track of encountered sequences. Do not want to generate non-unique entrapment sequences
    sequences_set = Set{String}()
    sizehint!(sequences_set, length(decoy_fasta_entries) + length(target_fasta_entries))
    union!(sequences_set, target_sequences)
    n = 1
    #For eacbh protein in the FASTA
    for target_entry in target_fasta_entries
        #Make a target and decoy for each peptide
        target_sequence = get_sequence(target_entry)
        decoy_sequence = reverse(target_sequence[1:(end -1)])*target_sequence[end]
        n_shuffle_attempts = 0
        #If reversal fails to generate a unique sequence, then shuffle 
        if decoy_sequence ∈ sequences_set
            while n_shuffle_attempts < max_shuffle_attempts
                decoy_sequence = shuffle_fast(get_sequence(target_entry))
                if decoy_sequence ∉ sequences_set
                    break
                end
                n_shuffle_attempts += 1
            end
        end
        if n_shuffle_attempts >= max_shuffle_attempts
            @warn "Exceeded max shuffle attempts for $target_entry"
        else
            decoy_fasta_entries[n] = FastaEntry(
                                            get_id(target_entry), 
                                            get_description(target_entry),# This justs wastes time and memory here 
                                            get_proteome(target_entry),
                                            decoy_sequence,
                                            get_base_pep_id(target_entry),
                                            get_entrapment_group_id(target_entry),
                                            true #Must be a decoy sequence
                                            )
            n += 1
            push!(sequences_set, decoy_sequence)
        end
    end
    #Sort the peptides
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
        if haskey(seq_to_fasta_entry, sequence)
            a += 1
            fasta_entry = seq_to_fasta_entry[sequence]
            accession = get_id(peptide)*";"*get_id(fasta_entry)
            seq_to_fasta_entry[sequence] = FastaEntry(accession, 
                                                        get_description(fasta_entry), 
                                                        get_proteome(fasta_entry),
                                                        get_sequence(fasta_entry),
                                                        get_base_pep_id(fasta_entry),
                                                        get_entrapment_group_id(fasta_entry), 
                                                        is_decoy(fasta_entry)
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