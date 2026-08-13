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

# src/fasta/fasta_digest.jl

# Set of valid amino acid characters for fast lookup (O(1) per character)
const VALID_AAS = Set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

const VALID_DIGEST_SPECIFICITIES = ("full", "semi", "semi-n", "semi-c")

"""Normalize and validate a FASTA digestion-specificity setting."""
function normalize_digest_specificity(specificity::AbstractString)::String
    normalized = replace(lowercase(strip(String(specificity))), '_' => '-')
    normalized in VALID_DIGEST_SPECIFICITIES || throw(ArgumentError(
        "specificity must be one of $(join(VALID_DIGEST_SPECIFICITIES, ", ")); got: $specificity"
    ))
    return normalized
end

function _digest_fully_specific_sequence(sequence::AbstractString,
                                         regex::Regex,
                                         max_length::Int,
                                         min_length::Int,
                                         missed_cleavages::Int
                                        )::Tuple{Vector{String}, Vector{UInt32}, Vector{UInt8}}
    
    function add_peptide!(peptides::Vector{String},
                          starts::Vector{UInt32},
                          n::Int,
                          sequence::AbstractString,
                          site::Int,
                          previous_sites::Vector{Int},
                          min_length::Int,
                          max_length::Int,
                          missed_cleavages::Int)
        
        for i in 1:min(n, missed_cleavages + 1)
            previous_site = previous_sites[end - i + 1]
            if ((site - previous_site) >= min_length) && 
                ((site - previous_site) <= max_length)
                push!(peptides, String(@view sequence[previous_site+1:site]))
                push!(starts, UInt32(previous_site + 1))
            end
        end

        for i in 1:length(previous_sites)-1
            previous_sites[i] = previous_sites[i+1]
        end
        previous_sites[end] = site
        return n + 1
    end

    peptides = String[]
    starts = Vector{UInt32}()
    previous_sites = zeros(Int, missed_cleavages + 1)
    previous_sites[1] = 0
    n = 1

    for site in eachmatch(regex, sequence, overlap = true)
        n = add_peptide!(peptides, starts, n, sequence, site.offset,
                        previous_sites, min_length, max_length,
                        missed_cleavages)
    end

    # Handle C-terminal peptides
    add_peptide!(peptides, starts, n, sequence, length(sequence),
                 previous_sites, min_length, max_length,
                 missed_cleavages)

    return peptides, starts, fill(UInt8(2), length(peptides))
end

"""
    digest_sequence(sequence, regex, max_length, min_length, missed_cleavages,
                    specificity) -> (peptides, starts, num_enzymatic_termini)

Digest a protein with configurable enzymatic specificity. `specificity` may be:

- `"full"`: both peptide termini must be enzymatic (the historical behavior),
- `"semi"`: at least one peptide terminus must be enzymatic,
- `"semi-n"`: the N terminus may be non-enzymatic; the C terminus must be enzymatic,
- `"semi-c"`: the C terminus may be non-enzymatic; the N terminus must be enzymatic.

Protein termini count as enzymatic.
"""
function digest_sequence(sequence::AbstractString,
                         regex::Regex,
                         max_length::Int,
                         min_length::Int,
                         missed_cleavages::Int,
                         specificity::AbstractString
                        )::Tuple{Vector{String}, Vector{UInt32}, Vector{UInt8}}
    normalized = normalize_digest_specificity(specificity)

    if normalized == "full"
        return _digest_fully_specific_sequence(
            sequence, regex, max_length, min_length, missed_cleavages
        )
    end

    isempty(sequence) && return String[], UInt32[], UInt8[]

    sequence_length = length(sequence)
    cleavage_mask = falses(sequence_length)
    for site in eachmatch(regex, sequence, overlap = true)
        1 <= site.offset <= sequence_length || continue
        cleavage_mask[site.offset] = true
    end

    # Prefix counts make the missed-cleavage test O(1) for every emitted peptide.
    cleavage_prefix = zeros(Int, sequence_length)
    running = 0
    @inbounds for idx in 1:sequence_length
        running += cleavage_mask[idx]
        cleavage_prefix[idx] = running
    end

    specific_ends = findall(cleavage_mask)
    if isempty(specific_ends) || specific_ends[end] != sequence_length
        push!(specific_ends, sequence_length)
    end

    @inline start_is_enzymatic(start_idx::Int) =
        start_idx == 1 || cleavage_mask[start_idx - 1]
    @inline end_is_enzymatic(end_idx::Int) =
        end_idx == sequence_length || cleavage_mask[end_idx]
    @inline function internal_cleavages(start_idx::Int, end_idx::Int)
        end_idx <= start_idx && return 0
        before_end = cleavage_prefix[end_idx - 1]
        before_start = start_idx > 1 ? cleavage_prefix[start_idx - 1] : 0
        return before_end - before_start
    end

    peptides = String[]
    starts = UInt32[]
    enzymatic_termini = UInt8[]

    function emit_candidate!(start_idx::Int, end_idx::Int, start_enzymatic::Bool)
        internal_cleavages(start_idx, end_idx) <= missed_cleavages || return
        end_enzymatic = end_is_enzymatic(end_idx)
        push!(peptides, String(@view sequence[start_idx:end_idx]))
        push!(starts, UInt32(start_idx))
        push!(enzymatic_termini,
              UInt8(start_enzymatic ? 1 : 0) + UInt8(end_enzymatic ? 1 : 0))
        return
    end

    for start_idx in 1:sequence_length
        min_end = start_idx + min_length - 1
        min_end > sequence_length && continue
        max_end = min(start_idx + max_length - 1, sequence_length)
        start_enzymatic = start_is_enzymatic(start_idx)

        # semi-c allows an arbitrary C terminus but requires an enzymatic N
        # terminus. General semi digestion does the same for enzymatic starts.
        if normalized == "semi-c" || (normalized == "semi" && start_enzymatic)
            start_enzymatic || continue
            for end_idx in min_end:max_end
                emit_candidate!(start_idx, end_idx, start_enzymatic)
            end
            continue
        end

        # semi-n, and general semi peptides with a non-enzymatic start, require
        # an enzymatic C terminus. Binary search avoids scanning every cleavage
        # site for every possible start.
        first_end = searchsortedfirst(specific_ends, min_end)
        last_end = searchsortedlast(specific_ends, max_end)
        first_end > last_end && continue
        for end_pos in first_end:last_end
            emit_candidate!(start_idx, specific_ends[end_pos], start_enzymatic)
        end
    end

    return peptides, starts, enzymatic_termini
end

"""
    digest_fasta(fasta::Vector{FastaEntry}, proteome_id::String; regex::Regex = r"[KR][^P|\$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1, specificity::AbstractString = "full")::Vector{FastaEntry}

Enzymatically digest protein sequences from FASTA entries into peptides.

# Parameters
- `fasta::Vector{FastaEntry}`: FASTA entries to digest
- `proteome_id::String`: Proteome identifier to assign to resulting peptides
- `regex::Regex`: Enzyme cleavage pattern (default: trypsin-like, cleaves after K or R except when followed by P)
- `max_length::Int`: Maximum peptide length to include (default: 40)
- `min_length::Int`: Minimum peptide length to include (default: 8)
- `missed_cleavages::Int`: Maximum missed cleavages allowed (default: 1)
- `specificity::AbstractString`: `"full"`, `"semi"`, `"semi-n"`, or `"semi-c"` (default: `"full"`)

# Returns
- `Vector{FastaEntry}`: Digested peptide entries as FastaEntry objects

# Details
For each protein in the input FASTA entries:
1. Digests the sequence using the specified enzyme pattern
2. Filters out peptides containing non-standard amino acids (only A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y are allowed)
3. Creates a new FastaEntry for each valid peptide
4. Assigns sequential base_pep_id values

The resulting FastaEntry objects contain:
- Original entry ID
- Empty description (to save memory)
- Specified proteome ID
- Digested peptide sequence
- No structural or isotopic modifications (missing)
- Zero charge state
- Number of enzymatic termini (one or two for the supported modes)
- Sequential base_pep_id
- Zero entrapment_group_id
- is_decoy set to false

# Examples
```julia
# Load FASTA entries
fasta_entries = parse_fasta("proteins.fasta", "human")

# Digest with default parameters (trypsin-like)
peptides = digest_fasta(fasta_entries, "human")

# Digest with custom parameters (AspN-like, cutting before D)
peptides = digest_fasta(
    fasta_entries, 
    "human",
    regex = r"[^D](?=[D])", 
    max_length = 30,
    min_length = 6,
    missed_cleavages = 2
)
```
"""
function digest_fasta(fasta::Vector{FastaEntry},
                     proteome_id::String;
                     regex::Regex = r"[KR][^P|$]",
                     max_length::Int = 40,
                     min_length::Int = 8,
                     missed_cleavages::Int = 1,
                     specificity::AbstractString = "full")::Vector{FastaEntry}

    peptides_fasta = Vector{FastaEntry}()
    base_pep_id = one(UInt32)
    for entry in fasta
        peptides, starts, enzymatic_termini = digest_sequence(
            get_sequence(entry),
            regex,
            max_length,
            min_length,
            missed_cleavages,
            specificity,
        )
        for (peptide, start_idx, num_enzymatic_termini) in
            zip(peptides, starts, enzymatic_termini)

            # Skip peptides containing non-standard amino acids
            if !all(aa -> aa ∈ VALID_AAS, peptide)
                continue
            end

            push!(peptides_fasta, FastaEntry(
                get_id(entry),
                "", # Skip description to save memory
                get_gene(entry),
                get_protein(entry),
                get_organism(entry),
                proteome_id,
                peptide,  # Now String instead of SubString
                start_idx,
                missing, #structural_mods 
                missing, #isotopic_mods 
                zero(UInt8),
                num_enzymatic_termini,
                zero(UInt32),  # base_target_id (will be assigned later)
                base_pep_id,
                zero(UInt32),  # entrapment_pair_id
                false          # is_decoy
            ))
            base_pep_id += one(UInt32)
        end
    end

    return peptides_fasta
end
