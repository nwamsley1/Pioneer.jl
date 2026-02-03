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
"""
    digest_sequence(sequence::AbstractString,
                    regex::Regex,
                    max_length::Int,
                    min_length::Int,
                    missed_cleavages::Int,
                    specificity::String = "full")::Tuple{Vector{String}, Vector{UInt32}, Vector{UInt8}}

Digest a protein sequence into peptides using the specified enzyme cleavage pattern.

# Parameters
- `sequence::AbstractString`: The protein sequence to digest
- `regex::Regex`: Enzyme cleavage pattern (matches cleavage sites)
- `max_length::Int`: Maximum peptide length to include
- `min_length::Int`: Minimum peptide length to include
- `missed_cleavages::Int`: Maximum number of internal cleavage sites allowed in a peptide
- `specificity::String`: Digestion specificity ("full", "semi-n", "semi-c", or "semi")

# Returns
- `Tuple{Vector{String}, Vector{UInt32}, Vector{UInt8}}`: Peptide sequences, start indices, and number of enzymatic termini

# Details
The function simulates enzymatic digestion by:
1. Finding all cleavage sites matching the regex pattern
2. Generating peptides between sites, including those with missed cleavages
3. Filtering peptides by length constraints
4. Converting all SubStrings to Strings to ensure memory stability

# Examples
```julia
# Trypsin-like digestion (cleaves after K or R except when followed by P)
peptides, _, _ = digest_sequence("MKVGPKAFRVLTEDEMAKR", r"[KR][^P]", 20, 5, 1)
# Returns ["MKVGPK", "VGPKAFR", "VLTEDEMAK", "VLTEDEMAKR", "AFRVLTEDEMAK"]

# No missed cleavages
peptides, _, _ = digest_sequence("MKVGPKAFRVLTEDEMAKR", r"[KR][^P]", 20, 5, 0)
# Returns ["VLTEDEMAK"]
```
"""
function digest_sequence(sequence::AbstractString,
                         regex::Regex,
                         max_length::Int,
                         min_length::Int,
                         missed_cleavages::Int,
                         specificity::String = "full")::Tuple{Vector{String}, Vector{UInt32}, Vector{UInt8}}
    if isempty(sequence)
        return String[], UInt32[], UInt8[]
    end

    normalized_specificity = lowercase(specificity)
    if !(normalized_specificity in ("full", "semi-n", "semi-c", "semi"))
        error("specificity must be one of \"full\", \"semi-n\", \"semi-c\", or \"semi\"; got: $specificity")
    end

    require_n_cleavage = normalized_specificity in ("full", "semi-c")
    require_c_cleavage = normalized_specificity in ("full", "semi-n")
    require_one_cleavage = normalized_specificity == "semi"

    sequence_length = lastindex(sequence)
    cleavage_mask = falses(sequence_length)
    for site in eachmatch(regex, sequence, overlap = true)
        if 1 <= site.offset <= sequence_length
            cleavage_mask[site.offset] = true
        end
    end

    cleavage_prefix = zeros(Int, sequence_length)
    running = 0
    for idx in 1:sequence_length
        if cleavage_mask[idx]
            running += 1
        end
        cleavage_prefix[idx] = running
    end

    function start_is_enzymatic(start_idx::Int)
        return start_idx == 1 || (start_idx > 1 && cleavage_mask[start_idx - 1])
    end

    function end_is_enzymatic(end_idx::Int)
        return end_idx == sequence_length || cleavage_mask[end_idx]
    end

    function internal_cleavages(start_idx::Int, end_idx::Int)
        if end_idx <= 1
            return 0
        end
        before_end = cleavage_prefix[end_idx - 1]
        before_start = start_idx > 1 ? cleavage_prefix[start_idx - 1] : 0
        return before_end - before_start
    end

    peptides = String[]
    starts = UInt32[]
    enzymatic_counts = UInt8[]
    for start_idx in 1:sequence_length
        start_enzymatic = start_is_enzymatic(start_idx)
        if require_n_cleavage && !start_enzymatic
            continue
        end

        min_end = start_idx + min_length - 1
        if min_end > sequence_length
            continue
        end
        max_end = min(start_idx + max_length - 1, sequence_length)
        for end_idx in min_end:max_end
            end_enzymatic = end_is_enzymatic(end_idx)
            if require_c_cleavage && !end_enzymatic
                continue
            end
            if require_one_cleavage && !(start_enzymatic || end_enzymatic)
                continue
            end
            if internal_cleavages(start_idx, end_idx) > missed_cleavages
                continue
            end

            push!(peptides, String(@view sequence[start_idx:end_idx]))
            push!(starts, UInt32(start_idx))
            num_enzymatic = UInt8((start_enzymatic ? 1 : 0) + (end_enzymatic ? 1 : 0))
            push!(enzymatic_counts, num_enzymatic)
        end
    end

    return peptides, starts, enzymatic_counts
end

"""
    digest_fasta(fasta::Vector{FastaEntry}, proteome_id::String; regex::Regex = r"[KR][^P|\$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1, specificity::String = "full")::Vector{FastaEntry}

Enzymatically digest protein sequences from FASTA entries into peptides.

# Parameters
- `fasta::Vector{FastaEntry}`: FASTA entries to digest
- `proteome_id::String`: Proteome identifier to assign to resulting peptides
- `regex::Regex`: Enzyme cleavage pattern (default: trypsin-like, cleaves after K or R except when followed by P)
- `max_length::Int`: Maximum peptide length to include (default: 40)
- `min_length::Int`: Minimum peptide length to include (default: 8)
- `missed_cleavages::Int`: Maximum missed cleavages allowed (default: 1)
- `specificity::String`: Digestion specificity ("full", "semi-n", "semi-c", or "semi")

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
                     specificity::String = "full")::Vector{FastaEntry}

    peptides_fasta = Vector{FastaEntry}()
    base_pep_id = one(UInt32)
    for entry in fasta
        peptides, starts, enzymatic_counts = digest_sequence(
            get_sequence(entry),
            regex,
            max_length,
            min_length,
            missed_cleavages,
            specificity,
        )
        for (peptide, start_idx, num_enzymatic_termini) in zip(peptides, starts, enzymatic_counts)

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
                zero(UInt32),  # base_target_id (will be assigned later)
                base_pep_id,
                zero(UInt8),   # entrapment_pair_id
                num_enzymatic_termini,
                false          # is_decoy
            ))
            base_pep_id += one(UInt32)
        end
    end

    return peptides_fasta
end
