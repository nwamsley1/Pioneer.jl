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
    digest_sequence(sequence::AbstractString, regex::Regex, max_length::Int, min_length::Int, missed_cleavages::Int)::Vector{String}

Digest a protein sequence into peptides using the specified enzyme cleavage pattern.

# Parameters
- `sequence::AbstractString`: The protein sequence to digest
- `regex::Regex`: Enzyme cleavage pattern (matches cleavage sites)
- `max_length::Int`: Maximum peptide length to include
- `min_length::Int`: Minimum peptide length to include
- `missed_cleavages::Int`: Maximum number of internal cleavage sites allowed in a peptide

# Returns
- `Vector{String}`: Vector of digested peptide sequences as Strings

# Details
The function simulates enzymatic digestion by:
1. Finding all cleavage sites matching the regex pattern
2. Generating peptides between sites, including those with missed cleavages
3. Filtering peptides by length constraints
4. Converting all SubStrings to Strings to ensure memory stability

# Examples
```julia
# Trypsin-like digestion (cleaves after K or R except when followed by P)
peptides = digest_sequence("MKVGPKAFRVLTEDEMAKR", r"[KR][^P]", 20, 5, 1)
# Returns ["MKVGPK", "VGPKAFR", "VLTEDEMAK", "VLTEDEMAKR", "AFRVLTEDEMAK"]

# No missed cleavages
peptides = digest_sequence("MKVGPKAFRVLTEDEMAKR", r"[KR][^P]", 20, 5, 0)
# Returns ["VLTEDEMAK"]
```
"""
function digest_sequence(sequence::AbstractString,
                        regex::Regex,
                        max_length::Int,
                        min_length::Int,
                        missed_cleavages::Int)::Tuple{Vector{String}, Vector{UInt32}}
    
    function add_peptide!(peptides::Vector{SubString{String}},
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
                # Convert SubString to String when adding to peptides
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

    peptides = Vector{SubString{String}}()
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
    n = add_peptide!(peptides, starts, n, sequence, length(sequence),
                    previous_sites, min_length, max_length,
                    missed_cleavages)

    return peptides, starts
end

"""
    digest_fasta(fasta::Vector{FastaEntry}, proteome_id::String; regex::Regex = r"[KR][^P|\$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)::Vector{FastaEntry}

Enzymatically digest protein sequences from FASTA entries into peptides.

# Parameters
- `fasta::Vector{FastaEntry}`: FASTA entries to digest
- `proteome_id::String`: Proteome identifier to assign to resulting peptides
- `regex::Regex`: Enzyme cleavage pattern (default: trypsin-like, cleaves after K or R except when followed by P)
- `max_length::Int`: Maximum peptide length to include (default: 40)
- `min_length::Int`: Minimum peptide length to include (default: 8)
- `missed_cleavages::Int`: Maximum missed cleavages allowed (default: 1)

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
                     missed_cleavages::Int = 1)::Vector{FastaEntry}

    peptides_fasta = Vector{FastaEntry}()
    base_pep_id = one(UInt32)
    for entry in fasta
        peptides, starts = digest_sequence(
            get_sequence(entry),
            regex,
            max_length,
            min_length,
            missed_cleavages,
        )
        for (peptide, start_idx) in zip(peptides, starts)

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
                false          # is_decoy
            ))
            base_pep_id += one(UInt32)
        end
    end

    return peptides_fasta
end
