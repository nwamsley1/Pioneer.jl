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
    PeptideSequenceSet

A data structure for efficiently storing and comparing peptide sequences with I/L equivalence.

The structure replaces all isoleucine (I) with leucine (L) in stored sequences to treat these
amino acids as equivalent for the purpose of comparisons.

# Fields
- `sequences::Set{String}`: Set of stored sequences with I replaced by L

# Methods
- `PeptideSequenceSet()`: Constructor to create an empty set
- `PeptideSequenceSet(fasta_entries::Vector{FastaEntry})`: Constructor to initialize from FASTA entries
- `push!(pss::PeptideSequenceSet, seq::AbstractString)`: Add a sequence to the set
- `in(seq::AbstractString, pss::PeptideSequenceSet)`: Check if a sequence exists in the set

# Examples
```julia
# Create an empty set
pss = PeptideSequenceSet()

# Add sequences
push!(pss, "PEPTIDE")
push!(pss, "PEPTLDE")  # Contains L instead of I

# Check membership (both return true due to I/L equivalence)
"PEPTIDE" in pss  # true
"PEPTLDE" in pss  # true

# Initialize from FASTA entries
pss = PeptideSequenceSet(fasta_entries)
```
"""
struct PeptideSequenceSet
    sequences::Set{Tuple{String, UInt8}}
    # Constructor
    function PeptideSequenceSet()
        new(Set{Tuple{String, UInt8}}())
    end
    function PeptideSequenceSet(fasta_entries::Vector{FastaEntry})
        pss = PeptideSequenceSet()
        # Add each sequence to the set
        for entry in fasta_entries
            push!(pss, get_sequence(entry), get_charge(entry))
        end
        return pss
    end
end

getSeqSet(s::PeptideSequenceSet) = s.sequences

import Base: push!
function push!(pss::PeptideSequenceSet, seq::AbstractString, charge::UInt8)
    # Replace 'I' with 'L' in the sequence and add to the set
    push!(pss.sequences, (replace(seq, 'I' => 'L'), charge))
    return pss
end

import Base: in
function in(seq_charge::Tuple{String, UInt8}, pss::PeptideSequenceSet)
    # Check if the modified sequence is in the set
    return (replace(first(seq_charge), 'I' => 'L'), last(seq_charge)) ∈ getSeqSet(pss)
end

const PROTEOME_INDEX_MIN_BLOCK_LENGTH = 3
const PROTEOME_INDEX_MAX_BLOCK_LENGTH = 20
const SHORT_EXACT_MAX_LENGTH = 7
const PACKED_AA_BITS = 5
const PACKED_AA_MASK = UInt64(0x1f)
const INVALID_AA_CODE = typemax(UInt8)
const RAW_SEQUENCE_AA_RESIDUES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
const RAW_SEQUENCE_AA_INDEX = Dict(aa => idx for (idx, aa) in enumerate(RAW_SEQUENCE_AA_RESIDUES))
const IL_COLLAPSED_AA_ORDER = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L',
                               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
const MAX_IL_COLLAPSED_AA_CODE = UInt8(length(IL_COLLAPSED_AA_ORDER) - 1)
const IL_COLLAPSED_AA_CODES = Dict{Char, UInt8}(
    'A' => UInt8(0),
    'C' => UInt8(1),
    'D' => UInt8(2),
    'E' => UInt8(3),
    'F' => UInt8(4),
    'G' => UInt8(5),
    'H' => UInt8(6),
    'K' => UInt8(7),
    'L' => UInt8(8),
    'M' => UInt8(9),
    'N' => UInt8(10),
    'P' => UInt8(11),
    'Q' => UInt8(12),
    'R' => UInt8(13),
    'S' => UInt8(14),
    'T' => UInt8(15),
    'V' => UInt8(16),
    'W' => UInt8(17),
    'Y' => UInt8(18),
    'I' => UInt8(8),
)
const TERMINAL_FRAGMENT_RESIDUE_MASS_KEYS = Dict{Char, Int32}(
    'A' => 7103711,
    'R' => 15610111,
    'N' => 11404293,
    'D' => 11502694,
    'C' => 10300919,
    'E' => 12904259,
    'Q' => 12805858,
    'G' => 5702146,
    'H' => 13705891,
    'I' => 11308406,
    'L' => 11308406,
    'K' => 12809496,
    'M' => 13104049,
    'F' => 14706841,
    'P' => 9705276,
    'S' => 8703203,
    'T' => 10104768,
    'W' => 18607931,
    'Y' => 16306333,
    'V' => 9906841,
)
const MUTATION_FALLBACK_RESTARTS = 100
const MUTATION_EXCLUSION_MAP = Dict{Char, Set{Char}}(
    'A' => Set(['R', 'K', 'H']), #Set(['V', 'Q', 'P']),
    'C' => Set(['R', 'K', 'H']), #Set(['A', 'S', 'G']),
    'D' => Set(['R', 'K', 'H']), #Set(['G', 'A', 'H', 'E']),
    'E' => Set(['R', 'K', 'H']), #Set(['D', 'A', 'W']),
    'F' => Set(['R', 'K', 'H']), #Set(['I', 'L', 'V', 'S']),
    'G' => Set(['R', 'K', 'H']), #Set(['N', 'A', 'S']),
    'H' => Set([]), #Set(['D', 'N', 'Q']),
    'I' => Set(['R', 'K', 'H', 'L']), #Set(['V', 'T', 'M', 'L']),
    'K' => Set([]), #Set(['G', 'R', 'Q']),
    'L' => Set(['R', 'K', 'H', 'I']), #Set(['V', 'P', 'F', 'I']),
    'M' => Set(['R', 'K', 'H']), #Set(['T', 'I', 'L', 'F']),
    'N' => Set(['R', 'K', 'H']), #Set(['S', 'G', 'E']),
    'P' => Set(['R', 'K', 'H']), #Set(['I', 'L', 'S', 'A']),
    'Q' => Set(['R', 'K', 'H']), #Set(['H', 'D', 'W']),
    'R' => Set([]), #Set(['H', 'K', 'Q']),
    'S' => Set(['R', 'K', 'H']), #Set(['V', 'N', 'A']),
    'T' => Set(['R', 'K', 'H']), #Set(['A', 'S', 'I', 'L']),
    'V' => Set(['R', 'K', 'H']), #Set(['I', 'L', 'A', 'M']),
    'W' => Set(['R', 'K', 'H']), #Set(['D', 'S', 'G']),
    'Y' => Set(['R', 'K', 'H']), #Set(['F', 'H', 'D']),
)
const KRH_RESTRICTED_MUTATION_CHOICES = Dict{Char, Vector{Char}}(
    'K' => ['R', 'H'],
    'R' => ['K', 'H'],
    'H' => ['R', 'K'],
)
const KRH_MUTATION_DESTINATION_BLOCKLIST = Set(keys(KRH_RESTRICTED_MUTATION_CHOICES))
const PROTEOME_INDEX_PROGRESS_CHUNK = 100_000

function collect_variable_mod_aas(
    variable_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
)
    variable_mod_aas = Set{Char}()
    isempty(variable_mod_names) && return variable_mod_aas

    for aa in RAW_SEQUENCE_AA_RESIDUES
        aa_string = string(aa)
        for mod in variable_mod_names
            if occursin(mod[:p], aa_string)
                push!(variable_mod_aas, aa)
                break
            end
        end
    end

    return variable_mod_aas
end

struct ProteomeDistanceIndex
    sequence_codes::Vector{UInt8}
    invalid_prefix::Vector{Int}
    block_indexes::Vector{Dict{UInt128, Vector{Int}}}
    short_exact_indexes::Vector{Set{UInt64}}
    amino_acid_counts::Vector{Int}
    min_block_length::Int
    max_block_length::Int
    short_exact_max_length::Int
end

function encode_sequence_codes(sequence::AbstractString)
    codes = Vector{UInt8}(undef, length(sequence))
    for (i, aa) in enumerate(sequence)
        codes[i] = get(IL_COLLAPSED_AA_CODES, aa, INVALID_AA_CODE)
    end
    return codes
end

function encode_block_key(codes::AbstractVector{UInt8}, start::Int, block_length::Int)
    key = zero(UInt128)
    @inbounds for idx in start:(start + block_length - 1)
        key = (key << 5) | UInt128(codes[idx])
    end
    return key
end

function encode_short_key(codes::AbstractVector{UInt8}, start::Int, block_length::Int)
    key = zero(UInt64)
    @inbounds for idx in start:(start + block_length - 1)
        key = (key << PACKED_AA_BITS) | UInt64(codes[idx])
    end
    return key
end

function ProteomeDistanceIndex(
    fasta_entries::Vector{FastaEntry};
    min_block_length::Int = PROTEOME_INDEX_MIN_BLOCK_LENGTH,
    max_block_length::Int = PROTEOME_INDEX_MAX_BLOCK_LENGTH,
    show_progress::Bool = false,
    log_progress::Bool = false,
)
    sequence_codes = UInt8[]
    amino_acid_counts = zeros(Int, length(RAW_SEQUENCE_AA_RESIDUES))
    total_length = sum(length(get_sequence(entry)) + 1 for entry in fasta_entries)
    total_residues = total_length - length(fasta_entries)
    total_block_lengths = max(0, max_block_length - min_block_length + 1)
    short_exact_max_length = SHORT_EXACT_MAX_LENGTH
    if log_progress
        @user_info "Building proteome distance index from $(length(fasta_entries)) proteins ($total_residues residues) with block lengths $min_block_length:$max_block_length"
        @user_info "Encoding proteome sequence buffer..."
    end
    encode_pbar = (show_progress && !isempty(fasta_entries)) ?
        ProgressBar(
            total=length(fasta_entries),
            leave=false,
            printing_delay=0.25,
        ) : nothing
    sizehint!(sequence_codes, total_length)

    for entry in fasta_entries
        sequence = get_sequence(entry)
        append!(sequence_codes, encode_sequence_codes(sequence))
        for aa in sequence
            aa_index = get(RAW_SEQUENCE_AA_INDEX, aa, 0)
            if aa_index > 0
                amino_acid_counts[aa_index] += 1
            end
        end
        push!(sequence_codes, INVALID_AA_CODE)
        if !isnothing(encode_pbar)
            update(encode_pbar)
        end
    end

    invalid_prefix = Vector{Int}(undef, length(sequence_codes) + 1)
    invalid_prefix[1] = 0
    for idx in eachindex(sequence_codes)
        invalid_prefix[idx + 1] = invalid_prefix[idx] + (sequence_codes[idx] == INVALID_AA_CODE)
    end

    block_indexes = [Dict{UInt128, Vector{Int}}() for _ in 1:max_block_length]
    short_exact_indexes = [Set{UInt64}() for _ in 1:short_exact_max_length]
    max_start = length(sequence_codes)
    total_windows = 0
    total_block_work = sum(max(max_start - block_length + 1, 0) for block_length in min_block_length:max_block_length)
    if log_progress && total_block_work > 0
        @user_info "Building proteome block indexes..."
    end
    block_pbar = (show_progress && total_block_work > 0) ?
        ProgressBar(
            total=total_block_work,
            leave=false,
            printing_delay=0.25,
        ) : nothing
    block_progress_pending = 0
    for block_length in min_block_length:max_block_length
        if max_start < block_length
            continue
        end
        index = block_indexes[block_length]
        for start in 1:(max_start - block_length + 1)
            if invalid_prefix[start + block_length] != invalid_prefix[start]
                if !isnothing(block_pbar)
                    block_progress_pending += 1
                    if block_progress_pending >= PROTEOME_INDEX_PROGRESS_CHUNK
                        update(block_pbar, block_progress_pending)
                        block_progress_pending = 0
                    end
                end
                continue
            end
            key = encode_block_key(sequence_codes, start, block_length)
            push!(get!(index, key, Int[]), start)
            total_windows += 1
            if !isnothing(block_pbar)
                block_progress_pending += 1
                if block_progress_pending >= PROTEOME_INDEX_PROGRESS_CHUNK
                    update(block_pbar, block_progress_pending)
                    block_progress_pending = 0
                end
            end
        end
    end
    if !isnothing(block_pbar) && block_progress_pending > 0
        update(block_pbar, block_progress_pending)
    end

    total_short_exact_work = sum(max(max_start - peptide_length + 1, 0) for peptide_length in 1:short_exact_max_length)
    if log_progress && total_short_exact_work > 0
        @user_info "Building short exact target indexes..."
    end
    short_exact_pbar = (show_progress && total_short_exact_work > 0) ?
        ProgressBar(
            total=total_short_exact_work,
            leave=false,
            printing_delay=0.25,
        ) : nothing
    short_exact_progress_pending = 0
    for peptide_length in 1:short_exact_max_length
        if max_start >= peptide_length
            exact_index = short_exact_indexes[peptide_length]
            for start in 1:(max_start - peptide_length + 1)
                if invalid_prefix[start + peptide_length] != invalid_prefix[start]
                    if !isnothing(short_exact_pbar)
                        short_exact_progress_pending += 1
                        if short_exact_progress_pending >= PROTEOME_INDEX_PROGRESS_CHUNK
                            update(short_exact_pbar, short_exact_progress_pending)
                            short_exact_progress_pending = 0
                        end
                    end
                    continue
                end
                push!(exact_index, encode_short_key(sequence_codes, start, peptide_length))
                if !isnothing(short_exact_pbar)
                    short_exact_progress_pending += 1
                    if short_exact_progress_pending >= PROTEOME_INDEX_PROGRESS_CHUNK
                        update(short_exact_pbar, short_exact_progress_pending)
                        short_exact_progress_pending = 0
                    end
                end
            end
        end
    end
    if !isnothing(short_exact_pbar) && short_exact_progress_pending > 0
        update(short_exact_pbar, short_exact_progress_pending)
    end

    if log_progress
        @user_info "Proteome distance index ready: $total_windows valid windows across $total_block_lengths block lengths"
    end

    return ProteomeDistanceIndex(
        sequence_codes,
        invalid_prefix,
        block_indexes,
        short_exact_indexes,
        amino_acid_counts,
        min_block_length,
        max_block_length,
        short_exact_max_length,
    )
end

function has_invalid_window(index::ProteomeDistanceIndex, start::Int, window_length::Int)
    stop = start + window_length - 1
    if start < 1 || stop > length(index.sequence_codes)
        return true
    end
    return index.invalid_prefix[stop + 1] != index.invalid_prefix[start]
end

function collect_candidate_starts!(
    starts::Vector{Int},
    seen::Set{Int},
    hits::Union{Nothing, Vector{Int}},
    offset::Int,
    candidate_length::Int,
    index::ProteomeDistanceIndex,
)
    isnothing(hits) && return nothing
    for hit in hits
        start = hit - offset
        if start < 1 || (start + candidate_length - 1) > length(index.sequence_codes)
            continue
        end
        if start ∉ seen
            push!(seen, start)
            push!(starts, start)
        end
    end
    return nothing
end

function has_exact_target_match(
    candidate_codes::AbstractVector{UInt8},
    index::ProteomeDistanceIndex,
)
    peptide_length = length(candidate_codes)
    if peptide_length <= index.short_exact_max_length
        exact_index = index.short_exact_indexes[peptide_length]
        isempty(exact_index) && return false
        return encode_short_key(candidate_codes, 1, peptide_length) ∈ exact_index
    end

    left_length = fld(peptide_length, 2)
    right_length = peptide_length - left_length

    if left_length < index.min_block_length || right_length < index.min_block_length
        error("ProteomeDistanceIndex was built for block lengths $(index.min_block_length):$(index.max_block_length), cannot query peptide length $peptide_length")
    end

    left_key = encode_block_key(candidate_codes, 1, left_length)
    right_key = encode_block_key(candidate_codes, left_length + 1, right_length)
    left_hits = get(index.block_indexes[left_length], left_key, nothing)
    right_hits = get(index.block_indexes[right_length], right_key, nothing)

    if isnothing(left_hits) || isnothing(right_hits)
        return false
    end

    use_left_hits = length(left_hits) <= length(right_hits)
    hits = use_left_hits ? left_hits : right_hits
    offset = use_left_hits ? 0 : left_length

    for hit in hits
        start = hit - offset
        if has_invalid_window(index, start, peptide_length)
            continue
        end

        exact_match = true
        @inbounds for peptide_idx in 1:peptide_length
            if index.sequence_codes[start + peptide_idx - 1] != candidate_codes[peptide_idx]
                exact_match = false
                break
            end
        end

        if exact_match
            return true
        end
    end

    return false
end

function has_target_adjacent_swap_match(
    candidate_codes::AbstractVector{UInt8},
    index::ProteomeDistanceIndex,
)
    peptide_length = length(candidate_codes)
    peptide_length <= 1 && return false

    swap_codes = copy(candidate_codes)
    for position in 1:(peptide_length - 1)
        if swap_codes[position] == swap_codes[position + 1]
            continue
        end

        swap_codes[position], swap_codes[position + 1] = swap_codes[position + 1], swap_codes[position]
        has_match = has_exact_target_match(swap_codes, index)
        swap_codes[position], swap_codes[position + 1] = swap_codes[position + 1], swap_codes[position]

        if has_match
            return true
        end
    end

    return false
end

function has_min_target_hamming_distance(
    candidate_codes::AbstractVector{UInt8},
    index::ProteomeDistanceIndex,
    min_target_hamming_distance::Int = 2,
)
    if min_target_hamming_distance < 1 || min_target_hamming_distance > 2
        error("min_target_hamming_distance must be 1 or 2, got $min_target_hamming_distance")
    end

    peptide_length = length(candidate_codes)
    if peptide_length <= index.short_exact_max_length
        exact_index = index.short_exact_indexes[peptide_length]
        isempty(exact_index) && return true

        base_key = encode_short_key(candidate_codes, 1, peptide_length)
        if base_key ∈ exact_index
            return false
        end

        if min_target_hamming_distance == 1
            return true
        end

        for position in 1:peptide_length
            bit_shift = PACKED_AA_BITS * (peptide_length - position)
            cleared_key = base_key & ~(PACKED_AA_MASK << bit_shift)
            original_code = candidate_codes[position]

            for alt_code in UInt8(0):MAX_IL_COLLAPSED_AA_CODE
                alt_code == original_code && continue
                if (cleared_key | (UInt64(alt_code) << bit_shift)) ∈ exact_index
                    return false
                end
            end
        end

        return true
    end

    left_length = fld(peptide_length, 2)
    right_length = peptide_length - left_length

    if left_length < index.min_block_length || right_length < index.min_block_length
        error("ProteomeDistanceIndex was built for block lengths $(index.min_block_length):$(index.max_block_length), cannot query peptide length $peptide_length")
    end

    left_key = encode_block_key(candidate_codes, 1, left_length)
    right_key = encode_block_key(candidate_codes, left_length + 1, right_length)
    left_hits = get(index.block_indexes[left_length], left_key, nothing)
    right_hits = get(index.block_indexes[right_length], right_key, nothing)

    if isnothing(left_hits) && isnothing(right_hits)
        return true
    end

    candidate_starts = Int[]
    seen = Set{Int}()
    collect_candidate_starts!(candidate_starts, seen, left_hits, 0, peptide_length, index)
    collect_candidate_starts!(candidate_starts, seen, right_hits, left_length, peptide_length, index)

    for start in candidate_starts
        if has_invalid_window(index, start, peptide_length)
            continue
        end

        mismatches = 0
        @inbounds for offset in 1:peptide_length
            if index.sequence_codes[start + offset - 1] != candidate_codes[offset]
                mismatches += 1
                mismatches >= min_target_hamming_distance && break
            end
        end

        if mismatches < min_target_hamming_distance
            return false
        end
    end

    return true
end

function has_min_target_hamming_distance(
    candidate::String,
    index::ProteomeDistanceIndex,
    min_target_hamming_distance::Int = 2,
)
    candidate_codes = encode_sequence_codes(candidate)
    if any(==(INVALID_AA_CODE), candidate_codes)
        return false
    end
    return has_min_target_hamming_distance(candidate_codes, index, min_target_hamming_distance)
end

function passes_target_distance_constraints(
    candidate::String,
    index::ProteomeDistanceIndex,
    min_target_hamming_distance::Int = 2,
)
    candidate_codes = encode_sequence_codes(candidate)
    if any(==(INVALID_AA_CODE), candidate_codes)
        return false
    end

    if !has_min_target_hamming_distance(candidate_codes, index, min_target_hamming_distance)
        return false
    end

    return !has_target_adjacent_swap_match(candidate_codes, index)
end

function terminal_fragment_length(peptide_length::Int)::Int
    return peptide_length <= 7 ? 3 : 4
end

function terminal_fragment_mass_key(
    sequence::String,
    start_idx::Int,
    stop_idx::Int,
)::Union{Nothing, Int64}
    mass_key = zero(Int64)
    @inbounds for idx in start_idx:stop_idx
        aa_mass_key = get(TERMINAL_FRAGMENT_RESIDUE_MASS_KEYS, sequence[idx], Int32(0))
        iszero(aa_mass_key) && return nothing
        mass_key += aa_mass_key
    end
    return mass_key
end

function has_distinct_terminal_fragment_masses(
    base_seq::String,
    candidate::String,
)::Bool
    peptide_length = length(base_seq)
    peptide_length == length(candidate) || return false

    fragment_length = terminal_fragment_length(peptide_length)
    peptide_length >= fragment_length || return false

    base_prefix = terminal_fragment_mass_key(base_seq, 1, fragment_length)
    candidate_prefix = terminal_fragment_mass_key(candidate, 1, fragment_length)
    if isnothing(base_prefix) || isnothing(candidate_prefix) || base_prefix == candidate_prefix
        return false
    end

    suffix_start = peptide_length - fragment_length + 1
    base_suffix = terminal_fragment_mass_key(base_seq, suffix_start, peptide_length)
    candidate_suffix = terminal_fragment_mass_key(candidate, suffix_start, peptide_length)
    if isnothing(base_suffix) || isnothing(candidate_suffix) || base_suffix == candidate_suffix
        return false
    end

    return true
end

function is_valid_decoy_candidate(
    candidate::String,
    base_seq::String,
    charges::Vector{UInt8},
    sequences_set::PeptideSequenceSet,
    proteome_index::Union{Nothing, ProteomeDistanceIndex},
    min_target_hamming_distance::Int = 2,
)
    if !has_distinct_terminal_fragment_masses(base_seq, candidate)
        return false
    end
    if any(((candidate, charge) ∈ sequences_set) for charge in charges)
        return false
    end
    return isnothing(proteome_index) || passes_target_distance_constraints(candidate, proteome_index, min_target_hamming_distance)
end

"""
    add_entrapment_sequences(
        target_fasta_entries::Vector{FastaEntry}, 
        entrapment_r::UInt8;
        max_shuffle_attempts::Int64 = 20
    )::Vector{FastaEntry}

Add entrapment sequences to a set of target peptides with modification handling.
Creates shuffled sequences while properly adjusting modification positions.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of target peptide entries (with modifications) to generate entrapment sequences for
- `entrapment_r::UInt8`: Number of entrapment sequences to generate per target
- `max_shuffle_attempts::Int64`: Maximum attempts to generate unique shuffled sequence (default: 20)

# Returns
- `Vector{FastaEntry}`: Combined vector of original entries and their entrapment sequences

# Details
For each target peptide:
1. Creates `entrapment_r` shuffled versions using `shuffle_sequence!()`
2. Adjusts modification positions to match shuffled sequence using `adjust_mod_positions`
3. Preserves C-terminal amino acid to maintain enzymatic cleavage properties  
4. Ensures each shuffled sequence is unique (I/L equivalence considered)
5. Sets entrapment_group_id to indicate the entrapment group
6. Maintains original metadata (base_target_id, base_pep_id, etc.) for tracking
7. Properly handles both structural and isotopic modifications

# Examples
```julia
# Create 2 entrapment sequences per target
entries_with_entrapment = add_entrapment_sequences(target_entries, UInt8(2))

# Create 3 entrapment sequences with more shuffle attempts for difficult sequences
entries_with_entrapment = add_entrapment_sequences(
    target_entries, 
    UInt8(3), 
    max_shuffle_attempts=50
)
```

# Notes
Entrapment sequences help assess false discovery rates in peptide identification.
The function uses I/L equivalence when checking for sequence uniqueness.
"""
function add_entrapment_sequences(
    target_fasta_entries::Vector{FastaEntry}, 
    entrapment_r::UInt8;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    entrapment_method::String = "shuffle"
)::Vector{FastaEntry}
    
    # Pre-allocate output vector
    entrapment_fasta_entries = Vector{FastaEntry}(
        undef, 
        length(target_fasta_entries) * entrapment_r
    )
    
    # Counters for tracking fallback to shuffle
    total_sequences = length(target_fasta_entries) * entrapment_r
    fallback_to_shuffle_count = 0
    
    # Track unique sequences
    sequences_set = PeptideSequenceSet(target_fasta_entries)#Set{String}()
    #sizehint!(sequences_set, length(entrapment_fasta_entries) + length(target_fasta_entries))
    #union!(sequences_set, target_sequences)
    
    n = 1

    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
       fixed_chars#['R','K']#Vector{Char}()
    )
    for target_entry in target_fasta_entries
        for entrapment_group_id in 1:entrapment_r
            n_shuffle_attempts = 0
            # Try the specified method first
            new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method=entrapment_method)
            
            # If it creates a duplicate and we're using reverse, fall back to shuffle
            if (new_sequence, get_charge(target_entry)) ∈ sequences_set && entrapment_method == "reverse"
                @debug_l2 "Reverse created duplicate for entrapment of $(get_sequence(target_entry)), falling back to shuffle"
                fallback_to_shuffle_count += 1
                # Fall back to shuffle since reverse is deterministic
                while n_shuffle_attempts < max_shuffle_attempts
                    new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method="shuffle")
                    
                    if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            elseif (new_sequence, get_charge(target_entry)) ∈ sequences_set
                # For shuffle method, keep trying with shuffle
                while n_shuffle_attempts < max_shuffle_attempts
                    new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method="shuffle")
                    
                    if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            end
            
            #Make sure the entrapment sequence is unique (I and L are equivalent)
            if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                    # Get sequence length for modification adjustment
                    seq_length = UInt8(length(get_sequence(target_entry)))
                    
                    # Adjust modification positions based on sequence shuffling
                    adjusted_structural_mods = adjust_mod_positions(
                        get_structural_mods(target_entry),
                        shuffle_seq.new_positions,
                        seq_length
                    )
                    
                    adjusted_isotopic_mods = adjust_mod_positions(
                        get_isotopic_mods(target_entry),
                        shuffle_seq.new_positions,
                        seq_length
                    )
                    
                    entrapment_fasta_entries[n] = FastaEntry(
                        get_id(target_entry),
                        get_description(target_entry),
                        get_gene(target_entry),
                        get_protein(target_entry),
                        get_organism(target_entry),
                        get_proteome(target_entry),
                        new_sequence,
                        get_start_idx(target_entry),
                        adjusted_structural_mods, #structural_mods - now properly adjusted
                        adjusted_isotopic_mods,   #isotopic_mods - now properly adjusted
                        get_charge(target_entry),
                        get_base_target_id(target_entry), # inherit base_target_id for tracking
                        get_base_pep_id(target_entry),
                        entrapment_group_id,
                        false
                    )
                    n += 1
                    push!(sequences_set, new_sequence, get_charge(target_entry))
            elseif n_shuffle_attempts >= max_shuffle_attempts
                @user_warn "Max shuffle attempts exceeded for $(get_sequence(target_entry))"
            end
        end
    end
    
    # Report statistics if using reverse method for entrapment
    #=
    if entrapment_method == "reverse" && total_sequences > 0
        if fallback_to_shuffle_count > 0
            @user_warn "Entrapment generation statistics for REVERSE method:"
            @user_warn "  Total entrapment sequences attempted: $total_sequences"
            @user_warn "  Sequences where reverse created duplicates: $fallback_to_shuffle_count"
            @user_warn "  Sequences successfully reversed: $(total_sequences - fallback_to_shuffle_count)"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_sequences, digits=1))%"
        end
    end
    =#
    
    return vcat(target_fasta_entries, entrapment_fasta_entries[1:n-1])
end

"""
    add_entrapment_sequences_grouped(
        target_fasta_entries::Vector{FastaEntry},
        entrapment_r::UInt8;
        max_shuffle_attempts::Int64 = 20,
        fixed_chars::Vector{Char} = Vector{Char}(),
        entrapment_method::String = "shuffle"
    )::Vector{FastaEntry}

Group-aware entrapment generation that ensures all modification variants of the
same base peptide sequence receive the same set of entrapment sequences.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of peptide entries (after modifications)
- `entrapment_r::UInt8`: Number of entrapment sequences to generate per base sequence
- `max_shuffle_attempts::Int64`: Max attempts to find a unique shuffled sequence
- `fixed_chars::Vector{Char}`: Optional characters to keep fixed when shuffling
- `fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Fixed modification patterns to reapply on final entrapment sequences
- `variable_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Variable modification patterns whose source and destination residues are excluded from mutation fallback
- `entrapment_method::String`: "shuffle" or "reverse" (reverse may fall back to shuffle)
- `proteome_index::Union{Nothing, ProteomeDistanceIndex}`: Optional target proteome index used to enforce decoy-style distance constraints
- `cleavage_regex::Regex`: Digest cleavage-site regex; internal matches are protected from mutation fallback
- `show_progress::Bool`: Show a progress bar while iterating grouped base sequences
- `log_progress::Bool`: Emit start/end summary logs for grouped generation

# Returns
- `Vector{FastaEntry}`: Combined vector of original entries and grouped entrapment entries

# Details
Algorithm:
1. Group entries by base sequence (ignoring modifications)
2. For each base sequence, generate `entrapment_r` unique entrapment sequences once
3. Reuse those sequences for all modification variants in the group, adjusting mod positions
4. Apply the same target-distance constraints and mutation fallback strategy used for decoys
5. Maintain I/L equivalence when checking uniqueness (via PeptideSequenceSet)
"""
function add_entrapment_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry},
    entrapment_r::UInt8;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(),
    variable_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(),
    entrapment_method::String = "shuffle",
    proteome_index::Union{Nothing, ProteomeDistanceIndex} = nothing,
    cleavage_regex::Regex = r"[KR][^P|$]",
    show_progress::Bool = true,
    log_progress::Bool = true,
)::Vector{FastaEntry}

    if iszero(entrapment_r)
        return target_fasta_entries
    end

    # Track existing sequences (I/L equivalence) including charges
    sequences_set = PeptideSequenceSet(target_fasta_entries)

    # Prepare shuffler
    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
        fixed_chars
    )
    variable_mod_aas = collect_variable_mod_aas(variable_mod_names)

    # Group entries by base sequence (ignore mods)
    groups = Dict{String, Vector{Int}}()
    for (idx, entry) in enumerate(target_fasta_entries)
        base_seq = get_sequence(entry)
        if !haskey(groups, base_seq)
            groups[base_seq] = Vector{Int}()
        end
        push!(groups[base_seq], idx)
    end

    # High-level diagnostics
    n_entries = length(target_fasta_entries)
    n_groups = length(groups)
    avg_variants = n_groups == 0 ? 0.0 : round(n_entries / n_groups, digits=2)

    entrapments_out = Vector{FastaEntry}()
    fallback_to_shuffle_count = 0
    mutation_fallback_count = 0
    total_attempted = length(groups) * Int(entrapment_r)
    pbar = (show_progress && n_groups > 0) ? ProgressBar(total=n_groups) : nothing
    if log_progress
        @user_info "Entrapment generation: $n_entries entries across $n_groups base sequence groups (avg variants/group=$(avg_variants), ratio=$(Int(entrapment_r)), method=$entrapment_method, min_target_hamming_distance=$(isnothing(proteome_index) ? 0 : 2))"
    end

    exhausted_groups = 0
    sample_logged = 0
    for (base_seq, idxs) in groups
        # Collect unique charges observed among variants (usually 0 at this stage)
        charges = unique([get_charge(target_fasta_entries[i]) for i in idxs])

        # Generate unique entrapment sequences for this base sequence
        entrap_seqs = Vector{String}()
        entrap_positions = Vector{Vector{UInt8}}()
        for entrapment_group_id in 1:entrapment_r
            new_sequence, positions_copy, used_shuffle_fallback, used_mutation_fallback, _ = find_group_decoy_candidate!(
                shuffle_seq,
                base_seq,
                idxs,
                target_fasta_entries,
                charges,
                sequences_set;
                max_shuffle_attempts = max_shuffle_attempts,
                fixed_chars = fixed_chars,
                variable_mod_aas = variable_mod_aas,
                decoy_method = entrapment_method,
                proteome_index = proteome_index,
                cleavage_regex = cleavage_regex,
            )

            if isnothing(new_sequence)
                exhausted_groups += 1
                continue
            end

            if used_shuffle_fallback
                fallback_to_shuffle_count += 1
            end
            if used_mutation_fallback
                mutation_fallback_count += 1
            end

            push!(entrap_seqs, new_sequence)
            push!(entrap_positions, positions_copy)

            # Reserve the sequence globally for all observed charges
            for c in charges
                push!(sequences_set, new_sequence, c)
            end
        end

        # Optional per-group diagnostics (debug level)
        if sample_logged < 5
            trunc_seq = length(base_seq) > 18 ? (base_seq[1:18] * "…") : base_seq
            sample_logged += 1
        end

        # Create entrapment entries for all variants using the same entrapment specs
        for idx in idxs
            target_entry = target_fasta_entries[idx]
            seq_length = UInt8(length(base_seq))
            for i in 1:length(entrap_seqs)
                adjusted_structural_mods = adjust_mod_positions(
                    get_structural_mods(target_entry),
                    entrap_positions[i],
                    seq_length
                )
                adjusted_structural_mods = merge_fixed_mods(
                    entrap_seqs[i],
                    adjusted_structural_mods,
                    fixed_mod_names,
                )
                adjusted_isotopic_mods = adjust_mod_positions(
                    get_isotopic_mods(target_entry),
                    entrap_positions[i],
                    seq_length
                )

                push!(entrapments_out, FastaEntry(
                    get_id(target_entry),
                    get_description(target_entry),
                    get_gene(target_entry),
                    get_protein(target_entry),
                    get_organism(target_entry),
                    get_proteome(target_entry),
                    entrap_seqs[i],
                    get_start_idx(target_entry),
                    adjusted_structural_mods,
                    adjusted_isotopic_mods,
                    get_charge(target_entry),
                    get_base_target_id(target_entry),
                    get_base_pep_id(target_entry),
                    UInt8(i),
                    false
                ))
            end
        end

        if !isnothing(pbar)
            update(pbar)
        end
    end

    if log_progress
        @user_info "Entrapment generation complete: $(length(entrapments_out)) entrapment entries created, exhausted_groups=$exhausted_groups, shuffle_fallbacks=$fallback_to_shuffle_count, mutation_fallbacks=$mutation_fallback_count"
    end

    # Report statistics if using reverse method for entrapment
    if entrapment_method == "reverse" && total_attempted > 0
        if fallback_to_shuffle_count > 0
            @user_warn "Entrapment generation (GROUPED) stats for REVERSE:"
            @user_warn "  Total entrapments attempted: $total_attempted"
            @user_warn "  Reverse duplicates: $fallback_to_shuffle_count"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_attempted, digits=1))%"
        end
    end

    # General summary

    return vcat(target_fasta_entries, entrapments_out)
end

mutable struct ShuffleSeq
    old_sequence::String 
    new_sequence::Vector{Char}
    new_positions::Vector{UInt8}
    movable_positions::Vector{UInt8}
    n_movable::UInt8
    sequence_length::UInt8
    fixed_chars::Vector{Char}
end

"""
    resetSequence!(shuffle_sequence::ShuffleSeq, sequence::String)

Resets the 'old_sequence' and 'new_sequence' attributes of the `ShuffleSeq` object.
and resets the 'sequence_length' attribute of the 'shuffle_sequence'
"""
function resetSequence!(shuffle_sequence::ShuffleSeq, sequence::String)

    #Fills the character string for the 'new_sequence' to match 'sequence'
    #and resets the 'sequence_length' attribute of the 'shuffle_sequence'
    shuffle_sequence.old_sequence = sequence 
    shuffle_sequence.sequence_length = length(sequence)
    for i in range(one(UInt8), UInt8(shuffle_sequence.sequence_length))
        shuffle_sequence.new_sequence[i] = sequence[i]
        shuffle_sequence.new_positions[i] = i
    end
    shuffle_sequence.n_movable = zero(UInt8) 
    return nothing
end

"""
    fillMovablePositions!(ss::ShuffleSeq, sequence) 
    
Fills the character string for the 'new_sequence' to match 'sequence'
and resets the 'sequence_length' attribute of the 'shuffle_sequence'
"""
function fillMovablePositions!(shuffle_sequence::ShuffleSeq)

    shuffle_sequence.n_movable = zero(UInt8)
    #shuffle_sequence.sequence_length-1 because the last amino-acid is fixed 
    for i in range(one(UInt8), shuffle_sequence.sequence_length-1)
        if shuffle_sequence.old_sequence[i] ∉ shuffle_sequence.fixed_chars
            shuffle_sequence.n_movable += one(UInt8)
            # If the character is fixed, keep it in the same position
            shuffle_sequence.movable_positions[shuffle_sequence.n_movable] = i
        end
    end
    return nothing
end

function permuteNewPositions!(shuffle_sequence::ShuffleSeq)
    perm = randperm(shuffle_sequence.n_movable)
    # Update new_positions based on the permutation
    for (new_idx, old_idx) in enumerate(perm)
        # Update sequence 
        shuffle_sequence.new_sequence[shuffle_sequence.movable_positions[new_idx]] = 
            shuffle_sequence.old_sequence[shuffle_sequence.movable_positions[old_idx]]
        # Update positions 
        shuffle_sequence.new_positions[shuffle_sequence.movable_positions[new_idx]] = 
            shuffle_sequence.movable_positions[old_idx]
    end
    return nothing
end

function reverseMovablePositions!(shuffle_sequence::ShuffleSeq)
    # Reverse the movable positions (all but the last amino acid)
    n = shuffle_sequence.n_movable
    
    # Create reversed mapping
    for i in 1:n
        old_pos = shuffle_sequence.movable_positions[i]
        new_pos = shuffle_sequence.movable_positions[n - i + 1]
        
        # Update sequence - place character from position i at position (n-i+1)
        shuffle_sequence.new_sequence[old_pos] = 
            shuffle_sequence.old_sequence[new_pos]
        
        # Update position mapping
        shuffle_sequence.new_positions[old_pos] = new_pos
    end
    return nothing
end

function shuffle_sequence!(
    shuffle_sequence::ShuffleSeq,
    sequence::String;
    method::String = "shuffle"
)
    # Reset the sequence and positions
    resetSequence!(shuffle_sequence, sequence)
    
    # Fill movable positions
    fillMovablePositions!(shuffle_sequence)
    
    # Apply the selected decoy generation method
    if method == "shuffle"
        permuteNewPositions!(shuffle_sequence)
    elseif method == "reverse"
        reverseMovablePositions!(shuffle_sequence)
    else
        error("Unknown decoy method: $method. Must be 'shuffle' or 'reverse'")
    end
    
    return String(shuffle_sequence.new_sequence[1:shuffle_sequence.sequence_length])
end

"""
Enhanced version of shuffle_fast_with_positions that keeps specified characters fixed.
Pre-allocated vectors are passed in to avoid allocations.
"""
function shuffle_fast_with_positions_and_fixed_chars!(
    s::String, 
    positions::Vector{UInt8}, 
    fixed_chars::Set{Char},
    fixed_positions::Vector{Int},  # Pre-allocated
    movable_positions::Vector{Int},  # Pre-allocated
    temp_positions::Vector{UInt8}  # Pre-allocated for temporary storage
)
    ss = sizeof(s)
    l = length(s)
    
    # Count fixed and movable positions
    n_fixed = 0
    n_movable = 0
    
    # Create indices vector for byte positions
    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end
    
    # Identify fixed and movable positions
    p = pointer(s)
    for j in 1:l
        c = Char(unsafe_load(p, v[j]))
        if j == l || c in fixed_chars  # Last position or fixed character
            n_fixed += 1
            fixed_positions[n_fixed] = j
        else
            n_movable += 1
            movable_positions[n_movable] = j
        end
    end
    
    # Generate permutation only for movable positions
    if n_movable > 0
        perm = randperm(n_movable)
        
        # Copy positions for movable characters to temp storage
        for idx in 1:n_movable
            temp_positions[idx] = positions[movable_positions[idx]]
        end
        
        # Apply permutation to positions
        for (new_idx, old_idx) in enumerate(perm)
            positions[movable_positions[new_idx]] = temp_positions[old_idx]
        end
        
        # Build the output string
        u = Vector{UInt8}(undef, ss)
        
        # Fill in the shuffled string
        for j in 1:l
            # Check if j is in fixed_positions (up to n_fixed)
            is_fixed = false
            for k in 1:n_fixed
                if fixed_positions[k] == j
                    is_fixed = true
                    break
                end
            end
            
            if is_fixed
                # Keep fixed characters in place
                u[v[j]] = unsafe_load(p, v[j])
            else
                # Find position in movable_positions
                idx = 0
                for k in 1:n_movable
                    if movable_positions[k] == j
                        idx = k
                        break
                    end
                end
                source_pos = movable_positions[perm[idx]]
                u[v[j]] = unsafe_load(p, v[source_pos])
            end
        end
        
        return String(u)
    else
        # All positions are fixed, return original string
        return s
    end
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

function build_fixed_mods_for_sequence(
    sequence::String,
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
)::Vector{PeptideMod}
    fixed_mods = PeptideMod[]
    for mod in fixed_mod_names
        getFixedMods!(fixed_mods, eachmatch(mod[:p], sequence), mod[:r])
    end
    sort!(fixed_mods)
    return fixed_mods
end

function merge_fixed_mods(
    sequence::String,
    adjusted_mods::Union{Missing, Vector{PeptideMod}},
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
)::Union{Missing, Vector{PeptideMod}}
    isempty(fixed_mod_names) && return adjusted_mods

    fixed_mod_names_set = Set(mod[:r] for mod in fixed_mod_names)
    merged_mods = PeptideMod[]

    if !ismissing(adjusted_mods) && !isempty(adjusted_mods)
        for mod in adjusted_mods
            if mod.mod_name ∉ fixed_mod_names_set
                push!(merged_mods, mod)
            end
        end
    end

    append!(merged_mods, build_fixed_mods_for_sequence(sequence, fixed_mod_names))

    if isempty(merged_mods)
        return ismissing(adjusted_mods) ? missing : PeptideMod[]
    end

    sort!(merged_mods)
    return merged_mods
end

function identity_positions(seq_length::Int)
    positions = Vector{UInt8}(undef, seq_length)
    for idx in 1:seq_length
        positions[idx] = UInt8(idx)
    end
    return positions
end

function collect_protected_positions(
    target_fasta_entries::Vector{FastaEntry},
    idxs::Vector{Int},
    base_seq::String,
    fixed_chars::Vector{Char},
    cleavage_regex::Regex = r"[KR][^P|$]",
)
    seq_length = length(base_seq)
    protected = falses(seq_length)
    protected[end] = true

    mark_internal_cleavage_positions!(protected, base_seq, cleavage_regex)

    for idx in 1:(seq_length - 1)
        if base_seq[idx] ∈ fixed_chars
            protected[idx] = true
        end
    end

    for entry_idx in idxs
        entry = target_fasta_entries[entry_idx]
        for mods in (get_structural_mods(entry), get_isotopic_mods(entry))
            if ismissing(mods) || isempty(mods)
                continue
            end
            for mod in mods
                if mod.aa == 'n'
                    protected[1] = true
                elseif mod.aa == 'c'
                    protected[end] = true
                elseif 1 <= mod.position <= seq_length
                    protected[Int(mod.position)] = true
                end
            end
        end
    end

    return protected
end

function mark_internal_cleavage_positions!(
    positions::AbstractVector{Bool},
    sequence::String,
    cleavage_regex::Regex,
)
    seq_length = length(sequence)
    for match in eachmatch(cleavage_regex, sequence, overlap = true)
        cleavage_position = match.offset
        if 1 <= cleavage_position < seq_length
            positions[cleavage_position] = true
        end
    end
    return positions
end

function internal_cleavage_positions(
    sequence::String,
    cleavage_regex::Regex,
)::BitVector
    positions = falses(length(sequence))
    mark_internal_cleavage_positions!(positions, sequence, cleavage_regex)
    return positions
end

function has_new_internal_cleavage_site(
    sequence::String,
    allowed_cleavage_positions::AbstractVector{Bool},
    cleavage_regex::Regex,
)::Bool
    seq_length = length(sequence)
    for match in eachmatch(cleavage_regex, sequence, overlap = true)
        cleavage_position = match.offset
        if 1 <= cleavage_position < seq_length && !allowed_cleavage_positions[cleavage_position]
            return true
        end
    end
    return false
end

function mutation_creates_new_internal_cleavage_site(
    candidate_chars::Vector{Char},
    position::Int,
    mutation_aa::Char,
    allowed_cleavage_positions::AbstractVector{Bool},
    cleavage_regex::Regex,
)::Bool
    original_aa = candidate_chars[position]
    candidate_chars[position] = mutation_aa
    creates_new = has_new_internal_cleavage_site(
        String(candidate_chars),
        allowed_cleavage_positions,
        cleavage_regex,
    )
    candidate_chars[position] = original_aa
    return creates_new
end

function inward_mutation_position_groups(seq_length::Int)
    position_groups = Vector{Vector{Int}}()
    left = 2
    right = seq_length - 1
    if seq_length <= 2 || left > right
        return position_groups
    end

    if isodd(seq_length)
        target_left = fld(seq_length + 1, 2)
        target_right = target_left
    else
        target_left = fld(seq_length, 2)
        target_right = target_left + 1
    end

    push!(position_groups, left == right ? [left] : [left, right])
    while left != target_left || right != target_right
        if left != target_left
            left += 1
            push!(position_groups, [left])
        end
        if right != target_right
            right -= 1
            if right != left
                push!(position_groups, [right])
            end
        end
    end

    return position_groups
end

function get_mutation_choices(
    source_aa::Char,
    variable_mod_aas::Set{Char},
)::Vector{Char}
    if source_aa ∈ variable_mod_aas
        return Char[]
    end

    source_code = get(IL_COLLAPSED_AA_CODES, source_aa, INVALID_AA_CODE)
    if haskey(KRH_RESTRICTED_MUTATION_CHOICES, source_aa)
        return [
            candidate_aa for candidate_aa in KRH_RESTRICTED_MUTATION_CHOICES[source_aa]
            if candidate_aa ∉ variable_mod_aas &&
               candidate_aa != source_aa &&
               get(IL_COLLAPSED_AA_CODES, candidate_aa, INVALID_AA_CODE) != source_code
        ]
    end

    excluded_aas = get(MUTATION_EXCLUSION_MAP, source_aa, Set{Char}())
    choices = Char[]
    for candidate_aa in RAW_SEQUENCE_AA_RESIDUES
        candidate_aa == source_aa && continue
        candidate_aa ∈ variable_mod_aas && continue
        candidate_aa ∈ KRH_MUTATION_DESTINATION_BLOCKLIST && continue
        candidate_aa ∈ excluded_aas && continue
        get(IL_COLLAPSED_AA_CODES, candidate_aa, INVALID_AA_CODE) == source_code && continue
        push!(choices, candidate_aa)
    end

    return choices
end

function sample_mutation_choice(
    choices::Vector{Char},
    amino_acid_counts::Union{Nothing,Vector{Int}},
)
    isempty(choices) && return nothing
    if isnothing(amino_acid_counts)
        return rand(choices)
    end

    total_weight = 0
    for aa in choices
        aa_index = get(RAW_SEQUENCE_AA_INDEX, aa, 0)
        total_weight += aa_index > 0 ? amino_acid_counts[aa_index] : 0
    end
    if total_weight <= 0
        return rand(choices)
    end

    threshold = rand(1:total_weight)
    running_weight = 0
    for aa in choices
        aa_index = get(RAW_SEQUENCE_AA_INDEX, aa, 0)
        running_weight += aa_index > 0 ? amino_acid_counts[aa_index] : 0
        if threshold <= running_weight
            return aa
        end
    end

    return choices[end]
end

function try_mutation_fallback(
    base_seq::String,
    idxs::Vector{Int},
    target_fasta_entries::Vector{FastaEntry},
    charges::Vector{UInt8},
    sequences_set::PeptideSequenceSet,
    proteome_index::Union{Nothing, ProteomeDistanceIndex},
    fixed_chars::Vector{Char},
    variable_mod_aas::Set{Char} = Set{Char}(),
    min_target_hamming_distance::Int = 2,
    shuffle_seq::Union{Nothing, ShuffleSeq} = nothing,
    cleavage_regex::Regex = r"[KR][^P|$]",
)
    protected = collect_protected_positions(target_fasta_entries, idxs, base_seq, fixed_chars, cleavage_regex)
    mutation_position_groups = inward_mutation_position_groups(length(base_seq))
    isempty(mutation_position_groups) && return nothing, nothing

    fallback_shuffle_seq = isnothing(shuffle_seq) ? ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
        fixed_chars
    ) : shuffle_seq
    amino_acid_counts = isnothing(proteome_index) ? nothing : proteome_index.amino_acid_counts
    seq_length = length(base_seq)

    for _ in 1:MUTATION_FALLBACK_RESTARTS
        shuffled_candidate = shuffle_sequence!(fallback_shuffle_seq, base_seq; method = "shuffle")
        candidate_chars = collect(shuffled_candidate)
        positions = Vector{UInt8}(fallback_shuffle_seq.new_positions[1:seq_length])
        allowed_cleavage_positions = internal_cleavage_positions(shuffled_candidate, cleavage_regex)
        for mutation_positions in mutation_position_groups
            mutated_in_group = false
            for position in mutation_positions
                protected[Int(positions[position])] && continue
                mutation_choices = get_mutation_choices(candidate_chars[position], variable_mod_aas)
                filter!(
                    mutation_aa -> !mutation_creates_new_internal_cleavage_site(
                        candidate_chars,
                        position,
                        mutation_aa,
                        allowed_cleavage_positions,
                        cleavage_regex,
                    ),
                    mutation_choices,
                )
                mutation_aa = sample_mutation_choice(mutation_choices, amino_acid_counts)
                isnothing(mutation_aa) && continue

                candidate_chars[position] = mutation_aa
                mutated_in_group = true
            end

            if mutated_in_group
                candidate = String(candidate_chars)
                if is_valid_decoy_candidate(candidate, base_seq, charges, sequences_set, proteome_index, min_target_hamming_distance)
                    return candidate, positions
                end
            end
        end
    end

    return nothing, nothing
end

function find_group_decoy_candidate_for_distance!(
    shuffle_seq::ShuffleSeq,
    base_seq::String,
    idxs::Vector{Int},
    target_fasta_entries::Vector{FastaEntry},
    charges::Vector{UInt8},
    sequences_set::PeptideSequenceSet;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    variable_mod_aas::Set{Char} = Set{Char}(),
    decoy_method::String = "shuffle",
    proteome_index::Union{Nothing, ProteomeDistanceIndex} = nothing,
    min_target_hamming_distance::Int = 2,
    cleavage_regex::Regex = r"[KR][^P|$]",
)
    candidate = shuffle_sequence!(shuffle_seq, base_seq; method = decoy_method)
    if is_valid_decoy_candidate(candidate, base_seq, charges, sequences_set, proteome_index, min_target_hamming_distance)
        return candidate, Vector{UInt8}(shuffle_seq.new_positions), false, false
    end

    shuffle_attempts = decoy_method == "shuffle" ? max(max_shuffle_attempts - 1, 0) : max_shuffle_attempts
    used_shuffle_fallback = decoy_method == "reverse"
    for _ in 1:shuffle_attempts
        candidate = shuffle_sequence!(shuffle_seq, base_seq; method = "shuffle")
        if is_valid_decoy_candidate(candidate, base_seq, charges, sequences_set, proteome_index, min_target_hamming_distance)
            return candidate, Vector{UInt8}(shuffle_seq.new_positions), used_shuffle_fallback, false
        end
    end

    candidate, positions = try_mutation_fallback(
        base_seq,
        idxs,
        target_fasta_entries,
        charges,
        sequences_set,
        proteome_index,
        fixed_chars,
        variable_mod_aas,
        min_target_hamming_distance,
        shuffle_seq,
        cleavage_regex,
    )
    if !isnothing(candidate)
        return candidate, positions, used_shuffle_fallback, true
    end

    return nothing, nothing, used_shuffle_fallback, false
end

function find_group_decoy_candidate!(
    shuffle_seq::ShuffleSeq,
    base_seq::String,
    idxs::Vector{Int},
    target_fasta_entries::Vector{FastaEntry},
    charges::Vector{UInt8},
    sequences_set::PeptideSequenceSet;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    variable_mod_aas::Set{Char} = Set{Char}(),
    decoy_method::String = "shuffle",
    proteome_index::Union{Nothing, ProteomeDistanceIndex} = nothing,
    cleavage_regex::Regex = r"[KR][^P|$]",
)
    candidate, positions, used_shuffle_fallback, used_mutation_fallback = find_group_decoy_candidate_for_distance!(
        shuffle_seq,
        base_seq,
        idxs,
        target_fasta_entries,
        charges,
        sequences_set;
        max_shuffle_attempts = max_shuffle_attempts,
        fixed_chars = fixed_chars,
        variable_mod_aas = variable_mod_aas,
        decoy_method = decoy_method,
        proteome_index = proteome_index,
        min_target_hamming_distance = 2,
        cleavage_regex = cleavage_regex,
    )
    return candidate, positions, used_shuffle_fallback, used_mutation_fallback, false
end


"""
    add_decoy_sequences(target_fasta_entries::Vector{FastaEntry}; max_shuffle_attempts::Int64 = 20)

Creates decoy sequences for target peptides by reversing all but the last amino acid.
If reversal creates a duplicate sequence, falls back to shuffling.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of target peptide entries to generate decoys for
- `max_shuffle_attempts::Int64`: Maximum attempts to generate unique shuffled sequence when reversal creates a duplicate (default: 20)

# Returns
- `Vector{FastaEntry}`: Sorted vector containing both original entries and their decoys

# Details
For each target peptide:
1. Reverses the sequence keeping the last amino acid fixed
2. If the resulting sequence already exists, tries shuffling instead
3. Updates modification positions to match the reversed/shuffled sequence
4. Sets is_decoy=true for decoy entries
5. Maintains original metadata (base_pep_id, entrapment_group_id) for tracking
6. Returns a combined list of target and decoy sequences, sorted by sequence

# Examples
```julia
# Add reverse decoys to a set of target entries
all_entries = add_decoy_sequences(target_entries)

# Add decoys with more shuffle attempts
all_entries = add_decoy_sequences(target_entries, max_shuffle_attempts=50)
```

# Notes
- Preserves C-terminal amino acid to maintain enzymatic cleavage properties
- Correctly handles modifications, updating their positions to match the reversed sequence
- Uses I/L equivalence when checking for sequence uniqueness
- Entries are sorted by sequence in the output for efficient lookup
"""
function add_decoy_sequences(
    target_fasta_entries::Vector{FastaEntry}; 
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    decoy_method::String = "shuffle"
    )
    # Pre-allocate space for decoy entries
    decoy_fasta_entries = Vector{FastaEntry}(undef, length(target_fasta_entries))
    
    # Set to track unique sequences
    sequences_set = PeptideSequenceSet(target_fasta_entries)
    
    # Counters for tracking fallback to shuffle
    total_sequences = length(target_fasta_entries)
    fallback_to_shuffle_count = 0
    
    # Initialize position tracking vector (max peptide length of 255 should be sufficient)
    #positions = Vector{UInt8}(undef, 255)
    
    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
        fixed_chars#['R','K']#Vector{Char}()
    )
    n = 1
    for target_entry in target_fasta_entries
        target_sequence = get_sequence(target_entry)
        charge = get_charge(target_entry)
        seq_length = UInt8(length(target_sequence))

        # Create decoy sequence using the specified method
        decoy_sequence = shuffle_sequence!(shuffle_seq, target_sequence; method=decoy_method)
                
        n_shuffle_attempts = 0
        
        # If the decoy creates a duplicate, need to handle differently based on method
        if (decoy_sequence, charge) ∈ sequences_set
            if decoy_method == "reverse"
                # If reverse creates a duplicate, fall back to shuffle
                # (reverse is deterministic, so retrying won't help)
                @debug_l2 "Reverse created duplicate for $target_sequence, falling back to shuffle"
                fallback_to_shuffle_count += 1
                while n_shuffle_attempts < max_shuffle_attempts
                    decoy_sequence = shuffle_sequence!(shuffle_seq, target_sequence; method="shuffle")
                    
                    if (decoy_sequence, charge) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            else
                # For shuffle, keep trying with shuffle
                while n_shuffle_attempts < max_shuffle_attempts
                    decoy_sequence = shuffle_sequence!(shuffle_seq, target_sequence; method="shuffle")
                    
                    if (decoy_sequence, charge) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            end
        end
        
        if n_shuffle_attempts >= max_shuffle_attempts
            @user_warn "Exceeded max shuffle attempts for $(get_sequence(target_entry))"
        else
            # Adjust modification positions based on sequence manipulation
            adjusted_structural_mods = adjust_mod_positions(
                get_structural_mods(target_entry),
                shuffle_seq.new_positions,
                seq_length
            )
            
            adjusted_isotopic_mods = adjust_mod_positions(
                get_isotopic_mods(target_entry),
                shuffle_seq.new_positions,
                seq_length
            )
            
            # Create decoy entry with adjusted modifications
            decoy_fasta_entries[n] = FastaEntry(
                get_id(target_entry),
                get_description(target_entry),
                get_gene(target_entry),
                get_protein(target_entry),
                get_organism(target_entry),
                get_proteome(target_entry),
                decoy_sequence,
                get_start_idx(target_entry),
                adjusted_structural_mods,
                adjusted_isotopic_mods,
                get_charge(target_entry),
                get_base_target_id(target_entry), # inherit base_target_id for tracking
                get_base_pep_id(target_entry),  # inherit base_pep_id for pairing
                get_entrapment_pair_id(target_entry),
                true  # This is a decoy sequence
            )
            
            n += 1
            push!(sequences_set, decoy_sequence, get_charge(target_entry))
        end
    end
    
    # Report statistics if using reverse method
    #=
    if decoy_method == "reverse"
        if fallback_to_shuffle_count > 0
            @user_warn "Decoy generation statistics for REVERSE method:"
            @user_warn "  Total sequences attempted: $total_sequences"
            @user_warn "  Sequences where reverse created duplicates: $fallback_to_shuffle_count"
            @user_warn "  Sequences successfully reversed: $(total_sequences - fallback_to_shuffle_count)"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_sequences, digits=1))%"
        else
            @user_info "Successfully reversed all $total_sequences sequences without duplicates"
        end
    end
    =#
    # Sort the peptides by sequence
    return sort(vcat(target_fasta_entries, decoy_fasta_entries[1:n-1]), by = x -> get_sequence(x))
end

"""
    add_decoy_sequences_grouped(
        target_fasta_entries::Vector{FastaEntry};
        max_shuffle_attempts::Int64 = 20,
        fixed_chars::Vector{Char} = Vector{Char}(),
        decoy_method::String = "shuffle"
    )::Vector{FastaEntry}

Group-aware decoy generation that ensures all modification variants of the same
base peptide sequence share a single decoy sequence and mod position mapping.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Peptide entries to generate decoys for (typically includes targets and entrapments)
- `max_shuffle_attempts::Int64`: Max attempts to find a unique shuffled sequence
- `fixed_chars::Vector{Char}`: Optional set of characters kept fixed when shuffling
- `fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Fixed modification patterns to reapply on the final decoy sequence
- `variable_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Variable modification patterns whose source and destination residues are excluded from mutation fallback
- `decoy_method::String`: "shuffle" or "reverse" (reverse may fall back to shuffle)
- `proteome_index::Union{Nothing, ProteomeDistanceIndex}`: Optional target proteome index used to enforce a minimum Hamming distance of 2 from target substrings
- `cleavage_regex::Regex`: Digest cleavage-site regex; internal matches are protected from mutation fallback
- `show_progress::Bool`: Show a progress bar while iterating grouped base sequences
- `log_progress::Bool`: Emit start/end summary logs for grouped generation

# Returns
- `Vector{FastaEntry}`: Sorted vector with both original entries and their decoys

# Details
Algorithm:
1. Group by base sequence (ignoring modifications)
2. For each base sequence, generate one decoy sequence once (respect I/L equivalence and charges)
3. If shuffling cannot satisfy the acceptance rules, fall back to an inward stochastic mutation strategy
4. Apply the same position mapping to all modification variants in the group
5. Preserve metadata and set `is_decoy = true`
"""
function add_decoy_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry};
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    fixed_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(),
    variable_mod_names::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(),
    decoy_method::String = "shuffle",
    proteome_index::Union{Nothing, ProteomeDistanceIndex} = nothing,
    cleavage_regex::Regex = r"[KR][^P|$]",
    show_progress::Bool = true,
    log_progress::Bool = true,
)::Vector{FastaEntry}

    # Track sequences (I/L equivalence) with charge awareness
    sequences_set = PeptideSequenceSet(target_fasta_entries)

    # Prepare shuffler/reverser
    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
        fixed_chars
    )
    variable_mod_aas = collect_variable_mod_aas(variable_mod_names)

    # Group entries by base sequence (sequence only, ignore mods)
    groups = Dict{String, Vector{Int}}()
    for (idx, entry) in enumerate(target_fasta_entries)
        base_seq = get_sequence(entry)
        if !haskey(groups, base_seq)
            groups[base_seq] = Vector{Int}()
        end
        push!(groups[base_seq], idx)
    end

    # High-level diagnostics
    n_entries = length(target_fasta_entries)
    n_groups = length(groups)
    avg_variants = n_groups == 0 ? 0.0 : round(n_entries / n_groups, digits=2)
    decoy_entries = Vector{FastaEntry}()
    fallback_to_shuffle_count = 0
    total_groups = length(groups)
    pbar = (show_progress && total_groups > 0) ? ProgressBar(total=total_groups) : nothing
    if log_progress
        @user_info "Decoy generation: $n_entries entries across $n_groups base sequence groups (avg variants/group=$(avg_variants), method=$decoy_method, min_target_hamming_distance=$(isnothing(proteome_index) ? 0 : 2))"
    end

    sample_logged = 0
    exhausted_groups = 0
    mutation_fallback_count = 0
    for (base_seq, idxs) in groups
        # Unique charges across variants in this group
        charges = unique([get_charge(target_fasta_entries[i]) for i in idxs])

        decoy_sequence, positions_copy, used_shuffle_fallback, used_mutation_fallback, _ = find_group_decoy_candidate!(
            shuffle_seq,
            base_seq,
            idxs,
            target_fasta_entries,
            charges,
            sequences_set;
            max_shuffle_attempts = max_shuffle_attempts,
            fixed_chars = fixed_chars,
            variable_mod_aas = variable_mod_aas,
            decoy_method = decoy_method,
            proteome_index = proteome_index,
            cleavage_regex = cleavage_regex,
        )

        if isnothing(decoy_sequence)
            exhausted_groups += 1
            if log_progress
                if !isnothing(pbar)
                    println()
                    flush(stdout)
                end
                @user_warn "Decoy generation exhausted group: original_sequence=$base_seq, variants=$(length(idxs)), charges=$charges"
                if !isnothing(pbar)
                    println()
                    flush(stdout)
                end
            end
            if !isnothing(pbar)
                update(pbar)
            end
            continue
        end

        if used_shuffle_fallback
            fallback_to_shuffle_count += 1
        end
        if used_mutation_fallback
            mutation_fallback_count += 1
        end

        seq_length = UInt8(length(base_seq))

        # Reserve the decoy sequence across all charges
        for c in charges
            push!(sequences_set, decoy_sequence, c)
        end

        # Build decoys for each variant in this group using the same mapping
        for idx in idxs
            target_entry = target_fasta_entries[idx]

            adjusted_structural_mods = adjust_mod_positions(
                get_structural_mods(target_entry),
                positions_copy,
                seq_length
            )
            adjusted_structural_mods = merge_fixed_mods(
                decoy_sequence,
                adjusted_structural_mods,
                fixed_mod_names,
            )
            adjusted_isotopic_mods = adjust_mod_positions(
                get_isotopic_mods(target_entry),
                positions_copy,
                seq_length
            )

            push!(decoy_entries, FastaEntry(
                get_id(target_entry),
                get_description(target_entry),
                get_gene(target_entry),
                get_protein(target_entry),
                get_organism(target_entry),
                get_proteome(target_entry),
                decoy_sequence,
                get_start_idx(target_entry),
                adjusted_structural_mods,
                adjusted_isotopic_mods,
                get_charge(target_entry),
                get_base_target_id(target_entry),
                get_base_pep_id(target_entry),
                get_entrapment_pair_id(target_entry),
                true
            ))
        end

        if !isnothing(pbar)
            update(pbar)
        end
    end

    if log_progress
        @user_info "Decoy generation complete: $(length(decoy_entries)) decoy entries created, exhausted_groups=$exhausted_groups, shuffle_fallbacks=$fallback_to_shuffle_count, mutation_fallbacks=$mutation_fallback_count"
    end

    # Report statistics if using reverse method
    #=
    if decoy_method == "reverse" && total_groups > 0
        if fallback_to_shuffle_count > 0
            @user_warn "Decoy generation (GROUPED) stats for REVERSE:"
            @user_warn "  Total base sequences attempted: $total_groups"
            @user_warn "  Reverse duplicates: $fallback_to_shuffle_count"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_groups, digits=1))%"
        else
            @user_info "Successfully reversed all $total_groups base sequences without duplicates"
        end
    end
    =#
    # General summary
    return sort(vcat(target_fasta_entries, decoy_entries), by = x -> get_sequence(x))
end

"""
    combine_shared_peptides(peptides::Vector{FastaEntry})::Vector{FastaEntry}

Combines entries that share identical peptide sequences by concatenating their protein accessions.

# Parameters
- `peptides::Vector{FastaEntry}`: Vector of peptide entries that may contain duplicates

# Returns
- `Vector{FastaEntry}`: Vector of unique peptide entries with concatenated metadata

# Details
This function:
1. Identifies peptides with identical sequences (considering I/L as equivalent)
2. For shared peptides, combines their protein accessions with semicolon separators
3. Also combines proteome identifiers and descriptions if multiple exist
4. Preserves other metadata from the first encountered instance of each sequence
5. Returns a vector containing only unique peptide sequences

# Examples
```julia
# Original entries with shared sequence
entries = [
    FastaEntry("P1", "desc1", "geneA", "protA", "human", "human", "PEPTIDE", 1, missing, missing, 0, 0, 0, 0, 0, 0, false),
    FastaEntry("P2", "desc2", "geneB", "protB", "human", "human", "PEPTIDE", 2, missing, missing, 0, 0, 0, 0, 0, 0, false),
    FastaEntry("P3", "desc3", "geneC", "protC", "human", "human", "UNIQUE", 3, missing, missing, 0, 0, 0, 0, 0, 0, false)
]

# Combine shared peptides
combined = combine_shared_peptides(entries)
# Results in 2 entries:
# 1. FastaEntry("P1;P2", "desc1;desc2", "geneA;geneB", "protA;protB", "human;human", "human;human", "PEPTIDE", 1, missing, missing, 0, 0, 0, 0, 0, 0, false)
# 2. FastaEntry("P3", "desc3", "geneC", "protC", "human", "human", "UNIQUE", 3, missing, missing, 0, 0, 0, 0, 0, 0, false)
```

# Notes
This function helps handle peptides that map to multiple proteins while maintaining
unique sequences in the library. It's particularly useful in bottom-up proteomics
where shared peptides are common.
"""
function combine_shared_peptides(peptides::Vector{FastaEntry})
    seq_to_fasta_entry = Dictionary{String, FastaEntry}()
    n = 0
    a = 0
    base_pep_id = one(UInt32)
    for peptide in peptides
        sequence = get_sequence(peptide)
        sequence_il_equiv = replace(sequence, 'I' => 'L')
        if haskey(seq_to_fasta_entry, sequence_il_equiv)
            a += 1
            fasta_entry = seq_to_fasta_entry[sequence_il_equiv]
            accession = get_id(peptide)*";"*get_id(fasta_entry)
            proteome = get_proteome(peptide)*";"*get_proteome(fasta_entry)
            description = get_description(peptide)*";"*get_description(fasta_entry)
            gene = get_gene(peptide)*";"*get_gene(fasta_entry)
            protein = get_protein(peptide)*";"*get_protein(fasta_entry)
            organism = get_organism(peptide)*";"*get_organism(fasta_entry)
            seq_to_fasta_entry[sequence_il_equiv] = FastaEntry(
                                                        accession,
                                                        description,
                                                        gene,
                                                        protein,
                                                        organism,
                                                        proteome,
                                                        get_sequence(fasta_entry),
                                                        get_start_idx(fasta_entry),
                                                        get_structural_mods(fasta_entry),
                                                        get_isotopic_mods(fasta_entry),
                                                        get_charge(fasta_entry),
                                                        get_base_target_id(fasta_entry), # preserve base_target_id
                                                        base_pep_id,
                                                        get_entrapment_pair_id(fasta_entry), 
                                                        is_decoy(fasta_entry)
                                                        )
            base_pep_id += one(UInt32)
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

function assign_base_pep_ids!(fasta_entries::Vector{FastaEntry})
    """
    Assign sequential base_pep_id values starting from 1.
    Called after add_mods to identify unique peptides (sequence + modifications).
    Each peptide variant gets a unique base_pep_id for tracking through charge variants.
    
    Returns:
    - Int: Number of entries processed
    """
    
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        
        # Create new FastaEntry with sequential base_pep_id
        fasta_entries[i] = FastaEntry(
            get_id(entry),
            get_description(entry),
            get_gene(entry),
            get_protein(entry),
            get_organism(entry),
            get_proteome(entry),
            get_sequence(entry),
            get_start_idx(entry),
            get_structural_mods(entry),
            get_isotopic_mods(entry),
            get_charge(entry),
            get_base_target_id(entry), # preserve base_target_id
            UInt32(i),               # base_pep_id - sequential assignment
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(fasta_entries)
end


function assign_base_target_ids!(fasta_entries::Vector{FastaEntry})
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        # Create new FastaEntry with assigned base_target_id
        fasta_entries[i] = FastaEntry(
            get_id(entry),
            get_description(entry),
            get_gene(entry),
            get_protein(entry),
            get_organism(entry),
            get_proteome(entry),
            get_sequence(entry),
            get_start_idx(entry),
            get_structural_mods(entry),
            get_isotopic_mods(entry),
            get_charge(entry),
            UInt32(i),   # assign grouped base_target_id
            get_base_pep_id(entry),    # preserve existing base_pep_id
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(fasta_entries)
end
