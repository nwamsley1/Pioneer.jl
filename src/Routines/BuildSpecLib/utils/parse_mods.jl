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
    fillVarModStrings!(
        all_mods::Vector{Vector{PeptideMod}},
        var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
        fixed_mods::Vector{PeptideMod},
        max_var_mods::Int)

Fill a vector with all valid combinations of variable modifications up to a specified maximum.

# Parameters
- `all_mods::Vector{Vector{PeptideMod}}`: Pre-allocated vector to store all modification combinations
- `var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}`: Potential variable modification sites
- `fixed_mods::Vector{PeptideMod}`: Fixed modifications to include in all combinations
- `max_var_mods::Int`: Maximum number of variable modifications to apply to a single peptide

# Effects
Fills `all_mods` with different combinations of modifications, each element containing 
a complete set of modifications (both fixed and variable) for a peptide variant

# Returns
`nothing` - Modifies `all_mods` in place

# Notes
The function generates combinations using the `combinations` function and ensures
that each combination includes all fixed modifications. The combinations are
sorted for consistent output.
"""
function fillVarModStrings!(
                            all_mods::Vector{Vector{PeptideMod}},
                            var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                            fixed_mods::Vector{PeptideMod},
                            max_var_mods::Int,
                            min_var_mods::Int = 0
                            )
    i = 1
    #from max(1,min_var_mods) to length(var_mod_matches) variable mods
    for N in max(1, min_var_mods):min(max_var_mods, length(var_mod_matches))
        #Each combination of 'N' var mods
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
    #Version with 0 variable mods — only when min_var_mods==0 (else phospho-only)
    if min_var_mods == 0
        all_mods[i] = fixed_mods
        i += 1
    end
    sort!(all_mods)
    return nothing
end

# Empirical-library parsing helpers are disabled with ParseSpecLib /
# BasicEmpiricalLibrary. BuildSpecLib's FASTA path does not use them.

"""
    getIsoModMasses!(iso_mod_masses::Vector{Float32},
                     structural_mods::String,
                     isotopic_mods::String,
                     iso_mods_dict::Dict{String, Dict{String, Float32}})

Maps isotopic modification masses to their positions in a peptide sequence by filling an array of masses.
Uses structural modifications to identify modification types and isotopic modifications to determine channels.

# Arguments
- `iso_mod_masses::Vector{Float32}`: Pre-allocated vector to store isotopic modification masses
- `structural_mods::String`: String containing structural modifications in format "(position,aa,mod_name)"
- `isotopic_mods::String`: String containing isotopic channels in format "(mod_index, channel)" where mod_index 
   refers to the position of the modification in structural_mods (1-based)
- `iso_mods_dict::Dict{String, Dict{String, Float32}}`: Nested dictionary mapping modification names and channels to masses

# Effects
- Fills `iso_mod_masses` with isotopic modification masses at sequence positions
- All other positions remain 0.0

# Example
```julia
iso_mods_dict = Dict(
    "exTag1" => Dict("d0" => 0.0f0, "d4" => 4.0f0),
    "exTag3" => Dict("d0" => 0.0f0, "d4" => 4.0f0)
)
# Third modification in structural_mods gets d4 channel
structural_mods = "(1,I,exTag1)(5,L,exTag3)(7,K,exTag1)"
isotopic_mods = "(3, d4)"  # Refers to third mod "(7,K,exTag1)"
iso_mod_masses = zeros(Float32, 255)
getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
# Results in iso_mod_masses[7] = 4.

iso_mods_dict = Dict(
    "exTag1" => Dict("d0" => 0.0f0, "d4" => 4.0f0),
    "exTag3" => Dict("d0" => 0.0f0, "d4" => 4.0f0)
)

# Test 1: Basic indexing
structural_mods = "(1,I,exTag1)(5,L,exTag3)(7,K,exTag1)"
isotopic_mods = "(3, d4)"  # Third mod gets d4
iso_mod_masses = zeros(Float32, 255)
getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
@assert iso_mod_masses[7] == 4.0f0 "Failed: Third mod should have mass 4.0"
@assert iso_mod_masses[1] == 0.0f0 "Failed: First mod should have mass 0.0"

# Test 2: Multiple mods
structural_mods = "(1,n,exTag1)(5,L,exTag3)(16,c,exTag1)"
isotopic_mods = "(1, d0)(2, d4)"  # First mod gets d0, second gets d4
iso_mod_masses = zeros(Float32, 255)
getIsoModMasses!(iso_mod_masses, structural_mods, isotopic_mods, iso_mods_dict)
@assert iso_mod_masses[1] == 0.0f0 "Failed: First mod should have mass 0.0"
@assert iso_mod_masses[5] == 4.0f0 "Failed: Second mod should have mass 4.0"

println("All tests passed!")
```
"""
function getIsoModMasses!(iso_mod_masses::Vector{Float32},
                         structural_mods::String,
                         isotopic_mods::String,
                         iso_mods_dict::Dict{String, Dict{String, Float32}})
    
    # Reset the masses array
    fill!(iso_mod_masses, 0.0f0)
    
    # Return early if no mods
    (isempty(structural_mods) || isempty(isotopic_mods)) && return
    
    # Parse structural mods into sequence of (position, mod_name) pairs
    struct_mod_regex = r"\((\d+),(?:[A-Z]|[nc]),([^,\)]+)\)"
    structural_mods_seq = Vector{Tuple{Int, String}}()
    
    for m in eachmatch(struct_mod_regex, structural_mods)
        pos = parse(Int, m.captures[1])
        mod_name = m.captures[2]
        if haskey(iso_mods_dict, mod_name)
            push!(structural_mods_seq, (pos, mod_name))
        end
    end
    
    # Parse isotopic mods and apply masses
    iso_mod_regex = r"\((\d+),\s*([^)]+)\)"
    
    for m in eachmatch(iso_mod_regex, isotopic_mods)
        mod_idx = parse(Int, m.captures[1])
        channel = strip(m.captures[2])
        
        # Check if this index exists in our structural mods sequence
        if 1 <= mod_idx <= length(structural_mods_seq)
            pos, mod_name = structural_mods_seq[mod_idx]
            if haskey(iso_mods_dict, mod_name) && haskey(iso_mods_dict[mod_name], channel)
                iso_mod_masses[pos] = iso_mods_dict[mod_name][channel]
            end
        end
    end
end

"""
    get_aa_masses!(aa_masses::Vector{T}, sequence::String) where {T<:AbstractFloat}

Fill a pre-allocated vector with amino acid masses for each position in the sequence.

# Arguments
- `aa_masses::Vector{T}`: Pre-allocated vector to store amino acid masses
- `sequence::String`: Peptide sequence

# Effects
Fills `aa_masses` with the mass of each amino acid at its corresponding position

# Examples
```julia
# Test 1: Simple tripeptide
sequence = "PAK"
aa_masses = zeros(Float32, 3)
get_aa_masses!(aa_masses, sequence)
aa_masses .≈ [97.05276,    # P
  71.03711,    # A
  128.09496]   # K

# Test 3: Sequence with all 20 standard amino acids
sequence = "ACDEFGHIKLMNPQRSTVWY"
aa_masses = zeros(Float32, 20)
get_aa_masses!(aa_masses, sequence)
@assert aa_masses[1] ≈ 71.03711f0  # A
@assert aa_masses[8] ≈ 113.08406f0 # I
@assert aa_masses[20] ≈ 163.06333f0 # Y
```

Note: All masses are monoisotopic masses in Daltons (Da).
See AA_to_mass dictionary for the complete mass mapping.
"""
function get_aa_masses!(aa_masses::Vector{T}, sequence::AbstractString) where {T<:AbstractFloat}
    fill!(aa_masses, zero(T))
    for (i, aa) in enumerate(sequence)
        aa_masses[i] = convert(T, AA_to_mass[aa])
    end
end

"""
    get_structural_mod_masses!(mod_masses::Vector{T}, 
                             structural_mods::String,
                             mod_to_mass::Dict{String, T}) where {T<:AbstractFloat}

Fill a pre-allocated vector with structural modification masses.

# Arguments
- `mod_masses::Vector{T}`: Pre-allocated vector to store modification masses
- `structural_mods::String`: String containing modifications in format "(index,aa,mod_name)"
- `mod_to_mass::Dict{String, T}`: Dictionary mapping modification names to masses

# Example
```julia
sequence = "PEPTIDE"
structural_mods = "(1,P,Phospho)(7,E,Acetyl)"
mod_masses = zeros(Float32, length(sequence))
mod_to_mass = Dict("Phospho" => 79.966331f0, "Acetyl" => 42.010565f0)
get_structural_mod_masses!(mod_masses, structural_mods, mod_to_mass)
```
"""
function get_structural_mod_masses!(mod_masses::Vector{T}, 
                                  structural_mods::AbstractString,
                                  mod_to_mass::Dict{String, T}) where {T<:AbstractFloat}
    fill!(mod_masses, zero(T))
    
    # If no mods, return early
    isempty(structural_mods) && return
    
    # Parse modification entries
    mod_regex = r"\((\d+),([A-Z]|[nc]),([^)]+)\)"
    
    for m in eachmatch(mod_regex, structural_mods)
        index = parse(Int, m.captures[1])
        mod_name = m.captures[3]
        
        if haskey(mod_to_mass, mod_name)
            mod_masses[index] = mod_to_mass[mod_name]
        end
    end
end

function get_fragment_mz(
                        start_idx::Integer, stop_idx::Integer, 
                        base_type::Char, charge::UInt8, 
                        aa_masses::Vector{T},
                        structural_mod_masses::Vector{T},
                        iso_mod_masses::Vector{T})::T where {T<:AbstractFloat}
    #H2O and PROTON are global constants defined in Pioneer.jl 
    # Get modification mass for this fragment range
    frag_mass = zero(T)
    for idx in range(start_idx, stop_idx)
        frag_mass += aa_masses[idx] + structural_mod_masses[idx] + iso_mod_masses[idx]
    end

    # Calculate base mass depending on ion type
    base_mass = if base_type == 'b'
        PROTON * charge  # b-ions just get protons
    else
        PROTON * charge + H2O  # y-ions get water plus protons
    end
    
    # Return final m/z
    return (frag_mass + base_mass) / charge
end

function get_precursor_mz(
                        seq_length::Integer, charge::UInt8, 
                        aa_masses::Vector{T},
                        structural_mod_masses::Vector{T},
                        iso_mod_masses::Vector{T})::T where {T<:AbstractFloat}
    #H2O and PROTON are global constants defined in Pioneer.jl 
    # Get modification mass for this fragment range
    precursor_mass = zero(T)
    charge = convert(T, charge)
    for idx in range(1, seq_length)
        precursor_mass += aa_masses[idx] + structural_mod_masses[idx] + iso_mod_masses[idx]
    end
    
    # Return final m/z
    return (precursor_mass + PROTON*charge + H2O) / charge
end

"""
    get_sulfur_counts!(sulfur_counts::Vector{Int8}, sequence::String, structural_mods::String, mods_to_sulfur_diff::Dict{String, Int8})

Fill a pre-allocated vector with sulfur counts for each position in the sequence,
including both amino acid sulfurs (C, M) and modification-based sulfurs.

# Arguments
- `sulfur_counts::Vector{Int8}`: Pre-allocated vector to store sulfur counts
- `sequence::String`: Peptide sequence
- `structural_mods::String`: String containing modifications in format "(pos,aa,mod_name)"
- `mods_to_sulfur_diff::Dict{String, Int8}`: Dictionary mapping modification names to sulfur count changes

# Effects
Fills `sulfur_counts` with the number of sulfurs at each position
"""
function get_sulfur_counts!(sulfur_counts::Vector{Int8}, 
                          sequence::AbstractString, 
                          structural_mods::AbstractString,
                          mods_to_sulfur_diff::Dict{String, Int8})
    fill!(sulfur_counts, zero(Int8))
    
    # Count sulfurs in sequence
    for (i, aa) in enumerate(sequence)
        sulfur_counts[i] = Int8((aa == 'C') | (aa == 'M'))
    end
    
    # Add sulfurs from modifications
    if !isempty(structural_mods)
        mod_regex = r"\((\d+),([A-Z]|[nc]),([^)]+)\)"
        for m in eachmatch(mod_regex, structural_mods)
            idx = parse(Int, m.captures[1])
            mod_name = m.captures[3]
            if haskey(mods_to_sulfur_diff, mod_name)
                sulfur_counts[idx] += mods_to_sulfur_diff[mod_name]
            end
        end
    end
end

"""
    calculate_reversed_mz!(df::DataFrame, iso_mod_masses::Vector{Float32})

Reverses the peptide sequence (except last AA) and recalculates fragment m/z values.
Terminal modifications stay in place and are included in m/z calculations.
"""
function calculate_mz_and_sulfur_count!(df::DataFrame, 
                    structural_mod_to_mass::Dict{String, T},
                    iso_mods_dict::Dict{String, Dict{String, T}},
                    mods_to_sulfur_diff::Dict{String, Int8}) where {T<:AbstractFloat}
    
    aa_masses = zeros(Float32, 255)
    structural_mod_masses = zeros(Float32, 255)
    iso_mod_masses = zeros(Float32, 255)
    sulfur_counts = zeros(Int8, 255)
    
    for row in eachrow(df)
        seq_length = length(row.sequence)
        base_type = first(row.frag_type)
        
        # Calculate masses
        get_aa_masses!(aa_masses, row.sequence)
        get_structural_mod_masses!(structural_mod_masses, row.structural_mods, structural_mod_to_mass)
        getIsoModMasses!(iso_mod_masses, row.structural_mods, row.isotopic_mods, iso_mods_dict)
        
        # Calculate sulfur counts
        get_sulfur_counts!(sulfur_counts, row.sequence, row.structural_mods, mods_to_sulfur_diff)
        
        # Get fragment indices
        start_idx, stop_idx = get_fragment_indices(
            base_type,
            row.frag_series_number,
            UInt8(seq_length)
        )
        
        # Calculate m/z values
        row.frag_mz = get_fragment_mz(
            start_idx, stop_idx,
            base_type, row.frag_charge,
            aa_masses,
            structural_mod_masses,
            iso_mod_masses
        )
        
        row.prec_mz = get_precursor_mz(
            seq_length, row.prec_charge,
            aa_masses,
            structural_mod_masses,
            iso_mod_masses
        )
        
        # Set sulfur counts
        row.prec_sulfur_count = sum(sulfur_counts[1:seq_length])
        row.frag_sulfur_count = sum(sulfur_counts[start_idx:stop_idx])
    end
end

"""
    reverseSequence(sequence::String, structural_mods::String)

Reverse a peptide sequence while preserving first and last position modifications.
Internal modifications are moved to correspond with their amino acids in the reversed sequence.

# Arguments
- `sequence::String`: Original peptide sequence
- `structural_mods::String`: String containing modifications in format "(pos,aa,mod_name)"

# Returns
- `Tuple{String, String}`: (reversed sequence, updated modification String)

# Examples
```julia
# Test 1: Basic reversal with terminal modifications
sequence = "PEPTIDE"
mods = "(4,T,Phospho)(7,E,Acetyl)"
rev_seq, rev_mods = reverseSequence(sequence, mods)
@assert rev_seq == "DITPEPE"
@assert rev_mods == "(3,T,Phospho)(7,E,Acetyl)"

# Test 2: Reversal with internal modification
sequence = "PEPTIDE"
mods = "(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)"
rev_seq, rev_mods = reverseSequence(sequence, mods)
@assert rev_seq == "DITPEPE"
@assert rev_mods == "(3,T,Phospho)(6,P,mymod)(7,E,Acetyl)"

# Test 2: Reversal with internal modification
sequence = "PEPTIDE"
mods = "(1,n,mymod-nterm)(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
rev_seq, rev_mods = reverseSequence(sequence, mods)
@assert rev_seq  == "DITPEPE"
@assert rev_mods == "(1,n,mymod-nterm)(3,T,Phospho)(6,P,mymod)(7,E,Acetyl)(7,c,mymod-cterm)"

```
"""
function reverseSequence(sequence::AbstractString, structural_mods::String)
    # Early return if no modifications
    if isempty(structural_mods)
        return reverse(sequence[1:end-1])*sequence[end], ""
    end
    
    # Pre-count number of modifications for array pre-allocation
    n_mods = count("(", structural_mods)
    
    # Pre-allocate arrays
    mod_positions = Vector{Int}(undef, n_mods)
    mod_aas = Vector{Char}(undef, n_mods)
    mod_names = Vector{String}(undef, n_mods)
    
    # Parse modifications into parallel arrays
    mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    i = 1
    for m in eachmatch(mod_regex, structural_mods)
        mod_positions[i] = parse(Int, m.captures[1])
        mod_aas[i] = first(m.captures[2])
        mod_names[i] = m.captures[3]
        i += 1
    end
    
    # Create reversed sequence (keeping last AA in place)
    seq_length = length(sequence)
    rev_seq = reverse(sequence[1:end-1])*sequence[end]
    
    # Calculate new positions and create tuples for sorting
    n_final_mods = n_mods
    final_mods = Vector{Tuple{Int, String}}(undef, n_final_mods)
    j = 1
    for i in 1:n_mods
        pos = mod_positions[i]
        aa = mod_aas[i]
        name = mod_names[i]
        
        if aa == 'n'
            final_mods[j] = (1, "(1,n,$name)")
        elseif pos == seq_length
            final_mods[j] = (seq_length, "($pos,$aa,$name)")
        else
            new_pos = seq_length - pos
            new_aa = rev_seq[new_pos]
            final_mods[j] = (new_pos, "($new_pos,$new_aa,$name)")
        end
        j += 1
    end
    
    # Sort by new positions and join
    sort!(final_mods, by=first)
    return rev_seq, join(last.(final_mods))
end

# COMMENTED OUT: Functions below depend on BasicEmpiricalLibrary which is disabled
# Uncomment when BasicEmpiricalLibrary is re-enabled

