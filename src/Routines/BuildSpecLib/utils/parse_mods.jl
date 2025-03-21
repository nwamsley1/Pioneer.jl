"""
    matchVarMods(sequence::String, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})

Find all positions in a peptide sequence where variable modifications can be applied based on specified patterns.

# Arguments
- `sequence::String`: The peptide sequence to search for modification sites
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Vector of modification patterns and their corresponding names
    - `:p`: Regex pattern that matches amino acids to be modified
    - `:r`: String identifier for the modification (e.g., "Unimod:35")

# Returns
- `Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}`: Vector of matches and their associated modification names
    - `:regex_match`: RegexMatch object containing the match position and matched text
    - `:name`: String identifier of the modification

# Examples
```julia
sequence = "PEPMTIDME"
var_mods = [(p=r"M", r="Unimod:35")]
matches = matchVarMods(sequence, var_mods)
# Returns matches for M at positions 4 and 8 with Unimod:35 modification

# Multiple modification patterns
var_mods = [(p=r"M", r="Unimod:35"), (p=r"[ST]", r="Unimod:21")]
matches = matchVarMods(sequence, var_mods)
# Returns matches for M (Unimod:35) and S/T (Unimod:21) positions
```
"""
function matchVarMods(sequence::String, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    var_mod_matches = Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}()
    for mod in var_mods
        for mod_match in eachmatch(mod[:p], sequence)
            push!(var_mod_matches, (regex_match=mod_match, name=mod[:r]))
        end
    end
    return var_mod_matches
end

"""
    countVarModCombinations(var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                           max_var_mods::Int)

Calculate the total number of possible variable modification combinations for a sequence, including the unmodified version.

# Arguments
- `var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}`: Vector of modification matches from matchVarMods()
- `max_var_mods::Int`: Maximum number of variable modifications allowed per sequence

# Returns
- `Int`: Total number of possible modification combinations, including the unmodified sequence

# Examples
```julia
sequence = "PEPMTIDME"
var_mods = [(p=r"M", r="Unimod:35")]
matches = matchVarMods(sequence, var_mods)
n_combinations = countVarModCombinations(matches, 2)
# Returns 4: one unmodified sequence + two single mods + one double mod
```

Note: The function uses the binomial coefficient to calculate combinations of modifications 
up to max_var_mods and adds 1 for the unmodified sequence.
"""
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

"""
    fillVarModStrings!(var_mod_strings::Vector{String},
                       var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
                       fixed_mods_string::String,
                       max_var_mods::Int)

Fill a pre-allocated vector with all possible combinations of variable modifications, including fixed modifications.

# Arguments
- `var_mod_strings::Vector{String}`: Pre-allocated vector to store modification strings
- `var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}`: Vector of modification matches from matchVarMods()
- `fixed_mods_string::String`: String containing any fixed modifications in format "(index,aa,mod_name)"
- `max_var_mods::Int`: Maximum number of variable modifications allowed per sequence

# Output Format
Each modification String follows the format: "(index,amino_acid,modification_name)", concatenated for multiple modifications.

# Examples
```julia
sequence = "PEPMTIDME"
var_mods = [(p=r"M", r="Unimod:35")]
matches = matchVarMods(sequence, var_mods)
n_combinations = countVarModCombinations(matches, 2)
var_mods_string = Vector{String}(undef, n_combinations)
fillVarModStrings!(var_mods_string, matches, "", 2)
# Results in:
# ["(4,M,Unimod:35)",
#  "(8,M,Unimod:35)",
#  "(4,M,Unimod:35)(8,M,Unimod:35)",
#  ""]
```

Note: The function uses Combinatorics.combinations to generate all possible modification combinations.
Requires `using Combinatorics`.
"""
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

"""
    parseEmpiricalLibraryMods(sequence::String)

Parse modifications embedded in a sequence String and return them in a standardized format.
Modifications in the sequence are in the format AA(mod_name).

# Arguments
- `sequence::String`: Peptide sequence with embedded modifications

# Returns
- `String`: Modified String containing modifications in format "(position,amino_acid,mod_name)"

# Examples
```julia
# Test case 1
sequence = "I(exTag1)LSISADI(Unimod:35)ETIGEILK(exTag1)"
result = parseEmpiricalLibraryMods(sequence)
@assert result == "(1,I,exTag1)(8,I,Unimod:35)(16,K,exTag1)" "Test 1 failed"

# Test case 2
sequence = "n(exTag1)PEPTIDE(exTag1)"
result = parseEmpiricalLibraryMods(sequence)
@assert result == "(1,n,exTag1)(7,E,exTag1)" "Test 2 failed"

# Test case 3
sequence = "n(exTag1)PEPTIDEc(exTag1)"
result = parseEmpiricalLibraryMods(sequence)
@assert result == "(1,n,exTag1)(7,c,exTag1)" "Test 3 failed"

println("All tests passed!")
```
"""
function parseEmpiricalLibraryMods(sequence::String)
    # Regex to find modifications in sequence
    mod_regex = r"([A-Z]|[nc])\(([^)]+)\)"
    
    # Check for terminal markers
    has_n_term = startswith(sequence, "n")
    
    # Track position adjustments due to removing modification text
    pos_adjust = 0
    if has_n_term
        pos_adjust += 1  # Start position counting after 'n'
    end
    mods = String[]
    
    # Find all modifications
    for m in eachmatch(mod_regex, sequence)
        aa = m.captures[1]
        mod_name = m.captures[2]
        if 'n' ∈ aa
            push!(mods, "(1,n,$mod_name)")
        elseif 'c' ∈ aa
            pos = m.offset - pos_adjust - 1
            push!(mods, "($pos,c,$mod_name)")
        else
            pos = m.offset - pos_adjust
            push!(mods, "($pos,$aa,$mod_name)")
        end
        # Adjust for removal of (mod_name)
        pos_adjust += length(m.match) - 1
    end
    
    # Sort by position and join
    sort!(mods, by=x->parse(Int, match(r"\((\d+)", x).captures[1]))
    return join(mods)
end

"""
    getMassOffset(start_idx::Integer, stop_idx::Integer, mod_masses::Vector{Float32})

Calculate the total modification mass for a sequence fragment.

# Arguments
- `start_idx::Integer`: Starting position in the sequence (inclusive)
- `stop_idx::Integer`: Ending position in the sequence (inclusive)
- `mod_masses::Vector{Float32}`: Array containing modification masses at each position

# Returns
- `Float32`: Total mass of modifications in the specified range

# Examples
```julia
mod_masses = [229.163f0, 0.0f0, 0.0f0, 0.0f0, 229.163f0]
mass_offset = getMassOffset(1, 3, mod_masses)  # Returns mass for positions 1-3
```
"""
function getMassOffset(start_idx::Integer, stop_idx::Integer, mod_masses::Vector{Float32})
    mass_offset = 0.0f0
    for i in start_idx:stop_idx
        mass_offset += mod_masses[i]
    end
    return mass_offset
end

"""
    addIsoMods(isotopic_mods::String, 
               structural_mods::String, 
               mod_key::String,
               mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}})

Creates a new String combining existing isotopic modifications with additional isotopic channel information 
for specified structural modifications. The function identifies all instances of the target modification in 
the structural modifications String and adds corresponding isotopic channel annotations.

# Arguments
- `isotopic_mods::String`: Existing isotopic modifications in format "(mod_index, channel)", where mod_index refers 
   to the position of the modification in the structural_mods String
- `structural_mods::String`: Structural modifications in format "(position,amino_acid,mod_name)"
- `mod_key::String`: Target modification to add isotopic channel for (e.g., "exTag1")
- `mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}}`: Isotopic channel information
    - `:channel`: String identifier for the isotopic channel (e.g., "d0")
    - `:mass`: Mass shift associated with the channel (unused in String generation)

# Returns
- `String`: New String containing all isotopic modifications in format "(mod_index, channel)", where mod_index 
   refers to the order of appearance in structural_mods (1-based indexing)

# Examples
```julia
# Basic case with structural mod indices
structural_mods = "(1,n,exTag1)(5,L,exTag3)(16,c,exTag1)"
isotopic_mods = ""
result = addIsoMods(isotopic_mods, structural_mods, "exTag1", (channel="d0", mass=0.0f0))
# Returns "(1, d0)(3, d0)"
# Here:
#   - 1 refers to first modification "(1,n,exTag1)"
#   - 3 refers to third modification "(16,c,exTag1)"

# With existing isotopic mods
structural_mods = "(1,n,exTag1)(5,L,exTag3)(16,c,exTag1)"
isotopic_mods = "(2, d4)"  # Refers to second structural mod "(5,L,exTag3)"
result = addIsoMods(isotopic_mods, structural_mods, "exTag1", (channel="d0", mass=0.0f0))
# Returns "(1, d0)(2, d4)(3, d0)"
# Here:
#   - 1, 3 refer to first and third exTag1 modifications
#   - 2 refers to existing isotopic mod on second structural modification
```

Note: Mod indices in the output are sorted numerically for consistent ordering. If no instances 
of the target modification are found, returns the original isotopic_mods String unchanged.
The indices in the output refer to the order of appearance in structural_mods (1-based), not 
to sequence positions.
"""
function addIsoMods(isotopic_mods::String, 
                    structural_mods::String, 
                    mod_key::String,
                    mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}})
    
    # Parse structural mods to find indices of target modification
    struct_mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    mod_indices = Int[]  # Store the index of the modification, not its position
    
    # Find all instances where the target modification occurs
    for (mod_idx, m) in enumerate(eachmatch(struct_mod_regex, structural_mods))
        mod_name = m.captures[3]
        if mod_name == mod_key
            push!(mod_indices, mod_idx)
        end
    end
    
    # If no target modifications found, return original isotopic mods
    isempty(mod_indices) && return isotopic_mods
    
    # Parse existing isotopic modifications
    iso_mod_regex = r"\((\d+),\s*([^)]+)\)"
    existing_mods = Dict{Int, String}()
    
    if !isempty(isotopic_mods)
        for m in eachmatch(iso_mod_regex, isotopic_mods)
            idx = parse(Int, m.captures[1])
            channel = m.captures[2]
            existing_mods[idx] = channel
        end
    end
    
    # Add new modifications
    for idx in mod_indices
        existing_mods[idx] = mod_channel.channel
    end
    
    # Create sorted output
    sorted_indices = sort(collect(keys(existing_mods)))
    return join("($idx, $(existing_mods[idx]))" for idx in sorted_indices)
end

function addIsoMods(isotopic_mods::Missing, 
    structural_mods::String, 
    mod_key::String,
    mod_channel::NamedTuple{(:channel, :mass), Tuple{String, Float32}})

    isotopic_mods = addIsoMods("", structural_mods, mod_key, mod_channel)
    if length(isotopic_mods) == 0
        return missing
    else
        return isotopic_mods
    end
end

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

"""
    getRevDecoys!(speclibdf::BasicEmpiricalLibrary)

Generates reversed decoy sequences for a spectral library while maintaining consistency
across fragments of the same precursor. Reverses each unique precursor sequence once
and applies the same reversed sequence to all its fragments.

# Arguments
- `speclibdf::BasicEmpiricalLibrary`: Spectral library containing peptide sequences,
   modifications, and precursor indices

# Effects
- Modifies `speclibdf` in place by appending non-duplicate reversed sequences
- Sets `is_decoy=true` for all appended reversed sequences
- Maintains consistency across fragments of the same precursor
- Discards all fragments of a precursor if its reversed sequence matches any forward sequence

# Run test if this file is executed directly
```julia
println("test")
```
"""
function getRevDecoys!(speclibdf::BasicEmpiricalLibrary)
    # Make a deep copy and sort by precursor_idx
    rev_libdf = deepcopy(speclibdf.libdf)
    sort!(rev_libdf, :precursor_idx)
    
    # Create a Set of forward sequences for O(1) lookup
    forward_seqs = Set(speclibdf.libdf.sequence)
    
    # Track which rows to keep
    keep_rows = trues(size(rev_libdf, 1))
    
    # Track current precursor and its reversed sequence
    current_precursor = nothing
    current_rev_seq = nothing
    current_rev_mods = nothing
    
    # Process rows maintaining precursor consistency
    for (idx, row) in enumerate(eachrow(rev_libdf))
        if row.precursor_idx != current_precursor
            # New precursor encountered - need to generate new reversal
            current_precursor = row.precursor_idx
            
            # Generate reversed sequence
            rev_sequence, rev_structural_mods = reverseSequence(row.sequence, row.structural_mods)
            
            # If reversed sequence exists in forward set, mark all rows of this precursor for removal
            if rev_sequence in forward_seqs
                idx_start = idx
                while idx_start <= size(rev_libdf, 1) && 
                      rev_libdf[idx_start, :precursor_idx] == current_precursor
                    keep_rows[idx_start] = false
                    idx_start += 1
                end
                continue
            end
            
            # Store reversed sequence for this precursor group
            current_rev_seq = rev_sequence
            current_rev_mods = rev_structural_mods
            push!(forward_seqs, rev_sequence)
        end
        
        # Apply current reversed sequence to this row
        row.sequence = current_rev_seq
        row.structural_mods = current_rev_mods
        row.is_decoy = true
    end
    
    # Filter out duplicates and append
    return rev_libdf[keep_rows, :]
    #append!(speclibdf.libdf, rev_libdf[keep_rows, :])
end

"""
    shuffleSequence(sequence::String, structural_mods::String)

Randomly shuffle a peptide sequence (keeping last AA fixed) and update modification positions accordingly.

# Arguments
- `sequence::String`: Original peptide sequence
- `structural_mods::String`: Modifications in format "(pos,aa,mod_name)"

# Returns
- `Tuple{String, String}`: (shuffled sequence, updated modification String)

# Examples
```julia
# Test 1: Shuffle with terminal modifications
Random.seed!(1844)
sequence = "PEPTIDE"
mods = "(1,n,mymod-nterm)(1,P,mymod)(4,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
shuffled_seq, shuffled_mods = shuffleSequence(sequence, mods)

# Test that:
# 1. Seed set so sequence should always be DPTPIEE
@assert shuffled_seq == "DPTPIEE"
@assert shuffled_mods == "(1,n,mymod-nterm)(2,P,mymod)(3,T,Phospho)(7,E,Acetyl)(7,c,mymod-cterm)"
```
"""
function shuffleSequence(sequence::AbstractString, structural_mods::String)
    # Early return if no modifications
    if isempty(structural_mods)
        # Shuffle all but last AA
        return join(Random.shuffle(collect(sequence[1:end-1])))*sequence[end], ""
    end
    
    # Pre-count number of modifications
    n_mods = count("(", structural_mods)
    
    # Pre-allocate arrays
    mod_positions = Vector{Int}(undef, n_mods)
    mod_aas = Vector{Char}(undef, n_mods)
    mod_names = Vector{String}(undef, n_mods)
    
    # Parse modifications
    mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    i = 1
    for m in eachmatch(mod_regex, structural_mods)
        mod_positions[i] = parse(Int, m.captures[1])
        mod_aas[i] = first(m.captures[2])
        mod_names[i] = m.captures[3]
        i += 1
    end
    
    # Create shuffled sequence (keeping last AA in place)
    seq_length = length(sequence)
    all_but_last = collect(sequence[1:end-1])
    shuffle!(all_but_last)
    shuffled_seq = String(all_but_last) * sequence[end]
    
    # Create position mapping from original to shuffled sequence
    pos_map = Dict{Int, Int}()
    for (old_pos, aa) in enumerate(sequence[1:end-1])
        # Find where this AA went in the shuffled sequence
        new_pos = findfirst(==(aa), shuffled_seq)
        pos_map[old_pos] = new_pos
    end
    pos_map[seq_length] = seq_length  # Last position stays the same
    
    # Calculate new positions and create tuples for sorting
    final_mods = Vector{Tuple{Int, String}}(undef, n_mods)
    for i in 1:n_mods
        pos = mod_positions[i]
        aa = mod_aas[i]
        name = mod_names[i]
        
        if aa == 'n'
            final_mods[i] = (1, "(1,n,$name)")
        elseif pos == seq_length
            final_mods[i] = (seq_length, "($pos,$aa,$name)")
        else
            new_pos = pos_map[pos]
            new_aa = shuffled_seq[new_pos]
            final_mods[i] = (new_pos, "($new_pos,$new_aa,$name)")
        end
    end
    
    # Sort by new positions and join
    sort!(final_mods, by=first)
    return shuffled_seq, join(last.(final_mods))
end

"""
    getShuffledEntrapmentSeqs!(speclibdf::BasicEmpiricalLibrary)

Generates shuffled decoy sequences for a spectral library while maintaining consistency
across fragments of the same precursor. Shuffles each unique precursor sequence once
and applies the same shuffled sequence to all its fragments.

# Arguments
- `speclibdf::BasicEmpiricalLibrary`: Spectral library containing peptide sequences,
   modifications, and precursor indices

# Effects
- Modifies `speclibdf` in place by appending non-duplicate shuffled sequences
- Sets `is_decoy=true` for all appended shuffled sequences
- Maintains consistency across fragments of the same precursor
- Discards all fragments of a precursor if shuffling fails after 20 attempts

#Tests 
```julia
# Run test if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    using DataFrames
    
    # Test setup
    struct BasicEmpiricalLibrary
        libdf::DataFrame
    end
    
    # Create test data with multiple fragments per precursor
    test_df = DataFrame(
        sequence = repeat(["PEPTIDE", "PROTEIN"], inner=3),
        structural_mods = repeat(["(2,P,Oxidation)", ""], inner=3),
        precursor_idx = repeat([1, 2], inner=3),
        fragment_idx = repeat(1:3, outer=2),
        is_decoy = falses(6)
    )
    test_lib = BasicEmpiricalLibrary(test_df)
    
    # Generate decoys
    getShuffledEntrapmentSeqs!(test_lib)
    
    # Verify results
    result_df = test_lib.libdf
    
    # Basic checks
    @assert size(result_df, 1) == 12 "Should have original + decoy rows"
    @assert all(result_df[7:end, :is_decoy]) "New rows should be marked as decoys"
    
    # Check precursor consistency
    for precursor_idx in unique(result_df.precursor_idx)
        decoy_rows = result_df[result_df.is_decoy .& (result_df.precursor_idx .== precursor_idx), :]
        if !isempty(decoy_rows)
            # All fragments of same precursor should have same sequence
            @assert length(unique(decoy_rows.sequence)) == 1 "All fragments of precursor precursor_idx should have same sequence"
            # Sequence should be different from original
            orig_seq = first(result_df[.!result_df.is_decoy .& (result_df.precursor_idx .== precursor_idx), :sequence])
            @assert first(decoy_rows.sequence) != orig_seq "Decoy sequence should differ from original"
        end
    end
    
    println("All tests passed!")
end
```
"""
function getShuffledEntrapmentSeqs!(speclibdf::BasicEmpiricalLibrary, entrapment_group_id::Integer)
    # Make a deep copy and sort by precursor_idx
    shuffle_libdf = deepcopy(speclibdf.libdf)
    sort!(shuffle_libdf, :precursor_idx)
    shuffle_libdf[!,:entrapment_group_id] .= UInt8(entrapment_group_id)
    # Create a Set of forward sequences for O(1) lookup
    forward_seqs = Set(speclibdf.libdf.sequence)
    
    # Track which rows to keep
    keep_rows = trues(size(shuffle_libdf, 1))
    
    # Track current precursor and its shuffled sequence
    current_precursor = nothing
    current_shuffled_seq = nothing
    current_shuffled_mods = nothing
    
    # Process rows maintaining precursor consistency
    for (idx, row) in enumerate(eachrow(shuffle_libdf))
        if row.precursor_idx != current_precursor
            # New precursor encountered - need to generate new shuffle
            current_precursor = row.precursor_idx
            
            # Try to find a valid shuffle
            found_valid_shuffle = false
            shuffle_attempts = 0
            while shuffle_attempts < 20 && !found_valid_shuffle
                shuffled_sequence, shuffled_mods = shuffleSequence(row.sequence, row.structural_mods)
                if !(shuffled_sequence in forward_seqs)
                    current_shuffled_seq = shuffled_sequence
                    current_shuffled_mods = shuffled_mods
                    push!(forward_seqs, shuffled_sequence)
                    found_valid_shuffle = true
                end
                shuffle_attempts += 1
            end
            
            # If no valid shuffle found, mark all rows of this precursor for removal
            if !found_valid_shuffle
                idx_start = idx
                while idx_start <= size(shuffle_libdf, 1) && 
                      shuffle_libdf[idx_start, :precursor_idx] == current_precursor
                    keep_rows[idx_start] = false
                    idx_start += 1
                end
                continue
            end
        end
        
        # Apply current shuffled sequence to this row
        row.sequence = current_shuffled_seq
        row.structural_mods = current_shuffled_mods
        row.is_decoy = true
    end
    
    # Filter out failed shuffles and append
    append!(speclibdf.libdf, shuffle_libdf[keep_rows, :])
end
