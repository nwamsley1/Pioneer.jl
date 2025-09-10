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

# Essential modification functions needed for BuildSpecLib
# Extracted from parse_mods.jl to avoid BasicEmpiricalLibrary dependency

using Combinatorics: combinations

"""
    getFixedMods!(
        fixed_mods::Vector{PeptideMod},
        mod_matches::Base.RegexMatchIterator,
        mod_name::String)

Add fixed modifications to the provided vector based on regex matches.

# Parameters
- `fixed_mods::Vector{PeptideMod}`: Vector to store fixed modifications
- `mod_matches::Base.RegexMatchIterator`: Iterator of regex matches for modification sites
- `mod_name::String`: Name of the modification to apply

# Effects
Adds `PeptideMod` objects to the `fixed_mods` vector for each regex match

# Returns
`nothing` - Modifies `fixed_mods` in place
"""
function getFixedMods!(
                        fixed_mods::Vector{PeptideMod},
                        mod_matches::Base.RegexMatchIterator,
                        mod_name::String)
    for mod_match in mod_matches
        index = UInt8(mod_match.offset)
        aa = mod_match.match
        push!(fixed_mods, PeptideMod(index, first(aa), mod_name))
    end
    return nothing
end

"""
    matchVarMods(sequence::String, var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})

Find all potential variable modification sites in a peptide sequence based on regex patterns.

# Arguments
- `sequence::String`: Peptide sequence to search for modification sites
- `var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}`: Vector of regex patterns and modification names

# Returns
- `Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}}`: Vector of matches with positions and names

# Examples
```julia
sequence = "PEPTIDE"
var_mods = [(p=r"M", r="Oxidation"), (p=r"[ST]", r="Phospho")]
matches = matchVarMods(sequence, var_mods)
# Returns matches for M (Oxidation) and S/T (Phospho) positions
```

# Notes
This function identifies all potential modification sites without applying any
combinatorial constraints. Use `countVarModCombinations` and `fillVarModStrings!`
to generate actual modification combinations.
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
- `Int`: Total number of modification combinations (including unmodified)

# Examples
```julia
sequence = "METHAMPHETAMINE"
var_mods = [(p=r"M", r="Unimod:35")]
matches = matchVarMods(sequence, var_mods)
n_combinations = countVarModCombinations(matches, 2)
# Returns 4: one unmodified sequence + two single mods + one double mod
```
"""
function countVarModCombinations(var_mod_matches::Vector{NamedTuple{(:regex_match, :name), Tuple{RegexMatch, String}}},
    max_var_mods::Int)
    n_var_mods = length(var_mod_matches)
    n_var_mod_combinations = 0
    for n_mods in 0:min(max_var_mods, n_var_mods)
        n_var_mod_combinations += binomial(n_var_mods, n_mods)
    end
    return n_var_mod_combinations
end

