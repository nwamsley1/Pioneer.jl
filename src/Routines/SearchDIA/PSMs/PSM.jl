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
    PSM

    Every concrete instance of type `T` that inherits from `PSM` (peptide spectrum match) should implement the following methods:

- ScoreFragmentMatches!(results::UnorderedDictionary{UInt32, `T`}, matches::Vector{FragmentMatch})). 
- ModifyFeatures!(score::`T`, transition::Transition, mass::Union{Missing, Float32}, intensity::Union{Missing, Float32})
- makePSMsDict((::`T`) - Creats the score dictionary 
- Score!(PSMs_dict::Dict, unscored_PSMs::Vector{`T`}) - Calculates features for each PSM  in Vector{`T`} and adds them to the score dictionary. (See examples)

### Examples
PSM_dict =  makeScoreDict(`T`)
Dict(
                :FeatureA => [],
                :FeatureB => [],
                ...
                )
Score!(PSM_dict, `Vector{T<:PSMFeatures}`)
            Dict(
                :FeatureA => [a, a, ...],
                :FeatureB => [b, b, ...],
                ...
                )
"""
abstract type PSM end