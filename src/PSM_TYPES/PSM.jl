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
export PSM