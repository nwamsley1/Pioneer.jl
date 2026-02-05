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

#=
Unified configuration for PSM scoring algorithm.
Combines all trait types into a single configuration.
=#

"""
    ScoringConfig{M,P,T,F,I,B}

Complete configuration for PSM scoring algorithm.
Combines all trait types into a single parameterized configuration.

# Type Parameters
- `M <: PSMScoringModel`: ML model (LightGBMScorer, ProbitScorer)
- `P <: PairingStrategy`: Pairing method (RandomPairing, NoPairing)
- `T <: TrainingDataStrategy`: Training data selection
- `F <: FeatureSelectionStrategy`: Feature selection method
- `I <: IterationScheme`: Iteration configuration
- `B <: MBRUpdateStrategy`: MBR feature computation
"""
struct ScoringConfig{M<:PSMScoringModel, P<:PairingStrategy, T<:TrainingDataStrategy,
                     F<:FeatureSelectionStrategy, I<:IterationScheme, B<:MBRUpdateStrategy}
    model::M
    pairing::P
    training_data::T
    feature_selection::F
    iteration_scheme::I
    mbr_update::B
end

"""
    default_scoring_config(; match_between_runs=true, features=Symbol[], hyperparams=Dict())

Create default LightGBM scoring configuration.

# Arguments
- `match_between_runs::Bool`: Whether to enable MBR features
- `features::Vector{Symbol}`: Base feature set
- `hyperparams::Dict{Symbol,Any}`: LightGBM hyperparameters

# Returns
- `ScoringConfig`: Configured scoring configuration
"""
function default_scoring_config(;
    match_between_runs::Bool = true,
    features::Vector{Symbol} = Symbol[],
    hyperparams::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    max_q_value::Float32 = 0.01f0,
    min_pep_threshold::Float32 = 0.90f0
)
    # Default MBR features
    mbr_features = match_between_runs ? [
        :MBR_max_pair_prob,
        :MBR_log2_weight_ratio,
        :MBR_log2_explained_ratio,
        :MBR_rv_coefficient,
        :MBR_is_missing
    ] : Symbol[]

    # Separate base features from MBR features
    base_features = [f for f in features if !startswith(String(f), "MBR_")]

    # Get iteration scheme from hyperparams
    iter_scheme = get(hyperparams, :iter_scheme, [100, 200, 200])

    return ScoringConfig(
        LightGBMScorer(hyperparams),
        RandomPairing(),
        QValueNegativeMining(max_q_value, min_pep_threshold),
        IterativeFeatureSelection(base_features, mbr_features, length(iter_scheme)),
        FixedIterationScheme(iter_scheme),
        match_between_runs ? PairBasedMBR(max_q_value) : NoMBR()
    )
end

"""
    probit_scoring_config(; features=Symbol[])

Create probit regression scoring configuration (no MBR support).
"""
function probit_scoring_config(; features::Vector{Symbol} = Symbol[])
    return ScoringConfig(
        ProbitScorer(),
        NoPairing(),
        AllDataSelection(),
        StaticFeatureSelection(features),
        SinglePassScheme(),
        NoMBR()
    )
end

# Note: build_scoring_config function is defined in ScoringSearch/model_config.jl
# because it depends on ModelConfig which is defined there
