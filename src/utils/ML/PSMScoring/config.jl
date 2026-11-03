# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#=
Unified configuration for PSM scoring algorithm.
Combines all trait types into a single configuration.
=#

"""
    ScoringConfig{M,T,F}

Complete configuration for PSM scoring algorithm.
Combines all trait types into a single parameterized configuration.

# Type Parameters
- `M <: PSMScoringModel`: ML model (LightGBMScorer, ProbitScorer)
- `T <: TrainingDataStrategy`: Training data selection
- `F <: FeatureSelectionStrategy`: Feature selection method
"""
struct ScoringConfig{M<:PSMScoringModel, T<:TrainingDataStrategy,
                     F<:FeatureSelectionStrategy}
    model::M
    training_data::T
    feature_selection::F
    num_rounds::Int
end

"""
    probit_scoring_config(; features=Symbol[])

Create probit regression scoring configuration.
"""
function probit_scoring_config(; features::Vector{Symbol} = Symbol[])
    return ScoringConfig(
        ProbitScorer(),
        AllDataSelection(),
        StaticFeatureSelection(features),
        30
    )
end

# Note: build_scoring_config function is defined in ScoringSearch/model_config.jl
# because it depends on ModelConfig which is defined there
