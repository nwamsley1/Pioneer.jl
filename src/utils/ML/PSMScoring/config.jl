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
