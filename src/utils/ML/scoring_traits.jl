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
Trait hierarchy for PSM scoring algorithms.
Each trait represents a single point of abstraction in the scoring pipeline.
=#

#############################################################################
# 1. PSMScoringModel - ML Algorithm Selection
#############################################################################

"""
Abstract type for PSM scoring models.
Concrete types define the ML algorithm used for semi-supervised learning.
"""
abstract type PSMScoringModel end

"""
    LightGBMScorer <: PSMScoringModel

LightGBM gradient boosting classifier for PSM scoring.
"""
struct LightGBMScorer <: PSMScoringModel
    hyperparams::Dict{Symbol, Any}
end

# Convenience constructor with defaults
function LightGBMScorer(;
    feature_fraction::Float64 = 0.5,
    learning_rate::Float64 = 0.05,
    min_data_in_leaf::Int = 500,
    bagging_fraction::Float64 = 0.25,
    min_gain_to_split::Float64 = 0.5,
    max_depth::Int = 10,
    num_leaves::Int = 63
)
    return LightGBMScorer(Dict{Symbol, Any}(
        :feature_fraction => feature_fraction,
        :learning_rate => learning_rate,
        :min_data_in_leaf => min_data_in_leaf,
        :bagging_fraction => bagging_fraction,
        :min_gain_to_split => min_gain_to_split,
        :max_depth => max_depth,
        :num_leaves => num_leaves
    ))
end

"""
    ProbitScorer <: PSMScoringModel

Probit regression classifier for PSM scoring.
"""
struct ProbitScorer <: PSMScoringModel
    max_iter::Int
    step_size::Float64
end

ProbitScorer() = ProbitScorer(30, 1.0)

#############################################################################
# 2. PairingStrategy - Target-Decoy Pairing
#############################################################################

"""
Abstract type for target-decoy pairing strategies.
Controls how precursors are grouped for semi-supervised learning.
"""
abstract type PairingStrategy end

"""
    RandomPairing <: PairingStrategy

Randomly pairs precursors within iRT bins regardless of target/decoy status.
Creates T-T, T-D, and D-D pairs.
"""
struct RandomPairing <: PairingStrategy
    seed::Int
    bin_size::Int
end

RandomPairing() = RandomPairing(1844, 1000)

"""
    NoPairing <: PairingStrategy

Each PSM gets a unique pair_id (no pairing).
Used for single-run analysis without MBR.
"""
struct NoPairing <: PairingStrategy end

#############################################################################
# 3. TrainingDataStrategy - Per-Iteration Selection
#############################################################################

"""
Abstract type for training data selection strategies.
Controls how positive/negative examples are selected at each iteration.
"""
abstract type TrainingDataStrategy end

"""
    QValueNegativeMining <: TrainingDataStrategy

Standard iterative training strategy:
- Iteration 1: Use all data
- Later iterations: Filter by q-value, convert worst targets to negatives via PEP
"""
struct QValueNegativeMining <: TrainingDataStrategy
    max_q_value::Float32
    min_pep_threshold::Float32
end

QValueNegativeMining() = QValueNegativeMining(0.01f0, 0.90f0)

"""
    AllDataSelection <: TrainingDataStrategy

Use all training data without filtering (used for first iteration).
"""
struct AllDataSelection <: TrainingDataStrategy end

#############################################################################
# 4. FeatureSelectionStrategy - Feature Selection
#############################################################################

"""
Abstract type for feature selection strategies.
Controls which features are used at each iteration.
"""
abstract type FeatureSelectionStrategy end

"""
    IterativeFeatureSelection <: FeatureSelectionStrategy

Uses base features early, adds MBR features at specified iteration.
"""
struct IterativeFeatureSelection <: FeatureSelectionStrategy
    base_features::Vector{Symbol}
    mbr_features::Vector{Symbol}
    mbr_start_iteration::Int
end

"""
    StaticFeatureSelection <: FeatureSelectionStrategy

Uses fixed feature set for all iterations.
"""
struct StaticFeatureSelection <: FeatureSelectionStrategy
    features::Vector{Symbol}
end

#############################################################################
# 5. IterationScheme - Training Iterations
#############################################################################

"""
Abstract type for iteration schemes.
Controls the number and configuration of training iterations.
"""
abstract type IterationScheme end

"""
    FixedIterationScheme <: IterationScheme

Fixed number of boosting rounds per iteration.
Default: [100, 200, 200] = 3 iterations with increasing rounds.
"""
struct FixedIterationScheme <: IterationScheme
    rounds::Vector{Int}
end

FixedIterationScheme() = FixedIterationScheme([100, 200, 200])

"""
    SinglePassScheme <: IterationScheme

Single training pass with fixed rounds (for probit regression).
"""
struct SinglePassScheme <: IterationScheme
    num_rounds::Int
end

SinglePassScheme() = SinglePassScheme(30)

#############################################################################
# 6. MBRUpdateStrategy - Match-Between-Runs
#############################################################################

"""
Abstract type for match-between-runs feature computation.
Controls how MBR features are computed and updated.
"""
abstract type MBRUpdateStrategy end

"""
    PairBasedMBR <: MBRUpdateStrategy

Computes MBR features based on paired precursors across runs.
"""
struct PairBasedMBR <: MBRUpdateStrategy
    q_cutoff::Float32
end

PairBasedMBR() = PairBasedMBR(0.01f0)

"""
    NoMBR <: MBRUpdateStrategy

MBR disabled - no MBR features computed.
"""
struct NoMBR <: MBRUpdateStrategy end

# Query functions for MBR support
has_mbr_support(::PairBasedMBR) = true
has_mbr_support(::NoMBR) = false

#############################################################################
# 7. ScoringPhase - Training vs Prediction Phase
#############################################################################

"""
Abstract type for scoring phases.
Controls whether training or prediction behavior is used in unified loop.
"""
abstract type ScoringPhase end

"""
    TrainingPhase <: ScoringPhase

Phase 1: Train models and predict on training data.
- Trains models and stores them
- Computes q-values in main loop (needed for select_training_data)
- Uses update_mbr_features_train_only!
"""
struct TrainingPhase <: ScoringPhase end

"""
    PredictionPhase <: ScoringPhase

Phase 2: Apply stored models to held-out test data.
- Retrieves pre-trained models (no training)
- Q-values computed inside update_mbr_features_test_only!
- Stores baseline at mbr_start_iter - 1
"""
struct PredictionPhase <: ScoringPhase end

#############################################################################
# Wrapper struct for trained probit models
#############################################################################

"""
Wrapper struct for trained probit models.
"""
struct ProbitModel
    beta::Vector{Float64}
    features::Vector{Symbol}
end
