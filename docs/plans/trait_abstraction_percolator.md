# Trait-Based Abstraction for Percolator Functions

## Overview

Abstract `sort_of_percolator!` and `train_lightgbm_model` (renamed from `*_in_memory`) using Julia traits (abstract types + multiple dispatch) to enable extensible PSM scoring with different ML models, pairing strategies, feature selection, and MBR support.

**Key Abstractions:**
1. **PSMContainer** - Abstract container for PSM data (DataFrame-backed, Arrow streaming, etc.)
2. **PSMScoringModel** - ML algorithm selection (LightGBM, Probit, etc.)
3. **PairingStrategy** - Target-decoy pairing method
4. **TrainingDataStrategy** - Per-iteration training data selection
5. **FeatureSelectionStrategy** - Feature selection per iteration
6. **IterationScheme** - Training iteration configuration
7. **MBRUpdateStrategy** - Match-between-runs feature computation

## Key Files

- `src/utils/ML/percolatorSortOf.jl` - Main algorithm (lines 242-431)
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl` - Entry point (lines 331-384)
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_config.jl` - Existing config
- `src/utils/ML/lightgbm_utils.jl` - LightGBM wrappers

## Function Renaming

| Old Name | New Name |
|----------|----------|
| `sort_of_percolator_in_memory!` | `sort_of_percolator!` |
| `train_lightgbm_model_in_memory` | `train_lightgbm_model` |
| `train_probit_model_in_memory` | `train_probit_model` |

---

## Part 0: PSMContainer Abstraction

### File: `src/utils/ML/psm_container.jl`

```julia
#=
PSMContainer abstraction for PSM data storage.
Provides a unified interface for different backing implementations.
=#

export AbstractPSMContainer, DataFramePSMContainer
export get_column, set_column!, has_column, nrows, get_view, copy_container
export get_cv_folds, get_fold_indices, get_train_indices

#############################################################################
# Abstract Interface
#############################################################################

"""
    AbstractPSMContainer

Abstract type for PSM data containers.
Implementations must provide access to PSM data with standard column operations.

# Required Interface
- `nrows(container)` - Number of PSMs
- `get_column(container, col::Symbol)` - Get column data
- `set_column!(container, col::Symbol, data)` - Set column data
- `has_column(container, col::Symbol)` - Check if column exists
- `get_view(container, indices)` - Get view of subset
- `copy_container(container)` - Create a copy
- `get_cv_folds(container)` - Get unique CV fold values
- `get_fold_indices(container, fold)` - Get indices for a CV fold
- `get_train_indices(container, test_fold)` - Get training indices (all except test_fold)
- `sort_container!(container, cols)` - Sort in place by columns
- `to_dataframe(container)` - Convert to DataFrame for ML training
"""
abstract type AbstractPSMContainer end

# Required columns for PSM scoring
const REQUIRED_PSM_COLUMNS = [
    :target,           # Bool - target/decoy label
    :cv_fold,          # UInt8 - cross-validation fold
    :precursor_idx,    # UInt32 - precursor identifier
    :ms_file_idx,      # UInt32 - MS file index
    :irt_pred,         # Float32 - predicted iRT
    :irt_obs,          # Float32 - observed iRT
    :isotopes_captured # Tuple{Int8,Int8} - isotope info
]

# Optional columns for MBR
const MBR_PSM_COLUMNS = [
    :weights,                  # Vector{Float32} - fragment weights
    :irts,                     # Vector{Float32} - fragment iRTs
    :weight,                   # Float32 - precursor weight
    :log2_intensity_explained  # Float32 - MS1 intensity metric
]

#############################################################################
# DataFrame-backed Implementation
#############################################################################

"""
    DataFramePSMContainer <: AbstractPSMContainer

PSM container backed by a DataFrame.
This is the standard in-memory implementation.
"""
struct DataFramePSMContainer <: AbstractPSMContainer
    data::DataFrame

    function DataFramePSMContainer(df::DataFrame)
        # Validate required columns
        for col in REQUIRED_PSM_COLUMNS
            if !hasproperty(df, col)
                error("Missing required column: $col")
            end
        end
        new(df)
    end
end

# Allow construction from DataFrame without validation (for internal use)
DataFramePSMContainer(df::DataFrame, ::Val{:unsafe}) = new(df)

# Core interface implementation
nrows(container::DataFramePSMContainer) = nrow(container.data)

function get_column(container::DataFramePSMContainer, col::Symbol)
    return container.data[!, col]
end

function set_column!(container::DataFramePSMContainer, col::Symbol, data)
    container.data[!, col] = data
    return nothing
end

function has_column(container::DataFramePSMContainer, col::Symbol)
    return hasproperty(container.data, col)
end

function get_view(container::DataFramePSMContainer, indices::AbstractVector{<:Integer})
    return DataFramePSMContainerView(container, indices)
end

function copy_container(container::DataFramePSMContainer)
    return DataFramePSMContainer(copy(container.data), Val(:unsafe))
end

function get_cv_folds(container::DataFramePSMContainer)
    return unique(container.data[!, :cv_fold])
end

function get_fold_indices(container::DataFramePSMContainer, fold)
    cv_col = container.data[!, :cv_fold]
    return findall(==(fold), cv_col)
end

function get_train_indices(container::DataFramePSMContainer, test_fold)
    cv_col = container.data[!, :cv_fold]
    return findall(!=(test_fold), cv_col)
end

function sort_container!(container::DataFramePSMContainer, cols::Vector{Symbol})
    sort!(container.data, cols)
    return nothing
end

function to_dataframe(container::DataFramePSMContainer)
    return container.data
end

# Access underlying DataFrame (for backward compatibility)
Base.getproperty(c::DataFramePSMContainer, s::Symbol) = s === :data ? getfield(c, :data) : get_column(c, s)

#############################################################################
# View Implementation (for train/test splits)
#############################################################################

"""
    DataFramePSMContainerView <: AbstractPSMContainer

A view into a DataFramePSMContainer for a subset of rows.
Used for train/test splits without copying data.
"""
struct DataFramePSMContainerView <: AbstractPSMContainer
    parent::DataFramePSMContainer
    indices::Vector{Int}
end

nrows(view::DataFramePSMContainerView) = length(view.indices)

function get_column(view::DataFramePSMContainerView, col::Symbol)
    return @view get_column(view.parent, col)[view.indices]
end

function set_column!(view::DataFramePSMContainerView, col::Symbol, data)
    parent_col = get_column(view.parent, col)
    parent_col[view.indices] = data
    return nothing
end

function has_column(view::DataFramePSMContainerView, col::Symbol)
    return has_column(view.parent, col)
end

function get_view(view::DataFramePSMContainerView, indices::AbstractVector{<:Integer})
    # Compose views
    return DataFramePSMContainerView(view.parent, view.indices[indices])
end

function copy_container(view::DataFramePSMContainerView)
    # Copy creates a new container with just the viewed rows
    return DataFramePSMContainer(copy(view.parent.data[view.indices, :]), Val(:unsafe))
end

function to_dataframe(view::DataFramePSMContainerView)
    return @view view.parent.data[view.indices, :]
end

# Convenience: iterate over rows
Base.eachindex(c::AbstractPSMContainer) = 1:nrows(c)

#############################################################################
# Utility Functions
#############################################################################

"""
    select_subset(container::AbstractPSMContainer, predicate::Function) -> AbstractPSMContainer

Select rows matching the predicate function.
Returns a copy with only matching rows.
"""
function select_subset(container::AbstractPSMContainer, predicate::Function)
    df = to_dataframe(container)
    mask = [predicate(df, i) for i in 1:nrows(container)]
    indices = findall(mask)
    return copy_container(get_view(container, indices))
end

"""
    get_feature_matrix(container::AbstractPSMContainer, features::Vector{Symbol}) -> Matrix{Float32}

Extract features as a matrix for ML training.
"""
function get_feature_matrix(container::AbstractPSMContainer, features::Vector{Symbol})
    n = nrows(container)
    m = length(features)
    X = Matrix{Float32}(undef, n, m)
    for (j, feat) in enumerate(features)
        col = get_column(container, feat)
        for i in 1:n
            X[i, j] = Float32(col[i])
        end
    end
    return X
end

"""
    get_labels(container::AbstractPSMContainer) -> Vector{Bool}

Get target labels for ML training.
"""
function get_labels(container::AbstractPSMContainer)
    return collect(get_column(container, :target))
end
```

---

## Part 1: Trait Hierarchy

### File: `src/utils/ML/scoring_traits.jl`

```julia
#=
Trait hierarchy for PSM scoring algorithms.
Each trait represents a single point of abstraction in the scoring pipeline.
=#

export PSMScoringModel, LightGBMScorer, ProbitScorer
export PairingStrategy, RandomPairing, NoPairing
export TrainingDataStrategy, QValueNegativeMining, AllDataSelection
export FeatureSelectionStrategy, IterativeFeatureSelection, StaticFeatureSelection
export IterationScheme, FixedIterationScheme, SinglePassScheme
export MBRUpdateStrategy, PairBasedMBR, NoMBR

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
struct NoPairing end

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
```

---

## Part 2: ScoringConfig Struct

### File: `src/utils/ML/scoring_config.jl`

```julia
#=
Unified configuration for PSM scoring algorithm.
Combines all trait types into a single configuration.
=#

export ScoringConfig, default_scoring_config

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
```

---

## Part 3: Trait Dispatch Functions

### File: `src/utils/ML/pairing.jl`

```julia
#=
Pairing strategy implementations.
Extracted from percolatorSortOf.jl lines 94-240.
Uses AbstractPSMContainer for data access.
=#

export assign_pairs!, getIrtBins

const PAIRING_RANDOM_SEED = 1844
const IRT_BIN_SIZE = 1000

"""
    assign_pairs!(psms::AbstractPSMContainer, strategy::PairingStrategy) -> Nothing

Assign pair IDs to PSMs based on the pairing strategy.
Modifies `psms` in-place by adding/updating :pair_id and :irt_bin_idx columns.
"""
function assign_pairs!(psms::AbstractPSMContainer, strategy::RandomPairing)
    n = nrows(psms)
    set_column!(psms, :pair_id, zeros(UInt32, n))

    # Compute iRT bins
    irt_pred = collect(get_column(psms, :irt_pred))
    irt_bin_idx = getIrtBins(irt_pred, strategy.bin_size)
    set_column!(psms, :irt_bin_idx, irt_bin_idx)

    # Get grouping columns
    cv_fold = collect(get_column(psms, :cv_fold))
    isotopes = collect(get_column(psms, :isotopes_captured))
    precursor_idx = collect(get_column(psms, :precursor_idx))
    pair_ids = collect(get_column(psms, :pair_id))

    # Group by (irt_bin_idx, cv_fold, isotopes_captured)
    # Use a dictionary to track groups
    groups = Dict{Tuple{UInt32, UInt8, Tuple{Int8,Int8}}, Vector{Int}}()
    for i in 1:n
        key = (irt_bin_idx[i], cv_fold[i], isotopes[i])
        if !haskey(groups, key)
            groups[key] = Int[]
        end
        push!(groups[key], i)
    end

    # Assign pair IDs within each group
    last_pair_id = zero(UInt32)
    for indices in values(groups)
        last_pair_id = assign_random_pair_ids_for_indices!(
            pair_ids, precursor_idx, indices, last_pair_id, strategy.seed
        )
    end

    set_column!(psms, :pair_id, pair_ids)
    return nothing
end

function assign_pairs!(psms::AbstractPSMContainer, ::NoPairing)
    n = nrows(psms)
    set_column!(psms, :pair_id, UInt32.(1:n))
    set_column!(psms, :irt_bin_idx, zeros(UInt32, n))
    return nothing
end

# Helper functions
function getIrtBins(irts::AbstractVector{R}, bin_size::Int = IRT_BIN_SIZE) where {R<:Real}
    sort_idx = sortperm(irts)
    bin_idx, bin_count = zero(UInt32), zero(UInt32)
    bin_idxs = similar(irts, UInt32, length(irts))
    for idx in sort_idx
        bin_count += one(UInt32)
        bin_idxs[idx] = bin_idx
        if bin_count >= bin_size
            bin_idx += one(UInt32)
            bin_count = zero(UInt32)
        end
    end
    return bin_idxs
end

function assign_random_pair_ids_for_indices!(
    pair_ids::Vector{UInt32},
    precursor_idx::Vector{UInt32},
    indices::Vector{Int},
    last_pair_id::UInt32,
    seed::Int
)
    # Get unique precursors in this group
    group_precursors = precursor_idx[indices]
    unique_precursors = unique(group_precursors)
    n_precursors = length(unique_precursors)

    # Shuffle precursors
    perm = randperm(MersenneTwister(seed), n_precursors)
    shuffled_precursors = unique_precursors[perm]

    # Create mapping
    precursor_to_pair = Dict{UInt32, UInt32}()
    n_pairs = n_precursors ÷ 2

    for i in 1:n_pairs
        last_pair_id += one(UInt32)
        precursor_to_pair[shuffled_precursors[2*i - 1]] = last_pair_id
        precursor_to_pair[shuffled_precursors[2*i]] = last_pair_id
    end

    # Handle odd number
    if isodd(n_precursors)
        last_pair_id += one(UInt32)
        precursor_to_pair[shuffled_precursors[end]] = last_pair_id
    end

    # Assign pair IDs for this group
    for idx in indices
        pair_ids[idx] = precursor_to_pair[precursor_idx[idx]]
    end

    return last_pair_id
end
```

### File: `src/utils/ML/model_training.jl`

```julia
#=
Model training dispatch functions.
Extracted from percolatorSortOf.jl lines 807-830.
Uses AbstractPSMContainer for data access.
=#

export train_model, predict_scores

"""
    train_model(model::PSMScoringModel, psms::AbstractPSMContainer, features, num_rounds) -> TrainedModel

Train a PSM scoring model on the given data.
Returns a trained model object that can be used for prediction.
"""
function train_model(model::LightGBMScorer, psms::AbstractPSMContainer,
                     features::Vector{Symbol}, num_rounds::Int)
    hp = model.hyperparams

    classifier = build_lightgbm_classifier(
        num_iterations = num_rounds,
        max_depth = get(hp, :max_depth, 10),
        learning_rate = get(hp, :learning_rate, 0.05),
        num_leaves = get(hp, :num_leaves, 63),
        feature_fraction = get(hp, :feature_fraction, 0.5),
        bagging_fraction = get(hp, :bagging_fraction, 0.25),
        bagging_freq = get(hp, :bagging_fraction, 0.25) < 1 ? 1 : 0,
        min_data_in_leaf = get(hp, :min_data_in_leaf, 500),
        min_gain_to_split = get(hp, :min_gain_to_split, 0.5),
        is_unbalance = true
    )

    # Use container interface for feature extraction
    feature_matrix = get_feature_matrix(psms, features)
    labels = get_labels(psms)
    return fit_lightgbm_model(classifier, feature_matrix, labels; positive_label=true)
end

function train_model(model::ProbitScorer, psms::AbstractPSMContainer,
                     features::Vector{Symbol}, num_rounds::Int)
    y = get_labels(psms)
    X = get_feature_matrix(psms, features)

    chunk_size = max(1, nrows(psms) ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:nrows(psms), chunk_size)

    beta = zeros(Float64, length(features))
    beta = ProbitRegression(beta, Matrix{Float64}(X), y, data_chunks;
                           max_iter=model.max_iter, step_size=model.step_size)

    return ProbitModel(beta, features)
end

"""
    predict_scores(model, psms::AbstractPSMContainer) -> Vector{Float32}

Get probability predictions from a trained model.
"""
function predict_scores(model::LightGBMModel, psms::AbstractPSMContainer)
    feature_matrix = get_feature_matrix(psms, model.features)
    return lightgbm_predict(model, feature_matrix; output_type=Float32)
end

function predict_scores(model::ProbitModel, psms::AbstractPSMContainer)
    scores = zeros(Float32, nrows(psms))
    chunk_size = max(1, nrows(psms) ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:nrows(psms), chunk_size)
    X = get_feature_matrix(psms, model.features)
    ModelPredictProbs!(scores, Matrix{Float64}(X), model.beta, data_chunks)
    return scores
end

"""
Wrapper struct for trained probit models.
"""
struct ProbitModel
    beta::Vector{Float64}
    features::Vector{Symbol}
end
```

### File: `src/utils/ML/training_selection.jl`

```julia
#=
Training data selection strategies.
Extracted from percolatorSortOf.jl lines 972-1009.
Uses AbstractPSMContainer for data access.
=#

export select_training_data

"""
    select_training_data(psms::AbstractPSMContainer, strategy::TrainingDataStrategy, iteration::Int) -> AbstractPSMContainer

Select training data based on the strategy and current iteration.
Returns a copy of the filtered training data.
"""
function select_training_data(psms::AbstractPSMContainer, ::AllDataSelection, ::Int)
    return copy_container(psms)
end

function select_training_data(psms::AbstractPSMContainer, strategy::QValueNegativeMining, iteration::Int)
    if iteration == 1
        return copy_container(psms)
    end

    # Deep copy to avoid modifying original labels
    psms_train = copy_container(psms)

    # Convert worst-scoring targets to negatives using PEP
    trace_probs = collect(get_column(psms_train, :trace_prob))
    targets = collect(get_column(psms_train, :target))

    order = sortperm(trace_probs, rev=true)
    sorted_scores = trace_probs[order]
    sorted_targets = targets[order]

    PEPs = Vector{Float32}(undef, length(order))
    get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

    idx_cutoff = findfirst(x -> x >= strategy.min_pep_threshold, PEPs)
    if !isnothing(idx_cutoff)
        worst_idxs = order[idx_cutoff:end]
        targets[worst_idxs] .= false
        set_column!(psms_train, :target, targets)
    end

    # Filter: keep all decoys + targets passing q-value threshold
    q_values = get_column(psms_train, :q_value)
    targets = get_column(psms_train, :target)
    keep_mask = [(!t) || (t && q <= strategy.max_q_value) for (t, q) in zip(targets, q_values)]
    keep_indices = findall(keep_mask)

    return copy_container(get_view(psms_train, keep_indices))
end
```

### File: `src/utils/ML/feature_selection.jl`

```julia
#=
Feature selection strategies.
Uses AbstractPSMContainer for data access.
=#

export get_features, filter_available_features

"""
    get_features(strategy::FeatureSelectionStrategy, iteration::Int, total_iterations::Int) -> Vector{Symbol}

Get the features to use for the current iteration.
"""
function get_features(strategy::IterativeFeatureSelection, iteration::Int, ::Int)
    if iteration < strategy.mbr_start_iteration
        return strategy.base_features
    else
        return vcat(strategy.base_features, strategy.mbr_features)
    end
end

function get_features(strategy::StaticFeatureSelection, ::Int, ::Int)
    return strategy.features
end

"""
    filter_available_features(features::Vector{Symbol}, psms::AbstractPSMContainer) -> Vector{Symbol}

Filter features to only those available in the container.
"""
function filter_available_features(features::Vector{Symbol}, psms::AbstractPSMContainer)
    return [f for f in features if has_column(psms, f)]
end
```

### File: `src/utils/ML/iteration_scheme.jl`

```julia
#=
Iteration scheme implementations.
=#

export get_iteration_rounds, get_iterations_count

"""
    get_iteration_rounds(scheme::IterationScheme) -> Vector{Int}

Get the number of boosting rounds for each iteration.
"""
function get_iteration_rounds(scheme::FixedIterationScheme)
    return scheme.rounds
end

function get_iteration_rounds(scheme::SinglePassScheme)
    return [scheme.num_rounds]
end

"""
    get_iterations_count(scheme::IterationScheme) -> Int

Get total number of training iterations.
"""
get_iterations_count(scheme::FixedIterationScheme) = length(scheme.rounds)
get_iterations_count(scheme::SinglePassScheme) = 1
```

### File: `src/utils/ML/mbr_update.jl`

```julia
#=
MBR feature update strategies.
Extracted from percolatorSortOf.jl lines 839-970.
Uses AbstractPSMContainer for data access.
=#

export initialize_mbr_columns!, update_mbr_features!

"""
    initialize_mbr_columns!(psms::AbstractPSMContainer, strategy::MBRUpdateStrategy) -> Nothing

Initialize MBR-related columns in the PSM container.
"""
function initialize_mbr_columns!(psms::AbstractPSMContainer, ::PairBasedMBR)
    n = nrows(psms)
    set_column!(psms, :trace_prob, zeros(Float32, n))
    set_column!(psms, :q_value, zeros(Float64, n))
    set_column!(psms, :MBR_max_pair_prob, zeros(Float32, n))
    set_column!(psms, :MBR_best_irt_diff, zeros(Float32, n))
    set_column!(psms, :MBR_log2_weight_ratio, zeros(Float32, n))
    set_column!(psms, :MBR_log2_explained_ratio, zeros(Float32, n))
    set_column!(psms, :MBR_rv_coefficient, zeros(Float32, n))
    set_column!(psms, :MBR_is_best_decoy, trues(n))
    set_column!(psms, :MBR_num_runs, zeros(Int32, n))
    set_column!(psms, :MBR_transfer_candidate, falses(n))
    set_column!(psms, :MBR_is_missing, falses(n))
    return nothing
end

function initialize_mbr_columns!(psms::AbstractPSMContainer, ::NoMBR)
    n = nrows(psms)
    set_column!(psms, :trace_prob, zeros(Float32, n))
    set_column!(psms, :q_value, zeros(Float64, n))
    return nothing
end

"""
    update_mbr_features!(psms_train, psms_test, iteration, mbr_start_iter, strategy) -> Nothing

Update MBR features based on the strategy.
"""
function update_mbr_features!(psms_train::AbstractPSMContainer, psms_test::AbstractPSMContainer,
                               iteration::Int, mbr_start_iter::Int, strategy::PairBasedMBR)
    if iteration >= mbr_start_iter - 1
        trace_probs = get_column(psms_test, :trace_prob)
        targets = get_column(psms_test, :target)
        q_values = collect(get_column(psms_test, :q_value))
        get_qvalues!(collect(trace_probs), collect(targets), q_values)
        set_column!(psms_test, :q_value, q_values)

        # summarize_precursors! needs DataFrame - convert temporarily
        summarize_precursors!(to_dataframe(psms_test), q_cutoff=strategy.q_cutoff)
        summarize_precursors!(to_dataframe(psms_train), q_cutoff=strategy.q_cutoff)
    end
    return nothing
end

function update_mbr_features!(::AbstractPSMContainer, ::AbstractPSMContainer,
                               ::Int, ::Int, ::NoMBR)
    return nothing
end

# summarize_precursors! stays in percolatorSortOf.jl - operates on DataFrame
# due to complex groupby operations that require DataFrame
```

---

## Part 4: Main Generic Function

### File: `src/utils/ML/percolator_generic.jl`

```julia
#=
Generic percolator scoring function using trait-based abstraction.
Uses AbstractPSMContainer for data access.
=#

export percolator_scoring!

"""
    percolator_scoring!(psms::AbstractPSMContainer, config::ScoringConfig;
                        show_progress::Bool=true, verbose::Bool=false) -> Dict{UInt8, Any}

Generic PSM scoring using configurable traits.

# Arguments
- `psms::AbstractPSMContainer`: PSMs to score (modified in-place)
- `config::ScoringConfig`: Configuration specifying all algorithm components
- `show_progress::Bool`: Whether to show progress bar
- `verbose::Bool`: Whether to print verbose logging

# Returns
- `Dict{UInt8, Vector{TrainedModel}}`: Dictionary mapping CV fold to trained models
"""
function percolator_scoring!(psms::AbstractPSMContainer, config::ScoringConfig;
                             show_progress::Bool = true, verbose::Bool = false)

    # Step 1: Apply pairing strategy
    assign_pairs!(psms, config.pairing)
    sort_container!(psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])

    # Log dataset info
    if verbose
        targets = get_column(psms, :target)
        n_targets = sum(targets)
        n_decoys = sum(.!targets)
        @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys (total: $(nrows(psms)) PSMs)"
    end

    # Step 2: Initialize columns
    initialize_mbr_columns!(psms, config.mbr_update)

    # Step 3: Get iteration scheme
    iteration_rounds = get_iteration_rounds(config.iteration_scheme)
    total_iterations = length(iteration_rounds)
    mbr_start_iter = total_iterations  # MBR features used in last iteration

    # Determine actual iterations to run
    iterations_to_run = has_mbr_support(config.mbr_update) ? total_iterations : max(mbr_start_iter - 1, 1)

    # Step 4: Setup cross-validation
    unique_cv_folds = get_cv_folds(psms)
    models = Dict{UInt8, Vector{Any}}()

    # Pre-compute fold indices
    fold_indices = Dict(fold => get_fold_indices(psms, fold) for fold in unique_cv_folds)
    train_indices = Dict(fold => get_train_indices(psms, fold) for fold in unique_cv_folds)

    # Step 5: Allocate output arrays
    n = nrows(psms)
    prob_test = zeros(Float32, n)
    prob_train = zeros(Float32, n)
    nonMBR_estimates = zeros(Float32, n)
    MBR_estimates = zeros(Float32, n)

    # Progress tracking
    total_progress_steps = length(unique_cv_folds) * iterations_to_run
    pbar = show_progress ? ProgressBar(total=total_progress_steps) : nothing

    Random.seed!(1776)

    # Step 6: Cross-validation training loop
    for test_fold_idx in unique_cv_folds
        train_idx = train_indices[test_fold_idx]
        test_idx = fold_indices[test_fold_idx]

        # Create views for train/test splits
        psms_train = get_view(psms, train_idx)
        psms_test = get_view(psms, test_idx)

        fold_models = Vector{Any}(undef, total_iterations)

        for (itr, num_round) in enumerate(iteration_rounds)
            # Get training data based on strategy
            training_strategy = itr == 1 ? AllDataSelection() : config.training_data
            psms_train_itr = select_training_data(psms_train, training_strategy, itr)

            # Get features for this iteration
            features = get_features(config.feature_selection, itr, total_iterations)
            features = filter_available_features(features, psms_train_itr)

            # Train model
            model = train_model(config.model, psms_train_itr, features, num_round)
            fold_models[itr] = model

            # Predict on training set
            prob_train[train_idx] = predict_scores(model, psms_train)
            set_column!(psms_train, :trace_prob, prob_train[train_idx])

            # Compute q-values on training set
            trace_probs_train = collect(get_column(psms_train, :trace_prob))
            targets_train = collect(get_column(psms_train, :target))
            q_values_train = zeros(Float64, length(train_idx))
            get_qvalues!(trace_probs_train, targets_train, q_values_train)
            set_column!(psms_train, :q_value, q_values_train)

            # Predict on test set
            prob_test[test_idx] = predict_scores(model, psms_test)
            set_column!(psms_test, :trace_prob, prob_test[test_idx])

            # Store non-MBR baseline before MBR iteration
            if itr == (mbr_start_iter - 1)
                nonMBR_estimates[test_idx] = prob_test[test_idx]
            end

            # Update MBR features
            update_mbr_features!(psms_train, psms_test, itr, mbr_start_iter, config.mbr_update)

            show_progress && update(pbar)

            # Early exit for non-MBR mode
            if !has_mbr_support(config.mbr_update) && itr == (mbr_start_iter - 1)
                break
            end
        end

        # Store final predictions
        if has_mbr_support(config.mbr_update)
            MBR_estimates[test_idx] = collect(get_column(psms_test, :trace_prob))
        else
            prob_test[test_idx] = collect(get_column(psms_test, :trace_prob))
        end

        models[test_fold_idx] = fold_models
    end

    # Step 7: Finalize scoring
    finalize_scoring!(psms, config.mbr_update, prob_test, nonMBR_estimates, MBR_estimates)

    return models
end

"""
    finalize_scoring!(psms, mbr_strategy, prob_test, nonMBR_estimates, MBR_estimates)

Finalize PSM scoring by setting appropriate probability columns and MBR transfer candidates.
"""
function finalize_scoring!(psms::AbstractPSMContainer, strategy::PairBasedMBR,
                           ::Vector{Float32}, nonMBR_estimates::Vector{Float32},
                           MBR_estimates::Vector{Float32})
    # Compute q-values on non-MBR baseline
    targets = collect(get_column(psms, :target))
    qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
    get_qvalues!(nonMBR_estimates, targets, qvals_prev)
    pass_mask = (qvals_prev .<= strategy.q_cutoff)

    # Find passing probability threshold
    has_passing = !isempty(pass_mask) && any(pass_mask)
    prob_thresh = has_passing ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)

    # Mark transfer candidates
    mbr_max_pair_prob = collect(get_column(psms, :MBR_max_pair_prob))
    transfer_candidates = .!pass_mask .& (mbr_max_pair_prob .>= prob_thresh)
    set_column!(psms, :MBR_transfer_candidate, transfer_candidates)

    # Store probabilities
    set_column!(psms, :trace_prob, nonMBR_estimates)
    set_column!(psms, :MBR_boosted_trace_prob, MBR_estimates)
end

function finalize_scoring!(psms::AbstractPSMContainer, ::NoMBR,
                           prob_test::Vector{Float32}, ::Vector{Float32},
                           ::Vector{Float32})
    set_column!(psms, :trace_prob, prob_test)
end

#############################################################################
# DataFrame convenience wrapper
#############################################################################

"""
    percolator_scoring!(psms::DataFrame, config::ScoringConfig; kwargs...) -> Dict

Convenience wrapper that wraps DataFrame in DataFramePSMContainer.
"""
function percolator_scoring!(psms::DataFrame, config::ScoringConfig;
                             show_progress::Bool = true, verbose::Bool = false)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    return percolator_scoring!(container, config; show_progress, verbose)
end
```

---

## Part 5: Updated Entry Points

### File: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl` (modified)

```julia
"""
    train_lightgbm_model(best_psms, file_paths, precursors, model_config,
                         match_between_runs, max_q_value_lightgbm_rescore,
                         max_q_value_mbr_itr, min_PEP_neg_threshold_itr;
                         show_progress=true) -> Models

Trains LightGBM model using configuration from ModelConfig.
Delegates to generic percolator_scoring! with appropriate config.

# Arguments
- `best_psms::DataFrame`: PSMs to score
- `file_paths::Vector{String}`: Paths to PSM files (for compatibility)
- `precursors::LibraryPrecursors`: Library precursor information
- `model_config::ModelConfig`: Model configuration
- `match_between_runs::Bool`: Whether to enable MBR
- `max_q_value_lightgbm_rescore::Float32`: Q-value threshold
- `max_q_value_mbr_itr::Float32`: MBR q-value threshold
- `min_PEP_neg_threshold_itr::Float32`: PEP threshold for negative mining
- `show_progress::Bool`: Whether to show progress bars

# Returns
- Dictionary of trained models per CV fold
"""
function train_lightgbm_model(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    show_progress::Bool = true
)
    # Add required columns
    best_psms[!, :accession_numbers] = [getAccessionNumbers(precursors)[pid]
                                         for pid in best_psms[!, :precursor_idx]]
    best_psms[!, :q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!, :decoy] = best_psms[!, :target] .== false

    # Build configuration from ModelConfig
    config = build_scoring_config(
        model_config,
        match_between_runs,
        max_q_value_lightgbm_rescore,
        min_PEP_neg_threshold_itr
    )

    # Delegate to generic function
    return percolator_scoring!(best_psms, config; show_progress=show_progress)
end

"""
    build_scoring_config(model_config, match_between_runs, max_q_value, min_pep_threshold) -> ScoringConfig

Build a ScoringConfig from a ModelConfig and parameters.
"""
function build_scoring_config(
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value::Float32,
    min_pep_threshold::Float32
)
    # Select model type
    model = if model_config.model_type == :lightgbm
        LightGBMScorer(model_config.hyperparams)
    elseif model_config.model_type == :probit
        ProbitScorer(get(model_config.hyperparams, :max_iter, 30), 1.0)
    else
        error("Unsupported model type: $(model_config.model_type)")
    end

    # Get features (filter MBR features from base)
    base_features = [f for f in model_config.features if !startswith(String(f), "MBR_")]
    mbr_features = match_between_runs ? [
        :MBR_max_pair_prob,
        :MBR_log2_weight_ratio,
        :MBR_log2_explained_ratio,
        :MBR_rv_coefficient,
        :MBR_is_missing
    ] : Symbol[]

    # Get iteration scheme
    iter_scheme = get(model_config.hyperparams, :iter_scheme, [100, 200, 200])

    return ScoringConfig(
        model,
        RandomPairing(),
        QValueNegativeMining(max_q_value, min_pep_threshold),
        IterativeFeatureSelection(base_features, mbr_features, length(iter_scheme)),
        FixedIterationScheme(iter_scheme),
        match_between_runs ? PairBasedMBR(max_q_value) : NoMBR()
    )
end

# Backward compatibility alias
const train_lightgbm_model_in_memory = train_lightgbm_model
```

### File: `src/utils/ML/percolatorSortOf.jl` (modified)

```julia
# Add backward compatibility wrapper at the end of file

"""
    sort_of_percolator!(psms, features, match_between_runs; kwargs...) -> Models

Generic PSM scoring function. This is the new name for sort_of_percolator_in_memory!.
See `percolator_scoring!` for the trait-based implementation.
"""
function sort_of_percolator!(psms::DataFrame,
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_lightgbm_rescore::Float32 = 0.01f0,
                  max_q_value_mbr_itr::Float32 = 0.20f0,
                  min_PEP_neg_threshold_itr = 0.90f0,
                  feature_fraction::Float64 = 0.5,
                  learning_rate::Float64 = 0.15,
                  min_data_in_leaf::Int = 1,
                  bagging_fraction::Float64 = 0.5,
                  min_gain_to_split::Float64 = 0.0,
                  max_depth::Int = 10,
                  num_leaves::Int = 63,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = false,
                  show_progress::Bool = true,
                  verbose_logging::Bool = false)

    # Build config from kwargs
    base_features = [f for f in features if !startswith(String(f), "MBR_")]
    mbr_features = [f for f in features if startswith(String(f), "MBR_")]

    config = ScoringConfig(
        LightGBMScorer(Dict{Symbol,Any}(
            :feature_fraction => feature_fraction,
            :learning_rate => learning_rate,
            :min_data_in_leaf => min_data_in_leaf,
            :bagging_fraction => bagging_fraction,
            :min_gain_to_split => min_gain_to_split,
            :max_depth => max_depth,
            :num_leaves => num_leaves
        )),
        RandomPairing(),
        QValueNegativeMining(max_q_value_lightgbm_rescore, Float32(min_PEP_neg_threshold_itr)),
        IterativeFeatureSelection(base_features, mbr_features, length(iter_scheme)),
        FixedIterationScheme(iter_scheme),
        match_between_runs ? PairBasedMBR(max_q_value_lightgbm_rescore) : NoMBR()
    )

    return percolator_scoring!(psms, config; show_progress=show_progress, verbose=verbose_logging)
end

# Backward compatibility alias
const sort_of_percolator_in_memory! = sort_of_percolator!
```

---

## Part 6: File Organization

### New Files to Create

```
src/utils/ML/
├── psm_container.jl         # AbstractPSMContainer abstraction (Part 0)
├── scoring_traits.jl        # Abstract type definitions (Part 1)
├── scoring_config.jl        # ScoringConfig struct + factories (Part 2)
├── pairing.jl               # PairingStrategy implementations (Part 3)
├── model_training.jl        # PSMScoringModel train/predict (Part 3)
├── training_selection.jl    # TrainingDataStrategy implementations (Part 3)
├── feature_selection.jl     # FeatureSelectionStrategy implementations (Part 3)
├── iteration_scheme.jl      # IterationScheme implementations (Part 3)
├── mbr_update.jl            # MBRUpdateStrategy implementations (Part 3)
└── percolator_generic.jl    # Main percolator_scoring! function (Part 4)
```

### Files to Modify

1. **`src/utils/ML/percolatorSortOf.jl`**
   - Add `sort_of_percolator!` as main function name
   - Keep `sort_of_percolator_in_memory!` as alias
   - Extract helper functions to new files

2. **`src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`**
   - Rename `train_lightgbm_model_in_memory` → `train_lightgbm_model`
   - Add `build_scoring_config` function
   - Keep old name as alias

3. **`src/utils/ML/ML.jl`** (or main module)
   - Add includes for new files
   - Export new types and functions

---

## Implementation Order

1. **Phase 1**: Create `psm_container.jl` with AbstractPSMContainer and DataFramePSMContainer
2. **Phase 2**: Create `scoring_traits.jl` with all abstract types
3. **Phase 3**: Create `scoring_config.jl` with ScoringConfig struct
4. **Phase 4**: Create dispatch function files:
   - `pairing.jl`
   - `model_training.jl`
   - `training_selection.jl`
   - `feature_selection.jl`
   - `iteration_scheme.jl`
   - `mbr_update.jl`
5. **Phase 5**: Create `percolator_generic.jl` with main function
6. **Phase 6**: Update `score_psms.jl` with new names and `build_scoring_config`
7. **Phase 7**: Update `percolatorSortOf.jl` with wrappers and aliases
8. **Phase 8**: Update module includes and exports
9. **Phase 9**: Add tests for new trait system

---

## Verification

1. Run existing tests: `julia --project -e 'using Pkg; Pkg.test()'`
2. Test backward compatibility: ensure old function names still work
3. Test new trait-based API with LightGBM config
4. Test new trait-based API with Probit config
5. Verify MBR features computed correctly
6. Compare FDR results between old and new implementations
7. Integration test: `SearchDIA("./data/ecoli_test/ecoli_test_params.json")`
