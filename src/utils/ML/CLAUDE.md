# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the ML utilities module.

## Overview

The ML module provides machine learning utilities for Pioneer.jl, including:
- **PSM Scoring**: Trait-based semi-supervised learning for peptide-spectrum match scoring
- **FDR Control**: Q-value and PEP estimation utilities
- **Regression Models**: Probit regression, spectral linear regression, spline fitting
- **Smoothing**: Whittaker-Henderson smoothing, B-splines

## Module Structure

```
src/utils/ML/
├── # Trait-Based PSM Scoring System (NEW)
├── scoring_traits.jl         # 6 trait types defining algorithm components
├── scoring_config.jl         # ScoringConfig combining traits
├── psm_container.jl          # AbstractPSMContainer data abstraction
├── percolator_generic.jl     # Main percolator_scoring!() entry point
├── pairing.jl                # Target-decoy pairing strategies
├── model_training.jl         # ML model training (LightGBM, Probit)
├── training_selection.jl     # Training data selection strategies
├── feature_selection.jl      # Feature selection strategies
├── iteration_scheme.jl       # Iteration configuration
├── mbr_update.jl             # Match-between-runs feature computation
│
├── # Legacy/Wrapper
├── percolatorSortOf.jl       # Legacy wrapper → delegates to percolator_generic.jl
│
├── # ML Utilities
├── lightgbm_utils.jl         # LightGBM API abstraction
├── probitRegression.jl       # Probit regression implementation
├── fdrUtilities.jl           # Q-value and PEP calculation
├── ftrUtilities.jl           # FTR utilities
│
├── # Regression & Splines
├── spectralLinearRegression.jl   # Huber-loss spectral deconvolution
├── libraryBSpline.jl             # B-spline library models
├── uniformBasisCubicSpline.jl    # Cubic spline implementation
├── piecewiseLinearFunction.jl    # Piecewise linear functions
├── wittakerHendersonSmoothing.jl # Whittaker-Henderson smoothing
│
└── CLAUDE.md                 # This file
```

## Trait-Based PSM Scoring Architecture

The PSM scoring system uses a trait-based design with 6 independent abstraction points:

### Trait Types (scoring_traits.jl)

```julia
# 1. PSMScoringModel - ML algorithm selection
abstract type PSMScoringModel end
├── LightGBMScorer       # Gradient boosting (default)
└── ProbitScorer         # Probit regression

# 2. PairingStrategy - Target-decoy pairing
abstract type PairingStrategy end
├── RandomPairing        # Random pairs within iRT bins
└── NoPairing            # No pairing (single-run mode)

# 3. TrainingDataStrategy - Per-iteration data selection
abstract type TrainingDataStrategy end
├── QValueNegativeMining # Filter by q-value, convert worst to negatives
└── AllDataSelection     # Use all data (iteration 1)

# 4. FeatureSelectionStrategy - Feature set selection
abstract type FeatureSelectionStrategy end
├── IterativeFeatureSelection  # Base features → add MBR features
└── StaticFeatureSelection     # Fixed feature set

# 5. IterationScheme - Training iterations
abstract type IterationScheme end
├── FixedIterationScheme  # Fixed rounds per iteration [100, 200, 200]
└── SinglePassScheme      # Single pass (for probit)

# 6. MBRUpdateStrategy - Match-between-runs
abstract type MBRUpdateStrategy end
├── PairBasedMBR         # MBR features from paired precursors
└── NoMBR                # MBR disabled

# 7. ScoringPhase - Training vs Prediction
abstract type ScoringPhase end
├── TrainingPhase        # Train models, compute q-values
└── PredictionPhase      # Apply stored models to test data

# 8. MemoryStrategy - Memory Processing
abstract type MemoryStrategy end
├── InMemoryProcessing   # All PSMs in memory (default)
└── OutOfMemoryProcessing # File-by-file processing for large datasets
```

### Configuration (scoring_config.jl)

```julia
# Combine all traits into single config
struct ScoringConfig{M,P,T,F,I,B}
    model::M              # PSMScoringModel
    pairing::P            # PairingStrategy
    training_data::T      # TrainingDataStrategy
    feature_selection::F  # FeatureSelectionStrategy
    iteration_scheme::I   # IterationScheme
    mbr_update::B         # MBRUpdateStrategy
end

# Create default LightGBM config
config = default_scoring_config(
    match_between_runs=true,
    features=base_features,
    max_q_value=0.01f0
)

# Create probit config
config = probit_scoring_config(features=features)
```

### PSM Container (psm_container.jl)

```julia
# Abstract interface for PSM data
abstract type AbstractPSMContainer end

# DataFrame-backed implementation
struct DataFramePSMContainer <: AbstractPSMContainer
    data::DataFrame
end

# View for train/test splits (no copy)
struct DataFramePSMContainerView <: AbstractPSMContainer
    parent::DataFramePSMContainer
    indices::Vector{Int}
end

# Required interface methods
nrows(container)
get_column(container, :col)
set_column!(container, :col, data)
has_column(container, :col)
get_view(container, indices)
get_cv_folds(container)
get_fold_indices(container, fold)
get_train_indices(container, test_fold)
to_dataframe(container)
```

### Main Entry Point (percolator_generic.jl)

```julia
# Generic scoring function (in-memory - default)
models = percolator_scoring!(psms, config; show_progress=true, verbose=false)

# Out-of-memory scoring for large datasets
memory = OutOfMemoryProcessing(
    max_training_psms = 5_000_000,  # Max PSMs to sample for training
    fold_file_paths = Dict(         # CV fold -> file paths
        0x00 => ["fold0_file1.arrow", "fold0_file2.arrow"],
        0x01 => ["fold1_file1.arrow", "fold1_file2.arrow"]
    )
)
models = percolator_scoring!(psms, config; memory=memory, show_progress=true)

# Algorithm flow:
# 1. assign_pairs!(psms, config.pairing)
# 2. initialize_mbr_columns!(psms, config.mbr_update)
# 3. For each CV fold:
#    a. For each iteration:
#       - select_training_data(psms_train, config.training_data, itr)
#       - get_features_for_iteration(config.feature_selection, itr)
#       - train_model!(config.model, train_data, features, num_rounds)
#       - predict!(model, test_data)
#       - update_mbr_features!(psms, config.mbr_update) [if applicable]
# 4. Return trained models
```

## Key Files Reference

### percolatorSortOf.jl

Legacy wrapper that constructs a `ScoringConfig` and delegates to `percolator_scoring!`:

```julia
function sort_of_percolator!(psms::AbstractPSMContainer, features, match_between_runs; kwargs...)
    # Build ScoringConfig from keyword arguments
    config = ScoringConfig(
        LightGBMScorer(hyperparams),
        RandomPairing(),
        QValueNegativeMining(max_q_value, min_pep_threshold),
        IterativeFeatureSelection(base_features, mbr_features, n_iters),
        FixedIterationScheme(iter_scheme),
        match_between_runs ? PairBasedMBR(max_q_value) : NoMBR()
    )
    return percolator_scoring!(psms, config; show_progress, verbose)
end
```

**Constants kept for backward compatibility:**
- `PAIRING_RANDOM_SEED = 1844` - Random seed for reproducible pairing
- `IRT_BIN_SIZE = 1000` - Size of iRT bins for pairing

**Helper functions:**
- `update_pair_statistics()` - Running statistics for MBR (used by OOM code)
- `irt_residual()` - Compute iRT residual for MBR features
- `summarize_precursors!()` - Compute MBR features per pair group
- `train_booster()` - LightGBM training helper
- MBR helper functions (pad_equal_length, MBR_rv_coefficient, etc.)

### pairing.jl

Implements `PairingStrategy` trait:

```julia
# Random pairing within iRT bins
function assign_pairs!(psms::AbstractPSMContainer, strategy::RandomPairing)
    # 1. Compute iRT bins
    # 2. Group by (irt_bin, cv_fold, isotopes_captured)
    # 3. Within each group, shuffle precursors and pair consecutively
    # Creates pair_id column
end

# No pairing (each PSM gets unique pair_id)
function assign_pairs!(psms::AbstractPSMContainer, ::NoPairing)
    # Assign sequential pair_ids
end

# Core iRT binning function
function getIrtBins(irts::AbstractVector{R}, bin_size::Int = 1000)
    # Returns bin indices based on sorted iRT order
end
```

### model_training.jl

Implements `PSMScoringModel` trait:

```julia
# LightGBM training
function train_model!(model::LightGBMScorer, data, features, num_rounds)
    # Uses lightgbm_utils.jl API
    return fit_lightgbm_model(classifier, features, labels)
end

# Probit training
function train_model!(model::ProbitScorer, data, features, num_rounds)
    # Uses probitRegression.jl
    return ProbitModel(beta, features)
end
```

### fdrUtilities.jl

FDR estimation utilities:

```julia
# Q-value calculation (in-place)
get_qvalues!(scores, targets, qvals; doSort=true)

# PEP calculation (in-place)
get_PEP!(scores, targets, PEPs; doSort=true)

# Combined wrapper
computeFDR!(scores, targets, qvals, PEPs)
```

## Common Development Patterns

### Adding a New Pairing Strategy

```julia
# 1. Define new strategy type
struct MyPairingStrategy <: PairingStrategy
    my_param::Float64
end

# 2. Implement assign_pairs! in pairing.jl
function assign_pairs!(psms::AbstractPSMContainer, strategy::MyPairingStrategy)
    # Your pairing logic
    set_column!(psms, :pair_id, pair_ids)
    set_column!(psms, :irt_bin_idx, bin_ids)
end
```

### Adding a New ML Model

```julia
# 1. Define new model type
struct MyMLModel <: PSMScoringModel
    hyperparams::Dict{Symbol, Any}
end

# 2. Implement train_model! in model_training.jl
function train_model!(model::MyMLModel, data, features, num_rounds)
    # Your training logic
    return trained_model
end

# 3. Implement predict!
function predict!(model::MyTrainedModel, data)
    # Your prediction logic
    return probabilities
end
```

### Debugging PSM Scoring

```julia
# Enable verbose logging
models = percolator_scoring!(psms, config; verbose=true)

# Check intermediate state
println("Unique pair_ids: ", length(unique(get_column(psms, :pair_id))))
println("CV folds: ", get_cv_folds(psms))

# Inspect model performance
for (fold, fold_models) in models
    println("Fold $fold: $(length(fold_models)) models")
end
```

## Memory Management

### In-Memory Processing (InMemoryProcessing)
- `DataFramePSMContainer` holds all PSMs in memory
- `DataFramePSMContainerView` provides zero-copy views for train/test splits
- Feature matrices are allocated per iteration
- Default behavior for datasets that fit in memory

### Out-of-Memory Processing (OutOfMemoryProcessing)
Implemented in `percolator_generic.jl` for datasets too large to fit in memory.

**Training Phase:**
- Samples pair_ids proportionally to meet `max_training_psms` budget
- All PSMs with same pair_id are either all sampled or all excluded (maintains pair integrity)
- Only sampled PSMs are loaded into memory for training
- Q-values computed on sampled training data only

**Prediction Phase (file-by-file):**
- Pass 1: Apply trained models to each file, write predictions back
- Pass 2: Compute global q-values across all files in CV fold
- MBR: Collect pair statistics across files, then apply features

**Key Functions:**
- `scan_pair_ids()` - Scan files to get pair_id counts
- `sample_pair_ids()` - Sample pair_ids proportionally
- `load_sampled_psms()` - Load only sampled PSMs from files
- `compute_global_qvalues!()` - Two-pass global q-value computation
- `collect_pair_statistics_oom()` - Stream MBR statistics across files
- `apply_mbr_features_from_stats!()` - Apply MBR features using statistics

## Recent Changes (2025-02)

### Trait-Based Refactoring
- **New files**: scoring_traits.jl, scoring_config.jl, psm_container.jl, percolator_generic.jl, pairing.jl, model_training.jl, training_selection.jl, feature_selection.jl, iteration_scheme.jl, mbr_update.jl
- **Removed duplicate functions** from percolatorSortOf.jl:
  - `getIrtBins()` → now in pairing.jl
  - `getIrtBins!()` → now in pairing.jl
  - `assign_random_target_decoy_pairs!()` → now in pairing.jl
  - `assignPairIds!()` → now in pairing.jl
  - `assign_pair_ids()` → now in pairing.jl
- **percolatorSortOf.jl** is now a thin wrapper that builds `ScoringConfig` and calls `percolator_scoring!()`

### Benefits of New Architecture
1. **Composability**: Mix and match algorithm components
2. **Testability**: Each trait can be unit tested independently
3. **Extensibility**: Add new strategies without modifying core code
4. **Type Safety**: Compile-time dispatch based on config types
