#=
ScoringWorkspace abstraction for CV fold setup and output array management.
Dispatches on container type to allow future OOM implementations.
=#

#############################################################################
# Abstract Type
#############################################################################

"""
    AbstractScoringWorkspace

Abstract type for scoring workspace management.
Bundles CV fold information, output arrays, and model storage.

Future OOM implementations can specialize this to map folds to file paths
and write output per-file rather than to global arrays.
"""
abstract type AbstractScoringWorkspace end

#############################################################################
# In-Memory Implementation
#############################################################################

"""
    InMemoryScoringWorkspace <: AbstractScoringWorkspace

In-memory workspace that pre-computes fold indices and allocates global output arrays.
"""
struct InMemoryScoringWorkspace <: AbstractScoringWorkspace
    psms::DataFramePSMContainer
    cv_folds::Vector{UInt8}
    fold_indices::Dict{UInt8, Vector{Int}}
    train_indices::Dict{UInt8, Vector{Int}}
    models::Dict{UInt8, Vector{Any}}
    prob_test::Vector{Float32}
    prob_train::Vector{Float32}
    nonMBR_estimates::Vector{Float32}
    MBR_estimates::Vector{Float32}
end

#############################################################################
# Factory Function (dispatch point for future OOM)
#############################################################################

"""
    setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)

Create a scoring workspace dispatched on container type.
For in-memory containers, pre-computes fold/train indices and allocates output arrays.
"""
function setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)
    cv_folds = get_cv_folds(psms)
    fold_idx = Dict(fold => get_fold_indices(psms, fold) for fold in cv_folds)
    train_idx = Dict(fold => get_train_indices(psms, fold) for fold in cv_folds)
    n = nrows(psms)
    return InMemoryScoringWorkspace(
        psms, cv_folds, fold_idx, train_idx,
        Dict{UInt8, Vector{Any}}(),
        zeros(Float32, n), zeros(Float32, n),
        zeros(Float32, n), zeros(Float32, n)
    )
end

#############################################################################
# Training Data Preparation (dispatch point for future OOM sampling)
#############################################################################

"""
    prepare_training_data(psms::DataFramePSMContainer, config::ScoringConfig)

Prepare training data from the PSM container.
In-memory: returns psms unchanged (all data already available).
Future OOM: would sample PSMs from files and return a DataFramePSMContainer.
"""
function prepare_training_data(psms::DataFramePSMContainer, ::ScoringConfig)
    return psms
end

#############################################################################
# Interface Functions (all overridable for future OOM)
#############################################################################

# CV fold iteration
get_cv_folds(ws::InMemoryScoringWorkspace) = ws.cv_folds

# Training phase access
get_training_view(ws::InMemoryScoringWorkspace, fold::UInt8) = get_view(ws.psms, ws.train_indices[fold])
get_train_indices(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.train_indices[fold]

# Prediction phase access
get_test_view(ws::InMemoryScoringWorkspace, fold::UInt8) = get_view(ws.psms, ws.fold_indices[fold])
get_test_indices(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.fold_indices[fold]

# Output arrays
get_train_output(ws::InMemoryScoringWorkspace) = ws.prob_train
get_test_output(ws::InMemoryScoringWorkspace) = ws.prob_test
get_baseline_output(ws::InMemoryScoringWorkspace) = ws.nonMBR_estimates
get_mbr_output(ws::InMemoryScoringWorkspace) = ws.MBR_estimates

# Model storage
store_fold_models!(ws::InMemoryScoringWorkspace, fold::UInt8, models::Vector{Any}) = (ws.models[fold] = models)
get_fold_models(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.models[fold]
get_all_models(ws::InMemoryScoringWorkspace) = ws.models

# Container access (for finalize_scoring! and other top-level ops)
get_psms(ws::InMemoryScoringWorkspace) = ws.psms
