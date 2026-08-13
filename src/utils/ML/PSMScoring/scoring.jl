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
Generic percolator scoring function using trait-based abstraction.
Uses AbstractPSMContainer for data access.
=#

#############################################################################
# Phase-dispatched helper functions
#############################################################################

"""
    get_or_train_model(phase, fold_models, itr, config, psms_itr, features, num_round) -> Model

Get or train a model depending on the phase.
- TrainingPhase: Train a new model and store it
- PredictionPhase: Retrieve pre-trained model
"""
function get_or_train_model(::TrainingPhase, fold_models::Vector{Any}, itr::Int,
                            config::ScoringConfig, psms_itr::AbstractPSMContainer,
                            features::Vector{Symbol}, num_round::Int)
    model = train_model(config.model, psms_itr, features, num_round)
    fold_models[itr] = model
    return model
end

function get_or_train_model(::PredictionPhase, fold_models::Vector{Any}, itr::Int,
                            ::ScoringConfig, ::AbstractPSMContainer,
                            ::Vector{Symbol}, ::Int)
    return fold_models[itr]
end

"""
    get_training_subset(phase, psms, config, itr, total_iterations) -> AbstractPSMContainer

Get training data subset for model training.
- TrainingPhase: Apply training data selection strategy
- PredictionPhase: Return all data (no subset needed)
"""
function get_training_subset(::TrainingPhase, psms::AbstractPSMContainer,
                             config::ScoringConfig, itr::Int, total_iterations::Int)
    psms_itr = select_training_data(psms, config.training_data, itr)

    # Get features and filter to available
    features = get_features(config.feature_selection, itr, total_iterations)
    features = filter_available_features(features, psms_itr)

    return psms_itr, features
end

function get_training_subset(::PredictionPhase, psms::AbstractPSMContainer,
                             config::ScoringConfig, itr::Int, total_iterations::Int)
    # For prediction, no training subset needed - features only matter for shape
    features = get_features(config.feature_selection, itr, total_iterations)
    features = filter_available_features(features, psms)
    return psms, features
end

#############################################################################
# Unified fold iteration processing
#############################################################################

"""
    process_fold!(phase, workspace, fold, num_rounds, config, pbar) -> Nothing

Train (or apply) a single model for one CV fold.
"""
function process_fold!(
    phase::ScoringPhase,
    workspace::AbstractScoringWorkspace,
    fold::UInt8,
    num_rounds::Int,
    config::ScoringConfig,
    pbar
)
    psms_view = get_phase_view(workspace, phase, fold)
    indices = get_phase_indices(workspace, phase, fold)
    fold_models = init_fold_models(workspace, phase, fold, 1)
    prob_output = get_phase_output(workspace, phase)

    # Get training subset (TrainingPhase) or full data (PredictionPhase)
    psms_itr, features = get_training_subset(phase, psms_view, config, 1, 1)

    # Train or retrieve model
    model = get_or_train_model(phase, fold_models, 1, config, psms_itr, features, num_rounds)

    # Predict on PSM view
    prob_output[indices] = predict_scores(model, psms_view)
    set_column!(psms_view, :trace_prob, prob_output[indices])

    !isnothing(pbar) && update(pbar)

    commit_fold!(workspace, phase, fold, fold_models)
    return nothing
end

#############################################################################
# Main percolator scoring function
#############################################################################

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
                             show_progress::Bool = true, verbose::Bool = false,
                             print_importances::Bool = false)

    # Step 1: Prepare training data (pairing, sorting, column init)
    prepare_training_data!(psms, config)

    # Step 2: Single training pass
    num_rounds = config.num_rounds

    # Steps 4-5: Setup workspace (CV folds + output arrays)
    workspace = setup_scoring_workspace(psms, config)

    # Log dataset info
    if verbose
        resolved_psms = get_psms(workspace)
        targets = collect(get_column(resolved_psms, :target))
        n_targets = sum(targets)
        n_decoys = sum(.!targets)
        @debug_l1 "ML Training Dataset: $n_targets targets, $n_decoys decoys (total: $(nrows(resolved_psms)) PSMs)"
    end

    # Progress tracking (2 phases × n_folds)
    total_progress_steps = 2 * length(get_cv_folds(workspace))
    pbar = show_progress ? ProgressBar(total=total_progress_steps) : nothing

    Random.seed!(1776)

    # Step 6: PHASE 1 - Train all models (predict on training data only)
    for fold in get_cv_folds(workspace)
        process_fold!(TrainingPhase(), workspace, fold, num_rounds, config, pbar)
    end

    # Print feature importances from the first fold
    if print_importances
        first_fold = first(get_cv_folds(workspace))
        fold_models = get_fold_models(workspace, first_fold)
        imp = importance(fold_models[1])
        if imp !== nothing
            sorted_imp = sort(imp, by=x -> x[2], rev=true)
            @debug_l1 "LightGBM Feature Importances (gain, fold $first_fold):"
            for (fname, gain) in sorted_imp
                @debug_l1 "  $(rpad(fname, 40)) $(round(gain, digits=2))"
            end
        end
    end

    # Step 7: PHASE 2 - Predict on held-out test data
    for fold in get_cv_folds(workspace)
        process_fold!(PredictionPhase(), workspace, fold, num_rounds, config, pbar)
        store_final_predictions!(workspace, fold)
    end

    # Step 8: Finalize scoring
    finalize_scoring!(workspace)

    return get_all_models(workspace)
end

#############################################################################
# Workspace-level finalize_scoring! dispatch
#############################################################################

"""
    finalize_scoring!(workspace)

Workspace-level dispatch for finalization. Sets trace_prob = prob_test.
"""
function finalize_scoring!(ws::InMemoryScoringWorkspace)
    set_column!(get_psms(ws), :trace_prob, get_test_output(ws))
end

#############################################################################
# DataFrame convenience wrapper
#############################################################################

"""
    percolator_scoring!(psms::DataFrame, config::ScoringConfig; kwargs...) -> Dict

Convenience wrapper that wraps DataFrame in DataFramePSMContainer.
"""
function percolator_scoring!(psms::DataFrame, config::ScoringConfig;
                             show_progress::Bool = true, verbose::Bool = false,
                             print_importances::Bool = false)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    return percolator_scoring!(container, config; show_progress, verbose, print_importances)
end
