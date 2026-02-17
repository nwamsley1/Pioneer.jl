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
    training_strategy = itr == 1 ? AllDataSelection() : config.training_data
    psms_itr = select_training_data(psms, training_strategy, itr)

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

"""
    compute_and_set_qvalues!(phase, psms) -> Nothing

Compute and set q-values on PSM data.
- TrainingPhase: Compute q-values (needed for select_training_data in next iteration)
- PredictionPhase: No-op (q-values computed inside update_mbr_features_test_only!)
"""
function compute_and_set_qvalues!(::TrainingPhase, psms::AbstractPSMContainer)
    trace_probs = collect(Float32, get_column(psms, :trace_prob))
    targets = collect(Bool, get_column(psms, :target))
    q_values = zeros(Float64, nrows(psms))
    get_qvalues!(trace_probs, targets, q_values)
    set_column!(psms, :q_value, q_values)
    return nothing
end

function compute_and_set_qvalues!(::PredictionPhase, ::AbstractPSMContainer)
    return nothing
end

"""
    update_mbr!(phase, psms, itr, mbr_start_iter, strategy) -> Nothing

Update MBR features based on phase.
- TrainingPhase: Uses update_mbr_features_train_only!
- PredictionPhase: Uses update_mbr_features_test_only!
"""
function update_mbr!(::TrainingPhase, psms::AbstractPSMContainer,
                     itr::Int, mbr_start_iter::Int, strategy::MBRUpdateStrategy)
    update_mbr_features_train_only!(psms, itr, mbr_start_iter, strategy)
    return nothing
end

function update_mbr!(::PredictionPhase, psms::AbstractPSMContainer,
                     itr::Int, mbr_start_iter::Int, strategy::MBRUpdateStrategy)
    update_mbr_features_test_only!(psms, itr, mbr_start_iter, strategy)
    return nothing
end

"""
    store_baseline!(phase, baseline, indices, probs, itr, mbr_start_iter) -> Nothing

Store baseline predictions before MBR iteration.
- TrainingPhase: No-op
- PredictionPhase: Store at mbr_start_iter - 1
"""
function store_baseline!(::TrainingPhase, ::Vector{Float32}, ::Vector{Int},
                         ::Vector{Float32}, ::Int, ::Int)
    return nothing
end

function store_baseline!(::PredictionPhase, baseline::Vector{Float32}, indices::Vector{Int},
                         probs::Vector{Float32}, itr::Int, mbr_start_iter::Int)
    if itr == mbr_start_iter - 1
        baseline[indices] .= probs[indices]
    end
    return nothing
end

#############################################################################
# Unified fold iteration processing
#############################################################################

"""
    process_fold_iterations!(phase, workspace, fold, iteration_rounds,
                             config, mbr_start_iter, pbar) -> Nothing

Process one fold through all iterations for either training or prediction phase.
Uses trait dispatch to control phase-specific behavior. Fold-specific data
(PSM views, indices, models, output arrays) is accessed via phase-dispatched
workspace accessors.
"""
function process_fold_iterations!(
    phase::ScoringPhase,
    workspace::AbstractScoringWorkspace,
    fold::UInt8,
    iteration_rounds::Vector{Int},
    config::ScoringConfig,
    mbr_start_iter::Int,
    pbar
)
    psms_view = get_phase_view(workspace, phase, fold)
    indices = get_phase_indices(workspace, phase, fold)
    fold_models = init_fold_models(workspace, phase, fold, length(iteration_rounds))
    prob_output = get_phase_output(workspace, phase)
    baseline_output = get_baseline_output(workspace)

    total_iterations = length(iteration_rounds)

    for (itr, num_round) in enumerate(iteration_rounds)
        # Get training subset (TrainingPhase) or full data (PredictionPhase)
        psms_itr, features = get_training_subset(phase, psms_view, config, itr, total_iterations)

        # Train or retrieve model
        model = get_or_train_model(phase, fold_models, itr, config, psms_itr, features, num_round)

        # Predict on PSM view
        prob_output[indices] = predict_scores(model, psms_view)
        set_column!(psms_view, :trace_prob, prob_output[indices])

        # Store baseline before MBR iteration (PredictionPhase only)
        store_baseline!(phase, baseline_output, indices, prob_output, itr, mbr_start_iter)

        # Compute q-values (TrainingPhase only - needed for next iteration's select_training_data)
        compute_and_set_qvalues!(phase, psms_view)

        # Update MBR features
        update_mbr!(phase, psms_view, itr, mbr_start_iter, config.mbr_update)

        # Update progress
        !isnothing(pbar) && update(pbar)

        # Early exit for non-MBR mode
        if !has_mbr_support(config.mbr_update) && itr == (mbr_start_iter - 1)
            break
        end
    end

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
                             show_progress::Bool = true, verbose::Bool = false)

    # Step 1: Prepare training data (pairing, sorting, column init — dispatched on container type)
    prepare_training_data!(psms, config)

    # Step 2: Get iteration scheme
    iteration_rounds = get_iteration_rounds(config.iteration_scheme)
    total_iterations = length(iteration_rounds)
    mbr_start_iter = total_iterations  # MBR features used in last iteration

    # Determine actual iterations to run
    iterations_to_run = has_mbr_support(config.mbr_update) ? total_iterations : max(mbr_start_iter - 1, 1)

    # Steps 4-5: Setup workspace (CV folds + output arrays, dispatched on container type)
    workspace = setup_scoring_workspace(psms, config)

    # Log dataset info (after workspace setup so we can use the resolved container)
    if verbose
        resolved_psms = get_psms(workspace)
        targets = collect(get_column(resolved_psms, :target))
        n_targets = sum(targets)
        n_decoys = sum(.!targets)
        @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys (total: $(nrows(resolved_psms)) PSMs)"
    end

    # Progress tracking (2 phases × n_folds × iterations_per_fold)
    total_progress_steps = 2 * length(get_cv_folds(workspace)) * iterations_to_run
    pbar = show_progress ? ProgressBar(total=total_progress_steps) : nothing

    Random.seed!(1776)

    # Step 6: PHASE 1 - Train all models (predict on training data only)
    for fold in get_cv_folds(workspace)
        process_fold_iterations!(TrainingPhase(), workspace, fold,
                                 iteration_rounds, config, mbr_start_iter, pbar)
    end

    # Step 7: PHASE 2 - Predict on held-out test data (dispatched on phase + workspace type)
    for fold in get_cv_folds(workspace)
        process_fold_iterations!(PredictionPhase(), workspace, fold,
                                 iteration_rounds, config, mbr_start_iter, pbar)
        store_final_predictions!(workspace, fold, config.mbr_update)
    end

    # Step 8: Finalize scoring (dispatched on workspace type)
    finalize_scoring!(workspace, config.mbr_update)

    return get_all_models(workspace)
end

"""
    finalize_scoring!(psms, mbr_strategy, prob_test, nonMBR_estimates, MBR_estimates)

Finalize PSM scoring by setting appropriate probability columns and MBR transfer candidates.
"""
function finalize_scoring!(psms::AbstractPSMContainer, strategy::PairBasedMBR,
                           ::Vector{Float32}, nonMBR_estimates::Vector{Float32},
                           MBR_estimates::Vector{Float32})
    # Compute q-values on non-MBR baseline
    targets = collect(Bool, get_column(psms, :target))
    qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
    get_qvalues!(nonMBR_estimates, targets, qvals_prev)
    pass_mask = (qvals_prev .<= strategy.q_cutoff)

    # Find passing probability threshold
    has_passing = !isempty(pass_mask) && any(pass_mask)
    prob_thresh = has_passing ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)

    # Mark transfer candidates
    mbr_max_pair_prob = collect(Float32, get_column(psms, :MBR_max_pair_prob))
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
# Workspace-level finalize_scoring! dispatch
#############################################################################

"""
    finalize_scoring!(workspace, mbr_strategy)

Workspace-level dispatch for finalization. InMemoryScoringWorkspace delegates
to the psms-level method; ArrowFileScoringWorkspace is defined in scoring_workspace.jl.
"""
function finalize_scoring!(ws::InMemoryScoringWorkspace, mbr_strategy)
    finalize_scoring!(get_psms(ws), mbr_strategy,
                      get_test_output(ws), get_baseline_output(ws),
                      get_mbr_output(ws))
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
