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
        targets = collect(get_column(psms, :target))
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
            trace_probs_train = collect(Float32, get_column(psms_train, :trace_prob))
            targets_train = collect(Bool, get_column(psms_train, :target))
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
            MBR_estimates[test_idx] = collect(Float32, get_column(psms_test, :trace_prob))
        else
            prob_test[test_idx] = collect(Float32, get_column(psms_test, :trace_prob))
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
