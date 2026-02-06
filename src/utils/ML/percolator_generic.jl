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
Supports both in-memory and out-of-memory processing via MemoryStrategy trait.
=#

#############################################################################
# Memory Strategy Helper Functions (In-Memory and OOM side-by-side)
#############################################################################

#= Pair ID Sampling - OOM only =#

"""
    scan_pair_ids(file_paths::Vector{String}) -> (pair_counts, pair_to_files)

Scan files to get pair_id -> count mapping and pair_id -> file paths mapping.
Used for proportional sampling in OOM mode.
"""
function scan_pair_ids(file_paths::Vector{String})
    pair_counts = Dict{UInt32, Int}()
    pair_to_files = Dict{UInt32, Vector{String}}()

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        pair_ids = tbl.pair_id
        for pid in pair_ids
            pair_counts[pid] = get(pair_counts, pid, 0) + 1
            files = get!(pair_to_files, pid, String[])
            if isempty(files) || files[end] != fpath
                push!(files, fpath)
            end
        end
    end

    return pair_counts, pair_to_files
end

"""
    sample_pair_ids(pair_counts::Dict{UInt32, Int}, max_psms::Int) -> Set{UInt32}

Sample pair_ids proportionally to achieve target PSM count.
Returns Set of sampled pair_ids. All PSMs with same pair_id are either all
sampled or all excluded (maintains pair integrity for MBR).
"""
function sample_pair_ids(pair_counts::Dict{UInt32, Int}, max_psms::Int)
    total_psms = sum(values(pair_counts))
    if total_psms <= max_psms
        return Set(keys(pair_counts))  # Use all
    end

    # Sample proportionally
    sample_rate = max_psms / total_psms
    sampled = Set{UInt32}()
    rng = Random.MersenneTwister(1776)

    for (pid, _) in pair_counts
        if rand(rng) < sample_rate
            push!(sampled, pid)
        end
    end

    return sampled
end

"""
    load_sampled_psms(file_paths::Vector{String}, sampled_pairs::Set{UInt32}) -> DataFrame

Load only PSMs with sampled pair_ids from files.
"""
function load_sampled_psms(file_paths::Vector{String}, sampled_pairs::Set{UInt32})
    dfs = DataFrame[]
    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        df = DataFrame(tbl)
        mask = [pid in sampled_pairs for pid in df.pair_id]
        if any(mask)
            push!(dfs, df[mask, :])
        end
    end
    return isempty(dfs) ? DataFrame() : vcat(dfs...)
end

#= Global Q-value Computation (OOM) =#

"""
    compute_global_qvalues!(file_paths::Vector{String})

Compute global q-values across all files in a CV fold and write back.
Two-pass approach: first collect all scores/targets, then compute and write back.
"""
function compute_global_qvalues!(file_paths::Vector{String})
    # Collect all scores and targets
    all_scores = Float32[]
    all_targets = Bool[]
    file_ranges = Tuple{Int,Int}[]  # (start, end) for each file

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        start_idx = length(all_scores) + 1
        append!(all_scores, tbl.trace_prob)
        append!(all_targets, tbl.target)
        push!(file_ranges, (start_idx, length(all_scores)))
    end

    # Compute global q-values
    q_values = zeros(Float64, length(all_scores))
    get_qvalues!(all_scores, all_targets, q_values)

    # Write back to each file
    for (i, fpath) in enumerate(file_paths)
        df = DataFrame(Arrow.Table(fpath); copycols=true)
        start_idx, end_idx = file_ranges[i]
        df.q_value = q_values[start_idx:end_idx]
        writeArrow(fpath, df)
    end
end

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
    process_fold_iterations!(phase, psms_view, indices, fold_models, iteration_rounds,
                             config, mbr_start_iter, prob_output, baseline_output, pbar) -> Nothing

Process one fold through all iterations for either training or prediction phase.
Uses trait dispatch to control phase-specific behavior.
"""
function process_fold_iterations!(
    phase::ScoringPhase,
    psms_view::AbstractPSMContainer,
    indices::Vector{Int},
    fold_models::Vector{Any},
    iteration_rounds::Vector{Int},
    config::ScoringConfig,
    mbr_start_iter::Int,
    prob_output::Vector{Float32},
    baseline_output::Vector{Float32},
    pbar
)
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

    return nothing
end

#############################################################################
# Main percolator scoring function
#############################################################################

"""
    percolator_scoring!(psms::AbstractPSMContainer, config::ScoringConfig;
                        memory::MemoryStrategy=InMemoryProcessing(),
                        show_progress::Bool=true, verbose::Bool=false) -> Dict{UInt8, Any}

Generic PSM scoring using configurable traits.

# Arguments
- `psms::AbstractPSMContainer`: PSMs to score (modified in-place for in-memory mode)
- `config::ScoringConfig`: Configuration specifying all algorithm components
- `memory::MemoryStrategy`: Memory strategy (InMemoryProcessing or OutOfMemoryProcessing)
- `show_progress::Bool`: Whether to show progress bar
- `verbose::Bool`: Whether to print verbose logging

# Returns
- `Dict{UInt8, Vector{TrainedModel}}`: Dictionary mapping CV fold to trained models
"""
function percolator_scoring!(psms::AbstractPSMContainer, config::ScoringConfig;
                             memory::MemoryStrategy = InMemoryProcessing(),
                             show_progress::Bool = true, verbose::Bool = false)
    return percolator_scoring_impl!(psms, config, memory; show_progress, verbose)
end

# In-memory implementation (existing code path)
function percolator_scoring_impl!(psms::AbstractPSMContainer, config::ScoringConfig,
                                   ::InMemoryProcessing;
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

    # Step 6: PHASE 1 - Train all models (predict on training data only)
    for test_fold_idx in unique_cv_folds
        train_idx = train_indices[test_fold_idx]
        psms_train = get_view(psms, train_idx)

        fold_models = Vector{Any}(undef, total_iterations)
        process_fold_iterations!(
            TrainingPhase(), psms_train, train_idx, fold_models,
            iteration_rounds, config, mbr_start_iter,
            prob_train, nonMBR_estimates, pbar
        )
        models[test_fold_idx] = fold_models
    end

    # Step 7: PHASE 2 - Apply stored models to held-out test data
    for test_fold_idx in unique_cv_folds
        test_idx = fold_indices[test_fold_idx]
        fold_models = models[test_fold_idx]
        psms_test = get_view(psms, test_idx)

        process_fold_iterations!(
            PredictionPhase(), psms_test, test_idx, fold_models,
            iteration_rounds, config, mbr_start_iter,
            prob_test, nonMBR_estimates, pbar
        )

        # Store final predictions
        if has_mbr_support(config.mbr_update)
            MBR_estimates[test_idx] = collect(Float32, get_column(psms_test, :trace_prob))
        else
            prob_test[test_idx] = collect(Float32, get_column(psms_test, :trace_prob))
        end
    end

    # Step 8: Finalize scoring
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
# Out-of-Memory Implementation
#############################################################################

"""
    percolator_scoring_impl!(psms, config, memory::OutOfMemoryProcessing; kwargs...)

Out-of-memory implementation that processes files one at a time.
Training uses sampled data, prediction is file-by-file with global q-values.
"""
function percolator_scoring_impl!(::AbstractPSMContainer, config::ScoringConfig,
                                   memory::OutOfMemoryProcessing;
                                   show_progress::Bool = true, verbose::Bool = false)
    unique_cv_folds = collect(keys(memory.fold_file_paths))
    models = Dict{UInt8, Vector{Any}}()

    iteration_rounds = get_iteration_rounds(config.iteration_scheme)
    total_iterations = length(iteration_rounds)
    mbr_start_iter = total_iterations

    iterations_to_run = has_mbr_support(config.mbr_update) ? total_iterations : max(mbr_start_iter - 1, 1)

    # Progress tracking - training + prediction for each fold
    total_progress_steps = length(unique_cv_folds) * iterations_to_run * 2
    pbar = show_progress ? ProgressBar(total=total_progress_steps) : nothing

    Random.seed!(1776)

    # Phase 1: Training (one fold at a time)
    for test_fold in unique_cv_folds
        # Get files for training (all folds except test fold)
        train_files = String[]
        for (f, paths) in memory.fold_file_paths
            if f != test_fold
                append!(train_files, paths)
            end
        end

        # Scan and sample pair_ids
        pair_counts, _ = scan_pair_ids(train_files)
        sampled_pairs = sample_pair_ids(pair_counts, memory.max_training_psms)

        # Load sampled PSMs
        df = load_sampled_psms(train_files, sampled_pairs)

        if verbose
            @user_info "OOM Training Fold $test_fold: Sampled $(nrow(df)) PSMs from $(length(train_files)) files"
        end

        if nrow(df) == 0
            @warn "No training data sampled for fold $test_fold"
            continue
        end

        training_psms = DataFramePSMContainer(df, Val(:unsafe))

        # Apply pairing to sampled data (pair_id already assigned during file split)
        # Initialize MBR columns
        initialize_mbr_columns!(training_psms, config.mbr_update)

        # Sort for MBR grouping
        sort_container!(training_psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])

        # Train models - dummy arrays since we don't need outputs in memory
        n_train = nrows(training_psms)
        prob_train = zeros(Float32, n_train)
        nonMBR_estimates = zeros(Float32, n_train)
        train_indices = collect(1:n_train)

        fold_models = Vector{Any}(undef, total_iterations)
        process_fold_iterations!(
            TrainingPhase(), training_psms, train_indices, fold_models,
            iteration_rounds, config, mbr_start_iter,
            prob_train, nonMBR_estimates, pbar
        )
        models[test_fold] = fold_models

        # Free memory
        training_psms = nothing
        df = nothing
        GC.gc()
    end

    # Phase 2: Prediction (file-by-file)
    for test_fold in unique_cv_folds
        fold_models = models[test_fold]
        test_files = memory.fold_file_paths[test_fold]

        if verbose
            @user_info "OOM Prediction Fold $test_fold: Processing $(length(test_files)) files"
        end

        # Pass 1: Predict on each file
        for fpath in test_files
            df = DataFrame(Arrow.Table(fpath); copycols=true)
            psms = DataFramePSMContainer(df, Val(:unsafe))

            # Initialize MBR columns if needed
            if !has_column(psms, :trace_prob)
                initialize_mbr_columns!(psms, config.mbr_update)
            end

            # Apply each iteration's model
            for (itr, _) in enumerate(iteration_rounds)
                if itr > (mbr_start_iter - 1)
                    break
                end

                model = fold_models[itr]
                probs = predict_scores(model, psms)
                set_column!(psms, :trace_prob, probs)
            end

            # Write predictions back (use writeArrow to avoid mmap deadlock)
            writeArrow(fpath, to_dataframe(psms))

            # Update progress
            !isnothing(pbar) && update(pbar)
        end

        # Pass 2: Compute global q-values across all test files
        @user_info "Computing global q-values for fold $(test_fold)..."
        compute_global_qvalues!(test_files)
        @user_info "  Global q-values complete"

        # MBR: Collect pair statistics, then apply
        if has_mbr_support(config.mbr_update)
            @user_info "Computing MBR features (this may take a while)..."
            compute_and_apply_mbr_features_oom!(test_files, config.mbr_update)
            @user_info "  MBR features complete"

            @user_info "Finalizing MBR scoring..."
            finalize_mbr_scoring_oom!(test_files, fold_models[mbr_start_iter], config.mbr_update)
            @user_info "  MBR scoring finalized"
        end
    end

    return models
end

"""
    compute_and_apply_mbr_features_oom!(file_paths::Vector{String}, strategy::PairBasedMBR)

Compute MBR features across files using streaming approach.
Two-pass: first collect pair statistics, then apply features.
"""
function compute_and_apply_mbr_features_oom!(file_paths::Vector{String}, strategy::PairBasedMBR)
    # Collect pair statistics across all files
    pair_stats = collect_pair_statistics_oom(file_paths, strategy.q_cutoff)

    # Apply MBR features to each file
    for (file_num, fpath) in enumerate(file_paths)
        @user_info "  Applying MBR features: file $file_num/$(length(file_paths))"
        df = DataFrame(Arrow.Table(fpath); copycols=true)
        apply_mbr_features_from_stats!(df, pair_stats, strategy.q_cutoff)
        writeArrow(fpath, df)
    end
end

function compute_and_apply_mbr_features_oom!(::Vector{String}, ::NoMBR)
    return nothing
end

"""
    finalize_mbr_scoring_oom!(file_paths, mbr_model, strategy::PairBasedMBR)

Finalize MBR scoring for OOM pipeline. Mirrors the in-memory `finalize_scoring!` logic:
1. Collect baseline `trace_prob` and `q_value` across all files (written by Pass 1 + Pass 2)
2. Compute pass mask and probability threshold from baseline q-values
3. Re-predict with MBR-aware model, write `MBR_boosted_trace_prob` and `MBR_transfer_candidate`
"""
function finalize_mbr_scoring_oom!(file_paths::Vector{String}, mbr_model, strategy::PairBasedMBR)
    # Collect baseline predictions and q-values across all files
    all_baseline_probs = Float32[]
    all_qvals = Float64[]
    file_ranges = Tuple{Int,Int}[]

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        start_idx = length(all_baseline_probs) + 1
        append!(all_baseline_probs, tbl.trace_prob)
        append!(all_qvals, tbl.q_value)
        push!(file_ranges, (start_idx, length(all_baseline_probs)))
    end

    # Determine passing threshold from baseline (mirrors finalize_scoring!)
    pass_mask = all_qvals .<= strategy.q_cutoff
    has_passing = any(pass_mask)
    prob_thresh = has_passing ? minimum(all_baseline_probs[pass_mask]) : typemax(Float32)

    # Re-predict with MBR model and write finalized columns
    for (i, fpath) in enumerate(file_paths)
        df = DataFrame(Arrow.Table(fpath); copycols=true)
        psms = DataFramePSMContainer(df, Val(:unsafe))

        # Re-predict using MBR-aware model (now with real MBR features)
        mbr_probs = predict_scores(mbr_model, psms)

        # Set MBR_boosted_trace_prob (trace_prob stays as baseline)
        df[!, :MBR_boosted_trace_prob] = mbr_probs

        # Compute transfer candidates for this file's rows
        start_idx, end_idx = file_ranges[i]
        file_pass_mask = pass_mask[start_idx:end_idx]
        transfer_candidates = .!file_pass_mask .& (df.MBR_max_pair_prob .>= prob_thresh)
        df[!, :MBR_transfer_candidate] = transfer_candidates

        writeArrow(fpath, df)
    end
end

function finalize_mbr_scoring_oom!(::Vector{String}, ::Any, ::NoMBR)
    return nothing
end

"""
    PairStatisticsOOM

Running statistics for MBR feature computation in OOM mode.
Tracks the best 2 PSMs per pair group for MBR transfer.
"""
mutable struct PairStatisticsOOM
    # Best PSM info
    best_prob_1::Float32
    best_ms_file_idx_1::UInt32
    best_log2_weights_1::Vector{Float32}
    best_irts_1::Vector{Float32}
    best_irt_residual_1::Float32
    best_weight_1::Float32
    best_log2_intensity_explained_1::Float32
    is_best_decoy_1::Bool

    # Second best PSM info
    best_prob_2::Float32
    best_ms_file_idx_2::UInt32
    best_log2_weights_2::Vector{Float32}
    best_irts_2::Vector{Float32}
    best_irt_residual_2::Float32
    best_weight_2::Float32
    best_log2_intensity_explained_2::Float32
    is_best_decoy_2::Bool

    # Aggregate stats
    num_runs_passing::Int32
    unique_passing_runs::Set{UInt32}
end

function PairStatisticsOOM()
    return PairStatisticsOOM(
        -Inf32, zero(UInt32), Float32[], Float32[], 0f0, 0f0, 0f0, true,
        -Inf32, zero(UInt32), Float32[], Float32[], 0f0, 0f0, 0f0, true,
        zero(Int32), Set{UInt32}()
    )
end

"""
    collect_pair_statistics_oom(file_paths::Vector{String}, q_cutoff::Float32)

Collect running statistics per (pair_id, isotopes_captured) across files.
"""
function collect_pair_statistics_oom(file_paths::Vector{String}, q_cutoff::Float32)
    # Key: (pair_id, isotopes_captured)
    pair_stats = Dict{Tuple{UInt32, Tuple{Int8, Int8}}, PairStatisticsOOM}()

    for (file_num, fpath) in enumerate(file_paths)
        @user_info "  Collecting pair statistics: file $file_num/$(length(file_paths))"
        tbl = Arrow.Table(fpath)
        n = length(tbl.pair_id)

        for i in 1:n
            pair_id = tbl.pair_id[i]
            isotopes = tbl.isotopes_captured[i]
            key = (pair_id, isotopes)

            stats = get!(pair_stats, key, PairStatisticsOOM())

            prob = tbl.trace_prob[i]
            ms_file_idx = tbl.ms_file_idx[i]
            q_val = tbl.q_value[i]

            # Track passing runs
            if q_val <= q_cutoff
                push!(stats.unique_passing_runs, ms_file_idx)
            end

            # Update best/second-best
            if prob > stats.best_prob_1
                # Shift best_1 to best_2
                stats.best_prob_2 = stats.best_prob_1
                stats.best_ms_file_idx_2 = stats.best_ms_file_idx_1
                stats.best_log2_weights_2 = stats.best_log2_weights_1
                stats.best_irts_2 = stats.best_irts_1
                stats.best_irt_residual_2 = stats.best_irt_residual_1
                stats.best_weight_2 = stats.best_weight_1
                stats.best_log2_intensity_explained_2 = stats.best_log2_intensity_explained_1
                stats.is_best_decoy_2 = stats.is_best_decoy_1

                # Update best_1
                stats.best_prob_1 = prob
                stats.best_ms_file_idx_1 = ms_file_idx
                stats.best_log2_weights_1 = log2.(Vector{Float32}(tbl.weights[i]))
                stats.best_irts_1 = Vector{Float32}(tbl.irts[i])
                stats.best_irt_residual_1 = Float32(tbl.irt_pred[i] - tbl.irt_obs[i])
                stats.best_weight_1 = tbl.weight[i]
                stats.best_log2_intensity_explained_1 = tbl.log2_intensity_explained[i]
                stats.is_best_decoy_1 = tbl.decoy[i]

            elseif prob > stats.best_prob_2
                # Update best_2
                stats.best_prob_2 = prob
                stats.best_ms_file_idx_2 = ms_file_idx
                stats.best_log2_weights_2 = log2.(Vector{Float32}(tbl.weights[i]))
                stats.best_irts_2 = Vector{Float32}(tbl.irts[i])
                stats.best_irt_residual_2 = Float32(tbl.irt_pred[i] - tbl.irt_obs[i])
                stats.best_weight_2 = tbl.weight[i]
                stats.best_log2_intensity_explained_2 = tbl.log2_intensity_explained[i]
                stats.is_best_decoy_2 = tbl.decoy[i]
            end
        end
    end

    # Compute num_runs_passing for each pair
    for (_, stats) in pair_stats
        stats.num_runs_passing = Int32(length(stats.unique_passing_runs))
    end

    return pair_stats
end

"""
    apply_mbr_features_from_stats!(df::DataFrame, pair_stats, q_cutoff::Float32)

Apply MBR features to a DataFrame using pre-computed pair statistics.
"""
function apply_mbr_features_from_stats!(df::DataFrame, pair_stats::Dict, q_cutoff::Float32)
    n = nrow(df)

    # Initialize MBR columns if missing
    if !hasproperty(df, :MBR_max_pair_prob)
        df.MBR_max_pair_prob = zeros(Float32, n)
        df.MBR_best_irt_diff = zeros(Float32, n)
        df.MBR_log2_weight_ratio = zeros(Float32, n)
        df.MBR_log2_explained_ratio = zeros(Float32, n)
        df.MBR_rv_coefficient = zeros(Float32, n)
        df.MBR_is_best_decoy = trues(n)
        df.MBR_num_runs = zeros(Int32, n)
        df.MBR_is_missing = falses(n)
    end

    for i in 1:n
        pair_id = df.pair_id[i]
        isotopes = df.isotopes_captured[i]
        key = (pair_id, isotopes)

        if !haskey(pair_stats, key)
            df.MBR_is_missing[i] = true
            continue
        end

        stats = pair_stats[key]
        ms_file_idx = df.ms_file_idx[i]

        # Count passing runs excluding current file
        current_passes = df.q_value[i] <= q_cutoff
        df.MBR_num_runs[i] = stats.num_runs_passing - (current_passes && ms_file_idx in stats.unique_passing_runs ? 1 : 0)

        # Get best PSM from a different run
        best_log2_weights = Float32[]
        best_irts = Float32[]
        best_weight = 0f0
        best_log2_ie = 0f0
        best_residual = 0f0
        best_prob = -1f0
        MBR_is_best_decoy = true

        if stats.best_ms_file_idx_1 != ms_file_idx && !isempty(stats.best_log2_weights_1)
            best_log2_weights = stats.best_log2_weights_1
            best_irts = stats.best_irts_1
            best_weight = stats.best_weight_1
            best_log2_ie = stats.best_log2_intensity_explained_1
            best_residual = stats.best_irt_residual_1
            best_prob = stats.best_prob_1
            MBR_is_best_decoy = stats.is_best_decoy_1
        elseif stats.best_ms_file_idx_2 != ms_file_idx && !isempty(stats.best_log2_weights_2)
            best_log2_weights = stats.best_log2_weights_2
            best_irts = stats.best_irts_2
            best_weight = stats.best_weight_2
            best_log2_ie = stats.best_log2_intensity_explained_2
            best_residual = stats.best_irt_residual_2
            best_prob = stats.best_prob_2
            MBR_is_best_decoy = stats.is_best_decoy_2
        else
            # No valid MBR match
            df.MBR_best_irt_diff[i] = -1f0
            df.MBR_rv_coefficient[i] = -1f0
            df.MBR_is_best_decoy[i] = true
            df.MBR_max_pair_prob[i] = -1f0
            df.MBR_log2_weight_ratio[i] = -1f0
            df.MBR_log2_explained_ratio[i] = -1f0
            df.MBR_is_missing[i] = true
            continue
        end

        df.MBR_max_pair_prob[i] = best_prob
        df.MBR_is_best_decoy[i] = MBR_is_best_decoy

        # Compute MBR features
        current_log2_weights = log2.(df.weights[i])
        current_irts = df.irts[i]
        current_residual = Float32(df.irt_pred[i] - df.irt_obs[i])

        df.MBR_best_irt_diff[i] = abs(best_residual - current_residual)
        df.MBR_log2_weight_ratio[i] = log2(df.weight[i] / best_weight)
        df.MBR_log2_explained_ratio[i] = df.log2_intensity_explained[i] - best_log2_ie

        # RV coefficient
        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, current_log2_weights)
        best_irts_padded, irts_padded = pad_rt_equal_length(best_irts, current_irts)
        df.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_irts_padded, weights_padded, irts_padded)
    end
end

#############################################################################
# DataFrame convenience wrapper
#############################################################################

"""
    percolator_scoring!(psms::DataFrame, config::ScoringConfig; kwargs...) -> Dict

Convenience wrapper that wraps DataFrame in DataFramePSMContainer.
"""
function percolator_scoring!(psms::DataFrame, config::ScoringConfig;
                             memory::MemoryStrategy = InMemoryProcessing(),
                             show_progress::Bool = true, verbose::Bool = false)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    return percolator_scoring!(container, config; memory, show_progress, verbose)
end
