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

"""
Machine learning-based global probability prediction for precursors.

This module replaces the simple logodds calculation with an adaptive ML approach
that trains LightGBM models to predict target/decoy labels using cross-run features.
"""

#==========================================================
Phase 1: Feature Calculation Helper Functions
==========================================================#

"""
    calculate_top_n_features(prec_probs::AbstractVector{Float32}, n_top_scores::Int) -> Vector{Float32}

Extract top N scores from a vector of probabilities, padding with 0.0 if needed.

# Arguments
- `prec_probs`: Vector of precursor probability scores across runs
- `n_top_scores`: Number of top scores to extract

# Returns
- Vector of length `n_top_scores` with top scores (descending), padded with 0.0 if needed
"""
function calculate_top_n_features(prec_probs::AbstractVector{Float32}, n_top_scores::Int)
    sorted_scores = sort(prec_probs, rev=true)
    top_n = Vector{Float32}(undef, n_top_scores)

    for i in 1:n_top_scores
        top_n[i] = i <= length(sorted_scores) ? sorted_scores[i] : 0.0f0
    end

    return top_n
end

"""
    calculate_precursor_features(prec_probs::AbstractVector{Float32},
                                 sqrt_n_runs::Int) -> NamedTuple

Calculate core features for a single precursor from its prec_prob values.

# Arguments
- `prec_probs`: Vector of precursor probability scores across runs
- `sqrt_n_runs`: Square root of number of runs (for logodds calculation)

# Returns
NamedTuple with 7 core features:
- `logodds_score`: Baseline logodds probability
- `n_observations`: Number of observations (runs) for this precursor
- `median_score`: Median probability (robust central tendency)
- `std_score`: Standard deviation (spread)
- `max_score`: Maximum probability (best evidence)
- `skewness`: Distribution asymmetry
- `kurtosis`: Distribution tail behavior
"""
function calculate_precursor_features(
    prec_probs::AbstractVector{Float32},
    sqrt_n_runs::Int
)
    # Core statistics
    n_obs = length(prec_probs)
    median_score = median(prec_probs)
    std_score = std(prec_probs)
    max_score = maximum(prec_probs)

    # Distribution shape features
    skew = skewness(prec_probs)
    kurt = kurtosis(prec_probs)

    # Logodds baseline
    logodds_score = logodds(prec_probs, sqrt_n_runs)

    return (
        logodds_score = logodds_score,
        n_observations = n_obs,
        median_score = median_score,
        std_score = std_score,
        max_score = max_score,
        skewness = skew,
        kurtosis = kurt
    )
end

#==========================================================
Phase 2: Feature Engineering Pipeline
==========================================================#

"""
    build_precursor_features_df(merged_df::DataFrame,
                                sqrt_n_runs::Int,
                                score_column::Symbol = :prec_prob) -> DataFrame

Build precursor-level feature DataFrame from merged PSM data using groupby.

# Arguments
- `merged_df`: DataFrame with columns [:precursor_idx, score_column, :target, :cv_fold]
- `sqrt_n_runs`: Square root of number of runs (for logodds calculation)
- `score_column`: Column name to use for scoring (:prec_prob or :MBR_boosted_prec_prob)

# Returns
DataFrame with one row per unique precursor, columns for 7 core features
"""
function build_precursor_features_df(
    merged_df::DataFrame,
    sqrt_n_runs::Int,
    score_column::Symbol = :prec_prob
)::DataFrame
    # Use groupby to get all score values for each precursor
    # This is similar to how current Step 5 works with transform!
    grouped = groupby(merged_df, :precursor_idx)

    # Pre-allocate result DataFrame with 7 core features
    n_precursors = length(grouped)
    precursor_features = DataFrame(
        precursor_idx = Vector{UInt32}(undef, n_precursors),
        target = Vector{Bool}(undef, n_precursors),
        cv_fold = Vector{UInt8}(undef, n_precursors),
        logodds_score = Vector{Float32}(undef, n_precursors),
        n_observations = Vector{Int}(undef, n_precursors),
        median_score = Vector{Float32}(undef, n_precursors),
        std_score = Vector{Float32}(undef, n_precursors),
        max_score = Vector{Float32}(undef, n_precursors),
        skewness = Vector{Float32}(undef, n_precursors),
        kurtosis = Vector{Float32}(undef, n_precursors)
    )

    # Process each group (precursor)
    for (idx, group) in enumerate(grouped)
        # Metadata (same for all rows in group)
        precursor_features.precursor_idx[idx] = first(group.precursor_idx)
        precursor_features.target[idx] = first(group.target)
        precursor_features.cv_fold[idx] = first(group.cv_fold)

        # Calculate features from score_column values (prec_prob or MBR_boosted_prec_prob)
        features = calculate_precursor_features(
            getproperty(group, score_column),
            sqrt_n_runs
        )

        # Assign 7 core features to DataFrame
        precursor_features.logodds_score[idx] = features.logodds_score
        precursor_features.n_observations[idx] = features.n_observations
        precursor_features.median_score[idx] = features.median_score
        precursor_features.std_score[idx] = features.std_score
        precursor_features.max_score[idx] = features.max_score
        precursor_features.skewness[idx] = features.skewness
        precursor_features.kurtosis[idx] = features.kurtosis
    end

    return precursor_features
end

"""
    get_training_data_for_global_prob_iteration!(
        precursor_features::DataFrame,
        itr::Int,
        global_prob_iter1::Union{Nothing, Vector{Float32}},
        max_global_qval_threshold::Float32,
        min_PEP_neg_threshold::Float32
    ) -> DataFrame

Filter and refine training data for iterative global probability training.

# Arguments
- `precursor_features`: DataFrame with precursor-level features and labels
- `itr`: Current iteration number (1 or 2)
- `global_prob_iter1`: Predictions from iteration 1 (nothing for iteration 1)
- `max_global_qval_threshold`: Max q-value for keeping targets (e.g., 0.01)
- `min_PEP_neg_threshold`: Min PEP for relabeling targets as negatives (e.g., 0.5)

# Returns
DataFrame with refined training set for current iteration

# Iteration 1
Returns copy of all data (no filtering)

# Iteration 2
1. Calculate PEP from global_prob_iter1
2. Relabel targets with PEP ≥ threshold as negatives (uncertain targets)
3. Filter to keep:
   - All original decoys
   - All relabeled uncertain targets (now negatives)
   - High-confidence targets (global_qval ≤ threshold)
"""
function get_training_data_for_global_prob_iteration!(
    precursor_features::DataFrame,
    itr::Int,
    global_prob_iter1::Union{Nothing, Vector{Float32}},
    max_global_qval_threshold::Float32,
    min_PEP_neg_threshold::Float32
)
    if itr == 1
        # Iteration 1: Train on all data
        return copy(precursor_features)
    else
        # Iteration 2: Refine training set
        @assert !isnothing(global_prob_iter1) "global_prob_iter1 required for iteration 2"

        # Create working copy to avoid modifying original
        features_refined = copy(precursor_features)

        # Add iteration 1 predictions to calculate PEP
        features_refined.global_prob_iter1 = global_prob_iter1

        # Calculate PEP from iteration 1 scores
        order = sortperm(features_refined.global_prob_iter1, rev=true)
        sorted_scores = features_refined.global_prob_iter1[order]
        sorted_targets = features_refined.target[order]

        PEPs = Vector{Float32}(undef, length(order))
        get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

        # Calculate global q-values from iteration 1 scores
        global_qvals = Vector{Float32}(undef, length(order))
        get_q_values!(sorted_scores, sorted_targets, global_qvals; doSort=false)

        # Map PEPs and q-values back to original order
        features_refined.PEP_iter1 = similar(PEPs)
        features_refined.global_qval_iter1 = similar(global_qvals)
        features_refined.PEP_iter1[order] = PEPs
        features_refined.global_qval_iter1[order] = global_qvals

        # Relabel uncertain targets (high PEP) as negatives
        uncertain_target_mask = features_refined.target .&&
                                (features_refined.PEP_iter1 .>= min_PEP_neg_threshold)
        features_refined.target[uncertain_target_mask] .= false

        @user_info "  Iteration 2: Relabeled $(sum(uncertain_target_mask)) uncertain targets as negatives"

        # Filter to keep:
        # - All negatives (original decoys + relabeled uncertain targets)
        # - High-confidence targets (q-value <= threshold)
        features_refined = subset(
            features_refined,
            [:target, :global_qval_iter1] => ByRow((t, q) ->
                (!t) || (t && q <= max_global_qval_threshold)
            )
        )

        @user_info "  Iteration 2: Training set size: $(nrow(features_refined)) precursors"
        @user_info "    - Targets: $(sum(features_refined.target))"
        @user_info "    - Negatives: $(sum(.!features_refined.target))"

        # Drop temporary columns
        select!(features_refined, Not([:global_prob_iter1, :PEP_iter1, :global_qval_iter1]))

        return features_refined
    end
end

#==========================================================
Phase 3: Model Training with CV Folds
==========================================================#

"""
    train_global_prob_models(precursor_features::DataFrame,
                            params,
                            n_iterations::Int,
                            max_global_qval_threshold::Float32,
                            min_PEP_neg_threshold::Float32) -> Tuple{Any, Any, Vector{Symbol}}

Train two LightGBM models (one per CV fold) to predict target vs decoy precursors using iterative refinement.

# Arguments
- `precursor_features`: DataFrame with precursor-level features
- `params`: ScoringSearchParameters with model hyperparameters
- `n_iterations`: Number of training iterations (default 2)
- `max_global_qval_threshold`: Max q-value for keeping targets in iteration 2
- `min_PEP_neg_threshold`: Min PEP for relabeling targets as negatives in iteration 2

# Returns
Tuple of (model_fold_0, model_fold_1, feature_cols)
"""
function train_global_prob_models(
    precursor_features::DataFrame,
    params,  # ScoringSearchParameters
    n_iterations::Int = 2,
    max_global_qval_threshold::Float32 = 0.01f0,
    min_PEP_neg_threshold::Float32 = 0.5f0
)::Tuple{Union{Nothing, Any}, Union{Nothing, Any}, Vector{Symbol}}
    # Explicitly construct feature columns (7 core features)
    feature_cols = Symbol[
        :logodds_score,
        :n_observations,
        :median_score,
        :std_score,
        :max_score,
        :skewness,
        :kurtosis
    ]

    @user_info "Training global probability models with $(n_iterations) iterations"
    @user_info "Features ($(length(feature_cols))): $(join(string.(feature_cols), ", "))"

    # Check CV folds
    unique_folds = unique(precursor_features.cv_fold)
    if length(unique_folds) != 2
        @user_warn "Expected 2 CV folds, got $(length(unique_folds)). Using logodds only."
        return (nothing, nothing, feature_cols)
    end

    # Storage for iteration 1 predictions (needed for iteration 2 filtering)
    global_prob_iter1 = nothing

    # Iterative training loop
    models = Vector{Any}(undef, 2)

    for itr in 1:n_iterations
        @user_info "Global Probability Training - Iteration $itr"

        # Get training data for this iteration
        features_itr = get_training_data_for_global_prob_iteration!(
            precursor_features,
            itr,
            global_prob_iter1,
            max_global_qval_threshold,
            min_PEP_neg_threshold
        )

        # Train models on current iteration's data
        for (fold_idx, fold) in enumerate(unique_folds)
            test_mask = features_itr.cv_fold .== fold
            train_mask = .!test_mask

            train_data = features_itr[train_mask, feature_cols]
            train_labels = features_itr.target[train_mask]

            # Train LightGBM classifier
            classifier = build_lightgbm_classifier(
                num_iterations = 200,
                max_depth = 10,
                learning_rate = 0.1,
                feature_fraction = 1.0,
                bagging_fraction = 0.5,
                bagging_freq = 1,
                min_data_in_leaf = 100,
                min_gain_to_split = 1.0
            )

            models[fold_idx] = fit_lightgbm_model(
                classifier, train_data, train_labels;
                positive_label=true
            )

            @user_info "  Fold $fold (iter $itr): trained on $(sum(train_mask)) precursors, test $(sum(test_mask))"
        end

        # Generate predictions for this iteration (on FULL dataset, not filtered)
        # This ensures we have predictions for all precursors for iteration 2 filtering
        if itr < n_iterations
            global_prob_iter1 = predict_global_probs_cv(
                precursor_features,  # Use FULL dataset for predictions
                models[1],
                models[2],
                feature_cols
            )
            @user_info "  Generated iteration $itr predictions for $(length(global_prob_iter1)) precursors"
        end
    end

    # Print feature importances after final iteration
    @user_info "Feature Importances (Final Iteration):"
    for (fold_idx, model) in enumerate(models)
        importances = lightgbm_feature_importances(model)
        if importances === nothing
            @user_warn "  Fold $fold_idx: No feature importances available"
        else
            feature_pairs = collect(zip(feature_cols, importances))
            sort!(feature_pairs, by=x->x[2], rev=true)
            @user_info "  Fold $fold_idx:"
            for (feat, score) in feature_pairs
                @user_info "    $feat: $(round(score, digits=3))"
            end
        end
    end

    return models[1], models[2], feature_cols
end

"""
    predict_global_probs_cv(precursor_features::DataFrame,
                           model_fold_0::Any,
                           model_fold_1::Any,
                           feature_cols::Vector{Symbol}) -> Vector{Float32}

Generate out-of-fold predictions using trained models.

# Arguments
- `precursor_features`: DataFrame with precursor-level features
- `model_fold_0`: Model trained on fold 0 data
- `model_fold_1`: Model trained on fold 1 data
- `feature_cols`: Feature column names

# Returns
Vector of predicted probabilities (out-of-fold)
"""
function predict_global_probs_cv(
    precursor_features::DataFrame,
    model_fold_0::Any,
    model_fold_1::Any,
    feature_cols::Vector{Symbol}
)::Vector{Float32}
    predictions = Vector{Float32}(undef, nrow(precursor_features))

    # Predict out-of-fold
    for row_idx in 1:nrow(precursor_features)
        fold = precursor_features.cv_fold[row_idx]
        model = fold == 0 ? model_fold_1 : model_fold_0

        # Extract features for this row
        features_df = precursor_features[row_idx:row_idx, feature_cols]
        pred = lightgbm_predict(model, features_df)
        predictions[row_idx] = pred[1]
    end

    return predictions
end

#==========================================================
Phase 4: Model Selection via AUC Comparison
==========================================================#

"""
    calculate_auc(scores::Vector{Float32}, labels::Vector{Bool}) -> Float64

Calculate area under ROC curve using trapezoidal rule.

# Arguments
- `scores`: Predicted scores (higher = more likely to be target)
- `labels`: True labels (true = target, false = decoy)

# Returns
AUC value between 0 and 1
"""
function calculate_auc(scores::Vector{Float32}, labels::Vector{Bool})::Float64
    # Sort by scores descending
    sorted_idx = sortperm(scores, rev=true)
    sorted_labels = labels[sorted_idx]

    n_pos = sum(sorted_labels)
    n_neg = length(sorted_labels) - n_pos

    if n_pos == 0 || n_neg == 0
        return 0.5  # Undefined AUC
    end

    # Calculate AUC using trapezoidal rule
    tpr = cumsum(sorted_labels) / n_pos
    fpr = cumsum(.!sorted_labels) / n_neg

    auc = 0.0
    for i in 2:length(fpr)
        auc += (fpr[i] - fpr[i-1]) * (tpr[i] + tpr[i-1]) / 2
    end

    return auc
end

"""
    select_best_global_prob_method(precursor_features::DataFrame,
                                   ml_predictions::Vector{Float32}) -> Symbol

Compare ML model vs logodds baseline using AUC on first CV fold.

# Arguments
- `precursor_features`: DataFrame with logodds_score and target columns
- `ml_predictions`: ML model predictions

# Returns
:ml if ML model has better AUC, :logodds otherwise
"""
function select_best_global_prob_method(
    precursor_features::DataFrame,
    ml_predictions::Vector{Float32}
)::Symbol
    # Use first CV fold for comparison (fold 0)
    fold_0_mask = precursor_features.cv_fold .== 0

    labels = precursor_features.target[fold_0_mask]
    logodds_scores = precursor_features.logodds_score[fold_0_mask]
    ml_scores = ml_predictions[fold_0_mask]

    auc_logodds = calculate_auc(logodds_scores, labels)
    auc_ml = calculate_auc(ml_scores, labels)

    @user_info "Global Probability Method Comparison (CV Fold 0):"
    @user_info "  Logodds AUC: $(round(auc_logodds, digits=4))"
    @user_info "  ML Model AUC: $(round(auc_ml, digits=4))"

    if auc_ml > auc_logodds
        @user_info "  ✓ Selected: ML Model (AUC improvement: +$(round(auc_ml - auc_logodds, digits=4)))"
        return :ml
    else
        @user_info "  ✓ Selected: Logodds (ML did not improve performance)"
        return :logodds
    end
end
