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
Streaming accumulator for building GlobalPrecFeatures per precursor.
Maintains bounded top-N scores buffer and OnlineStats for statistics.
"""
mutable struct PrecursorAccumulator{N}
    stats::Series  # OnlineStats: Mean, Variance, Skewness
    top_scores::Vector{Float32}  # Bounded buffer for top-N (sorted descending)
    n_above_thresh::Int32
    max_observed::Float32  # Track max for efficiency
end

function PrecursorAccumulator{N}() where N
    PrecursorAccumulator{N}(
        Series(Mean(), Variance(), Moments()),
        Float32[],
        Int32(0),
        -Inf32
    )
end

"""
    update_accumulator!(acc::PrecursorAccumulator{N}, prob::Float32,
                       prob_thresh::Float32) where N

Update accumulator with a new probability observation.
Maintains top-N scores using bounded insertion sort.
"""
function update_accumulator!(acc::PrecursorAccumulator{N}, prob::Float32, prob_thresh::Float32) where N
    # Update streaming statistics
    OnlineStats.fit!(acc.stats, prob)

    # Track max
    if prob > acc.max_observed
        acc.max_observed = prob
    end

    # Maintain top-N scores using bounded insertion
    if length(acc.top_scores) < N
        push!(acc.top_scores, prob)
        sort!(acc.top_scores; rev=true)
    elseif prob > acc.top_scores[end]
        acc.top_scores[end] = prob
        # Use insertion sort for small N (typically 5)
        for i in (N-1):-1:1
            if acc.top_scores[i+1] > acc.top_scores[i]
                acc.top_scores[i], acc.top_scores[i+1] = acc.top_scores[i+1], acc.top_scores[i]
            else
                break
            end
        end
    end

    # Count above threshold
    if prob >= prob_thresh
        acc.n_above_thresh += 1
    end

    return nothing
end

"""
    build_global_prec_features(merged_df::AbstractDataFrame, precursors, n_runs::Int;
                               top_n::Int=5, prob_thresh::Float32=0.95f0,
                               prob_column::Symbol=:prec_prob) ->
                               (Dict{UInt32, GlobalPrecFeatures{top_n}}, Dict{UInt32, Bool})

Build per-precursor feature summaries for global probability model from merged search results.
Uses single-pass streaming approach with accumulators to minimize memory and avoid groupby.

# Arguments
- `merged_df`: DataFrame with columns :precursor_idx and prob_column
- `precursors`: Precursor library for target/decoy labels
- `n_runs`: Total number of MS files in experiment
- `top_n`: Number of top scores to retain (default 5)
- `prob_thresh`: Threshold for n_above_thresh count (default 0.95)
- `prob_column`: Column to use for probability values (default :prec_prob, use :infold_prec_prob to avoid data leakage)

# Returns
- Dict mapping precursor_idx to GlobalPrecFeatures{top_n}
- Dict mapping precursor_idx to target label (Bool)
"""
function build_global_prec_features(
    merged_df::AbstractDataFrame,
    precursors,
    n_runs::Int;
    top_n::Int=5,
    prob_thresh::Float32=0.95f0,
    prob_column::Symbol=:prec_prob
)
    sqrt_n_runs = floor(Int, sqrt(n_runs))

    # Dictionary of streaming accumulators
    if top_n == 5
        accumulators = Dict{UInt32, PrecursorAccumulator{5}}()
    else
        accumulators = Dict{UInt32, PrecursorAccumulator{top_n}}()
    end

    # Single-pass streaming through DataFrame
    n_rows = nrow(merged_df)
    for row_idx in 1:n_rows
        pid = UInt32(merged_df.precursor_idx[row_idx])
        prob = Float32(merged_df[row_idx, prob_column])

        # Get or create accumulator
        if !haskey(accumulators, pid)
            if top_n == 5
                accumulators[pid] = PrecursorAccumulator{5}()
            else
                accumulators[pid] = PrecursorAccumulator{top_n}()
            end
        end

        update_accumulator!(accumulators[pid], prob, prob_thresh)
    end

    # Finalize: convert accumulators to GlobalPrecFeatures
    dict = Dict{UInt32, GlobalPrecFeatures{top_n}}()
    labels = Dict{UInt32, Bool}()

    # Get all is_decoy flags once (returns vector)
    is_decoy_array = getIsDecoy(precursors)

    for (pid, acc) in accumulators
        n_available = nobs(acc.stats[1])  # Observation count from Mean()

        # Pad top scores to fixed size
        top_scores = ntuple(i -> i <= length(acc.top_scores) ? acc.top_scores[i] : -1.0f0, top_n)

        # Extract statistics
        if n_available > 0
            meanp = Float32(value(acc.stats[1]))
            maxp = acc.max_observed
            # Min not tracked in streaming (would need Extrema() which doesn't support skewness in same Series)
            minp = -1.0f0
            stdp = n_available > 1 ? Float32(sqrt(value(acc.stats[2]))) : -1.0f0
            skewp = if n_available >= 3
                try
                    Float32(skewness(acc.stats[3]))
                catch
                    -1.0f0  # Return sentinel if skewness fails (near-zero variance)
                end
            else
                -1.0f0
            end
        else
            meanp = maxp = minp = stdp = skewp = -1.0f0
        end

        # Delta features
        n_top = length(acc.top_scores)
        del12 = (n_top >= 2) ? (top_scores[1] - top_scores[2]) : -1.0f0
        del23 = (n_top >= 3) ? (top_scores[2] - top_scores[3]) : -1.0f0

        # Baseline: logodds on top scores (approximation, but avoids storing all scores)
        lodds_baseline = if n_top > 0
            Float32(logodds(acc.top_scores, sqrt_n_runs))
        else
            -1.0f0
        end

        dict[pid] = GlobalPrecFeatures{top_n}(
            top_scores,
            Int32(n_available),
            (meanp, maxp, minp, stdp, skewp),
            (Int32(n_available), Int32(n_runs), acc.n_above_thresh),
            (del12, del23),
            lodds_baseline
        )

        # Get target label from library (true=target, false=decoy)
        labels[pid] = !is_decoy_array[pid]
    end

    return dict, labels
end

"""
    features_to_dataframe(features_dict::Dict{UInt32, GlobalPrecFeatures{N}},
                         labels_dict::Dict{UInt32, Bool}) -> DataFrame where N

Convert GlobalPrecFeatures dictionary to a DataFrame suitable for model training.
Each feature field is expanded into separate columns.
"""
function features_to_dataframe(
    features_dict::Dict{UInt32, GlobalPrecFeatures{N}},
    labels_dict::Dict{UInt32, Bool}
) where N
    n_prec = length(features_dict)

    # Pre-allocate vectors
    precursor_idx = Vector{UInt32}(undef, n_prec)
    target = Vector{Bool}(undef, n_prec)

    # Top scores
    top_scores = [Vector{Float32}(undef, n_prec) for _ in 1:N]

    # Stats
    mean_score = Vector{Float32}(undef, n_prec)
    max_score = Vector{Float32}(undef, n_prec)
    min_score = Vector{Float32}(undef, n_prec)
    std_score = Vector{Float32}(undef, n_prec)
    skew_score = Vector{Float32}(undef, n_prec)

    # Counts
    n_scores_available = Vector{Int32}(undef, n_prec)
    n_runs_with_score = Vector{Int32}(undef, n_prec)
    n_runs_total = Vector{Int32}(undef, n_prec)
    n_above_thresh = Vector{Int32}(undef, n_prec)

    # Deltas
    delta_12 = Vector{Float32}(undef, n_prec)
    delta_23 = Vector{Float32}(undef, n_prec)

    # Baseline
    logodds_baseline = Vector{Float32}(undef, n_prec)

    idx = 1
    for (pid, feat) in features_dict
        precursor_idx[idx] = pid
        target[idx] = labels_dict[pid]

        # Unpack top scores
        for i in 1:N
            top_scores[i][idx] = feat.top_scores[i]
        end

        n_scores_available[idx] = feat.n_scores_available

        # Unpack stats
        mean_score[idx] = feat.stats[1]
        max_score[idx] = feat.stats[2]
        min_score[idx] = feat.stats[3]
        std_score[idx] = feat.stats[4]
        skew_score[idx] = feat.stats[5]

        # Unpack counts
        n_runs_with_score[idx] = feat.counts[1]
        n_runs_total[idx] = feat.counts[2]
        n_above_thresh[idx] = feat.counts[3]

        # Unpack deltas
        delta_12[idx] = feat.deltas[1]
        delta_23[idx] = feat.deltas[2]

        logodds_baseline[idx] = feat.logodds_baseline

        idx += 1
    end

    # Build DataFrame
    df = DataFrame(
        precursor_idx = precursor_idx,
        target = target
    )

    # Add top scores
    for i in 1:N
        df[!, Symbol("top_", i)] = top_scores[i]
    end

    df[!, :n_scores_available] = n_scores_available
    df[!, :mean_score] = mean_score
    df[!, :max_score] = max_score
    df[!, :min_score] = min_score
    df[!, :std_score] = std_score
    df[!, :skew_score] = skew_score
    df[!, :n_runs_with_score] = n_runs_with_score
    df[!, :n_runs_total] = n_runs_total
    df[!, :n_above_thresh] = n_above_thresh
    df[!, :delta_12] = delta_12
    df[!, :delta_23] = delta_23
    df[!, :logodds_baseline] = logodds_baseline

    return df
end

"""
    train_global_prob_model(feat_df::DataFrame, folds::Vector{Int};
                           num_iterations::Int=100, learning_rate::Float64=0.1,
                           num_leaves::Int=8, max_depth::Int=3,
                           feature_fraction::Float64=0.8, bagging_fraction::Float64=0.5,
                           min_data_in_leaf::Int=1, min_gain_to_split::Float64=1.0) ->
                           Vector{Float32}

Train LightGBM classifier for global precursor probability with CV folds.
Returns out-of-fold predictions without data leakage.

# Arguments
- `feat_df`: Feature DataFrame with :target column and feature columns
- `folds`: Vector of fold assignments (one per row in feat_df)
- Hyperparameters matching MBR model from scoring_interface.jl

# Returns
- Vector of out-of-fold predictions (one per row in feat_df)
"""
function train_global_prob_model(
    feat_df::DataFrame,
    folds::Vector{UInt8};
    num_iterations::Int=100,
    learning_rate::Float64=0.1,
    num_leaves::Int=8,
    max_depth::Int=3,
    feature_fraction::Float64=0.8,
    bagging_fraction::Float64=0.5,
    min_data_in_leaf::Int=1,
    min_gain_to_split::Float64=1.0
)
    # Diagnostic logging for feature distributions
    @info "Global prob model diagnostics:" n_precursors=nrow(feat_df) n_folds=length(unique(folds))

    # Log target/decoy balance
    n_targets = count(feat_df.target)
    n_decoys = nrow(feat_df) - n_targets
    @info "Target/Decoy balance:" n_targets n_decoys ratio=n_targets/n_decoys

    # Log feature ranges for targets vs decoys
    target_mask = feat_df.target
    if :top_1 in propertynames(feat_df)
        target_top1_range = (minimum(feat_df.top_1[target_mask]), maximum(feat_df.top_1[target_mask]))
        decoy_top1_range = (minimum(feat_df.top_1[.!target_mask]), maximum(feat_df.top_1[.!target_mask]))
        @info "Top-1 score ranges:" targets=target_top1_range decoys=decoy_top1_range
    end

    if :logodds_baseline in propertynames(feat_df)
        target_baseline_range = (minimum(feat_df.logodds_baseline[target_mask]), maximum(feat_df.logodds_baseline[target_mask]))
        decoy_baseline_range = (minimum(feat_df.logodds_baseline[.!target_mask]), maximum(feat_df.logodds_baseline[.!target_mask]))
        @info "Baseline logodds ranges:" targets=target_baseline_range decoys=decoy_baseline_range
    end

    # Separate features and labels
    # Note: names() returns Strings, so use String literals not Symbols
    feature_cols = setdiff(names(feat_df), ["precursor_idx", "target"])
    X = feat_df[!, feature_cols]
    y = feat_df.target

    unique_folds = sort(unique(folds))

    # Out-of-fold predictions for diagnostics
    oof_preds = Vector{Float32}(undef, nrow(feat_df))
    fill!(oof_preds, -1.0f0)

    # Track per-fold AUCs and feature importances
    fold_aucs = Float32[]
    all_importances = Dict{Symbol, Vector{Float32}}()

    # Cross-validation
    for test_fold in unique_folds
        test_mask = folds .== test_fold
        train_mask = .!test_mask

        X_train = X[train_mask, :]
        y_train = y[train_mask]
        X_test = X[test_mask, :]
        y_test = y[test_mask]

        # Build and fit model
        lgbm = build_lightgbm_classifier(;
            num_iterations=num_iterations,
            learning_rate=learning_rate,
            num_leaves=num_leaves,
            max_depth=max_depth,
            feature_fraction=feature_fraction,
            bagging_fraction=bagging_fraction,
            bagging_freq=1,  # Enable bagging
            min_data_in_leaf=min_data_in_leaf,
            min_gain_to_split=min_gain_to_split,
            objective="binary",
            metric=["binary_logloss", "auc"],
            verbosity=-1
        )

        model = fit_lightgbm_model(lgbm, X_train, y_train)

        # Extract feature importances
        importances = lightgbm_feature_importances(model)
        if importances !== nothing
            for (i, feat) in enumerate(feature_cols)
                if !haskey(all_importances, Symbol(feat))
                    all_importances[Symbol(feat)] = Float32[]
                end
                push!(all_importances[Symbol(feat)], importances[i])
            end
        end

        # Predict on held-out fold
        fold_preds = lightgbm_predict(model, X_test; output_type=Float32)
        oof_preds[test_mask] = fold_preds

        # Diagnostic: check prediction distribution
        n_unique_preds = length(unique(fold_preds))
        pred_range = (minimum(fold_preds), maximum(fold_preds))
        @info "Fold $test_fold predictions:" n_unique=n_unique_preds range=pred_range

        # Calculate per-fold AUC
        fold_auc = calculate_auc(fold_preds, y_test)
        push!(fold_aucs, fold_auc)
        @info "Fold $test_fold AUC:" fold_auc n_test=sum(test_mask) n_train=sum(train_mask)
    end

    # Log overall fold AUC statistics
    @info "Per-fold AUC statistics:" mean=mean(fold_aucs) min=minimum(fold_aucs) max=maximum(fold_aucs)

    # Log feature importances (averaged across folds)
    if !isempty(all_importances)
        avg_importances = [(feat, mean(imps)) for (feat, imps) in all_importances]
        sort!(avg_importances, by=x->x[2], rev=true)
        @info "Top 10 feature importances (averaged across folds):"
        for (i, (feat, imp)) in enumerate(avg_importances[1:min(10, length(avg_importances))])
            @info "  $i. $feat: $(round(imp, digits=2))"
        end
    end

    return oof_preds
end

"""
    predict_global_prob(model::LightGBMModel, feat_df::DataFrame) -> Dict{UInt32, Float32}

Generate per-precursor probability predictions from trained model.
"""
function predict_global_prob(model::LightGBMModel, feat_df::DataFrame)
    # Note: names() returns Strings, so use String literals not Symbols
    feature_cols = setdiff(names(feat_df), ["precursor_idx", "target"])
    X = feat_df[!, feature_cols]

    preds = lightgbm_predict(model, X; output_type=Float32)

    prob_map = Dict{UInt32, Float32}()
    for (i, pid) in enumerate(feat_df.precursor_idx)
        prob_map[pid] = preds[i]
    end

    return prob_map
end

"""
    compare_global_prob_methods(model::LightGBMModel,
                               feat_df::DataFrame) ->
                               (Bool, Float32, Float32)

Compare model-based and baseline global probability methods using ROC AUC.
Compares model predictions against the logodds_baseline column in feat_df.
Returns (use_model, model_auc, baseline_auc) where use_model indicates whether
the ML model achieved better AUC than the baseline method.
"""
function compare_global_prob_methods(
    model::LightGBMModel,
    feat_df::DataFrame
)
    # Get model predictions (use same features as training)
    # Note: names() returns Strings, so use String literals not Symbols
    feature_cols = setdiff(names(feat_df), ["precursor_idx", "target"])
    X = feat_df[!, feature_cols]
    model_preds = lightgbm_predict(model, X; output_type=Float32)

    # Extract baseline scores and target flags
    baseline_preds = Vector{Float32}(feat_df.logodds_baseline)
    is_target = Vector{Bool}(feat_df.target)

    # Calculate AUC for both methods
    model_auc = calculate_auc(model_preds, is_target)
    baseline_auc = calculate_auc(baseline_preds, is_target)

    # Choose method with better AUC
    use_model = model_auc >= baseline_auc

    return use_model, model_auc, baseline_auc
end

"""
    compare_global_prob_methods_oof(oof_preds::Vector{Float32},
                                   feat_df::DataFrame) ->
                                   (Bool, Float32, Float32)

Compare out-of-fold model predictions against baseline global probability using ROC AUC.
This version uses pre-computed OOF predictions to avoid data leakage.
Returns (use_model, model_auc, baseline_auc) where use_model indicates whether
the ML model achieved better AUC than the baseline method.
"""
function compare_global_prob_methods_oof(
    oof_preds::Vector{Float32},
    feat_df::DataFrame
)
    # Extract baseline scores and target flags
    baseline_preds = Vector{Float32}(feat_df.logodds_baseline)
    is_target = Vector{Bool}(feat_df.target)

    # Calculate AUC for both methods with verbose diagnostics for model
    @info "Calculating model AUC with diagnostics..."
    model_auc = calculate_auc(oof_preds, is_target; verbose=true)

    @info "Calculating baseline AUC..."
    baseline_auc = calculate_auc(baseline_preds, is_target; verbose=false)

    # Choose method with better AUC
    use_model = model_auc >= baseline_auc

    return use_model, model_auc, baseline_auc
end

"""
    calculate_auc(scores::Vector{Float32}, is_target::Vector{Bool}; verbose::Bool=false)

Calculate ROC AUC using Mann-Whitney U statistic.
Higher scores should indicate higher probability of being a target.
"""
function calculate_auc(scores::Vector{Float32}, is_target::Vector{Bool}; verbose::Bool=false)
    n = length(scores)
    n == length(is_target) || throw(ArgumentError("scores and is_target must have same length"))

    # Diagnostic: score distributions
    if verbose
        target_scores = scores[is_target]
        decoy_scores = scores[.!is_target]
        @info "Score distributions:" target_range=(minimum(target_scores), maximum(target_scores)) decoy_range=(minimum(decoy_scores), maximum(decoy_scores))
        @info "Score quantiles (targets):" q25=quantile(target_scores, 0.25) q50=quantile(target_scores, 0.5) q75=quantile(target_scores, 0.75)
        @info "Score quantiles (decoys):" q25=quantile(decoy_scores, 0.25) q50=quantile(decoy_scores, 0.5) q75=quantile(decoy_scores, 0.75)

        # Check if all targets > all decoys (which would give AUC = 1.0)
        min_target = minimum(target_scores)
        max_decoy = maximum(decoy_scores)
        if min_target > max_decoy
            @warn "Perfect separation detected: minimum target score ($min_target) > maximum decoy score ($max_decoy)"
        end
    end

    # Sort by score (descending) and get ranks
    perm = sortperm(scores, rev=true)
    sorted_is_target = is_target[perm]

    # Calculate sum of ranks for targets (1-indexed ranks)
    n_targets = count(sorted_is_target)
    n_decoys = n - n_targets

    # Edge cases
    n_targets == 0 && return 0.0f0
    n_decoys == 0 && return 1.0f0

    # Sum of ranks for targets (ranks are 1-indexed)
    rank_sum = sum(i for (i, is_tgt) in enumerate(sorted_is_target) if is_tgt)

    # Mann-Whitney U statistic (adjusted for descending sort)
    # Since we sort descending (high scores first), we need the complement
    # U = number of (target, decoy) pairs where target score > decoy score
    u_standard = rank_sum - n_targets * (n_targets + 1) / 2
    u = n_targets * n_decoys - u_standard

    # AUC = U / (n_targets * n_decoys)
    auc = Float32(u / (n_targets * n_decoys))

    if verbose
        # Count correctly ordered pairs
        n_correct_pairs = Int(u)
        n_total_pairs = n_targets * n_decoys
        @info "Pair ordering:" n_correct=n_correct_pairs n_total=n_total_pairs pct_correct=100*n_correct_pairs/n_total_pairs
    end

    return auc
end
