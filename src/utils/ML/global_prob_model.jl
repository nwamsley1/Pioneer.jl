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
        Series(Mean(), Variance(), Skewness()),
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
    fit!(acc.stats, prob)

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
    build_global_prec_features(merged_path::String, precursors, n_runs::Int;
                               top_n::Int=5, prob_thresh::Float32=0.95f0) ->
                               (Dict{UInt32, GlobalPrecFeatures{top_n}}, Dict{UInt32, Bool})

Build per-precursor feature summaries for global probability model from merged search results.
Uses single-pass streaming approach with accumulators to minimize memory and avoid groupby.

# Arguments
- `merged_path`: Path to Arrow file with columns :precursor_idx, :prec_prob
- `precursors`: Precursor library for target/decoy labels
- `n_runs`: Total number of MS files in experiment
- `top_n`: Number of top scores to retain (default 5)
- `prob_thresh`: Threshold for n_above_thresh count (default 0.95)

# Returns
- Dict mapping precursor_idx to GlobalPrecFeatures{top_n}
- Dict mapping precursor_idx to target label (Bool)
"""
function build_global_prec_features(
    merged_path::String,
    precursors,
    n_runs::Int;
    top_n::Int=5,
    prob_thresh::Float32=0.95f0
)
    sqrt_n_runs = floor(Int, sqrt(n_runs))

    # Load Arrow table (memory-mapped, no full DataFrame materialization)
    merged_table = Arrow.Table(merged_path)

    # Dictionary of streaming accumulators
    if top_n == 5
        accumulators = Dict{UInt32, PrecursorAccumulator{5}}()
    else
        accumulators = Dict{UInt32, PrecursorAccumulator{top_n}}()
    end

    # Single-pass streaming through Arrow table
    n_rows = length(merged_table.precursor_idx)
    for row_idx in 1:n_rows
        pid = UInt32(merged_table.precursor_idx[row_idx])
        prob = Float32(merged_table.prec_prob[row_idx])

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
            skewp = n_available >= 3 ? Float32(value(acc.stats[3])) : -1.0f0
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

        # Get target label from library
        labels[pid] = getIsDecoy(precursors[pid])
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
                           num_iterations::Int=200, learning_rate::Float64=0.15,
                           num_leaves::Int=63, max_depth::Int=10,
                           feature_fraction::Float64=0.5, bagging_fraction::Float64=0.5,
                           min_data_in_leaf::Int=1, min_gain_to_split::Float64=0.0) ->
                           (LightGBMModel, Vector{Float32})

Train LightGBM classifier for global precursor probability with CV folds.
Returns fitted model and out-of-fold predictions for diagnostics.

# Arguments
- `feat_df`: Feature DataFrame with :target column and feature columns
- `folds`: Vector of fold assignments (one per row in feat_df)
- Hyperparameters matching MBR scoring defaults

# Returns
- Fitted LightGBMModel on all data
- Vector of out-of-fold predictions (for CV diagnostics)
"""
function train_global_prob_model(
    feat_df::DataFrame,
    folds::Vector{Int};
    num_iterations::Int=200,
    learning_rate::Float64=0.15,
    num_leaves::Int=63,
    max_depth::Int=10,
    feature_fraction::Float64=0.5,
    bagging_fraction::Float64=0.5,
    min_data_in_leaf::Int=1,
    min_gain_to_split::Float64=0.0
)
    # Separate features and labels
    feature_cols = setdiff(names(feat_df), [:precursor_idx, :target])
    X = feat_df[!, feature_cols]
    y = feat_df.target

    unique_folds = sort(unique(folds))

    # Out-of-fold predictions for diagnostics
    oof_preds = Vector{Float32}(undef, nrow(feat_df))
    fill!(oof_preds, -1.0f0)

    # Cross-validation
    for test_fold in unique_folds
        test_mask = folds .== test_fold
        train_mask = .!test_mask

        X_train = X[train_mask, :]
        y_train = y[train_mask]
        X_test = X[test_mask, :]

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

        # Predict on held-out fold
        oof_preds[test_mask] = lightgbm_predict(model, X_test; output_type=Float32)
    end

    # Train final model on all data
    lgbm_final = build_lightgbm_classifier(;
        num_iterations=num_iterations,
        learning_rate=learning_rate,
        num_leaves=num_leaves,
        max_depth=max_depth,
        feature_fraction=feature_fraction,
        bagging_fraction=bagging_fraction,
        bagging_freq=1,
        min_data_in_leaf=min_data_in_leaf,
        min_gain_to_split=min_gain_to_split,
        objective="binary",
        metric=["binary_logloss", "auc"],
        verbosity=-1
    )

    final_model = fit_lightgbm_model(lgbm_final, X, y)

    return final_model, oof_preds
end

"""
    predict_global_prob(model::LightGBMModel, feat_df::DataFrame) -> Dict{UInt32, Float32}

Generate per-precursor probability predictions from trained model.
"""
function predict_global_prob(model::LightGBMModel, feat_df::DataFrame)
    feature_cols = setdiff(names(feat_df), [:precursor_idx, :target])
    X = feat_df[!, feature_cols]

    preds = lightgbm_predict(model, X; output_type=Float32)

    prob_map = Dict{UInt32, Float32}()
    for (i, pid) in enumerate(feat_df.precursor_idx)
        prob_map[pid] = preds[i]
    end

    return prob_map
end

"""
    compare_global_prob_methods(model_scores::Dict{UInt32, Float32},
                               baseline_scores::Dict{UInt32, Float32},
                               merged_path::String,
                               params,
                               search_context) ->
                               (Bool, Int, Int)

Compare model-based and baseline global probability methods by computing q-values
and counting passing precursors. Returns (use_model, model_passing, baseline_passing).

Uses existing Pioneer q-value pipeline: assigns scores to :global_prob, writes to
temporary file, runs stream_sorted_merge + get_precursor_global_qval_spline.
"""
function compare_global_prob_methods(
    model_scores::Dict{UInt32, Float32},
    baseline_scores::Dict{UInt32, Float32},
    merged_path::String,
    params,
    search_context
)
    # Load original merged data
    merged_table = Arrow.Table(merged_path)

    temp_dir = mktempdir()

    try
        # Test model scores
        model_df = DataFrame(merged_table)
        model_df.global_prob = [model_scores[UInt32(pid)] for pid in model_df.precursor_idx]
        model_path = joinpath(temp_dir, "model_merged.arrow")
        Arrow.write(model_path, model_df)

        # Sort and compute q-values for model
        model_sorted_path = joinpath(temp_dir, "model_sorted.arrow")
        model_ref = FileReference(model_path, :arrow)
        stream_sorted_merge([model_ref], model_sorted_path, :global_prob, :target;
                          batch_size=10_000_000, reverse=[true, true])

        model_qval_interp = get_precursor_global_qval_spline(
            model_sorted_path, params, search_context
        )

        # Count passing precursors for model
        model_sorted = DataFrame(Arrow.Table(model_sorted_path))
        model_qvals = [model_qval_interp(p) for p in model_sorted.global_prob]
        model_passing = count(<=(params.precursor_global_qvalue_threshold), model_qvals)

        # Test baseline scores
        baseline_df = DataFrame(merged_table)
        baseline_df.global_prob = [baseline_scores[UInt32(pid)] for pid in baseline_df.precursor_idx]
        baseline_path = joinpath(temp_dir, "baseline_merged.arrow")
        Arrow.write(baseline_path, baseline_df)

        # Sort and compute q-values for baseline
        baseline_sorted_path = joinpath(temp_dir, "baseline_sorted.arrow")
        baseline_ref = FileReference(baseline_path, :arrow)
        stream_sorted_merge([baseline_ref], baseline_sorted_path, :global_prob, :target;
                          batch_size=10_000_000, reverse=[true, true])

        baseline_qval_interp = get_precursor_global_qval_spline(
            baseline_sorted_path, params, search_context
        )

        # Count passing precursors for baseline
        baseline_sorted = DataFrame(Arrow.Table(baseline_sorted_path))
        baseline_qvals = [baseline_qval_interp(p) for p in baseline_sorted.global_prob]
        baseline_passing = count(<=(params.precursor_global_qvalue_threshold), baseline_qvals)

        # Choose method with more passing precursors
        use_model = model_passing >= baseline_passing

        return use_model, model_passing, baseline_passing

    finally
        # Clean up temp directory
        rm(temp_dir; recursive=true, force=true)
    end
end
