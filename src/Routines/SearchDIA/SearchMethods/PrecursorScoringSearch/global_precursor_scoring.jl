# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pioneer.jl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

const GLOBAL_PRECURSOR_SCORE_FEATURES = Symbol[
    :empirical_global_score,
    :top1_prec_prob,
    :top2_prec_prob,
    :top2_logodds_score,
    :top3_logodds_score,
    :mean_prec_prob,
    :n_runs_observed,
    :n_prob_gt_0_9,
    :n_prob_gt_0_99,
]

const GLOBAL_PRECURSOR_SMALL_SCORE_FEATURES = Symbol[
    :empirical_global_score,
    :top2_logodds_score,
    :top3_logodds_score,
]

const GLOBAL_PRECURSOR_SMALL_DATASET_MAX_CANDIDATES = 300_000

const GLOBAL_PRECURSOR_LGBM_HP = (
    num_iterations = 100,
    learning_rate = 0.05,
    max_depth = 4,
    num_leaves = 15,
    min_data_in_leaf = 100,
    min_gain_to_split = 0.0,
    feature_fraction = 1.0,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    is_unbalance = true,
    max_bin = 255,
    lambda_l1 = 1.0,
    lambda_l2 = 1.0,
)

const GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT = 100
const GLOBAL_PRECURSOR_MAX_TRAIN = 1_000_000
const GLOBAL_PRECURSOR_CV_FOLDS = (UInt8(0), UInt8(1))

struct GlobalPrecursorInputs
    probabilities::Dict{UInt32, Vector{Float32}}
    targets::Dict{UInt32, Bool}
    folds::Dict{UInt32, UInt8}
end

function _select_global_precursor_score_features(n_candidates::Int)
    if n_candidates < GLOBAL_PRECURSOR_SMALL_DATASET_MAX_CANDIDATES
        return GLOBAL_PRECURSOR_SMALL_SCORE_FEATURES
    end
    return GLOBAL_PRECURSOR_SCORE_FEATURES
end

function _logodds_from_sorted(
    sorted_probabilities::AbstractVector{T},
    top_n::Int,
) where {T<:AbstractFloat}
    n_selected = min(length(sorted_probabilities), top_n)
    selected = @view sorted_probabilities[1:n_selected]
    probability_floor = 0.1f0
    probability_ceiling = 1.0f0 - 1.0f-6
    logits = log.(
        clamp.(selected, probability_floor, probability_ceiling) ./
        (1.0f0 .- clamp.(selected, probability_floor, probability_ceiling)),
    )
    return 1.0f0 / (1.0f0 + exp(-sum(logits) / n_selected))
end

"""
    logodds(probabilities, top_n)

Combine the highest `top_n` probabilities by averaging their log odds.
"""
function logodds(
    probabilities::AbstractVector{T},
    top_n::Int,
) where {T<:AbstractFloat}
    return _logodds_from_sorted(sort(probabilities; rev = true), top_n)
end

function _collect_global_precursor_inputs(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
)
    probabilities = Dict{UInt32, Vector{Float32}}()
    targets = Dict{UInt32, Bool}()
    folds = Dict{UInt32, UInt8}()
    sizehint!(probabilities, n_precursors)
    sizehint!(targets, n_precursors)
    sizehint!(folds, n_precursors)

    for ref in refs
        columns = materialize_columns(
            ref,
            [:precursor_idx, :prec_prob, :target, :cv_fold],
        )
        @inbounds for row in axes(columns, 1)
            precursor_idx = UInt32(columns.precursor_idx[row])
            precursor_probabilities = get!(probabilities, precursor_idx) do
                Float32[]
            end
            push!(precursor_probabilities, Float32(columns.prec_prob[row]))
            targets[precursor_idx] = Bool(columns.target[row])
            folds[precursor_idx] = UInt8(columns.cv_fold[row])
        end
    end

    isempty(probabilities) &&
        throw(ArgumentError("Global precursor scoring requires observed precursors"))
    return GlobalPrecursorInputs(probabilities, targets, folds)
end

function _build_global_precursor_feature_table(
    inputs::GlobalPrecursorInputs,
    n_runs_total::Int,
)
    precursor_ids = sort!(collect(keys(inputs.probabilities)))
    n_precursors = length(precursor_ids)
    top_run_count = isqrt(n_runs_total)
    feature_columns = Dict(
        feature => Vector{Float32}(undef, n_precursors)
        for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
    )

    @inbounds for (row, precursor_idx) in enumerate(precursor_ids)
        probabilities = inputs.probabilities[precursor_idx]
        sorted_probabilities = sort(probabilities; rev = true)
        n_observed = length(sorted_probabilities)
        top1 = sorted_probabilities[1]
        top2 = n_observed >= 2 ? sorted_probabilities[2] : 0.0f0

        feature_columns[:empirical_global_score][row] =
            _logodds_from_sorted(sorted_probabilities, top_run_count)
        feature_columns[:top1_prec_prob][row] = top1
        feature_columns[:top2_prec_prob][row] = top2
        feature_columns[:top2_logodds_score][row] =
            _logodds_from_sorted(sorted_probabilities, 2)
        feature_columns[:top3_logodds_score][row] =
            _logodds_from_sorted(sorted_probabilities, 3)
        feature_columns[:mean_prec_prob][row] = Float32(mean(probabilities))
        feature_columns[:n_runs_observed][row] = Float32(n_observed)
        feature_columns[:n_prob_gt_0_9][row] = Float32(count(>(0.9f0), probabilities))
        feature_columns[:n_prob_gt_0_99][row] = Float32(count(>(0.99f0), probabilities))
    end

    table = DataFrame(
        precursor_idx = precursor_ids,
        target = Bool[inputs.targets[precursor_idx] for precursor_idx in precursor_ids],
        cv_fold = UInt8[inputs.folds[precursor_idx] for precursor_idx in precursor_ids],
    )
    for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
        table[!, feature] = feature_columns[feature]
    end
    return table
end

function _log_global_precursor_feature_importances(
    models::Vector{LightGBMModel},
    features::Vector{Symbol},
)
    total_gains = Dict(feature => 0.0 for feature in features)
    for model in models
        for (feature, gain) in importance(model)
            total_gains[feature] += Float64(gain)
        end
    end

    sorted_gains = sort!(collect(total_gains); by = pair -> -last(pair))
    lines = ["Global precursor LightGBM feature gains (summed across folds):"]
    for (feature, gain) in sorted_gains
        push!(lines, "    $(rpad(string(feature), 52)) $(round(Int, gain))")
    end
    @debug_l1 join(lines, "\n")
    return nothing
end

function _score_global_precursor_features_oof(
    table::DataFrame,
    features::Vector{Symbol};
    lgbm_hp::NamedTuple = GLOBAL_PRECURSOR_LGBM_HP,
    min_training_class_count::Int = GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT,
    max_train::Int = GLOBAL_PRECURSOR_MAX_TRAIN,
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    observed_folds = sort!(unique(table.cv_fold))
    observed_folds == collect(GLOBAL_PRECURSOR_CV_FOLDS) ||
        throw(ArgumentError("Global precursor scoring requires CV folds 0 and 1"))

    targets = Vector{Bool}(table.target)
    training_mask = nothing
    previous_target_q01 = -1
    best_state = nothing

    for iteration in 1:max_iterations
        scores = Vector{Float32}(undef, nrow(table))
        models = LightGBMModel[]
        n_train_targets = 0
        n_train_decoys = 0
        iteration_valid = true

        for test_fold in GLOBAL_PRECURSOR_CV_FOLDS
            test_indices = findall(==(test_fold), table.cv_fold)
            train_indices = Int[
                row for row in axes(table, 1)
                if table.cv_fold[row] != test_fold &&
                   (training_mask === nothing || training_mask[row])
            ]

            if length(train_indices) > max_train
                rng = Random.MersenneTwister(1776 + 1000 * iteration + Int(test_fold))
                Random.shuffle!(rng, train_indices)
                resize!(train_indices, max_train)
            end

            train_targets = targets[train_indices]
            fold_target_count = count(train_targets)
            fold_decoy_count = length(train_targets) - fold_target_count
            if min(fold_target_count, fold_decoy_count) < min_training_class_count
                if iteration == 1
                    throw(ArgumentError(
                        "Global precursor scoring requires at least " *
                        "$min_training_class_count targets and decoys in each training fold",
                    ))
                end
                iteration_valid = false
                break
            end

            classifier = build_lightgbm_classifier(; lgbm_hp...)
            model = fit_lightgbm_model(
                classifier,
                view(table, train_indices, features),
                train_targets;
                positive_label = true,
            )
            scores[test_indices] .= lightgbm_predict(
                model,
                view(table, test_indices, features);
                output_type = Float32,
            )
            n_train_targets += fold_target_count
            n_train_decoys += fold_decoy_count
            push!(models, model)
        end

        if !iteration_valid
            @debug_l1 "Global precursor semi-supervised stopping: iteration " *
                       "$iteration does not retain enough targets and decoys; " *
                       "using iteration $(best_state.iter)"
            break
        end

        scores .= Float32.(clamp.(scores, 1.0f-6, 1.0f0 - 1.0f-4))
        metrics = _scoring_semisupervised_metrics_and_mask(
            scores,
            targets;
            train_q_threshold = train_q_threshold,
            stop_q_threshold = stop_q_threshold,
        )
        current_state = (
            scores = scores,
            models = models,
            target_q01 = metrics.target_q01,
            decoy_q01 = metrics.decoy_q01,
            iter = iteration,
        )
        @debug_l1 "Global precursor semi-supervised iteration $iteration: " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "q≤.01 targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        if iteration > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01,
            metrics.target_q01;
            min_fraction = min_gain,
        )
            best_state = _scoring_better_iteration_state(best_state, current_state)
            @debug_l1 "Global precursor semi-supervised stopping: iteration " *
                       "$iteration did not improve q≤.01 targets by " *
                       "$(round(100 * min_gain, digits = 2))%; using iteration " *
                       "$(best_state.iter)"
            break
        end

        best_state = _scoring_better_iteration_state(best_state, current_state)
        iteration == max_iterations && break
        previous_target_q01 = metrics.target_q01
        training_mask = metrics.training_mask
    end

    _log_global_precursor_feature_importances(best_state.models, features)
    return best_state
end

"""
    build_global_precursor_score_dicts(refs, n_precursors, n_runs_total)

Build one score-distribution feature row per precursor and return out-of-fold
global LightGBM scores with the corresponding target labels. For a single-run
search, return the individual precursor probabilities without training a model.
"""
function _build_single_run_precursor_score_dicts(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
)
    score_dict = Dict{UInt32, Float32}()
    target_dict = Dict{UInt32, Bool}()
    sizehint!(score_dict, n_precursors)
    sizehint!(target_dict, n_precursors)

    for ref in refs
        columns = materialize_columns(ref, [:precursor_idx, :prec_prob, :target])
        @inbounds for row in axes(columns, 1)
            precursor_idx = UInt32(columns.precursor_idx[row])
            score_dict[precursor_idx] = Float32(columns.prec_prob[row])
            target_dict[precursor_idx] = Bool(columns.target[row])
        end
    end
    return score_dict, target_dict
end

function build_global_precursor_score_dicts(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
    n_runs_total::Int,
)
    n_runs_total == 1 &&
        return _build_single_run_precursor_score_dicts(refs, n_precursors)

    inputs = _collect_global_precursor_inputs(refs, n_precursors)
    table = _build_global_precursor_feature_table(inputs, n_runs_total)
    features = _select_global_precursor_score_features(length(inputs.probabilities))
    scored = _score_global_precursor_features_oof(
        table,
        features,
    )

    score_dict = Dict{UInt32, Float32}()
    sizehint!(score_dict, nrow(table))
    @inbounds for row in axes(table, 1)
        score_dict[table.precursor_idx[row]] = scored.scores[row]
    end
    @debug_l1 "Global precursor LightGBM scored $(length(score_dict)) precursors " *
              "with $(length(features)) features; " *
              "selected semi-supervised iteration $(scored.iter)"
    return score_dict, inputs.targets
end
