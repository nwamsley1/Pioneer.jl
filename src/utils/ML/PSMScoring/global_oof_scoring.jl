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

const SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD = 0.03f0
const SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD = 0.01f0
const SCORING_SEMISUPERVISED_MIN_TARGET_GAIN = 0.01f0
const SCORING_SEMISUPERVISED_MAX_ITERATIONS = 8

function _count_targets_at_qvalue(
    scores::AbstractVector{<:AbstractFloat},
    targets::AbstractVector{Bool};
    q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    fdr_scale_factor::Float32 = 1.0f0,
)
    q_values = Vector{Float32}(undef, length(scores))
    get_qvalues!(
        scores,
        targets,
        q_values;
        fdr_scale_factor = fdr_scale_factor,
    )
    return count(
        target && q_value <= q_threshold
        for (target, q_value) in zip(targets, q_values)
    )
end

"""
    _select_global_scores(model_scores, empirical_scores, targets; kwargs...)

Select model scores only when more targets pass the configured global q-value
threshold than with the empirical scores. Experiment-wide q-values are not
used. Ties retain the empirical scores.
"""
function _select_global_scores(
    model_scores::Vector{Float32},
    empirical_scores::AbstractVector{Float32},
    targets::AbstractVector{Bool};
    scoring_name::AbstractString,
    q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    fdr_scale_factor::Float32 = 1.0f0,
)
    model_target_count = _count_targets_at_qvalue(
        model_scores,
        targets;
        q_threshold = q_threshold,
        fdr_scale_factor = fdr_scale_factor,
    )
    empirical_target_count = _count_targets_at_qvalue(
        empirical_scores,
        targets;
        q_threshold = q_threshold,
        fdr_scale_factor = fdr_scale_factor,
    )
    use_model = model_target_count > empirical_target_count
    source = use_model ? :lightgbm : :empirical
    selected_scores = use_model ? model_scores : empirical_scores

    @debug_l1 "$scoring_name score selection at global q≤$q_threshold: " *
              "LightGBM targets=$model_target_count, " *
              "empirical targets=$empirical_target_count; selected $source"
    return (
        scores = selected_scores,
        source = source,
        model_target_count = model_target_count,
        empirical_target_count = empirical_target_count,
    )
end

function _scoring_semisupervised_train_mask(
    targets::AbstractVector{Bool},
    q_values::AbstractVector{<:Real};
    q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
)
    n = length(targets)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = !targets[i] || Float32(q_values[i]) <= q_threshold
    end
    return mask
end

function _scoring_semisupervised_metrics_and_mask(
    scores::AbstractVector{<:Real},
    targets::AbstractVector{Bool};
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
)
    n = length(scores)
    order = sortperm(scores; rev = true, alg = QuickSort)
    training_mask = BitVector(undef, n)
    total_targets = count(targets)
    total_decoys = n - total_targets
    suffix_targets = 0
    suffix_decoys = 0
    min_q = Inf32
    target_q01 = 0
    decoy_q01 = 0

    @inbounds for sorted_pos in n:-1:1
        i = Int(order[sorted_pos])
        prefix_targets = total_targets - suffix_targets
        prefix_decoys = total_decoys - suffix_decoys
        raw_q = prefix_targets == 0 ? Inf32 : Float32(prefix_decoys) / Float32(prefix_targets)
        min_q = min(min_q, raw_q)
        train_q_pass = min_q <= train_q_threshold
        stop_q_pass = min_q <= stop_q_threshold
        is_target = targets[i]

        if stop_q_pass
            if is_target
                target_q01 += 1
            else
                decoy_q01 += 1
            end
        end
        training_mask[i] = !is_target || train_q_pass

        if is_target
            suffix_targets += 1
        else
            suffix_decoys += 1
        end
    end

    return (
        training_mask = training_mask,
        target_q01 = target_q01,
        decoy_q01 = decoy_q01,
    )
end

function _scoring_training_target_decoy_counts(
    targets::AbstractVector{Bool},
    training_mask::Union{Nothing, AbstractVector{Bool}},
)
    if training_mask === nothing
        n_targets = count(targets)
        return n_targets, length(targets) - n_targets
    end
    n_targets = 0
    n_decoys = 0
    @inbounds for i in eachindex(targets)
        training_mask[i] || continue
        if targets[i]
            n_targets += 1
        else
            n_decoys += 1
        end
    end
    return n_targets, n_decoys
end

function _scoring_target_gain_sufficient(
    previous_targets::Integer,
    current_targets::Integer;
    min_fraction::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
)
    previous_targets <= 0 && return current_targets > previous_targets
    return Float64(current_targets) >= Float64(previous_targets) * (1.0 + Float64(min_fraction))
end

function _scoring_better_iteration_state(previous_state, current_state)
    previous_state === nothing && return current_state
    current_state === nothing && return previous_state
    return current_state.target_q01 >= previous_state.target_q01 ? current_state : previous_state
end

function _split_scoring_train_masks(
    row_counts::AbstractVector{<:Integer},
    training_mask::AbstractVector{Bool},
)
    masks = Vector{BitVector}(undef, length(row_counts))
    offset = 0
    for (file_idx, n_integer) in enumerate(row_counts)
        n = Int(n_integer)
        masks[file_idx] = BitVector(@view training_mask[(offset + 1):(offset + n)])
        offset += n
    end
    return masks
end

function _log_global_feature_importances(
    models::Vector{LightGBMModel},
    features::Vector{Symbol},
    scoring_name::AbstractString,
)
    total_gains = Dict(feature => 0.0 for feature in features)
    for model in models
        for (feature, gain) in importance(model)
            total_gains[feature] += Float64(gain)
        end
    end

    sorted_gains = sort!(collect(total_gains); by = pair -> -last(pair))
    lines = ["$scoring_name LightGBM feature gains (summed across folds):"]
    for (feature, gain) in sorted_gains
        push!(lines, "    $(rpad(string(feature), 52)) $(round(Int, gain))")
    end
    @debug_l1 join(lines, "\n")
    return nothing
end

"""
    _score_global_features_oof(table, features, cv_folds; kwargs...)

Train one deterministic LightGBM model per held-out CV fold and return
out-of-fold scores for every row. Semi-supervised iterations retain all decoys
and targets passing the configured training q-value threshold.
"""
function _score_global_features_oof(
    table::DataFrame,
    features::Vector{Symbol},
    cv_folds::Vector{UInt8};
    scoring_name::AbstractString,
    lgbm_hp::NamedTuple,
    min_training_class_count::Int,
    max_train::Int,
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    observed_folds = sort!(unique(Vector{UInt8}(table.cv_fold)))
    observed_folds == cv_folds || throw(ArgumentError(
        "$scoring_name scoring requires CV folds $(join(cv_folds, ", "))",
    ))
    length(cv_folds) >= 2 || throw(ArgumentError(
        "$scoring_name scoring requires at least two CV folds",
    ))

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

        for test_fold in cv_folds
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
                        "$scoring_name scoring requires at least " *
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
            @debug_l1 "$scoring_name semi-supervised stopping: iteration " *
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
        @debug_l1 "$scoring_name semi-supervised iteration $iteration: " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "q≤.01 targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        if iteration > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01,
            metrics.target_q01;
            min_fraction = min_gain,
        )
            best_state = _scoring_better_iteration_state(best_state, current_state)
            @debug_l1 "$scoring_name semi-supervised stopping: iteration " *
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

    _log_global_feature_importances(best_state.models, features, scoring_name)
    return best_state
end
