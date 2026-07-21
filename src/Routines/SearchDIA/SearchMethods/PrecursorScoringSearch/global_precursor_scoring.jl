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
    :top3_prec_prob,
    :top2_logodds_score,
    :top3_logodds_score,
    :mean_prec_prob,
    :median_prec_prob,
    :std_prec_prob,
    :min_prec_prob,
    :top1_top2_gap,
    :top2_top3_gap,
    :n_runs_observed,
    :n_prob_gt_0_5,
    :n_prob_gt_0_9,
    :n_prob_gt_0_99,
]

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
        sort!(probabilities; rev = true)
        sorted_probabilities = probabilities
        n_observed = length(sorted_probabilities)
        top1 = sorted_probabilities[1]
        top2 = n_observed >= 2 ? sorted_probabilities[2] : 0.0f0
        top3 = n_observed >= 3 ? sorted_probabilities[3] : 0.0f0

        feature_columns[:empirical_global_score][row] =
            _logodds_from_sorted(sorted_probabilities, top_run_count)
        feature_columns[:top1_prec_prob][row] = top1
        feature_columns[:top2_prec_prob][row] = top2
        feature_columns[:top3_prec_prob][row] = top3
        feature_columns[:top2_logodds_score][row] =
            _logodds_from_sorted(sorted_probabilities, 2)
        feature_columns[:top3_logodds_score][row] =
            _logodds_from_sorted(sorted_probabilities, 3)
        feature_columns[:mean_prec_prob][row] = Float32(mean(probabilities))
        midpoint = (n_observed + 1) ÷ 2
        feature_columns[:median_prec_prob][row] = isodd(n_observed) ?
            sorted_probabilities[midpoint] :
            (sorted_probabilities[midpoint] + sorted_probabilities[midpoint + 1]) / 2.0f0
        feature_columns[:std_prec_prob][row] =
            n_observed > 1 ? Float32(std(probabilities)) : 0.0f0
        feature_columns[:min_prec_prob][row] = sorted_probabilities[end]
        feature_columns[:top1_top2_gap][row] = n_observed >= 2 ? top1 - top2 : 0.0f0
        feature_columns[:top2_top3_gap][row] = n_observed >= 3 ? top2 - top3 : 0.0f0
        feature_columns[:n_runs_observed][row] = Float32(n_observed)
        feature_columns[:n_prob_gt_0_5][row] = Float32(count(>(0.5f0), probabilities))
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

"""
    build_global_precursor_score_dicts(
        refs,
        n_precursors,
        n_runs_total,
        fdr_scale_factor,
    )

Build one score-distribution feature row per precursor and return out-of-fold
global scores with the corresponding target labels. Select the LightGBM scores
only when they identify more targets at 1% FDR than the empirical global score.
For a single-run search, return the individual precursor probabilities without
training a model.
"""
function build_global_precursor_score_dicts(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
    n_runs_total::Int,
    fdr_scale_factor::Float32,
)
    n_runs_total == 1 &&
        return _build_single_run_precursor_score_dicts(refs, n_precursors)

    inputs = _collect_global_precursor_inputs(refs, n_precursors)
    table = _build_global_precursor_feature_table(inputs, n_runs_total)
    scored = _score_global_features_oof(
        table,
        GLOBAL_PRECURSOR_SCORE_FEATURES,
        collect(GLOBAL_PRECURSOR_CV_FOLDS);
        scoring_name = "Global precursor",
        lgbm_hp = GLOBAL_PRECURSOR_LGBM_HP,
        min_training_class_count = GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT,
        max_train = GLOBAL_PRECURSOR_MAX_TRAIN,
    )
    selected = _select_global_scores(
        scored.scores,
        table.empirical_global_score,
        table.target;
        scoring_name = "Global precursor",
        fdr_scale_factor = fdr_scale_factor,
    )

    score_dict = Dict{UInt32, Float32}()
    sizehint!(score_dict, nrow(table))
    @inbounds for row in axes(table, 1)
        score_dict[table.precursor_idx[row]] = selected.scores[row]
    end
    @debug_l1 "Global precursor scoring scored $(length(score_dict)) precursors " *
              "with $(length(GLOBAL_PRECURSOR_SCORE_FEATURES)) features; " *
              "selected $(selected.source) after semi-supervised iteration $(scored.iter)"
    return score_dict, inputs.targets
end
