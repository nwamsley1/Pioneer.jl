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
    :n_passing_runs,
    :n_prob_gt_0_5,
    :n_prob_gt_0_9,
    :n_prob_gt_0_99,
    :observed_run_centrality_mean,
    :observed_run_centrality_max,
    :missing_run_similarity_mass_approx,
]

const GLOBAL_PRECURSOR_MONOTONE_INCREASING_FEATURES = (
    :empirical_global_score,
    :top1_prec_prob,
    :top2_prec_prob,
    :top3_prec_prob,
    :top2_logodds_score,
    :top3_logodds_score,
    :mean_prec_prob,
    :median_prec_prob,
    :min_prec_prob,
    :n_passing_runs,
    :n_prob_gt_0_5,
    :n_prob_gt_0_9,
    :n_prob_gt_0_99,
)

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
    monotone_constraints_method = "intermediate",
)

const GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT = 100
const GLOBAL_PRECURSOR_MAX_TRAIN = 1_000_000
const GLOBAL_PRECURSOR_CV_FOLDS = (UInt8(0), UInt8(1))

struct GlobalPrecursorRunScore
    ms_file_idx::UInt32
    probability::Float32
end

struct GlobalPrecursorInputs
    run_scores::Dict{UInt32, Vector{GlobalPrecursorRunScore}}
    targets::Dict{UInt32, Bool}
    folds::Dict{UInt32, UInt8}
end

"""
    build_precursor_run_similarity(refs, score_floor, n_runs_total)

Build the shared run-similarity atlas from precursors that pass the frozen
pre-global score threshold in each run. Target and decoy rows are treated
identically.
"""
function build_precursor_run_similarity(
    refs::Vector{PSMFileReference},
    score_floor::Float32,
    n_runs_total::Int,
)::RunSimilarityAtlas
    observed_precursors_by_run = Dict(
        UInt32(file_idx) => UInt32[]
        for file_idx in 1:max(n_runs_total, 0)
    )

    for (ref_position, ref) in enumerate(refs)
        ref.row_count > 0 || continue
        has_ms_file_idx = :ms_file_idx in column_names(ref)
        requested = Symbol[:precursor_idx, :prec_prob]
        has_ms_file_idx && push!(requested, :ms_file_idx)
        columns = materialize_columns(ref, requested)

        @inbounds for row in axes(columns, 1)
            file_idx = has_ms_file_idx ?
                UInt32(columns.ms_file_idx[row]) : UInt32(ref_position)
            observed_precursors = get!(
                observed_precursors_by_run,
                file_idx,
                UInt32[],
            )
            Float32(columns.prec_prob[row]) >= score_floor || continue
            push!(observed_precursors, UInt32(columns.precursor_idx[row]))
        end
    end

    return build_run_similarity(observed_precursors_by_run)
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
    run_scores = Dict{UInt32, Vector{GlobalPrecursorRunScore}}()
    targets = Dict{UInt32, Bool}()
    folds = Dict{UInt32, UInt8}()
    sizehint!(run_scores, n_precursors)
    sizehint!(targets, n_precursors)
    sizehint!(folds, n_precursors)

    for (ref_position, ref) in enumerate(refs)
        has_ms_file_idx = :ms_file_idx in column_names(ref)
        requested = Symbol[:precursor_idx, :prec_prob, :target, :cv_fold]
        has_ms_file_idx && push!(requested, :ms_file_idx)
        columns = materialize_columns(
            ref,
            requested,
        )
        @inbounds for row in axes(columns, 1)
            precursor_idx = UInt32(columns.precursor_idx[row])
            precursor_run_scores = get!(run_scores, precursor_idx) do
                GlobalPrecursorRunScore[]
            end
            push!(
                precursor_run_scores,
                GlobalPrecursorRunScore(
                    has_ms_file_idx ?
                        UInt32(columns.ms_file_idx[row]) : UInt32(ref_position),
                    Float32(columns.prec_prob[row]),
                ),
            )
            targets[precursor_idx] = Bool(columns.target[row])
            folds[precursor_idx] = UInt8(columns.cv_fold[row])
        end
    end

    isempty(run_scores) &&
        throw(ArgumentError("Global precursor scoring requires observed precursors"))
    return GlobalPrecursorInputs(run_scores, targets, folds)
end

function _build_global_precursor_feature_table(
    inputs::GlobalPrecursorInputs,
    n_runs_total::Int,
    ;
    run_similarity::Union{Nothing, RunSimilarityAtlas} = nothing,
    run_score_floor::Float32 = 0.0f0,
)
    precursor_ids = sort!(collect(keys(inputs.run_scores)))
    n_precursors = length(precursor_ids)
    top_run_count = isqrt(n_runs_total)
    feature_columns = Dict(
        feature => Vector{Float32}(undef, n_precursors)
        for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
    )
    @inbounds for (row, precursor_idx) in enumerate(precursor_ids)
        run_scores = inputs.run_scores[precursor_idx]
        sorted_probabilities = sort!(
            Float32[run_score.probability for run_score in run_scores];
            rev = true,
        )
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
        feature_columns[:mean_prec_prob][row] = Float32(mean(sorted_probabilities))
        midpoint = (n_observed + 1) ÷ 2
        feature_columns[:median_prec_prob][row] = isodd(n_observed) ?
            sorted_probabilities[midpoint] :
            (sorted_probabilities[midpoint] + sorted_probabilities[midpoint + 1]) / 2.0f0
        feature_columns[:std_prec_prob][row] =
            n_observed > 1 ? Float32(std(sorted_probabilities)) : 0.0f0
        feature_columns[:min_prec_prob][row] = sorted_probabilities[end]
        feature_columns[:top1_top2_gap][row] = n_observed >= 2 ? top1 - top2 : 0.0f0
        feature_columns[:top2_top3_gap][row] = n_observed >= 3 ? top2 - top3 : 0.0f0
        feature_columns[:n_prob_gt_0_5][row] =
            Float32(count(>(0.5f0), sorted_probabilities))
        feature_columns[:n_prob_gt_0_9][row] =
            Float32(count(>(0.9f0), sorted_probabilities))
        feature_columns[:n_prob_gt_0_99][row] =
            Float32(count(>(0.99f0), sorted_probabilities))

        centrality_sum = 0.0f0
        centrality_max = 0.0f0
        n_passing_runs = 0
        for run_score in run_scores
            run_score.probability >= run_score_floor || continue
            n_passing_runs += 1
            centrality = run_similarity === nothing ?
                0.0f0 :
                run_centrality(run_similarity, run_score.ms_file_idx)
            centrality_sum += centrality
            centrality_max = max(centrality_max, centrality)
        end
        centrality_mean = n_passing_runs > 0 ?
            centrality_sum / Float32(n_passing_runs) : 0.0f0
        n_missing_runs = max(n_runs_total - n_passing_runs, 0)
        feature_columns[:n_passing_runs][row] = Float32(n_passing_runs)
        feature_columns[:observed_run_centrality_mean][row] = centrality_mean
        feature_columns[:observed_run_centrality_max][row] = centrality_max
        feature_columns[:missing_run_similarity_mass_approx][row] =
            centrality_mean * Float32(n_missing_runs)
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

"""
    build_global_precursor_score_dicts(refs, n_precursors, n_runs_total; kwargs...)

Build one score-distribution and run-context feature row per precursor and
select out-of-fold global LightGBM scores only when they retain more total
precursor-run target IDs after applying both the configured global and
experiment-wide q-value thresholds than the empirical top-run log-odds scores.
For a single-run search, return the individual precursor probabilities without
training a model.
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
    ;
    run_similarity::Union{Nothing, RunSimilarityAtlas} = nothing,
    run_score_floor::Float32 = 0.0f0,
    q_value_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    fdr_scale_factor::Float32 = 1.0f0,
)
    n_runs_total == 1 &&
        return _build_single_run_precursor_score_dicts(refs, n_precursors)

    inputs = _collect_global_precursor_inputs(refs, n_precursors)
    table = _build_global_precursor_feature_table(
        inputs,
        n_runs_total;
        run_similarity = run_similarity,
        run_score_floor = run_score_floor,
    )
    scored = _score_global_features_oof(
        table,
        GLOBAL_PRECURSOR_SCORE_FEATURES,
        collect(GLOBAL_PRECURSOR_CV_FOLDS);
        scoring_name = "Global precursor",
        lgbm_hp = GLOBAL_PRECURSOR_LGBM_HP,
        min_training_class_count = GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT,
        max_train = GLOBAL_PRECURSOR_MAX_TRAIN,
        monotone_increasing_features =
            GLOBAL_PRECURSOR_MONOTONE_INCREASING_FEATURES,
    )
    selected = _select_global_scores(
        scored.scores,
        table.empirical_global_score,
        table.target;
        scoring_name = "Global precursor",
        experiment_wide_id_counts = Int.(table.n_passing_runs),
        q_threshold = q_value_threshold,
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
