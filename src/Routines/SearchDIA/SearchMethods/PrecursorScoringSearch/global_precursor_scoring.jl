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
    :std_prec_prob,
    :min_prec_prob,
    :top1_top2_gap,
    :top2_top3_gap,
    :n_runs_observed,
    :n_prob_gt_0_5,
    :n_prob_gt_0_9,
    :n_prob_gt_0_99,
]

# Per-precursor cross-run features are accumulated ONLINE in a single streaming pass into a
# dense `Vector{GlobalAcc}` indexed by precursor_idx — no `Dict{pid → Vector}` and no per-precursor
# heap vectors. `median_prec_prob` was dropped precisely because it is the only feature that needs
# the full per-precursor distribution; every remaining feature is a bounded online reduction
# (top-K tuple, Welford mean/variance, min, threshold counts). The top-K logodds uses the top
# `min(isqrt(n_runs), GLOBAL_TOPK_CAP)` scores, so the accumulator keeps exactly the top
# GLOBAL_TOPK_CAP. For cohorts with isqrt(n_runs) > GLOBAL_TOPK_CAP (>100 runs) this caps the
# logodds window; ≤100-run cohorts are unaffected.
const GLOBAL_TOPK_CAP = 10

struct GlobalAcc
    top::NTuple{GLOBAL_TOPK_CAP, Float32}   # top-K prec_prob, descending; trailing 0f0 for empty slots
    mean::Float32                            # Welford running mean
    m2::Float32                              # Welford running sum of squared deviations
    minv::Float32
    count::Int32
    n_gt50::Int32
    n_gt90::Int32
    n_gt99::Int32
end

const EMPTY_GLOBAL_ACC = GlobalAcc(
    ntuple(_ -> 0.0f0, GLOBAL_TOPK_CAP), 0.0f0, 0.0f0, Inf32,
    Int32(0), Int32(0), Int32(0), Int32(0),
)

# Insert `x` into a descending NTuple keeping the top-N, dropping the smallest.
# `@generated` so each branch returns a tuple built from COMPILE-TIME-literal positions —
# a runtime insertion index would make the tuple construction dynamic and heap-allocate per
# call (measured: ~112 B/call). This version is 0 bytes/call. Same idea as `TopK4`, generalized.
@generated function _topk_insert(t::NTuple{N, Float32}, x::Float32) where {N}
    branches = Expr[]
    for p in 1:N
        elems = Any[i < p ? :(t[$i]) : (i == p ? :x : :(t[$(i - 1)])) for i in 1:N]
        push!(branches, :(x > t[$p] && return ($(elems...),)))
    end
    quote
        @inbounds begin
            $(branches...)
            return t
        end
    end
end

@inline function _update_global_acc(a::GlobalAcc, p::Float32)
    count = a.count + Int32(1)
    delta = p - a.mean
    mean = a.mean + delta / count
    m2 = a.m2 + delta * (p - mean)
    return GlobalAcc(
        _topk_insert(a.top, p),
        mean, m2,
        min(a.minv, p),
        count,
        a.n_gt50 + (p > 0.5f0  ? Int32(1) : Int32(0)),
        a.n_gt90 + (p > 0.9f0  ? Int32(1) : Int32(0)),
        a.n_gt99 + (p > 0.99f0 ? Int32(1) : Int32(0)),
    )
end

# Average log-odds of the top `top_n` (capped by the observed count and the stored K), then sigmoid.
# Bit-identical to `_logodds_from_sorted` over the same top elements.
@inline function _logodds_from_topk(top::NTuple{K, Float32}, n_observed::Int, top_n::Int) where {K}
    m = min(n_observed, top_n, K)
    floor_val = 0.1f0
    ceil_val = 1.0f0 - 1.0f-6
    s = 0.0f0
    @inbounds for i in 1:m
        c = clamp(top[i], floor_val, ceil_val)
        s += log(c / (1.0f0 - c))
    end
    return 1.0f0 / (1.0f0 + exp(-s / m))
end

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
    acc::Vector{GlobalAcc}      # indexed by precursor_idx + 1
    targets::Vector{Bool}       # indexed by precursor_idx + 1
    folds::Vector{UInt8}        # indexed by precursor_idx + 1
    n_precursors::Int
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

# Single streaming pass: fold each row into the dense per-precursor accumulator. `prec_prob`
# is a main column of the merged file (written in Step 1b), so one Arrow.Table open reads all
# four columns — no sidecar join, no column materialization. Isolated as a function barrier so
# the concrete Arrow column types specialize and element access does not box.
function _accumulate_global_inputs!(
    acc::Vector{GlobalAcc},
    targets::Vector{Bool},
    folds::Vector{UInt8},
    pid_col::AbstractVector,
    prob_col::AbstractVector,
    target_col::AbstractVector,
    fold_col::AbstractVector,
)
    # Plain integer indices — NOT `eachindex(pid_col)`. For a multi-batch Arrow file the columns
    # are `ChainedVector`s and `eachindex` returns a column-specific `IndexIterator`; using pid's
    # iterator to index the OTHER columns mis-reads across chunk boundaries (a merged file written
    # as fold0+fold1 record batches hit exactly this). Integer indexing is correct for every Arrow
    # column type since all four share the file's row count.
    @inbounds for i in 1:length(pid_col)
        pid1 = Int(pid_col[i]) + 1
        acc[pid1] = _update_global_acc(acc[pid1], Float32(prob_col[i]))
        targets[pid1] = Bool(target_col[i])
        folds[pid1] = UInt8(fold_col[i])
    end
    return nothing
end

function _collect_global_precursor_inputs(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
)
    acc = fill(EMPTY_GLOBAL_ACC, n_precursors + 1)
    targets = zeros(Bool, n_precursors + 1)
    folds = zeros(UInt8, n_precursors + 1)

    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        hasproperty(tbl, :prec_prob) ||
            error("_collect_global_precursor_inputs: :prec_prob missing from $(file_path(ref)); " *
                  "it should be a main column after Step 1b")
        _accumulate_global_inputs!(
            acc, targets, folds,
            tbl.precursor_idx, tbl.prec_prob, tbl.target, tbl.cv_fold,
        )
    end

    any(a -> a.count > 0, acc) ||
        throw(ArgumentError("Global precursor scoring requires observed precursors"))
    return GlobalPrecursorInputs(acc, targets, folds, n_precursors)
end

function _build_global_precursor_feature_table(
    inputs::GlobalPrecursorInputs,
    n_runs_total::Int,
)
    acc = inputs.acc
    # Observed precursors, ascending precursor_idx (iterating pid order gives them pre-sorted,
    # matching the old `sort!(collect(keys(...)))`). Count first, then fill — no push!.
    n_obs = 0
    @inbounds for pid1 in 2:(inputs.n_precursors + 1)
        acc[pid1].count > 0 && (n_obs += 1)
    end
    precursor_ids = Vector{UInt32}(undef, n_obs)
    row_targets = Vector{Bool}(undef, n_obs)
    row_folds = Vector{UInt8}(undef, n_obs)
    k = 0
    @inbounds for pid1 in 2:(inputs.n_precursors + 1)
        acc[pid1].count > 0 || continue
        k += 1
        precursor_ids[k] = UInt32(pid1 - 1)
        row_targets[k] = inputs.targets[pid1]
        row_folds[k] = inputs.folds[pid1]
    end

    top_run_count = min(isqrt(n_runs_total), GLOBAL_TOPK_CAP)
    feature_columns = _compute_global_feature_columns(acc, precursor_ids, top_run_count)

    table = DataFrame(
        precursor_idx = precursor_ids,
        target = row_targets,
        cv_fold = row_folds,
    )
    for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
        table[!, feature] = feature_columns[feature]
    end
    return table
end

# Function barrier: everything enters as a typed argument, so nothing is captured from an
# outer scope. Building the feature vectors inline in the caller boxed `n_obs` (loop-assigned
# AND captured by the column-allocating comprehension), which deoptimized the whole per-precursor
# loop into ~1.3 KB/precursor of boxing. Here `n_obs` is a fresh single-assigned local.
function _compute_global_feature_columns(
    acc::Vector{GlobalAcc},
    precursor_ids::Vector{UInt32},
    top_run_count::Int,
)
    n_obs = length(precursor_ids)
    feature_columns = Dict(
        feature => Vector{Float32}(undef, n_obs)
        for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
    )
    empirical  = feature_columns[:empirical_global_score]
    c_top1     = feature_columns[:top1_prec_prob]
    c_top2     = feature_columns[:top2_prec_prob]
    c_top3     = feature_columns[:top3_prec_prob]
    c_lo2      = feature_columns[:top2_logodds_score]
    c_lo3      = feature_columns[:top3_logodds_score]
    c_mean     = feature_columns[:mean_prec_prob]
    c_std      = feature_columns[:std_prec_prob]
    c_min      = feature_columns[:min_prec_prob]
    c_gap12    = feature_columns[:top1_top2_gap]
    c_gap23    = feature_columns[:top2_top3_gap]
    c_nobs     = feature_columns[:n_runs_observed]
    c_gt50     = feature_columns[:n_prob_gt_0_5]
    c_gt90     = feature_columns[:n_prob_gt_0_9]
    c_gt99     = feature_columns[:n_prob_gt_0_99]

    @inbounds for row in 1:n_obs
        a = acc[Int(precursor_ids[row]) + 1]
        top = a.top
        n_observed = Int(a.count)
        top1 = top[1]
        top2 = n_observed >= 2 ? top[2] : 0.0f0
        top3 = n_observed >= 3 ? top[3] : 0.0f0

        empirical[row] = _logodds_from_topk(top, n_observed, top_run_count)
        c_top1[row] = top1
        c_top2[row] = top2
        c_top3[row] = top3
        c_lo2[row] = _logodds_from_topk(top, n_observed, 2)
        c_lo3[row] = _logodds_from_topk(top, n_observed, 3)
        c_mean[row] = a.mean
        c_std[row] = n_observed > 1 ? sqrt(a.m2 / (n_observed - 1)) : 0.0f0
        c_min[row] = a.minv
        c_gap12[row] = n_observed >= 2 ? top1 - top2 : 0.0f0
        c_gap23[row] = n_observed >= 3 ? top2 - top3 : 0.0f0
        c_nobs[row] = Float32(n_observed)
        c_gt50[row] = Float32(a.n_gt50)
        c_gt90[row] = Float32(a.n_gt90)
        c_gt99[row] = Float32(a.n_gt99)
    end
    return feature_columns
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
    scored = _score_global_precursor_features_oof(
        table,
        GLOBAL_PRECURSOR_SCORE_FEATURES,
    )

    score_dict = Dict{UInt32, Float32}()
    target_dict = Dict{UInt32, Bool}()
    sizehint!(score_dict, nrow(table))
    sizehint!(target_dict, nrow(table))
    @inbounds for row in axes(table, 1)
        pid = table.precursor_idx[row]
        score_dict[pid] = scored.scores[row]
        target_dict[pid] = table.target[row]
    end
    @debug_l1 "Global precursor LightGBM scored $(length(score_dict)) precursors " *
              "with $(length(GLOBAL_PRECURSOR_SCORE_FEATURES)) features; " *
              "selected semi-supervised iteration $(scored.iter)"
    return score_dict, target_dict
end
