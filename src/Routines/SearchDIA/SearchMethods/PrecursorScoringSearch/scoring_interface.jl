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
Scoring interface for precursor probability aggregation and dictionary construction.
"""

#==========================================================
Per-File Precursor Aggregation
==========================================================#

"""
    _aggregate_trace_to_precursor_probs!(df)

Per-file Bayesian aggregation of trace-level → precursor-level probabilities.
Groups by (precursor_idx, ms_file_idx). Since ms_file_idx is constant within
a single file, this is effectively grouping by precursor_idx alone.
"""
function _aggregate_trace_to_precursor_probs!(df::DataFrame)
    prob_agg = p -> begin
        trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
        clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
    end
    transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => prob_agg => :prec_prob)
    return df
end

"""
    aggregate_per_file!(refs)

Per-file precursor probability aggregation (no MBR filtering).
"""
function aggregate_per_file!(refs::Vector{PSMFileReference})
    for ref in refs
        df = load_with_sidecars(ref)
        _aggregate_trace_to_precursor_probs!(df)
        # Write only :prec_prob as a row-aligned sidecar instead of rewriting
        # the entire main file. Downstream readers locate :prec_prob via the
        # PSMFileReference's sidecar registry.
        side_path = file_path(ref) * ".prec_prob.sidecar.arrow"
        writeArrow(side_path, DataFrame(prec_prob = df.prec_prob))
        register_sidecar!(ref, side_path, [:prec_prob])
    end
    return nothing
end

#==========================================================
Additional Interface Functions
==========================================================#

"""
    logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}

Combine probabilities using a log-odds average.
The final value is converted back to a probability via the logistic function.
"""
function _logodds_from_sorted(sorted_probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    isempty(sorted_probs) && return 0.0f0
    n = min(length(sorted_probs), top_n)
    selected = @view sorted_probs[1:n]
    eps = 1f-6
    # Convert to log-odds, clip to avoid Inf or negative contribution
    logodds = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(logodds) / n
    return 1.0f0 / (1 + exp(-avg))
end

function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    sorted_probs = sort(probs; rev = true)
    return _logodds_from_sorted(sorted_probs, top_n)
end

#==========================================================
Dictionary + Sidecar Helper Functions for OOM Scoring Pipeline
==========================================================#

"""
    build_precursor_global_prob_dicts(refs, sqrt_n_runs, n_precursors; kwargs...)
    → (global_prob_dict, target_dict)

Build one score-distribution feature row per precursor, train a dedicated
global LightGBM with the existing protein-group two-fold split, and return
out-of-fold precursor scores. Features summarize the existing per-run
`prec_prob` values, including the exact empirical global score being replaced.

The previous top-`sqrt(n_runs)` log-odds equation is retained as a fallback for
small or legacy inputs that do not have both CV folds or both target classes.
"""
function _collect_global_precursor_inputs(
    refs::Vector{PSMFileReference},
    n_precursors::Int,
)
    prob_acc = Dict{UInt32, Vector{Float32}}()
    sizehint!(prob_acc, n_precursors)
    target_dict = Dict{UInt32, Bool}()
    sizehint!(target_dict, n_precursors)
    fold_dict = Dict{UInt32, UInt8}()
    sizehint!(fold_dict, n_precursors)

    nonempty_refs = filter(ref -> ref.row_count > 0, refs)
    has_cv_fold = !isempty(nonempty_refs) &&
                  all(ref -> :cv_fold in column_names(ref), nonempty_refs)

    for ref in nonempty_refs
        requested = Symbol[:precursor_idx, :prec_prob, :target]
        has_cv_fold && push!(requested, :cv_fold)
        cols_df = materialize_columns(ref, requested)
        n = nrow(cols_df)
        n == 0 && continue
        prec_ids = cols_df.precursor_idx
        prec_probs = cols_df.prec_prob
        targets = cols_df.target
        folds = has_cv_fold ? cols_df.cv_fold : nothing

        @inbounds for i in 1:n
            pid = UInt32(prec_ids[i])
            if !haskey(prob_acc, pid)
                prob_acc[pid] = Float32[]
                target_dict[pid] = Bool(targets[i])
                has_cv_fold && (fold_dict[pid] = UInt8(folds[i]))
            else
                target_dict[pid] == Bool(targets[i]) ||
                    error("Inconsistent target label for precursor $pid")
                if has_cv_fold
                    fold_dict[pid] == UInt8(folds[i]) ||
                        error("Inconsistent cv_fold for precursor $pid")
                end
            end
            push!(prob_acc[pid], Float32(prec_probs[i]))
        end
    end

    return (
        prob_acc = prob_acc,
        target_dict = target_dict,
        fold_dict = fold_dict,
        has_cv_fold = has_cv_fold,
    )
end

function _empirical_global_prob_dict(
    prob_acc::Dict{UInt32, Vector{Float32}},
    sqrt_n_runs::Int,
)
    global_prob_dict = Dict{UInt32, Float32}()
    sizehint!(global_prob_dict, length(prob_acc))
    for (pid, probs) in prob_acc
        global_prob_dict[pid] = logodds(probs, sqrt_n_runs)
    end
    return global_prob_dict
end

"""
    _build_global_precursor_feature_table(inputs, sqrt_n_runs, n_runs_total)

Summarize the existing per-run precursor probabilities for every observed
precursor. Returns `(table, features)` with exactly one row per precursor.
"""
function _build_global_precursor_feature_table(
    inputs;
    sqrt_n_runs::Int,
    n_runs_total::Int,
)
    pids = sort!(collect(keys(inputs.target_dict)))
    n_precursors = length(pids)
    table = DataFrame(
        precursor_idx = pids,
        target = Bool[inputs.target_dict[pid] for pid in pids],
        cv_fold = UInt8[inputs.fold_dict[pid] for pid in pids],
    )
    feature_columns = Dict(
        feature => Vector{Float32}(undef, n_precursors)
        for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
    )
    denominator = Float32(max(n_runs_total, 1))

    @inbounds for (row, pid) in enumerate(pids)
        probs = inputs.prob_acc[pid]
        sorted_probs = sort(probs; rev = true)
        n_observed = length(sorted_probs)
        top1 = sorted_probs[1]
        top2 = n_observed >= 2 ? sorted_probs[2] : 0.0f0
        top3 = n_observed >= 3 ? sorted_probs[3] : 0.0f0

        feature_columns[:empirical_global_score][row] =
            _logodds_from_sorted(sorted_probs, sqrt_n_runs)
        feature_columns[:top1_prec_prob][row] = top1
        feature_columns[:top2_prec_prob][row] = top2
        feature_columns[:top3_prec_prob][row] = top3
        feature_columns[:top2_logodds_score][row] = _logodds_from_sorted(sorted_probs, 2)
        feature_columns[:top3_logodds_score][row] = _logodds_from_sorted(sorted_probs, 3)
        feature_columns[:mean_prec_prob][row] = Float32(mean(probs))
        feature_columns[:median_prec_prob][row] = Float32(median(probs))
        feature_columns[:std_prec_prob][row] = n_observed > 1 ? Float32(std(probs)) : 0.0f0
        feature_columns[:min_prec_prob][row] = minimum(probs)
        feature_columns[:top1_top2_gap][row] = n_observed >= 2 ? top1 - top2 : 0.0f0
        feature_columns[:top2_top3_gap][row] = n_observed >= 3 ? top2 - top3 : 0.0f0
        feature_columns[:n_runs_observed][row] = Float32(n_observed)
        feature_columns[:observed_run_fraction][row] = Float32(n_observed) / denominator
        feature_columns[:n_prob_gt_0_5][row] = Float32(count(>(0.5f0), probs))
        feature_columns[:n_prob_gt_0_9][row] = Float32(count(>(0.9f0), probs))
        feature_columns[:n_prob_gt_0_99][row] = Float32(count(>(0.99f0), probs))
    end

    for feature in GLOBAL_PRECURSOR_SCORE_FEATURES
        table[!, feature] = feature_columns[feature]
    end
    return (table = table, features = copy(GLOBAL_PRECURSOR_SCORE_FEATURES))
end

function _log_global_precursor_feature_importances(
    models::Vector{LightGBMModel},
    features::Vector{Symbol},
)
    total_gains = Dict(feature => 0.0 for feature in features)
    any_importance = false
    for model in models
        model_importance = importance(model)
        model_importance === nothing && continue
        any_importance = true
        for (feature, gain) in model_importance
            total_gains[feature] += Float64(gain)
        end
    end
    any_importance || return nothing

    sorted_gains = sort!(collect(total_gains); by = pair -> -last(pair))
    lines = ["Global precursor LightGBM feature gains (summed across folds):"]
    for (feature, gain) in sorted_gains
        push!(lines, "    $(rpad(string(feature), 52)) $(round(Int, gain))")
    end
    @debug_l1 join(lines, "\n")
    return nothing
end

"""
    _score_global_precursor_features_oof(table, features; ...)

Fit one LightGBM per protein-group fold and return out-of-fold precursor
scores. After the initial fit with all targets and decoys, retain all decoys but
only q-value-passing targets for subsequent iterations. Returns `nothing` when
the initial training fold lacks enough examples from either class.
"""
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
    isempty(features) && return nothing
    folds = sort!(unique(table.cv_fold))
    folds == UInt8[0, 1] || return nothing

    targets = Vector{Bool}(table.target)
    training_mask = nothing
    previous_target_q01 = -1
    best_state = nothing

    for iter_idx in 1:max_iterations
        scores = Vector{Float32}(undef, nrow(table))
        models = LightGBMModel[]
        n_train_targets = 0
        n_train_decoys = 0
        iteration_valid = true

        for test_fold in folds
            test_idx = findall(==(test_fold), table.cv_fold)
            train_idx = Int[
                row for row in axes(table, 1)
                if table.cv_fold[row] != test_fold &&
                   (training_mask === nothing || training_mask[row])
            ]
            if isempty(test_idx)
                iteration_valid = false
                break
            end
            train_targets = targets[train_idx]
            fold_targets = count(train_targets)
            fold_decoys = length(train_targets) - fold_targets
            if min(fold_targets, fold_decoys) < min_training_class_count
                iteration_valid = false
                break
            end

            if length(train_idx) > max_train
                rng = Random.MersenneTwister(1776 + 1000 * iter_idx + Int(test_fold))
                Random.shuffle!(rng, train_idx)
                resize!(train_idx, max_train)
                train_targets = targets[train_idx]
                fold_targets = count(train_targets)
                fold_decoys = length(train_targets) - fold_targets
                if min(fold_targets, fold_decoys) < min_training_class_count
                    iteration_valid = false
                    break
                end
            end

            classifier = build_lightgbm_classifier(; lgbm_hp...)
            model = fit_lightgbm_model(
                classifier,
                view(table, train_idx, features),
                targets[train_idx];
                positive_label = true,
            )
            scores[test_idx] .= lightgbm_predict(
                model,
                view(table, test_idx, features);
                output_type = Float32,
            )
            n_train_targets += fold_targets
            n_train_decoys += fold_decoys
            push!(models, model)
        end

        if !iteration_valid
            if iter_idx == 1
                return nothing
            end
            @debug_l1 "Global precursor semi-supervised stopping: iter $iter_idx " *
                       "does not retain enough targets and decoys in both folds; " *
                       "using iter $(best_state.iter)"
            break
        end

        scores .= Float32.(clamp.(scores, 1f-6, 1f0 - 1f-4))
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
            iter = iter_idx,
        )
        @debug_l1 "Global precursor semi-supervised iter $iter_idx: " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "q≤.01 targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        if iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01,
            metrics.target_q01;
            min_fraction = min_gain,
        )
            best_state = _scoring_better_iteration_state(best_state, current_state)
            @debug_l1 "Global precursor semi-supervised stopping: iter $iter_idx " *
                       "q≤.01 targets=$(metrics.target_q01) did not improve by " *
                       "$(round(100 * min_gain, digits=2))% over $previous_target_q01; " *
                       "using iter $(best_state.iter) with q≤.01 targets=$(best_state.target_q01)"
            break
        end

        best_state = _scoring_better_iteration_state(best_state, current_state)
        if iter_idx == max_iterations
            @debug_l1 "Global precursor semi-supervised stopping: hit max iterations " *
                       "$max_iterations; using iter $(best_state.iter) with " *
                       "q≤.01 targets=$(best_state.target_q01)"
            break
        end

        previous_target_q01 = metrics.target_q01
        training_mask = metrics.training_mask
    end

    best_state === nothing && return nothing
    _log_global_precursor_feature_importances(best_state.models, features)
    return best_state
end

function build_precursor_global_prob_dicts(
    refs::Vector{PSMFileReference},
    sqrt_n_runs::Int,
    n_precursors::Int;
    n_runs_total::Int = max(length(refs), 1),
    lgbm_hp::NamedTuple = GLOBAL_PRECURSOR_LGBM_HP,
    min_training_class_count::Int = GLOBAL_PRECURSOR_MIN_TRAINING_CLASS_COUNT,
    max_train::Int = GLOBAL_PRECURSOR_MAX_TRAIN,
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    inputs = _collect_global_precursor_inputs(refs, n_precursors)

    fallback(reason) = begin
        @debug_l1 "Global precursor LightGBM fallback to empirical log-odds: $reason"
        return _empirical_global_prob_dict(inputs.prob_acc, sqrt_n_runs), inputs.target_dict
    end

    isempty(inputs.prob_acc) && return fallback("no observed precursors")
    inputs.has_cv_fold || return fallback("cv_fold column unavailable")
    unique_folds = sort!(unique(collect(values(inputs.fold_dict))))
    unique_folds == UInt8[0, 1] || return fallback("requires both CV folds")

    for test_fold in unique_folds
        train_pids = UInt32[pid for pid in keys(inputs.target_dict)
                            if inputs.fold_dict[pid] != test_fold]
        n_targets = count(pid -> inputs.target_dict[pid], train_pids)
        n_decoys = length(train_pids) - n_targets
        min(n_targets, n_decoys) >= min_training_class_count ||
            return fallback("insufficient target/decoy examples in a training fold")
    end

    built = _build_global_precursor_feature_table(
        inputs;
        sqrt_n_runs = sqrt_n_runs,
        n_runs_total = n_runs_total,
    )

    scored = _score_global_precursor_features_oof(
        built.table,
        built.features;
        lgbm_hp = lgbm_hp,
        min_training_class_count = min_training_class_count,
        max_train = max_train,
        train_q_threshold = train_q_threshold,
        stop_q_threshold = stop_q_threshold,
        min_gain = min_gain,
        max_iterations = max_iterations,
    )
    scored === nothing && return fallback("global LightGBM training requirements not met")

    global_prob_dict = Dict{UInt32, Float32}()
    sizehint!(global_prob_dict, nrow(built.table))
    @inbounds for row in 1:nrow(built.table)
        global_prob_dict[built.table.precursor_idx[row]] = scored.scores[row]
    end
    @debug_l1 "Global precursor LightGBM scored $(length(global_prob_dict)) precursors " *
              "with $(length(built.features)) score-distribution features; " *
              "selected semi-supervised iter $(scored.iter)"
    return global_prob_dict, inputs.target_dict
end

"""
    build_global_qval_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global q-values from a score dictionary without any file I/O.
"""
function build_global_qval_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[get(target_dict, pid, false) for pid in pids]

    # Sort descending by score
    perm = sortperm(scores; rev=true)
    permute!(pids, perm)
    permute!(scores, perm)
    permute!(targets, perm)

    # Compute q-values
    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals; fdr_scale_factor=fdr_scale)

    # Build dictionary
    qval_dict = Dict{UInt32, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[pids[i]] = qvals[i]
    end
    return qval_dict
end

"""
    build_global_pep_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global posterior error probabilities (local FDR) from a score dictionary
without any file I/O. Parallel to `build_global_qval_dict_from_scores` but uses
`get_PEP!` (PAVA-fit) instead of cumulative q-values.
"""
function build_global_pep_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[get(target_dict, pid, false) for pid in pids]

    peps = Vector{Float32}(undef, n)
    get_PEP!(scores, targets, peps; doSort=true, fdr_scale_factor=fdr_scale)

    pep_dict = Dict{UInt32, Float32}()
    sizehint!(pep_dict, n)
    for i in 1:n
        pep_dict[pids[i]] = peps[i]
    end
    return pep_dict
end


"""
    write_score_sidecars(refs, columns; temp_prefix) → Vector{PSMFileReference}

Extract only the named columns from each file into a temporary Arrow sidecar file.
"""
function write_score_sidecars(
    refs::Vector{<:FileReference},
    columns::Vector{Symbol};
    temp_prefix::String = "sidecar"
)
    sidecar_refs = PSMFileReference[]
    for ref in refs
        # Use materialize_columns so columns are pulled from main OR any
        # registered sidecar (e.g. :prec_prob now lives in a sidecar after
        # aggregate_per_file!).
        df = ref isa PSMFileReference ? materialize_columns(ref, columns) :
             DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))[!, columns]
        nrow(df) == 0 && continue
        temp_path = tempname() * "_$(temp_prefix).arrow"
        writeArrow(temp_path, df)
        push!(sidecar_refs, PSMFileReference(temp_path))
    end
    return sidecar_refs
end

"""
    build_qvalue_spline_from_refs(refs, score_col, merged_path; ...) → Union{Nothing, NamedTuple}

Encapsulates the full sidecar lifecycle: write → sort → merge → cleanup → spline computation.
"""
function build_qvalue_spline_from_refs(
    refs::Vector{<:FileReference},
    score_col::Symbol,
    merged_path::String;
    batch_size::Int = 10_000_000,
    compute_pep::Bool = false,
    min_pep_points_per_bin::Int = 100,
    fdr_scale_factor::Float32 = 1.0f0,
    temp_prefix::String = "sidecar"
)
    sidecar_refs = write_score_sidecars(refs, [score_col, :target]; temp_prefix=temp_prefix)
    isempty(sidecar_refs) && return nothing

    try
        sort_file_by_keys!(sidecar_refs, score_col, :target; reverse=[true, true])
        stream_sorted_merge(sidecar_refs, merged_path, score_col, :target;
                           batch_size=batch_size, reverse=[true, true])
    finally
        GC.gc(false)
        for ref in sidecar_refs
            safeRm(file_path(ref), nothing; force=true)
        end
    end

    qval_spline = get_qvalue_spline(merged_path, score_col, false;
        min_pep_points_per_bin=min_pep_points_per_bin,
        fdr_scale_factor=fdr_scale_factor)

    pep_interp = if compute_pep
        get_pep_interpolation(merged_path, score_col;
            fdr_scale_factor=fdr_scale_factor)
    else
        nothing
    end

    return (; qval_spline, pep_interp)
end
