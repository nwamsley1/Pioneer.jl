# Learned chromatogram-boundary ranking utilities.
#
# The learner is deliberately scoped to candidate-level peak-shape features and
# self-supervised quant consistency targets. Context features such as charge,
# precursor m/z, file intensity percentile, and transmission state are avoided so
# the model learns how good the proposed bounds look rather than which precursor
# class it came from.

const BOUNDARY_CANDIDATE_FEATURES = Symbol[
    :candidate_width,
    :endpoint_height_fraction,
    :peak_prominence_score,
    :log2_smoothed_apex_weight,
    :asymmetry_penalty,
    :irt_asymmetry_delta,
    :baseline_internal_trough_score,
    :baseline_internal_trough_log2_ratio,
    :left_excluded_signal_fraction,
    :right_excluded_signal_fraction,
    :left_boundary_recovery_fraction,
    :right_boundary_recovery_fraction,
    :left_outside_peak_fraction,
    :right_outside_peak_fraction,
    :internal_dip_recovery_score,
    :left_edge_valley_log2_ratio,
    :right_edge_valley_log2_ratio,
    :included_nonapex_max_log2_ratio,
    :included_nonapex_max_irt_distance,
    :mainsearch_left_bound_delta,
    :mainsearch_right_bound_delta,
]
const BOUNDARY_RANKING_MAX_LABEL = 30
const DEBUG_BOUNDARY_CANDIDATE_TARGETS = Ref{Set{Tuple{UInt32,UInt16}}}(
    Set([
        (UInt32(370714), UInt16(1)),
        (UInt32(909488), UInt16(3)),
        (UInt32(1047223), UInt16(3)),
    ])
)

@inline function _boundary_log_area(area)
    area_f = Float32(area)
    return isfinite(area_f) && area_f > 0.0f0 ? log2(area_f) : Inf32
end

function _median_or_nothing(values::Vector{Float32})
    isempty(values) && return nothing
    return Float32(median(values))
end

function _median_excluding_group(entries::Vector{Tuple{UInt64, Float32}}, excluded_group::UInt64)
    vals = Float32[]
    sizehint!(vals, length(entries))
    for (group_id, value) in entries
        group_id == excluded_group && continue
        isfinite(value) || continue
        push!(vals, value)
    end
    return _median_or_nothing(vals)
end

function _get_column_or_fill(df::AbstractDataFrame, col::Symbol, default)
    hasproperty(df, col) && return df[!, col]
    return fill(default, nrow(df))
end

"""
    sample_boundary_candidate_groups(candidates; max_groups, rng, group_col)

Randomly sample whole boundary-candidate groups. Keeping complete groups avoids
training examples where one candidate for a chromatogram is present but its
alternatives were discarded.
"""
function sample_boundary_candidate_groups(
    candidates::DataFrame;
    max_groups::Integer,
    rng::AbstractRNG = Random.default_rng(),
    group_col::Symbol = :boundary_group_id,
)
    max_groups <= 0 && return candidates[1:0, :]
    nrow(candidates) == 0 && return copy(candidates)

    groups = unique(candidates[!, group_col])
    length(groups) <= max_groups && return copy(candidates)

    shuffled = collect(groups)
    Random.shuffle!(rng, shuffled)
    keep = Set(shuffled[1:Int(max_groups)])
    mask = [group_id in keep for group_id in candidates[!, group_col]]
    return candidates[mask, :]
end

function _fallback_rows_by_group(
    candidates::AbstractDataFrame,
    group_col::Symbol,
    fallback_col::Symbol,
)
    groups = candidates[!, group_col]
    fallback = _get_column_or_fill(candidates, fallback_col, false)
    rows_by_group = Dict{UInt64, Vector{Int}}()

    for i in 1:nrow(candidates)
        group_id = UInt64(groups[i])
        push!(get!(() -> Int[], rows_by_group, group_id), i)
    end

    fallback_row = Dict{UInt64, Int}()
    for (group_id, rows) in rows_by_group
        chosen = rows[1]
        for row in rows
            if Bool(fallback[row])
                chosen = row
                break
            end
        end
        fallback_row[group_id] = chosen
    end

    return rows_by_group, fallback_row
end

function _fallback_log_entries(
    candidates::AbstractDataFrame,
    fallback_row::Dict{UInt64, Int},
    group_col::Symbol,
    precursor_col::Symbol,
    file_col::Symbol,
    protein_col::Symbol,
    area_col::Symbol,
    isotope_col::Union{Nothing, Symbol},
)
    pid_col = candidates[!, precursor_col]
    file_ids = candidates[!, file_col]
    protein_keys = candidates[!, protein_col]
    isotope_keys = isotope_col === nothing ? nothing : candidates[!, isotope_col]
    quant_selected = _get_column_or_fill(candidates, :quant_trace_selected, true)
    areas = candidates[!, area_col]

    by_precursor = Dict{UInt32, Vector{Tuple{UInt64, Float32}}}()
    by_precursor_run_isotope = Dict{Tuple{UInt32, UInt16}, Vector{Tuple{Any, UInt64, Float32}}}()
    entries = NamedTuple[]

    for (group_id, row) in fallback_row
        log_area = _boundary_log_area(areas[row])
        isfinite(log_area) || continue

        pid = UInt32(pid_col[row])
        file_id = UInt16(file_ids[row])
        protein_key = protein_keys[row]
        isotope_key = isotope_col === nothing ? nothing : isotope_keys[row]
        is_quant_selected = Bool(quant_selected[row])

        if is_quant_selected
            push!(get!(() -> Tuple{UInt64, Float32}[], by_precursor, pid), (group_id, log_area))
        end
        if isotope_col !== nothing
            push!(
                get!(() -> Tuple{Any, UInt64, Float32}[], by_precursor_run_isotope, (pid, file_id)),
                (isotope_key, group_id, log_area),
            )
        end
        push!(entries, (
            group_id = UInt64(group_id),
            row = row,
            precursor_idx = pid,
            ms_file_idx = file_id,
            protein_key = protein_key,
            isotope_key = isotope_key,
            quant_trace_selected = is_quant_selected,
            log_area = log_area,
        ))
    end

    return entries, by_precursor, by_precursor_run_isotope
end

function _crossrun_reference_residuals(
    entries::Vector{NamedTuple},
    by_precursor::Dict{UInt32, Vector{Tuple{UInt64, Float32}}},
)
    residuals = Dict{Tuple{Any, UInt16}, Vector{Tuple{UInt32, UInt64, Float32}}}()
    for entry in entries
        entry.quant_trace_selected || continue
        center = _median_excluding_group(by_precursor[entry.precursor_idx], entry.group_id)
        center === nothing && continue
        residual = Float32(entry.log_area - center)
        push!(
            get!(() -> Tuple{UInt32, UInt64, Float32}[], residuals,
                (entry.protein_key, entry.ms_file_idx)),
            (entry.precursor_idx, entry.group_id, residual),
        )
    end
    return residuals
end

function _crossrun_loss_for_candidate(
    row::Int,
    candidates::AbstractDataFrame,
    by_precursor::Dict{UInt32, Vector{Tuple{UInt64, Float32}}},
    crossrun_residuals::Dict{Tuple{Any, UInt16}, Vector{Tuple{UInt32, UInt64, Float32}}},
    group_col::Symbol,
    precursor_col::Symbol,
    file_col::Symbol,
    protein_col::Symbol,
    area_col::Symbol,
    min_crossrun_refs::Integer,
)
    quant_selected = _get_column_or_fill(candidates, :quant_trace_selected, true)
    Bool(quant_selected[row]) || return Inf32

    group_id = UInt64(candidates[row, group_col])
    pid = UInt32(candidates[row, precursor_col])
    center = get(by_precursor, pid, nothing)
    center === nothing && return Inf32
    precursor_center = _median_excluding_group(center, group_id)
    precursor_center === nothing && return Inf32

    refs = get(crossrun_residuals, (candidates[row, protein_col], UInt16(candidates[row, file_col])), nothing)
    refs === nothing && return Inf32

    vals = Float32[]
    for (ref_pid, ref_group, residual) in refs
        ref_pid == pid && continue
        ref_group == group_id && continue
        isfinite(residual) || continue
        push!(vals, residual)
    end
    length(vals) < min_crossrun_refs && return Inf32

    expected = Float32(precursor_center + median(vals))
    log_area = _boundary_log_area(candidates[row, area_col])
    isfinite(log_area) || return Inf32
    return abs(log_area - expected)
end

function _isotope_loss_for_candidate(
    row::Int,
    candidates::AbstractDataFrame,
    by_precursor_run_isotope::Dict{Tuple{UInt32, UInt16}, Vector{Tuple{Any, UInt64, Float32}}},
    group_col::Symbol,
    precursor_col::Symbol,
    file_col::Symbol,
    area_col::Symbol,
    isotope_col::Union{Nothing, Symbol},
    min_isotope_refs::Integer,
)
    isotope_col === nothing && return Inf32

    group_id = UInt64(candidates[row, group_col])
    pid = UInt32(candidates[row, precursor_col])
    isotope_key = candidates[row, isotope_col]
    refs = get(by_precursor_run_isotope, (pid, UInt16(candidates[row, file_col])), nothing)
    refs === nothing && return Inf32

    vals = Float32[]
    for (ref_isotope, ref_group, log_area) in refs
        ref_isotope == isotope_key && continue
        ref_group == group_id && continue
        isfinite(log_area) || continue
        push!(vals, log_area)
    end
    length(vals) < min_isotope_refs && return Inf32

    expected = Float32(median(vals))
    log_area = _boundary_log_area(candidates[row, area_col])
    isfinite(log_area) || return Inf32
    return abs(log_area - expected)
end

function _combined_boundary_loss(
    crossrun_loss::Float32,
    isotope_loss::Float32,
    is_fallback::Bool,
    crossrun_weight::Float32,
    isotope_weight::Float32,
)
    total = 0.0f0
    weight = 0.0f0
    if isfinite(crossrun_loss)
        total += crossrun_weight * crossrun_loss
        weight += crossrun_weight
    end
    if isfinite(isotope_loss)
        total += isotope_weight * isotope_loss
        weight += isotope_weight
    end
    weight > 0.0f0 && return total / weight

    # No self-supervised objective is available for this candidate group.
    # Leave it unlabeled so it does not teach the model that fallback bounds
    # are correct just because they are the fallback.
    return Inf32
end

"""
    label_boundary_candidate_targets!(candidates; ...)

Add self-supervised training labels for candidate boundary ranking.

Targets are selected within each `boundary_group_id`. If enough related
precursors exist, candidates are ranked by how well their log area agrees with
the leave-one-group-out run shift of other precursors from the same protein.
When isotope traces are present, a second loss compares each candidate against
the fallback areas of other isotope traces from the same precursor and run.
Groups without either signal remain unlabeled; the second-derivative candidate
is still used as the runtime fallback when the model cannot be trained.
"""
function label_boundary_candidate_targets!(
    candidates::DataFrame;
    group_col::Symbol = :boundary_group_id,
    precursor_col::Symbol = :precursor_idx,
    protein_col::Symbol = :protein_key,
    file_col::Symbol = :ms_file_idx,
    area_col::Symbol = :peak_area,
    fallback_col::Symbol = :is_fallback,
    isotope_col::Union{Nothing, Symbol} = hasproperty(candidates, :isotope_key) ? :isotope_key : nothing,
    min_crossrun_refs::Integer = 2,
    min_isotope_refs::Integer = 1,
    crossrun_weight::Real = 1.0f0,
    isotope_weight::Real = 1.0f0,
)
    n = nrow(candidates)
    candidates[!, :boundary_consistency_loss] = fill(Inf32, n)
    candidates[!, :boundary_isotope_loss] = fill(Inf32, n)
    candidates[!, :boundary_training_loss] = fill(Inf32, n)
    candidates[!, :boundary_candidate_target] = fill(false, n)
    n == 0 && return candidates

    rows_by_group, fallback_row = _fallback_rows_by_group(candidates, group_col, fallback_col)
    entries, by_precursor, by_precursor_run_isotope = _fallback_log_entries(
        candidates,
        fallback_row,
        group_col,
        precursor_col,
        file_col,
        protein_col,
        area_col,
        isotope_col,
    )
    crossrun_residuals = _crossrun_reference_residuals(entries, by_precursor)
    fallback = _get_column_or_fill(candidates, fallback_col, false)

    for row in 1:n
        crossrun_loss = _crossrun_loss_for_candidate(
            row,
            candidates,
            by_precursor,
            crossrun_residuals,
            group_col,
            precursor_col,
            file_col,
            protein_col,
            area_col,
            min_crossrun_refs,
        )
        isotope_loss = _isotope_loss_for_candidate(
            row,
            candidates,
            by_precursor_run_isotope,
            group_col,
            precursor_col,
            file_col,
            area_col,
            isotope_col,
            min_isotope_refs,
        )
        candidates[row, :boundary_consistency_loss] = crossrun_loss
        candidates[row, :boundary_isotope_loss] = isotope_loss
        candidates[row, :boundary_training_loss] = _combined_boundary_loss(
            crossrun_loss,
            isotope_loss,
            Bool(fallback[row]),
            Float32(crossrun_weight),
            Float32(isotope_weight),
        )
    end

    for (_, rows) in rows_by_group
        best_row = rows[1]
        best_loss = candidates[best_row, :boundary_training_loss]
        for row in rows
            loss = candidates[row, :boundary_training_loss]
            if loss < best_loss ||
               (loss == best_loss && Bool(fallback[row]) && !Bool(fallback[best_row]))
                best_row = row
                best_loss = loss
            end
        end
        if isfinite(best_loss)
            candidates[best_row, :boundary_candidate_target] = true
        end
    end

    return candidates
end

function _available_boundary_features(candidates::AbstractDataFrame, features::Vector{Symbol})
    return [feature for feature in features if hasproperty(candidates, feature)]
end

function _boundary_feature_frame(candidates::AbstractDataFrame, features::Vector{Symbol})
    feature_data = DataFrame()
    for feature in features
        vals = Vector{Float32}(undef, nrow(candidates))
        col = candidates[!, feature]
        for i in eachindex(vals)
            value = Float32(coalesce(col[i], 0.0f0))
            vals[i] = isfinite(value) ? value : 0.0f0
        end
        feature_data[!, feature] = vals
    end
    return feature_data
end

function _boundary_relevance_from_losses(losses::AbstractVector)
    n = length(losses)
    finite_losses = sort(unique(Float32(loss) for loss in losses if isfinite(Float32(loss))))
    isempty(finite_losses) && return zeros(Int32, n)

    loss_rank = Dict{Float32, Int}()
    for (rank, loss) in enumerate(finite_losses)
        loss_rank[loss] = rank
    end

    relevance = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        loss = Float32(losses[i])
        relevance[i] = isfinite(loss) ?
            Int32(n - loss_rank[loss]) :
            Int32(0)
    end
    max_relevance = maximum(relevance)
    if max_relevance > BOUNDARY_RANKING_MAX_LABEL
        @inbounds for i in eachindex(relevance)
            raw = relevance[i]
            if raw > 0
                scaled = ceil(
                    Int,
                    Float64(raw) * Float64(BOUNDARY_RANKING_MAX_LABEL) / Float64(max_relevance),
                )
                relevance[i] = Int32(clamp(scaled, 1, BOUNDARY_RANKING_MAX_LABEL))
            end
        end
    end
    return relevance
end

function boundary_ranking_training_data(
    candidates::AbstractDataFrame,
    features::Vector{Symbol};
    group_col::Symbol = :boundary_group_id,
    loss_col::Symbol = :boundary_training_loss,
)
    nrow(candidates) == 0 && return DataFrame(), Int32[], Int[]

    groups = sort(unique(candidates[!, group_col]))
    sorted_rows = Int[]
    group_sizes = Int[]
    relevance = Int32[]
    sizehint!(sorted_rows, nrow(candidates))
    sizehint!(relevance, nrow(candidates))

    for group_id in groups
        rows = findall(==(group_id), candidates[!, group_col])
        isempty(rows) && continue
        any(isfinite(Float32(candidates[row, loss_col])) for row in rows) || continue
        append!(sorted_rows, rows)
        push!(group_sizes, length(rows))
        append!(relevance, _boundary_relevance_from_losses(candidates[rows, loss_col]))
    end

    training = candidates[sorted_rows, :]
    return _boundary_feature_frame(training, features), relevance, group_sizes
end

function _build_boundary_ranker(;
    num_iterations::Integer = 100,
    max_depth::Integer = 4,
    num_leaves::Integer = 15,
    learning_rate::Real = 0.05,
    feature_fraction::Real = 0.85,
    bagging_fraction::Real = 0.85,
    bagging_freq::Integer = 1,
    min_data_in_leaf::Integer = 5,
    min_gain_to_split::Real = 0.0,
    lambda_l2::Real = 1.0,
)
    return LightGBM.LGBMRanking(
        objective = "lambdarank",
        metric = ["ndcg"],
        eval_at = [1, 3, 5],
        learning_rate = float(learning_rate),
        num_iterations = Int(num_iterations),
        num_leaves = Int(num_leaves),
        max_depth = Int(max_depth),
        feature_fraction = float(feature_fraction),
        bagging_fraction = float(bagging_fraction),
        bagging_freq = Int(bagging_freq),
        min_data_in_leaf = Int(min_data_in_leaf),
        min_gain_to_split = float(min_gain_to_split),
        lambda_l2 = float(lambda_l2),
        num_threads = Threads.nthreads(),
        verbosity = -1,
        seed = 1776,
        deterministic = true,
        force_row_wise = true,
    )
end

"""
    train_boundary_candidate_model(candidates; ...)

Train a small LightGBM ranker to order boundary candidates. Returns `nothing`
when there is not enough relevance diversity, preserving the second-derivative
fallback.
"""
function train_boundary_candidate_model(
    candidates::DataFrame;
    features::Vector{Symbol} = BOUNDARY_CANDIDATE_FEATURES,
    max_groups::Integer = 20_000,
    min_positive::Integer = 25,
    min_negative::Integer = 25,
    rng::AbstractRNG = Random.default_rng(),
)
    nrow(candidates) == 0 && return nothing

    training = sample_boundary_candidate_groups(candidates; max_groups = max_groups, rng = rng)
    hasproperty(training, :boundary_candidate_target) ||
        label_boundary_candidate_targets!(training)

    available = _available_boundary_features(training, features)
    isempty(available) && return nothing

    feature_frame, relevance, group_sizes = boundary_ranking_training_data(training, available)
    isempty(relevance) && return nothing

    n_positive = count(>(0), relevance)
    n_negative = count(==(0), relevance)
    (n_positive < min_positive || n_negative < min_negative) && return nothing

    model = _build_boundary_ranker(
        num_iterations = 100,
        learning_rate = 0.05,
        max_depth = 4,
        num_leaves = 15,
        feature_fraction = 0.85,
        bagging_fraction = 0.85,
        bagging_freq = 1,
        min_data_in_leaf = max(5, min(100, length(relevance) ÷ 50)),
        min_gain_to_split = 0.0,
        lambda_l2 = 1.0,
    )
    LightGBM.fit!(model, feature_matrix(feature_frame, available), relevance; group = group_sizes, verbosity = -1)
    return LightGBMModel(model, available, nothing)
end

function boundary_cv_folds(candidates::AbstractDataFrame)
    return sort(unique(UInt8.(candidates[!, :cv_fold])))
end

function boundary_training_candidates_for_cv_fold(
    candidates::AbstractDataFrame,
    heldout_fold::UInt8,
)
    mask = UInt8.(candidates[!, :cv_fold]) .!= heldout_fold
    return candidates[mask, :]
end

function train_boundary_candidate_model_for_cv_fold(
    candidates::DataFrame,
    heldout_fold::UInt8;
    features::Vector{Symbol} = BOUNDARY_CANDIDATE_FEATURES,
    max_groups::Integer = 20_000,
    min_positive::Integer = 25,
    min_negative::Integer = 25,
    min_crossrun_refs::Integer = 2,
    rng::AbstractRNG = Random.default_rng(),
)
    training = DataFrame(boundary_training_candidates_for_cv_fold(candidates, heldout_fold))
    nrow(training) == 0 && return nothing
    if !hasproperty(training, :boundary_candidate_target) ||
       !hasproperty(training, :boundary_training_loss)
        label_boundary_candidate_targets!(
            training;
            min_crossrun_refs = min_crossrun_refs,
        )
    end

    return train_boundary_candidate_model(
        training;
        features = features,
        max_groups = max_groups,
        min_positive = min_positive,
        min_negative = min_negative,
        rng = rng,
    )
end

function train_boundary_candidate_models_by_cv_fold(
    candidates::DataFrame;
    features::Vector{Symbol} = BOUNDARY_CANDIDATE_FEATURES,
    max_groups::Integer = 20_000,
    min_positive::Integer = 25,
    min_negative::Integer = 25,
    min_crossrun_refs::Integer = 2,
    rng::AbstractRNG = Random.default_rng(),
)
    models = Dict{UInt8, Union{LightGBMModel, Nothing}}()
    for heldout_fold in boundary_cv_folds(candidates)
        models[heldout_fold] = train_boundary_candidate_model_for_cv_fold(
            candidates,
            heldout_fold;
            features = features,
            max_groups = max_groups,
            min_positive = min_positive,
            min_negative = min_negative,
            min_crossrun_refs = min_crossrun_refs,
            rng = rng,
        )
    end
    return models
end

function score_boundary_candidates!(
    candidates::DataFrame,
    model::Union{LightGBMModel, Nothing};
    score_col::Symbol = :boundary_model_score,
)
    if model === nothing
        fallback = _get_column_or_fill(candidates, :is_fallback, false)
        candidates[!, score_col] = Float32[Bool(x) ? 1.0f0 : 0.0f0 for x in fallback]
        return candidates
    end

    available = [feature for feature in model.features if hasproperty(candidates, feature)]
    if length(available) != length(model.features)
        fallback = _get_column_or_fill(candidates, :is_fallback, false)
        candidates[!, score_col] = Float32[Bool(x) ? 1.0f0 : 0.0f0 for x in fallback]
        return candidates
    end

    candidates[!, score_col] = Float32.(predict(model, _boundary_feature_frame(candidates, model.features)))
    return candidates
end

function score_boundary_candidates_crossfold!(
    candidates::DataFrame,
    models::AbstractDict{UInt8, <:Union{LightGBMModel, Nothing}};
    score_col::Symbol = :boundary_model_score,
)
    candidates[!, score_col] = zeros(Float32, nrow(candidates))
    for cv_fold in boundary_cv_folds(candidates)
        rows = findall(==(cv_fold), UInt8.(candidates[!, :cv_fold]))
        isempty(rows) && continue
        fold_candidates = candidates[rows, :]
        score_boundary_candidates!(
            fold_candidates,
            get(models, cv_fold, nothing);
            score_col = score_col,
        )
        candidates[rows, score_col] = fold_candidates[!, score_col]
    end
    return candidates
end

function boundary_model_feature_importance_lines(importances)
    sorted_imp = sort(collect(importances), by = x -> -x[2])
    lines = ["Learned chromatogram boundary LGBM feature gains (all $(length(sorted_imp))):"]
    for (fname, gain) in sorted_imp
        push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
    end
    return lines
end

function log_boundary_model_feature_importances(importances)
    lines = boundary_model_feature_importance_lines(importances)
    @debug_l1 join(lines, "\n")
    return nothing
end

function log_boundary_model_feature_importances(model::LightGBMModel)
    imp = importance(model)
    imp === nothing && return nothing
    log_boundary_model_feature_importances(imp)
    return nothing
end

function log_boundary_cv_model_feature_importances(
    models::AbstractDict{UInt8, <:Union{LightGBMModel, Nothing}},
)
    for fold in sort(collect(keys(models)))
        model = models[fold]
        if model === nothing
            @debug_l1 "Learned chromatogram boundary LGBM feature gains fold=$(fold): skipped=true reason=insufficient_training_data"
            continue
        end
        imp = importance(model)
        imp === nothing && continue
        lines = boundary_model_feature_importance_lines(imp)
        lines[1] = lines[1] * " fold=$(fold)"
        @debug_l1 join(lines, "\n")
    end
    return nothing
end

function boundary_candidate_category_counts(selected_candidates::AbstractDataFrame)
    counts = Dict{String, Int}(label => 0 for label in BOUNDARY_CANDIDATE_CATEGORY_LABELS)
    hasproperty(selected_candidates, :candidate_category) || return counts

    for category in selected_candidates[!, :candidate_category]
        label = String(category)
        counts[label] = get(counts, label, 0) + 1
    end
    return counts
end

function boundary_candidate_category_tally_lines(selected_candidates::AbstractDataFrame)
    total = nrow(selected_candidates)
    counts = boundary_candidate_category_counts(selected_candidates)
    lines = ["Learned chromatogram boundary selected candidate categories (total $(total)):"]

    for label in BOUNDARY_CANDIDATE_CATEGORY_LABELS
        push!(lines, "    $(rpad(label, 24)) $(counts[label])")
    end

    extra_labels = sort(setdiff(collect(keys(counts)), collect(BOUNDARY_CANDIDATE_CATEGORY_LABELS)))
    for label in extra_labels
        push!(lines, "    $(rpad(label, 24)) $(counts[label])")
    end

    return lines
end

function log_boundary_candidate_category_tally(selected_candidates::AbstractDataFrame)
    lines = boundary_candidate_category_tally_lines(selected_candidates)
    @debug_l1 join(lines, "\n")
    return nothing
end

function quant_boundary_candidate_rows(candidates::AbstractDataFrame)
    hasproperty(candidates, :quant_trace_selected) || return copy(candidates)
    mask = Bool.(candidates[!, :quant_trace_selected])
    return candidates[mask, :]
end

function _boundary_debug_value(df::AbstractDataFrame, row::Integer, col::Symbol)
    hasproperty(df, col) || return "NA"
    value = df[Int(row), col]
    if value isa AbstractFloat
        isfinite(value) || return string(value)
        return string(round(Float64(value), sigdigits = 5))
    end
    return string(value)
end

function _boundary_debug_selected_keys(selected_candidates::AbstractDataFrame)
    keys = Set{Tuple{UInt64, UInt16}}()
    if !hasproperty(selected_candidates, :boundary_group_id) ||
       !hasproperty(selected_candidates, :candidate_index)
        return keys
    end

    for row in eachrow(selected_candidates)
        push!(keys, (UInt64(row.boundary_group_id), UInt16(row.candidate_index)))
    end
    return keys
end

function boundary_candidate_debug_lines(
    candidates::AbstractDataFrame,
    selected_candidates::AbstractDataFrame;
    target_precursor_idx::UInt32 = UInt32(370714),
    target_ms_file_idx::UInt16 = UInt16(1),
)
    hasproperty(candidates, :precursor_idx) ||
        return ["Boundary candidate debug precursor_idx=$(target_precursor_idx) ms_file_idx=$(target_ms_file_idx) rows=0 reason=missing_precursor_idx"]
    hasproperty(candidates, :ms_file_idx) ||
        return ["Boundary candidate debug precursor_idx=$(target_precursor_idx) ms_file_idx=$(target_ms_file_idx) rows=0 reason=missing_ms_file_idx"]

    rows = findall(
        i -> UInt32(candidates[i, :precursor_idx]) == target_precursor_idx &&
             UInt16(candidates[i, :ms_file_idx]) == target_ms_file_idx,
        axes(candidates, 1),
    )
    if hasproperty(candidates, :boundary_group_id) && hasproperty(candidates, :candidate_index)
        sort!(rows, by = i -> (UInt64(candidates[i, :boundary_group_id]), UInt16(candidates[i, :candidate_index])))
    end

    selected_keys = _boundary_debug_selected_keys(selected_candidates)
    lines = [
        "Boundary candidate debug precursor_idx=$(target_precursor_idx) " *
        "ms_file_idx=$(target_ms_file_idx) rows=$(length(rows))",
    ]

    for row in rows
        group_id = hasproperty(candidates, :boundary_group_id) ?
            UInt64(candidates[row, :boundary_group_id]) : UInt64(0)
        candidate_index = hasproperty(candidates, :candidate_index) ?
            UInt16(candidates[row, :candidate_index]) : UInt16(0)
        selected = (group_id, candidate_index) in selected_keys
        range_idx = _boundary_debug_value(candidates, row, :candidate_start_idx) *
            ":" * _boundary_debug_value(candidates, row, :candidate_stop_idx)
        range_scan = _boundary_debug_value(candidates, row, :candidate_start_scan) *
            ":" * _boundary_debug_value(candidates, row, :candidate_stop_scan)

        push!(lines, join((
            "  group=$(group_id)",
            "candidate_index=$(candidate_index)",
            "selected=$(selected)",
            "category=$(_boundary_debug_value(candidates, row, :candidate_category))",
            "range_idx=$(range_idx)",
            "range_scan=$(range_scan)",
            "area=$(_boundary_debug_value(candidates, row, :peak_area))",
            "points=$(_boundary_debug_value(candidates, row, :points_integrated))",
            "score=$(_boundary_debug_value(candidates, row, :boundary_model_score))",
            "loss=$(_boundary_debug_value(candidates, row, :boundary_training_loss))",
            "crossrun_loss=$(_boundary_debug_value(candidates, row, :boundary_crossrun_loss))",
            "isotope_loss=$(_boundary_debug_value(candidates, row, :boundary_isotope_loss))",
            "fallback=$(_boundary_debug_value(candidates, row, :is_fallback))",
            "quant_trace_selected=$(_boundary_debug_value(candidates, row, :quant_trace_selected))",
            "isotope_key=$(_boundary_debug_value(candidates, row, :isotope_key))",
            "width=$(_boundary_debug_value(candidates, row, :candidate_width))",
            "endpoint_valley=$(_boundary_debug_value(candidates, row, :endpoint_valley_score))",
            "endpoint_height=$(_boundary_debug_value(candidates, row, :endpoint_height_fraction))",
            "secondary_peak=$(_boundary_debug_value(candidates, row, :secondary_peak_penalty))",
            "irt_asymmetry=$(_boundary_debug_value(candidates, row, :irt_asymmetry_delta))",
            "baseline_disconnect=$(_boundary_debug_value(candidates, row, :baseline_disconnected_signal_fraction))",
            "baseline_nonapex_lobe=$(_boundary_debug_value(candidates, row, :baseline_largest_nonapex_lobe_log2_ratio))",
            "baseline_trough=$(_boundary_debug_value(candidates, row, :baseline_internal_trough_score)):$(_boundary_debug_value(candidates, row, :baseline_internal_trough_log2_ratio))",
            "dip_recovery=$(_boundary_debug_value(candidates, row, :internal_dip_recovery_score))",
            "edge_valley_ratio=$(_boundary_debug_value(candidates, row, :left_edge_valley_log2_ratio)):$(_boundary_debug_value(candidates, row, :right_edge_valley_log2_ratio))",
            "included_nonapex=$(_boundary_debug_value(candidates, row, :included_nonapex_max_log2_ratio))",
            "included_nonapex_irt=$(_boundary_debug_value(candidates, row, :included_nonapex_max_irt_distance))",
            "right_excluded=$(_boundary_debug_value(candidates, row, :right_excluded_signal_fraction))",
            "right_recovery=$(_boundary_debug_value(candidates, row, :right_boundary_recovery_fraction))",
            "right_outside_peak=$(_boundary_debug_value(candidates, row, :right_outside_peak_fraction))",
            "mainsearch_delta=$(_boundary_debug_value(candidates, row, :mainsearch_left_bound_delta)):$(_boundary_debug_value(candidates, row, :mainsearch_right_bound_delta))",
        ), " "))
    end

    return lines
end

function log_boundary_candidate_debug(
    candidates::AbstractDataFrame,
    selected_candidates::AbstractDataFrame;
    target_precursor_idx::Union{Nothing,UInt32} = nothing,
    target_ms_file_idx::Union{Nothing,UInt16} = nothing,
)
    DEBUG_CONSOLE_LEVEL[] >= 1 || return nothing
    targets = if target_precursor_idx !== nothing && target_ms_file_idx !== nothing
        [(target_precursor_idx, target_ms_file_idx)]
    else
        sort!(collect(DEBUG_BOUNDARY_CANDIDATE_TARGETS[]))
    end

    for (precursor_idx, ms_file_idx) in targets
        lines = boundary_candidate_debug_lines(
            candidates,
            selected_candidates;
            target_precursor_idx = precursor_idx,
            target_ms_file_idx = ms_file_idx,
        )
        @debug_l1 join(lines, "\n")
    end
    return nothing
end

function select_boundary_candidate_rows(
    candidates::DataFrame,
    model::Union{LightGBMModel, Nothing};
    group_col::Symbol = :boundary_group_id,
    score_col::Symbol = :boundary_model_score,
)
    selection_candidates = quant_boundary_candidate_rows(candidates)
    nrow(selection_candidates) == 0 && return copy(selection_candidates)
    scored = copy(selection_candidates)
    score_boundary_candidates!(scored, model; score_col = score_col)

    rows_by_group, _ = _fallback_rows_by_group(scored, group_col, :is_fallback)
    fallback = _get_column_or_fill(scored, :is_fallback, false)
    selected_rows = Int[]
    sizehint!(selected_rows, length(rows_by_group))

    for (_, rows) in rows_by_group
        best_row = rows[1]
        best_score = scored[best_row, score_col]
        for row in rows
            score = scored[row, score_col]
            if score > best_score ||
               (score == best_score && Bool(fallback[row]) && !Bool(fallback[best_row]))
                best_row = row
                best_score = score
            end
        end
        push!(selected_rows, best_row)
    end

    sort!(selected_rows)
    return scored[selected_rows, :]
end

function select_boundary_candidate_rows_crossfold(
    candidates::DataFrame,
    models::AbstractDict{UInt8, <:Union{LightGBMModel, Nothing}};
)
    selection_candidates = quant_boundary_candidate_rows(candidates)
    nrow(selection_candidates) == 0 && return copy(selection_candidates)

    selected_by_fold = DataFrame[]
    for cv_fold in boundary_cv_folds(selection_candidates)
        rows = findall(==(cv_fold), UInt8.(selection_candidates[!, :cv_fold]))
        isempty(rows) && continue
        fold_candidates = DataFrame(selection_candidates[rows, :])
        push!(
            selected_by_fold,
            select_boundary_candidate_rows(
                fold_candidates,
                get(models, cv_fold, nothing),
            ),
        )
    end

    isempty(selected_by_fold) && return selection_candidates[1:0, :]
    return reduce(vcat, selected_by_fold; cols = :union)
end

function append_boundary_candidate_rows!(
    rows::Vector{NamedTuple},
    candidates,
    boundary_group_id::UInt64,
    ms_file_idx::Integer,
    precursor_idx::UInt32,
    cv_fold::UInt8,
    protein_key,
    sequence_key,
    isotope_key,
    target::Bool,
    quant_trace_selected::Bool,
)
    candidates === nothing && return rows

    for candidate in candidates
        push!(rows, (;
            boundary_group_id = boundary_group_id,
            ms_file_idx = UInt16(ms_file_idx),
            precursor_idx = precursor_idx,
            cv_fold = cv_fold,
            protein_key = protein_key,
            sequence_key = sequence_key,
            isotope_key = isotope_key,
            target = target,
            quant_trace_selected = quant_trace_selected,
            candidate...
        ))
    end

    return rows
end

function boundary_candidate_table(
    boundary_candidate_data::AbstractVector,
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Integer;
    extra_rows::Vector{NamedTuple} = NamedTuple[],
)
    precursors = getPrecursors(getSpecLib(search_context))
    has_isotope_key = hasproperty(psms, :isotopes_captured)
    has_target = hasproperty(psms, :target)
    rows = copy(extra_rows)

    for i in eachindex(boundary_candidate_data)
        candidates = boundary_candidate_data[i]
        candidates === nothing && continue

        pid = UInt32(psms[i, :precursor_idx])
        boundary_group_id = (UInt64(ms_file_idx) << 32) | UInt64(i)
        protein_key = getAccessionNumbers(precursors)[pid]
        sequence_key = getSequence(precursors)[pid]
        cv_fold = getCvFold(precursors, pid)
        isotope_key = has_isotope_key ? psms[i, :isotopes_captured] : (Int8(0), Int8(0))
        target = has_target ? Bool(psms[i, :target]) : true

        append_boundary_candidate_rows!(
            rows,
            candidates,
            boundary_group_id,
            ms_file_idx,
            pid,
            cv_fold,
            protein_key,
            sequence_key,
            isotope_key,
            target,
            true,
        )
    end

    return DataFrame(rows)
end

function boundary_candidate_dir(search_context::SearchContext)
    return joinpath(getDataOutDir(search_context), "temp_data", "boundary_candidates")
end

function boundary_candidate_path(search_context::SearchContext, ms_file_idx::Integer)
    return joinpath(boundary_candidate_dir(search_context), "file_$(Int(ms_file_idx)).arrow")
end

function write_boundary_candidate_table!(
    boundary_candidate_data::AbstractVector,
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Integer,
    extra_rows::Vector{NamedTuple} = NamedTuple[],
)
    candidates = boundary_candidate_table(
        boundary_candidate_data,
        psms,
        search_context,
        ms_file_idx;
        extra_rows = extra_rows,
    )
    nrow(candidates) == 0 && return nothing

    out_dir = boundary_candidate_dir(search_context)
    isdir(out_dir) || mkpath(out_dir)
    writeArrow(boundary_candidate_path(search_context, ms_file_idx), candidates)
    return nothing
end

function collect_isotope_boundary_candidate_rows(
    chromatograms::DataFrame,
    passing_psms::DataFrame,
    selected_quant_trace::Union{Nothing, Dict{UInt32, Tuple{Tuple{Int8, Int8}, Float32}}},
    search_context::SearchContext,
    ms_file_idx::Integer,
    min_fraction_transmitted::Float32,
    λ::Float32;
    max_groups::Integer = typemax(Int),
    rng::AbstractRNG = Random.default_rng(),
)
    nrow(chromatograms) == 0 && return NamedTuple[]
    hasproperty(chromatograms, :isotopes_captured) || return NamedTuple[]

    # The quant path may be combined by precursor only. For learner-only isotope
    # traces, regroup the already extracted points by captured isotope set after
    # the fallback quantification has finished.
    sort_chromatograms_for_integration!(chromatograms, SeperateTraces())

    precursors = getPrecursors(getSpecLib(search_context))
    rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
    chrom_index, max_chrom_len = build_chrom_index(chromatograms, SeperateTraces())
    max_chrom_len <= 0 && return NamedTuple[]

    rt_all = chromatograms[!, :rt]::AbstractVector{Float32}
    scan_idx_all = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    intensity_all = chromatograms[!, :intensity]::AbstractVector{Float32}
    fraction_all = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}

    psm_row_by_pid = Dict{UInt32, Int}()
    for i in 1:nrow(passing_psms)
        psm_row_by_pid[UInt32(passing_psms[i, :precursor_idx])] = i
    end

    ws = WHWorkspace(max_chrom_len)
    state = Chromatogram(zeros(Float32, max_chrom_len), zeros(Float32, max_chrom_len), 0)
    rows = NamedTuple[]
    extra_group_idx = 0
    has_target = hasproperty(passing_psms, :target)

    chrom_keys = collect(keys(chrom_index))
    Random.shuffle!(rng, chrom_keys)
    max_iter = min(length(chrom_keys), max(0, Int(max_groups)))

    for chrom_key in @view chrom_keys[1:max_iter]
        chrom_range = chrom_index[chrom_key]
        pid, isotope_key = chrom_key
        psm_row = get(psm_row_by_pid, pid, 0)
        psm_row == 0 && continue

        selected = selected_quant_trace === nothing ? nothing : get(selected_quant_trace, pid, nothing)
        selected !== nothing && isotope_key == selected[1] && continue
        isnothing(findfirst(x -> x > 0.0f0, @view intensity_all[chrom_range])) && continue

        avg_cycle_time = (rt_all[last(chrom_range)] - rt_all[first(chrom_range)]) / length(chrom_range)
        apex_scan = find_nearest_scan(
            @view(scan_idx_all[chrom_range]),
            UInt32(passing_psms[psm_row, :scan_idx]),
        )

        candidate_ref = Ref{Any}(nothing)
        integrate_chrom(
            @view(rt_all[chrom_range]),
            @view(scan_idx_all[chrom_range]),
            @view(intensity_all[chrom_range]),
            @view(fraction_all[chrom_range]),
            apex_scan,
            ws,
            state,
            avg_cycle_time,
            λ;
            min_fraction_transmitted = min_fraction_transmitted,
            mainsearch_1pct_start_scan = hasproperty(passing_psms, :mainsearch_1pct_start_scan) ?
                UInt32(passing_psms[psm_row, :mainsearch_1pct_start_scan]) :
                UInt32(0),
            mainsearch_1pct_stop_scan = hasproperty(passing_psms, :mainsearch_1pct_stop_scan) ?
                UInt32(passing_psms[psm_row, :mainsearch_1pct_stop_scan]) :
                UInt32(0),
            rt_to_irt_model = rt_to_irt_model,
            boundary_candidate_data = candidate_ref,
        )
        reset!(state)
        candidate_ref[] === nothing && continue

        extra_group_idx += 1
        boundary_group_id = (UInt64(ms_file_idx) << 32) |
            UInt64(nrow(passing_psms) + extra_group_idx)
        target = has_target ? Bool(passing_psms[psm_row, :target]) : true
        append_boundary_candidate_rows!(
            rows,
            candidate_ref[],
            boundary_group_id,
            ms_file_idx,
            pid,
            getCvFold(precursors, pid),
            getAccessionNumbers(precursors)[pid],
            getSequence(precursors)[pid],
            isotope_key,
            target,
            false,
        )
    end

    return rows
end

function load_boundary_candidate_tables(search_context::SearchContext)
    candidate_dir = boundary_candidate_dir(search_context)
    isdir(candidate_dir) || return DataFrame()

    paths = filter(path -> endswith(path, ".arrow"), readdir(candidate_dir; join = true))
    isempty(paths) && return DataFrame()

    tables = DataFrame[]
    sizehint!(tables, length(paths))
    for path in paths
        push!(tables, DataFrame(Tables.columntable(Arrow.Table(path))))
    end
    isempty(tables) && return DataFrame()
    return reduce(vcat, tables; cols = :union)
end

function _candidate_lookup_key(precursor_idx, isotope_key)
    return (UInt32(precursor_idx), isotope_key)
end

function apply_selected_boundary_candidates!(
    selected_candidates::DataFrame,
    search_context::SearchContext,
)
    nrow(selected_candidates) == 0 && return 0

    updated = 0
    ms_data = getMSData(search_context)
    for group in groupby(selected_candidates, :ms_file_idx)
        ms_file_idx = Int(first(group.ms_file_idx))
        psm_path = getPassingPsms(ms_data)[ms_file_idx]
        (isempty(psm_path) || !isfile(psm_path)) && continue

        lookup = Dict{Tuple{UInt32, Any}, Tuple{Float32, UInt32, UInt32}}()
        for row in eachrow(group)
            area = Float32(row.peak_area)
            area > 0.0f0 && isfinite(area) || continue
            lookup[_candidate_lookup_key(row.precursor_idx, row.isotope_key)] = (
                area,
                UInt32(row.new_best_scan),
                UInt32(row.points_integrated),
            )
        end
        isempty(lookup) && continue

        psms = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
        nrow(psms) == 0 && continue
        has_isotope_key = hasproperty(psms, :isotopes_captured)

        for i in 1:nrow(psms)
            isotope_key = has_isotope_key ? psms[i, :isotopes_captured] : (Int8(0), Int8(0))
            value = get(lookup, _candidate_lookup_key(psms[i, :precursor_idx], isotope_key), nothing)
            value === nothing && continue

            psms[i, :peak_area] = value[1]
            psms[i, :new_best_scan] = value[2]
            psms[i, :points_integrated] = value[3]
            updated += 1
        end

        writeArrow(psm_path, psms)
    end

    return updated
end
