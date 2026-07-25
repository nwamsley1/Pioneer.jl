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
    log_probit_feature_importance(feature_names, β_fitted, X; context = "protein_probit")

Log fitted protein probit coefficients together with one-sigma effects.
"""
function log_probit_feature_importance(
    feature_names::Vector{Symbol},
    β_fitted::AbstractVector{<:Real},
    X::Matrix{Float64};
    context::String = "protein_probit"
)
    n_coefficients = max(length(β_fitted) - 1, 0)
    n_features = min(length(feature_names), size(X, 2), n_coefficients)

    if n_features == 0
        @user_warn "Protein probit fit produced no feature coefficients to log"
        return
    end

    intercept = Float64(β_fitted[1])
    @debug_l1 "Protein probit model coefficients context=$(context) intercept=$(round(intercept, digits=6)) n_features=$(n_features) importance_metric=abs(coefficient*feature_std)"

    feature_summaries = NamedTuple[]
    for j in 1:n_features
        feature_std = Float64(std(view(X, :, j)))
        coefficient = Float64(β_fitted[j + 1])
        one_sigma_effect = coefficient * feature_std
        push!(
            feature_summaries,
            (
                feature = String(feature_names[j]),
                coefficient = coefficient,
                abs_coefficient = abs(coefficient),
                feature_std = feature_std,
                one_sigma_effect = one_sigma_effect,
                abs_one_sigma_effect = abs(one_sigma_effect)
            )
        )
    end

    sort!(feature_summaries, by = x -> x.abs_one_sigma_effect, rev = true)

    for (rank, summary) in enumerate(feature_summaries)
        @debug_l1 "Protein probit feature importance context=$(context) rank=$(rank) feature=$(summary.feature) coefficient=$(round(summary.coefficient, digits=6)) feature_std=$(round(summary.feature_std, digits=6)) one_sigma_effect=$(round(summary.one_sigma_effect, digits=6)) abs_one_sigma_effect=$(round(summary.abs_one_sigma_effect, digits=6))"
    end

    return
end

function make_protein_feature_importance_logger(
    feature_names::Vector{Symbol};
    context::AbstractString
)
    return function(β_fitted::AbstractVector{<:Real}, X::Matrix{Float64}, iteration::Int)
        log_probit_feature_importance(
            feature_names,
            β_fitted,
            X;
            context = context * "_iter_$(iteration)"
        )
        return nothing
    end
end

"""
    build_protein_semisupervised_training_set(scores, targets, prefix_shape, n_non_mbr_peptides; q_value_threshold = 0.01f0, max_positive_pep_threshold = 1.0f0, mined_negative_prefix_shape_threshold = -0.20f0, mined_negative_pep_threshold = 0.90f0, keep_non_mined_targets_as_positive = true)

Build labels for semi-supervised protein probit training from a score vector.
Targets are mined as negatives when `PEP >= mined_negative_pep_threshold` or
when a target with at most one non-MBR unique peptide has
`prefix_shape <= mined_negative_prefix_shape_threshold`. Shared-peptide scores,
counts, and prefix shapes do not affect negative mining.
If `keep_non_mined_targets_as_positive=true`, remaining targets stay positive
for training; otherwise only q-value-passing targets stay positive and the rest
are dropped.
"""
function build_protein_semisupervised_training_set(
    scores::AbstractVector{<:Real},
    targets::AbstractVector{Bool},
    prefix_shape::AbstractVector{<:Real},
    n_non_mbr_peptides::AbstractVector{<:Integer};
    q_value_threshold::Float32 = 0.01f0,
    max_positive_pep_threshold::Float32 = 1.0f0,
    mined_negative_prefix_shape_threshold::Float32 = -0.20f0,
    mined_negative_pep_threshold::Float32 = 0.90f0,
    keep_non_mined_targets_as_positive::Bool = true
)
    n = length(scores)
    length(targets) == n || throw(ArgumentError("targets must have the same length as scores"))
    length(prefix_shape) == n || throw(ArgumentError("prefix_shape must have the same length as scores"))
    length(n_non_mbr_peptides) == n ||
        throw(ArgumentError("n_non_mbr_peptides must have the same length as scores"))
    qvals = Vector{Float32}(undef, n)
    peps = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals)
    get_PEP!(scores, targets, peps)
    
    confident_positive_mask = BitVector(undef, n)
    candidate_confident_positive_mask = BitVector(undef, n)
    mined_negative_mask = BitVector(undef, n)
    positive_mask = BitVector(undef, n)
    keep_mask = BitVector(undef, n)

    @inbounds for i in eachindex(
        scores,
        targets,
        prefix_shape,
        n_non_mbr_peptides,
        qvals,
        peps,
    )
        low_shape_singleton_target = targets[i] &&
                                     (n_non_mbr_peptides[i] <= 1) &&
                                     (Float32(prefix_shape[i]) <= mined_negative_prefix_shape_threshold)
        candidate_confident_positive_mask[i] = targets[i] &&
                                               (qvals[i] <= q_value_threshold) &&
                                               (peps[i] <= max_positive_pep_threshold)
        mined_negative_mask[i] = targets[i] &&
                                 ((peps[i] >= mined_negative_pep_threshold) ||
                                  low_shape_singleton_target)
        confident_positive_mask[i] = candidate_confident_positive_mask[i] &&
                                     !mined_negative_mask[i]
        positive_mask[i] = targets[i] &&
                           !mined_negative_mask[i] &&
                           (keep_non_mined_targets_as_positive || confident_positive_mask[i])
        keep_mask[i] = !targets[i] || positive_mask[i] || mined_negative_mask[i]
    end

    second_pass_labels = Vector{Bool}(positive_mask[keep_mask])

    return (
        keep_mask = keep_mask,
        labels = second_pass_labels,
        positive_mask = positive_mask,
        confident_positive_mask = confident_positive_mask,
        mined_negative_mask = mined_negative_mask,
        requested_mined_negative_count = sum(mined_negative_mask),
        mined_negative_prefix_shape_threshold = mined_negative_prefix_shape_threshold,
        mined_negative_pep_threshold = mined_negative_pep_threshold,
        qvals = qvals,
        peps = peps
    )
end

@inline function _count_mask_changes(
    current_mask::AbstractVector{Bool},
    previous_mask::AbstractVector{Bool}
)::Int
    change_count = 0
    @inbounds for i in eachindex(current_mask, previous_mask)
        change_count += current_mask[i] != previous_mask[i]
    end
    return change_count
end

@inline function _count_training_label_state_changes(
    keep_mask::AbstractVector{Bool},
    positive_mask::AbstractVector{Bool},
    previous_keep_mask::AbstractVector{Bool},
    previous_positive_mask::AbstractVector{Bool}
)::Int
    change_count = 0
    @inbounds for i in eachindex(keep_mask, positive_mask, previous_keep_mask, previous_positive_mask)
        current_state = positive_mask[i] ? Int8(1) : (keep_mask[i] ? Int8(-1) : Int8(0))
        previous_state = previous_positive_mask[i] ? Int8(1) : (previous_keep_mask[i] ? Int8(-1) : Int8(0))
        change_count += current_state != previous_state
    end
    return change_count
end


"""
    fit_probit_model(X::Matrix{Float64}, y::Vector{Bool})

Fit a probit regression model for protein group classification.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `y::Vector{Bool}`: Target labels (true for targets, false for decoys)

# Returns
- `β_fitted`: Fitted coefficients
- `X_mean`: Feature means for standardization
- `X_std`: Feature standard deviations for standardization
"""
function fit_probit_model(X::Matrix{Float64}, y::AbstractVector{Bool})
    size(X, 2) > 0 || throw(ArgumentError("No valid features remaining for probit regression after variance filtering"))

    # Standardize features
    #X_mean = mean(X, dims=1)
    #X_std = std(X, dims=1)
    #X_std[X_std .== 0] .= 1.0  # Should not happen after filtering, but defensive
    X_standardized = X #(X .- X_mean) ./ X_std

    # Add intercept column
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])

    # Initialize coefficients
    β = zeros(Float64, size(X_with_intercept, 2))

    # Create data chunks for parallel processing
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, length(y) / n_chunks))
    data_chunks = Iterators.partition(1:length(y), chunk_size)

    # Fit probit model
    β_fitted = Pioneer.ProbitRegression(β, X_df, y, data_chunks, max_iter=30)

    return β_fitted#, vec(X_mean), vec(X_std)
end

"""
    fit_probit_model_semisupervised(X, y, initial_scores, prefix_shape, n_non_mbr_peptides; q_value_threshold = 0.01f0, min_prefix_shape_neg_threshold = -0.20f0, min_pep_neg_threshold = 0.90f0, max_positive_pep_threshold = 1.0f0, n_iterations = 10, context = "protein_probit", iteration_debug_callback = nothing)

Fit the protein probit model by seeding iteration 1 labels from raw initial scores,
then fitting and refining with the full feature set. Later outer iterations stop
early when label churn stays small for two consecutive rounds.
"""
function fit_probit_model_semisupervised(
    X::Matrix{Float64},
    y::AbstractVector{Bool},
    initial_scores::AbstractVector{<:Real},
    prefix_shape::AbstractVector{<:Real},
    n_non_mbr_peptides::AbstractVector{<:Integer};
    q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold::Float32 = -0.20f0,
    min_pep_neg_threshold::Float32 = 0.90f0,
    max_positive_pep_threshold::Float32 = 1.0f0,
    n_iterations::Int = 10,
    context::AbstractString = "protein_probit",
    iteration_debug_callback = nothing,
    feature_importance_logger = nothing
)
    size(X, 1) == length(y) || throw(ArgumentError("X must have the same number of rows as y"))
    length(initial_scores) == length(y) || throw(ArgumentError("initial_scores must have the same length as y"))
    length(prefix_shape) == length(y) || throw(ArgumentError("prefix_shape must have the same length as y"))
    length(n_non_mbr_peptides) == length(y) ||
        throw(ArgumentError("n_non_mbr_peptides must have the same length as y"))

    row_change_threshold = max(5, ceil(Int, 0.001 * length(y)))
    near_convergence_check_start_iteration = 4
    near_convergence_patience = 2
    near_convergence_streak = 0
    last_plot_iteration = 1

    # Iteration 1 is seeded from the raw protein-group score. At this stage we
    # keep non-mined target rows as positives so the model starts from a broad,
    # conservative training set rather than an aggressively filtered one.
    ss_initial = build_protein_semisupervised_training_set(
        initial_scores,
        y,
        prefix_shape,
        n_non_mbr_peptides;
        q_value_threshold = q_value_threshold,
        max_positive_pep_threshold = max_positive_pep_threshold,
        mined_negative_prefix_shape_threshold = min_prefix_shape_neg_threshold,
        mined_negative_pep_threshold = min_pep_neg_threshold,
        keep_non_mined_targets_as_positive = true
    )
    initial_keep_mask = ss_initial.keep_mask
    initial_labels = ss_initial.labels
    last_plot_state = (
        positive_mask = ss_initial.positive_mask,
        confident_positive_mask = ss_initial.confident_positive_mask,
        mined_negative_mask = ss_initial.mined_negative_mask
    )

    initial_target_count = sum(initial_labels)
    initial_decoy_count = sum(initial_keep_mask) - initial_target_count
    if initial_target_count < 10 || initial_decoy_count < 10
        @user_warn "Protein probit semi-supervised iteration has insufficient data; falling back to all rows for initial pg_score fit" context = context iteration = 1 selected_targets = initial_target_count selected_decoys = initial_decoy_count
        initial_keep_mask = trues(length(y))
        initial_labels = y
        last_plot_state = (
            positive_mask = copy(y),
            confident_positive_mask = ss_initial.confident_positive_mask,
            mined_negative_mask = falses(length(y))
        )
    end

    # Fit the first model on the seeded labels.
    X_initial_iteration = X[initial_keep_mask, :]
    β_current = fit_probit_model(X_initial_iteration, initial_labels)
    if !isnothing(feature_importance_logger)
        feature_importance_logger(β_current, X_initial_iteration, 1)
    end
    if !isnothing(iteration_debug_callback)
        iteration_debug_callback(1, last_plot_state)
    end

    if n_iterations <= 1
        return β_current
    end

    previous_keep_mask = copy(initial_keep_mask)
    previous_positive_mask = copy(last_plot_state.positive_mask)

    for iteration in 2:n_iterations
        # Score all rows with the current model, then rebuild the training set
        # using only q-passing targets as positives. This is the semi-supervised
        # refinement step.
        iteration_scores = calculate_probit_scores(X, β_current)
        ss = build_protein_semisupervised_training_set(
            iteration_scores,
            y,
            prefix_shape,
            n_non_mbr_peptides;
            q_value_threshold = q_value_threshold,
            max_positive_pep_threshold = max_positive_pep_threshold,
            mined_negative_prefix_shape_threshold = min_prefix_shape_neg_threshold,
            mined_negative_pep_threshold = min_pep_neg_threshold,
            keep_non_mined_targets_as_positive = false
        )

        iteration_target_count = sum(ss.labels)
        iteration_decoy_count = sum(ss.keep_mask) - iteration_target_count
        if iteration_target_count < 10 || iteration_decoy_count < 10
            @user_warn "Protein probit semi-supervised iteration has insufficient data; using previous iteration model" context = context iteration = iteration selected_targets = iteration_target_count selected_decoys = iteration_decoy_count
            if !isnothing(iteration_debug_callback) && last_plot_iteration > 1
                iteration_debug_callback(last_plot_iteration, last_plot_state)
            end
            return β_current
        end

        # Early stopping is based on label churn, not score drift. Once the set
        # of kept rows and their effective positive/negative states stop changing
        # much, additional outer iterations are usually wasted.
        changed_keep_rows = _count_mask_changes(ss.keep_mask, previous_keep_mask)
        changed_label_rows = _count_training_label_state_changes(
            ss.keep_mask,
            ss.positive_mask,
            previous_keep_mask,
            previous_positive_mask
        )
        near_converged =
            (iteration >= near_convergence_check_start_iteration) &&
            (changed_keep_rows <= row_change_threshold) &&
            (changed_label_rows <= row_change_threshold)
        near_convergence_streak = near_converged ? (near_convergence_streak + 1) : 0
        last_plot_iteration = iteration
        last_plot_state = (
            positive_mask = ss.positive_mask,
            confident_positive_mask = ss.confident_positive_mask,
            mined_negative_mask = ss.mined_negative_mask
        )

        if near_convergence_streak >= near_convergence_patience
            if !isnothing(iteration_debug_callback) && last_plot_iteration > 1
                iteration_debug_callback(last_plot_iteration, last_plot_state)
            end
            return β_current
        end

        # Refit on the updated semi-supervised labels and continue.
        X_iteration = X[ss.keep_mask, :]
        β_current = fit_probit_model(X_iteration, ss.labels)
        if !isnothing(feature_importance_logger)
            feature_importance_logger(β_current, X_iteration, iteration)
        end
        previous_keep_mask = copy(ss.keep_mask)
        previous_positive_mask = copy(ss.positive_mask)
    end

    if !isnothing(iteration_debug_callback) && last_plot_iteration > 1
        iteration_debug_callback(last_plot_iteration, last_plot_state)
    end

    return β_current
end

const RUN_LEVEL_PROTEIN_LGBM_HP = (
    num_iterations = 100,
    learning_rate = 0.05,
    max_depth = 4,
    num_leaves = 15,
    min_data_in_leaf = 10,
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

const PROTEIN_MONOTONE_INCREASING_FEATURES = (
    :pg_score,
    :ambiguous_pg_score,
    :shared_peptide_coverage_logit,
    :shared_coverage_log_ratio,
    :peptide_coverage_logit,
    :any_common_peps,
    :coverage_log_ratio,
    :precursor_consensus_prefix_shape,
    :shared_precursor_consensus_prefix_shape,
    :mbr_recovered_peptides,
)

const PROTEIN_MONOTONE_DECREASING_FEATURES = (
    :mbr_only_protein,
)

"""
    _protein_lightgbm_monotone_constraints(feature_names)

Build the LightGBM monotonic-constraint vector in the exact order of the
features supplied to the fitted model. A value of `1` requires predictions to
be nondecreasing as that feature increases, `-1` requires predictions to be
nonincreasing, and `0` leaves the feature unconstrained.
"""
function _protein_lightgbm_monotone_constraints(
    feature_names::AbstractVector{Symbol},
)
    return Int[
        if feature in PROTEIN_MONOTONE_INCREASING_FEATURES
            1
        elseif feature in PROTEIN_MONOTONE_DECREASING_FEATURES
            -1
        else
            0
        end
        for feature in feature_names
    ]
end

function _log_protein_lightgbm_feature_importance(
    model::LightGBMModel,
    iteration::Int;
    context::AbstractString,
)
    DEBUG_CONSOLE_LEVEL[] >= 1 || return nothing
    model_importance = importance(model)
    model_importance === nothing && return nothing

    sorted_gains = sort!(
        collect(model_importance);
        by = pair -> -last(pair),
    )
    lines = [
        "Protein LightGBM feature gains context=$(context)_iter_$(iteration):",
    ]
    for (feature, gain) in sorted_gains
        push!(
            lines,
            "    $(rpad(string(feature), 44)) $(round(Int, gain))",
        )
    end
    @debug_l1 join(lines, "\n")
    return nothing
end

function _protein_training_composition(
    keep_mask::AbstractVector{Bool},
    positive_mask::AbstractVector{Bool},
    mined_negative_mask::AbstractVector{Bool},
    targets::AbstractVector{Bool},
)
    target_positives = 0
    mined_target_negatives = 0
    decoys = 0
    dropped_targets = 0

    @inbounds for row in eachindex(
        keep_mask,
        positive_mask,
        mined_negative_mask,
        targets,
    )
        if targets[row]
            if positive_mask[row]
                target_positives += 1
            elseif mined_negative_mask[row]
                mined_target_negatives += 1
            elseif !keep_mask[row]
                dropped_targets += 1
            end
        elseif keep_mask[row]
            decoys += 1
        end
    end

    return (
        target_positives = target_positives,
        mined_target_negatives = mined_target_negatives,
        decoys = decoys,
        dropped_targets = dropped_targets,
        negative_labels = mined_target_negatives + decoys,
        training_rows =
            target_positives + mined_target_negatives + decoys,
    )
end

function _protein_qvalue_counts(
    qvals::AbstractVector{<:Real},
    targets::AbstractVector{Bool},
    q_value_threshold::Float32,
)
    target_ids = 0
    decoys = 0
    @inbounds for row in eachindex(qvals, targets)
        qvals[row] <= q_value_threshold || continue
        if targets[row]
            target_ids += 1
        else
            decoys += 1
        end
    end
    return (target_ids = target_ids, decoys = decoys)
end

function _log_protein_lightgbm_iteration(
    iteration::Int,
    context::AbstractString,
    current_training,
    qvalue_counts,
    next_training;
    q_value_threshold::Float32,
)
    q_threshold_label = round(Float64(q_value_threshold), digits = 4)
    @debug_l1 "Protein LightGBM semi-supervised iter $iteration context=$context:\n" *
              "    IDs at q≤$q_threshold_label: targets=$(qvalue_counts.target_ids), " *
              "decoys=$(qvalue_counts.decoys)\n" *
              "    current training: target positives=$(current_training.target_positives), " *
              "mined target negatives=$(current_training.mined_target_negatives), " *
              "decoys=$(current_training.decoys), rows=$(current_training.training_rows)\n" *
              "    next training: target positives=$(next_training.target_positives), " *
              "mined target negatives=$(next_training.mined_target_negatives), " *
              "decoys=$(next_training.decoys), dropped targets=$(next_training.dropped_targets)"
    return nothing
end

"""
    fit_protein_lightgbm_semisupervised(feature_data, y, initial_scores, prefix_shape, n_non_mbr_peptides; kwargs...)

Fit a run-level protein LightGBM model using the same protein-specific
semi-supervised labels and negative-mining rules as the probit model. As in
precursor scoring, stop when a new iteration fails to improve target IDs at
the configured q-value threshold by at least 1%, retain the iteration with the
most passing targets, and cap training at eight iterations by default.
"""
function fit_protein_lightgbm_semisupervised(
    feature_data::AbstractDataFrame,
    y::AbstractVector{Bool},
    initial_scores::AbstractVector{<:Real},
    prefix_shape::AbstractVector{<:Real},
    n_non_mbr_peptides::AbstractVector{<:Integer};
    q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold::Float32 = -0.20f0,
    min_pep_neg_threshold::Float32 = 0.90f0,
    max_positive_pep_threshold::Float32 = 1.0f0,
    min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    n_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
    context::AbstractString = "protein_lightgbm",
    iteration_debug_callback = nothing,
    lgbm_hp::NamedTuple = RUN_LEVEL_PROTEIN_LGBM_HP,
)
    nrow(feature_data) == length(y) ||
        throw(ArgumentError("feature_data must have the same number of rows as y"))
    length(initial_scores) == length(y) ||
        throw(ArgumentError("initial_scores must have the same length as y"))
    length(prefix_shape) == length(y) ||
        throw(ArgumentError("prefix_shape must have the same length as y"))
    length(n_non_mbr_peptides) == length(y) ||
        throw(ArgumentError("n_non_mbr_peptides must have the same length as y"))

    feature_names = Symbol.(names(feature_data))
    monotone_constraints =
        _protein_lightgbm_monotone_constraints(feature_names)
    increasing_features = feature_names[monotone_constraints .== 1]
    decreasing_features = feature_names[monotone_constraints .== -1]
    monotone_constraints_method =
        get(lgbm_hp, :monotone_constraints_method, "basic")
    @debug_l1 "Protein LightGBM monotone constraints context=$context: " *
              "method=$monotone_constraints_method; " *
              "nondecreasing features=$(join(increasing_features, ", ")); " *
              "nonincreasing features=$(join(decreasing_features, ", "))"

    last_plot_iteration = 1

    ss_initial = build_protein_semisupervised_training_set(
        initial_scores,
        y,
        prefix_shape,
        n_non_mbr_peptides;
        q_value_threshold = q_value_threshold,
        max_positive_pep_threshold = max_positive_pep_threshold,
        mined_negative_prefix_shape_threshold =
            min_prefix_shape_neg_threshold,
        mined_negative_pep_threshold = min_pep_neg_threshold,
        keep_non_mined_targets_as_positive = true,
    )
    initial_keep_mask = ss_initial.keep_mask
    initial_labels = ss_initial.labels
    last_plot_state = (
        positive_mask = ss_initial.positive_mask,
        confident_positive_mask = ss_initial.confident_positive_mask,
        mined_negative_mask = ss_initial.mined_negative_mask,
    )

    initial_target_count = sum(initial_labels)
    initial_mined_target_negative_count =
        sum(ss_initial.mined_negative_mask)
    initial_decoy_count = count(.!y .& initial_keep_mask)
    initial_negative_label_count =
        initial_mined_target_negative_count + initial_decoy_count
    if initial_target_count < 10 || initial_negative_label_count < 10
        @user_warn "Protein LightGBM semi-supervised iteration has insufficient data; falling back to all rows for initial pg_score fit" context = context iteration = 1 selected_target_positives = initial_target_count selected_mined_target_negatives = initial_mined_target_negative_count selected_decoys = initial_decoy_count selected_negative_labels = initial_negative_label_count
        initial_keep_mask = trues(length(y))
        initial_labels = y
        last_plot_state = (
            positive_mask = copy(y),
            confident_positive_mask = ss_initial.confident_positive_mask,
            mined_negative_mask = falses(length(y)),
        )
    end

    model_current = fit_lightgbm_model(
        build_lightgbm_classifier(
            ;
            lgbm_hp...,
            monotone_constraints = monotone_constraints,
        ),
        view(feature_data, initial_keep_mask, :),
        initial_labels;
        positive_label = true,
    )
    _log_protein_lightgbm_feature_importance(
        model_current,
        1;
        context = context,
    )
    if !isnothing(iteration_debug_callback)
        iteration_debug_callback(1, last_plot_state)
    end

    current_training = _protein_training_composition(
        initial_keep_mask,
        last_plot_state.positive_mask,
        last_plot_state.mined_negative_mask,
        y,
    )

    # Evaluate each fitted model before deciding whether to train the next one.
    # As in precursor scoring, continuation depends only on target IDs passing
    # the q-value threshold, not on training-label or row churn.
    previous_target_ids = -1
    best_state = nothing
    for model_iteration in 1:n_iterations
        iteration_scores = lightgbm_predict(
            model_current,
            feature_data;
            output_type = Float32,
        )
        ss = build_protein_semisupervised_training_set(
            iteration_scores,
            y,
            prefix_shape,
            n_non_mbr_peptides;
            q_value_threshold = q_value_threshold,
            max_positive_pep_threshold = max_positive_pep_threshold,
            mined_negative_prefix_shape_threshold =
                min_prefix_shape_neg_threshold,
            mined_negative_pep_threshold = min_pep_neg_threshold,
            keep_non_mined_targets_as_positive = false,
        )

        next_training = _protein_training_composition(
            ss.keep_mask,
            ss.positive_mask,
            ss.mined_negative_mask,
            y,
        )
        qvalue_counts =
            _protein_qvalue_counts(ss.qvals, y, q_value_threshold)
        next_iteration = model_iteration + 1

        _log_protein_lightgbm_iteration(
            model_iteration,
            context,
            current_training,
            qvalue_counts,
            next_training;
            q_value_threshold = q_value_threshold,
        )

        current_state = (
            model = model_current,
            target_q01 = qvalue_counts.target_ids,
            decoy_q01 = qvalue_counts.decoys,
            iter = model_iteration,
            plot_iteration = last_plot_iteration,
            plot_state = last_plot_state,
        )
        if model_iteration > 1 && !_scoring_target_gain_sufficient(
            previous_target_ids,
            qvalue_counts.target_ids;
            min_fraction = min_gain,
        )
            best_state =
                _scoring_better_iteration_state(best_state, current_state)
            @debug_l1 "Protein LightGBM semi-supervised stopping context=$context after iter $model_iteration: " *
                      "q≤$(round(Float64(q_value_threshold), digits = 4)) target IDs=$(qvalue_counts.target_ids) " *
                      "did not improve by $(round(100 * min_gain, digits = 2))% over $previous_target_ids; " *
                      "using iter $(best_state.iter) with target IDs=$(best_state.target_q01)"
            if !isnothing(iteration_debug_callback) &&
               best_state.plot_iteration > 1
                iteration_debug_callback(
                    best_state.plot_iteration,
                    best_state.plot_state,
                )
            end
            return best_state.model
        end

        best_state = _scoring_better_iteration_state(
            best_state,
            current_state,
        )
        if model_iteration == n_iterations
            @debug_l1 "Protein LightGBM semi-supervised stopping context=$context: " *
                      "hit max iterations $n_iterations; using iter $(best_state.iter) " *
                      "with q≤$(round(Float64(q_value_threshold), digits = 4)) " *
                      "target IDs=$(best_state.target_q01)"
            if !isnothing(iteration_debug_callback) &&
               best_state.plot_iteration > 1
                iteration_debug_callback(
                    best_state.plot_iteration,
                    best_state.plot_state,
                )
            end
            return best_state.model
        end

        if next_training.target_positives < 10 ||
           next_training.negative_labels < 10
            @debug_l1 "Protein LightGBM semi-supervised stopping context=$context after iter $model_iteration: " *
                      "next iteration has insufficient class support " *
                      "(target positives=$(next_training.target_positives), " *
                      "negative labels=$(next_training.negative_labels)); " *
                      "using iter $(best_state.iter) with q≤$(round(Float64(q_value_threshold), digits = 4)) " *
                      "target IDs=$(best_state.target_q01)"
            @user_warn "Protein LightGBM semi-supervised iteration has insufficient data; using best iteration model" context = context iteration = next_iteration selected_target_positives = next_training.target_positives selected_mined_target_negatives = next_training.mined_target_negatives selected_decoys = next_training.decoys selected_negative_labels = next_training.negative_labels best_iteration = best_state.iter best_target_ids = best_state.target_q01
            if !isnothing(iteration_debug_callback) &&
               best_state.plot_iteration > 1
                iteration_debug_callback(
                    best_state.plot_iteration,
                    best_state.plot_state,
                )
            end
            return best_state.model
        end

        continue_reason = if model_iteration == 1
            "established target-ID baseline=$(qvalue_counts.target_ids)"
        else
            "target IDs improved from $previous_target_ids to $(qvalue_counts.target_ids), meeting the $(round(100 * min_gain, digits = 2))% minimum gain"
        end
        @debug_l1 "Protein LightGBM semi-supervised continuing context=$context after iter $model_iteration: " *
                  "$continue_reason; fitting iter $next_iteration"

        previous_target_ids = qvalue_counts.target_ids
        last_plot_iteration = next_iteration
        last_plot_state = (
            positive_mask = ss.positive_mask,
            confident_positive_mask = ss.confident_positive_mask,
            mined_negative_mask = ss.mined_negative_mask,
        )
        model_current = fit_lightgbm_model(
            build_lightgbm_classifier(
                ;
                lgbm_hp...,
                monotone_constraints = monotone_constraints,
            ),
            view(feature_data, ss.keep_mask, :),
            ss.labels;
            positive_label = true,
        )
        _log_protein_lightgbm_feature_importance(
            model_current,
            next_iteration;
            context = context,
        )
        current_training = next_training
    end

    return model_current
end

"""
    calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}, X_mean::Vector{Float64}, X_std::Vector{Float64})

Calculate probit probability scores for new data.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `β::Vector{Float64}`: Fitted coefficients
- `X_mean::Vector{Float64}`: Feature means from training
- `X_std::Vector{Float64}`: Feature standard deviations from training

# Returns
- `Vector{Float64}`: Probability scores
"""
function calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}
    #, X_mean::Vector{Float64}, X_std::Vector{Float64}
    )
    # Standardize using training statistics
    X_standardized = X#(X .- X_mean') ./ X_std'
    
    # Add intercept
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])
    
    # Create data chunks
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, size(X, 1) / n_chunks))
    data_chunks = Iterators.partition(1:size(X, 1), chunk_size)
    
    # Calculate probabilities
    prob_scores = zeros(Float64, size(X, 1))
    Pioneer.ModelPredictProbs!(prob_scores, X_df, β, data_chunks)
    
    return prob_scores
end

# ==========================================================
# Multi-fold Cross-Validation Probit Analysis
# ==========================================================

"""
    detect_unique_cv_folds(precursors::LibraryPrecursors) -> Vector{UInt8}

Detect unique CV fold values from the library precursors.

# Returns
- Sorted vector of unique CV fold values
"""
function detect_unique_cv_folds(precursors::LibraryPrecursors)
    return sort(unique(precursors.pid_to_cv_fold))
end

"""
    assign_protein_group_cv_folds!(all_protein_groups::DataFrame, 
                                  protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}})

Assign CV fold to each protein group based on a pre-built mapping.

# Arguments
- `all_protein_groups`: DataFrame of protein groups to update (modified in-place)
- `protein_to_cv_fold`: Pre-built dictionary mapping protein names to cv_fold assignments

# Process
Adds cv_fold column to protein groups DataFrame based on the mapping
"""
function assign_protein_group_cv_folds!(
    all_protein_groups::DataFrame,
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}
)
    # Assign cv_fold to protein groups using the pre-built mapping
    cv_folds = Vector{UInt8}(undef, nrow(all_protein_groups))
    missing_count = 0
    
    for (i, protein_name) in enumerate(all_protein_groups.protein_name)
        if haskey(protein_to_cv_fold, protein_name)
            cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
        else
            missing_count += 1
            # Default to fold 0 if no matching peptide found
            cv_folds[i] = UInt8(0)
        end
    end
    
    if missing_count > 0
        @user_warn "There were $missing_count protein groups without matching peptides, assigned to fold 0"
    end

    all_protein_groups[!, :cv_fold] = cv_folds
    return protein_to_cv_fold
end

"""
    apply_probit_scores_multifold!(pg_refs::Vector{ProteinGroupFileReference},
                                  protein_to_cv_fold::Dict{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
                                  models::Dict{UInt8, Vector{Float64}},
                                  feature_names::Vector{Symbol})

Apply probit models to protein group files based on their CV fold. When probit
scoring is unavailable, convert every initial log-sum score to a probability.

# Arguments
- `pg_refs`: Vector of protein group file references to update
- `protein_to_cv_fold`: Pre-built mapping of protein names to cv_fold
- `models`: Dictionary mapping CV fold to fitted model coefficients
- `feature_names`: Feature names used in the model
- `use_probit_scores`: Whether every fold has a trained probit model
"""
function apply_probit_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    models::Dict{UInt8, LightGBMModel},
    feature_names::Vector{Symbol};
    use_probit_scores::Bool,
)
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Assign CV folds using pre-built mapping
            cv_folds = Vector{UInt8}(undef, nrow(df))
            for (i, protein_name) in enumerate(df.protein_name)
                if haskey(protein_to_cv_fold, protein_name)
                    cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
                else
                    # Default to fold 0 if not found
                    cv_folds[i] = UInt8(0)
                end
            end
            df[!, :cv_fold] = cv_folds
            
            # Save original scores for comparison
            df[!, :old_pg_score] = copy(df.pg_score)
            
            if use_probit_scores
                for (fold, model) in models
                    mask = df.cv_fold .== fold
                    if sum(mask) > 0
                        scoring_df = df[mask, :]
                        df[mask, :pg_score] = lightgbm_predict(
                            model,
                            view(scoring_df, :, feature_names);
                            output_type = Float32,
                        )
                    end
                end
            else
                df[!, :pg_score] = _initial_protein_probabilities(df.old_pg_score)
            end
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Remove temporary training-only columns before returning.
            columns_to_remove = Symbol[:cv_fold]
            hasproperty(df, :_ambiguous_peptide_count) &&
                push!(columns_to_remove, :_ambiguous_peptide_count)
            select!(df, Not(columns_to_remove))
            
            return df
        end
    end
end

function _initial_protein_probabilities(scores::AbstractVector{<:Real})
    return Float32.(clamp.(
        1.0f0 .- exp.(-scores),
        1.0f-6,
        1.0f0 - 1.0f-6,
    ))
end

function _protein_probit_training_data_sufficient(
    all_protein_groups::DataFrame,
    cv_folds::Vector{UInt8},
)
    targets = all_protein_groups.target
    folds = all_protein_groups.cv_fold
    for test_fold in cv_folds
        n_train_targets = 0
        n_train_decoys = 0
        @inbounds for row in eachindex(targets, folds)
            folds[row] == test_fold && continue
            if targets[row]
                n_train_targets += 1
            else
                n_train_decoys += 1
            end
        end
        min(n_train_targets, n_train_decoys) >= 10 || return false
    end
    return true
end

"""
    perform_probit_analysis_multifold(all_protein_groups::DataFrame,
                                     qc_folder::String,
                                     pg_refs::Vector{ProteinGroupFileReference},
                                     precursors::LibraryPrecursors;
                                     protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
                                     show_improvement = true,
                                     skip_scoring = false)

Perform probit regression analysis with automatic CV fold detection from library.

# Arguments
- `all_protein_groups`: DataFrame with protein group data
- `qc_folder`: Folder for QC plots
- `pg_refs`: Protein group file references
- `precursors`: Library precursors containing CV fold information
- `protein_to_cv_fold`: Pre-built mapping of proteins to CV folds
- `show_improvement`: Whether to report improvement metrics

# Process
1. Detects unique CV folds from library
2. Assigns CV folds to protein groups based on the pre-built mapping
3. Trains separate probit models only when every fold has sufficient data
4. Applies every model, or converts every initial score to a probability
5. Updates protein group files if provided

Returns `true` only when probit models were trained and applied to every fold.
"""
function perform_probit_analysis_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    show_improvement = true,
    skip_scoring = false,
    write_qc_plots::Bool = true,
    train_q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold_itr::Float32 = -0.20f0,
    min_pep_neg_threshold_itr::Float32 = 0.90f0
)
    #skip_scoring = true 
    #@user_info "Skipped scoring!!!"
    # 1. Detect unique CV folds from library
    unique_cv_folds = detect_unique_cv_folds(precursors)

    # 2. Assign CV folds to protein groups based on the mapping
    assign_protein_group_cv_folds!(all_protein_groups, protein_to_cv_fold)

    feature_names = protein_probit_feature_names()
    feature_filter_df = all_protein_groups
    remove_zero_variance_columns!(feature_names, feature_filter_df)

    features_available = !isempty(feature_names)
    if !features_available
        @user_warn "No valid features remaining for multifold probit regression after filtering"
    end

    rows_available = nrow(feature_filter_df) > 0
    if !rows_available
        @user_warn "No protein groups are available for protein probit regression"
    end

    models = Dict{UInt8, LightGBMModel}()
    use_probit_scores = !skip_scoring &&
                        features_available &&
                        rows_available &&
                        _protein_probit_training_data_sufficient(
                            all_protein_groups,
                            unique_cv_folds,
                        )

    if use_probit_scores
        for test_fold in unique_cv_folds
            train_mask = all_protein_groups.cv_fold .!= test_fold
            train_df = all_protein_groups[train_mask, :]

            y_train = train_df.target
            initial_scores_train = train_df.pg_score
            context = "protein_lightgbm_multifold_fold_$(test_fold)"
            iteration_debug_callback = write_qc_plots ?
                make_protein_iteration_plot_callback(
                    train_df,
                    test_fold,
                    qc_folder;
                    context = context,
                    file_idx_to_name = file_idx_to_name
                ) :
                nothing
            # Fit model
            fitted_model = fit_protein_lightgbm_semisupervised(
                view(train_df, :, feature_names),
                y_train,
                initial_scores_train,
                train_df.precursor_consensus_prefix_shape,
                train_df.n_non_mbr_peptides;
                q_value_threshold = train_q_value_threshold,
                min_prefix_shape_neg_threshold = min_prefix_shape_neg_threshold_itr,
                min_pep_neg_threshold = min_pep_neg_threshold_itr,
                context = context,
                iteration_debug_callback = iteration_debug_callback
            )
            models[test_fold] = fitted_model
        end
    end
    
    all_protein_groups[!, :old_pg_score] = copy(all_protein_groups[!, :pg_score])  # Save for comparison
    
    if use_probit_scores
        for test_fold in unique_cv_folds
            test_mask = all_protein_groups.cv_fold .== test_fold
            n_test = sum(test_mask)
            
            if n_test > 0
                scoring_df =
                    view(all_protein_groups, test_mask, feature_names)
                all_protein_groups[test_mask, :pg_score] =
                    lightgbm_predict(
                        models[test_fold],
                        scoring_df;
                        output_type = Float32,
                    )
            end
        end
    else
        all_protein_groups[!, :pg_score] =
            _initial_protein_probabilities(all_protein_groups.old_pg_score)
    end
    
    if show_improvement && use_probit_scores
        # Calculate improvement at 1% FDR
        old_qvalues = zeros(Float32, nrow(all_protein_groups))
        new_qvalues = zeros(Float32, nrow(all_protein_groups))
        
        get_qvalues!(all_protein_groups[!, :old_pg_score], all_protein_groups[!, :target], old_qvalues)
        get_qvalues!(all_protein_groups[!, :pg_score], all_protein_groups[!, :target], new_qvalues)
        
        old_passing = sum((old_qvalues .<= 0.01f0) .& all_protein_groups.target)
        new_passing = sum((new_qvalues .<= 0.01f0) .& all_protein_groups.target)
        
        percent_improvement = 0.0
        if old_passing > 0
            percent_improvement = round(100.0 * (new_passing - old_passing) / old_passing, digits=1)
        end
    end

    # 8. Update protein group files if provided
    if !isempty(pg_refs)
        apply_probit_scores_multifold!(
            pg_refs,
            protein_to_cv_fold,
            models,
            feature_names;
            use_probit_scores = use_probit_scores,
        )
    end
    
    # Clean up temporary column
    select!(all_protein_groups, Not(:cv_fold))
    return use_probit_scores
end
