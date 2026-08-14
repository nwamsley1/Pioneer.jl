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
    build_protein_semisupervised_training_set(scores, targets, prefix_shape, n_non_mbr_peptides; q_value_threshold = 0.01f0, max_positive_pep_threshold = 1.0f0, mined_negative_prefix_shape_threshold = -0.20f0, mined_negative_pep_threshold = 0.90f0, keep_non_mined_targets_as_positive = true)

Build labels for semi-supervised run-level protein training from a score vector.
Targets are mined as negatives when `PEP >= mined_negative_pep_threshold` or
when `prefix_shape <= mined_negative_prefix_shape_threshold` and the target has
at most one non-MBR unique peptide.
Shared-peptide scores, counts, and prefix shapes do not affect negative mining.
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
        low_shape_limited_support_target =
            targets[i] &&
            (n_non_mbr_peptides[i] <= 1) &&
            (Float32(prefix_shape[i]) <= mined_negative_prefix_shape_threshold)
        candidate_confident_positive_mask[i] = targets[i] &&
                                               (qvals[i] <= q_value_threshold) &&
                                               (peps[i] <= max_positive_pep_threshold)
        mined_negative_mask[i] = targets[i] &&
                                 ((peps[i] >= mined_negative_pep_threshold) ||
                                  low_shape_limited_support_target)
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
    :ambiguous_pg_score,
    :shared_peptide_coverage_logit,
    :shared_coverage_log_ratio,
    :peptide_coverage_logit,
    :all_peptide_coverage_logit,
    :any_common_peps,
    :coverage_log_ratio,
    :precursor_consensus_prefix_shape,
    :shared_precursor_consensus_prefix_shape,
    :single_non_mbr_prefix_shape,
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

Fit a run-level protein LightGBM model using protein-specific semi-supervised
labels and negative-mining rules. As in precursor scoring, stop when a new
iteration fails to improve target IDs at the configured q-value threshold by
at least 1%, retain the iteration with the most passing targets, and cap
training at eight iterations by default.
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

# ==========================================================
# Multi-fold run-level protein scoring
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
    apply_protein_scores_multifold!(pg_refs::Vector{ProteinGroupFileReference},
                                    protein_to_cv_fold::Dict{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
                                    models::Dict{UInt8, LightGBMModel},
                                    feature_names::Vector{Symbol})

Apply run-level protein models to protein-group files based on their CV fold.
When model scoring is unavailable, convert every initial log-sum score to a
probability.

# Arguments
- `pg_refs`: Vector of protein group file references to update
- `protein_to_cv_fold`: Pre-built mapping of protein names to cv_fold
- `models`: Dictionary mapping CV fold to fitted LightGBM models
- `feature_names`: Feature names used in the model
- `use_model_scores`: Whether every fold has a trained model
"""
function apply_protein_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    models::Dict{UInt8, LightGBMModel},
    feature_names::Vector{Symbol};
    use_model_scores::Bool,
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
            
            if use_model_scores
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

function _protein_model_training_data_sufficient(
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
    perform_protein_scoring_multifold(all_protein_groups::DataFrame,
                                     qc_folder::String,
                                     pg_refs::Vector{ProteinGroupFileReference},
                                     precursors::LibraryPrecursors;
                                     protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
                                     skip_scoring = false)

Perform run-level protein LightGBM scoring with library-derived CV folds.

# Arguments
- `all_protein_groups`: DataFrame with protein group data
- `qc_folder`: Folder for QC plots
- `pg_refs`: Protein group file references
- `precursors`: Library precursors containing CV fold information
- `protein_to_cv_fold`: Pre-built mapping of proteins to CV folds

# Process
1. Detects unique CV folds from library
2. Assigns CV folds to protein groups based on the pre-built mapping
3. Trains separate LightGBM models only when every fold has sufficient data
4. Applies every model, or converts every initial score to a probability
5. Updates protein group files if provided

Returns `true` only when models were trained and applied to every fold.
"""
function perform_protein_scoring_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    skip_scoring = false,
    write_qc_plots::Bool = true,
    train_q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold_itr::Float32 = -0.20f0,
    min_pep_neg_threshold_itr::Float32 = 0.90f0
)
    unique_cv_folds = detect_unique_cv_folds(precursors)

    assign_protein_group_cv_folds!(all_protein_groups, protein_to_cv_fold)

    feature_names = run_level_protein_feature_names()
    feature_filter_df = all_protein_groups
    remove_zero_variance_columns!(feature_names, feature_filter_df)

    features_available = !isempty(feature_names)
    if !features_available
        @user_warn "No valid features remaining for run-level protein scoring after filtering"
    end

    rows_available = nrow(feature_filter_df) > 0
    if !rows_available
        @user_warn "No protein groups are available for run-level protein scoring"
    end

    models = Dict{UInt8, LightGBMModel}()
    use_model_scores = !skip_scoring &&
                       features_available &&
                       rows_available &&
                       _protein_model_training_data_sufficient(
                           all_protein_groups,
                           unique_cv_folds,
                       )

    if use_model_scores
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
    
    initial_pg_scores = copy(all_protein_groups.pg_score)
    
    if use_model_scores
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
            _initial_protein_probabilities(initial_pg_scores)
    end
    
    if !isempty(pg_refs)
        apply_protein_scores_multifold!(
            pg_refs,
            protein_to_cv_fold,
            models,
            feature_names;
            use_model_scores = use_model_scores,
        )
    end
    
    # Clean up temporary column
    select!(all_protein_groups, Not(:cv_fold))
    return use_model_scores
end
