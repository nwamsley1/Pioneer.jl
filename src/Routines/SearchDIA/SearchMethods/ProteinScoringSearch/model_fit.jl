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
    @user_info "Protein probit model coefficients context=$(context) intercept=$(round(intercept, digits=6)) n_features=$(n_features) importance_metric=abs(coefficient*feature_std)"

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
        @user_info "Protein probit feature importance context=$(context) rank=$(rank) feature=$(summary.feature) coefficient=$(round(summary.coefficient, digits=6)) feature_std=$(round(summary.feature_std, digits=6)) one_sigma_effect=$(round(summary.one_sigma_effect, digits=6)) abs_one_sigma_effect=$(round(summary.abs_one_sigma_effect, digits=6))"
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
    build_protein_semisupervised_training_set(scores, targets, prefix_shape, n_peptides; q_value_threshold = 0.01f0, max_positive_pep_threshold = 1.0f0, mined_negative_prefix_shape_threshold = -0.20f0, mined_negative_pep_threshold = 0.90f0, keep_non_mined_targets_as_positive = true)

Build labels for semi-supervised protein probit training from a score vector.
Targets with `PEP >= mined_negative_pep_threshold` are mined as negatives.
Singleton targets with `prefix_shape <= mined_negative_prefix_shape_threshold`
are mined as negatives regardless of q-value. If
`keep_non_mined_targets_as_positive=true`, remaining targets stay positive for
training; otherwise only q-value-passing targets stay positive and the rest are
dropped.
"""
function build_protein_semisupervised_training_set(
    scores::AbstractVector{<:Real},
    targets::AbstractVector{Bool},
    prefix_shape::AbstractVector{<:Real},
    n_peptides::AbstractVector{<:Integer};
    q_value_threshold::Float32 = 0.01f0,
    max_positive_pep_threshold::Float32 = 1.0f0,
    mined_negative_prefix_shape_threshold::Float32 = -0.20f0,
    mined_negative_pep_threshold::Float32 = 0.90f0,
    keep_non_mined_targets_as_positive::Bool = true
)
    n = length(scores)
    length(targets) == n || throw(ArgumentError("targets must have the same length as scores"))
    length(prefix_shape) == n || throw(ArgumentError("prefix_shape must have the same length as scores"))
    length(n_peptides) == n || throw(ArgumentError("n_peptides must have the same length as scores"))
    qvals = Vector{Float32}(undef, n)
    peps = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals)
    get_PEP!(scores, targets, peps)
    
    confident_positive_mask = BitVector(undef, n)
    candidate_confident_positive_mask = BitVector(undef, n)
    mined_negative_mask = BitVector(undef, n)
    positive_mask = BitVector(undef, n)
    keep_mask = BitVector(undef, n)

    @inbounds for i in eachindex(scores, targets, prefix_shape, n_peptides, qvals, peps)
        low_shape_singleton_target = targets[i] &&
                                     (n_peptides[i] == 1) &&
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
    fit_probit_model_semisupervised(X, y, initial_scores, prefix_shape, n_peptides; q_value_threshold = 0.01f0, min_prefix_shape_neg_threshold = -0.20f0, min_pep_neg_threshold = 0.90f0, max_positive_pep_threshold = 1.0f0, n_iterations = 10, context = "protein_probit", iteration_debug_callback = nothing)

Fit the protein probit model by seeding iteration 1 labels from raw initial scores,
then fitting and refining with the full feature set. Later outer iterations stop
early when label churn stays small for two consecutive rounds.
"""
function fit_probit_model_semisupervised(
    X::Matrix{Float64},
    y::AbstractVector{Bool},
    initial_scores::AbstractVector{<:Real},
    prefix_shape::AbstractVector{<:Real},
    n_peptides::AbstractVector{<:Integer};
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
    length(n_peptides) == length(y) || throw(ArgumentError("n_peptides must have the same length as y"))

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
        n_peptides;
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
            n_peptides;
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

Apply probit models to protein group files based on their CV fold.

# Arguments
- `pg_refs`: Vector of protein group file references to update
- `protein_to_cv_fold`: Pre-built mapping of protein names to cv_fold
- `models`: Dictionary mapping CV fold to fitted model coefficients
- `feature_names`: Feature names used in the model
"""
function apply_probit_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    models::Dict{UInt8, Vector{Float64}},
    feature_names::Vector{Symbol};
    skip_scoring = false
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
            
            # Apply appropriate model to each fold (skip if skip_scoring = true)
            if !skip_scoring
                for (fold, model) in models
                    mask = df.cv_fold .== fold
                    if sum(mask) > 0
                        scoring_df = df[mask, :]
                        X = Matrix{Float64}(scoring_df[:, feature_names])
                        df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
                    end
                end
            else
                # When skip_scoring is true, pg_score contains log-sum scores: -sum(log1p(-p))
                # Convert to probabilities: p = 1 - exp(-pg_score)
                df[!, :pg_score] = 1.0f0 .- exp.(-df.pg_score)
                # Clamp to avoid numerical issues with logodds function
                df[!, :pg_score] = clamp.(df.pg_score, 1f-6, 1f0 - 1f-6)
            end
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Remove temporary cv_fold column before returning
            select!(df, Not(:cv_fold))
            
            return df
        end
    end
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
3. Trains separate probit model for each fold
4. Applies models to held-out data
5. Updates protein group files if provided
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
    log_feature_importance::Bool = true,
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

    # 3. Check distribution
    fold_counts = Dict{UInt8, Int}()
    for fold in all_protein_groups.cv_fold
        fold_counts[fold] = get(fold_counts, fold, 0) + 1
    end

    feature_names = protein_probit_feature_names()
    feature_filter_df = all_protein_groups
    remove_zero_variance_columns!(feature_names, feature_filter_df)

    if isempty(feature_names)
        @user_warn "No valid features remaining for multifold probit regression after filtering"
        return
    end

    if nrow(feature_filter_df) == 0
        @user_warn "No protein groups are available for protein probit regression"
        return
    end

    # 5. Train probit model for each fold (skip if skip_scoring = true)
    models = Dict{UInt8, Vector{Float64}}()

    if !skip_scoring
        for test_fold in unique_cv_folds
            # Get training data (all folds except test_fold)
            train_mask = all_protein_groups.cv_fold .!= test_fold
            
            # Check if we have sufficient data
            n_train_targets = sum(all_protein_groups[train_mask, :target])
            n_train_decoys = sum(.!all_protein_groups[train_mask, :target])
            
            if n_train_targets < 10 || n_train_decoys < 10
                @user_warn "Insufficient training data for fold $test_fold" targets=n_train_targets decoys=n_train_decoys
                continue
            end

            train_df = all_protein_groups[train_mask, :]

            X_train = Matrix{Float64}(train_df[:, feature_names])
            y_train = train_df.target
            initial_scores_train = train_df.pg_score
            context = "protein_probit_multifold_fold_$(test_fold)"
            iteration_debug_callback = write_qc_plots ?
                make_protein_iteration_plot_callback(
                    train_df,
                    test_fold,
                    qc_folder;
                    context = context,
                    file_idx_to_name = file_idx_to_name
                ) :
                nothing
            feature_importance_logger = log_feature_importance ?
                make_protein_feature_importance_logger(feature_names; context = context) :
                nothing
            
            # Fit model
            β_fitted = fit_probit_model_semisupervised(
                X_train,
                y_train,
                initial_scores_train,
                train_df.precursor_consensus_prefix_shape,
                train_df.n_peptides;
                q_value_threshold = train_q_value_threshold,
                min_prefix_shape_neg_threshold = min_prefix_shape_neg_threshold_itr,
                min_pep_neg_threshold = min_pep_neg_threshold_itr,
                context = context,
                iteration_debug_callback = iteration_debug_callback,
                feature_importance_logger = feature_importance_logger
            )
            models[test_fold] = β_fitted
        end
    end
    
    # 6. Apply models to their respective test folds (skip if skip_scoring = true)
    all_protein_groups[!, :old_pg_score] = copy(all_protein_groups[!, :pg_score])  # Save for comparison
    
    if !skip_scoring
        for test_fold in unique_cv_folds
            if !haskey(models, test_fold)
                @user_warn "No model available for fold $test_fold, keeping original scores"
                continue
            end
            
            test_mask = all_protein_groups.cv_fold .== test_fold
            n_test = sum(test_mask)
            
            if n_test > 0
                X_test = Matrix{Float64}(all_protein_groups[test_mask, feature_names])
                prob_scores = calculate_probit_scores(X_test, models[test_fold])
                all_protein_groups[test_mask, :pg_score] = Float32.(prob_scores)
            end
        end
    end
    
    # 7. Report improvement if requested (skip if skip_scoring = true)
    if show_improvement && !skip_scoring
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
            skip_scoring = skip_scoring
        )
    end
    
    # Clean up temporary column
    select!(all_protein_groups, Not(:cv_fold))
end
