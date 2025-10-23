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

#############################################################################
# Random Precursor Pairing Implementation
#############################################################################
#
# PAIRING STRATEGY:
# Precursors are randomly paired within iRT bins (size 1000) regardless of
# target/decoy status. This allows:
#   - Target-Target pairs
#   - Target-Decoy pairs
#   - Decoy-Decoy pairs
#
# ALGORITHM:
# 1. Within each iRT bin, collect all unique precursor_idx values
# 2. Randomly shuffle using fixed seed (PAIRING_RANDOM_SEED = 1844)
# 3. Pair consecutive elements: (1,2), (3,4), (5,6), etc.
# 4. Both members of a pair share the same pair_id
# 5. If odd number of precursors, last one becomes singleton pair
#
# BENEFITS:
# - Every precursor has a pair (except 0-1 singleton per bin)
# - More diverse training data for ML models
# - Simpler implementation than target-decoy constrained pairing
# - ML models learn from all pairing types (T-T, T-D, D-D)
#
#############################################################################

const PAIRING_RANDOM_SEED = 1844  # Fixed seed for reproducible pairing
const IRT_BIN_SIZE = 1000

#############################################################################
# Running Statistics Helper Functions
#############################################################################

"""
    update_pair_statistics(current_stats, new_prob::Float32)

Updates running statistics for MBR pair probabilities in a memory-efficient manner.
Maintains best 2, worst 2, running mean, and count without storing all values.
"""
function update_pair_statistics(current_stats, new_prob::Float32)
    # Update count and running mean
    new_count = current_stats.count_pairs + 1
    new_mean = (current_stats.mean_prob * current_stats.count_pairs + new_prob) / new_count

    # Update best probabilities (existing logic enhanced)
    new_best_1, new_best_2 = if new_prob > current_stats.best_prob_1
        (new_prob, current_stats.best_prob_1)
    elseif new_prob > current_stats.best_prob_2
        (current_stats.best_prob_1, new_prob)
    else
        (current_stats.best_prob_1, current_stats.best_prob_2)
    end

    # Update worst probabilities (new logic)
    new_worst_1, new_worst_2 = if new_count == 1
        (new_prob, zero(Float32))  # First probability
    elseif new_count == 2
        (min(current_stats.worst_prob_1, new_prob), max(current_stats.worst_prob_1, new_prob))
    elseif new_prob < current_stats.worst_prob_1
        (new_prob, current_stats.worst_prob_1)  # New minimum
    elseif new_prob < current_stats.worst_prob_2
        (current_stats.worst_prob_1, new_prob)  # New second minimum
    else
        (current_stats.worst_prob_1, current_stats.worst_prob_2)  # No change
    end

    return merge(current_stats, (
        best_prob_1 = new_best_1,
        best_prob_2 = new_best_2,
        worst_prob_1 = new_worst_1,
        worst_prob_2 = new_worst_2,
        mean_prob = new_mean,
        count_pairs = new_count
    ))
end

function getIrtBins(irts::AbstractVector{R}) where {R <:Real}
    sort_idx = sortperm(irts)
    bin_idx, bin_count = zero(UInt32), zero(UInt32)
    bin_idxs = similar(irts, UInt32, length(irts))
    for idx in sort_idx
        bin_count += one(UInt32)
        bin_idxs[idx] = bin_idx
        if bin_count >= IRT_BIN_SIZE
            bin_idx += one(UInt32)
            bin_count = zero(UInt32)
        end
    end
    return bin_idxs 
end


function getIrtBins!(psms::AbstractDataFrame)
    psms[!, :irt_bin_idx] = getIrtBins(psms.irt_pred)
    return psms
end


@inline function irt_residual(psms::AbstractDataFrame, idx::Integer)
    return Float32(psms.irt_pred[idx] - psms.irt_obs[idx])
end


function assign_random_target_decoy_pairs!(psms::DataFrame)
    last_pair_id = zero(UInt32)
    psms[!,:pair_id] = zeros(UInt32, nrow(psms))  # Initialize pair_id column
    psms[!,:irt_bin_idx] = getIrtBins(psms.irt_pred)  # Ensure irt_bin_idx column exists

    irt_bin_groups = groupby(psms, :irt_bin_idx)
    for (irt_bin_idx, sub_psms) in pairs(irt_bin_groups)
        last_pair_id = assignPairIds!(sub_psms, last_pair_id)
    end
end


function assignPairIds!(psms::AbstractDataFrame, last_pair_id::UInt32)
    psms[!,:pair_id], last_pair_id = assign_pair_ids(
        psms.target, psms.decoy, psms.precursor_idx, psms.irt_bin_idx, last_pair_id
    )
    return last_pair_id
end

"""
    assign_pair_ids(target, decoy, precursor_idx, irt_bin_idx, last_pair_id)

Randomly assign pair IDs to precursors within an iRT bin, allowing all pairing types.

# Arguments
- `target::AbstractVector{Bool}`: Target/decoy labels (used only for statistics)
- `decoy::AbstractVector{Bool}`: Decoy flags (used only for statistics)
- `precursor_idx::AbstractVector{UInt32}`: Unique precursor identifiers
- `irt_bin_idx::AbstractVector{UInt32}`: iRT bin identifiers
- `last_pair_id::UInt32`: Last assigned pair_id (for continuity across bins)

# Returns
- `pair_ids::Vector{UInt32}`: Pair ID for each row
- `last_pair_id::UInt32`: Updated last pair_id value

# Pairing Strategy
Within each iRT bin, all unique precursors are randomly shuffled and paired
consecutively: (1,2), (3,4), (5,6), etc. This creates:
- **Target-Target pairs**: Both members are targets
- **Target-Decoy pairs**: One target, one decoy
- **Decoy-Decoy pairs**: Both members are decoys
- **Singleton pairs**: If odd number, last precursor pairs with itself

# Notes
- Uses fixed random seed (`PAIRING_RANDOM_SEED`) for reproducibility
- Both members of a pair receive the same `pair_id`
- Pairing is independent of target/decoy status
- Reports pairing statistics to @debug_l2 log level
"""
function assign_pair_ids(
    target::AbstractVector{Bool}, decoy::AbstractVector{Bool},
    precursor_idx::AbstractVector{UInt32}, irt_bin_idx::AbstractVector{UInt32},
    last_pair_id::UInt32
)
    # Get all unique precursors in this iRT bin (regardless of target/decoy status)
    unique_precursors = unique(precursor_idx)
    n_precursors = length(unique_precursors)

    # Randomly shuffle all precursors using fixed seed for reproducibility
    perm = randperm(MersenneTwister(PAIRING_RANDOM_SEED), n_precursors)
    shuffled_precursors = unique_precursors[perm]

    # Create mapping from precursor_idx to pair_id
    precursor_idx_to_pair_id = Dict{UInt32, UInt32}()

    # Pair consecutive elements: (1,2), (3,4), (5,6), etc.
    n_pairs = n_precursors ÷ 2
    for i in 1:n_pairs
        last_pair_id += one(UInt32)
        precursor_idx_to_pair_id[shuffled_precursors[2*i - 1]] = last_pair_id
        precursor_idx_to_pair_id[shuffled_precursors[2*i]] = last_pair_id
    end

    # Handle odd number - last precursor gets singleton pair_id
    if isodd(n_precursors)
        last_pair_id += one(UInt32)
        precursor_idx_to_pair_id[shuffled_precursors[end]] = last_pair_id
    end

    # Count pairing types for reporting
    n_target_target = 0
    n_target_decoy = 0
    n_decoy_decoy = 0
    n_singleton = isodd(n_precursors) ? 1 : 0

    # Build target/decoy lookup for statistics
    precursor_is_target = Dict{UInt32, Bool}()
    for (i, pid) in enumerate(precursor_idx)
        if !haskey(precursor_is_target, pid)
            precursor_is_target[pid] = target[i]
        end
    end

    # Count pair types
    for i in 1:n_pairs
        p1 = shuffled_precursors[2*i - 1]
        p2 = shuffled_precursors[2*i]
        is_t1 = precursor_is_target[p1]
        is_t2 = precursor_is_target[p2]

        if is_t1 && is_t2
            n_target_target += 1
        elseif !is_t1 && !is_t2
            n_decoy_decoy += 1
        else
            n_target_decoy += 1
        end
    end

    # Report pairing statistics for this bin
    @debug_l2 "iRT bin $(first(irt_bin_idx)): $(n_pairs) pairs (T-T: $n_target_target, T-D: $n_target_decoy, D-D: $n_decoy_decoy), Singletons: $n_singleton"

    # Map all rows to their pair_ids
    pair_ids = similar(precursor_idx, UInt32)
    for (row_idx, pid) in enumerate(precursor_idx)
        pair_ids[row_idx] = precursor_idx_to_pair_id[pid]
    end

    return pair_ids, last_pair_id
end

function sort_of_percolator_in_memory!(psms::DataFrame,
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_lightgbm_rescore::Float32 = 0.01f0,
                  max_q_value_mbr_itr::Float32 = 0.20f0,
                  min_PEP_neg_threshold_itr = 0.90f0,
                  feature_fraction::Float64 = 0.5,
                  learning_rate::Float64 = 0.15,
                  min_data_in_leaf::Int = 1,
                  bagging_fraction::Float64 = 0.5,
                  min_gain_to_split::Float64 = 0.0,
                  max_depth::Int = 10,
                  num_leaves::Int = 63,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = false,
                  show_progress::Bool = true,
                  verbose_logging::Bool = false)
    
    # Apply random target-decoy pairing before ML training
    assign_random_target_decoy_pairs!(psms)
    #Faster if sorted first (handle missing pair_id values)
    sort!(psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])
    # Display target/decoy/entrapment counts for training dataset
    if verbose_logging
        n_targets = sum(psms.target)
        n_decoys = sum(.!psms.target)
        n_entrapments = hasproperty(psms, :entrapment) ? sum(psms.entrapment) : 0
        n_total = nrow(psms)

        if n_entrapments > 0
            @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys, $n_entrapments entrapments (total: $n_total PSMs)"
        else
            @user_info "ML Training Dataset: $n_targets targets, $n_decoys decoys (total: $n_total PSMs)"
        end
    end

    prob_test   = zeros(Float32, nrow(psms))  # final CV predictions
    prob_train  = zeros(Float32, nrow(psms))  # temporary, used during training
    #first_pass_estimates = zeros(Float32, nrow(psms)) # store first iteration estimates
    MBR_estimates = zeros(Float32, nrow(psms)) # optional MBR layer
    nonMBR_estimates  = zeros(Float32, nrow(psms)) # keep track of last nonMBR test scores

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, LightGBMModelVector}()
    mbr_start_iter = length(iter_scheme)
    iterations_per_fold = match_between_runs ? length(iter_scheme) : max(mbr_start_iter - 1, 1)

    cv_fold_col = psms[!, :cv_fold]
    fold_indices = Dict(fold => findall(==(fold), cv_fold_col) for fold in unique_cv_folds)
    train_indices = Dict(fold => findall(!=(fold), cv_fold_col) for fold in unique_cv_folds)

    Random.seed!(1776)
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    total_progress_steps = length(unique_cv_folds) * iterations_per_fold
    pbar = show_progress ? ProgressBar(total=total_progress_steps) : nothing

    # Collect final-iteration training sets if requested
    final_train_parts = Vector{DataFrame}()

    for test_fold_idx in unique_cv_folds

        initialize_prob_group_features!(psms, match_between_runs)

        train_idx = train_indices[test_fold_idx]
        test_idx  = fold_indices[test_fold_idx]
        
        psms_train = @view psms[train_idx, :]
        psms_test  = @view psms[test_idx, :]

        # Display counts for this CV fold
        if verbose_logging
            n_train_targets = sum(psms_train.target)
            n_train_decoys = sum(.!psms_train.target)
            n_test_targets = sum(psms_test.target)
            n_test_decoys = sum(.!psms_test.target)

            if hasproperty(psms, :entrapment)
                n_train_entrapments = sum(psms_train.entrapment)
                n_test_entrapments = sum(psms_test.entrapment)
                @user_info "Fold $test_fold_idx - Train: $n_train_targets targets, $n_train_decoys decoys, $n_train_entrapments entrapments | Test: $n_test_targets targets, $n_test_decoys decoys, $n_test_entrapments entrapments"
            else
                @user_info "Fold $test_fold_idx - Train: $n_train_targets targets, $n_train_decoys decoys | Test: $n_test_targets targets, $n_test_decoys decoys"
            end
        end

        fold_models = LightGBMModelVector(undef, length(iter_scheme))

        for (itr, num_round) in enumerate(iter_scheme)
            psms_train_itr = get_training_data_for_iteration!(psms_train,
                                                              itr,
                                                              match_between_runs,
                                                              max_q_value_lightgbm_rescore,
                                                              max_q_value_mbr_itr,
                                                              min_PEP_neg_threshold_itr,
                                                              itr >= mbr_start_iter)

            train_feats = itr < mbr_start_iter ? non_mbr_features : features

            bst = train_booster(psms_train_itr, train_feats, num_round;
                               feature_fraction=feature_fraction,
                               learning_rate=learning_rate,
                               min_data_in_leaf=min_data_in_leaf,
                               bagging_fraction=bagging_fraction,
                               min_gain_to_split=min_gain_to_split,
                               max_depth=max_depth,
                               num_leaves=num_leaves)
                               
            fold_models[itr] = bst

            # Print feature importances for each iteration and fold
            if print_importance
                importances = lightgbm_feature_importances(bst)
                if importances === nothing
                    @user_warn "LightGBM backend did not provide feature importances for iteration $(itr)."
                else
                    feature_pairs = collect(zip(bst.features, importances))
                    # Sort by importance in descending order
                    sort!(feature_pairs, by=x->x[2], rev=true)
                    @user_info "Feature Importances - Fold $(test_fold_idx), Iteration $(itr) ($(length(feature_pairs)) features):"
                    for i in 1:10:length(feature_pairs)
                        chunk = feature_pairs[i:min(i+9, end)]
                        feat_strs = ["$(feat):$(round(score, digits=3))" for (feat, score) in chunk]
                        @user_info "  " * join(feat_strs, " | ")
                    end
                end
            end

            #predict_fold!(bst, psms_train, psms_test, train_feats)
            # **temporary predictions for training only**
            prob_train[train_idx] = predict(bst, psms_train)
            psms_train[!,:prob] = prob_train[train_idx]
            get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)

            # **predict held-out fold**
            prob_test[test_idx] = predict(bst, psms_test)
            psms_test[!,:prob] = prob_test[test_idx]

            #if itr == 1
            #    first_pass_estimates[test_idx] = prob_test[test_idx]
            #end

            if itr == (mbr_start_iter - 1)
			    nonMBR_estimates[test_idx] = prob_test[test_idx]
            end

            if match_between_runs
                update_mbr_features!(psms_train, psms_test, prob_test,
                                     test_idx, itr, mbr_start_iter,
                                     max_q_value_lightgbm_rescore)
            end

            show_progress && update(pbar)

            if (!match_between_runs) && itr == (mbr_start_iter - 1)
                break
            end
        end
        # Make predictions on hold out data.
        if match_between_runs
            MBR_estimates[test_idx] = psms_test.prob
        else
            prob_test[test_idx] = psms_test.prob
        end
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end

    if match_between_runs
        # Determine which precursors failed the q-value cutoff prior to MBR
        qvals_prev = Vector{Float32}(undef, length(nonMBR_estimates))
        get_qvalues!(nonMBR_estimates, psms.target, qvals_prev)
        pass_mask = (qvals_prev .<= max_q_value_lightgbm_rescore)
        has_passing_psms = !isempty(pass_mask) && any(pass_mask)
        prob_thresh = has_passing_psms ? minimum(nonMBR_estimates[pass_mask]) : typemax(Float32)
        # Label as transfer candidates only those failing the q-value cutoff but
        # whose best matched pair surpassed the passing probability threshold.
        psms[!, :MBR_transfer_candidate] .= .!pass_mask .&
                                            (psms.MBR_max_pair_prob .>= prob_thresh)

        # Store base trace probabilities (non-MBR)
        psms[!, :trace_prob] = nonMBR_estimates
        # Store MBR-enhanced trace probabilities
        psms[!, :MBR_boosted_trace_prob] = MBR_estimates
        #psms[!,:first_pass_prob] = first_pass_estimates
    else
        # Only store base trace probabilities
        psms[!, :trace_prob] = prob_test
        # Do NOT create MBR_boosted_trace_prob column
        #psms[!,:first_pass_prob] = first_pass_estimates
    end

    return models
end

function sort_of_percolator_out_of_memory!(psms::DataFrame,
                    file_paths::Vector{String},
                    features::Vector{Symbol},
                    match_between_runs::Bool = true;
                    max_q_value_lightgbm_rescore::Float32 = 0.01f0,
                    max_q_value_mbr_itr::Float32 = 0.20f0,
                    min_PEP_neg_threshold_itr::Float32 = 0.90f0,
                    feature_fraction::Float64 = 0.5,
                    learning_rate::Float64 = 0.15,
                    min_data_in_leaf::Int = 1,
                    bagging_fraction::Float64 = 0.5,
                    min_gain_to_split::Float64 = 0.0,
                    max_depth::Int = 10,
                    num_leaves::Int = 63,
                    iter_scheme::Vector{Int} = [100, 200, 200],
                    print_importance::Bool = false)

    function getBestScorePerPrec!(
        prec_to_best_score_new::Dictionary,
        file_paths::Vector{String},
        models::Dictionary{UInt8,LightGBMModel},
        non_mbr_models::Union{Dictionary{UInt8,LightGBMModel}, Nothing},
        features::Vector{Symbol},
        non_mbr_features::Vector{Symbol},
        match_between_runs::Bool;
        is_last_iteration::Bool = false)
    
        # Reset counts for new scores
        reset_precursor_scores!(prec_to_best_score_new)
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            
            probs = predict_cv_models(models, psms_subset, features)
            
            if match_between_runs && !is_last_iteration
                #Update maximum probabilities for tracked precursors 
                qvals = zeros(Float32, nrow(psms_subset))
                get_qvalues!(probs, psms_subset.target, qvals)

                for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                    prob = probs[i]
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        # Update running statistics with new probability
                        updated_stats = update_pair_statistics(scores, prob)

                        if prob > scores.best_prob_1
                           new_scores = merge(updated_stats, (
                                # replace best_prob_2 with best_prob_1
                                best_prob_2                     = scores.best_prob_1,
                                best_log2_weights_2             = scores.best_log2_weights_1,
                                best_irts_2                     = scores.best_irts_1,
                                best_irt_residual_2             = scores.best_irt_residual_1,
                                best_weight_2                   = scores.best_weight_1,
                                best_log2_intensity_explained_2 = scores.best_log2_intensity_explained_1,
                                best_ms_file_idx_2              = scores.best_ms_file_idx_1,
                                is_best_decoy_2                 = scores.is_best_decoy_1,
                                # overwrite best_prob_1
                                best_prob_1                     = prob,
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_irts_1                     = psms_subset.irts[i],
                                best_irt_residual_1             = irt_residual(psms_subset, i),
                                best_weight_1                   = psms_subset.weight[i],
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                is_best_decoy_1                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores

                        elseif prob > scores.best_prob_2
                            # overwrite best_prob_2
                            new_scores = merge(updated_stats, (
                                best_prob_2                     = prob,
                                best_log2_weights_2             = log2.(psms_subset.weights[i]),
                                best_irts_2                     = psms_subset.irts[i],
                                best_irt_residual_2             = irt_residual(psms_subset, i),
                                best_weight_2                   = psms_subset.weight[i],
                                best_log2_intensity_explained_2 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_2              = psms_subset.ms_file_idx[i],
                                is_best_decoy_2                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores
                        else
                            # No change to best/second best, but update running stats
                            prec_to_best_score_new[key] = updated_stats
                        end

                        if qvals[i] <= max_q_value_lightgbm_rescore
                            push!(scores.unique_passing_runs, psms_subset.ms_file_idx[i])
                        end

                    else
                        insert!(prec_to_best_score_new, key, (
                                best_prob_1                     = prob,
                                best_prob_2                     = zero(Float32),
                                worst_prob_1                    = prob,
                                worst_prob_2                    = zero(Float32),
                                mean_prob                       = prob,
                                count_pairs                     = Int32(1),
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_log2_weights_2             = Vector{Float32}(),
                                best_irts_1                     = psms_subset.irts[i],
                                best_irts_2                     = Vector{Float32}(),
                                best_irt_residual_1             = irt_residual(psms_subset, i),
                                best_irt_residual_2             = zero(Float32),
                                best_weight_1                   = psms_subset.weight[i],
                                best_weight_2                   = zero(Float32),
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_log2_intensity_explained_2 = zero(Float32),
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                best_ms_file_idx_2              = zero(UInt32),
                                is_best_decoy_1                 = psms_subset.decoy[i],
                                is_best_decoy_2                 = false,
                                unique_passing_runs             = ( qvals[i] <= max_q_value_lightgbm_rescore ?
                                                                    Set{UInt16}([psms_subset.ms_file_idx[i]]) :
                                                                    Set{UInt16}() )
                            ))
                    end
                end
            end


            if is_last_iteration
                if match_between_runs
                    update_mbr_probs!(psms_subset, probs, max_q_value_lightgbm_rescore)
                else
                    psms_subset.prob = probs
                end

            end
        end
    

        # Compute probs and features for next round
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = predict_cv_models(models, psms_subset, features)

            # On last iteration with MBR, also get non-MBR predictions
            non_mbr_probs = if is_last_iteration && match_between_runs && non_mbr_models !== nothing
                predict_cv_models(non_mbr_models, psms_subset, non_mbr_features)
            else
                nothing
            end

            for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                psms_subset[i,:prob] = probs[i]


                if match_between_runs && !is_last_iteration
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        psms_subset.MBR_num_runs[i] = length(scores.unique_passing_runs)

                        best_log2_weights = Float32[]
                        best_irts = Float32[]
                        best_weight = zero(Float32)
                        best_log2_ie = zero(Float32)
                        best_residual = zero(Float32)

                        if (scores.best_ms_file_idx_1 != psms_subset.ms_file_idx[i]) &&
                           (!isempty(scores.best_log2_weights_1))
                            best_log2_weights                   = scores.best_log2_weights_1
                            best_irts                           = scores.best_irts_1
                            best_weight                         = scores.best_weight_1
                            best_log2_ie                        = scores.best_log2_intensity_explained_1
                            best_residual                       = scores.best_irt_residual_1
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_1
                            MBR_is_best_decoy                   = scores.is_best_decoy_1
                        elseif (scores.best_ms_file_idx_2 != psms_subset.ms_file_idx[i]) &&
                               (!isempty(scores.best_log2_weights_2))
                            best_log2_weights                   = scores.best_log2_weights_2
                            best_irts                           = scores.best_irts_2
                            best_weight                         = scores.best_weight_2
                            best_log2_ie                        = scores.best_log2_intensity_explained_2
                            best_residual                       = scores.best_irt_residual_2
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_2
                            MBR_is_best_decoy                   = scores.is_best_decoy_2
                        else
                            psms_subset.MBR_best_irt_diff[i]        = -1.0f0
                            psms_subset.MBR_rv_coefficient[i]       = -1.0f0
                            psms_subset.MBR_is_best_decoy[i]        = true
                            psms_subset.MBR_max_pair_prob[i]        = -1.0f0
                            psms_subset.MBR_log2_weight_ratio[i]    = -1.0f0
                            psms_subset.MBR_log2_explained_ratio[i] = -1.0f0
                            psms_subset.MBR_is_missing[i]           = true
                            continue
                        end

                        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(psms_subset.weights[i]))
                        best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_irts, psms_subset.irts[i])

                        current_residual = irt_residual(psms_subset, i)
                        psms_subset.MBR_best_irt_diff[i] = abs(best_residual - current_residual)
                        psms_subset.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
                        psms_subset.MBR_log2_weight_ratio[i] = log2(psms_subset.weight[i] / best_weight)
                        psms_subset.MBR_log2_explained_ratio[i] = psms_subset.log2_intensity_explained[i] - best_log2_ie
                        psms_subset.MBR_is_best_decoy[i] = MBR_is_best_decoy
                    end
                end
            end

            # On last iteration with MBR: probs=MBR-boosted, non_mbr_probs=base
            # Otherwise: probs=base, non_mbr_probs=nothing
            mbr_probs_arg = if is_last_iteration && match_between_runs
                probs  # MBR-boosted predictions
            else
                nothing
            end
            probs_arg = if is_last_iteration && match_between_runs && non_mbr_probs !== nothing
                non_mbr_probs  # Non-MBR predictions
            else
                probs  # Regular predictions
            end

            write_subset(
                file_path,
                psms_subset,
                probs_arg,
                mbr_probs_arg,
                match_between_runs,
                max_q_value_lightgbm_rescore;
                dropVectors = is_last_iteration,
            )
        end
        
        return prec_to_best_score_new
    end


    unique_cv_folds = unique(psms[!, :cv_fold])
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, LightGBMModelVector}()
    mbr_start_iter = length(iter_scheme)
    iterations_per_fold = match_between_runs ? length(iter_scheme) : max(mbr_start_iter - 1, 1)
    total_progress_steps = length(unique_cv_folds) * iterations_per_fold
    pbar = ProgressBar(total=total_progress_steps)
    Random.seed!(1776);
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        #Clear prob stats 
        initialize_prob_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])

        for (itr, num_round) in enumerate(iter_scheme)

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_lightgbm_rescore,
                                                                max_q_value_mbr_itr,
                                                                min_PEP_neg_threshold_itr,
                                                                itr >= length(iter_scheme))
            ###################
            #Train a model on the n-1 training folds.
            train_feats = itr < length(iter_scheme) ? non_mbr_features : features
            bst = train_booster(psms_train_itr, train_feats, num_round;
                               feature_fraction=feature_fraction,
                               learning_rate=learning_rate,
                               min_data_in_leaf=min_data_in_leaf,
                               bagging_fraction=bagging_fraction,
                               min_gain_to_split=min_gain_to_split,
                               max_depth=max_depth,
                               num_leaves=num_leaves)
            if !haskey(models, test_fold_idx)
                insert!(
                    models,
                    test_fold_idx,
                    LightGBMModelVector([bst])
                )
            else
                push!(models[test_fold_idx], bst)
            end
            # Print feature importances
            if print_importance
                importances = lightgbm_feature_importances(bst)
                if importances === nothing
                    @user_warn "LightGBM backend did not provide feature importances for iteration $(itr)."
                else
                    feature_pairs = collect(zip(bst.features, importances))
                    @user_info "Feature Importances ($(length(feature_pairs)) features):"
                    for i in 1:10:length(feature_pairs)
                        chunk = feature_pairs[i:min(i+9, end)]
                        feat_strs = ["$(feat):$(round(score, digits=3))" for (feat, score) in chunk]
                        @user_info "  " * join(feat_strs, " | ")
                    end
                end
            end

            # Get probabilities for training sample so we can get q-values
            psms_train[!,:prob] = lightgbm_predict(bst, psms_train; output_type=Float32)
            
            if match_between_runs
                summarize_precursors!(psms_train, q_cutoff = max_q_value_lightgbm_rescore)
            end

            show_progress && update(pbar)
        end
    end
    
    pbar = ProgressBar(total=length(iter_scheme))
    prec_to_best_score = Dictionary{@NamedTuple{pair_id::UInt32,
                                                isotopes::Tuple{Int8,Int8}},
                                    @NamedTuple{best_prob_1::Float32,
                                                best_prob_2::Float32,
                                                worst_prob_1::Float32,
                                                worst_prob_2::Float32,
                                                mean_prob::Float32,
                                                count_pairs::Int32,
                                                best_log2_weights_1::Vector{Float32},
                                                best_log2_weights_2::Vector{Float32},
                                                best_irts_1::Vector{Float32},
                                                best_irts_2::Vector{Float32},
                                                best_irt_residual_1::Float32,
                                                best_irt_residual_2::Float32,
                                                best_weight_1::Float32,
                                                best_weight_2::Float32,
                                                best_log2_intensity_explained_1::Float32,
                                                best_log2_intensity_explained_2::Float32,
                                                best_ms_file_idx_1::UInt32,
                                                best_ms_file_idx_2::UInt32,
                                                is_best_decoy_1::Bool,
                                                is_best_decoy_2::Bool,
                                                unique_passing_runs::Set{UInt16}}}()

    for (train_iter, num_round) in enumerate(iter_scheme)
        models_for_iter = Dictionary{UInt8,LightGBMModel}()
        for test_fold_idx in unique_cv_folds
            insert!(models_for_iter, test_fold_idx, models[test_fold_idx][train_iter])
        end

        # On last iteration with MBR, also get previous (non-MBR) models
        non_mbr_models_for_iter = if (train_iter == length(iter_scheme)) && match_between_runs && (train_iter > 1)
            non_mbr_models = Dictionary{UInt8,LightGBMModel}()
            for test_fold_idx in unique_cv_folds
                insert!(non_mbr_models, test_fold_idx, models[test_fold_idx][train_iter - 1])
            end
            non_mbr_models
        else
            nothing
        end

        # Determine which features to use
        current_features = (train_iter == length(iter_scheme)) ? features : non_mbr_features

        prec_to_best_score = getBestScorePerPrec!(
            prec_to_best_score,
            file_paths,
            models_for_iter,
            non_mbr_models_for_iter,
            current_features,
            non_mbr_features,
            match_between_runs;
            is_last_iteration = (train_iter == length(iter_scheme))
        )
        update(pbar)
    end

    return models
end

function train_booster(psms::AbstractDataFrame, features, num_round;
                       feature_fraction::Float64,
                       learning_rate::Float64,
                       min_data_in_leaf::Int,
                       bagging_fraction::Float64,
                       min_gain_to_split::Float64,
                       max_depth::Int,
                       num_leaves::Int)

    classifier = build_lightgbm_classifier(
        num_iterations = num_round,
        max_depth = max_depth,
        learning_rate = learning_rate,
        num_leaves = num_leaves,
        feature_fraction = feature_fraction,
        bagging_fraction = bagging_fraction,
        bagging_freq = bagging_fraction < 1 ? 1 : 0,
        min_data_in_leaf = min_data_in_leaf,
        min_gain_to_split = min_gain_to_split,
    )
    feature_frame = psms[:, features]
    return fit_lightgbm_model(classifier, feature_frame, psms.target; positive_label=true)
end

function predict_fold!(bst, psms_train::AbstractDataFrame,
                       psms_test::AbstractDataFrame, features)
    psms_test[!, :prob] = lightgbm_predict(bst, psms_test; output_type=Float32)
    psms_train[!, :prob] = lightgbm_predict(bst, psms_train; output_type=Float32)
    get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)
end

function update_mbr_features!(psms_train::AbstractDataFrame,
                              psms_test::AbstractDataFrame,
                              prob_test::Vector{Float32},
                              test_fold_idxs,
                              itr::Int,
                              mbr_start_iter::Int,
                              max_q_value_lightgbm_rescore::Float32)
    if itr >= mbr_start_iter - 1
        get_qvalues!(psms_test.prob, psms_test.target, psms_test.q_value)
        summarize_precursors!(psms_test, q_cutoff = max_q_value_lightgbm_rescore)
        summarize_precursors!(psms_train, q_cutoff = max_q_value_lightgbm_rescore)
    end
    if itr == mbr_start_iter - 1
        prob_test[test_fold_idxs] = psms_test.prob
    end
end

function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)
    # Diagnostic: Show isotope and pairing interaction
    n_unique_pairs = length(unique(psms.pair_id))
    unique_isotopes = unique(psms.isotopes_captured)
    n_unique_isotopes = length(unique_isotopes)

    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
    n_pair_isotope_groups = length(pair_groups)

    @debug_l2 "MBR Feature Computation: $n_unique_pairs unique pair_ids × $n_unique_isotopes isotope combinations = $n_pair_isotope_groups groups"
    @debug_l2 "Isotope combinations present: $unique_isotopes"

    Threads.@threads for idx in eachindex(pair_groups)
        _, sub_psms = pair_groups[idx]
        
        # Efficient way to find the top 2 precursors so we can do MBR on the 
        # best precursor match that isn't itself. It's always one of the top 2.

        # single pass: record the best PSM index & prob per run
        offset = Int(minimum(sub_psms.ms_file_idx))
        range_len = Int(maximum(sub_psms.ms_file_idx)) - offset + 1
        best_i = zeros(Int, range_len)
        best_p = fill(-Inf, range_len)
        for (i, run) in enumerate(sub_psms.ms_file_idx)
            idx = Int(run) - offset + 1
            p = sub_psms.prob[i]
            if p > best_p[idx]
                best_p[idx] = p
                best_i[idx] = i
            end
        end

        # if more than one run, find the global top-2 runs by their best-PSM prob
        run_best_indices = zeros(Int, range_len)
        runs = findall(!=(0), best_i)
        if length(runs) > 1
            # track top two runs (r1 > r2)
            r1 = 0; p1 = -Inf
            r2 = 0; p2 = -Inf
            for r in runs
                p = best_p[r]
                if p > p1
                    r2, p2 = r1, p1
                    r1, p1 = r, p
                elseif p > p2
                    r2, p2 = r, p
                end
            end

            # assign, for each run, the best index in “any other” run
            for r in runs
                run_best_indices[r] = (r == r1 ? best_i[r2] : best_i[r1])
            end
        end

        # Compute MBR features
        num_runs_passing = length(sub_psms.ms_file_idx[sub_psms.q_value .<= q_cutoff])
        for i in 1:nrow(sub_psms)
            sub_psms.MBR_num_runs[i] = num_runs_passing - (sub_psms.q_value[i] .<= q_cutoff)

            idx = Int(sub_psms.ms_file_idx[i]) - offset + 1
            best_idx = run_best_indices[idx]
            if best_idx == 0 || sub_psms.MBR_num_runs[i] == 0
                sub_psms.MBR_best_irt_diff[i]           = -1.0f0
                sub_psms.MBR_rv_coefficient[i]          = -1.0f0
                sub_psms.MBR_is_best_decoy[i]           = true
                sub_psms.MBR_log2_weight_ratio[i]       = -1.0f0
                sub_psms.MBR_log2_explained_ratio[i]    = -1.0f0
                sub_psms.MBR_max_pair_prob[i]           = -1.0f0
                sub_psms.MBR_is_missing[i]              = true
                continue
            end

            best_log2_weights = log2.(sub_psms.weights[best_idx])
            best_iRTs = sub_psms.irts[best_idx]
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(sub_psms.weights[i]))
            best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_iRTs, sub_psms.irts[i])

            sub_psms.MBR_max_pair_prob[i] = sub_psms.prob[best_idx]
            best_residual = irt_residual(sub_psms, best_idx)
            current_residual = irt_residual(sub_psms, i)
            sub_psms.MBR_best_irt_diff[i] = abs(best_residual - current_residual)
            sub_psms.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
            sub_psms.MBR_log2_weight_ratio[i] = log2(sub_psms.weight[i] / sub_psms.weight[best_idx])
            sub_psms.MBR_log2_explained_ratio[i] = sub_psms.log2_intensity_explained[i] - sub_psms.log2_intensity_explained[best_idx]
            sub_psms.MBR_is_best_decoy[i] = sub_psms.decoy[best_idx]
        end
    end
end

function initialize_prob_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    n = nrow(psms)
    psms[!, :prob]      = zeros(Float32, n)
    psms[!, :q_value]   = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
        psms[!, :MBR_best_irt_diff]             = zeros(Float32, n)
        psms[!, :MBR_log2_weight_ratio]         = zeros(Float32, n)
        psms[!, :MBR_log2_explained_ratio]      = zeros(Float32, n)
        psms[!, :MBR_rv_coefficient]            = zeros(Float32, n)
        psms[!, :MBR_is_best_decoy]             = trues(n)
        psms[!, :MBR_num_runs]                  = zeros(Int32, n)
        psms[!, :MBR_transfer_candidate]        = falses(n)
        psms[!, :MBR_is_missing]                = falses(n)
    end

    return psms
end

function get_training_data_for_iteration!(
    psms_train::AbstractDataFrame,
    itr::Int,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    last_iter::Bool
)
   
    if itr == 1
        # Train on all precursors during first iteration. 
        return copy(psms_train)
    else
        # Do a shallow copy to avoid overwriting target/decoy labels
        psms_train_itr = copy(psms_train)

        # Convert the worst-scoring targets to negatives using PEP estimate
        order = sortperm(psms_train_itr.prob, rev=true)
        sorted_scores  = psms_train_itr.prob[order]
        sorted_targets = psms_train_itr.target[order]
        PEPs = Vector{Float32}(undef, length(order))
        get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)
        idx_cutoff = findfirst(x -> x >= min_PEP_neg_threshold_itr, PEPs)
        if !isnothing(idx_cutoff)
            worst_idxs = order[idx_cutoff:end]
            psms_train_itr.target[worst_idxs] .= false
        end

        # Take all decoys and targets passing q_thresh
        psms_train_itr = subset(
            psms_train_itr,
            [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_lightgbm_rescore))
        )

        return psms_train_itr
    end
end

function dropVectorColumns!(df)
    to_drop = String[]
    for col in names(df)
        if eltype(df[!, col]) <: AbstractVector
            push!(to_drop, col)
        end
    end
    # 2) Drop those columns in place
    select!(df, Not(to_drop))
end

"""
    reset_precursor_scores!(dict)

Set all values of `dict` to an empty precursor score tuple.  This helps reuse
the same dictionary between iterations without reallocating memory.
"""
function reset_precursor_scores!(dict)
    for key in keys(dict)
        dict[key] = (
            best_prob_1 = zero(Float32),
            best_prob_2 = zero(Float32),
            best_log2_weights_1 = Vector{Float32}(),
            best_log2_weights_2 = Vector{Float32}(),
            best_irts_1 = Vector{Float32}(),
            best_irts_2 = Vector{Float32}(),
            best_irt_residual_1 = zero(Float32),
            best_irt_residual_2 = zero(Float32),
            best_weight_1 = zero(Float32),
            best_weight_2 = zero(Float32),
            best_log2_intensity_explained_1 = zero(Float32),
            best_log2_intensity_explained_2 = zero(Float32),
            best_ms_file_idx_1 = zero(UInt32),
            best_ms_file_idx_2 = zero(UInt32),
            is_best_decoy_1 = false,
            is_best_decoy_2 = false,
            unique_passing_runs = Set{UInt16}(),
        )
    end
    return dict
end

"""
    predict_cv_models(models, df, features)

Return a vector of probabilities for `df` using the cross validation `models`.
"""
function predict_cv_models(models::Dictionary{UInt8,LightGBMModel},
                           df::AbstractDataFrame,
                           features::Vector{Symbol})
    probs = zeros(Float32, nrow(df))
    for (fold_idx, bst) in pairs(models)
        fold_rows = findall(==(fold_idx), df[!, :cv_fold])
        if !isempty(fold_rows)
            probs[fold_rows] = lightgbm_predict(bst, df[fold_rows, :]; output_type=Float32)
        end
    end
    return probs
end

"""
    update_mbr_probs!(df, probs, qval_thresh)

Store final MBR probabilities and mark transfer candidates as those
failing the pre-MBR q-value threshold but whose best matched pair passed
the corresponding probability cutoff.
"""
function update_mbr_probs!(
    df::AbstractDataFrame,
    probs::AbstractVector{Float32},
    mbr_probs::AbstractVector{Float32},
    qval_thresh::Float32,
)
    prev_qvals = similar(probs)
    get_qvalues!(probs, df.target, prev_qvals)
    pass_mask = (prev_qvals .<= qval_thresh) .& df.target
    prob_thresh = any(pass_mask) ? minimum(probs[pass_mask]) : typemax(Float32)
    df[!, :MBR_transfer_candidate] = (prev_qvals .> qval_thresh) .&
                                     (df.MBR_max_pair_prob .>= prob_thresh)
    df[!, :trace_prob] = probs
    df[!, :MBR_boosted_trace_prob] = mbr_probs
    return df
end

"""
    write_subset(file_path, df, probs, match_between_runs, qval_thresh; dropVectors=false)

Write the updated subset to disk, optionally dropping vector columns.
The `qval_thresh` argument is used to mark transfer candidates when
`match_between_runs` is true and `dropVectors` is set.
"""
function write_subset(
    file_path::String,
    df::DataFrame,
    probs::AbstractVector{Float32},
    mbr_probs::Union{AbstractVector{Float32}, Nothing},
    match_between_runs::Bool,
    qval_thresh::Float32;
    dropVectors::Bool=false,
)
    if dropVectors
        if match_between_runs
            update_mbr_probs!(df, probs, mbr_probs, qval_thresh)
        else
            df[!, :trace_prob] = probs
        end
        writeArrow(file_path, dropVectorColumns!(df))
    else
        if match_between_runs && mbr_probs !== nothing
            df[!, :trace_prob] = probs
            df[!, :MBR_boosted_trace_prob] = mbr_probs
        else
            df[!, :trace_prob] = probs
        end
        writeArrow(file_path, convert_subarrays(df))
    end
end

function MBR_rv_coefficient(weights_A::AbstractVector{<:Real},
    times_A::AbstractVector{<:Real},
    weights_B::AbstractVector{<:Real},
    times_B::AbstractVector{<:Real})

    # Construct two Nx2 matrices, each row is (weight, time)
    X = hcat(collect(weights_A), collect(times_A))
    Y = hcat(collect(weights_B), collect(times_B))

    # Compute cross-products (Gram matrices)
    Sx = X' * X
    Sy = Y' * Y

    # Numerator: trace(Sx * Sy)
    numerator = tr(Sx * Sy)

    # Denominator: sqrt( trace(Sx*Sx)* trace(Sy*Sy) )
    denominator = sqrt(tr(Sx * Sx) * tr(Sy * Sy))

    # Protect against zero in denominator (e.g. if X or Y is all zeros)
    if denominator == 0
        return 0.0
    end

    return numerator / denominator
end

function pad_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, d)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        return padded_data
    end

    d = length(y) - length(x)
    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, d), y)
    else
        # Pad y
        return (x, pad(y, abs(d)))
    end
end


function pad_rt_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, n, d, rt_step)
        left  = div(d, 2) 
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        @inbounds @fastmath for i in range(1, left)
            padded_data[i] = padded_data[left+1] - (rt_step * ((left+1) - i))
        end 
        @inbounds @fastmath for i in range(1, right)
            padded_data[i+left+n] = padded_data[left+n] + (rt_step * i)
        end 
        return padded_data
    end

    nx = length(x)
    ny = length(y)
    d = ny - nx
    rt_step_x = (x[end] - x[1]) / nx
    rt_step_y = (y[end] - y[1]) / ny
    rt_step = nx > 1 ? rt_step_x : rt_step_y

    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, nx, d, rt_step), y)
    else
        # Pad y
        return (x, pad(y, ny, abs(d), rt_step))
    end
end

function summarize_prob(probs::AbstractVector{Float32})
    minimum(probs), maximum(probs), mean(probs)
end

function convert_subarrays(df::DataFrame)
    for col_name in names(df)
        col_data = df[!, col_name]
        # If it's a SubArray, let's convert it to a plain vector:
        if eltype(col_data) <: AbstractVector
           df[!, col_name] = [Vector(col_data[i]) for i in eachindex(col_data)]
        end
    end
    return df
end
