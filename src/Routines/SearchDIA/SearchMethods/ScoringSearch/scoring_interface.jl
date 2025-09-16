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
Interface functions for ScoringSearch that work exclusively with file references.

These functions provide a clean abstraction layer between ScoringSearch and file operations,
ensuring all file access goes through the reference system.
"""
#==========================================================
Reference-Based Wrappers for Direct File Operations
==========================================================#

"""
    calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    
Calculate global protein scores and add them to files via references.
Returns the score dictionary for downstream use.
"""
function calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    sqrt_n_runs = max(1, floor(Int, sqrt(length(pg_refs))))
    acc_to_scores = Dict{ProteinKey, Vector{Float32}}()

    # First pass: collect scores per protein across all files
    for ref in pg_refs
        process_with_memory_limit(ref,
            batch -> begin
                for row in eachrow(batch)
                    key = ProteinKey(row.protein_name, row.target, row.entrap_id)
                    push!(get!(acc_to_scores, key, Float32[]), row.pg_score)
                end
            end
        )
    end

    # Compute global score using log-odds combination
    acc_to_global_score = Dict{ProteinKey, Float32}()
    for (key, scores) in acc_to_scores
        acc_to_global_score[key] = logodds(scores, sqrt_n_runs)
    end
    
    # Second pass: add global_pg_score column and sort
    for ref in pg_refs
        add_column_and_sort!(ref, :global_pg_score, 
            batch -> begin
                scores = Vector{Float32}(undef, nrow(batch))
                for i in 1:nrow(batch)
                    key = ProteinKey(
                        batch.protein_name[i],
                        batch.target[i],
                        batch.entrap_id[i]
                    )
                    scores[i] = get(acc_to_global_score, key, batch.pg_score[i])
                end
                scores
            end,
            :global_pg_score, :target;  # sort keys
            reverse=true
        )
    end
    
    return acc_to_global_score
end

"""
    add_trace_qvalues(fdr_scale_factor::Float32)

Add a column `:trace_qval` based on the `:prob` column using target/decoy q-values.
"""
function add_trace_qvalues(fdr_scale_factor::Float32)
    op = function(df)
        qvals = Vector{Float32}(undef, nrow(df))
        get_qvalues!(df.prob, df.target, qvals; fdr_scale_factor=fdr_scale_factor)
        df[!, :trace_qval] = qvals
        return df
    end
    return "add_trace_qvalues" => op
end



"""
    add_prec_prob(prob_col::Symbol)

Compute run-specific precursor probabilities from the given probability column.
"""
function add_prec_prob(prob_col::Symbol)
    op = function(df)
        transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
                   prob_col => (p -> 1.0f0 - 0.000001f0 - exp(sum(log1p.(-p)))) => :prec_prob)
        return df
    end
    return "add_prec_prob" => op
end

"""
    logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}

Combine probabilities using a log-odds average. 
The final value is converted back to a probability via the logistic function.
"""
function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    # Sort descending and select the top n probabilities
    sorted = sort(probs; rev=true)
    selected = sorted[1:n]
    eps = 1f-6
    # Convert to log-odds, clip to avoid Inf or negative contribution
    logodds = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(logodds) / n
    return 1.0f0 / (1 + exp(-avg))
end


#==========================================================
Reference-based Sort and Filter Functions
==========================================================#





"""
    apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference}, 
                        β_fitted::Vector{Float64}, feature_names::Vector{Symbol})
    
Apply probit regression scores to protein group files.
Note: This function is called from utils.jl and needs access to calculate_probit_scores.
"""
function apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference},
                             β_fitted::Vector{Float64},
                             feature_names::Vector{Symbol})
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Calculate probit scores (function from utils.jl)
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted)
            
            # Update scores
            df[!, :old_pg_score] = copy(df.pg_score)
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by new scores
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            return df
        end
    end
end

#==========================================================
Scoring-Specific Pipeline Operations
==========================================================#

function apply_mbr_filter!(
    merged_df::DataFrame,
    params,
    fdr_scale_factor::Float32,
)
    # 1) identify transfer candidates
    candidate_mask = merged_df.MBR_transfer_candidate

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| # T->D
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false)) # D->T
        #merged_df.decoy # D->D or T->D
    )

    # 3) compute threshold using the local bad_mask
    @user_info "merged_df columns:"
    for col_name in names(merged_df)
        @user_info "  $col_name"
    end
    
    # NEW: ML-based filtering
    MIN_MBR_CANDIDATES_FOR_ML = 0#20000  # Hardcoded threshold
    n_candidates = sum(candidate_mask)
    
    if n_candidates >= MIN_MBR_CANDIDATES_FOR_ML
        @user_info "Using ML-based MBR filtering ($n_candidates candidates >= $MIN_MBR_CANDIDATES_FOR_ML threshold)"
        filtered_prob_col = apply_ml_mbr_filter!(
            merged_df, candidate_mask, is_bad_transfer, params
        )
    else
        @user_info "Using threshold-based MBR filtering ($n_candidates candidates < $MIN_MBR_CANDIDATES_FOR_ML threshold)"
        filtered_prob_col = apply_threshold_mbr_filter!(
            merged_df, candidate_mask, is_bad_transfer, params
        )
    end

    merged_df._filtered_prob = filtered_prob_col
    return :_filtered_prob
end

function apply_threshold_mbr_filter!(
    merged_df::DataFrame,
    candidate_mask::AbstractVector{Bool},
    is_bad_transfer::AbstractVector{Bool},
    params
)
    # Original threshold-based method
    τ = get_ftr_threshold(
        merged_df.prob,
        merged_df.target,
        is_bad_transfer,
        params.max_MBR_false_transfer_rate;
        mask = candidate_mask,
    )

    @user_warn "τ is $τ at α = $(params.max_MBR_false_transfer_rate)"
    
    # Apply threshold filtering
    filtered_probs = ifelse.(
        candidate_mask .& (merged_df.prob .< τ),
        0.0f0,
        merged_df.prob,
    )
    
    return filtered_probs
end

function apply_ml_mbr_filter!(
    merged_df::DataFrame,
    candidate_mask::AbstractVector{Bool},
    is_bad_transfer::AbstractVector{Bool},
    params
)
    # Extract training data (only candidates)
    training_data = merged_df[candidate_mask, :]
    y_train = is_bad_transfer[candidate_mask]
    
    @user_info "ML MBR Training: $(nrow(training_data)) candidates, $(sum(y_train)) bad transfers ($(round(mean(y_train)*100, digits=1))%)"
    
    # Feature selection and preprocessing
    feature_cols = select_mbr_features(training_data)
    X_train, all_feature_names = prepare_mbr_features(training_data[:, feature_cols])
    
    # Check if cv_fold column exists
    if !hasproperty(training_data, :cv_fold)
        @user_warn "No cv_fold column found - falling back to threshold method"
        return apply_threshold_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)
    end
    
    # Cross-validation using existing cv_fold assignments
    unique_folds = unique(training_data.cv_fold)
    @user_info "Using $(length(unique_folds)) CV folds for MBR ML training: $unique_folds"
    
    bad_transfer_scores = zeros(Float64, nrow(training_data))
    
    for test_fold in unique_folds
        # Split data by fold
        test_mask = training_data.cv_fold .== test_fold
        train_mask = .!test_mask
        
        X_train_fold = X_train[train_mask, :]
        y_train_fold = y_train[train_mask]
        X_test_fold = X_train[test_mask, :]
        
        @user_info "  Fold $test_fold: $(sum(train_mask)) train, $(sum(test_mask)) test, $(sum(y_train_fold)) train bad transfers"
        
        # Train probit regression model on this fold
        model = train_probit_mbr_model(X_train_fold, y_train_fold, all_feature_names)
        
        # Predict on test fold (out-of-fold predictions)
        test_predictions = predict_probit_mbr(model, X_test_fold, all_feature_names)
        bad_transfer_scores[test_mask] = test_predictions
    end
    
    # Determine threshold based on desired FTR using out-of-fold predictions
    τ = calibrate_ml_threshold(
        bad_transfer_scores, y_train, Float64(params.max_MBR_false_transfer_rate)
    )
    
    n_filtered = sum(bad_transfer_scores .>= τ)
    @user_info "ML MBR Filter: τ = $(round(τ, digits=4)), $n_filtered/$(length(bad_transfer_scores)) candidates filtered ($(round(n_filtered/length(bad_transfer_scores)*100, digits=1))%)"
    
    # Apply filtering to all rows
    filtered_probs = copy(merged_df.prob)
    candidate_indices = findall(candidate_mask)
    
    for (i, idx) in enumerate(candidate_indices)
        if bad_transfer_scores[i] >= τ
            filtered_probs[idx] = 0.0f0
        end
    end
    
    return filtered_probs
end

function select_mbr_features(df::DataFrame)
    # Focused feature set: prob + irt_error + MBR features
    core_features = [:prob]
    
    # iRT error feature
    irt_features = [:irt_error]
    
    # MBR-specific features (if available)
    mbr_features = [
        :MBR_max_pair_prob, :MBR_best_irt_diff, :MBR_rv_coefficient,
        :MBR_log2_weight_ratio, :MBR_log2_explained_ratio, :MBR_num_runs
    ]
    
    # Filter to available columns
    available_features = Symbol[]
    for feature_set in [core_features, irt_features, mbr_features]
        for feature in feature_set
            if hasproperty(df, feature) && !all(ismissing, df[!, feature])
                push!(available_features, feature)
            end
        end
    end
    
    @user_info "Selected $(length(available_features)) focused features for MBR ML filtering: $(available_features)"
    return available_features
end

function prepare_mbr_features(df::DataFrame)
    # Handle missing values and add intercept
    processed_df = copy(df)
    
    # Add intercept column (constant 1.0 for all rows)
    processed_df[!, :intercept] = ones(Float64, nrow(processed_df))
    
    for col in names(processed_df)
        col_data = processed_df[!, col]
        
        if eltype(col_data) <: Union{Missing, <:Number}
            # Replace missing with median for numeric
            non_missing = collect(skipmissing(col_data))
            if !isempty(non_missing)
                median_val = median(non_missing)
                processed_df[!, col] = coalesce.(col_data, median_val)
            else
                # All missing - replace with zero
                processed_df[!, col] = zeros(eltype(non_missing), nrow(processed_df))
            end
        end
    end
    
    # Convert to Matrix for ML training
    return Matrix{Float64}(processed_df), Symbol[:intercept; Symbol.(names(df))]
end

function train_mbr_filter_model(X::Matrix, y::AbstractVector{Bool}, feature_names::Vector{Symbol}, params)
    # Create DataFrame with proper feature names for EvoTrees API
    training_df = DataFrame(X, feature_names)
    training_df[!, :target] = y
    
    # Simple XGBoost configuration for fast training
    config = EvoTreeClassifier(
        loss = :logloss,
        nrounds = 50,        # Fast training
        max_depth = 4,       # Prevent overfitting
        eta = 0.1,           # Conservative learning rate
        rowsample = 0.8,     # Bootstrap for robustness (was subsample)
        colsample = 0.8,     # Feature bagging
        gamma = 1.0          # Regularization
    )
    
    # Handle extreme class imbalance
    pos_rate = mean(y)
    if pos_rate < 0.01 || pos_rate > 0.99
        @user_warn "Extreme class imbalance in MBR training: $(round(pos_rate*100, digits=2))% positive - results may be unreliable"
    end
    
    model = fit(config, training_df; target_name = :target, feature_names = feature_names, verbosity = 0)
    
    # Feature importance logging
    try
        importance = EvoTrees.importance(model)
        top_features = first(importance, min(5, length(importance)))
        @user_info "Top MBR filter features: $top_features"
    catch e
        @user_warn "Could not extract feature importance: $e"
    end
    
    return model
end

function calibrate_ml_threshold(scores::AbstractVector, y_true::AbstractVector{Bool}, target_ftr::Float64)
    # Find threshold that achieves target false transfer rate
    # FTR = bad_transfers_flagged / total_flagged
    
    if isempty(scores) || all(y_true .== false)
        @user_warn "No bad transfers in training data - using maximum score as threshold"
        return isempty(scores) ? 1.0 : maximum(scores)
    end
    
    sorted_indices = sortperm(scores, rev=true)
    sorted_scores = scores[sorted_indices]
    sorted_labels = y_true[sorted_indices]
    
    cumulative_bad = cumsum(sorted_labels)
    cumulative_total = 1:length(sorted_labels)
    
    # False transfer rate = bad_transfers / total_flagged
    ftr = cumulative_bad ./ cumulative_total
    
    # Find first point where FTR <= target
    valid_idx = findfirst(x -> x <= target_ftr, ftr)
    
    if isnothing(valid_idx)
        @user_warn "Could not achieve target FTR $(target_ftr), using most restrictive threshold"
        return maximum(sorted_scores)
    end
    
    threshold = sorted_scores[valid_idx]
    achieved_ftr = ftr[valid_idx]
    
    @user_info "MBR ML threshold calibration: threshold=$(round(threshold, digits=4)), achieved FTR=$(round(achieved_ftr, digits=4)), target FTR=$(target_ftr)"
    return threshold
end

function train_probit_mbr_model(X::Matrix, y::AbstractVector{Bool}, feature_names::Vector{Symbol})
    # Use the existing probit regression implementation
    # X already includes intercept column from prepare_mbr_features
    
    # Create DataFrame for logging (feature names help with debugging)
    X_df = DataFrame(X, feature_names)
    
    @user_info "Training probit model with $(size(X, 1)) samples, $(size(X, 2)) features"
    @user_info "Feature names: $feature_names"
    @user_info "Class balance: $(sum(y))/$(length(y)) positive ($(round(mean(y)*100, digits=1))%)"
    
    # Initialize coefficients
    β = zeros(Float64, size(X, 2))
    
    # Create data chunks for parallel processing (following FirstPassSearch pattern)
    tasks_per_thread = 10
    M = length(y)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = Iterators.partition(1:M, chunk_size)
    
    # Fit probit model using existing implementation (max_iter=30 as in FirstPassSearch)
    β_fitted = ProbitRegression(β, X_df, y, data_chunks, max_iter=30)
    
    @user_info "Probit coefficients fitted: $(round.(β_fitted, digits=4))"
    
    return β_fitted
end

function predict_probit_mbr(β::Vector{Float64}, X::Matrix, feature_names::Vector{Symbol})
    # Create DataFrame from feature matrix (following FirstPassSearch pattern)
    X_df = DataFrame(X, feature_names)
    
    # Create scores vector
    scores = zeros(Float32, size(X, 1))
    
    # Create data chunks for parallel processing (following FirstPassSearch pattern)
    tasks_per_thread = 10
    M = size(X, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = Iterators.partition(1:M, chunk_size)
    
    # Use ModelPredict! to get Z-scores (following FirstPassSearch pattern)
    ModelPredict!(scores, X_df, β, data_chunks)
    
    # Convert Z-scores to probabilities using erf (following FirstPassSearch pattern)
    probabilities = zeros(Float32, length(scores))
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                probabilities[i] = (1 + Pioneer.SpecialFunctions.erf(scores[i]/sqrt(2)))/2
            end
        end
    end
    fetch.(tasks)
    
    return probabilities
end

"""
    add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)

Add best trace indicator based on isotope trace type.
"""
function add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)
    op = function(df)  # df is passed by transform_and_write!
        if seperateTraces(isotope_type)
            # Extract columns with type assertions for performance
            precursor_idx_col = df.precursor_idx::AbstractVector{UInt32}
            isotopes_captured_col = df.isotopes_captured::AbstractVector{Tuple{Int8, Int8}}   
            
            # Efficient vectorized operation for separate traces
            df[!,:best_trace] = [
                (precursor_idx=precursor_idx_col[i], 
                 isotopes_captured=isotopes_captured_col[i]) ∈ best_traces
                for i in eachindex(precursor_idx_col)
            ]
        else
            # Group-based operation for combined traces
            transform!(groupby(df, :precursor_idx),
                      :prob => (p -> begin
                          best_idx = argmax(p)
                          result = falses(length(p))
                          result[best_idx] = true
                          result
                      end) => :best_trace)
        end
        return df
    end
    return "" => op
end

"""
    get_quant_necessary_columns() -> Vector{Symbol}

Get the standard columns needed for quantification analysis.
"""
function get_quant_necessary_columns()
    return [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :global_qval,
        :run_specific_qval,
        :prec_mz,
        :pep,
        :weight,
        :target,
        :rt,
        :irt_obs,
        :missed_cleavage,
        :Mox,
        :isotopes_captured,
        :scan_idx,
        :entrapment_group_id,
        :ms_file_idx
    ]
end
