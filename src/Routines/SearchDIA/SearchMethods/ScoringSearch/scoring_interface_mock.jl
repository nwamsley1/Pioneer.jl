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
Clean refactored MBR filtering interface using trait-based design.
"""

#==========================================================
MBR Filter Method Traits
==========================================================#

abstract type MBRFilterMethod end

struct ThresholdFilter <: MBRFilterMethod end
struct ProbitFilter <: MBRFilterMethod end 
struct XGBoostFilter <: MBRFilterMethod end

struct FilterResult
    method_name::String
    scores::Vector{Float64}
    threshold::Float64
    n_passing::Int
end

#==========================================================
Main MBR Filtering Interface
==========================================================#

"""
    apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

Apply MBR filtering using automatic method selection.
Tests threshold, probit, and XGBoost methods, selects the one that passes the most candidates.
"""
function apply_mbr_filter!(
    merged_df::DataFrame,
    candidate_mask::AbstractVector{Bool},
    is_bad_transfer::AbstractVector{Bool},
    params
)
    # Extract candidate data once
    candidate_data = merged_df[candidate_mask, :]
    candidate_labels = is_bad_transfer[candidate_mask]
    
    n_candidates = length(candidate_labels)
    n_bad = sum(candidate_labels)
    @user_info "MBR Filtering: $n_candidates candidates, $n_bad bad transfers ($(round(100*n_bad/n_candidates, digits=1))%)"
    
    # Test all methods and store results
    methods = [ThresholdFilter(), ProbitFilter(), XGBoostFilter()]
    results = FilterResult[]
    
    for method in methods
        result = train_and_evaluate(method, candidate_data, candidate_labels, params)
        if result !== nothing
            push!(results, result)
        end
    end
    
    # Select best method (most candidates passing)
    if isempty(results)
        error("No MBR filtering methods succeeded")
    end
    
    best_result = results[argmax([r.n_passing for r in results])]
    
    @user_info "MBR Method Selection:"
    for result in results
        marker = result === best_result ? " ✓" : ""
        @user_info "  $(result.method_name): $(result.n_passing)/$n_candidates pass ($(round(100*result.n_passing/n_candidates, digits=1))%)$marker"
    end
    
    # Apply best method's filtering
    return apply_filtering(best_result, merged_df, candidate_mask, params)
end

#==========================================================
Method-Specific Training and Evaluation
==========================================================#

"""
    train_and_evaluate(method, candidate_data, candidate_labels, params) -> FilterResult

Train a filtering method and evaluate performance. Returns FilterResult with scores and threshold.
"""
function train_and_evaluate(method::ThresholdFilter, candidate_data::DataFrame, candidate_labels::AbstractVector{Bool}, params)
    try
        τ = get_ftr_threshold(
            candidate_data.prob,
            candidate_labels,
            params.max_MBR_false_transfer_rate
        )
        n_passing = sum(candidate_data.prob .>= τ)
        
        return FilterResult("Threshold", candidate_data.prob, τ, n_passing)
    catch e
        @user_warn "Threshold method failed: $e"
        return nothing
    end
end

function train_and_evaluate(method::ProbitFilter, candidate_data::DataFrame, candidate_labels::AbstractVector{Bool}, params)
    try
        # Check CV fold availability
        if !hasproperty(candidate_data, :cv_fold)
            @user_warn "No cv_fold column - skipping probit method"
            return nothing
        end
        
        # Feature preparation
        feature_cols = select_mbr_features(candidate_data)
        X, feature_names = prepare_mbr_features(candidate_data[:, feature_cols])
        
        # Cross-validation training
        scores = run_cv_training(method, X, candidate_labels, candidate_data.cv_fold, feature_names, params)
        
        # Calibrate threshold
        τ = calibrate_ml_threshold(scores, candidate_labels, Float64(params.max_MBR_false_transfer_rate))
        n_passing = sum(scores .< τ)  # Lower score = better for probit
        
        # Optional: Save to Arrow
        save_to_arrow(method, candidate_data, X, feature_names, candidate_labels, scores)
        
        return FilterResult("Probit", scores, τ, n_passing)
    catch e
        @user_warn "Probit method failed: $e"
        return nothing
    end
end

function train_and_evaluate(method::XGBoostFilter, candidate_data::DataFrame, candidate_labels::AbstractVector{Bool}, params)
    try
        # Check CV fold availability  
        if !hasproperty(candidate_data, :cv_fold)
            @user_warn "No cv_fold column - skipping XGBoost method"
            return nothing
        end
        
        # Feature preparation
        feature_cols = select_mbr_features(candidate_data)
        X, feature_names = prepare_mbr_features(candidate_data[:, feature_cols])
        
        # Cross-validation training
        scores = run_cv_training(method, X, candidate_labels, candidate_data.cv_fold, feature_names, params)
        
        # Calibrate threshold
        τ = calibrate_ml_threshold(scores, candidate_labels, Float64(params.max_MBR_false_transfer_rate))
        n_passing = sum(scores .< τ)  # Lower score = better for XGBoost
        
        # Optional: Save to Arrow
        save_to_arrow(method, candidate_data, X, feature_names, candidate_labels, scores)
        
        return FilterResult("XGBoost", scores, τ, n_passing)
    catch e
        @user_warn "XGBoost method failed: $e"
        return nothing
    end
end

#==========================================================
Cross-Validation Training
==========================================================#

"""
    run_cv_training(method, X, labels, cv_folds, feature_names, params) -> Vector{Float64}

Perform cross-validation training and return out-of-fold scores.
"""
function run_cv_training(method::ProbitFilter, X::Matrix, labels::AbstractVector{Bool}, cv_folds::AbstractVector, feature_names::Vector{Symbol}, params)
    out_of_fold_scores = zeros(Float64, length(labels))
    
    for fold in unique(cv_folds)
        test_mask = cv_folds .== fold
        train_mask = .!test_mask
        
        # Train model on fold
        model = train_probit_model(X[train_mask, :], labels[train_mask], feature_names)
        
        # Predict on test fold
        test_scores = predict_probit_model(model, X[test_mask, :], feature_names)
        out_of_fold_scores[test_mask] = test_scores
    end
    
    return out_of_fold_scores
end

function run_cv_training(method::XGBoostFilter, X::Matrix, labels::AbstractVector{Bool}, cv_folds::AbstractVector, feature_names::Vector{Symbol}, params)
    out_of_fold_scores = zeros(Float64, length(labels))
    
    for fold in unique(cv_folds)
        test_mask = cv_folds .== fold
        train_mask = .!test_mask
        
        # Train model on fold
        model = train_xgboost_model(X[train_mask, :], labels[train_mask], feature_names, params)
        
        # Predict on test fold
        test_predictions = EvoTrees.predict(model, X[test_mask, :])
        out_of_fold_scores[test_mask] = test_predictions[:, 2]  # Probability of positive class
    end
    
    return out_of_fold_scores
end

#==========================================================
Model Training Functions
==========================================================#

function train_probit_model(X::Matrix, y::AbstractVector{Bool}, feature_names::Vector{Symbol})
    # Convert to expected format for probit regression
    y_probit = convert(Vector{Bool}, y .== false)  # Invert labels for probit
    
    # Train probit model
    β = ProbitRegression(X, y_probit)
    
    return β
end

function predict_probit_model(β::Vector{Float64}, X::Matrix, feature_names::Vector{Symbol})
    scores = zeros(Float64, size(X, 1))
    ModelPredict!(β, X, scores)
    return scores
end

function train_xgboost_model(X::Matrix, y::AbstractVector{Bool}, feature_names::Vector{Symbol}, params)
    # Create training DataFrame
    training_df = DataFrame(X, feature_names)
    training_df[!, :target] = y .== false  # Invert labels for XGBoost
    
    # Configure XGBoost
    config = EvoTreeClassifier(
        loss = :logloss,
        nrounds = 50,
        max_depth = 4,
        eta = 0.1,
        rowsample = 0.5,
        colsample = 0.8,
        gamma = 1.0
    )
    
    # Train model
    model = EvoTrees.fit_evotree(config, training_df; target_name=:target)
    
    return model
end

#==========================================================
Filtering Application
==========================================================#

"""
    apply_filtering(result, merged_df, candidate_mask, params) -> Vector{Float32}

Apply the filtering result to the full dataframe.
"""
function apply_filtering(result::FilterResult, merged_df::DataFrame, candidate_mask::AbstractVector{Bool}, params)
    filtered_probs = copy(merged_df.prob)
    candidate_indices = findall(candidate_mask)
    
    if result.method_name == "Threshold"
        # Simple threshold on probability
        for idx in candidate_indices
            if merged_df.prob[idx] < result.threshold
                filtered_probs[idx] = 0.0f0
            end
        end
    else
        # ML-based filtering using scores
        for (i, idx) in enumerate(candidate_indices)
            if result.scores[i] >= result.threshold  # Higher score = worse candidate
                filtered_probs[idx] = 0.0f0
            end
        end
    end
    
    return filtered_probs
end

#==========================================================
Feature Processing (Simplified)
==========================================================#

function select_mbr_features(df::DataFrame)
    # Core features for MBR filtering
    candidate_features = [:prob, :irt_error, :MBR_max_pair_prob, :MBR_best_irt_diff, 
                         :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio, :MBR_num_runs]
    
    # Filter to available columns
    available_features = Symbol[]
    for feature in candidate_features
        if hasproperty(df, feature) && !all(ismissing, df[!, feature])
            push!(available_features, feature)
        end
    end
    
    return available_features
end

function prepare_mbr_features(df::DataFrame)
    # Simple preprocessing: handle missing values and add intercept
    processed_df = copy(df)
    
    # Replace missing with median
    for col in names(processed_df)
        col_data = processed_df[!, col]
        if any(ismissing, col_data)
            non_missing = collect(skipmissing(col_data))
            if !isempty(non_missing)
                median_val = median(non_missing)
                processed_df[!, col] = coalesce.(col_data, median_val)
            end
        end
    end
    
    # Build feature matrix with intercept
    n_rows = nrow(processed_df)
    feature_names = [:intercept; Symbol.(names(processed_df))]
    
    # Create matrix
    X = hcat(ones(Float64, n_rows), Matrix{Float64}(processed_df))
    
    return X, feature_names
end

function calibrate_ml_threshold(scores::AbstractVector, y_true::AbstractVector{Bool}, target_ftr::Float64)
    """Find score threshold that achieves target FTR."""
    return get_ftr_threshold(scores, y_true, target_ftr)
end

#==========================================================
Optional Arrow Export
==========================================================#

function save_to_arrow(method::ProbitFilter, candidate_data::DataFrame, X::Matrix, feature_names::Vector{Symbol}, labels::AbstractVector{Bool}, scores::Vector{Float64})
    try
        export_df = copy(candidate_data)
        
        # Add processed features
        for (i, name) in enumerate(feature_names)
            export_df[!, Symbol("processed_", name)] = X[:, i]
        end
        
        # Add results
        export_df[!, :is_bad_transfer] = labels
        export_df[!, :probit_score] = scores
        
        # Save
        arrow_path = "/Users/nathanwamsley/Desktop/mbr_probit_features_scores.arrow"
        Arrow.write(arrow_path, export_df)
        @user_info "Saved probit data: $(nrow(export_df)) rows, $(ncol(export_df)) columns"
    catch e
        @user_warn "Arrow export failed: $e"
    end
end

function save_to_arrow(method::XGBoostFilter, candidate_data::DataFrame, X::Matrix, feature_names::Vector{Symbol}, labels::AbstractVector{Bool}, scores::Vector{Float64})
    try
        export_df = copy(candidate_data)
        
        # Add processed features
        for (i, name) in enumerate(feature_names)
            export_df[!, Symbol("processed_", name)] = X[:, i]
        end
        
        # Add results
        export_df[!, :is_bad_transfer] = labels
        export_df[!, :xgboost_score] = scores
        
        # Save
        arrow_path = "/Users/nathanwamsley/Desktop/mbr_xgboost_features_scores.arrow"
        Arrow.write(arrow_path, export_df)
        @user_info "Saved XGBoost data: $(nrow(export_df)) rows, $(ncol(export_df)) columns"
    catch e
        @user_warn "Arrow export failed: $e"
    end
end

save_to_arrow(method::ThresholdFilter, args...) = nothing  # No export for threshold method

#==========================================================
Additional Interface Functions (Preserved from Original)
==========================================================#

"""
    calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    
Calculate global protein scores and add them to files via references.
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
            :global_pg_score, :target;
            reverse=true
        )
    end
    
    return acc_to_global_score
end

function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    if length(probs) == 1
        return first(probs)
    end
    
    # Use log-odds transformation for score combination
    # This handles probabilities close to 0 and 1 more stably
    log_odds_sum = zero(T)
    count = 0
    
    # Take top N scores to avoid being dominated by many weak scores
    sorted_probs = sort(probs, rev=true)
    n_scores = min(top_n, length(sorted_probs))
    
    for i in 1:n_scores
        p = clamp(sorted_probs[i], T(1e-7), T(1-1e-7))  # Avoid log(0)
        log_odds_sum += log(p / (1 - p))
        count += 1
    end
    
    if count == 0
        return T(0.5)  # Default neutral score
    end
    
    # Convert back to probability
    avg_log_odds = log_odds_sum / count
    return T(1 / (1 + exp(-avg_log_odds)))
end

# Add other essential functions as needed...
function add_trace_qvalues(fdr_scale_factor::Float32)
    op = function(df)
        qvals = Vector{Float32}(undef, nrow(df))
        get_qvalues!(df.prob, df.target, qvals; fdr_scale_factor=fdr_scale_factor)
        df[!, :trace_qval] = qvals
        return df
    end
    return "add_trace_qvalues" => op
end

function add_prec_prob(prob_col::Symbol)
    op = function(df)
        gdf = groupby(df, [:precursor_idx, :filename])
        combined = combine(gdf) do chunk
            chunk[!, :prec_prob] .= 1.0f0 .- reduce(*, 1.0f0 .- chunk[!, prob_col])
            return chunk
        end
        return combined
    end
    return "add_prec_prob" => op  
end