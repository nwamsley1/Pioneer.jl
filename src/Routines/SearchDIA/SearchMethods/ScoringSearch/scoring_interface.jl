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
        @debug "ProbitFilter: Selected features: $feature_cols"
        @debug "ProbitFilter: Candidate data size: $(size(candidate_data))"
        @debug "ProbitFilter: Labels length: $(length(candidate_labels))"
        
        # Work directly with DataFrame like FirstPassSearch
        feature_data = candidate_data[:, feature_cols]
        @debug "ProbitFilter: Feature data size: $(size(feature_data))"
        
        # Cross-validation training using DataFrame
        scores = run_cv_training(method, feature_data, candidate_labels, candidate_data.cv_fold, params)
        
        # Calibrate threshold
        τ = calibrate_ml_threshold(scores, candidate_labels, Float64(params.max_MBR_false_transfer_rate))
        n_passing = sum(scores .>= τ)  # Higher score = better for probit
        
        
        return FilterResult("Probit", scores, τ, n_passing)
    catch e
        @user_warn "Probit method failed with error: $(typeof(e)): $e"
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
        @debug "XGBoostFilter: Selected features: $feature_cols"
        @debug "XGBoostFilter: Candidate data size: $(size(candidate_data))"
        @debug "XGBoostFilter: Labels length: $(length(candidate_labels))"

        # Work directly with DataFrame like FirstPassSearch
        feature_data = candidate_data[:, feature_cols]
        @debug "XGBoostFilter: Feature data size: $(size(feature_data))"

        # Cross-validation training using DataFrame
        scores = run_cv_training(method, feature_data, candidate_labels, candidate_data.cv_fold, params)
        
        # Calibrate threshold
        τ = calibrate_ml_threshold(scores, candidate_labels, Float64(params.max_MBR_false_transfer_rate))
        n_passing = sum(scores .>= τ)  # Higher score = better for XGBoost
        
        
        return FilterResult("XGBoost", scores, τ, n_passing)
    catch e
        @user_warn "XGBoost method failed with error: $(typeof(e)): $e"
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
function run_cv_training(method::ProbitFilter, feature_data::DataFrame, labels::AbstractVector{Bool}, cv_folds::AbstractVector, params)
    out_of_fold_scores = zeros(Float64, length(labels))
    
    for fold in unique(cv_folds)
        test_mask = cv_folds .== fold
        train_mask = .!test_mask
        
        # Train model on fold using DataFrame slices (like FirstPassSearch)
        train_data = feature_data[train_mask, :]
        train_labels = labels[train_mask]
        test_data = feature_data[test_mask, :]
        
        # Train and predict using DataFrames directly
        model, valid_cols = train_probit_model_df(train_data, train_labels, params)
        test_scores = predict_probit_model_df(model, test_data, valid_cols)
        out_of_fold_scores[test_mask] = test_scores
    end
    
    return out_of_fold_scores
end

function run_cv_training(method::XGBoostFilter, feature_data::DataFrame, labels::AbstractVector{Bool}, cv_folds::AbstractVector, params)
    out_of_fold_scores = zeros(Float64, length(labels))

    for fold in unique(cv_folds)
        test_mask = cv_folds .== fold
        train_mask = .!test_mask

        # Train model on fold using DataFrame slices (like FirstPassSearch)
        train_data = feature_data[train_mask, :]
        train_labels = labels[train_mask]
        test_data = feature_data[test_mask, :]

        # Train and predict using DataFrames directly
        model = train_xgboost_model_df(train_data, train_labels, params)
        test_predictions = EvoTrees.predict(model, test_data)
        @debug "XGBoost predictions shape: $(size(test_predictions)), type: $(typeof(test_predictions))"

        # Handle different prediction formats
        if isa(test_predictions, Vector)
            out_of_fold_scores[test_mask] = test_predictions
        elseif size(test_predictions, 2) >= 2
            out_of_fold_scores[test_mask] = test_predictions[:, 2]  # Probability of positive class
        else
            @warn "Unexpected XGBoost prediction format: $(size(test_predictions))"
            out_of_fold_scores[test_mask] = vec(test_predictions)
        end
    end

    return out_of_fold_scores
end

#==========================================================
Model Training Functions
==========================================================#

function train_probit_model_df(feature_data::DataFrame, y::AbstractVector{Bool}, params)
    # Convert to expected format for probit regression
    y_probit = convert(Vector{Bool}, y .== false)  # Invert labels for probit

    # Check for problematic columns (zero variance or containing Inf/NaN)
    valid_cols = Symbol[]
    for col in names(feature_data)
        col_data = feature_data[!, col]

        # Check for Inf or NaN values
        if any(isinf, col_data) || any(isnan, col_data)
            @user_warn "Probit MBR: Skipping column $col (contains Inf/NaN values)"
            continue
        end

        # Check for zero variance (constant columns)
        if length(unique(col_data)) <= 1 || var(col_data) ≈ 0.0
            @user_warn "Probit MBR: Skipping column $col (zero variance)"
            continue
        end

        push!(valid_cols, Symbol(col))
    end

    if isempty(valid_cols)
        @user_warn "Probit MBR: No valid columns for training - all have zero variance or Inf/NaN"
        throw(ArgumentError("No valid features for probit regression"))
    end

    # Use only valid columns
    filtered_data = feature_data[:, valid_cols]
    @user_info "Probit MBR: Using $(length(valid_cols))/$(ncol(feature_data)) valid features: $(join(valid_cols, ", "))"

    # Initialize coefficients for filtered data
    β = zeros(Float64, size(filtered_data, 2))

    # Create data chunks for parallel processing
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, length(y_probit) / n_chunks))
    data_chunks = Iterators.partition(1:length(y_probit), chunk_size)

    # Train probit model directly with DataFrame (like FirstPassSearch)
    β_fitted = ProbitRegression(β, filtered_data, y_probit, data_chunks, max_iter=30)

    return β_fitted, valid_cols  # Return both model and valid column names
end

function predict_probit_model_df(β::Vector{Float64}, feature_data::DataFrame, valid_cols::Vector{Symbol})
    scores = zeros(Float64, size(feature_data, 1))

    # Use only the valid columns that were used for training
    filtered_data = feature_data[:, valid_cols]

    # Create data chunks for parallel processing
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, size(filtered_data, 1) / n_chunks))
    data_chunks = Iterators.partition(1:size(filtered_data, 1), chunk_size)

    # Call ModelPredict! directly with DataFrame (like FirstPassSearch)
    ModelPredict!(scores, filtered_data, β, data_chunks)
    return scores
end

function train_xgboost_model_df(feature_data::DataFrame, y::AbstractVector{Bool}, params)
    # Create training DataFrame with target column
    training_df = copy(feature_data)
    training_df[!, :target] = y .== false  # Invert labels for XGBoost

    # Configure XGBoost
    config = EvoTreeClassifier(
        nrounds = 100,
        max_depth = 3,
        eta = 0.1,
        rowsample = 0.5,
        colsample = 0.8,
        gamma = 1.0
    )

    # Train model directly with DataFrame (no Matrix conversion)
    model = EvoTrees.fit(config, training_df; target_name=:target)

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
            if result.scores[i] < result.threshold  # Lower score = worse candidate (bad transfer)
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
    candidate_features = [:prob, :irt_error, :rt_diff, :MBR_max_pair_prob, :MBR_best_irt_diff,
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

