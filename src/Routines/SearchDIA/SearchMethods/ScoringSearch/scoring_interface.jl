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
struct LightGBMFilter <: MBRFilterMethod end

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
    apply_mbr_filter!(merged_df, params, fdr_scale_factor)

Wrapper function for ScoringSearch compatibility.
Generates candidate mask and bad transfer labels automatically, then applies MBR filtering.
Returns column name for filtered probabilities.
"""
function apply_mbr_filter!(
    merged_df::DataFrame,
    params
)
    # 1) identify transfer candidates
    candidate_mask = merged_df.MBR_transfer_candidate
    n_candidates = sum(candidate_mask)

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| # T->D
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false)) # D->T
    )

    # 3) Apply the main filtering function
    filtered_probs = apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

    # 4) Add filtered probabilities as a new column and return column name
    merged_df[!, :MBR_filtered_prob] = filtered_probs
    return :MBR_filtered_prob
end

"""
    apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

Apply MBR filtering using automatic method selection.
Tests threshold, probit, and LightGBM methods, selects the one that passes the most candidates.
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
    
    # Handle case with no MBR candidates
    if n_candidates == 0
        @user_warn "No MBR transfer candidates found - returning original probabilities unchanged"
        return merged_df.prob
    end

    # Test all methods and store results
    methods = [ThresholdFilter(), ProbitFilter(), LightGBMFilter()]
    results = FilterResult[]
    
    for method in methods
        result = train_and_evaluate(method, candidate_data, candidate_labels, params)
        if result !== nothing
            push!(results, result)
        end
    end
    
    # Select best method (most candidates passing)
    if isempty(results)
        @user_warn "No MBR filtering methods succeeded - returning original probabilities unchanged"
        return merged_df.prob
    end
    
    best_result = results[argmax([r.n_passing for r in results])]
    
    @user_info "MBR Method Selection:"
    for result in results
        marker = result === best_result ? " ✓" : ""
        @user_info "  $(result.method_name): $(result.n_passing)/$n_candidates pass ($(round(100*result.n_passing/n_candidates, digits=1))%)$marker"
    end
    
    # Compute per-candidate MBR q-values (FTR) based on the selected method's scores
    try
        # Scores aligned to candidate_data row order
        scores = best_result.scores
        # Use is_bad_transfer as numerator flag; denominator counts all candidates at threshold
        candidate_qvals = Vector{Float32}(undef, length(scores))
        # get_ftr! expects placeholders for targets; it accumulates total candidates internally
        get_ftr!(scores, trues(length(scores)), candidate_labels, candidate_qvals)

        # Map back to full dataframe rows (only when candidates exist)
        full_qvals = Vector{Union{Missing, Float32}}(missing, nrow(merged_df))
        cand_indices = findall(candidate_mask)
        @inbounds for (i, idx) in enumerate(cand_indices)
            full_qvals[idx] = candidate_qvals[i]
        end

        # Add columns for downstream pipelines and outputs
        merged_df[!, :MBR_candidate] = merged_df.MBR_transfer_candidate
        merged_df[!, :MBR_transfer_q_value] = full_qvals
    catch e
        @user_warn "Failed to compute per-candidate MBR q-values: $(typeof(e)) — $(e)"
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
    # Handle empty candidate data
    if isempty(candidate_data) || !hasproperty(candidate_data, :prob)
        return nothing
    end

    # candidate_labels represents bad transfer flags
    τ = get_ftr_threshold(
        candidate_data.prob,
        candidate_labels,
        params.max_MBR_false_transfer_rate
    )

    # Handle edge case where threshold is infinite (no valid threshold found)
    if isinf(τ)
        n_passing = 0
    else
        n_passing = sum(candidate_data.prob .>= τ)
    end

    return FilterResult("Threshold", candidate_data.prob, τ, n_passing)
end

function train_and_evaluate(method::ProbitFilter, candidate_data::DataFrame, candidate_labels::AbstractVector{Bool}, params)
    try
        # Handle empty candidate data
        if isempty(candidate_data)
            return nothing
        end

        # Check CV fold availability
        if !hasproperty(candidate_data, :cv_fold)
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
        # Handle probit training errors gracefully in main search/tests
        @user_warn "Probit method failed: $(typeof(e)) — $(e)"
        return nothing
    end
end

function train_and_evaluate(method::LightGBMFilter, candidate_data::DataFrame, candidate_labels::AbstractVector{Bool}, params)
    try
        # Handle empty candidate data
        if isempty(candidate_data)
            return nothing
        end

        # Check CV fold availability
        if !hasproperty(candidate_data, :cv_fold)
            return nothing
        end

        # Feature preparation
        feature_cols = select_mbr_features(candidate_data)
        @debug "LightGBMFilter: Selected features: $feature_cols"
        @debug "LightGBMFilter: Candidate data size: $(size(candidate_data))"
        @debug "LightGBMFilter: Labels length: $(length(candidate_labels))"

        # Work directly with DataFrame like FirstPassSearch
        feature_data = candidate_data[:, feature_cols]
        @debug "LightGBMFilter: Feature data size: $(size(feature_data))"

        # Cross-validation training using DataFrame
        scores = run_cv_training(method, feature_data, candidate_labels, candidate_data.cv_fold, params)

        # Calibrate threshold
        τ = calibrate_ml_threshold(scores, candidate_labels, Float64(params.max_MBR_false_transfer_rate))
        n_passing = sum(scores .>= τ)  # Higher score = better for LightGBM


        return FilterResult("LightGBM", scores, τ, n_passing)
    catch e
        @user_warn "LightGBM method failed with error: $(typeof(e)): $e"
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

function run_cv_training(method::LightGBMFilter, feature_data::DataFrame, labels::AbstractVector{Bool}, cv_folds::AbstractVector, params)
    out_of_fold_scores = zeros(Float64, length(labels))

    for fold in unique(cv_folds)
        test_mask = cv_folds .== fold
        train_mask = .!test_mask

        # Train model on fold using DataFrame slices (like FirstPassSearch)
        train_data = feature_data[train_mask, :]
        train_labels = labels[train_mask]
        test_data = feature_data[test_mask, :]

        # Train and predict using DataFrames directly
        model = train_lightgbm_model_df(train_data, train_labels, params)
        test_predictions = lightgbm_predict(model, test_data)
        @debug "LightGBM predictions shape: $(size(test_predictions)), type: $(typeof(test_predictions))"
        out_of_fold_scores[test_mask] = test_predictions
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
            continue
        end

        # Check for zero variance (constant columns)
        if length(unique(col_data)) <= 1 || var(col_data) ≈ 0.0
            continue
        end

        push!(valid_cols, Symbol(col))
    end

    if isempty(valid_cols)
        throw(ArgumentError("No valid features for probit regression"))
    end

    # Use only valid columns
    filtered_data = feature_data[:, valid_cols]
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

function train_lightgbm_model_df(feature_data::DataFrame, y::AbstractVector{Bool}, params)
    labels = y .== false  # Invert labels so true indicates a good transfer
    classifier = build_lightgbm_classifier(
        num_iterations = 100,      # matches EvoTrees nrounds
        max_depth = 3,             # matches EvoTrees max_depth
        num_leaves = 8,            # 2^3 for max_depth=3
        learning_rate = 0.1,       # matches EvoTrees eta
        feature_fraction = 0.8,    # matches EvoTrees colsample
        bagging_fraction = 0.5,    # matches EvoTrees rowsample
        bagging_freq = 1,
        min_data_in_leaf = 1,      # matches EvoTrees min_child_weight
        min_gain_to_split = 1.0,   # matches EvoTrees gamma
    )
    return fit_lightgbm_model(classifier, feature_data, labels; positive_label=true)
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
    # Core features for MBR filtering (including MS1-MS2 RT difference feature)
    candidate_features = [:prob, :irt_error, :ms1_ms2_rt_diff, :MBR_max_pair_prob, :MBR_best_irt_diff,
                         :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio]#, :MBR_num_runs]
    
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

function calibrate_ml_threshold(scores::AbstractVector, is_bad_transfer::AbstractVector{Bool}, target_ftr::Float64)
    """Find score threshold that achieves target FTR."""
    return get_ftr_threshold(scores, is_bad_transfer, target_ftr)
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
        :MBR_candidate,
        :MBR_transfer_q_value,
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
