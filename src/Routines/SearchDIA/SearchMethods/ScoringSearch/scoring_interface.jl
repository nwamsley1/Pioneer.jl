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
    # 1) identify transfer candidates based on MBR_boosted_trace_prob
    candidate_mask = merged_df.MBR_transfer_candidate
    n_candidates = sum(candidate_mask)

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| # T->D
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false)) # D->T
    )

    # 3) Apply filtering and get filtered probabilities
    filtered_probs = apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

    # 4) Modify both trace_prob and MBR_boosted_trace_prob columns IN-PLACE
    merged_df[!, :MBR_boosted_trace_prob] = filtered_probs.MBR_boosted_trace_prob
    merged_df[!, :trace_prob] = filtered_probs.trace_prob

    # 5) No return value needed (modifies columns in place)
    return nothing
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
        return (MBR_boosted_trace_prob = merged_df.MBR_boosted_trace_prob,
                trace_prob = merged_df.trace_prob)
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
        return (MBR_boosted_trace_prob = merged_df.MBR_boosted_trace_prob,
                trace_prob = merged_df.trace_prob)
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
    return apply_filtering(best_result, merged_df, candidate_mask, is_bad_transfer, params)
end

#==========================================================
OOM Streaming MBR Filter + Per-File Precursor Aggregation
==========================================================#

"""
    _aggregate_trace_to_precursor_probs!(df, has_mbr)

Per-file Bayesian aggregation of trace-level → precursor-level probabilities.
Groups by (precursor_idx, ms_file_idx). Since ms_file_idx is constant within
a single file, this is effectively grouping by precursor_idx alone.
"""
function _aggregate_trace_to_precursor_probs!(df::DataFrame, has_mbr::Bool)
    prob_agg = p -> begin
        trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
        clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
    end
    if has_mbr && hasproperty(df, :MBR_boosted_trace_prob)
        transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
                   :MBR_boosted_trace_prob => prob_agg => :MBR_boosted_prec_prob)
    end
    transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => prob_agg => :prec_prob)
end

"""
    _compute_ftr_threshold_streaming(merged_path, target_ftr) -> (threshold, n_passing)

Stream through a merged Arrow file (sorted by score desc) computing FTR.
Returns the minimum score where FTR ≤ target_ftr and the count of candidates
passing at that threshold. Cleans up the merged file after reading.

Same pattern as _compute_prob_threshold_from_merged (scoring_workspace.jl).
"""
function _compute_ftr_threshold_streaming(merged_path::String, target_ftr::Float64)
    tbl = Arrow.Table(merged_path)
    scores = tbl.score
    bad_flags = tbl.is_bad_transfer
    n = length(scores)

    n_total = 0
    n_bad = 0
    threshold = typemax(Float64)
    n_passing = 0

    @inbounds for i in 1:n
        n_total += 1
        if bad_flags[i]
            n_bad += 1
        end
        if n_total > 0 && n_bad / n_total <= target_ftr
            threshold = Float64(scores[i])
            n_passing = n_total
        end
    end

    # Release mmap before deleting (same pattern as scoring_workspace.jl)
    tbl = nothing
    safeRm(merged_path, nothing; force=true)
    return threshold, n_passing
end

"""
    apply_mbr_filter_and_aggregate_per_file!(refs, valid_file_indices, params)

OOM-safe MBR filtering and precursor probability aggregation.
Uses a 4-phase streaming approach:
  Phase 1: SAMPLE — reservoir-sample candidates from files for training
  Phase 2: TRAIN  — train per-fold Probit + LightGBM on bounded sample
  Phase 3: EVAL   — predict all candidates streaming, compute FTR thresholds, select best method
  Phase 4: APPLY  — re-predict + apply threshold + aggregate precursor probs per file

At no point is more than one file's data in memory (beyond the bounded training sample).
"""
function apply_mbr_filter_and_aggregate_per_file!(
    refs::Vector{PSMFileReference},
    valid_file_indices::Vector{Int64},
    params
)
    has_mbr = false
    # Quick check: does the first non-empty file have MBR columns?
    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        if length(tbl) > 0
            has_mbr = hasproperty(tbl, :MBR_boosted_trace_prob)
            break
        end
    end

    if !has_mbr
        # No MBR: just aggregate per-file and return
        for ref in refs
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            _aggregate_trace_to_precursor_probs!(df, false)
            write_arrow_file(ref, df)
        end
        return nothing
    end

    # Training columns — only what select_mbr_features + model training actually need
    training_columns = [
        :trace_prob, :irt_error, :ms1_ms2_rt_diff,
        :MBR_max_pair_prob, :MBR_best_irt_diff, :MBR_rv_coefficient,
        :MBR_log2_weight_ratio, :MBR_log2_explained_ratio,
        :MBR_boosted_trace_prob,
        :cv_fold,
        :target, :decoy, :MBR_is_best_decoy
    ]

    max_candidates = params.max_mbr_training_candidates

    #--------------------------------------------------------------
    # Phase 1: SAMPLE — stream files, reservoir-sample candidates
    #--------------------------------------------------------------
    sample_df = DataFrame()
    n_candidates_total = 0
    rng = MersenneTwister(1844)

    for ref in refs
        df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
        mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
        n_cand = sum(mask)
        n_cand == 0 && continue
        n_candidates_total += n_cand

        cols_available = filter(c -> hasproperty(df, c), training_columns)
        file_candidates = df[mask, cols_available]

        if nrow(sample_df) == 0
            sample_df = file_candidates
        else
            append!(sample_df, file_candidates)
        end
        # Downsample if over budget
        if nrow(sample_df) > max_candidates
            keep = randperm(rng, nrow(sample_df))[1:max_candidates]
            sample_df = sample_df[sort(keep), :]
        end
    end

    # Handle case with no MBR candidates at all
    if n_candidates_total == 0
        @user_warn "No MBR transfer candidates found - returning original probabilities unchanged"
        for ref in refs
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            df[!, :MBR_candidate] = falses(nrow(df))
            df[!, :MBR_transfer_q_value] = Vector{Union{Missing, Float32}}(missing, nrow(df))
            _aggregate_trace_to_precursor_probs!(df, true)
            write_arrow_file(ref, df)
        end
        return nothing
    end

    #--------------------------------------------------------------
    # Phase 2: TRAIN — per-fold models on bounded sample
    #--------------------------------------------------------------
    sample_labels = (
        (sample_df.target .& coalesce.(sample_df.MBR_is_best_decoy, false)) .|
        (sample_df.decoy .& .!coalesce.(sample_df.MBR_is_best_decoy, false))
    )

    feature_cols = select_mbr_features(sample_df)
    feature_data = sample_df[:, feature_cols]
    folds = sort(unique(sample_df.cv_fold))

    # Train Probit: test_fold => (model, valid_cols)
    probit_models = Dict{eltype(folds), Any}()
    probit_valid_cols = Dict{eltype(folds), Vector{Symbol}}()
    probit_ok = true
    try
        for fold in folds
            train_mask = sample_df.cv_fold .!= fold
            model, vcols = train_probit_model_df(feature_data[train_mask, :], sample_labels[train_mask], params)
            probit_models[fold] = model
            probit_valid_cols[fold] = vcols
        end
    catch e
        @user_warn "Probit training failed: $(typeof(e)) — $(e)"
        probit_ok = false
    end

    # Train LightGBM: test_fold => booster
    lgbm_models = Dict{eltype(folds), Any}()
    lgbm_ok = true
    try
        for fold in folds
            train_mask = sample_df.cv_fold .!= fold
            lgbm_models[fold] = train_lightgbm_model_df(feature_data[train_mask, :], sample_labels[train_mask], params)
        end
    catch e
        @user_warn "LightGBM training failed: $(typeof(e)) — $(e)"
        lgbm_ok = false
    end

    sample_df = nothing; feature_data = nothing  # free training data

    #--------------------------------------------------------------
    # Phase 3: EVAL — predict all candidates, compute FTR thresholds
    #--------------------------------------------------------------
    method_names = String["Threshold"]
    probit_ok && push!(method_names, "Probit")
    lgbm_ok && push!(method_names, "LightGBM")
    method_temp_refs = Dict(m => PSMFileReference[] for m in method_names)

    # Single pass through all files — predict candidates with all methods
    for ref in refs
        df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
        mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
        n_cand = sum(mask)
        n_cand == 0 && continue

        cand_indices = findall(mask)
        cand_features = df[cand_indices, feature_cols]

        is_bad = (
            (df.target[cand_indices] .& coalesce.(df.MBR_is_best_decoy[cand_indices], false)) .|
            (df.decoy[cand_indices] .& .!coalesce.(df.MBR_is_best_decoy[cand_indices], false))
        )

        for method_name in method_names
            scores = if method_name == "Threshold"
                Float64.(df.MBR_boosted_trace_prob[cand_indices])
            elseif method_name == "Probit"
                out = zeros(Float64, n_cand)
                for fold in folds
                    fm = df.cv_fold[cand_indices] .== fold
                    any(fm) || continue
                    out[fm] = predict_probit_model_df(probit_models[fold], cand_features[fm, :], probit_valid_cols[fold])
                end
                out
            else  # LightGBM
                out = zeros(Float64, n_cand)
                for fold in folds
                    fm = df.cv_fold[cand_indices] .== fold
                    any(fm) || continue
                    out[fm] = lightgbm_predict(lgbm_models[fold], cand_features[fm, :])
                end
                out
            end

            # Sort by score desc, write temp Arrow
            perm = sortperm(scores; rev=true)
            temp_path = tempname() * "_mbr_ftr_$(method_name).arrow"
            writeArrow(temp_path, DataFrame(score=scores[perm], is_bad_transfer=is_bad[perm]))
            tref = PSMFileReference(temp_path)
            mark_sorted!(tref, :score)
            push!(method_temp_refs[method_name], tref)
        end
    end

    # For each method: merge → stream FTR → threshold + n_passing
    method_results = Dict{String, Tuple{Float64, Int}}()
    for method_name in method_names
        trefs = method_temp_refs[method_name]
        if isempty(trefs)
            method_results[method_name] = (typemax(Float64), 0)
            continue
        end
        merged_path = tempname() * "_mbr_ftr_merged_$(method_name).arrow"
        stream_sorted_merge(trefs, merged_path, :score; reverse=true)
        GC.gc(false)
        for tref in trefs
            safeRm(file_path(tref), nothing; force=true)
        end
        threshold, n_passing = _compute_ftr_threshold_streaming(
            merged_path, Float64(params.max_MBR_false_transfer_rate))
        method_results[method_name] = (threshold, n_passing)
    end

    # Select best method (most passing candidates)
    best_method = method_names[argmax([method_results[m][2] for m in method_names])]
    best_threshold = method_results[best_method][1]

    @user_info "MBR Method Selection:"
    for method_name in method_names
        _, n_pass = method_results[method_name]
        marker = method_name == best_method ? " ✓" : ""
        @user_info "  $(method_name): $(n_pass)/$(n_candidates_total) pass ($(round(100*n_pass/n_candidates_total, digits=1))%)$marker"
    end

    #--------------------------------------------------------------
    # Phase 4: APPLY — re-predict + threshold + aggregate per file
    #--------------------------------------------------------------
    for ref in refs
        df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
        cand_mask = hasproperty(df, :MBR_transfer_candidate) ? df.MBR_transfer_candidate : falses(nrow(df))
        n_cand = sum(cand_mask)

        if n_cand > 0
            cand_indices = findall(cand_mask)

            # Re-predict with best method
            if best_method == "Threshold"
                for row in cand_indices
                    if Float64(df.MBR_boosted_trace_prob[row]) < best_threshold
                        df.MBR_boosted_trace_prob[row] = 0.0f0
                        df.trace_prob[row] = 0.0f0
                    end
                end
            else
                cand_features = df[cand_indices, feature_cols]
                scores = zeros(Float64, n_cand)
                for fold in folds
                    fm = df.cv_fold[cand_indices] .== fold
                    any(fm) || continue
                    if best_method == "Probit"
                        scores[fm] = predict_probit_model_df(
                            probit_models[fold], cand_features[fm, :], probit_valid_cols[fold])
                    else
                        scores[fm] = lightgbm_predict(lgbm_models[fold], cand_features[fm, :])
                    end
                end
                for (i, row) in enumerate(cand_indices)
                    if scores[i] < best_threshold
                        df.MBR_boosted_trace_prob[row] = 0.0f0
                        df.trace_prob[row] = 0.0f0
                    end
                end
            end

            # Zero out bad transfers regardless of score
            for row in cand_indices
                is_bad = (df.target[row] && coalesce(df.MBR_is_best_decoy[row], false)) ||
                         (df.decoy[row] && !coalesce(df.MBR_is_best_decoy[row], false))
                if is_bad
                    df.MBR_boosted_trace_prob[row] = 0.0f0
                    df.trace_prob[row] = 0.0f0
                end
            end

            # Diagnostic columns
            df[!, :MBR_candidate] = df.MBR_transfer_candidate
            df[!, :MBR_transfer_q_value] = Vector{Union{Missing, Float32}}(missing, nrow(df))
        else
            df[!, :MBR_candidate] = falses(nrow(df))
            df[!, :MBR_transfer_q_value] = Vector{Union{Missing, Float32}}(missing, nrow(df))
        end

        # Precursor aggregation (per-file, no merge needed)
        _aggregate_trace_to_precursor_probs!(df, true)
        write_arrow_file(ref, df)
    end

    return nothing
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
    if isempty(candidate_data)
        return nothing
    end

    # Use MBR-boosted scores if available, otherwise base scores
    score_column = hasproperty(candidate_data, :MBR_boosted_trace_prob) ?
                   :MBR_boosted_trace_prob : :trace_prob

    if !hasproperty(candidate_data, score_column)
        return nothing
    end

    # candidate_labels represents bad transfer flags
    τ = get_ftr_threshold(
        candidate_data[!, score_column],
        candidate_labels,
        params.max_MBR_false_transfer_rate
    )

    # Handle edge case where threshold is infinite (no valid threshold found)
    if isinf(τ)
        n_passing = 0
    else
        n_passing = sum(candidate_data[!, score_column] .>= τ)
    end

    return FilterResult("Threshold", candidate_data[!, score_column], τ, n_passing)
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
        is_unbalance = true
    )
    return fit_lightgbm_model(classifier, feature_data, labels; positive_label=true)
end

#==========================================================
Filtering Application
==========================================================#

"""
    apply_filtering(result, merged_df, candidate_mask, is_bad_transfer, params) -> NamedTuple

Apply the filtering result to the full dataframe, zeroing out candidates below the
selected threshold **and** any bad transfers regardless of their score.
Returns a NamedTuple with both MBR_boosted_trace_prob and trace_prob filtered.
"""
function apply_filtering(result::FilterResult, 
    merged_df::DataFrame, 
    candidate_mask::AbstractVector{Bool}, 
    is_bad_transfer::AbstractVector{Bool},
    params)
    
    filtered_MBR_boosted_trace_probs = copy(merged_df.MBR_boosted_trace_prob)
    filtered_trace_probs = copy(merged_df.trace_prob)
    candidate_indices = findall(candidate_mask)

    if result.method_name == "Threshold"
        # Simple threshold on probability
        for idx in candidate_indices
            if merged_df.MBR_boosted_trace_prob[idx] < result.threshold
                filtered_MBR_boosted_trace_probs[idx] = 0.0f0
                filtered_trace_probs[idx] = 0.0f0
            end
        end
    else
        # ML-based filtering using scores
        for (i, idx) in enumerate(candidate_indices)
            if result.scores[i] < result.threshold  # Lower score = worse candidate (bad transfer)
                filtered_MBR_boosted_trace_probs[idx] = 0.0f0
                filtered_trace_probs[idx] = 0.0f0
            end
        end
    end

    # Explicitly remove bad transfers even if they cleared the threshold
    for idx in candidate_indices
        if is_bad_transfer[idx]
            filtered_MBR_boosted_trace_probs[idx] = 0.0f0
            filtered_trace_probs[idx] = 0.0f0
        end
    end

    return (MBR_boosted_trace_prob = filtered_MBR_boosted_trace_probs,
            trace_prob = filtered_trace_probs)
end

#==========================================================
Feature Processing (Simplified)
==========================================================#

function select_mbr_features(df::DataFrame)
    # Core features for MBR filtering (including MS1-MS2 RT difference feature)
    candidate_features = [
                        :trace_prob,
                        :irt_error, :ms1_ms2_rt_diff, :MBR_max_pair_prob, :MBR_best_irt_diff,
                        :MBR_rv_coefficient, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio
                        #, :MBR_num_runs
                        ]
    
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
function get_quant_necessary_columns(match_between_runs::Bool)
    # Get conditional q-value column names
    qval_cols = get_qval_column_names(match_between_runs)

    base_cols = [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        qval_cols.global_qval,  # Conditional: :MBR_boosted_global_qval or :global_qval
        :run_specific_qval,
        :prec_mz,
        qval_cols.pep,  # Always :pep but included for consistency
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

    if match_between_runs
        # Add MBR-specific columns including MBR_boosted_qval
        return vcat(base_cols, [
            :MBR_boosted_global_prob,
            :MBR_boosted_prec_prob,
            :MBR_boosted_trace_prob,
            :MBR_candidate,
            :MBR_transfer_q_value,
            qval_cols.qval  # Add :MBR_boosted_qval
        ])
    else
        # Add standard qval for non-MBR mode
        return vcat(base_cols, [qval_cols.qval])  # Add :qval
    end
end

"""
    get_qval_column_names(match_between_runs::Bool) -> NamedTuple

Returns the appropriate q-value column names based on whether MBR is enabled.

# Arguments
- `match_between_runs::Bool`: Whether match-between-runs (MBR) is enabled

# Returns
A NamedTuple with fields:
- `qval::Symbol`: Precursor q-value column name
- `global_qval::Symbol`: Global precursor q-value column name
- `pep::Symbol`: PEP column name (always `:pep`)

# Examples
```julia
qval_cols = get_qval_column_names(true)  # MBR enabled
# Returns: (qval = :MBR_boosted_qval, global_qval = :MBR_boosted_global_qval, pep = :pep)

qval_cols = get_qval_column_names(false)  # MBR disabled
# Returns: (qval = :qval, global_qval = :global_qval, pep = :pep)

# Usage:
df_filtered = filter(row -> row[qval_cols.qval] < 0.01, df)
```
"""
function get_qval_column_names(match_between_runs::Bool)
    if match_between_runs
        return (
            qval = :MBR_boosted_qval,
            global_qval = :MBR_boosted_global_qval,
            pep = :pep
        )
    else
        return (
            qval = :qval,
            global_qval = :global_qval,
            pep = :pep
        )
    end
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
                      :trace_prob => (p -> begin
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

#==========================================================
Dictionary + Sidecar Helper Functions for OOM Scoring Pipeline
==========================================================#

"""
    build_precursor_global_prob_dicts(refs, sqrt_n_runs, has_mbr, n_precursors)
    → (global_prob_dict, mbr_global_prob_dict, target_dict)

Stream per-file to build precursor_idx → global_prob dictionaries.
Reads only (precursor_idx, prec_prob, [MBR_boosted_prec_prob], target) via mmap
instead of loading all 20+ columns.

`n_precursors` is used for `sizehint!()` pre-allocation of all dictionaries,
obtained via `length(getPrecursors(getSpecLib(search_context)))`.
"""
function build_precursor_global_prob_dicts(
    refs::Vector{PSMFileReference},
    sqrt_n_runs::Int,
    has_mbr::Bool,
    n_precursors::Int
)
    # Pre-allocate accumulation dictionaries with known upper bound
    prob_acc = Dict{UInt32, Vector{Float32}}()
    sizehint!(prob_acc, n_precursors)
    target_dict = Dict{UInt32, Bool}()
    sizehint!(target_dict, n_precursors)
    mbr_acc = if has_mbr
        d = Dict{UInt32, Vector{Float32}}()
        sizehint!(d, n_precursors)
        d
    else
        nothing
    end

    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        n = length(tbl.precursor_idx)
        n == 0 && continue
        prec_ids = tbl.precursor_idx
        prec_probs = tbl.prec_prob
        targets = tbl.target
        mbr_probs = has_mbr && hasproperty(tbl, :MBR_boosted_prec_prob) ? tbl.MBR_boosted_prec_prob : nothing

        @inbounds for i in 1:n
            pid = prec_ids[i]
            if !haskey(prob_acc, pid)
                prob_acc[pid] = Float32[]
                target_dict[pid] = targets[i]
                has_mbr && mbr_probs !== nothing && (mbr_acc[pid] = Float32[])
            end
            push!(prob_acc[pid], prec_probs[i])
            has_mbr && mbr_probs !== nothing && push!(mbr_acc[pid], mbr_probs[i])
        end
    end

    # Compute logodds per precursor
    global_prob_dict = Dict{UInt32, Float32}()
    sizehint!(global_prob_dict, length(prob_acc))
    for (pid, probs) in prob_acc
        global_prob_dict[pid] = logodds(probs, sqrt_n_runs)
    end

    mbr_global_prob_dict = Dict{UInt32, Float32}()
    if has_mbr && mbr_acc !== nothing
        sizehint!(mbr_global_prob_dict, length(mbr_acc))
        for (pid, probs) in mbr_acc
            mbr_global_prob_dict[pid] = logodds(probs, sqrt_n_runs)
        end
    end

    return global_prob_dict, mbr_global_prob_dict, target_dict
end

"""
    build_global_qval_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global q-values from a score dictionary without any file I/O.
Replaces get_precursor_global_qval_dict (which loaded a full merged DataFrame).
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
    write_score_sidecars(refs, columns; temp_prefix) → Vector{PSMFileReference}

Extract only the named columns from each file into a temporary Arrow sidecar file.
Uses Arrow mmap to avoid loading all columns. Caller must cleanup returned refs.
"""
function write_score_sidecars(
    refs::Vector{<:FileReference},
    columns::Vector{Symbol};
    temp_prefix::String = "sidecar"
)
    sidecar_refs = PSMFileReference[]
    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        n = length(Tables.getcolumn(tbl, first(columns)))
        n == 0 && continue

        # Extract only needed columns (collect to materialize from mmap)
        col_data = NamedTuple{Tuple(columns)}(Tuple(collect(Tables.getcolumn(tbl, c)) for c in columns))
        temp_path = tempname() * "_$(temp_prefix).arrow"
        writeArrow(temp_path, DataFrame(; col_data...))
        push!(sidecar_refs, PSMFileReference(temp_path))
    end
    return sidecar_refs
end

"""
    build_qvalue_spline_from_refs(refs, score_col, merged_path; ...) → Union{Nothing, NamedTuple}

Encapsulates the full sidecar lifecycle: write → sort → merge → cleanup → spline computation.
Returns `nothing` if all files are empty, otherwise returns `(; qval_spline, pep_interp)`.
`pep_interp` is `nothing` when `compute_pep=false`.

Uses `try/finally` to guarantee sidecar cleanup even on error.
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

"""
    build_protein_global_score_dicts(pg_refs, sqrt_n_runs, n_proteins)
    → (global_pg_score_dict, pg_name_to_global_pg_score)

Stream per-file protein group files reading only (protein_name, target, entrap_id, pg_score).
Returns composite-key dictionaries for protein global score computation.

`n_proteins` is used for `sizehint!()` pre-allocation,
obtained via `length(getProteins(getSpecLib(search_context)))`.
"""
function build_protein_global_score_dicts(
    pg_refs::Vector{ProteinGroupFileReference},
    sqrt_n_runs::Int,
    n_proteins::Int
)
    # Pre-allocate accumulation dictionary with known upper bound
    score_acc = Dict{Tuple{String,Bool,UInt8}, Vector{Float32}}()
    sizehint!(score_acc, n_proteins)

    for ref in pg_refs
        tbl = Arrow.Table(file_path(ref))
        n = length(tbl.protein_name)
        n == 0 && continue
        @inbounds for i in 1:n
            key = (tbl.protein_name[i], tbl.target[i], tbl.entrap_id[i])
            if !haskey(score_acc, key)
                score_acc[key] = Float32[]
            end
            push!(score_acc[key], tbl.pg_score[i])
        end
    end

    # Compute global scores via logodds
    n_observed = length(score_acc)
    global_pg_score_dict = Dict{Tuple{String,Bool,UInt8}, Float32}()
    sizehint!(global_pg_score_dict, n_observed)
    pg_name_to_global_pg_score = Dict{ProteinKey, Float32}()
    sizehint!(pg_name_to_global_pg_score, n_observed)

    for (key, scores) in score_acc
        gs = logodds(scores, sqrt_n_runs)
        global_pg_score_dict[key] = gs
        pg_name_to_global_pg_score[ProteinKey(key[1], key[2], key[3])] = gs
    end

    return global_pg_score_dict, pg_name_to_global_pg_score
end

"""
    build_protein_global_qval_dict(global_pg_score_dict)
    → Dict{Tuple{String,Bool,UInt8}, Float32}

Compute protein global q-values directly from score dictionary.
Replaces get_protein_global_qval_dict (which loaded a full merged DataFrame).
"""
function build_protein_global_qval_dict(
    global_pg_score_dict::Dict{Tuple{String,Bool,UInt8}, Float32}
)
    n = length(global_pg_score_dict)
    keys_vec = collect(keys(global_pg_score_dict))
    scores = Float32[global_pg_score_dict[k] for k in keys_vec]
    targets = Bool[k[2] for k in keys_vec]

    # Sort by (score desc, target desc) for proper FDR
    perm = sortperm(collect(zip(scores, targets)); by=x->(-x[1], -x[2]))
    permute!(keys_vec, perm)
    permute!(scores, perm)
    permute!(targets, perm)

    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals)

    qval_dict = Dict{Tuple{String,Bool,UInt8}, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[keys_vec[i]] = qvals[i]
    end
    return qval_dict
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
        get_qvalues!(df.trace_prob, df.target, qvals; fdr_scale_factor=fdr_scale_factor)
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
