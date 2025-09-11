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

#==========================================================
PSM sampling and scoring 
==========================================================#

# Constant for model selection threshold
const MAX_FOR_MODEL_SELECTION = 200_000
# Backward-compatible alias used in discussions/documentation
const MAX_FOR_MODEL_SELECTION_PSMS = MAX_FOR_MODEL_SELECTION

"""
    score_precursor_isotope_traces(second_pass_folder::String, 
                                  file_paths::Vector{String},
                                  precursors::LibraryPrecursors,
                                  match_between_runs::Bool,
                                  max_q_value_xgboost_rescore::Float32,
                                  max_q_value_xgboost_mbr_rescore::Float32,
                                  min_PEP_neg_threshold_xgboost_rescore::Float32,
                                  max_psms_in_memory::Int64)

Main entry point for PSM scoring with automatic model selection based on dataset size.

# Three-Case Logic
1. PSMs ≥ max_psms_in_memory: Out-of-memory processing with default XGBoost
2. PSMs < max_psms_in_memory AND ≥ 100K: In-memory with default/advanced XGBoost
3. PSMs < 100K: In-memory with automatic model comparison

# Arguments
- `second_pass_folder`: Folder containing second pass PSM files
- `file_paths`: Vector of PSM file paths
- `precursors`: Library precursors
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_xgboost_rescore`: Max q-value for EvoTrees/XGBoost rescoring
- `max_q_value_xgboost_mbr_rescore`: Max q-value for MBR rescoring
- `min_PEP_neg_threshold_xgboost_rescore`: Min PEP for negative relabeling
- `max_psms_in_memory`: Maximum PSMs to keep in memory

# Returns
- Trained EvoTrees/XGBoost models or nothing for probit regression
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    max_psms_in_memory::Int64,
    q_value_threshold::Float32 = 0.01f0  # Default to 1% if not specified
)
    # Step 1: Count PSMs and determine processing approach
    psms_count = get_psms_count(file_paths)
    
    if psms_count >= max_psms_in_memory
        # Case 1: Out-of-memory processing with default XGBoost
        @user_info "Using out-of-memory processing for $psms_count PSMs (≥ $max_psms_in_memory)"
        best_psms = sample_psms_for_xgboost(second_pass_folder, psms_count, max_psms_in_memory)
        # Use a ModelConfig (AdvancedXGBoost by default) for OOM path
        model_config = create_default_advanced_xgboost_config()
        models = score_precursor_isotope_traces_out_of_memory!(
            best_psms,
            file_paths,
            precursors,
            model_config,
            match_between_runs,
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore
        )
    else
        # In-memory processing - load PSMs first
        best_psms = load_psms_for_xgboost(second_pass_folder)
        
        if psms_count >= MAX_FOR_MODEL_SELECTION  # 100K
            # Case 2: In-memory with default/advanced XGBoost (no comparison)
            @user_info "Using in-memory advanced XGBoost for $psms_count PSMs (< $max_psms_in_memory but ≥ 100K)"
            model_config = create_default_advanced_xgboost_config()
        else
            # Case 3: In-memory with automatic model comparison (<100K)
            model_config = select_psm_scoring_model(
                best_psms, file_paths, precursors, match_between_runs,
                max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
                min_PEP_neg_threshold_xgboost_rescore, q_value_threshold
            )
        end
        
        # Execute selected model using unified function
        @user_info "Training final model: $(model_config.name)"
        models = score_precursor_isotope_traces_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
        )
        
        # Write scored PSMs to files
        write_scored_psms_to_files!(best_psms, file_paths)
    end
    
    # Clean up
    best_psms = nothing
    GC.gc()
    
    return models
end

"""
    select_psm_scoring_model(best_psms::DataFrame, ...) -> ModelConfig

Selects the appropriate PSM scoring model based on dataset size and characteristics.

# Model Selection Logic
- ≥100K PSMs: Returns default/advanced XGBoost configuration (no comparison)
- <100K PSMs: Trains each model and selects based on training performance

# Returns
- ModelConfig object specifying the selected model and its hyperparameters
"""
function select_psm_scoring_model(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    q_value_threshold::Float32
)
    psms_count = size(best_psms, 1)
    
    if psms_count >= MAX_FOR_MODEL_SELECTION
        @user_info "Using default advanced XGBoost for $psms_count PSMs (≥ 100K)"
        return create_default_advanced_xgboost_config()
    else
        # Get model configurations from model_config.jl
        model_configs = create_model_configurations()
        best_model_config = nothing
        best_target_count = 0
        
        @user_info "Model comparison: Testing $(length(model_configs)) models on $psms_count PSMs"
        @user_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        for config in model_configs
            # Create deepcopy to avoid side effects between models
            psms_copy = deepcopy(best_psms)
            
            try
                # Train model with suppressed output during comparison
                # Redirect both stdout AND stderr to suppress all progress bars
                score_precursor_isotope_traces_in_memory(
                    psms_copy, file_paths, precursors, config,
                    match_between_runs, max_q_value_xgboost_rescore,
                    max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore,
                    false  # show_progress = false during comparison
                )
                
                # Count passing targets from DataFrame
                target_count = count_passing_targets(psms_copy, q_value_threshold)
                @user_info "  $(config.name): $(target_count) IDs at q ≤ $q_value_threshold"
                
                # Update best model if this one is better
                if target_count > best_target_count
                    best_target_count = target_count
                    best_model_config = config
                end
            catch e
                @user_warn "Failed to train $(config.name):"
                # Show error type and message without data
                if isa(e, MethodError)
                    @user_warn "  MethodError: $(e.f) with $(length(e.args)) arguments"
                    @user_warn "  Argument types: $(typeof.(e.args))"
                else
                    @user_warn "  $(typeof(e)): $(sprint(showerror, e))"
                end
                # Show just the top of stack trace
                st = stacktrace(catch_backtrace())
                if length(st) > 0
                    @user_warn "  at $(st[1].file):$(st[1].line)"
                end
                # Continue with next model
            end
        end
        
        @user_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        # If no model succeeded, fall back to SimpleXGBoost
        if best_model_config === nothing
            @user_warn "All models failed, defaulting to SimpleXGBoost"
            best_model_config = model_configs[1]  # SimpleXGBoost is first
        else
            @user_info "✓ Selected: $(best_model_config.name) ($(best_target_count) IDs)"
        end
        
        return best_model_config
    end
end

"""
    count_passing_targets(scored_psms::DataFrame, qvalue_threshold::Float64) -> Int

Counts the number of target PSMs passing the q-value threshold.

# Arguments
- `scored_psms`: DataFrame containing scored PSMs with :prob and :target columns
- `qvalue_threshold`: Q-value threshold for counting (typically 0.01)

# Returns
- Total count of target PSMs with q_value ≤ threshold
"""
function count_passing_targets(scored_psms::DataFrame, qvalue_threshold::Float32)
    # Use existing q-value calculation logic
    if :prob in propertynames(scored_psms) && :target in propertynames(scored_psms)
        # Calculate q-values from probabilities
        qvals = Vector{Float32}(undef, nrow(scored_psms))
        get_qvalues!(scored_psms.prob, scored_psms.target, qvals)
        
        # Count targets passing threshold
        return sum((scored_psms.target .== true) .& (qvals .<= qvalue_threshold))
    else
        error("DataFrame must contain :prob and :target columns")
    end
end

"""
    create_default_advanced_xgboost_config() -> ModelConfig

Creates the default advanced XGBoost configuration for large datasets.
"""
function create_default_advanced_xgboost_config()
    return ModelConfig(
        "AdvancedXGBoost",
        :xgboost,
        ADVANCED_FEATURE_SET,
        Dict(
            :colsample_bytree => 0.5,
            :colsample_bynode => 0.5,
            :min_child_weight => 5,
            :gamma => 1.0,
            :subsample => 0.25,
            :max_depth => 10,
            :eta => 0.05,
            :iter_scheme => [100, 200, 200]
        )
    )
end

"""
    score_precursor_isotope_traces_in_memory(psms::DataFrame, ...) -> Models

Unified in-memory PSM scoring function that executes the specified model.

# Arguments
- `model_config`: ModelConfig specifying which model to use and its hyperparameters

# Supported Models
- SimpleXGBoost: Full feature set, standard hyperparameters
- AdvancedXGBoost: Full feature set, advanced hyperparameters
- ProbitRegression: Linear probit model with CV folds from library
- SuperSimplified: Minimal 5-feature XGBoost model
"""
function score_precursor_isotope_traces_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    show_progress::Bool = true
)
    
    if model_config.model_type == :xgboost
        return train_xgboost_model_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore,
            show_progress
        )
    elseif model_config.model_type == :probit
        return train_probit_model_in_memory(
            best_psms, file_paths, precursors, model_config, match_between_runs,
            min_PEP_neg_threshold_xgboost_rescore
        )
    else
        error("Unsupported model type: $(model_config.model_type)")
    end
end

"""
    train_xgboost_model_in_memory(...) -> Models

Trains XGBoost model using configuration from ModelConfig.
"""
function train_xgboost_model_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    show_progress::Bool = true
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features and hyperparams from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    if match_between_runs
        append!(features, [
            :MBR_num_runs, 
            :MBR_max_pair_prob,
            :MBR_log2_weight_ratio, 
            :MBR_log2_explained_ratio,
            :MBR_rv_coefficient, 
            :MBR_best_irt_diff,
            :MBR_is_missing
        ])
    end
    
    hp = model_config.hyperparams
    return sort_of_percolator_in_memory!(
        best_psms, features, match_between_runs;
        max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
        min_PEP_neg_threshold_xgboost_rescore,
        colsample_bytree = get(hp, :colsample_bytree, 0.5),
        min_child_weight = get(hp, :min_child_weight, 5),
        gamma = get(hp, :gamma, 1.0),
        subsample = get(hp, :subsample, 0.25),
        max_depth = get(hp, :max_depth, 10),
        eta = get(hp, :eta, 0.05),
        iter_scheme = get(hp, :iter_scheme, [100, 200, 200]),
        show_progress = show_progress
    )
end

"""
    train_probit_model_in_memory(...) -> Nothing

Trains probit regression model using configuration from ModelConfig.
"""
function train_probit_model_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    
    probit_regression_scoring_cv!(
        best_psms,
        file_paths,
        features,
        match_between_runs;
        neg_mining_pep_threshold = min_PEP_neg_threshold_xgboost_rescore
    )
    
    # File writing removed - will be done at higher level
    
    return nothing  # Probit doesn't return models
end

"""
     get_psms_count(quant_psms_folder::String)::Integer

Sample PSMs from multiple files for EvoTrees/XGBoost model training.

# Arguments
- `quant_psms_folder`: Folder containing PSM Arrow files

# Process
1. Counts total PSMs across .arrow files in the directory
"""
function get_psms_count(file_paths::Vector{String})

    #file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]

    psms_count = 0

    for file_path in file_paths
        psms_count += length(Arrow.Table(file_path)[1])
    end

    return psms_count
end


"""
    sample_psms_for_xgboost(quant_psms_folder::String, psms_count::Integer max_psms::Integer) -> DataFrame

Sample PSMs from multiple files for EvoTrees/XGBoost model training.

# Arguments
- `quant_psms_folder`: Folder containing PSM Arrow files
- `psms_count`: number of psms across all the arrow files
- `max_psms`: Maximum number of PSMs to sample for training

# Process
1. Proportionally samples from each file
2. Combines samples into single DataFrame
"""
function sample_psms_for_xgboost(quant_psms_folder::String, psms_count::Integer, max_psms::Integer)

    file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]

    # Initialize an empty DataFrame to store the results
    result_df = DataFrame()

    for file_path in file_paths
        # Read the Arrow table
        arrow_table = Arrow.Table(file_path)
        
        # Get the number of rows
        num_rows = length(arrow_table[1])
        
        # Calculate the number of rows to sample (1/N'th of the total)
        sample_size = min(ceil(Int, (num_rows/psms_count)*max_psms), num_rows) #ceil(Int, num_rows / N)

        # Generate sorted random indices for sampling
        sampled_indices = sort!(sample(MersenneTwister(1776), 1:num_rows, sample_size, replace=false))
        
        # Sample the rows and convert to DataFrame
        sampled_df = DataFrame(arrow_table)[sampled_indices, :]
        
        # Append to the result DataFrame
        append!(result_df, sampled_df)
    end

    return result_df
end

"""
     get_psms_count(quant_psms_folder::String)::Integer

Loads all PSMs from multiple files for EvoTrees/XGBoost model training.
"""
function load_psms_for_xgboost(quant_psms_folder::String)
    file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]
    return DataFrame(Tables.columntable(Arrow.Table(file_paths)))
end

"""
    score_precursor_isotope_traces_in_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::LibraryPrecursors) -> EvoTreesModels

Train EvoTrees/XGBoost models for PSM scoring. All psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained EvoTrees/XGBoost models or simplified model if insufficient PSMs.
"""

"""
    probit_regression_scoring_cv!(psms::DataFrame,
                                  file_paths::Vector{String},
                                  features::Vector{Symbol},
                                  match_between_runs::Bool;
                                  n_folds::Int64 = 3,
                                  neg_mining_pep_threshold::Float32 = 0.90f0)

Alternative PSM scoring using probit regression with cross-validation.

This is a simpler alternative to XGBoost/EvoTrees for small datasets (<100k PSMs).
Uses linear probit model with cross-validation, similar to FirstPassSearch but with CV folds.
No iterative refinement or max_prob updates - single pass training only.

# Arguments
- `psms`: DataFrame containing PSMs to score
- `file_paths`: Vector of file paths for CV fold assignment
- `features`: Feature columns to use for scoring (same as XGBoost)
- `match_between_runs`: Whether MBR was performed
- `n_folds`: Number of cross-validation folds (default: 3)

# Modifies
- Adds columns to psms: `:prob`, `:q_value`, `:best_psm`, `:cv_fold`
"""
function probit_regression_scoring_cv!(
    psms::DataFrame,
    file_paths::Vector{String},
    features::Vector{Symbol},
    match_between_runs::Bool;
    neg_mining_pep_threshold::Float32 = 0.90f0
)
    # Step 1: Validate CV fold column exists (from library assignments)
    if !(:cv_fold in propertynames(psms))
        error("PSMs must have :cv_fold column from library assignments")
    end
    
    # Detect unique CV folds from data (should be [0, 1] from library)
    unique_cv_folds = sort(unique(psms[!, :cv_fold]))
    n_folds = length(unique_cv_folds)
    
    # Validate we have expected 2 folds with values 0 and 1
    if unique_cv_folds != UInt8[0, 1]
        @user_warn "Unexpected CV folds: $unique_cv_folds (expected [0, 1] from library)"
    end
    
    # Step 2: Initialize probability array (like XGBoost does)
    # Use a separate array instead of directly modifying DataFrame
    prob_estimates = zeros(Float32, size(psms, 1))
    
    # Step 3: Use the exact features passed in (already filtered by the caller)
    # Add intercept column if not present
    if !(:intercept in propertynames(psms))
        psms[!, :intercept] = ones(Float32, size(psms, 1))
    end
    
    # Filter features to those that exist in the DataFrame
    available_features = filter(col -> col in propertynames(psms), features)
    if :intercept ∉ available_features
        push!(available_features, :intercept)
    end
    
    # Step 4: Train probit model per CV fold (two-pass with 1% q-value and optional negative mining)
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = Iterators.partition(1:size(psms, 1), chunk_size)
    
    for fold_idx in unique_cv_folds
        # Get train/test masks (using actual fold values from library)
        test_mask = psms[!, :cv_fold] .== fold_idx
        train_mask = .!test_mask
        
        # Check for sufficient training data
        train_targets = train_mask .& psms[!, :target]
        train_decoys = train_mask .& .!psms[!, :target]
        
        n_train_targets = sum(train_targets)
        n_train_decoys = sum(train_decoys)
        
        if n_train_targets < 100 || n_train_decoys < 100
            @user_warn "Insufficient training data for fold $fold_idx (targets: $n_train_targets, decoys: $n_train_decoys)"
            psms[test_mask, :prob] .= 0.5f0
            continue
        end
        
        # Get training data
        train_data = psms[train_mask, :]
        train_targets_bool = train_data[!, :target]
        
        # Prepare chunks for training data
        train_chunk_size = max(1, size(train_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
        train_chunks = Iterators.partition(1:size(train_data, 1), train_chunk_size)
        
        # Initialize coefficients
        β = zeros(Float64, length(available_features))
        
        # Pass 1: Fit probit model using existing ProbitRegression on full training set
        β = Pioneer.ProbitRegression(
            β, 
            train_data[!, available_features], 
            train_targets_bool, 
            train_chunks, 
            max_iter = 30
        )
        # Coefficient diagnostics removed per request
        
        # Predict probabilities on train and test sets (needed for selection and scoring)
        train_probs = zeros(Float32, size(train_data, 1))
        Pioneer.ModelPredictProbs!(train_probs, train_data[!, available_features], β, train_chunks)
        
        # Predict probabilities on test set
        test_data = psms[test_mask, :]
        test_probs = zeros(Float32, size(test_data, 1))
        
        test_chunk_size = max(1, size(test_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
        test_chunks = Iterators.partition(1:size(test_data, 1), test_chunk_size)
        
        # Use ModelPredictProbs! to get probabilities
        Pioneer.ModelPredictProbs!(
            test_probs,
            test_data[!, available_features],
            β,
            test_chunks
        )
        prob_estimates[test_mask] = test_probs

        # Compute train q-values and PEPs for selection
        qvals_train = Vector{Float32}(undef, size(train_data, 1))
        get_qvalues!(train_probs, train_targets_bool, qvals_train)
        pep_train = Vector{Float32}(undef, size(train_data, 1))
        get_PEP!(train_probs, train_targets_bool, pep_train)

        # Build second-pass training set
        # Hard-coded 1% q-value threshold for selecting positive targets
        q_thresh = 0.01f0
        confident_pos = (train_targets_bool .== true) .& (qvals_train .<= q_thresh)
        decoys_train = (train_targets_bool .== false)
        # Negative mining: re-label worst targets as negatives using provided PEP threshold
        mined_negs = (train_targets_bool .== true) .& (pep_train .>= neg_mining_pep_threshold)

        second_mask_rel = decoys_train .| confident_pos .| mined_negs
        if any(second_mask_rel)
            second_train_data = train_data[second_mask_rel, :]
            second_labels = copy(train_targets_bool[second_mask_rel])
            # Relabel mined negatives to false
            relabel_idx = findall(mined_negs[second_mask_rel])
            for i in relabel_idx
                second_labels[i] = false
            end

            # Fit pass 2 model
            second_chunk_size = max(1, size(second_train_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
            second_chunks = Iterators.partition(1:size(second_train_data, 1), second_chunk_size)
            β2 = zeros(Float64, length(available_features))
            β2 = Pioneer.ProbitRegression(
                β2,
                second_train_data[!, available_features],
                second_labels,
                second_chunks,
                max_iter = 30
            )
            # Coefficient diagnostics removed per request
            # Re-predict on test set with pass 2 model
            Pioneer.ModelPredictProbs!(
                test_probs,
                test_data[!, available_features],
                β2,
                test_chunks
            )
            prob_estimates[test_mask] = test_probs
        end
    end
    
    # Step 5: Assign probabilities to DataFrame (like XGBoost does at line 105 of percolatorSortOf.jl)
    psms[!, :prob] = prob_estimates
    
    # If MBR is enabled, create MBR columns (probit doesn't do separate MBR scoring)
    # This is needed for compatibility with apply_mbr_filter! in ScoringSearch
    if match_between_runs
        psms[!, :MBR_prob] = copy(prob_estimates)
        # MBR_is_best_decoy is required by apply_mbr_filter! 
        # For probit, we don't have MBR transfer info, so set conservatively
        psms[!, :MBR_is_best_decoy] = fill(missing, size(psms, 1))  # missing means not an MBR transfer
    end
    
    # Ensure prob column has proper type and no NaN/Inf values
    # Check for any NaN or Inf values and replace them
    n_nan = sum(isnan.(psms.prob))
    n_inf = sum(isinf.(psms.prob))
    if n_nan > 0 || n_inf > 0
        @user_warn "Found $n_nan NaN and $n_inf Inf values in probabilities, replacing with 0.5"
        psms.prob[isnan.(psms.prob) .| isinf.(psms.prob)] .= 0.5f0
    end
    
    # Ensure values are in valid range but avoid values too close to 0 or 1
    # The downstream aggregation formula breaks with probabilities > 0.999999
    # Use more conservative clamping to ensure we stay below the threshold
    # Force Float32 to ensure proper clamping
    psms.prob = Float32.(clamp.(psms.prob, Float32(1e-6), Float32(0.9999)))
    
    # Clean up columns added during probit regression that shouldn't be persisted
    if :intercept in propertynames(psms)
        select!(psms, Not(:intercept))
    end
    if :cv_fold in propertynames(psms)
        select!(psms, Not(:cv_fold))
    end
    
    
    return nothing
end


"""
    score_precursor_isotope_traces_out_of_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::LibraryPrecursors) -> EvoTreesModels

Train EvoTrees/XGBoost models for PSM scoring. Only a subset of psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained EvoTrees/XGBoost models or simplified model if insufficient PSMs.
"""
function score_precursor_isotope_traces_out_of_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    # Features from model_config; do not include :target
    features = [f for f in model_config.features if hasproperty(best_psms, f)];
    if match_between_runs
        append!(features, [
            :MBR_rv_coefficient,
            :MBR_best_irt_diff,
            :MBR_num_runs,
            :MBR_max_pair_prob,
            :MBR_log2_weight_ratio,
            :MBR_log2_explained_ratio,
            :MBR_is_missing
            ])
    end

    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
    best_psms[!,:decoy] = best_psms[!,:target].==false;

    # Hyperparameters from model_config
    hp = model_config.hyperparams
    models = sort_of_percolator_out_of_memory!(
                            best_psms, 
                            file_paths,
                            features,
                            match_between_runs;
                            max_q_value_xgboost_rescore,
                            max_q_value_xgboost_mbr_rescore,
                            min_PEP_neg_threshold_xgboost_rescore,
                            colsample_bytree = get(hp, :colsample_bytree, 0.5), 
                            min_child_weight = get(hp, :min_child_weight, 5), 
                            gamma = get(hp, :gamma, 1.0),
                            subsample = get(hp, :subsample, 0.25), 
                            max_depth = get(hp, :max_depth, 10),
                            eta = get(hp, :eta, 0.05), 
                            iter_scheme = get(hp, :iter_scheme, [100, 200, 200]),
                            print_importance = false);
    return models;#best_psms
end

"""
    write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})

Write scored PSMs back to Arrow files, grouped by ms_file_idx.
This function is separated from scoring to allow model comparison without file I/O.

# Arguments
- `psms`: DataFrame containing scored PSMs with ms_file_idx column
- `file_paths`: Vector of file paths for valid files only
"""
function write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})
    dropVectorColumns!(psms) # avoids writing issues
    
    # Create mapping from unique ms_file_idx values to file paths
    unique_file_indices = unique(psms[:, :ms_file_idx])
    sort!(unique_file_indices)
    
    # Check that we have enough file paths for all file indices
    if length(file_paths) != length(unique_file_indices)
        error("Mismatch: $(length(file_paths)) file paths provided but $(length(unique_file_indices)) unique file indices found in PSM data")
    end
    
    # Create mapping: original file index → output file path
    index_to_path = Dict(zip(unique_file_indices, file_paths))
    
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        file_idx = ms_file_idx[:ms_file_idx]
        if haskey(index_to_path, file_idx)
            fpath = index_to_path[file_idx]
            writeArrow(fpath, gpsms)
        else
            @warn "No output path found for file index $file_idx, skipping"
        end
    end
end
