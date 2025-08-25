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
const MAX_FOR_MODEL_SELECTION = 100_000

"""
    score_precursor_isotope_traces(second_pass_folder::String, 
                                  file_paths::Vector{String},
                                  precursors::LibraryPrecursors,
                                  match_between_runs::Bool,
                                  max_q_value_xgboost_rescore::Float32,
                                  max_q_value_xgboost_mbr_rescore::Float32,
                                  min_PEP_neg_threshold_xgboost_rescore::Float32,
                                  q_value_threshold::Float32,
                                  max_MBR_false_transfer_rate::Float32,
                                  max_psms_in_memory::Int64)

Main entry point for PSM scoring. Automatically chooses between in-memory and out-of-memory
processing based on data size.

# Arguments
- `second_pass_folder`: Folder containing second pass PSM files
- `file_paths`: Vector of PSM file paths
- `precursors`: Library precursors
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_xgboost_rescore`: Max q-value for EvoTrees/XGBoost rescoring
- `max_q_value_xgboost_mbr_rescore`: Max q-value for MBR rescoring
- `min_PEP_neg_threshold_xgboost_rescore`: Min PEP for negative relabeling
- `q_value_threshold`: Max q-value for final output
- `max_MBR_false_transfer_rate`: Max FTR for MBR
- `max_psms_in_memory`: Maximum PSMs to keep in memory

# Returns
- Trained EvoTrees/XGBoost models
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    max_psms_in_memory::Int64;
    enable_model_comparison::Bool = false,
    validation_split_ratio::Float64 = 0.2,
    qvalue_threshold::Float64 = 0.01,
    min_psms_for_comparison::Int = 1000,
    max_psms_for_comparison::Int = 100000,
    output_dir::String = "."
)
    # Count total PSMs
    psms_count = get_psms_count(file_paths)
    
    if psms_count > max_psms_in_memory
        # Use out-of-memory algorithm
        best_psms = sample_psms_for_xgboost(second_pass_folder, psms_count, max_psms_in_memory)
        models = score_precursor_isotope_traces_out_of_memory!(
            best_psms,
            file_paths,
            precursors,
            match_between_runs,
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore
        )
    else
        # Use in-memory algorithm
        best_psms = load_psms_for_xgboost(second_pass_folder)
        
        # Check if model comparison should be enabled
        if enable_model_comparison && 
           psms_count >= min_psms_for_comparison && 
           psms_count <= max_psms_for_comparison
            
            models = score_precursor_isotope_traces_in_memory_with_comparison!(
                best_psms,
                file_paths,
                precursors,
                match_between_runs,
                max_q_value_xgboost_rescore,
                max_q_value_xgboost_mbr_rescore,
                min_PEP_neg_threshold_xgboost_rescore,
                enable_model_comparison,
                validation_split_ratio,
                qvalue_threshold,
                output_dir
            )
            
            # If model comparison fails, fall back to existing method
            if models === nothing
                @user_warn "Model comparison failed, falling back to standard approach"
                models = score_precursor_isotope_traces_in_memory!(
                    best_psms,
                    file_paths,
                    precursors,
                    match_between_runs,
                    max_q_value_xgboost_rescore,
                    max_q_value_xgboost_mbr_rescore,
                    min_PEP_neg_threshold_xgboost_rescore
                )
            end
        else
            # Use standard approach
            models = score_precursor_isotope_traces_in_memory!(
                best_psms,
                file_paths,
                precursors,
                match_between_runs,
                max_q_value_xgboost_rescore,
                max_q_value_xgboost_mbr_rescore,
                min_PEP_neg_threshold_xgboost_rescore
            )
        end
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
- <100K PSMs: Runs model comparison and returns best performing model

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
    validation_split_ratio::Float64,
    qvalue_threshold::Float64,
    output_dir::String
)
    psms_count = size(best_psms, 1)
    
    if psms_count >= MAX_FOR_MODEL_SELECTION
        @user_info "Using default advanced XGBoost for $psms_count PSMs (≥ 100K)"
        return create_default_advanced_xgboost_config()
    else
        @user_info "Running model comparison for $psms_count PSMs (< 100K)"
        
        # Use the existing model comparison framework from model_comparison.jl
        models = score_precursor_isotope_traces_in_memory_with_comparison!(
            best_psms,
            file_paths,
            precursors,
            match_between_runs,
            max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore,
            min_PEP_neg_threshold_xgboost_rescore,
            true,  # enable_model_comparison
            validation_split_ratio,
            qvalue_threshold,
            output_dir
        )
        
        # For now, return the SimpleXGBoost config as default
        # In future, this should extract the best model config from the comparison
        return ModelConfig(
            "SimpleXGBoost",
            :xgboost,
            REDUCED_FEATURE_SET,
            Dict(
                :colsample_bytree => 0.8,
                :colsample_bynode => 0.8,
                :min_child_weight => 20,
                :gamma => 0.1,
                :subsample => 0.8,
                :max_depth => 4,
                :eta => 0.1,
                :iter_scheme => [150, 300, 300]
            )
        )
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
        REDUCED_FEATURE_SET,
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
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    @user_info "Running $(model_config.name) model for PSM scoring"
    
    if model_config.model_type == :xgboost
        return train_xgboost_model_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_xgboost_rescore,
            max_q_value_xgboost_mbr_rescore, min_PEP_neg_threshold_xgboost_rescore
        )
    elseif model_config.model_type == :probit
        return train_probit_model_in_memory(
            best_psms, file_paths, precursors, model_config, match_between_runs
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
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features and hyperparams from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    if match_between_runs
        append!(features, [
            :MBR_rv_coefficient, :MBR_best_irt_diff, :MBR_num_runs,
            :MBR_max_pair_prob, :MBR_log2_weight_ratio, :MBR_log2_explained_ratio
        ])
    end
    
    hp = model_config.hyperparams
    return sort_of_percolator_in_memory!(
        best_psms, file_paths, features, match_between_runs;
        max_q_value_xgboost_rescore, max_q_value_xgboost_mbr_rescore,
        min_PEP_neg_threshold_xgboost_rescore,
        colsample_bytree = get(hp, :colsample_bytree, 0.5),
        colsample_bynode = get(hp, :colsample_bynode, 0.5),
        min_child_weight = get(hp, :min_child_weight, 5),
        gamma = get(hp, :gamma, 1.0),
        subsample = get(hp, :subsample, 0.25),
        max_depth = get(hp, :max_depth, 10),
        eta = get(hp, :eta, 0.05),
        iter_scheme = get(hp, :iter_scheme, [100, 200, 200])
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
    match_between_runs::Bool
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
        match_between_runs
    )
    
    # Write results back to files
    dropVectorColumns!(best_psms)
    for (ms_file_idx, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    
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
                                  n_folds::Int64 = 3)

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
    match_between_runs::Bool
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
    
    # Verify fold distribution
    for fold_idx in unique_cv_folds
        n_in_fold = sum(psms.cv_fold .== fold_idx)
        @debug_l1 "CV fold $fold_idx: $n_in_fold PSMs"
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
    push!(available_features, :intercept)
    
    @debug_l1 "Probit regression using $(length(available_features)-1) features (plus intercept)"
    
    # Step 4: Train probit model per CV fold
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
        
        @debug_l1 "Training probit model for fold $fold_idx with $n_train_targets targets and $n_train_decoys decoys"
        
        # Get training data
        train_data = psms[train_mask, :]
        train_targets_bool = train_data[!, :target]
        
        # Prepare chunks for training data
        train_chunk_size = max(1, size(train_data, 1) ÷ (tasks_per_thread * Threads.nthreads()))
        train_chunks = Iterators.partition(1:size(train_data, 1), train_chunk_size)
        
        # Initialize coefficients
        β = zeros(Float64, length(available_features))
        
        # Fit probit model using existing ProbitRegression
        β = Pioneer.ProbitRegression(
            β, 
            train_data[!, available_features], 
            train_targets_bool, 
            train_chunks, 
            max_iter = 30
        )
        
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
        @debug_l1 "Created MBR columns for compatibility (MBR_prob same as prob for probit)"
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
    
    # Debug: Verify clamping worked
    @debug_l1 "After clamping: prob range = $(minimum(psms.prob)) to $(maximum(psms.prob))"
    
    Arrow.write("/Users/nathanwamsley/Desktop/test_arrow_psms.arrow", psms)
    @debug_l1 "Probit regression scoring complete. Probabilities assigned to $(size(psms, 1)) PSMs"
    
    return nothing
end

function score_precursor_isotope_traces_in_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    if size(best_psms, 1) > 1#100_000
        file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
        features = [
            :missed_cleavage,
            :Mox,
            :prec_mz,
            :sequence_length,
            :charge,
            :irt_pred,
            :irt_error,
            :irt_diff,
            :max_y_ions,
            :y_ions_sum,
            :longest_y,
            :y_count,
            :b_count,
            :isotope_count,
            :total_ions,
            :best_rank,
            :best_rank_iso,
            :topn,
            :topn_iso,
            :gof,
            :max_fitted_manhattan_distance,
            :max_fitted_spectral_contrast,
            :max_matched_residual,
            :max_unmatched_residual,
            :max_gof,
            :fitted_spectral_contrast,
            :spectral_contrast,
            :max_matched_ratio,
            :err_norm,
            :poisson,
            :weight,
            :log2_intensity_explained,
            :tic,
            :num_scans,
            :smoothness,
            :rt_diff,
            :ms1_irt_diff,
            :weight_ms1,
            :gof_ms1,
            :max_matched_residual_ms1,
            :max_unmatched_residual_ms1,
            :fitted_spectral_contrast_ms1,
            :error_ms1,
            :m0_error_ms1,
            :n_iso_ms1,
            :big_iso_ms1,
            :ms1_features_missing,
            :percent_theoretical_ignored,
            :scribe,
            :max_scribe
        ];
        features = [f for f in features if hasproperty(best_psms, f)];
        if match_between_runs
            append!(features, [
                :MBR_rv_coefficient,
                :MBR_best_irt_diff,
                :MBR_num_runs,
                :MBR_max_pair_prob,
                :MBR_log2_weight_ratio,
                :MBR_log2_explained_ratio
                ])
        end

        
        best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
        best_psms[!,:decoy] = best_psms[!,:target].==false;

        models = sort_of_percolator_in_memory!(
                                best_psms, 
                                file_paths,
                                features,
                                match_between_runs;
                                max_q_value_xgboost_rescore,
                                max_q_value_xgboost_mbr_rescore,
                                min_PEP_neg_threshold_xgboost_rescore,
                                colsample_bytree = 0.5, 
                                colsample_bynode = 0.5,
                                min_child_weight = 5, 
                                gamma = 1.0,
                                subsample = 0.25, 
                                max_depth = 10,
                                eta = 0.05, 
                                iter_scheme = [100, 200, 200],
                                print_importance = false);
        return models;#best_psms
    else
        @user_warn "Less than 100,000 psms: $(size(best_psms, 1)). Training with simplified target-decoy discrimination model..."
        if size(best_psms, 1) > 1000
            @user_warn "Less than 100,000 psms. Training with simplified target-decoy discrimination model..."
            features = [ 
                :missed_cleavage,
                :Mox,
                :sequence_length,
                :charge,
                :irt_error,
                :irt_diff,
                :y_count,
                :max_fitted_manhattan_distance,
                :max_matched_residual,
                :max_unmatched_residual,
                :max_gof,
                :fitted_spectral_contrast,
                :spectral_contrast,
                :err_norm,
                :weight,
                :log2_intensity_explained,
                :percent_theoretical_ignored,
                :scribe,
            ];

                    features = [
            :missed_cleavage,
            :Mox,
            :prec_mz,
            :sequence_length,
            :charge,
            :irt_pred,
            :irt_error,
            :irt_diff,
            :max_y_ions,
            :y_ions_sum,
            :longest_y,
            :y_count,
            :b_count,
            :isotope_count,
            :total_ions,
            :best_rank,
            :best_rank_iso,
            :topn,
            :topn_iso,
            :gof,
            :max_fitted_manhattan_distance,
            :max_fitted_spectral_contrast,
            :max_matched_residual,
            :max_unmatched_residual,
            :max_gof,
            :fitted_spectral_contrast,
            :spectral_contrast,
            :max_matched_ratio,
            :err_norm,
            :poisson,
            :weight,
            :log2_intensity_explained,
            :tic,
            :num_scans,
            :smoothness,
            :rt_diff,
            :ms1_irt_diff,
            :weight_ms1,
            :gof_ms1,
            :max_matched_residual_ms1,
            :max_unmatched_residual_ms1,
            :fitted_spectral_contrast_ms1,
            :error_ms1,
            :m0_error_ms1,
            :n_iso_ms1,
            :big_iso_ms1,
            :ms1_features_missing,
            :percent_theoretical_ignored,
            :scribe,
            :max_scribe
        ];

        else
             @user_warn "Less than 1,000 psms. Training with super simplified target-decoy discrimination model..."
             features = [ 
                :fitted_spectral_contrast,
                :max_matched_residual,
                :max_unmatched_residual,
                :err_norm,
                :log2_intensity_explained,
            ];
        end
        file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
        
        best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
        best_psms[!,:decoy] = best_psms[!,:target].==false;
        
        # OPTION 1: Probit regression (SIMPLE, NO ITERATIVE REFINEMENT)
        # @info "Using probit regression for small dataset (<100k PSMs)"
        #=
        probit_regression_scoring_cv!(
             best_psms,
             file_paths,
             features,
             match_between_runs
             # No n_folds parameter - uses library-assigned CV folds
        )
        # Debug: Check prob values BEFORE dropVectorColumns
        @debug_l1 "Probit: BEFORE dropVectorColumns - prob range: $(minimum(best_psms.prob)) to $(maximum(best_psms.prob))"
        @debug_l1 "  DataFrame has $(nrow(best_psms)) rows, $(ncol(best_psms)) columns"
        if match_between_runs
            @debug_l1 "  MBR_prob exists: $(:MBR_prob in names(best_psms)), MBR_is_best_decoy exists: $(:MBR_is_best_decoy in names(best_psms))"
        end
        
        # Write the scored PSMs back to their original files, just like XGBoost does
        # CRITICAL: Drop vector columns before writing to avoid Arrow serialization issues
        dropVectorColumns!(best_psms)
        
        # Debug: Check prob values AFTER dropVectorColumns
        @debug_l1 "Probit: AFTER dropVectorColumns - prob range: $(minimum(best_psms.prob)) to $(maximum(best_psms.prob))"
        @debug_l1 "  DataFrame has $(nrow(best_psms)) rows, $(ncol(best_psms)) columns"
        @debug_l1 "  Columns: prob=$(:prob in names(best_psms)), decoy=$(:decoy in names(best_psms))"
        if match_between_runs
            @debug_l1 "  MBR columns after drop: MBR_prob=$(:MBR_prob in names(best_psms)), MBR_is_best_decoy=$(:MBR_is_best_decoy in names(best_psms))"
        end
        
        for (ms_file_idx, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
            # Debug: Check each group before writing
            @debug_l1 "  Writing group $ms_file_idx: $(nrow(gpsms)) rows, prob range: $(minimum(gpsms.prob)) to $(maximum(gpsms.prob))"
            
            fpath = file_paths[ms_file_idx[:ms_file_idx]]
            writeArrow(fpath, gpsms)
        end
        
        # Debug: Check what was written to files
        @debug_l1 "Probit: Verifying written files"
        for (i, fpath) in enumerate(file_paths)
            df_check = DataFrame(Arrow.Table(fpath))
            @debug_l1 "  File $i after write: $(nrow(df_check)) rows, prob range: $(minimum(df_check.prob)) to $(maximum(df_check.prob))"
            if match_between_runs && :MBR_prob in names(df_check)
                @debug_l1 "    MBR_prob range: $(minimum(df_check.MBR_prob)) to $(maximum(df_check.MBR_prob))"
            end
            @debug_l1 "    Columns in file: $(names(df_check)[1:min(10, length(names(df_check)))])..."
        end
        
        models = nothing  # Probit doesn't return models
        
        # OPTION 2: XGBoost/EvoTrees (WITH ITERATIVE REFINEMENT)
        #see src/utils/ML/percolatorSortOf.jl
        #Train EvoTrees/XGBoost model to score each precursor trace. Target-decoy descrimination
        =#
        @user_warn "TEST TEST TEST XGboost"
        models = sort_of_percolator_in_memory!(
                                best_psms, 

                                file_paths,
                                features,
                                match_between_runs;
                                max_q_value_xgboost_rescore,
                                max_q_value_xgboost_mbr_rescore,
                                min_PEP_neg_threshold_xgboost_rescore,
                                colsample_bytree = 0.8, 
                                colsample_bynode = 0.8,
                                min_child_weight = 20, 
                                gamma = 0.1,
                                subsample = 0.8, 
                                max_depth = 4,
                                eta = 0.1,
                                iter_scheme = [150, 300, 300],
                                print_importance = true);
        
        return models;
    end
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
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32
)
    if size(best_psms, 1) > 100000
        file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
        features = [ 
            :missed_cleavage,
            :Mox,
            :prec_mz,
            :sequence_length,
            :charge,
            :irt_pred,
            :irt_error,
            :irt_diff,
            :max_y_ions,
            :y_ions_sum,
            :longest_y,
            :y_count,
            :b_count,
            :isotope_count,
            :total_ions,
            :best_rank,
            :best_rank_iso,
            :topn,
            :topn_iso,
            :gof,
            :max_fitted_manhattan_distance,
            :max_fitted_spectral_contrast,
            :max_matched_residual,
            :max_unmatched_residual,
            :max_gof,
            :fitted_spectral_contrast,
            :spectral_contrast,
            :max_matched_ratio,
            :err_norm,
            :poisson,
            :weight,
            :log2_intensity_explained,
            :tic,
            :num_scans,
            :smoothness,
            :rt_diff,
            :ms1_irt_diff,
            :weight_ms1,
            :gof_ms1,
            :max_matched_residual_ms1,
            :max_unmatched_residual_ms1,
            :fitted_spectral_contrast_ms1,
            :error_ms1,
            :m0_error_ms1,
            :n_iso_ms1,
            :big_iso_ms1,
            :ms1_features_missing,
            :percent_theoretical_ignored,
            :scribe,
            :max_scribe,
            :target
        ];
        features = [f for f in features if hasproperty(best_psms, f)];
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

        models = sort_of_percolator_out_of_memory!(
                                best_psms, 
                                file_paths,
                                features,
                                match_between_runs;
                                max_q_value_xgboost_rescore,
                                max_q_value_xgboost_mbr_rescore,
                                min_PEP_neg_threshold_xgboost_rescore,
                                colsample_bytree = 0.5, 
                                colsample_bynode = 0.5,
                                min_child_weight = 5, 
                                gamma = 1,
                                subsample = 0.25, 
                                max_depth = 10,
                                eta = 0.05, 
                                iter_scheme = [100, 200, 200],
                                print_importance = false);
        return models;#best_psms
    else
        @user_warn "Less than 100,000 psms. Training with simplified target-decoy discrimination model..."
        file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
        features = [ 
            :missed_cleavage,
            :Mox,
            :sequence_length,
            :charge,
            :irt_error,
            :irt_diff,
            :y_count,
            :max_fitted_manhattan_distance,
            :max_matched_residual,
            :max_unmatched_residual,
            :max_gof,
            :err_norm,
            :weight,
            :log2_intensity_explained,
        ];
        best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
        best_psms[!,:decoy] = best_psms[!,:target].==false;
        #see src/utils/ML/percolatorSortOf.jl
        #Train EvoTrees/XGBoost model to score each precursor trace. Target-decoy descrimination
        models = sort_of_percolator_out_of_memory!(
                                best_psms, 
                                file_paths,
                                features,
                                match_between_runs;
                                max_q_value_xgboost_rescore,
                                max_q_value_xgboost_mbr_rescore,
                                min_PEP_neg_threshold_xgboost_rescore,
                                colsample_bytree = 1.0, 
                                colsample_bynode = 1.0,
                                min_child_weight = 100, 
                                gamma = 0,
                                subsample = 1.0, 
                                max_depth = 3,
                                eta = 0.01, 
                                iter_scheme = [200],
                                print_importance = false);
        return models;#best_psms
    end
end