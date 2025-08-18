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
    max_psms_in_memory::Int64
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
    
    # Clean up
    best_psms = nothing
    GC.gc()
    
    return models
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
    match_between_runs::Bool;
    n_folds::Int64 = 3
)
    # Step 1: Add CV fold column based on file index (same as XGBoost)
    psms[!, :cv_fold] = zeros(UInt8, size(psms, 1))
    
    # Always use :ms_file_idx column - it's created by SecondPassSearch
    ms_file_idx_col = :ms_file_idx
    
    # Check if the column exists
    if !(ms_file_idx_col in propertynames(psms))
        @error "Column :ms_file_idx not found in PSMs. Available columns: $(propertynames(psms))"
        error("Required column :ms_file_idx is missing from PSMs DataFrame")
    end
    
    #Will need to determine this on the protien/accession_numbers level from the spectral library
    #This is fine for testing atm. 
    for fold_idx in 1:n_folds
        fold_mask = (psms[!, ms_file_idx_col] .% n_folds) .== (fold_idx - 1)
        psms[fold_mask, :cv_fold] .= UInt8(fold_idx)
    end
    
    # Step 2: Initialize score columns
    psms[!, :prob] = zeros(Float32, size(psms, 1))
    
    # Step 3: Use the exact features passed in (already filtered by the caller)
    # Add intercept column if not present
    if !(:intercept in propertynames(psms))
        psms[!, :intercept] = ones(Float32, size(psms, 1))
    end
    
    # Filter features to those that exist in the DataFrame
    available_features = filter(col -> col in propertynames(psms), features)
    push!(available_features, :intercept)
    
    @info "Probit regression using $(length(available_features)-1) features (plus intercept)"
    
    # Step 4: Train probit model per CV fold
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = Iterators.partition(1:size(psms, 1), chunk_size)
    
    for fold_idx in 1:n_folds
        # Get train/test masks
        test_mask = psms[!, :cv_fold] .== fold_idx
        train_mask = .!test_mask
        
        # Check for sufficient training data
        train_targets = train_mask .& psms[!, :target]
        train_decoys = train_mask .& .!psms[!, :target]
        
        n_train_targets = sum(train_targets)
        n_train_decoys = sum(train_decoys)
        
        if n_train_targets < 100 || n_train_decoys < 100
            @warn "Insufficient training data for fold $fold_idx (targets: $n_train_targets, decoys: $n_train_decoys)"
            psms[test_mask, :prob] .= 0.5f0
            continue
        end
        
        @info "Training probit model for fold $fold_idx with $n_train_targets targets and $n_train_decoys decoys"
        
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
        
        psms[test_mask, :prob] = test_probs
    end
    
    # Match XGBoost behavior - just fill in prob column, no q-values or best_psm selection
    
    # Ensure prob column has proper type and no NaN/Inf values
    # Check for any NaN or Inf values and replace them
    n_nan = sum(isnan.(psms.prob))
    n_inf = sum(isinf.(psms.prob))
    if n_nan > 0 || n_inf > 0
        @warn "Found $n_nan NaN and $n_inf Inf values in probabilities, replacing with 0.5"
        psms.prob[isnan.(psms.prob) .| isinf.(psms.prob)] .= 0.5f0
    end
    
    # Ensure values are in valid range [0,1] but avoid exact 0 and 1
    # XGBoost/EvoTrees never produces exact 0 or 1, so match that behavior
    # This is important for downstream log operations
    psms.prob = clamp.(psms.prob, eps(Float32), 1.0f0 - eps(Float32))
    
    # Clean up columns added during probit regression that shouldn't be persisted
    if :intercept in propertynames(psms)
        select!(psms, Not(:intercept))
    end
    if :cv_fold in propertynames(psms)
        select!(psms, Not(:cv_fold))
    end
    
    Arrow.write("/Users/nathanwamsley/Desktop/test_arrow_psms.arrow", psms)
    @info "Probit regression scoring complete. Probabilities assigned to $(size(psms, 1)) PSMs (range: $(minimum(psms.prob)) to $(maximum(psms.prob)))"
    
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
    if size(best_psms, 1) > 100_000
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
        @warn "Less than 100,000 psms: $(size(best_psms, 1)). Training with simplified target-decoy discrimination model..."
        if size(best_psms, 1) > 1000
            @warn "Less than 100,000 psms. Training with simplified target-decoy discrimination model..."
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
        else
             @warn "Less than 1,000 psms. Training with super simplified target-decoy discrimination model..."
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
        
        probit_regression_scoring_cv!(
             best_psms,
             file_paths,
             features,
             match_between_runs;
             n_folds = 3
        )
        
        # Write the scored PSMs back to their original files, just like XGBoost does
        # CRITICAL: Drop vector columns before writing to avoid Arrow serialization issues
        dropVectorColumns!(best_psms)
        
        for (ms_file_idx, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
            fpath = file_paths[ms_file_idx[:ms_file_idx]]
            writeArrow(fpath, gpsms)
        end
        
        # Debug: Check what was written to files
        @info "Probit: Written PSMs to files" n_files=length(file_paths)
        for (i, fpath) in enumerate(file_paths)
            df_check = DataFrame(Arrow.Table(fpath))
            @info "  File $i: $(nrow(df_check)) rows, prob range: $(minimum(df_check.prob)) to $(maximum(df_check.prob))"
        end
        
        models = nothing  # Probit doesn't return models
        
        # OPTION 2: XGBoost/EvoTrees (WITH ITERATIVE REFINEMENT)
        #see src/utils/ML/percolatorSortOf.jl
        #Train EvoTrees/XGBoost model to score each precursor trace. Target-decoy descrimination
        #=
        @warn "TEST TEST TEST XGboost"
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
                                iter_scheme = [300],
                                print_importance = true);
        =#
        return models
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
        @warn "Less than 100,000 psms. Training with simplified target-decoy discrimination model..."
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