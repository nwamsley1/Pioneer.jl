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
                                  max_q_value_lightgbm_rescore::Float32,
                                  max_q_value_mbr_itr::Float32,
                                  min_PEP_neg_threshold_itr::Float32,
                                  max_psms_in_memory::Int64)

Main entry point for PSM scoring with automatic model selection based on dataset size.

# Three-Case Logic
1. PSMs ≥ max_psms_in_memory: Out-of-memory processing with default LightGBM
2. PSMs < max_psms_in_memory AND ≥ 100K: In-memory with default/advanced LightGBM
3. PSMs < 100K: In-memory with automatic model comparison

# Arguments
- `second_pass_folder`: Folder containing second pass PSM files
- `file_paths`: Vector of PSM file paths
- `precursors`: Library precursors
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_lightgbm_rescore`: Max q-value for LightGBM rescoring
- `max_q_value_mbr_itr`: Max q-value for MBR transfers retained during iterative training (ITR)
- `min_PEP_neg_threshold_itr`: Min PEP threshold for relabeling weak targets as negatives during ITR
- `max_psms_in_memory`: Maximum PSMs to keep in memory
- `n_quantile_bins`: Number of quantile bins for feature discretization (1-65535)
- `q_value_threshold`: Q-value threshold for model comparison (default: 0.01)
- `ms1_scoring`: Whether to include MS1 scoring features (default: true)

# Returns
- Trained LightGBM models or nothing for probit regression
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    max_psms_in_memory::Int64,
    max_psm_memory_mb::Int64,  # User-specified memory limit for PSMs (MB)
    n_quantile_bins::Int64,
    q_value_threshold::Float32 = 0.01f0,  # Default to 1% if not specified
    ms1_scoring::Bool = true
)
    # Step 1: Count PSMs and determine processing approach
    psms_count = get_psms_count(file_paths)

    # Calculate memory-based threshold
    # Estimated bytes per PSM including all features and overhead
    # Note: max_psm_memory_mb is validated to be >= 50 MB at parameter loading
    psm_size_bytes = 500
    memory_based_threshold = floor(Int64, max_psm_memory_mb * 1e6 / psm_size_bytes)
    effective_threshold = min(max_psms_in_memory, memory_based_threshold)

    if psms_count >= effective_threshold
        # Case 1: Out-of-memory processing with default LightGBM
        @user_info "Using OUT-OF-MEMORY mode: $psms_count PSMs exceeds $effective_threshold threshold"
        @debug_l1 "\n[OOM] Out-of-memory percolator activated"
        @debug_l1 "\n[OOM]   PSM count: $psms_count"
        @debug_l1 "\n[OOM]   Memory threshold: $max_psm_memory_mb MB → $memory_based_threshold PSMs"
        @debug_l1 "\n[OOM]   Count threshold: $max_psms_in_memory PSMs"
        @debug_l1 "\n[OOM]   Effective threshold: $effective_threshold PSMs"
        @debug_l1 "\n[OOM]   Files: $(length(file_paths))"

        # Phase 0: Assign pair_ids across all files
        @debug_l1 "\n[OOM] === PHASE 0: Pair ID Assignment ==="
        phase0_timing = @timed assign_pair_ids_oom!(file_paths)
        @debug_l1 "\n[OOM] Phase 0 complete: $(round(phase0_timing.time, digits=2))s"

        # Phase 1: Sample complete pair groups for training (up to effective_threshold PSMs)
        @debug_l1 "\n[OOM] === PHASE 1: Sampling PSMs for Training ==="
        phase1_timing = @timed begin
            best_psms = sample_complete_pairs_for_training(
                file_paths,
                effective_threshold
            )
        end
        @debug_l1 "\n[OOM] Phase 1 complete: $(round(phase1_timing.time, digits=2))s, $(nrow(best_psms)) PSMs sampled"

        # Phase 2: Add quantile-binned features and save bin edges for OOM prediction
        @debug_l1 "\n[OOM] === PHASE 2: Feature Engineering ==="
        phase2_timing = @timed begin
            features_to_bin = [:prec_mz, :irt_pred, :weight, :tic]
            # Compute bin edges BEFORE adding binned features (so we capture edges from training data)
            bin_edges = compute_quantile_bin_edges(best_psms, features_to_bin, n_quantile_bins)
            add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)
        end
        @debug_l1 "\n[OOM] Phase 2 complete: $(round(phase2_timing.time, digits=2))s"
        @debug_l1 "\n[OOM]   Computed bin edges for $(length(bin_edges)) features"

        # Use AdvancedLightGBM config for OOM
        model_config = create_default_advanced_lightgbm_config(ms1_scoring)

        # Phase 3: Train models and apply to ALL files (this is the expensive part)
        @debug_l1 "\n[OOM] === PHASE 3: Training & Scoring ALL Files ==="
        @debug_l1 "\n[OOM]   Training on $(nrow(best_psms)) sampled PSMs"
        @debug_l1 "\n[OOM]   Will apply models to $(length(file_paths)) files with $psms_count total PSMs"
        phase3_timing = @timed begin
            models = score_precursor_isotope_traces_out_of_memory!(
                best_psms,
                file_paths,
                precursors,
                model_config,
                match_between_runs,
                max_q_value_lightgbm_rescore,
                max_q_value_mbr_itr,
                min_PEP_neg_threshold_itr;
                bin_edges = bin_edges
            )
        end
        @debug_l1 "\n[OOM] Phase 3 complete: $(round(phase3_timing.time, digits=2))s"
        @debug_l1 "\n[OOM] === OOM Pipeline Complete ==="
        total_time = phase0_timing.time + phase1_timing.time + phase2_timing.time + phase3_timing.time
        @debug_l1 "\n[OOM] Total OOM time: $(round(total_time, digits=2))s"
    else
        # In-memory processing - load PSMs first
        best_psms = load_psms_for_lightgbm(second_pass_folder)

        # Add quantile-binned features before training
        features_to_bin = [:prec_mz, :irt_pred, :weight, :tic]
        add_quantile_binned_features!(best_psms, features_to_bin, n_quantile_bins)

        if psms_count >= MAX_FOR_MODEL_SELECTION  # 100K
            # Case 2: In-memory with default/advanced LightGBM (no comparison)
            @user_info "Using in-memory advanced LightGBM for $psms_count PSMs (< $max_psms_in_memory but ≥ 100K)"
            model_config = create_default_advanced_lightgbm_config(ms1_scoring)
        else
            # Case 3: In-memory with automatic model comparison (<100K)
            model_config = select_psm_scoring_model(
                best_psms, file_paths, precursors, match_between_runs,
                max_q_value_lightgbm_rescore, max_q_value_mbr_itr,
                min_PEP_neg_threshold_itr, q_value_threshold, ms1_scoring
            )
        end
        
        # Execute selected model using unified function
        @user_info "Training final model: $(model_config.name)"
        models = score_precursor_isotope_traces_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_lightgbm_rescore,
            max_q_value_mbr_itr, min_PEP_neg_threshold_itr
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
- ≥100K PSMs: Returns default/advanced LightGBM configuration (no comparison)
- <100K PSMs: Trains each model and selects based on training performance

# Returns
- ModelConfig object specifying the selected model and its hyperparameters
"""
function select_psm_scoring_model(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    q_value_threshold::Float32,
    ms1_scoring::Bool = true
)
    psms_count = size(best_psms, 1)
    
    if psms_count >= MAX_FOR_MODEL_SELECTION
        @user_info "Using default advanced LightGBM for $psms_count PSMs (≥ 100K)"
        return create_default_advanced_lightgbm_config(ms1_scoring)
    else
        # Get model configurations from model_config.jl
        model_configs = create_model_configurations(ms1_scoring)
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
                    match_between_runs, max_q_value_lightgbm_rescore,
                    max_q_value_mbr_itr, min_PEP_neg_threshold_itr,
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
                throw(e)
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
        
        # If no model succeeded, fall back to SimpleLightGBM
        if best_model_config === nothing
            @user_warn "All models failed, defaulting to SimpleLightGBM"
            best_model_config = model_configs[1]  # SimpleLightGBM is first
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
- `scored_psms`: DataFrame containing scored PSMs with :trace_prob and :target columns
- `qvalue_threshold`: Q-value threshold for counting (typically 0.01)

# Returns
- Total count of target PSMs with q_value ≤ threshold
"""
function count_passing_targets(scored_psms::DataFrame, qvalue_threshold::Float32)
    # Use existing q-value calculation logic
    if :trace_prob in propertynames(scored_psms) && :target in propertynames(scored_psms)
        # Calculate q-values from probabilities
        qvals = Vector{Float32}(undef, nrow(scored_psms))
        get_qvalues!(scored_psms.trace_prob, scored_psms.target, qvals)

        # Count targets passing threshold
        return sum((scored_psms.target .== true) .& (qvals .<= qvalue_threshold))
    else
        error("DataFrame must contain :trace_prob and :target columns")
    end
end

"""
    create_default_advanced_lightgbm_config(ms1_scoring::Bool = true) -> ModelConfig

Creates the default advanced LightGBM configuration for large datasets.
"""
function create_default_advanced_lightgbm_config(ms1_scoring::Bool = true)
    features = copy(ADVANCED_FEATURE_SET)
    apply_ms1_filtering!(features, ms1_scoring)

    return ModelConfig(
        "AdvancedLightGBM",
        :lightgbm,
        features,
        Dict(
            :feature_fraction => 0.5,
            :min_data_in_leaf => 500,
            :min_gain_to_split => 0.5,
            :bagging_fraction => 0.25,
            :max_depth => 10,
            :num_leaves => 63,
            :learning_rate => 0.05,
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
- SimpleLightGBM: Full feature set, standard hyperparameters
- AdvancedLightGBM: Full feature set, advanced hyperparameters
- ProbitRegression: Linear probit model with CV folds from library
- SuperSimplified: Minimal 5-feature LightGBM model
"""
function score_precursor_isotope_traces_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    show_progress::Bool = true
)
    
    if model_config.model_type == :lightgbm
        return train_lightgbm_model_in_memory(
            best_psms, file_paths, precursors, model_config,
            match_between_runs, max_q_value_lightgbm_rescore,
            max_q_value_mbr_itr, min_PEP_neg_threshold_itr,
            show_progress
        )
    elseif model_config.model_type == :probit
        return train_probit_model_in_memory(
            best_psms, file_paths, precursors, model_config, match_between_runs,
            min_PEP_neg_threshold_itr
        )
    else
        error("Unsupported model type: $(model_config.model_type)")
    end
end

"""
    train_lightgbm_model_in_memory(...) -> Models

Trains LightGBM model using configuration from ModelConfig.
"""
function train_lightgbm_model_in_memory(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
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
            #:MBR_num_runs,
            :MBR_max_pair_prob,
            :MBR_log2_weight_ratio,
            :MBR_log2_explained_ratio,
            :MBR_rv_coefficient,
            #:MBR_best_irt_diff,
            :MBR_is_missing
        ])
    end

    # Diagnostic: Report which quantile-binned features are being used
    qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
    if !isempty(qbin_features)
        for qbin_feat in qbin_features
            n_unique = length(unique(skipmissing(best_psms[!, qbin_feat])))
        end
    end

    hp = model_config.hyperparams
    return sort_of_percolator_in_memory!(
        best_psms, features, match_between_runs;
        max_q_value_lightgbm_rescore, max_q_value_mbr_itr,
        min_PEP_neg_threshold_itr,
        feature_fraction = get(hp, :feature_fraction, 0.5),
        min_data_in_leaf = get(hp, :min_data_in_leaf, 500),
        min_gain_to_split = get(hp, :min_gain_to_split, 0.5),
        bagging_fraction = get(hp, :bagging_fraction, 0.25),
        max_depth = get(hp, :max_depth, 10),
        num_leaves = get(hp, :num_leaves, 63),
        learning_rate = get(hp, :learning_rate, 0.05),
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
    min_PEP_neg_threshold_itr::Float32
)
    # Add required columns
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!,:decoy] = best_psms[!,:target].==false
    
    # Get features from config
    features = [f for f in model_config.features if hasproperty(best_psms, f)]

    # Diagnostic: Report which quantile-binned features are being used
    qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
    if !isempty(qbin_features)
        for qbin_feat in qbin_features
            n_unique = length(unique(skipmissing(best_psms[!, qbin_feat])))
            col_type = eltype(best_psms[!, qbin_feat])
        end
    end

    probit_regression_scoring_cv!(
        best_psms,
        file_paths,
        features,
        match_between_runs;
        neg_mining_pep_threshold = min_PEP_neg_threshold_itr
    )
    
    # File writing removed - will be done at higher level
    
    return nothing  # Probit doesn't return models
end

"""
     get_psms_count(quant_psms_folder::String)::Integer

Sample PSMs from multiple files for LightGBM model training.

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
    sample_psms_for_lightgbm(quant_psms_folder::String, psms_count::Integer max_psms::Integer) -> DataFrame

Sample PSMs from multiple files for LightGBM model training.

# Arguments
- `quant_psms_folder`: Folder containing PSM Arrow files
- `psms_count`: number of psms across all the arrow files
- `max_psms`: Maximum number of PSMs to sample for training

# Process
1. Proportionally samples from each file
2. Combines samples into single DataFrame
"""
function sample_psms_for_lightgbm(quant_psms_folder::String, psms_count::Integer, max_psms::Integer)

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

Loads all PSMs from multiple files for LightGBM model training.
"""
function load_psms_for_lightgbm(quant_psms_folder::String)
    file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]
    return DataFrame(Tables.columntable(Arrow.Table(file_paths)))
end

"""
    score_precursor_isotope_traces_in_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::LibraryPrecursors) -> Dictionary{UInt8, LightGBMModel}

Train LightGBM models for PSM scoring. All psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained LightGBM models or simplified model if insufficient PSMs.
"""

"""
    probit_regression_scoring_cv!(psms::DataFrame,
                                  file_paths::Vector{String},
                                  features::Vector{Symbol},
                                  match_between_runs::Bool;
                                  n_folds::Int64 = 3,
                                  neg_mining_pep_threshold::Float32 = 0.90f0)

Alternative PSM scoring using probit regression with cross-validation.

This is a simpler alternative to LightGBM for small datasets (<100k PSMs).
Uses linear probit model with cross-validation, similar to FirstPassSearch but with CV folds.
No iterative refinement or max_prob updates - single pass training only.

# Arguments
- `psms`: DataFrame containing PSMs to score
- `file_paths`: Vector of file paths for CV fold assignment
- `features`: Feature columns to use for scoring (same as LightGBM)
- `match_between_runs`: Whether MBR was performed
- `n_folds`: Number of cross-validation folds (default: 3)

# Modifies
- Adds columns to psms: `:trace_prob`, `:q_value`, `:best_psm`, `:cv_fold`
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
    
    # Step 2: Initialize probability array (like LightGBM does)
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

    # Check for potential zero-variance issues
    problematic_features = Symbol[]
    for feature in available_features
        if feature != :intercept  # Skip intercept column
            col_data = psms[!, feature]
            if length(unique(col_data)) <= 1
                push!(problematic_features, feature)
            end
        end
    end

    if !isempty(problematic_features)
        @user_warn "Found $(length(problematic_features)) zero-variance features: $(join(string.(problematic_features), ", "))"

        # Remove zero-variance features from available_features
        filter!(feature -> !(feature in problematic_features), available_features)
        @user_info "After filtering: using $(length(available_features)) features"

        if length(available_features) <= 1  # Only intercept left
            @user_warn "No valid features remaining for probit regression after zero-variance filtering"
            return  # Exit early to avoid singular matrix
        end
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
            psms[test_mask, :trace_prob] .= 0.5f0
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
    
    # Step 5: Assign probabilities to DataFrame (like LightGBM does at line 105 of percolatorSortOf.jl)
    psms[!, :trace_prob] = prob_estimates

    # If MBR is enabled, create MBR columns (probit doesn't do separate MBR scoring)
    # This is needed for compatibility with apply_mbr_filter! in ScoringSearch
    if match_between_runs
        psms[!, :MBR_boosted_trace_prob] = copy(prob_estimates)
        # MBR_is_best_decoy is required by apply_mbr_filter!
        # For probit, we don't have MBR transfer info, so set conservatively
        psms[!, :MBR_is_best_decoy] = fill(missing, size(psms, 1))  # missing means not an MBR transfer
    end

    # Ensure trace_prob column has proper type and no NaN/Inf values
    # Check for any NaN or Inf values and replace them
    n_nan = sum(isnan.(psms.trace_prob))
    n_inf = sum(isinf.(psms.trace_prob))
    if n_nan > 0 || n_inf > 0
        @user_warn "Found $n_nan NaN and $n_inf Inf values in probabilities, replacing with 0.5"
        psms.trace_prob[isnan.(psms.trace_prob) .| isinf.(psms.trace_prob)] .= 0.5f0
    end

    # Ensure values are in valid range but avoid values too close to 0 or 1
    # The downstream aggregation formula breaks with probabilities > 0.999999
    # Use more conservative clamping to ensure we stay below the threshold
    # Force Float32 to ensure proper clamping
    psms.trace_prob = Float32.(clamp.(psms.trace_prob, Float32(1e-6), Float32(0.9999)))
    
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
    assign_pair_ids_oom!(file_paths::Vector{String})

Assign pair_ids to ALL PSMs across ALL files before sampling.
This ensures pair groups are never split between training sample and OOM data.

# Process
1. Pass 1: Collect all unique (precursor_idx, irt_pred, cv_fold, isotopes_captured, target) from all files
2. Pass 2: Assign pair_ids using same algorithm as in-memory (via assignPairIds!)
3. Pass 3: Write pair_ids back to all files

# Returns
- `pair_lookup`: Dictionary mapping (precursor_idx, cv_fold, isotopes_captured) -> pair_id
"""
function assign_pair_ids_oom!(file_paths::Vector{String})
    @debug_l1 "\n[OOM] === Pass 1: Collecting precursor info from $(length(file_paths)) files ==="

    # Pass 1: Collect all unique (precursor_idx, irt_pred, cv_fold, isotopes_captured, target)
    pass1_timing = @timed begin
        precursor_info = DataFrame(
            precursor_idx = UInt32[],
            irt_pred = Float32[],
            cv_fold = UInt8[],
            isotopes_captured = Tuple{Int8,Int8}[],
            target = Bool[]
        )

        for (i, file_path) in enumerate(file_paths)
            df = DataFrame(Arrow.Table(file_path))
            # Get unique precursors from this file
            unique_precs = unique(df, [:precursor_idx, :cv_fold, :isotopes_captured])
            append!(precursor_info, select(unique_precs,
                :precursor_idx, :irt_pred, :cv_fold, :isotopes_captured, :target))
            if i % 10 == 0 || i == length(file_paths)
                @debug_l2 "\n[OOM]   Read $i/$(length(file_paths)) files, $(nrow(precursor_info)) precursors collected..."
            end
        end

        # Deduplicate across files
        unique!(precursor_info, [:precursor_idx, :cv_fold, :isotopes_captured])
    end
    @debug_l1 "\n[OOM]   Pass 1 complete: $(round(pass1_timing.time, digits=2))s"
    @debug_l1 "\n[OOM]   Found $(nrow(precursor_info)) unique precursor-isotope combinations"

    # Pass 2: Assign pair_ids using same algorithm as in-memory
    @debug_l1 "\n[OOM] === Pass 2: Assigning pair_ids ==="
    pass2_timing = @timed begin
        # Sort for deterministic ordering regardless of file processing order
        sort!(precursor_info, [:cv_fold, :isotopes_captured, :irt_pred, :precursor_idx])

        # Add irt_bin_idx using the existing function
        precursor_info[!, :irt_bin_idx] = getIrtBins(precursor_info.irt_pred)

        # Add decoy column (required by assignPairIds! which calls assign_pair_ids)
        precursor_info[!, :decoy] = .!precursor_info.target

        # Assign pair_ids within groups using existing assignPairIds! function
        last_pair_id = zero(UInt32)
        precursor_info[!, :pair_id] = zeros(UInt32, nrow(precursor_info))

        for group in groupby(precursor_info, [:irt_bin_idx, :cv_fold, :isotopes_captured])
            last_pair_id = assignPairIds!(group, last_pair_id)
        end

        # Create lookup: (precursor_idx, cv_fold, isotopes) -> pair_id
        # Use column iteration instead of eachrow() for performance
        pair_lookup = Dict{Tuple{UInt32, UInt8, Tuple{Int8,Int8}}, UInt32}()
        for (pid, cv, iso, pair) in zip(precursor_info.precursor_idx,
                                         precursor_info.cv_fold,
                                         precursor_info.isotopes_captured,
                                         precursor_info.pair_id)
            key = (pid, cv, iso)
            pair_lookup[key] = pair
        end
    end
    @debug_l1 "\n[OOM]   Pass 2 complete: $(round(pass2_timing.time, digits=2))s"
    @debug_l1 "\n[OOM]   Assigned $(length(pair_lookup)) pair_ids"

    # Pass 3: Write pair_ids back to all files
    @debug_l1 "\n[OOM] === Pass 3: Writing pair_ids to $(length(file_paths)) files ==="
    pass3_timing = @timed begin
        for (file_idx, file_path) in enumerate(file_paths)
            @debug_l2 "\n[OOM]   File $file_idx/$(length(file_paths)): $file_path"

            # Create file reference for streaming column operation
            ref = PSMFileReference(file_path)

            # Define compute function using pair_lookup closure
            read_time = @elapsed begin
                compute_pair_id = function(df_batch)
                    n = nrow(df_batch)
                    result = Vector{UInt32}(undef, n)
                    pids = df_batch.precursor_idx
                    cvs = df_batch.cv_fold
                    isos = df_batch.isotopes_captured
                    for j in 1:n
                        key = (pids[j], cvs[j], isos[j])
                        result[j] = get(pair_lookup, key, zero(UInt32))
                    end
                    return result
                end
            end
            @debug_l2 "\n[OOM]     Function setup: $(round(read_time * 1000, digits=1))ms"

            # Use existing streaming infrastructure to add column
            write_time = @elapsed add_column_to_file!(ref, :pair_id, compute_pair_id)
            @debug_l2 "\n[OOM]     add_column_to_file!: $(round(write_time, digits=2))s"
        end
    end
    @debug_l1 "\n[OOM]   Pass 3 complete: $(round(pass3_timing.time, digits=2))s"

    return pair_lookup
end

"""
    sample_complete_pairs_for_training(file_paths, max_psms)

Sample COMPLETE pair groups for training. Never splits a pair between
training sample and OOM data.

# Process
1. Collect all unique pair_ids and count their PSMs across all files
2. Randomly select pair_ids until we reach target PSM count
3. Load PSMs belonging to selected pairs from all files

# Arguments
- `file_paths`: Vector of PSM file paths
- `max_psms`: Maximum number of PSMs to sample for training

# Returns
- DataFrame containing sampled PSMs (complete pairs only)
"""
function sample_complete_pairs_for_training(
    file_paths::Vector{String},
    max_psms::Int64
)
    @debug_l1 "\n[OOM] Sampling complete pairs for training (target: $max_psms PSMs)..."
    @debug_l1 "\n[OOM]   Reading pair counts from $(length(file_paths)) files..."

    # Collect all unique pair_ids and their PSM counts
    count_timing = @timed begin
        pair_counts = Dict{UInt32, Int}()
        for (i, file_path) in enumerate(file_paths)
            df = DataFrame(Arrow.Table(file_path))
            # Iterate directly over column vector (fast)
            for pair_id in df.pair_id
                pair_counts[pair_id] = get(pair_counts, pair_id, 0) + 1
            end
            if i % 10 == 0 || i == length(file_paths)
                @debug_l2 "\n[OOM]     Processed $i/$(length(file_paths)) files..."
            end
        end
    end
    @debug_l1 "\n[OOM]   Pair counting: $(round(count_timing.time, digits=2))s"

    all_pair_ids = collect(keys(pair_counts))
    @debug_l1 "\n[OOM]   Found $(length(all_pair_ids)) unique pair_ids"

    # Shuffle and select pairs until we reach target PSM count
    Random.seed!(1776)
    shuffled_pairs = shuffle(all_pair_ids)

    selected_pairs = Set{UInt32}()
    total_psms = 0
    for pair_id in shuffled_pairs
        if total_psms + pair_counts[pair_id] <= max_psms
            push!(selected_pairs, pair_id)
            total_psms += pair_counts[pair_id]
        end
        if total_psms >= max_psms * 0.9  # Stop at 90% to avoid overshoot
            break
        end
    end

    @debug_l1 "\n[OOM]   Selected $(length(selected_pairs)) pairs with $total_psms PSMs"

    # Load PSMs belonging to selected pairs
    @debug_l1 "\n[OOM]   Loading selected PSMs from files..."
    load_timing = @timed begin
        sampled_dfs = DataFrame[]
        for (i, file_path) in enumerate(file_paths)
            df = DataFrame(Arrow.Table(file_path))
            mask = in.(df.pair_id, Ref(selected_pairs))
            if any(mask)
                push!(sampled_dfs, df[mask, :])
            end
            if i % 10 == 0 || i == length(file_paths)
                @debug_l2 "\n[OOM]     Loaded $i/$(length(file_paths)) files..."
            end
        end
    end
    @debug_l1 "\n[OOM]   File loading: $(round(load_timing.time, digits=2))s"

    concat_timing = @timed begin
        sampled_psms = vcat(sampled_dfs...)
    end
    @debug_l1 "\n[OOM]   DataFrame concatenation: $(round(concat_timing.time, digits=2))s"
    @debug_l1 "\n[OOM]   Loaded $(nrow(sampled_psms)) PSMs for training"

    return sampled_psms
end

"""
    score_precursor_isotope_traces_out_of_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::LibraryPrecursors) -> Dictionary{UInt8, LightGBMModel}

Train LightGBM models for PSM scoring using out-of-memory processing.
Only a subset of psms (complete pair groups) are kept in memory for training.

# Arguments
- `best_psms`: Sampled PSMs (complete pair groups) for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information
- `model_config`: Model configuration with features and hyperparameters
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_lightgbm_rescore`: Max q-value for LightGBM rescoring
- `max_q_value_mbr_itr`: Max q-value for MBR transfers during iterative training
- `min_PEP_neg_threshold_itr`: Min PEP threshold for relabeling weak targets

# Returns
Trained LightGBM models dictionary.
"""
function score_precursor_isotope_traces_out_of_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32;
    bin_edges::Dict{Symbol, Vector{Float64}} = Dict{Symbol, Vector{Float64}}()
)
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    # Features from model_config; do not include :target
    features = [f for f in model_config.features if hasproperty(best_psms, f)]
    if match_between_runs
        append!(features, [
            :MBR_rv_coefficient,
            #:MBR_best_irt_diff,
            #:MBR_num_runs,
            :MBR_max_pair_prob,
            :MBR_log2_weight_ratio,
            :MBR_log2_explained_ratio,
            :MBR_is_missing
        ])
    end

    # Diagnostic: Report which quantile-binned features are being used
    qbin_features = filter(f -> endswith(string(f), "_qbin"), features)
    if !isempty(qbin_features)
        @debug_l2 "\n[OOM] LightGBM using $(length(qbin_features)) quantile-binned features"
        for qbin_feat in qbin_features
            n_unique = length(unique(skipmissing(best_psms[!, qbin_feat])))
            @debug_l2 "\n[OOM]   $qbin_feat: $n_unique unique values in training sample"
        end
    end

    best_psms[!, :accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!, :precursor_idx]]
    best_psms[!, :q_value] = zeros(Float32, size(best_psms, 1))
    best_psms[!, :decoy] = best_psms[!, :target] .== false

    # Hyperparameters from model_config
    hp = model_config.hyperparams
    models = sort_of_percolator_out_of_memory!(
        best_psms,
        file_paths,
        features,
        match_between_runs;
        max_q_value_lightgbm_rescore,
        max_q_value_mbr_itr,
        min_PEP_neg_threshold_itr,
        feature_fraction = get(hp, :feature_fraction, 0.5),
        min_data_in_leaf = get(hp, :min_data_in_leaf, 500),
        min_gain_to_split = get(hp, :min_gain_to_split, 0.5),
        bagging_fraction = get(hp, :bagging_fraction, 0.25),
        max_depth = get(hp, :max_depth, 10),
        num_leaves = get(hp, :num_leaves, 63),
        learning_rate = get(hp, :learning_rate, 0.05),
        iter_scheme = get(hp, :iter_scheme, [100, 200, 200]),
        print_importance = false,
        bin_edges = bin_edges
    )
    return models
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
