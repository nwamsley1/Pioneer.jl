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

#==========================================================
OOM (Out-of-Memory) Scoring Infrastructure
==========================================================#

# OOM feature set - excludes _qbin features since global quantile computation
# requires loading all values into memory. For OOM processing, we use raw features.
const OOM_FEATURE_SET = [
    :missed_cleavage,
    :Mox,
    :prec_mz,          # Use raw instead of :prec_mz_qbin
    :sequence_length,
    :charge,
    :irt_pred,         # Use raw instead of :irt_pred_qbin
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
    :weight,           # Use raw instead of :weight_qbin
    :log2_intensity_explained,
    :tic,              # Use raw instead of :tic_qbin
    :num_scans,
    :smoothness,
    :ms1_ms2_rt_diff,
    :gof_ms1,
    :max_matched_residual_ms1,
    :max_unmatched_residual_ms1,
    :fitted_spectral_contrast_ms1,
    :error_ms1,
    :m0_error_ms1,
    :n_iso_ms1,
    :big_iso_ms1,
    :rt_max_intensity_ms1,
    :rt_diff_max_intensity_ms1,
    :ms1_features_missing,
    :percent_theoretical_ignored,
    :scribe,
    :max_scribe
]

"""
    sizeof_arrow_column_type(col_type) -> Int

Returns the estimated size in bytes per element for an Arrow column type.
Used for memory estimation without loading data.
"""
function sizeof_arrow_column_type(col_type)
    # Ensure we're working with a Type, not a value
    # Some Arrow column eltypes can return unexpected results
    try
        # Handle common primitive types using pattern matching
        if col_type === Float32 || col_type === Int32 || col_type === UInt32
            return 4
        elseif col_type === Float64 || col_type === Int64 || col_type === UInt64
            return 8
        elseif col_type === Float16 || col_type === Int16 || col_type === UInt16
            return 2
        elseif col_type === Int8 || col_type === UInt8 || col_type === Bool
            return 1
        elseif col_type isa DataType
            # Try subtype checks for more complex types
            if col_type <: AbstractVector
                # For vector columns, estimate based on typical length
                # PSM vectors (weights, irts) are typically ~10-20 elements of Float32
                return 60  # Conservative estimate: 15 elements × 4 bytes
            elseif col_type <: AbstractString
                return 20  # Estimate average string length
            elseif col_type <: Tuple
                # For tuples like Tuple{Int8, Int8}, sum up component sizes
                if hasproperty(col_type, :parameters) && !isempty(col_type.parameters)
                    return sum(sizeof_arrow_column_type(t) for t in col_type.parameters)
                else
                    return 8
                end
            elseif col_type <: Real
                # Any other numeric type
                return sizeof(col_type)
            else
                return 8  # Default for unknown DataTypes
            end
        elseif col_type isa UnionAll
            # Handle Union types like Union{Missing, T}
            # Just estimate conservatively
            return 8
        else
            return 8  # Default fallback for any other case
        end
    catch
        # If anything fails, return a safe default
        return 8
    end
end

"""
    estimate_arrow_dataframe_memory_mb(file_paths::Vector{String})::Float64

Estimate the total memory required to load all PSM files as DataFrames.
Uses schema-based estimation without loading actual data.

# Arguments
- `file_paths`: Vector of Arrow file paths

# Returns
- Estimated memory in megabytes
"""
function estimate_arrow_dataframe_memory_mb(file_paths::Vector{String})::Float64
    if isempty(file_paths)
        return 0.0
    end

    # Read schema from first file
    first_table = Arrow.Table(file_paths[1])
    schema_names = Tables.columnnames(first_table)

    # Calculate bytes per PSM based on column types
    bytes_per_psm = 0
    for col_name in schema_names
        col = Tables.getcolumn(first_table, col_name)
        col_type = eltype(col)
        bytes_per_psm += sizeof_arrow_column_type(col_type)
    end

    # Count total rows across all files (reads only metadata)
    total_rows = 0
    for fp in file_paths
        table = Arrow.Table(fp)
        total_rows += length(Tables.getcolumn(table, first(schema_names)))
    end

    # Add overhead for DataFrame (column vectors, indices, etc.) - ~20%
    overhead_factor = 1.2

    return (bytes_per_psm * total_rows * overhead_factor) / (1024 * 1024)
end

"""
    calculate_max_sampleable_psms(file_paths::Vector{String}, max_mb::Int)::Int

Calculate the maximum number of PSMs that can be sampled within a memory budget.

# Arguments
- `file_paths`: Vector of Arrow file paths
- `max_mb`: Maximum memory budget in megabytes

# Returns
- Maximum number of PSMs that can fit in the memory budget
"""
function calculate_max_sampleable_psms(file_paths::Vector{String}, max_mb::Int)::Int
    if isempty(file_paths)
        return 0
    end

    # Read schema from first file
    first_table = Arrow.Table(file_paths[1])
    schema_names = Tables.columnnames(first_table)

    # Calculate bytes per PSM
    bytes_per_psm = 0
    for col_name in schema_names
        col = Tables.getcolumn(first_table, col_name)
        col_type = eltype(col)
        bytes_per_psm += sizeof_arrow_column_type(col_type)
    end

    # Account for overhead
    overhead_factor = 1.2
    effective_bytes_per_psm = bytes_per_psm * overhead_factor

    return floor(Int, (max_mb * 1024 * 1024) / effective_bytes_per_psm)
end

"""
    should_use_oom_processing(file_paths::Vector{String}, max_psm_memory_mb::Int,
                               max_psms_in_memory::Int, force_oom::Bool)::Bool

Determine whether to use out-of-memory processing based on memory estimation.

# Arguments
- `file_paths`: Vector of Arrow file paths
- `max_psm_memory_mb`: Memory-based threshold in MB
- `max_psms_in_memory`: Count-based threshold (fallback)
- `force_oom`: If true, always use OOM processing (useful for testing/debugging)

# Returns
- true if OOM processing should be used, false otherwise
"""
function should_use_oom_processing(file_paths::Vector{String}, max_psm_memory_mb::Int,
                                    max_psms_in_memory::Int, force_oom::Bool)::Bool
    # Force OOM for testing/debugging
    if force_oom
        return true
    end

    # Memory-based decision
    estimated_mb = estimate_arrow_dataframe_memory_mb(file_paths)
    if estimated_mb > max_psm_memory_mb
        return true
    end

    # Count-based decision (fallback)
    psms_count = get_psms_count(file_paths)
    return psms_count >= max_psms_in_memory
end

"""
    score_precursor_isotope_traces(second_pass_folder::String,
                                  file_paths::Vector{String},
                                  precursors::LibraryPrecursors,
                                  match_between_runs::Bool,
                                  max_q_value_lightgbm_rescore::Float32,
                                  max_q_value_mbr_itr::Float32,
                                  min_PEP_neg_threshold_itr::Float32,
                                  max_psms_in_memory::Int64,
                                  max_psm_memory_mb::Int64,
                                  force_oom::Bool)

Main entry point for PSM scoring with automatic model selection based on dataset size.

# Processing Strategy
1. If force_oom=true OR memory/count threshold exceeded: Out-of-memory streaming processing
2. PSMs < threshold AND ≥ 200K: In-memory with default/advanced LightGBM
3. PSMs < 200K: In-memory with automatic model comparison

# Arguments
- `second_pass_folder`: Folder containing second pass PSM files
- `file_paths`: Vector of PSM file paths
- `precursors`: Library precursors
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_lightgbm_rescore`: Max q-value for LightGBM rescoring
- `max_q_value_mbr_itr`: Max q-value for MBR transfers retained during iterative training (ITR)
- `min_PEP_neg_threshold_itr`: Min PEP threshold for relabeling weak targets as negatives during ITR
- `max_psms_in_memory`: Maximum PSMs to keep in memory (count-based threshold)
- `max_psm_memory_mb`: Memory-based threshold in MB (default: 5000 = 5GB)
- `force_oom`: If true, always use OOM processing (useful for testing/debugging)
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
    max_psm_memory_mb::Int64,
    force_oom::Bool,
    n_quantile_bins::Int64,
    q_value_threshold::Float32 = 0.01f0,  # Default to 1% if not specified
    ms1_scoring::Bool = true
)
    # Step 1: Count PSMs and determine processing approach
    psms_count = get_psms_count(file_paths)

    # Determine if OOM processing is needed
    use_oom = should_use_oom_processing(file_paths, Int(max_psm_memory_mb), Int(max_psms_in_memory), force_oom)

    if use_oom
        # Case 1: Out-of-memory streaming processing
        estimated_mb = estimate_arrow_dataframe_memory_mb(file_paths)
        if force_oom
            @user_info "Using out-of-memory streaming processing for $psms_count PSMs (FORCED, estimated $(round(estimated_mb, digits=1)) MB)"
        else
            @user_info "Using out-of-memory streaming processing for $psms_count PSMs (estimated $(round(estimated_mb, digits=1)) MB exceeds $(max_psm_memory_mb) MB threshold)"
        end

        # Create OOM model config (uses raw features, no _qbin)
        model_config = create_oom_lightgbm_config(ms1_scoring)

        models = score_precursor_isotope_traces_oom!(
            second_pass_folder,
            file_paths,
            precursors,
            model_config,
            match_between_runs,
            max_q_value_lightgbm_rescore,
            max_q_value_mbr_itr,
            min_PEP_neg_threshold_itr,
            Int(max_psm_memory_mb),
            n_quantile_bins
        )
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
    create_oom_lightgbm_config(ms1_scoring::Bool = true) -> ModelConfig

Creates the LightGBM configuration for out-of-memory processing.
Uses OOM_FEATURE_SET (excludes _qbin features since quantile computation
requires loading all values into memory).
"""
function create_oom_lightgbm_config(ms1_scoring::Bool = true)
    features = copy(OOM_FEATURE_SET)
    apply_ms1_filtering!(features, ms1_scoring)

    return ModelConfig(
        "OOMLightGBM",
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

# DISABLED: Out-of-memory processing - hardcoded to always use in-memory approach
# This function is commented out in favor of always using in-memory processing.
# Preserved for potential future use if needed for extremely large datasets.
#=
"""
    score_precursor_isotope_traces_out_of_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::LibraryPrecursors) -> Dictionary{UInt8, LightGBMModel}

Train LightGBM models for PSM scoring. Only a subset of psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained LightGBM models or simplified model if insufficient PSMs.
"""
function score_precursor_isotope_traces_out_of_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32
)
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    # Features from model_config; do not include :target
    features = [f for f in model_config.features if hasproperty(best_psms, f)];
    if match_between_runs
        append!(features, [
            :MBR_rv_coefficient,
            :MBR_best_irt_diff,
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
        @user_info "OOM LightGBM using $(length(qbin_features)) quantile-binned features: $(join(string.(qbin_features), ", "))"
        for qbin_feat in qbin_features
            n_unique = length(unique(skipmissing(best_psms[!, qbin_feat])))
            @user_info "  $qbin_feat: $n_unique unique values in training sample"
        end
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
                            print_importance = false);
    return models;#best_psms
end
=#

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

#==========================================================
OOM (Out-of-Memory) Scoring Functions
==========================================================#

"""
    score_precursor_isotope_traces_oom!(second_pass_folder, file_paths, precursors,
                                        model_config, match_between_runs, ...)

Out-of-memory scoring pipeline for large datasets. Processes files via streaming,
calculates global q-values across all files, and trains models on sampled data.

# Key Design Decisions
1. **Global Q-Values**: After each iteration, apply model to ALL files, calculate
   GLOBAL q-values, then re-sample based on global q-values
2. **Sampling Unit = pair_id**: Never split pairs between sampled/unsampled data
3. **Training-eligible PSMs only**: For iterations 2+, only load decoys + targets
   passing q-value threshold
4. **Quantile features DISABLED**: OOM version uses raw features (no _qbin)
5. **Arrow.jl approach**: Process files ONE AT A TIME with `Arrow.Table(single_file)`
   for predictable memory control

# Arguments
- `second_pass_folder`: Folder containing second pass PSM files
- `file_paths`: Vector of PSM file paths
- `precursors`: Library precursors
- `model_config`: OOM LightGBM configuration
- `match_between_runs`: Whether to perform match between runs
- `max_q_value_lightgbm_rescore`: Max q-value for training eligibility
- `max_q_value_mbr_itr`: Max q-value for MBR transfers retained during ITR
- `min_PEP_neg_threshold_itr`: Min PEP threshold for negative mining
- `max_memory_mb`: Maximum memory budget in MB
- `n_quantile_bins`: Number of quantile bins (not used in OOM, for API compatibility)

# Returns
- Dict{UInt8, Vector{LightGBM.Booster}} - trained models per CV fold per iteration
"""
function score_precursor_isotope_traces_oom!(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    model_config::ModelConfig,
    match_between_runs::Bool,
    max_q_value_lightgbm_rescore::Float32,
    max_q_value_mbr_itr::Float32,
    min_PEP_neg_threshold_itr::Float32,
    max_memory_mb::Int,
    n_quantile_bins::Int64
)
    # Get hyperparameters from config
    hp = model_config.hyperparams
    iter_scheme = get(hp, :iter_scheme, [100, 200, 200])

    # Calculate target sample size based on memory budget
    target_psm_count = if max_memory_mb > 0
        calculate_max_sampleable_psms(file_paths, max_memory_mb)
    else
        # Default: sample ~500K PSMs
        500_000
    end
    @user_info "OOM: Target sample size: $target_psm_count PSMs"

    # ========================================
    # Step 1: Assign pair_ids to all files (streaming, one-time)
    # ========================================
    @user_info "OOM Step 1: Assigning pair_ids via streaming..."
    precursor_to_pair_id = assign_pair_ids_streaming!(file_paths, precursors)
    @user_info "OOM Step 1 complete: $(length(precursor_to_pair_id)) precursors assigned to pairs"

    # Get features from config (OOM features - no _qbin)
    features = [f for f in model_config.features]
    if match_between_runs
        append!(features, [
            :MBR_max_pair_prob,
            :MBR_log2_weight_ratio,
            :MBR_log2_explained_ratio,
            :MBR_rv_coefficient,
            :MBR_is_missing
        ])
    end
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    # Store models per CV fold - each fold gets a vector of models (one per iteration)
    models = Dict{UInt8, Vector{Any}}()

    # Determine iteration structure
    mbr_start_iter = length(iter_scheme)
    iterations_per_fold = match_between_runs ? length(iter_scheme) : max(mbr_start_iter - 1, 1)

    # Get unique CV folds from first file
    first_table = Arrow.Table(file_paths[1])
    unique_cv_folds = unique(Tables.getcolumn(first_table, :cv_fold))
    first_table = nothing

    # ========================================
    # Iterative Training Loop
    # ========================================
    for (itr, num_rounds) in enumerate(iter_scheme)
        @user_info "OOM Iteration $itr/$iterations_per_fold: Training with $num_rounds rounds"

        # ========================================
        # Step 2: Sample PSMs for this iteration
        # ========================================
        if itr == 1
            # Iteration 1: Sample from ALL PSMs (no q-value filtering)
            @user_info "OOM Step 2 (Iter $itr): Sampling from all PSMs..."
            sampled_psms = sample_psms_by_pairs(
                file_paths, precursor_to_pair_id, target_psm_count, precursors
            )
        else
            # Iterations 2+: Sample from training-eligible PSMs only
            # (decoys + targets passing GLOBAL q-value threshold)
            @user_info "OOM Step 2 (Iter $itr): Sampling from q-value filtered PSMs..."
            sampled_psms = sample_psms_by_pairs_with_qvalue_filter(
                file_paths, precursor_to_pair_id, target_psm_count,
                max_q_value_lightgbm_rescore, precursors
            )
        end
        @user_info "OOM Step 2 (Iter $itr) complete: Sampled $(nrow(sampled_psms)) PSMs"

        # Add required columns
        sampled_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in sampled_psms[!,:precursor_idx]]
        sampled_psms[!,:q_value] = zeros(Float32, nrow(sampled_psms))
        sampled_psms[!,:decoy] = sampled_psms[!,:target] .== false

        # Initialize MBR columns if needed
        if match_between_runs
            initialize_prob_group_features!(sampled_psms, true)
        end

        # ========================================
        # Step 3: Train single iteration model on sampled data
        # ========================================
        @user_info "OOM Step 3 (Iter $itr): Training models..."
        train_feats = itr < mbr_start_iter ? non_mbr_features : features
        train_feats = [f for f in train_feats if hasproperty(sampled_psms, f)]

        for test_fold_idx in unique_cv_folds
            # Get training data (all folds except test fold)
            train_mask = sampled_psms.cv_fold .!= test_fold_idx
            psms_train = sampled_psms[train_mask, :]

            # Apply negative mining for iteration 2+
            if itr >= 2
                psms_train = apply_negative_mining_oom!(
                    psms_train, min_PEP_neg_threshold_itr
                )
            end

            # Train booster
            bst = train_booster(psms_train, train_feats, num_rounds;
                               feature_fraction = get(hp, :feature_fraction, 0.5),
                               learning_rate = get(hp, :learning_rate, 0.05),
                               min_data_in_leaf = get(hp, :min_data_in_leaf, 500),
                               bagging_fraction = get(hp, :bagging_fraction, 0.25),
                               min_gain_to_split = get(hp, :min_gain_to_split, 0.5),
                               max_depth = get(hp, :max_depth, 10),
                               num_leaves = get(hp, :num_leaves, 63))

            # Store model
            if !haskey(models, test_fold_idx)
                models[test_fold_idx] = Vector{Any}()
            end
            push!(models[test_fold_idx], bst)

            # Predict on test fold of sampled data for q-value computation
            test_mask = sampled_psms.cv_fold .== test_fold_idx
            if any(test_mask)
                sampled_psms[test_mask, :trace_prob] = predict(bst, sampled_psms[test_mask, :])
            end
        end

        # Calculate q-values on sampled data (for next iteration's training data selection)
        get_qvalues!(sampled_psms.trace_prob, sampled_psms.target, sampled_psms.q_value)

        # ========================================
        # Step 4: Apply model to ALL files (streaming)
        # ========================================
        @user_info "OOM Step 4 (Iter $itr): Applying models to all files..."
        itr_models = Dict(fold => models[fold][itr] for fold in unique_cv_folds)
        apply_models_streaming!(file_paths, itr_models, train_feats)

        # ========================================
        # Step 5: Calculate GLOBAL q-values across ALL files
        # ========================================
        @user_info "OOM Step 5 (Iter $itr): Computing global q-values..."
        calculate_global_qvalues_streaming!(file_paths)

        # ========================================
        # Step 6: Compute MBR features for next iteration (if needed)
        # ========================================
        if match_between_runs && itr >= mbr_start_iter - 1
            @user_info "OOM Step 6 (Iter $itr): Computing MBR features..."
            compute_mbr_features_streaming!(file_paths, precursor_to_pair_id,
                                            max_q_value_lightgbm_rescore)
        end

        # Clean up sampled data
        sampled_psms = nothing
        GC.gc()

        # Break early if not doing MBR and we've completed non-MBR iterations
        if !match_between_runs && itr == (mbr_start_iter - 1)
            break
        end
    end

    @user_info "OOM scoring complete"
    return models
end

"""
    sample_psms_by_pairs(file_paths, precursor_to_pair_id, target_count, precursors)

Sample PSMs by complete pairs (never split pairs between sampled/unsampled).
For iteration 1, samples from all PSMs.

# Returns
- DataFrame containing sampled PSMs with pair integrity preserved
"""
function sample_psms_by_pairs(
    file_paths::Vector{String},
    precursor_to_pair_id::Dict{UInt32, UInt32},
    target_count::Int,
    precursors::LibraryPrecursors
)
    # Pass 1: Collect all unique pair_ids and count PSMs per pair
    pair_psm_counts = Dict{UInt32, Int}()
    total_psms = 0

    for file_path in file_paths
        table = Arrow.Table(file_path)
        prec_idx_col = Tables.getcolumn(table, :precursor_idx)

        for i in 1:length(prec_idx_col)
            prec_idx = prec_idx_col[i]
            if haskey(precursor_to_pair_id, prec_idx)
                pair_id = precursor_to_pair_id[prec_idx]
                pair_psm_counts[pair_id] = get(pair_psm_counts, pair_id, 0) + 1
                total_psms += 1
            end
        end
    end

    # Calculate sampling probability
    n_pairs = length(pair_psm_counts)
    avg_psms_per_pair = total_psms / max(n_pairs, 1)
    target_pairs = ceil(Int, target_count / max(avg_psms_per_pair, 1))

    # Sample pairs
    all_pair_ids = collect(keys(pair_psm_counts))
    rng = MersenneTwister(1776)
    n_sample = min(target_pairs, length(all_pair_ids))
    sampled_pairs = Set(sample(rng, all_pair_ids, n_sample, replace=false))

    @user_info "OOM Sampling: Selected $(length(sampled_pairs)) pairs from $n_pairs total pairs"

    # Pass 2: Load PSMs belonging to sampled pairs
    result_dfs = DataFrame[]

    for file_path in file_paths
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Filter to sampled pairs
        mask = [haskey(precursor_to_pair_id, pid) &&
                precursor_to_pair_id[pid] in sampled_pairs
                for pid in df.precursor_idx]

        if any(mask)
            filtered_df = df[mask, :]
            # Add pair_id column
            filtered_df[!, :pair_id] = [precursor_to_pair_id[pid] for pid in filtered_df.precursor_idx]
            push!(result_dfs, filtered_df)
        end
    end

    if isempty(result_dfs)
        error("No PSMs sampled - check pair_id assignment")
    end

    return vcat(result_dfs...)
end

"""
    sample_psms_by_pairs_with_qvalue_filter(file_paths, precursor_to_pair_id,
                                             target_count, max_q_value, precursors)

Sample PSMs by complete pairs, but only from training-eligible PSMs.
Training-eligible = decoys OR targets passing q-value threshold.

# Returns
- DataFrame containing sampled training-eligible PSMs with pair integrity preserved
"""
function sample_psms_by_pairs_with_qvalue_filter(
    file_paths::Vector{String},
    precursor_to_pair_id::Dict{UInt32, UInt32},
    target_count::Int,
    max_q_value::Float32,
    precursors::LibraryPrecursors
)
    # Pass 1: Find pair_ids that have training-eligible PSMs
    eligible_pairs = Dict{UInt32, Int}()  # pair_id -> count of eligible PSMs

    for file_path in file_paths
        table = Arrow.Table(file_path)
        prec_idx_col = Tables.getcolumn(table, :precursor_idx)
        target_col = Tables.getcolumn(table, :target)
        qvalue_col = Tables.getcolumn(table, :q_value)

        for i in 1:length(prec_idx_col)
            prec_idx = prec_idx_col[i]
            if !haskey(precursor_to_pair_id, prec_idx)
                continue
            end

            # Training-eligible = decoys OR targets passing q-value
            is_eligible = !target_col[i] || qvalue_col[i] <= max_q_value
            if is_eligible
                pair_id = precursor_to_pair_id[prec_idx]
                eligible_pairs[pair_id] = get(eligible_pairs, pair_id, 0) + 1
            end
        end
    end

    # Calculate sampling
    total_eligible = sum(values(eligible_pairs))
    n_pairs = length(eligible_pairs)
    avg_psms_per_pair = total_eligible / max(n_pairs, 1)
    target_pairs = ceil(Int, target_count / max(avg_psms_per_pair, 1))

    # Sample pairs
    all_pair_ids = collect(keys(eligible_pairs))
    rng = MersenneTwister(1776)
    n_sample = min(target_pairs, length(all_pair_ids))
    sampled_pairs = Set(sample(rng, all_pair_ids, n_sample, replace=false))

    @user_info "OOM Q-filtered Sampling: Selected $(length(sampled_pairs)) pairs from $n_pairs eligible pairs"

    # Pass 2: Load training-eligible PSMs from sampled pairs
    result_dfs = DataFrame[]

    for file_path in file_paths
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Filter to: sampled pairs AND training-eligible
        mask = [
            haskey(precursor_to_pair_id, df.precursor_idx[i]) &&
            precursor_to_pair_id[df.precursor_idx[i]] in sampled_pairs &&
            (!df.target[i] || df.q_value[i] <= max_q_value)
            for i in 1:nrow(df)
        ]

        if any(mask)
            filtered_df = df[mask, :]
            # Add pair_id column
            filtered_df[!, :pair_id] = [precursor_to_pair_id[pid] for pid in filtered_df.precursor_idx]
            push!(result_dfs, filtered_df)
        end
    end

    if isempty(result_dfs)
        error("No training-eligible PSMs sampled - check q-value threshold")
    end

    return vcat(result_dfs...)
end

"""
    apply_negative_mining_oom!(psms_train, min_PEP_threshold)

Apply negative mining to training data: relabel worst-scoring targets as negatives
based on PEP threshold. Must match in-memory rules exactly.

# Returns
- Modified DataFrame with weak targets relabeled as negatives
"""
function apply_negative_mining_oom!(
    psms_train::DataFrame,
    min_PEP_threshold::Float32
)
    # Calculate PEP from current scores
    order = sortperm(psms_train.trace_prob, rev=true)
    sorted_scores = psms_train.trace_prob[order]
    sorted_targets = psms_train.target[order]
    PEPs = Vector{Float32}(undef, length(order))
    get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

    # Find cutoff and relabel
    idx_cutoff = findfirst(x -> x >= min_PEP_threshold, PEPs)
    if !isnothing(idx_cutoff)
        worst_idxs = order[idx_cutoff:end]
        # Only relabel targets, not decoys
        for idx in worst_idxs
            if psms_train.target[idx]
                psms_train.target[idx] = false
            end
        end
        n_relabeled = sum(psms_train.target[worst_idxs] .== false)
        @user_info "OOM Negative mining: Relabeled $n_relabeled weak targets as negatives (PEP ≥ $min_PEP_threshold)"
    end

    return psms_train
end

"""
    apply_models_streaming!(file_paths, models, features)

Apply CV fold models to all files via streaming. Loads each file, applies
appropriate model, writes back trace_prob column.

# Arguments
- `file_paths`: Vector of Arrow file paths
- `models`: Dict{UInt8, Booster} mapping cv_fold -> model
- `features`: Feature columns used for prediction
"""
function apply_models_streaming!(
    file_paths::Vector{String},
    models::Dict{UInt8, T},
    features::Vector{Symbol}
) where T
    for file_path in file_paths
        # Load full DataFrame for editing
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Ensure trace_prob column exists
        if !hasproperty(df, :trace_prob)
            df[!, :trace_prob] = zeros(Float32, nrow(df))
        end

        # Apply model based on CV fold
        for (fold, model) in models
            mask = df.cv_fold .== fold
            if any(mask)
                # Filter features to those present in DataFrame
                available_features = [f for f in features if hasproperty(df, f)]
                df[mask, :trace_prob] = predict(model, df[mask, available_features])
            end
        end

        writeArrow(file_path, df)
    end
end

"""
    collect_score_distribution_streaming(file_paths)

Collect all trace_prob scores from files for global q-value calculation.

# Returns
- (target_scores::Vector{Float32}, decoy_scores::Vector{Float32}) - sorted descending
"""
function collect_score_distribution_streaming(file_paths::Vector{String})
    target_scores = Float32[]
    decoy_scores = Float32[]

    for file_path in file_paths
        table = Arrow.Table(file_path)
        trace_prob_col = Tables.getcolumn(table, :trace_prob)
        target_col = Tables.getcolumn(table, :target)

        for i in 1:length(trace_prob_col)
            if target_col[i]
                push!(target_scores, trace_prob_col[i])
            else
                push!(decoy_scores, trace_prob_col[i])
            end
        end
    end

    sort!(target_scores, rev=true)
    sort!(decoy_scores, rev=true)

    return target_scores, decoy_scores
end

"""
    calculate_global_qvalues_streaming!(file_paths)

Calculate global q-values across all files and write back to each file.
Uses the same algorithm as in-memory get_qvalues! to ensure identical behavior.
"""
function calculate_global_qvalues_streaming!(file_paths::Vector{String})
    # Pass 1: Collect global score distribution
    target_scores, decoy_scores = collect_score_distribution_streaming(file_paths)
    @user_info "OOM Q-values: Collected $(length(target_scores)) target and $(length(decoy_scores)) decoy scores"

    # Pass 2: Write q-values to all files
    for file_path in file_paths
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Calculate q-values using existing function
        qvals = Vector{Float32}(undef, nrow(df))
        get_qvalues!(df.trace_prob, df.target, qvals)
        df[!, :q_value] = qvals

        writeArrow(file_path, df)
    end
end

"""
    compute_mbr_features_streaming!(file_paths, precursor_to_pair_id, max_q_value)

Compute MBR features across all files via streaming.
Two-pass approach: collect pair statistics, then apply features.
"""
function compute_mbr_features_streaming!(
    file_paths::Vector{String},
    precursor_to_pair_id::Dict{UInt32, UInt32},
    max_q_value::Float32
)
    # Pass 1: Collect pair MBR data
    pair_mbr_data = collect_pair_mbr_data_streaming(file_paths, precursor_to_pair_id, max_q_value)

    # Pass 2: Apply MBR features
    apply_mbr_features_streaming!(file_paths, pair_mbr_data, precursor_to_pair_id)
end

"""
    collect_pair_mbr_data_streaming(file_paths, precursor_to_pair_id, max_q_value)

Collect best PSM data per pair per run for MBR feature computation.

# Returns
- Dict mapping (pair_id, isotopes) -> Dict{ms_file_idx -> RunBestPSM}
"""
function collect_pair_mbr_data_streaming(
    file_paths::Vector{String},
    precursor_to_pair_id::Dict{UInt32, UInt32},
    max_q_value::Float32
)
    # Structure: (pair_id, isotopes) -> (ms_file_idx -> best PSM data)
    pair_mbr = Dict{
        @NamedTuple{pair_id::UInt32, isotopes::Tuple{Int8,Int8}},
        Dict{UInt32, @NamedTuple{
            trace_prob::Float32,
            weights::Vector{Float32},
            irts::Vector{Float32},
            weight::Float32,
            log2_intensity_explained::Float32,
            irt_residual::Float32,
            is_decoy::Bool,
            passes_qvalue::Bool
        }}
    }()

    for file_path in file_paths
        table = Arrow.Table(file_path)

        # Get column references
        prec_idx_col = Tables.getcolumn(table, :precursor_idx)
        isotopes_col = Tables.getcolumn(table, :isotopes_captured)
        ms_file_idx_col = Tables.getcolumn(table, :ms_file_idx)
        trace_prob_col = Tables.getcolumn(table, :trace_prob)
        target_col = Tables.getcolumn(table, :target)
        qvalue_col = Tables.getcolumn(table, :q_value)
        weights_col = Tables.getcolumn(table, :weights)
        irts_col = Tables.getcolumn(table, :irts)
        weight_col = Tables.getcolumn(table, :weight)
        l2ie_col = Tables.getcolumn(table, :log2_intensity_explained)
        irt_pred_col = Tables.getcolumn(table, :irt_pred)
        irt_obs_col = Tables.getcolumn(table, :irt_obs)

        for i in 1:length(prec_idx_col)
            prec_idx = prec_idx_col[i]
            if !haskey(precursor_to_pair_id, prec_idx)
                continue
            end

            pair_id = precursor_to_pair_id[prec_idx]
            key = (pair_id = pair_id, isotopes = isotopes_col[i])
            ms_file_idx = ms_file_idx_col[i]
            trace_prob = trace_prob_col[i]

            # Initialize pair entry if needed
            if !haskey(pair_mbr, key)
                pair_mbr[key] = Dict{UInt32, @NamedTuple{
                    trace_prob::Float32, weights::Vector{Float32}, irts::Vector{Float32},
                    weight::Float32, log2_intensity_explained::Float32, irt_residual::Float32,
                    is_decoy::Bool, passes_qvalue::Bool
                }}()
            end

            run_data = pair_mbr[key]

            # Update if this is the best PSM for this run
            if !haskey(run_data, ms_file_idx) || trace_prob > run_data[ms_file_idx].trace_prob
                run_data[ms_file_idx] = (
                    trace_prob = trace_prob,
                    weights = Vector{Float32}(weights_col[i]),
                    irts = Vector{Float32}(irts_col[i]),
                    weight = weight_col[i],
                    log2_intensity_explained = l2ie_col[i],
                    irt_residual = Float32(irt_pred_col[i] - irt_obs_col[i]),
                    is_decoy = !target_col[i],
                    passes_qvalue = qvalue_col[i] <= max_q_value
                )
            end
        end
    end

    return pair_mbr
end

"""
    apply_mbr_features_streaming!(file_paths, pair_mbr, precursor_to_pair_id)

Apply MBR features to all files using collected pair statistics.
"""
function apply_mbr_features_streaming!(
    file_paths::Vector{String},
    pair_mbr::Dict,
    precursor_to_pair_id::Dict{UInt32, UInt32}
)
    for file_path in file_paths
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Initialize MBR columns
        n = nrow(df)
        df[!, :MBR_max_pair_prob] = zeros(Float32, n)
        df[!, :MBR_best_irt_diff] = zeros(Float32, n)
        df[!, :MBR_log2_weight_ratio] = zeros(Float32, n)
        df[!, :MBR_log2_explained_ratio] = zeros(Float32, n)
        df[!, :MBR_rv_coefficient] = zeros(Float32, n)
        df[!, :MBR_is_best_decoy] = trues(n)
        df[!, :MBR_num_runs] = zeros(Int32, n)
        df[!, :MBR_is_missing] = falses(n)

        for i in 1:nrow(df)
            prec_idx = df.precursor_idx[i]
            if !haskey(precursor_to_pair_id, prec_idx)
                df.MBR_is_missing[i] = true
                continue
            end

            pair_id = precursor_to_pair_id[prec_idx]
            key = (pair_id = pair_id, isotopes = df.isotopes_captured[i])

            if !haskey(pair_mbr, key)
                df.MBR_is_missing[i] = true
                continue
            end

            run_data = pair_mbr[key]
            current_run = df.ms_file_idx[i]

            # Count runs with passing PSMs (excluding current run)
            n_passing = sum(data.passes_qvalue for (run_idx, data) in run_data
                           if run_idx != current_run; init=0)
            df.MBR_num_runs[i] = n_passing

            # Find best PSM from a DIFFERENT run
            best_other = nothing
            best_other_prob = -Inf32
            for (run_idx, data) in run_data
                if run_idx != current_run && data.trace_prob > best_other_prob
                    best_other_prob = data.trace_prob
                    best_other = data
                end
            end

            if isnothing(best_other) || n_passing == 0
                df.MBR_is_missing[i] = true
                df.MBR_best_irt_diff[i] = -1.0f0
                df.MBR_rv_coefficient[i] = -1.0f0
                df.MBR_max_pair_prob[i] = -1.0f0
                df.MBR_log2_weight_ratio[i] = -1.0f0
                df.MBR_log2_explained_ratio[i] = -1.0f0
                continue
            end

            # Compute MBR features
            df.MBR_max_pair_prob[i] = best_other.trace_prob
            df.MBR_is_best_decoy[i] = best_other.is_decoy

            # iRT residual difference
            current_residual = Float32(df.irt_pred[i] - df.irt_obs[i])
            df.MBR_best_irt_diff[i] = abs(best_other.irt_residual - current_residual)

            # Weight ratio
            df.MBR_log2_weight_ratio[i] = log2(df.weight[i] / best_other.weight)

            # Intensity explained ratio
            df.MBR_log2_explained_ratio[i] = df.log2_intensity_explained[i] - best_other.log2_intensity_explained

            # RV coefficient
            best_log2_weights = log2.(best_other.weights)
            current_log2_weights = log2.(df.weights[i])
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, current_log2_weights)
            best_irts_padded, irts_padded = pad_rt_equal_length(best_other.irts, df.irts[i])
            df.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_irts_padded,
                                                          weights_padded, irts_padded)
        end

        writeArrow(file_path, df)
    end
end
