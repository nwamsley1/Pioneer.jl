#==========================================================
PSM sampling and scoring 
==========================================================#
"""
     get_psms_count(quant_psms_folder::String)::Integer

Sample PSMs from multiple files for XGBoost model training.

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

Sample PSMs from multiple files for XGBoost model training.

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

Loads all PSMs from multiple files for XGBoost model training.
"""
function load_psms_for_xgboost(quant_psms_folder::String)
    file_paths = [fpath for fpath in readdir(quant_psms_folder, join=true) if endswith(fpath,".arrow")]
    return DataFrame(Tables.columntable(Arrow.Table(file_paths)))
end

"""
    score_precursor_isotope_traces_in_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::BasicLibraryPrecursors) -> XGBoostModels

Train XGBoost models for PSM scoring. All psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained XGBoost models or simplified model if insufficient PSMs.
"""
function score_precursor_isotope_traces_in_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::BasicLibraryPrecursors
)
    if size(best_psms, 1) > 100000
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    features = [ 
        :max_prob,
        :mean_prob,
        :min_prob,
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
    ];
    
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
    best_psms[!,:decoy] = best_psms[!,:target].==false;

    models = sort_of_percolator_in_memory!(
                            best_psms, 
                            file_paths,
                            features,
                            colsample_bytree = 0.5, 
                            colsample_bynode = 0.5,
                            min_child_weight = 5, 
                            gamma = 1,
                            subsample = 0.25, 
                            max_depth = 10,
                            eta = 0.05, 
                            iter_scheme = [100, 100, 200],
                            print_importance = false);
    return models;#best_psms
    else
        @warn "Less than 1,000,000 psms. Training with simplified target-decoy discrimination model..."
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
        #Train XGBoost model to score each precursor trace. Target-decoy descrimination
        models = sort_of_percolator_in_memory!(
                                best_psms, 
                                file_paths,
                                features,
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

"""
    score_precursor_isotope_traces_out_of_memory!(best_psms::DataFrame, file_paths::Vector{String},
                                  precursors::BasicLibraryPrecursors) -> XGBoostModels

Train XGBoost models for PSM scoring. Only a subset of psms are kept in memory

# Arguments
- `best_psms`: Sample of high-quality PSMs for training
- `file_paths`: Paths to PSM files
- `precursors`: Library precursor information

# Returns
Trained XGBoost models or simplified model if insufficient PSMs.
"""
function score_precursor_isotope_traces_out_of_memory!(
    best_psms::DataFrame,
    file_paths::Vector{String},
    precursors::BasicLibraryPrecursors
)
    if size(best_psms, 1) > 100000
    file_paths = [fpath for fpath in file_paths if endswith(fpath,".arrow")]
    features = [ 
        :max_prob,
        :mean_prob,
        :min_prob,
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
    ];
    
    best_psms[!,:accession_numbers] = [getAccessionNumbers(precursors)[pid] for pid in best_psms[!,:precursor_idx]]
    best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
    best_psms[!,:decoy] = best_psms[!,:target].==false;

    models = sort_of_percolator_out_of_memory!(
                            best_psms, 
                            file_paths,
                            features,
                            colsample_bytree = 0.5, 
                            colsample_bynode = 0.5,
                            min_child_weight = 5, 
                            gamma = 1,
                            subsample = 0.25, 
                            max_depth = 10,
                            eta = 0.05, 
                            iter_scheme = [100, 100, 200],
                            print_importance = false);
    return models;#best_psms
    else
        @warn "Less than 1,000,000 psms. Training with simplified target-decoy discrimination model..."
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
        #Train XGBoost model to score each precursor trace. Target-decoy descrimination
        models = sort_of_percolator_out_of_memory!(
                                best_psms, 
                                max_train_psms,
                                file_paths,
                                features,
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