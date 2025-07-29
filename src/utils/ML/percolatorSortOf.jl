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

function sort_of_percolator_in_memory!(psms::DataFrame, 
                  file_paths::Vector{String},
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_xgboost_rescore::Float32 = 0.01f0,
                  max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
                  min_PEP_neg_threshold_xgboost_rescore = 0.90f0,
                  colsample_bytree::Float64 = 0.5,
                  colsample_bynode::Float64 = 0.5,
                  eta::Float64 = 0.15,
                  min_child_weight::Int = 1,
                  subsample::Float64 = 0.5,
                  gamma::Int = 0,
                  max_depth::Int = 10,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = true)

    
    #Faster if sorted first
    sort!(psms, [:pair_id, :isotopes_captured])

    prob_estimates = zeros(Float32, nrow(psms))
    MBR_estimates  = zeros(Float32, nrow(psms))

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{EvoTrees.EvoTree}}()
    mbr_start_iter = length(iter_scheme)

    cv_fold_col = psms[!, :cv_fold]
    fold_indices = Dict(fold => findall(==(fold), cv_fold_col) for fold in unique_cv_folds)
    train_indices = Dict(fold => findall(!=(fold), cv_fold_col) for fold in unique_cv_folds)

    Random.seed!(1776)
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    pbar = ProgressBar(total=length(unique_cv_folds) * length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        initialize_prob_group_features!(psms, match_between_runs)
        psms_train = @view psms[train_indices[test_fold_idx], :]
        test_fold_idxs = fold_indices[test_fold_idx]
        test_fold_psms = @view psms[test_fold_idxs, :]
        fold_models = Vector{EvoTrees.EvoTree}(undef, length(iter_scheme))

        for (itr, num_round) in enumerate(iter_scheme)
            psms_train_itr = get_training_data_for_iteration!(psms_train,
                                                              itr,
                                                              match_between_runs,
                                                              max_q_value_xgboost_rescore,
                                                              max_q_value_xgboost_mbr_rescore,
                                                              min_PEP_neg_threshold_xgboost_rescore,
                                                              itr >= mbr_start_iter)

            train_feats = itr < mbr_start_iter ? non_mbr_features : features
            
            bst = train_booster(psms_train_itr, train_feats, num_round;
                               colsample_bytree=colsample_bytree,
                               colsample_bynode=colsample_bynode,
                               eta=eta,
                               min_child_weight=min_child_weight,
                               subsample=subsample,
                               gamma=gamma,
                               max_depth=max_depth)
            fold_models[itr] = bst

            predict_fold!(bst, psms_train, test_fold_psms, train_feats)

            if match_between_runs
                update_mbr_features!(psms_train, test_fold_psms, prob_estimates,
                                     test_fold_idxs, itr, mbr_start_iter,
                                     max_q_value_xgboost_rescore)
            end

            update(pbar)

            if (!match_between_runs) && itr == (mbr_start_iter - 1)
                break
            end
        end
        # Make predictions on hold out data.
        if match_between_runs
            MBR_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        else
            prob_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        end
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end
    psms[!,:prob] = prob_estimates

    if match_between_runs
        psms[!, :MBR_prob] = MBR_estimates
        @. psms.MBR_prob = clamp(
            psms.MBR_prob,
            0f0,
            max(psms.prob, coalesce(psms.MBR_max_pair_prob, psms.prob))
        )
    end
    
    dropVectorColumns!(psms) # avoids writing issues
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    return models
end




function sort_of_percolator_out_of_memory!(psms::DataFrame, 
                    file_paths::Vector{String},
                    features::Vector{Symbol},
                    match_between_runs::Bool = true; 
                    max_q_value_xgboost_rescore::Float32 = 0.01f0,
                    max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
                    min_PEP_neg_threshold_xgboost_rescore::Float32 = 0.90f0,
                    colsample_bytree::Float64 = 0.5, 
                    colsample_bynode::Float64 = 0.5,
                    eta::Float64 = 0.15, 
                    min_child_weight::Int = 1, 
                    subsample::Float64 = 0.5, 
                    gamma::Int = 0, 
                    max_depth::Int = 10,
                    iter_scheme::Vector{Int} = [100, 200, 200],
                    print_importance::Bool = true)

    function getBestScorePerPrec!(
        prec_to_best_score_new::Dictionary,
        file_paths::Vector{String},
        models::Dictionary{UInt8,EvoTrees.EvoTree},
        features::Vector{Symbol},
        match_between_runs::Bool;
        is_last_iteration::Bool = false)
    
        # Reset counts for new scores
        reset_precursor_scores!(prec_to_best_score_new)
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            
            probs = predict_cv_models(models, psms_subset, features)
            
            if match_between_runs && !is_last_iteration
                #Update maximum probabilities for tracked precursors 
                qvals = zeros(Float32, nrow(psms_subset))
                get_qvalues!(probs, psms_subset.target, qvals)

                for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                    prob = probs[i]
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]
    
                        if prob > scores.best_prob_1
                           new_scores = merge(scores, (
                                # replace best_prob_2 with best_prob_1
                                best_prob_2                     = scores.best_prob_1,   
                                best_log2_weights_2             = scores.best_log2_weights_1,
                                best_irts_2                     = scores.best_irts_1,
                                best_weight_2                   = scores.best_weight_1,
                                best_log2_intensity_explained_2 = scores.best_log2_intensity_explained_1,
                                best_ms_file_idx_2              = scores.best_ms_file_idx_1,
                                is_best_decoy_2                 = scores.is_best_decoy_1,
                                # overwrite best_prob_1
                                best_prob_1                     = prob,                
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_irts_1                     = psms_subset.irts[i],
                                best_weight_1                   = psms_subset.weight[i],
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                is_best_decoy_1                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores

                        elseif prob > scores.best_prob_2
                            # overwrite best_prob_2
                            new_scores = merge(scores, (
                                best_prob_2                     = prob,
                                best_log2_weights_2             = log2.(psms_subset.weights[i]),
                                best_irts_2                     = psms_subset.irts[i],
                                best_weight_2                   = psms_subset.weight[i],
                                best_log2_intensity_explained_2 = psms_subset.log2_intensity_explained[i],
                                best_ms_file_idx_2              = psms_subset.ms_file_idx[i],
                                is_best_decoy_2                 = psms_subset.decoy[i]
                            ))
                            prec_to_best_score_new[key] = new_scores
                        end

                        if qvals[i] <= max_q_value_xgboost_rescore
                            push!(scores.unique_passing_runs, psms_subset.ms_file_idx[i])
                        end

                    else
                        insert!(prec_to_best_score_new, key, (
                                best_prob_1                     = prob,
                                best_prob_2                     = zero(Float32),
                                best_log2_weights_1             = log2.(psms_subset.weights[i]),
                                best_log2_weights_2             = Vector{Float32}(),
                                best_irts_1                     = psms_subset.irts[i],
                                best_irts_2                     = Vector{Float32}(),
                                best_weight_1                   = psms_subset.weight[i],
                                best_weight_2                   = zero(Float32),
                                best_log2_intensity_explained_1 = psms_subset.log2_intensity_explained[i],
                                best_log2_intensity_explained_2 = zero(Float32),
                                best_ms_file_idx_1              = psms_subset.ms_file_idx[i],
                                best_ms_file_idx_2              = zero(UInt32),
                                is_best_decoy_1                 = psms_subset.decoy[i],
                                is_best_decoy_2                 = false,
                                unique_passing_runs             = ( qvals[i] <= max_q_value_xgboost_rescore ?
                                                                    Set{UInt16}([psms_subset.ms_file_idx[i]]) :
                                                                    Set{UInt16}() )
                            ))
                    end
                end
            end


            if is_last_iteration
                if match_between_runs
                    clamp_mbr_probs!(psms_subset, probs)
                else
                    psms_subset.prob = probs
                end

            end
        end
    

        # Compute probs and features for next round
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = predict_cv_models(models, psms_subset, features)

            for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                psms_subset[i,:prob] = probs[i]
                
                    
                if match_between_runs && !is_last_iteration
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        psms_subset.MBR_num_runs[i] = length(scores.unique_passing_runs)

                        best_log2_weights = Float32[]
                        best_irts = Float32[]
                        best_weight = zero(Float32)
                        best_log2_ie = zero(Float32)

                        if (scores.best_ms_file_idx_1 != psms_subset.ms_file_idx[i]) &&
                           (!isempty(scores.best_log2_weights_1))
                            best_log2_weights                   = scores.best_log2_weights_1
                            best_irts                           = scores.best_irts_1
                            best_weight                         = scores.best_weight_1
                            best_log2_ie                        = scores.best_log2_intensity_explained_1
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_1
                            MBR_is_best_decoy                   = scores.is_best_decoy_1
                        elseif (scores.best_ms_file_idx_2 != psms_subset.ms_file_idx[i]) &&
                               (!isempty(scores.best_log2_weights_2))
                            best_log2_weights                   = scores.best_log2_weights_2
                            best_irts                           = scores.best_irts_2
                            best_weight                         = scores.best_weight_2
                            best_log2_ie                        = scores.best_log2_intensity_explained_2
                            psms_subset.MBR_max_pair_prob[i]    = scores.best_prob_2
                            MBR_is_best_decoy                   = scores.is_best_decoy_2
                        else
                            psms_subset.MBR_best_irt_diff[i]        = -1.0f0
                            psms_subset.MBR_rv_coefficient[i]       = -1.0f0
                            psms_subset.MBR_is_best_decoy[i]        = false
                            psms_subset.MBR_max_pair_prob[i]        = -1.0f0
                            psms_subset.MBR_log2_weight_ratio[i]    = -1.0f0
                            psms_subset.MBR_log2_explained_ratio[i] = -1.0f0
                            psms_subset.sub_psms.MBR_is_missing[i]  = true
                            continue
                        end

                        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(psms_subset.weights[i]))
                        best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_irts, psms_subset.irts[i])

                        best_irt_at_apex = best_irts[argmax(best_log2_weights)]
                        psms_subset.MBR_best_irt_diff[i] = abs(best_irt_at_apex - psms_subset.irts[i][argmax(psms_subset.weights[i])])
                        psms_subset.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
                        psms_subset.MBR_log2_weight_ratio[i] = log2(psms_subset.weight[i] / best_weight)
                        psms_subset.MBR_log2_explained_ratio[i] = psms_subset.log2_intensity_explained[i] - best_log2_ie
                        psms_subset.MBR_is_best_decoy[i] = MBR_is_best_decoy
                    end
                end
            end

            write_subset(file_path, psms_subset, probs, match_between_runs; dropVectors=is_last_iteration)
        end
        
        return prec_to_best_score_new
    end


    unique_cv_folds = unique(psms[!, :cv_fold])
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{EvoTrees.EvoTree}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    Random.seed!(1776);
    non_mbr_features = [f for f in features if !startswith(String(f), "MBR_")]

    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        #Clear prob stats 
        initialize_prob_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])

        for (itr, num_round) in enumerate(iter_scheme)

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_xgboost_rescore,
                                                                max_q_value_xgboost_mbr_rescore,
                                                                min_PEP_neg_threshold_xgboost_rescore,
                                                                itr >= length(iter_scheme))
            ###################
            #Train a model on the n-1 training folds.
            train_feats = itr < length(iter_scheme) ? non_mbr_features : features
            bst = train_booster(psms_train_itr, train_feats, num_round;
                               colsample_bytree=colsample_bytree,
                               colsample_bynode=colsample_bynode,
                               eta=eta,
                               min_child_weight=min_child_weight,
                               subsample=subsample,
                               gamma=gamma,
                               max_depth=max_depth)
            if !haskey(models, test_fold_idx)
                insert!(
                    models,
                    test_fold_idx,
                    Vector{EvoTrees.EvoTree}([bst])
                )
            else
                push!(models[test_fold_idx], bst)
            end
            #print_importance = true
            print_importance ? println(collect(zip(feature_importance(bst)))[1:30]) : nothing

            # Get probabilities for training sample so we can get q-values
            psms_train[!,:prob] = EvoTrees.predict(bst, psms_train)
            
            if match_between_runs
                summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
            end

            update(pbar)
        end
    end
    
    @info "Scoring precursors with trained models"
    pbar = ProgressBar(total=length(iter_scheme))
    prec_to_best_score = Dictionary{@NamedTuple{pair_id::UInt32,
                                                isotopes::Tuple{Int8,Int8}},
                                    @NamedTuple{best_prob_1::Float32,
                                                best_prob_2::Float32,
                                                best_log2_weights_1::Vector{Float32},
                                                best_log2_weights_2::Vector{Float32},
                                                best_irts_1::Vector{Float32},
                                                best_irts_2::Vector{Float32},
                                                best_weight_1::Float32,
                                                best_weight_2::Float32,
                                                best_log2_intensity_explained_1::Float32,
                                                best_log2_intensity_explained_2::Float32,
                                                best_ms_file_idx_1::UInt32,
                                                best_ms_file_idx_2::UInt32,
                                                is_best_decoy_1::Bool,
                                                is_best_decoy_2::Bool,
                                                unique_passing_runs::Set{UInt16}}}()

    for (train_iter, num_round) in enumerate(iter_scheme)
        models_for_iter = Dictionary{UInt8,EvoTrees.EvoTree}()
        for test_fold_idx in unique_cv_folds
            insert!(models_for_iter, test_fold_idx, models[test_fold_idx][train_iter])
        end
        prec_to_best_score = getBestScorePerPrec!(
            prec_to_best_score,
            file_paths,
            models_for_iter,
            features,
            match_between_runs;
            is_last_iteration = (train_iter == length(iter_scheme))
        )
        update(pbar)
    end

    return models
end














function train_booster(psms::AbstractDataFrame, features, num_round;
                       colsample_bytree::Float64,
                       colsample_bynode::Float64,
                       eta::Float64,
                       min_child_weight::Int,
                       subsample::Float64,
                       gamma::Int,
                       max_depth::Int)

    config = EvoTreeRegressor(
        loss=:logloss,
        nrounds = num_round,
        max_depth = max_depth,
        min_weight = min_child_weight,
        rowsample = subsample,
        colsample = colsample_bytree,
        eta = eta,
        gamma = gamma
    )
    model = fit(config, psms; target_name = :target, feature_names = features)
    return model
end

function predict_fold!(bst, psms_train::AbstractDataFrame,
                       test_fold_psms::AbstractDataFrame, features)
    test_fold_psms[!, :prob] = predict(bst, test_fold_psms)
    psms_train[!, :prob] = predict(bst, psms_train)
    get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)
end

function update_mbr_features!(psms_train::AbstractDataFrame,
                              test_fold_psms::AbstractDataFrame,
                              prob_estimates::Vector{Float32},
                              test_fold_idxs,
                              itr::Int,
                              mbr_start_iter::Int,
                              max_q_value_xgboost_rescore::Float32)
    if itr >= mbr_start_iter - 1
        get_qvalues!(test_fold_psms.prob, test_fold_psms.target, test_fold_psms.q_value)
        summarize_precursors!(test_fold_psms, q_cutoff = max_q_value_xgboost_rescore)
        summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
    end
    if itr == mbr_start_iter - 1
        prob_estimates[test_fold_idxs] = test_fold_psms.prob
    end
end

function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)
    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
    Threads.@threads for idx in eachindex(pair_groups)
        _, sub_psms = pair_groups[idx]
        
        # Efficient way to find the top 2 precursors so we can do MBR on the 
        # best precursor match that isn't itself. It's always one of the top 2.

        # single pass: record the best PSM index & prob per run
        offset = Int(minimum(sub_psms.ms_file_idx))
        range_len = Int(maximum(sub_psms.ms_file_idx)) - offset + 1
        best_i = zeros(Int, range_len)
        best_p = fill(-Inf, range_len)
        for (i, run) in enumerate(sub_psms.ms_file_idx)
            idx = Int(run) - offset + 1
            p = sub_psms.prob[i]
            if p > best_p[idx]
                best_p[idx] = p
                best_i[idx] = i
            end
        end

        # if more than one run, find the global top-2 runs by their best-PSM prob
        run_best_indices = zeros(Int, range_len)
        runs = findall(!=(0), best_i)
        if length(runs) > 1
            # track top two runs (r1 > r2)
            r1 = 0; p1 = -Inf
            r2 = 0; p2 = -Inf
            for r in runs
                p = best_p[r]
                if p > p1
                    r2, p2 = r1, p1
                    r1, p1 = r, p
                elseif p > p2
                    r2, p2 = r, p
                end
            end

            # assign, for each run, the best index in “any other” run
            for r in runs
                run_best_indices[r] = (r == r1 ? best_i[r2] : best_i[r1])
            end
        end

        # Compute MBR features
        for i in 1:nrow(sub_psms)
            sub_psms.MBR_num_runs[i] = length(unique(sub_psms.ms_file_idx[sub_psms.q_value .<= q_cutoff]))

            idx = Int(sub_psms.ms_file_idx[i]) - offset + 1
            best_idx = run_best_indices[idx]
            if best_idx == 0
                sub_psms.MBR_best_irt_diff[i]           = -1.0f0
                sub_psms.MBR_rv_coefficient[i]          = -1.0f0
                sub_psms.MBR_is_best_decoy[i]           = false
                sub_psms.MBR_log2_weight_ratio[i]       = -1.0f0
                sub_psms.MBR_log2_explained_ratio[i]    = -1.0f0
                sub_psms.MBR_max_pair_prob[i]           = -1.0f0
                sub_psms.MBR_is_missing[i]              = true
                continue
            end

            best_log2_weights = log2.(sub_psms.weights[best_idx])
            best_iRTs = sub_psms.irts[best_idx]
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(sub_psms.weights[i]))
            best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_iRTs, sub_psms.irts[i])
            
            best_irt_at_apex = sub_psms.irts[best_idx][argmax(best_log2_weights)]
            sub_psms.MBR_max_pair_prob[i] = sub_psms.prob[best_idx]
            sub_psms.MBR_best_irt_diff[i] = abs(best_irt_at_apex - sub_psms.irts[i][argmax(sub_psms.weights[i])])
            sub_psms.MBR_rv_coefficient[i] = MBR_rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
            sub_psms.MBR_log2_weight_ratio[i] = log2(sub_psms.weight[i] / sub_psms.weight[best_idx])
            sub_psms.MBR_log2_explained_ratio[i] = sub_psms.log2_intensity_explained[i] - sub_psms.log2_intensity_explained[best_idx]
            sub_psms.MBR_is_best_decoy[i] = sub_psms.decoy[best_idx]
        end
    end

end




function initialize_prob_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    n = nrow(psms)
    psms[!, :prob]      = zeros(Float32, n)
    psms[!, :q_value]   = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_prob]                      = zeros(Float32, n)
        psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
        psms[!, :MBR_best_irt_diff]             = zeros(Float32, n)
        psms[!, :MBR_log2_weight_ratio]         = zeros(Float32, n)
        psms[!, :MBR_log2_explained_ratio]      = zeros(Float32, n)
        psms[!, :MBR_rv_coefficient]            = zeros(Float32, n)
        psms[!, :MBR_is_best_decoy]             = trues(n)
        psms[!, :MBR_num_runs]                  = zeros(Int32, n)
        psms[!, :MBR_transfer_candidate]        = falses(n)
        psms[!, :MBR_is_missing]                = falses(n)
    end

    return psms
end

function get_training_data_for_iteration!(
    psms_train::AbstractDataFrame,
    itr::Int,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32,
    min_PEP_neg_threshold_xgboost_rescore::Float32,
    last_iter::Bool
)
   
    if itr == 1
        # Train on all precursors during first iteration. 
        return psms_train
    else
        # Do a shallow copy to avoid overwriting target/decoy labels
        psms_train_itr = copy(psms_train)

        # Convert the worst-scoring targets to negatives using PEP estimate
        order = sortperm(psms_train_itr.prob, rev=true)
        sorted_scores  = psms_train_itr.prob[order]
        sorted_targets = psms_train_itr.target[order]
        PEPs = Vector{Float32}(undef, length(order))
        get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

        idx_cutoff = findfirst(x -> x >= min_PEP_neg_threshold_xgboost_rescore, PEPs)
        if !isnothing(idx_cutoff)
            worst_idxs = order[idx_cutoff:end]
            psms_train_itr.target[worst_idxs] .= false
        end

        # Also train on top scoring MBR candidates if requested
        if match_between_runs && last_iter
            # Determine prob threshold for precursors passing the q-value threshold
            max_prob_threshold = minimum(
                psms_train_itr.prob[
                    psms_train_itr.target .& (psms_train_itr.q_value .<= max_q_value_xgboost_rescore)
                ]
            )

            # Hacky way to ensure anything passing the initial q-value threshold
            # will pass the next q-value threshold
            psms_train_itr.q_value[psms_train_itr.q_value .<= max_q_value_xgboost_rescore] .= 0.0
            psms_train_itr.q_value[psms_train_itr.q_value .> max_q_value_xgboost_rescore]  .= 1.0

            # Must have at least one precursor passing the q-value threshold,
            # and the best precursor can't be a decoy
            psms_train_mbr = subset(
                psms_train_itr,
                [:MBR_is_best_decoy, :MBR_max_pair_prob, :prob, :MBR_is_missing] => ByRow((d, mp, p, im) ->
                    (!im && !d && mp >= max_prob_threshold && p < max_prob_threshold)
                );
                view = true
            )

            # Compute MBR q-values.
            get_qvalues!(psms_train_mbr[!,:prob], psms_train_mbr[!,:target], psms_train_mbr[!,:q_value])

            # Take all decoys and targets passing q_thresh (all 0's now) or mbr_q_thresh
            psms_train_itr = subset(
                psms_train_itr,
                [:target, :q_value, :MBR_is_best_decoy, :MBR_is_missing] => ByRow((t, q, MBR_d, im) -> (!t) || (t && !im && !MBR_d && q <= max_q_value_xgboost_mbr_rescore))
            )
        else
            # Take all decoys and targets passing q_thresh
            psms_train_itr = subset(
                psms_train_itr,
                [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_rescore))
            )
        end

        return psms_train_itr
    end
end



function dropVectorColumns!(df)
    to_drop = String[]
    for col in names(df)
        if eltype(df[!, col]) <: AbstractVector
            push!(to_drop, col)
        end
    end
    # 2) Drop those columns in place
    select!(df, Not(to_drop))
end

"""
    reset_precursor_scores!(dict)

Set all values of `dict` to an empty precursor score tuple.  This helps reuse
the same dictionary between iterations without reallocating memory.
"""
function reset_precursor_scores!(dict)
    for key in keys(dict)
        dict[key] = (
            best_prob_1 = zero(Float32),
            best_prob_2 = zero(Float32),
            best_log2_weights_1 = Vector{Float32}(),
            best_log2_weights_2 = Vector{Float32}(),
            best_irts_1 = Vector{Float32}(),
            best_irts_2 = Vector{Float32}(),
            best_weight_1 = zero(Float32),
            best_weight_2 = zero(Float32),
            best_log2_intensity_explained_1 = zero(Float32),
            best_log2_intensity_explained_2 = zero(Float32),
            best_ms_file_idx_1 = zero(UInt32),
            best_ms_file_idx_2 = zero(UInt32),
            is_best_decoy_1 = false,
            is_best_decoy_2 = false,
            unique_passing_runs = Set{UInt16}(),
        )
    end
    return dict
end

"""
    predict_cv_models(models, df, features)

Return a vector of probabilities for `df` using the cross validation `models`.
"""
function predict_cv_models(models::Dictionary{UInt8,EvoTrees.EvoTree},
                           df::AbstractDataFrame,
                           features::Vector{Symbol})
    probs = zeros(Float32, nrow(df))
    for (fold_idx, bst) in pairs(models)
        fold_rows = findall(==(fold_idx), df[!, :cv_fold])
        if !isempty(fold_rows)
            probs[fold_rows] = predict(bst, df[fold_rows, :])
        end
    end
    return probs
end

"""
    clamp_mbr_probs!(df, probs)

Clamp probabilities for MBR search so they never exceed the current best
probabilities for each PSM.
"""
function clamp_mbr_probs!(df::AbstractDataFrame, probs::AbstractVector{Float32})
    df[!, :MBR_prob] = probs
    @. df.MBR_prob = clamp(
        df.MBR_prob,
        0f0,
        max(df.prob, coalesce(df.MBR_max_pair_prob, df.prob)),
    )
    return df
end

"""
    write_subset(file_path, df, probs, match_between_runs; dropVectors=false)

Write the updated subset to disk, optionally dropping vector columns.
"""
function write_subset(file_path::String,
                      df::DataFrame,
                      probs::AbstractVector{Float32},
                      match_between_runs::Bool;
                      dropVectors::Bool=false)
    if dropVectors
        if match_between_runs
            clamp_mbr_probs!(df, probs)
        else
            df[!, :prob] = probs
        end
        writeArrow(file_path, dropVectorColumns!(df))
    else
        df[!, :prob] = probs
        writeArrow(file_path, convert_subarrays(df))
    end
end



function MBR_rv_coefficient(weights_A::AbstractVector{<:Real},
    times_A::AbstractVector{<:Real},
    weights_B::AbstractVector{<:Real},
    times_B::AbstractVector{<:Real})

    # Construct two Nx2 matrices, each row is (weight, time)
    X = hcat(collect(weights_A), collect(times_A))
    Y = hcat(collect(weights_B), collect(times_B))

    # Compute cross-products (Gram matrices)
    Sx = X' * X
    Sy = Y' * Y

    # Numerator: trace(Sx * Sy)
    numerator = tr(Sx * Sy)

    # Denominator: sqrt( trace(Sx*Sx)* trace(Sy*Sy) )
    denominator = sqrt(tr(Sx * Sx) * tr(Sy * Sy))

    # Protect against zero in denominator (e.g. if X or Y is all zeros)
    if denominator == 0
        return 0.0
    end

    return numerator / denominator
end

function pad_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, d)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        return padded_data
    end

    d = length(y) - length(x)
    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, d), y)
    else
        # Pad y
        return (x, pad(y, abs(d)))
    end
end


function pad_rt_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, n, d, rt_step)
        left  = div(d, 2) 
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        @inbounds @fastmath for i in range(1, left)
            padded_data[i] = padded_data[left+1] - (rt_step * ((left+1) - i))
        end 
        @inbounds @fastmath for i in range(1, right)
            padded_data[i+left+n] = padded_data[left+n] + (rt_step * i)
        end 
        return padded_data
    end

    nx = length(x)
    ny = length(y)
    d = ny - nx
    rt_step_x = (x[end] - x[1]) / nx
    rt_step_y = (y[end] - y[1]) / ny
    rt_step = nx > 1 ? rt_step_x : rt_step_y

    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, nx, d, rt_step), y)
    else
        # Pad y
        return (x, pad(y, ny, abs(d), rt_step))
    end
end

function summarize_prob(probs::AbstractVector{Float32})
    minimum(probs), maximum(probs), mean(probs)
end

function convert_subarrays(df::DataFrame)
    for col_name in names(df)
        col_data = df[!, col_name]
        # If it's a SubArray, let's convert it to a plain vector:
        if eltype(col_data) <: AbstractVector
           df[!, col_name] = [Vector(col_data[i]) for i in eachindex(col_data)]
        end
    end
    return df
end
