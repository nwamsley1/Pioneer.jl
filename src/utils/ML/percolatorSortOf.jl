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
                  q_value_threshold::Float32 = 0.01f0,
                  max_ftr::Float32 = 0.01f0,
                  colsample_bytree::Float64 = 0.5,
                  colsample_bynode::Float64 = 0.5,
                  eta::Float64 = 0.15,
                  min_child_weight::Int = 1,
                  subsample::Float64 = 0.5,
                  gamma::Int = 0,
                  max_depth::Int = 10,
                  iter_scheme::Vector{Int} = [100, 200, 200],
                  print_importance::Bool = true)
    
    
    # Helper functions for fold selection
    selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold == test_fold
    excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold != test_fold

    #Faster if sorted first 
    @info "sort time..."
    sort!(psms, [:pair_id, :isotopes_captured])

    #Final prob estimates
    prob_estimates = zeros(Float32, size(psms, 1))
    prob_estimates_stage1 = zeros(Float32, size(psms, 1))
    

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{Booster}}()
    # Always doing MBR on the last iteration
    mbr_start_iter = length(iter_scheme)
    
    # Train models for each fold
    Random.seed!(1776)
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        # Clear prob stats 
        initialize_prob_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])
        # Train models for each iteration
        fold_models = Vector{Booster}()

        for (itr, num_round) in enumerate(iter_scheme)

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_xgboost_rescore,
                                                                max_q_value_xgboost_mbr_rescore,
                                                                min_PEP_neg_threshold_xgboost_rescore,
                                                                itr >= mbr_start_iter)
                                         
            bst = xgboost(
                (psms_train_itr[!, features], psms_train_itr[!, :target]),
                num_round = num_round,
                colsample_bytree = colsample_bytree,
                colsample_bynode = colsample_bynode,
                scale_pos_weight = sum(psms_train_itr.decoy) / sum(psms_train_itr.target),
                gamma = gamma,
                max_depth = max_depth,
                eta = eta,
                min_child_weight = min_child_weight,
                subsample = subsample,
                objective = "binary:logistic",
                seed = rand(UInt32),
                watchlist = (;)
            )
            # Store feature names and print importance if requested
            bst.feature_names = string.(features)
            #print_importance = true
            if print_importance
                println(collect(zip(importance(bst))))
            end
            
            push!(fold_models, bst)
            test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
            test_fold_psms = @view(psms[test_fold_idxs,:])

            test_fold_psms[!,:prob] = XGBoost.predict(bst, test_fold_psms[!,features])
            psms_train[!,:prob] =  XGBoost.predict(bst, psms_train[!, features])

            # calculate training q-values to use for filtering the next iteration
            get_qvalues!(psms_train.prob, psms_train.target, psms_train.q_value)

            if match_between_runs
                # Compute the MBR features if we're going to use them on the next iteration
                if itr >= mbr_start_iter - 1
                    # calculate q-values on hold-out to properly summarize precursors
                    get_qvalues!(test_fold_psms.prob, test_fold_psms.target, test_fold_psms.q_value)
                    summarize_precursors!(test_fold_psms, q_cutoff = max_q_value_xgboost_rescore)
                    summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
                end
                # Keep track of the last non-MBR probabilites so we know which precursors are transfer candidates
                if itr == mbr_start_iter - 1
                    prob_estimates_stage1[test_fold_idxs] = test_fold_psms.prob
                end
            end

            update(pbar)

            # If we're not doing MBR, then skip the last iteration
            if (!match_between_runs) && itr == (mbr_start_iter - 1)
                break
            end
        end
        # Make predictions on hold out data.
        test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
        prob_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end
    psms[!,:prob] = prob_estimates
    
    if match_between_runs
        # Label transfer candidates using the probabilities prior to MBR features
        stage1_qvals = zeros(Float64, length(prob_estimates_stage1))
        get_qvalues!(prob_estimates_stage1, psms.target, stage1_qvals)
        psms.MBR_transfer_candidate .= (stage1_qvals .> q_value_threshold) .& !ismissing(psms.MBR_is_best_decoy)

        # Estimate FTR and filter
        is_bad_transfer = (psms.target .& [(!ismissing(d) && d) for d in psms.MBR_is_best_decoy]) .| (psms.decoy .& [(!ismissing(d) && !d) for d in psms.MBR_is_best_decoy])
        τ = get_ftr_threshold(psms.prob, psms.target, is_bad_transfer, max_ftr; mask=psms.MBR_transfer_candidate)
        filter!(
            row -> !(row.MBR_transfer_candidate && row.prob < τ),
            psms,
        )
    end

    transform!(
        groupby(psms, :precursor_idx),
        :prob => (p -> maximum(p)) => :global_prob
    )

    transform!(
        groupby(psms, [:precursor_idx, :ms_file_idx]),
        :prob => (p -> 1.0f0-0.000001f0-exp(sum(log1p.(-p)))) => :prec_prob
    )

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
                    q_value_threshold::Float32 = 0.01f0,
                    max_ftr::Float32 = 0.01f0,
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
        test_fold_filter::Function,
        test_fold_idx::UInt8,
        file_paths::Vector{String},
        bst::Booster,
        features::Vector{Symbol},
        match_between_runs::Bool;
        dropVectorCols::Bool = false,
        save_stage1_prob::Bool = false)

        prec_to_global_score = Dictionary{UInt32, Float32}()
    
        #Reset counts for new scores
        for (key, value) in pairs(prec_to_best_score_new)
            prec_to_best_score_new[key] = (
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
                unique_passing_runs = Set{UInt16}()
            )
        end
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            probs = XGBoost.predict(bst, psms_subset[!,features])
            
            if match_between_runs
                #Update maximum probabilities for tracked precursors 
                qvals = fill(Inf, nrow(psms_subset))
                idxs = [j for j in 1:nrow(psms_subset) if test_fold_filter(psms_subset.cv_fold[j], test_fold_idx)]
                
                if !isempty(idxs)
                    get_qvalues!(view(probs, idxs), view(psms_subset.target, idxs), view(qvals, idxs))
                end


                for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                    if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
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

                        if qvals[i] <= q_value_threshold
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

            # update global prob
            for (i, precursor_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                prob = probs[i]
                if haskey(prec_to_global_score, precursor_idx)
                    prec_to_global_score[precursor_idx] = max(prob, prec_to_global_score[precursor_idx])
                else
                    insert!(prec_to_global_score, precursor_idx, prob)
                end
            end
        end
    
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = XGBoost.predict(bst, psms_subset[!,features])

            for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                psms_subset[i,:prob] = probs[i]
                #if save_stage1_prob
                #    psms_subset[i,:prob_stage1] = probs[i]
                #end
                if match_between_runs
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        scores = prec_to_best_score_new[key]

                        psms_subset.MBR_num_runs[i] = length(scores.unique_passing_runs)

                        best_log2_weights = Float32[]
                        best_irts = Float32[]
                        best_weight = zero(Float32)
                        best_log2_ie = zero(Float32)
                        MBR_is_best_decoy = missing

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
                            psms_subset.MBR_best_irt_diff[i]        = missing
                            psms_subset.MBR_rv_coefficient[i]       = missing
                            psms_subset.MBR_is_best_decoy[i]        = missing
                            psms_subset.MBR_max_pair_prob[i]        = missing
                            psms_subset.MBR_log2_weight_ratio[i]    = missing
                            psms_subset.MBR_log2_explained_ratio[i] = missing
                            psms_subset.MBR_num_runs[i]             = length(scores.unique_passing_runs)
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

            # update global prob
            if !("global_prob" in names(psms_subset))
                psms_subset.global_prob = zeros(Float32, size(psms_subset, 1))
                psms_subset.prec_prob = zeros(Float32, size(psms_subset, 1))
            end

            for (i, precursor_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                if haskey(prec_to_global_score, precursor_idx)
                    psms_subset[i,:global_prob] = prec_to_global_score[precursor_idx]
                end
            end

            transform!(
                groupby(psms_subset, [:precursor_idx]),
                :prob => (p -> 1.0f0-0.000001f0-exp(sum(log1p.(-p)))) => :prec_prob
            )

            if dropVectorCols # last iteration
                Arrow.write(file_path, dropVectorColumns!(psms_subset))
            else
                Arrow.write(file_path, convert_subarrays(psms_subset))
            end
        end
        return prec_to_best_score_new
    end
                    
    function selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold == test_fold
    end

    function excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold != test_fold
    end

    unique_cv_folds = unique(psms[!, :cv_fold])
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{Booster}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    Random.seed!(1776);
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
            _seed_ = rand(UInt32)
            bst = xgboost((psms_train_itr[!,features], psms_train_itr[!,:target]), 
                            num_round=num_round, 
                            #monotone_constraints = monotone_constraints,
                            colsample_bytree = colsample_bytree, 
                            colsample_bynode = colsample_bynode,
                            scale_pos_weight = sum(psms_train_itr.decoy) / sum(psms_train_itr.target),
                            gamma = gamma, 
                            max_depth=max_depth, 
                            eta = eta, 
                            min_child_weight = min_child_weight, 
                            subsample = subsample, 
                            objective="binary:logistic",
                            seed = _seed_,
                            #max_bin = 128,
                            watchlist=(;)
                            )
            if !haskey(models, test_fold_idx)
                insert!(
                    models,
                    test_fold_idx,
                    Vector{Booster}([bst])
                )
            else
                push!(models[test_fold_idx], bst)
            end
            bst.feature_names = [string(x) for x in features]
            #print_importance = true
            print_importance ? println(collect(zip(importance(bst)))[1:30]) : nothing

            # Get probabilities for training sample so we can get q-values
            psms_train[!,:prob] = XGBoost.predict(bst, psms_train[!, features])
            
            if match_between_runs
                summarize_precursors!(psms_train, q_cutoff = max_q_value_xgboost_rescore)
            end

            update(pbar)
        end
    end
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
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
                                                     unique_passing_runs::Set{UInt16},
                                                     }}()

        for (train_iter, num_round) in enumerate(iter_scheme)
            model = models[test_fold_idx][train_iter]
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                selectTestFold,
                test_fold_idx,
                file_paths,
                model,
                features,
                match_between_runs;
                dropVectorCols = (train_iter == length(iter_scheme)) && (test_fold_idx == unique_cv_folds[end]), # remove on last iteration
                save_stage1_prob = (train_iter == length(iter_scheme) - 1)
            )
            update(pbar)
        end
    end

    # sub_psms.MBR_best_irt_diff[i] = missing
    # sub_psms.MBR_rv_coefficient[i] = missing
    # sub_psms.MBR_is_best_decoy[i] = missing
    # sub_psms.MBR_log2_weight_ratio[i] = missing
    # sub_psms.MBR_log2_explained_ratio[i] = missing
    # sub_psms.MBR_max_pair_prob[i] = missing



    
    # if match_between_runs
    #     all_stage1_probs = Float32[]
    #     all_final_probs = Float32[]
    #     all_targets = Bool[]
    #     all_decoys = Bool[]
    #     all_best_decoy = Vector{Union{Missing,Bool}}()

    #     for file_path in file_paths
    #         df = DataFrame(Arrow.Table(file_path))
    #         append!(all_stage1_probs, df.prob_stage1)
    #         append!(all_final_probs, df.prob)
    #         append!(all_targets, df.target)
    #         append!(all_decoys, df.decoy)
    #         append!(all_best_decoy, df.MBR_is_best_decoy)
    #     end

    #     stage1_qvals = Vector{Float64}(undef, length(all_stage1_probs))
    #     get_qvalues!(all_stage1_probs, all_targets, stage1_qvals)
    #     transfer_mask = (stage1_qvals .> q_value_threshold) .& .!ismissing.(all_best_decoy)
    #     bad_transfer = (all_targets .& map(d -> !ismissing(d) && d, all_best_decoy)) .|
    #                    (all_decoys .& map(d -> !ismissing(d) && !d, all_best_decoy))

    #     τ = get_ftr_threshold(all_final_probs, all_targets, bad_transfer, max_ftr;
    #                           mask = transfer_mask)

    #     offset = 1
    #     for file_path in file_paths
    #         df = DataFrame(Arrow.Table(file_path))
    #         n = nrow(df)
    #         idxs = offset:(offset + n - 1)
    #         cand = (stage1_qvals[idxs] .> q_value_threshold) .& .!ismissing.(df.MBR_is_best_decoy)
    #         df.MBR_transfer_candidate .= cand
    #         filter!(row -> !(row.MBR_transfer_candidate && row.prob < τ), df)

    #         transform!(groupby(df, :precursor_idx),
    #                    :prob => (p -> maximum(p)) => :global_prob)

    #         transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
    #                    :prob => (p -> 1.0f0-0.000001f0-exp(sum(log1p.(-p)))) => :prec_prob)

    #         if hasproperty(df, :prob_stage1)
    #             select!(df, Not(:prob_stage1))
    #         end
    #         Arrow.write(file_path, df)
    #         offset += n
    #     end
    # end

    return models
end














function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)   
    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
    Threads.@threads for idx in eachindex(pair_groups)
        _, sub_psms = pair_groups[idx]
        
        # Efficient way to find the top 2 precursors so we can do MBR on the 
        # best precursor match that isn't itself. It's always one of the top 2.

        # single pass: record the best PSM index & prob per run
        best_i = Dict{eltype(sub_psms.ms_file_idx), Int}()
        best_p = Dict{eltype(sub_psms.ms_file_idx), eltype(sub_psms.prob)}()
        for (i, run) in enumerate(sub_psms.ms_file_idx)
            p = sub_psms.prob[i]
            if !haskey(best_p, run) || p > best_p[run]
                best_p[run] = p
                best_i[run] = i
            end
        end

        # if more than one run, find the global top-2 runs by their best-PSM prob
        run_best_indices = Dict{eltype(sub_psms.ms_file_idx), Int}()
        runs = collect(keys(best_i))
        if length(runs) > 1
            # track top two runs (r1 > r2)
            r1 = nothing; p1 = -Inf
            r2 = nothing; p2 = -Inf
            for (run, p) in best_p
                if p > p1
                    r2, p2 = r1, p1
                    r1, p1 = run, p
                elseif p > p2
                    r2, p2 = run, p
                end
            end

            # assign, for each run, the best index in “any other” run
            for run in runs
                run_best_indices[run] = (run == r1 ? best_i[r2] : best_i[r1])
            end
        end

        # Compute MBR features
        for i in 1:nrow(sub_psms)
            sub_psms.MBR_num_runs[i] = length(unique(sub_psms.ms_file_idx[sub_psms.q_value .<= q_cutoff]))

            if !haskey(run_best_indices, sub_psms.ms_file_idx[i])
                sub_psms.MBR_best_irt_diff[i] = missing
                sub_psms.MBR_rv_coefficient[i] = missing
                sub_psms.MBR_is_best_decoy[i] = missing
                sub_psms.MBR_log2_weight_ratio[i] = missing
                sub_psms.MBR_log2_explained_ratio[i] = missing
                sub_psms.MBR_max_pair_prob[i] = missing
                continue
            end

            best_idx = run_best_indices[sub_psms.ms_file_idx[i]]
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
        psms[!, :MBR_max_pair_prob]             = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_best_irt_diff]             = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_log2_weight_ratio]         = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_log2_weight_ratio_min]     = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_log2_explained_ratio]      = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_rv_coefficient]            = Vector{Union{Missing, Float32}}(missing, n)
        psms[!, :MBR_num_runs]                  = zeros(Int32, n)
        psms[!, :MBR_is_best_decoy]             = Vector{Union{Missing, Bool}}(missing, n)
        psms[!, :MBR_transfer_candidate]        = falses(n)
        allowmissing!(psms, [:MBR_max_pair_prob, :MBR_rv_coefficient,
                              :MBR_best_irt_diff, :MBR_is_best_decoy, :MBR_log2_weight_ratio])
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
                [:MBR_is_best_decoy, :MBR_max_pair_prob, :prob] => ByRow((d, mp, p) ->
                    (!ismissing(d) && !d && !ismissing(mp) && !ismissing(p) && mp >= max_prob_threshold && p < max_prob_threshold)
                );
                view = true
            )

            # Compute MBR q-values.
            get_qvalues!(psms_train_mbr[!,:prob], psms_train_mbr[!,:target], psms_train_mbr[!,:q_value])

            # Take all decoys and targets passing q_thresh (all 0's now) or mbr_q_thresh
            psms_train_itr = subset(
                psms_train_itr,
                [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_mbr_rescore))
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