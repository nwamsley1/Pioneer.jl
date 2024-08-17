function getQvalues!(probs::Vector{U}, labels::Vector{Bool}, qvals::Vector{T}) where {T,U<:AbstractFloat}

    order = sortperm(probs, rev = true,alg=QuickSort) #Sort class probabilities
    targets = 0
    decoys = 0
    @inbounds @fastmath for i in order
            targets += labels[i]
            decoys += (1 - labels[i])
            qvals[i] = decoys/(targets + decoys)
    end
    #PSMs[:,:q_value] = q_values;
end

getQvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool}) = getQvalues!(PSMs, allowmissing(probs), allowmissing(labels))


#=
function getProbQuantiles!(psms::DataFrame, prob_column::Symbol, features::Vector{Symbol}, bst::Booster)
    psms[!,prob_column] = XGBoost.predict(bst, psms[!,features])
    grouped_psms = groupby(psms, [:precursor_idx,:isotopes_captured]); #
    for i in range(1, length(grouped_psms))
        median_prob = median(grouped_psms[i][!,prob_column])
        max_prob = minimum(grouped_psms[i][!,prob_column])
        q90_prob = quantile(grouped_psms[i][!,prob_column], 0.9)
        for j in range(1, size(grouped_psms[i], 1))
            grouped_psms[i][j,:median_prob] = median_prob
            grouped_psms[i][j,:max_prob] = max_prob
            grouped_psms[i][j,:q90_prob] = q90_prob
        end
    end
    gpsms = groupby(psms, [:file_name, :target, :accession_numbers])
    for i in range(1, length(gpsms))
        prec_count = sum(gpsms[i][!,prob_column].>0.9)
        max_pg_score = maximum(gpsms[i][!,prob_column])
        for j in range(1, size(gpsms[i], 1))
            gpsms[i][j,:pg_count] = prec_count
            gpsms[i][j,:max_pg_score] = max_pg_score
        end
    end 
    gpsms = groupby(psms, [:file_name, :sequence])
    for i in range(1, length(gpsms))
        prec_count = sum(gpsms[i][!,prob_column].>0.9)
        for j in range(1, size(gpsms[i], 1))
            gpsms[i][j,:pepgroup_count] = prec_count
        end
    end 
end
=#

function getBestScorePerPrec!(
    prec_to_best_score_old::Dictionary{UInt32, Float32},
    prec_to_best_score_new::Dictionary{UInt32, Float32},
    cv_fold::UInt8,
    file_paths::Vector{String},
    bst::Booster,
    features::Vector{Symbol}
)
    for file_path in file_paths
        arrow_table = Arrow.Table(file_path)
        #max_probs = zeros(Float32, length(arrow_table[1]))
        psms = DataFrame(Tables.columntable(arrow_table))
        for (i, prec_idx) in enumerate(arrow_table[:precursor_idx])
            if psms[i,:cv_fold] != cv_fold
                continue
            end
            if haskey(prec_to_best_score_old, prec_idx)
                psms[i,:max_prob] = prec_to_best_score_old[prec_idx]
            else
                psms[i,:max_prob] = zero(Float32)
                insert!(prec_to_best_score_old, prec_idx, zero(Float32))
            end
        end
        #Predict probabilites 
        #don't overwrite psms because that will affect other cv folds. 
        probs = XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        for (i, prec_idx) in enumerate(arrow_table[:precursor_idx])
            if psms[i,:cv_fold] != cv_fold
                continue
            end
            prob = probs[i]
            if haskey(prec_to_best_score_new, prec_idx)
                if prec_to_best_score_new[prec_idx] < prob
                     prec_to_best_score_new[prec_idx] = prob
                end
            else
                insert!(prec_to_best_score_new, prec_idx, prob)
            end
            psms[i,:prob] = probs[i]
        end
        #Update table with new max_prob column
        Arrow.write(file_path, psms)
    end
end

function getBestScorePerPrec!(
    prec_to_best_score_old::Dictionary{UInt32, Float32},
    prec_to_best_score_new::Dictionary{UInt32, Float32},
    psms::DataFrame,
    bst::Booster,
    features::Vector{Symbol}
)
    for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
        if haskey(prec_to_best_score_old, prec_idx)
            #max_depthprobs[i] = prec_to_best_score_old[prec_idx]
            psms[i,:max_prob] = prec_to_best_score_old[prec_idx]
        end
    end
    #Predict probabilites 
    probs = XGBoost.predict(bst, psms[!,features])
    #Update maximum probabilities for tracked precursors 
    for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
        if haskey(prec_to_best_score_new, prec_idx)
            if prec_to_best_score_new[prec_idx] < probs[i]
                prec_to_best_score_new[prec_idx] = probs[i]
            end
        end
    end
end



function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; 
                    colsample_bytree::Float64 = 0.5, 
                    colsample_bynode::Float64 = 0.5,
                    eta::Float64 = 0.15, 
                    min_child_weight::Int = 1, 
                    subsample::Float64 = 0.5, 
                    gamma::Int = 0, 
                    max_depth::Int = 10,
                    train_fraction::Float64 = 1.0,
                    iter_scheme::Vector{Int} = [100, 100, 200], 
                    print_importance::Bool = false)
    PSMs[!,:max_prob], PSMs[!,:median_prob], PSMs[!,:q90_prob] = zeros(Float32, size(PSMs)[1]), zeros(Float32, size(PSMs)[1]), zeros(Float32, size(PSMs)[1])
    PSMs[!,:prob_temp],PSMs[!,:prob] = zeros(Float32, size(PSMs)[1]),zeros(Float32, size(PSMs)[1])
    PSMs[!,:max_pg_score] = zeros(Float32, size(PSMs, 1))
    PSMs[!,:pg_count] = zeros(UInt16, size(PSMs, 1))
    PSMs[!,:pepgroup_count] = zeros(UInt16, size(PSMs, 1))
    folds = PSMs[!,:cv_fold]
    unique_cv_folds = unique(folds)
    bst = ""
    Random.seed!(128);
    #pbar = ProgressBar(total = sum(iter_scheme)*length(unique_cv_folds))
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    psms_test_folds = []
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        PSMs[!,:max_prob] .= zero(Float64)
        PSMs[!,:median_prob] .= zero(Float64)
        #row idxs for a
        test_fold_idxs = findall(x->x==test_fold_idx, folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        train_fold_idxs_sub = Random.randsubseq(train_fold_idxs_all[:], train_fraction)
        #psms_train = PSMs[train_fold_idxs_all,:]
        psms_train_sub = PSMs[train_fold_idxs_sub,:]
        psms_test_sub = PSMs[test_fold_idxs,:]

        for (train_iter, num_round) in enumerate(iter_scheme)
            #println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            ###################
            #Train
            #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
            bst = xgboost((psms_train_sub[!,features], psms_train_sub[!,:target]), 
                            num_round=num_round, 
                            #monotone_constraints = monotone_constraints,
                            colsample_bytree = colsample_bytree, 
                            colsample_bynode = colsample_bynode,
                            gamma = gamma, 
                            max_depth=max_depth, 
                            eta = eta, 
                            min_child_weight = min_child_weight, 
                            subsample = subsample, 
                            objective="binary:logistic",
                            seed = rand(UInt32),
                            #max_bin = 128,
                            watchlist=(;)
                            )
            bst.feature_names = [string(x) for x in features]
            print_importance ? println(collect(zip(importance(bst)))[1:30]) : nothing
            ###################
            #Apply to CV fold 
            getProbQuantiles!(psms_train_sub, :prob_temp, features, bst)
            #Apply to held-out data
            getProbQuantiles!(psms_test_sub, :prob, features, bst)
            #update(pbar, num_round)
        end
        push!(psms_test_folds, psms_test_sub)
  
    end
    bst.feature_names = [string(x) for x in features]
    return bst, vcat(psms_test_folds...)
end

function rankPSMs!(psms::DataFrame, 
                    file_paths::Vector{String},
                    features::Vector{Symbol}; 
                    colsample_bytree::Float64 = 0.5, 
                    colsample_bynode::Float64 = 0.5,
                    eta::Float64 = 0.15, 
                    min_child_weight::Int = 1, 
                    subsample::Float64 = 0.5, 
                    gamma::Int = 0, 
                    max_depth::Int = 10,
                    train_fraction::Float64 = 1.0,
                    iter_scheme::Vector{Int} = [100, 100, 200], 
                    print_importance::Bool = true)

    psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    folds = psms[!,:cv_fold]
    unique_cv_folds = unique(folds)
    bst = ""
    Random.seed!(128);
    #pbar = ProgressBar(total = sum(iter_scheme)*length(unique_cv_folds))
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    psms_test_folds = []
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        cv_psms = psms[train_fold_idxs_all,:]
        precs_to_score = unique(cv_psms[!,:precursor_idx])
        prec_to_best_score_new = Dictionary(
                                precs_to_score,
                                zeros(Float32, length(precs_to_score))
                            )
        prec_to_best_score_old = copy(prec_to_best_score_new)

        for (train_iter, num_round) in enumerate(iter_scheme)
            #println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            ###################
            #Train
            #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
            bst = xgboost((cv_psms[!,features], cv_psms[!,:target]), 
                            num_round=num_round, 
                            #monotone_constraints = monotone_constraints,
                            colsample_bytree = colsample_bytree, 
                            colsample_bynode = colsample_bynode,
                            gamma = gamma, 
                            max_depth=max_depth, 
                            eta = eta, 
                            min_child_weight = min_child_weight, 
                            subsample = subsample, 
                            objective="binary:logistic",
                            seed = rand(UInt32),
                            #max_bin = 128,
                            watchlist=(;)
                            )
            bst.feature_names = [string(x) for x in features]
            print_importance ? println(collect(zip(importance(bst)))[1:30]) : nothing

            ###################
            #Apply to held out data on disc. Get best/mean/min score accross
            #all runs for each precursor in this CV fold. 
            @time begin
            getBestScorePerPrec!(
                prec_to_best_score_old,
                prec_to_best_score_new,
                cv_psms,
                bst,
                features
            )

            getBestScorePerPrec!(
                prec_to_best_score_old,
                prec_to_best_score_new,
                test_fold_idx,
                file_paths,
                bst,
                features
            )
            
            prec_to_best_score_old = copy(prec_to_best_score_new)
            end
        end
        #push!(psms_test_folds, psms_test_sub)
  
    end
    bst.feature_names = [string(x) for x in features]
    return bst#bst, vcat(psms_test_folds...)
end