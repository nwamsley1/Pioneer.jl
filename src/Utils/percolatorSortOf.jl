#=
function getQvalues!(PSMs::DataFrame, probs::Vector{Union{Missing, Float64}}, labels::Vector{Union{Missing, Bool}})
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (length(probs),))
    order = reverse(sortperm(probs)) #Sort class probabilities
    targets = 0
    decoys = 0
    for i in order
        if labels[i] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_value] = q_values;
end
=#
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
    folds = best_psms[!,:cv_fold]
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
            #Apply to held-out data
            getProbQuantiles!(psms_train_sub, :prob_temp, features, bst)
            getProbQuantiles!(psms_test_sub, :prob, features, bst)
            #update(pbar, num_round)
        end
        push!(psms_test_folds, psms_test_sub)
  
    end
    bst.feature_names = [string(x) for x in features]
    return bst, vcat(psms_test_folds...)
end