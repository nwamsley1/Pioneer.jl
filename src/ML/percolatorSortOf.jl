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

function rankPSMs2!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, train_fraction::Float64 = 1.0, n_iters::Int = 2, print_importance::Bool = false)
    

    PSMs[!,:row] = collect(range(1, size(PSMs, 1)))
    PSMs[!,:max_prob] = zeros(Float64, size(PSMs)[1]) #Maximum probability observed per precursor trace
    PSMs[:,:prob_temp] = zeros(Float64, size(PSMs)[1]) #
    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    bst = ""
    println("features $features")
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    psms_test_folds = []
    for test_fold_idx in unique(PSMs[!,:cv_fold])
        
        #row idxs for a
        test_fold_idxs = findall(x->x==test_fold_idx, PSMs[!,:cv_fold])
        train_fold_idxs_all = findall(x->x!=test_fold_idx, PSMs[!,:cv_fold])
        train_fold_idxs = Random.randsubseq(train_fold_idxs_all[:], train_fraction)
        #Iterative training 
        train_features = nothing
        train_classes = nothing
        psms_train_sub = PSMs[train_fold_idxs_all,:]
        psms_test_sub = PSMs[test_fold_idxs,:]
        for train_iter in range(1, 3)
            println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            if train_iter != 3
                num_round = 100
            else
                num_round = 200
            end
            ###################
            #Train
            train_features = Matrix(psms_train_sub[:,features])#X[train_fold_idxs,:]
            train_classes = psms_train_sub[:,:target]
            #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
            bst = xgboost((train_features, train_classes), 
                            num_round=num_round, 
                            colsample_bytree = colsample_bytree, 
                            gamma = gamma, 
                            max_depth=max_depth, 
                            eta = eta, 
                            min_child_weight = min_child_weight, 
                            subsample = subsample, 
                            objective="binary:logistic",
                            watchlist=(;)
                            )
            ###################
            #Apply to held-out data
            psms_train_sub[i,:prob_temp] = XGBoost.predict(bst, train_features)
            grouped_psms = groupby(psms_train_sub, [:precursor_idx,:iso_rank]); 
            for i in ProgressBar(range(1, length(grouped_psms)))
                median_prob = median(grouped_psms[i][!,:prob_temp])
                for j in range(1, size(grouped_psms[i], 1))
                    grouped_psms[i][j,:max_prob] = median_prob
                end
            end

            psms_test_sub[i,:prob]  = XGBoost.predict(bst, psms_test_sub[:,features])
            grouped_psms = groupby(psms_test_sub, [:precursor_idx,:iso_rank]); #
            for i in ProgressBar(range(1, length(grouped_psms)))
                median_prob = median(grouped_psms[i][!,:prob])
                for j in range(1, size(grouped_psms[i], 1))
                    grouped_psms[i][j,:max_prob] = median_prob
                end
            end
        end
        push!(psms_test_folds, psms_test_sub)
  
    end
    bst.feature_names = [string(x) for x in features]
    return bst, folds, vcat(psms_test_folds...)
end

