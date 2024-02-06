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
function getQvalues!(probs::Vector{T}, labels::Vector{Bool}, qvals::Vector{T}) where {T<:AbstractFloat}

    order = sortperm(probs, rev = true,alg=QuickSort) #Sort class probabilities
    targets = 0
    decoys = 0
    @inbounds @fastmath for i in order
            decoys += labels[i]
            targets += (1 - labels[i])
            qvals[i] = decoys/(targets + decoys)
    end
    #PSMs[:,:q_value] = q_values;
end

getQvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool}) = getQvalues!(PSMs, allowmissing(probs), allowmissing(labels))

function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, train_fraction::Float64 = 1.0, n_iters::Int = 2, print_importance::Bool = false)

    #Initialize class probabilisites
    PSMs[:,:prob_temp] = zeros(Float64, size(PSMs)[1])
    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    #Sort by sequence. All psms corresponding to a given sequence
    #should be in the same cross-validation fold
    sort!(PSMs, :sequence);
     
    #Data structure for training. 
    #Must initialize AFTER sorting. 
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]
    
    #Assign psms to cross-validation folds
    folds = ones(Int64, size(PSMs)[1])
    fold = 1
    fold_ids = Tuple([_ for _ in range(1, n_folds)])
    Random.seed!(1234)
    for i in range(2, length(folds))
        if PSMs[i,:sequence] == PSMs[i-1,:sequence]
            folds[i] = fold
        else
            fold = rand(fold_ids)
            folds[i] = fold
        end
    end
    #XGBoost model
    println(unique(folds))
    bst = ""
    PSMs[:,:K_FOLD] .= 0
    
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    for test_fold_idx in range(1, n_folds)
        
        #row idxs for a
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        train_fold_idxs = Random.randsubseq(train_fold_idxs_all[:], train_fraction)
        #Iterative training 
        #for train_iter in range(1, n_iters)
        for train_iter in range(1, 5)
            println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            ###################
            #Train
            if train_iter != 5
                num_round = 50
            else
                num_round = 200
            end
            train_features = X[train_fold_idxs,:]
            train_classes = X_labels[train_fold_idxs,1]
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
                            #watchlist=(;)
                            )

            ###################
            #Apply to held-out data
            ŷ = XGBoost.predict(bst, X[train_fold_idxs_all,:])
            for i in eachindex(train_fold_idxs_all)
                PSMs[train_fold_idxs_all[i],:prob_temp] = (1 - ŷ[i])
            end
            ###################
            #Get Optimized Training Set
            order = reverse(sortperm(PSMs[train_fold_idxs_all,:prob_temp]))
            targets = 0
            decoys = 0
            new_train_fold_idxs = Vector{Int64}()
            for i in order
                if PSMs[train_fold_idxs_all[i],:decoy] == true
                    decoys += 1
                    if rand() <= train_fraction
                        push!(new_train_fold_idxs, train_fold_idxs_all[i])
                    end
                else
                    targets += 1
                    if decoys/(targets + decoys) <= 0.01
                        if rand() <= train_fraction
                            push!(new_train_fold_idxs, train_fold_idxs_all[i])
                        end
                    end
                end
            end
            train_fold_idxs = new_train_fold_idxs
        end
        test_fold_idxs = findall(x->x==test_fold_idx, folds)
        ŷ = XGBoost.predict(bst, X[test_fold_idxs,:])
        for i in eachindex(test_fold_idxs)
            PSMs[test_fold_idxs[i],:prob] = (1 - ŷ[i])
            PSMs[test_fold_idxs[i],:K_FOLD] = test_fold_idxs[i]
        end
    end
    bst.feature_names = [string(x) for x in features]
    return bst, folds
end


function rankPSMs2!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, train_fraction::Float64 = 1.0, n_iters::Int = 2, print_importance::Bool = false)
    #if Symbol("max_prob") ∉ features
    #    push!(features, Symbol("max_prob"))
    #end
    println("urmom")
    PSMs[!,:row] = collect(range(1, size(PSMs, 1)))
    PSMs[!,:max_prob] = zeros(Float64, size(PSMs)[1])
    #Initialize class probabilisites
    PSMs[:,:prob_temp] = zeros(Float64, size(PSMs)[1])
    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    #Sort by sequence. All psms corresponding to a given sequence
    #should be in the same cross-validation fold
    sort!(PSMs, :sequence);
    #Data structure for training. 
    #Must initialize AFTER sorting. 
    #X = Matrix(PSMs[:,features])
    #X_labels = PSMs[:, :decoy]
    #Assign psms to cross-validation folds
    folds = ones(Int64, size(PSMs)[1])
    fold = 1
    fold_ids = Tuple([_ for _ in range(1, n_folds)])
    Random.seed!(1234)
    for i in range(2, length(folds))
        if PSMs[i,:sequence] == PSMs[i-1,:sequence]
            folds[i] = fold
        else
            fold = rand(fold_ids)
            folds[i] = fold
        end
    end
    #XGBoost model
    println(unique(folds))
    bst = ""
    PSMs[:,:K_FOLD] .= 0
    println("features $features")
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    psms_test_folds = []
    for test_fold_idx in range(1, n_folds)
        
        #row idxs for a
        test_fold_idxs = findall(x->x==test_fold_idx, folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        train_fold_idxs = Random.randsubseq(train_fold_idxs_all[:], train_fraction)
        #Iterative training 
        train_features = nothing
        train_classes = nothing
        psms_train_sub = PSMs[train_fold_idxs_all,:]
        psms_test_sub = PSMs[test_fold_idxs,:]
        prec_to_prob = Dictionary{UInt32, Float32}()
        #for train_iter in range(1, n_iters)
        #    println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
        for train_iter in range(1, 3)
            println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            if train_iter != 3
                num_round = 100
            else
                num_round = 200
            end
            ###################
            #Train
            println("sum probs ", sum(psms_train_sub[:,:max_prob]))
            train_features = Matrix(psms_train_sub[:,features])#X[train_fold_idxs,:]
            train_classes = psms_train_sub[:,:decoy]
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
                            #watchlist=(;)
                            )

            ###################
            #Apply to held-out data
            ŷ = XGBoost.predict(bst, train_features)
            for i in range(1, size(psms_train_sub, 1))#eachindex(train_fold_idxs_all)
                psms_train_sub[i,:prob_temp] = (1 - ŷ[i])
            end
            println("start grouping")
            grouped_psms = groupby(psms_train_sub, :precursor_idx); #
            println("stop grouping")
            for i in ProgressBar(range(1, length(grouped_psms)))
                median_prob = median(grouped_psms[i][!,:prob_temp])
                if train_iter == 1
                    insert!(prec_to_prob, first(grouped_psms[i].precursor_idx), median_prob)
                end
                for j in range(1, size(grouped_psms[i], 1))
                    grouped_psms[i][j,:max_prob] = median_prob
                end
            end

            ŷ = XGBoost.predict(bst, psms_test_sub[:,features])
            for i in range(1, size(psms_test_sub, 1))
                psms_test_sub[i,:prob] = (1 - ŷ[i])
                #PSMs[test_fold_idxs[i],:K_FOLD] = test_fold_idxs[i]
            end

            println("start grouping")
            grouped_psms = groupby(psms_test_sub, :precursor_idx); #
            println("stop grouping")
            for i in ProgressBar(range(1, length(grouped_psms)))
                median_prob = median(grouped_psms[i][!,:prob])
                if train_iter == 1
                    insert!(prec_to_prob, first(grouped_psms[i].precursor_idx), median_prob)
                end
                for j in range(1, size(grouped_psms[i], 1))
                    grouped_psms[i][j,:max_prob] = median_prob
                end
            end

            ###################
            #Get Optimized Training Set
            #=
            order = reverse(sortperm(PSMs[train_fold_idxs_all,:prob_temp]))
            targets = 0
            decoys = 0
            new_train_fold_idxs = Vector{Int64}()
            for i in order
                if PSMs[train_fold_idxs_all[i],:decoy] == true
                    decoys += 1
                    if rand() <= train_fraction
                        push!(new_train_fold_idxs, train_fold_idxs_all[i])
                    end
                else
                    targets += 1
                    if decoys/(targets + decoys) <= 0.01
                        if rand() <= train_fraction
                            push!(new_train_fold_idxs, train_fold_idxs_all[i])
                        end
                    end
                end
            end
            train_fold_idxs = new_train_fold_idxs
            =#
        end
        push!(psms_test_folds, psms_test_sub)
  
    end
    bst.feature_names = [string(x) for x in features]
    return bst, folds, vcat(psms_test_folds...)
end

