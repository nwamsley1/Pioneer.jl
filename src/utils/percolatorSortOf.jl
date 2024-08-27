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
function getBestScorePerPrecTest!(
    #prec_to_best_score_old::Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    prec_to_best_score_new::Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    cv_fold::UInt8,
    file_paths::Vector{String},
    bst::Booster,
    features::Vector{Symbol})

    println("A")
    #Reset counts for new scores

    for (prec_idx, value) in pairs(prec_to_best_score_new)
        max_prob, mean_prob, min_prob, n = value
        prec_to_best_score_new[prec_idx] = (max_prob = zero(Float32), mean_prob = zero(Float32), min_prob = one(Float32), n = zero(UInt16))
    end
     
    for file_path in file_paths
        psms = DataFrame(Arrow.Table(file_path))
        #Apply old scores to each row 
        #Predict probabilites 
        #based on old scores 
        probs = psms[!,:prob]#XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        test_idx = 1
        m0,m1 = 0,0
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if psms[i,:cv_fold] != cv_fold
                continue
            end
            prob = probs[i]
            if haskey(prec_to_best_score_new, prec_idx)
                m0 += 1
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[prec_idx]
                if max_prob < prob
                    max_prob = prob
                end
                if min_prob > prob
                    min_prob = prob
                end
                mean_prob += prob
                n += one(UInt16)
                test_idx = prec_idx
                prec_to_best_score_new[prec_idx] = (max_prob = max_prob, mean_prob = mean_prob, min_prob = min_prob, n = n)
            else
                m1 += 1
                insert!(prec_to_best_score_new, prec_idx, 
                (max_prob = prob, mean_prob = prob, min_prob = prob, n = one(UInt16)))
            end
        end
    end

    for file_path in file_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        probs = XGBoost.predict(bst, psms[!,features])
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if psms[i,:cv_fold] != cv_fold
                continue
            end
            psms[i,:prob] = probs[i]
            if haskey(prec_to_best_score_new, prec_idx)
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[prec_idx]
                psms[i,:max_prob] = max_prob
                psms[i,:min_prob] = min_prob
                if n > 0
                    psms[i,:mean_prob] = mean_prob/n
                else
                    psms[i,:mean_prob] = zero(Float32)
                end
            end
        end
        Arrow.write(file_path, psms)
    end
    return prec_to_best_score_new
end

function getBestScorePerPrec!(
    #prec_to_best_score_old::Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    prec_to_best_score_new::Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    cv_fold::UInt8,
    file_paths::Vector{String},
    bst::Booster,
    features::Vector{Symbol})

    println("A")
    #Reset counts for new scores

    for (prec_idx, value) in pairs(prec_to_best_score_new)
        max_prob, mean_prob, min_prob, n = value
        prec_to_best_score_new[prec_idx] = (max_prob = zero(Float32), mean_prob = zero(Float32), min_prob = one(Float32), n = zero(UInt16))
    end
     
    for file_path in file_paths
        psms = DataFrame(Arrow.Table(file_path))
        #Apply old scores to each row 
        #Predict probabilites 
        #based on old scores 
        probs = psms[!,:prob]#XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        test_idx = 1
        m0,m1 = 0,0
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if psms[i,:cv_fold] == cv_fold
                continue
            end
            prob = probs[i]
            if haskey(prec_to_best_score_new, prec_idx)
                m0 += 1
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[prec_idx]
                if max_prob < prob
                    max_prob = prob
                end
                if min_prob > prob
                    min_prob = prob
                end
                mean_prob += prob
                n += one(UInt16)
                test_idx = prec_idx
                prec_to_best_score_new[prec_idx] = (max_prob = max_prob, mean_prob = mean_prob, min_prob = min_prob, n = n)
            else
                m1 += 1
                insert!(prec_to_best_score_new, prec_idx, 
                (max_prob = prob, mean_prob = prob, min_prob = prob, n = one(UInt16)))
            end
        end
    end

    for file_path in file_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        probs = XGBoost.predict(bst, psms[!,features])
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if psms[i,:cv_fold] == cv_fold
                continue
            end
            psms[i,:prob] = probs[i]
            if haskey(prec_to_best_score_new, prec_idx)
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[prec_idx]
                psms[i,:max_prob] = max_prob
                psms[i,:min_prob] = min_prob
                if n > 0
                    psms[i,:mean_prob] = mean_prob/n
                else
                    psms[i,:mean_prob] = zero(Float32)
                end
            end
        end
        Arrow.write(file_path, psms)
    end
    return prec_to_best_score_new
end

function getBestScorePerPrec!(
    prec_to_best_score_new::Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    psms::DataFrame,
    bst::Booster,
    features::Vector{Symbol}
)
    println("C")
    #println("is it working?")
    m = 0
    for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
        if haskey(prec_to_best_score_new, prec_idx)
            max_prob, mean_prob, min_prob, n = prec_to_best_score_new[prec_idx]
            psms[i,:max_prob] = max_prob
            psms[i,:min_prob] = min_prob
            psms[i,:mean_prob] = mean_prob/n
            m += 1
        end
    end
    #println("mean(psms[!,:max_prob]) ", mean(psms[!,:max_prob]))
    #println("m $m")
    #psms[!,:train_prob] = XGBoost.predict(bst, psms[!,features])
    #println("m $m")
    #display(describe(psms[!,[:max_prob,:min_prob,:mean_prob,:prob]]))
    #println("test describe ", describe(psms[!,:max_prob]))
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
    #println("TEST")
    println("TEST2")
    psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    folds = psms[!,:cv_fold]
    unique_cv_folds = unique(folds)
    bst = ""
    Random.seed!(128);
    #pbar = ProgressBar(total = sum(iter_scheme)*length(unique_cv_folds))
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    #println("features $features")
    models = Dictionary{UInt8, Vector{Booster}}()
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        psms_train = psms[train_fold_idxs_all,:]
        prec_to_best_score = Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)

            #println("Size of train_fold_idxs for $train_iter iteration for $test_fold_idx fold ,", length(train_fold_idxs))
            ###################
            #Train
            #Train a model on the n-1 training folds. Then apply it to get class probabilities for the test-fold. 
            bst = xgboost((psms_train[!,features], psms_train[!,:target]), 
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
            print_importance ? println(collect(zip(importance(bst)))[1:30]) : nothing
            println(collect(zip(importance(bst)))[1:5])
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                test_fold_idx,
                file_paths,
                bst,
                features
            )
            getBestScorePerPrec!(
                prec_to_best_score,
                psms_train,
                bst,
                features
            )
            println("length(keys(prec_to_best_score_new): ", length(keys(prec_to_best_score)))
        end
    end
    i = 1
    for test_fold_idx in unique_cv_folds
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        psms_train = psms[train_fold_idxs_all,:]
        prec_to_best_score = Dictionary{UInt32, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)
            println("i")
            i += 1
            model = models[test_fold_idx][train_iter]
            prec_to_best_score = getBestScorePerPrecTest!(
                prec_to_best_score,
                test_fold_idx,
                file_paths,
                bst,
                features
            )
        end
    end


    bst.feature_names = [string(x) for x in features]
    println("length(models) ", length(models))
    return models#bst, vcat(psms_test_folds...)
end