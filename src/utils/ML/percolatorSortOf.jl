function getBestScorePerPrec!(
    prec_to_best_score_new::Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    test_fold_filter::Function,
    test_fold_idx::UInt8,
    file_paths::Vector{String},
    bst::Booster,
    features::Vector{Symbol})

    #Reset counts for new scores
    for (key, value) in pairs(prec_to_best_score_new)
        max_prob, mean_prob, min_prob, n = value
        prec_to_best_score_new[key] = (max_prob = zero(Float32), mean_prob = zero(Float32), min_prob = one(Float32), n = zero(UInt16))
    end
     
    for file_path in file_paths
        psms = DataFrame(Arrow.Table(file_path))
        probs = XGBoost.predict(bst, psms[!,features])
        #Update maximum probabilities for tracked precursors 
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if !test_fold_filter(psms[i,:cv_fold], test_fold_idx)::Bool continue end
            prob = probs[i]
            key = (prec_idx = prec_idx, isotopes = psms[i,:isotopes_captured])
            if haskey(prec_to_best_score_new, key)
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
                if max_prob < prob
                    max_prob = prob
                end
                if min_prob > prob
                    min_prob = prob
                end
                mean_prob += prob
                n += one(UInt16)
                prec_to_best_score_new[key] = (max_prob = max_prob, mean_prob = mean_prob, min_prob = min_prob, n = n)
            else
                insert!(prec_to_best_score_new, key, 
                (max_prob = prob, mean_prob = prob, min_prob = prob, n = one(UInt16)))
            end
        end
    end

    for file_path in file_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(file_path)))
        probs = XGBoost.predict(bst, psms[!,features])
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            if !test_fold_filter(psms[i,:cv_fold], test_fold_idx)::Bool continue end
            psms[i,:prob] = probs[i]
            key = (prec_idx = prec_idx, isotopes = psms[i,:isotopes_captured])
            if haskey(prec_to_best_score_new, key)
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
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
    prec_to_best_score_new::Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
    psms::DataFrame)
    m = 0
    for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
        key = (prec_idx = prec_idx, isotopes = psms[i,:isotopes_captured])
        if haskey(prec_to_best_score_new, key)
            max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
            psms[i,:max_prob] = max_prob
            psms[i,:min_prob] = min_prob
            psms[i,:mean_prob] = mean_prob/n
            m += 1
        end
    end
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
                    iter_scheme::Vector{Int} = [100, 100, 200], 
                    print_importance::Bool = true)

    function selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold == test_fold
    end

    function excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold != test_fold
    end

    psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    folds = psms[!,:cv_fold]
    unique_cv_folds = unique(folds)
    Random.seed!(128);
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{Booster}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        psms_train = psms[train_fold_idxs_all,:]
        prec_to_best_score = Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)
            ###################
            #Train a model on the n-1 training folds.
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
            #println(collect(zip(importance(bst)))[1:5])
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                excludeTestFold,
                test_fold_idx,
                file_paths,
                bst,
                features
            )
            getBestScorePerPrec!(
                prec_to_best_score,
                psms_train
            )
            update(pbar)
        end
    end
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        prec_to_best_score = Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)
            model = models[test_fold_idx][train_iter]
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                selectTestFold,
                test_fold_idx,
                file_paths,
                model,
                features
            )
            update(pbar)
        end
    end
    
    return models#bst, vcat(psms_test_folds...)
end